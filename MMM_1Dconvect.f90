program MMM_1Dconvect
   use, intrinsic :: iso_fortran_env, only: real64
   use MMM_module_TimeINT 
   use MMM_module_spatial 
   use module_SchemeSetup
   use VTK_MultiMoment_Writer
   use module_physics
   
   use matlab_io_mod, only:save_bin_matlab
!!   use froto_vortex_functions, only: get_myU
   
   
   implicit none

   type(MeshDomain) :: Grid01               !Physical Grid
   type(scalar_spatial_data) :: Phi, Phi_init	    !Scalar Phi
   type(scalar_spatial_data) :: DUDX, DUDY, DVDX, DVDY   ! VELOCITY GRADIENT
   type(velocity_spatial_data) :: Velocity  !Velocity 
   type(TimeIntArray) :: TMP                !Work array
   type(Coefficient_Differentiation)  :: COEF   
    
   integer :: ncell, n_steps, npx, npy, unit, ios, iMesh, nx, ny, i
   real(wp) :: l, dx, dy, u, t_final, cfl,cfl0, l0, l1, dt, rho, c0, c1,X0,Y0,XL,YL
   real(wp) :: start_time, end_time, solve_time, pv_linf, cv_linf, l2_ratioPV, l2_ratioCV
   real(wp), parameter :: pi = acos(-1.0_wp)
   character(len=16) :: scase, method, scheme, METHNCs,    sSptScheme, sTimeINTmethod, sNcellX, sCFL


   character(len=256) :: filename
   integer :: NcellX, NcellY,  j , istart,iend,  NPxtot, NPytot
   
   character(len=256) :: char_buffer   
 
    REAL(wp) :: iDX, iDY, TimeNow 
    REAL(wp) :: FpERR=0.0_wp, FyzERR=0.0_wp, FxzERR=0.0_wp, FxyzERR=0.0_wp
    integer :: file_unit
   
   REAL(wp) :: T0=0.0 !initial time
 

 
   !read input file
   call read_setup('setup.txt', ng, ncell, cfl, u, Kdiff, t_final, scase, method, scheme, METHNCs)   


 
 
   NcellX = ncell
   NcellY = ncell
   NpX =  NcellX + 1
   NpY =  NcellY + 1
   
   
   
   iMesh = 1
   X0 = -5.0_wp ; XL = 5.0_wp
   Y0 = -5.0_wp ; YL = 5.0_wp


   !setup the method and other fundamental variables 
   call set_TimeINTmethod(  method, TimeINTmethod   )
   call set_MMMscheme(  scheme, METHNCs, SptScheme , isMM, METHNC   ) 
   !report setup of the case 
   print *, " ng     ncell          cfl      u           Kdiff      t_final  scase     method         scheme          isMM"   
   print '(I4, I10,E13.5, F11.3, E12.5, F12.3, XXX,A10,XXX,A10,XXX,A10, L9)',   &
   ng, ncell, cfl, u, Kdiff, t_final, scase, method, scheme, isMM   

	call setBandWidth_coeff(ng)
  !Initialize data	
	call Grid01%init( iMesh, ng, NpX, NpY, X0, XL, Y0, YL)
	call  PHI%init( Grid01 )
	call DUDX%init( Grid01 ) 
	call DUDY%init( Grid01 ) 
	call DVDX%init( Grid01 ) 
	call DVDY%init( Grid01 )	
	call Velocity%init( Grid01 )
	call TMP%init( TimeINTmethod, Grid01 )	
	call COEF%init( Grid01 )


! Initialize Velocity and Scalar field
        call get_fronto_final_phi( PHI, T0 ) 
	call get_fronto_initial_Vel( Velocity )


	DX = Grid01%DX
	DY = Grid01%DY	
	iDX = 1.0_wp / DX
	iDY = 1.0_wp / DY

	print*,"DX = ", DX

	PHI_init = PHI
	NPxtot= Grid01%NPxtot
	NPytot= Grid01%NPytot
     
associate( 	U => Velocity%U, V => Velocity%V )


!!!!!!!!!!!!!!!!!!!!!!!  VELOCITY GRADIENT FOR FDM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  IF( .NOT. isMM  ) THEN
  !------------------dbg  
 	 CALL Compute_FirstDerivativeX_FD( DUDX%p, U%p, Grid01, COEF, NPxtot, NPytot ) 
	 CALL Compute_FirstDerivativeY_FD( DVDY%p, V%p, Grid01, COEF, NPxtot, NPytot )
 

 	 CALL Compute_FirstDerivativeX_FD( DUDX%p, PHI%p, Grid01, COEF, NPxtot, NPytot ) 
	 CALL Compute_FirstDerivativeY_FD( DUDY%p, PHI%p, Grid01, COEF, NPxtot, NPytot )

!! DIAGNOSTIC CODE
!!	call save_paraview_point_binary_2d('DFDXp.vts', Grid01, DUDX, 'DFDXp')
!!	call save_paraview_point_binary_2d('DFDYp.vts', Grid01, DUDY, 'DFDYp')  
 
  
 	 CALL Compute_FirstDerivativeX_FD( DUDX%p, U%p, Grid01, COEF, NPxtot, NPytot ) 
	 CALL Compute_FirstDerivativeY_FD( DUDY%p, U%p, Grid01, COEF, NPxtot, NPytot )   
	 CALL Compute_FirstDerivativeX_FD( DVDX%p, V%p, Grid01, COEF, NPxtot, NPytot ) 	 
	 CALL Compute_FirstDerivativeY_FD( DVDY%p, V%p, Grid01, COEF, NPxtot, NPytot )  
  ELSE
!!!!!!!!!!!!!!!!!!!!!!!  VELOCITY GRADIENT FOR MMM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1  
	 CALL Compute_FirstDerivativeX_YZ( DUDX%p, U%xz, U%p ,Grid01, COEF, NPxtot, NPytot )
	 CALL Compute_FirstDerivativeX_YZ( DVDX%p, V%xz, V%p ,Grid01, COEF, NPxtot, NPytot ) 
	 
	 CALL Compute_FirstDerivativeY_XZ( DUDY%p, U%yz, U%p ,Grid01, COEF, NPxtot, NPytot ) 		 
	 CALL Compute_FirstDerivativeY_XZ( DVDY%p, V%yz, V%p ,Grid01, COEF, NPxtot, NPytot ) 		 
	 
	 	 
  	 CALL Compute_FirstDerivativeX_YZ( DUDX%yz, U%xyz, U%yz ,Grid01, COEF, NPxtot, NPytot -1 )  	 
  	 CALL Compute_FirstDerivativeX_YZ( DVDX%yz, V%xyz, V%yz ,Grid01, COEF, NPxtot, NPytot -1 )  	


	 CALL Compute_FirstDerivativeY_XZ( DVDY%xz, V%xyz, V%xz ,Grid01, COEF, NPxtot -1, NPytot ) 	  
	 CALL Compute_FirstDerivativeY_XZ( DUDY%xz, U%xyz, U%xz ,Grid01, COEF, NPxtot -1, NPytot ) 	  	 

 
!!!! ANALYTICAL CELL- and LINE-AVERAGE	 	 	 
     DUDX%xyz(:,1:(Grid01%IfEND-1) )=( Velocity%U%yz(:, 2:Grid01%IfEND) - Velocity%U%yz(  :,  1:(Grid01%IfEND - 1) ) )*iDX
     DVDX%xyz(:,1:(Grid01%IfEND-1))=( Velocity%V%yz(:, 2:Grid01%IfEND  ) - Velocity%V%yz(  :,  1:Grid01%IfEND - 1 ) )*iDX
     
     DVDY%xyz(1:(Grid01%JfEND - 1),:)=( Velocity%V%xz(   2:Grid01%JfEND, : ) - Velocity%V%xz(   1:Grid01%JfEND - 1, :) )*iDY 
     DUDY%xyz(1:(Grid01%JfEND - 1),:)=( Velocity%U%xz(   2:Grid01%JfEND, : ) - Velocity%U%xz(   1:Grid01%JfEND - 1, :) )*iDY
 
     
     DUDX%xz(:, 1:(Grid01%IfEND - 1))=( Velocity%U%p(:, 2:Grid01%IfEND  ) - Velocity%U%p(  :,  1:(Grid01%IfEND - 1) ) )*iDX     
     DVDX%xz(:, 1:(Grid01%IfEND - 1))=( Velocity%V%p(:, 2:Grid01%IfEND  ) - Velocity%V%p(  :,  1:(Grid01%IfEND - 1) ) )*iDX        
    

     DVDY%yz(1:(Grid01%JfEND - 1),:)=( Velocity%V%p(   2:Grid01%JfEND, : ) - Velocity%V%p(   1:Grid01%JfEND - 1, :) )*iDY 
     DUDY%yz(1:(Grid01%JfEND - 1),:)=( Velocity%U%p(   2:Grid01%JfEND, : ) - Velocity%U%p(   1:Grid01%JfEND - 1, :) )*iDY     
 

  ENDIF	 
end associate  

 
 
!!CALL save_paraview_staggered_2d( 'DUDXall' , Grid01, DUDX )
!!CALL save_paraview_staggered_2d( 'DVDXall' , Grid01, DVDY )


 
!call Write_MultiMoment_VTM( Velocity%U%p, Velocity%U%yz, Velocity%U%xz, Velocity%U%xyz, &
!                            Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'U' )

if( isMM) then
	call Write_MultiMoment_VTM( Velocity%U%p, Velocity%U%yz, Velocity%U%xz, Velocity%U%xyz, &
		                        Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'U' )
		                        
	call Write_MultiMoment_VTM( DUDX%p, DUDX%yz, DUDX%xz, DUDX%xyz, &
		                        Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'DUDX' )         
		                        
	call Write_MultiMoment_VTM( DVDY%p, DVDY%yz, DVDY%xz, DVDY%xyz, &
		                      Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'DVDY' )                            
endif  

call save_paraview_point_binary_2d('Fp_initial.vts', Grid01, PHI, 'PHI')     
call save_paraview_point_binary_2d('Up_initial.vts', Grid01, Velocity%U, 'U')
call save_paraview_point_binary_2d('Vp_initial.vts', Grid01, Velocity%V, 'V')




!call save_paraview_point_binary_2d('DUDX.vts', Grid01, Velocity%U, 'DUDX%p')
!call save_paraview_point_binary_2d('DUDY.vts', Grid01, Velocity%V, 'DUDY%p')




istart = Velocity%Mesh%IcSTART
iend   = Velocity%Mesh%IcEND

IF( isMM) THEN
	call save_paraview_cell_2d('PHIxyz_initial.vts', Grid01, Velocity%U%xyz, 'PHIxyz')
	call save_paraview_cell_2d('Uxyz_initial.vts', Grid01, Velocity%U%xyz, 'Uxyz')
	call save_paraview_cell_2d('Vxyz_initial.vts', Grid01, Velocity%U%xyz, 'Vxyz')
 
ENDIF

 
   


 


 

 
 

   if (u > 1e-10) then
              dt = cfl * (dx / u)
   else
              dt = cfl * (dx / 1.0_wp )	    
   endif

 
 

   if( dt .GT. t_final) DT=1.0_wp*t_final

   n_steps = int(t_final / dt)
   dt = t_final / real(n_steps, wp)
    

 

   cfl0 = cfl 
   TMP%DT = DT



 
!CALL set_fronto_timedep_BC_simp( PHI, t_final )

!call save_paraview_point_binary_2d_w_GC('PHIp_fin_ghost.vts', Grid01, PHI, 'PHIgh') 

!	call Write_MultiMoment_VTM( PHI%p, PHI%yz, PHI%xz, PHI%xyz, &
!		                        Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'FwGH' )
! 
! 
!stop "check ghostcells"












TimeNOW = T0
!--------------- Time Integration loop --------------------------------------------------
   call cpu_time(start_time)
   i = 1
   do i = 1, n_steps

      call TimeIntegration(   Grid01, Velocity, DUDX, DUDY, DVDX, DVDY, PHI , TMP, COEF , TimeNOW )
      TimeNOW = TimeNOW + DT      
     IF( mod(i, n_steps/10) .EQ.0) print *,   i, " / ", n_steps
   end do
   call cpu_time(end_time)
   solve_time = end_time - start_time
   
   
   
!      write(6,*) "  t_final ", t_final   
!      write(6,*) " Final Time step ,  DT,  Total number of Time step, Final Time", i-1, dt, n_steps, t_final
!      write(6,*) " SOLUTION TIME = ", solve_time
      
CALL set_fronto_timedep_BC_simp( PHI, t_final  + T0  )      
      
!--------------- Time Integration loop ----------------------------------------------------
call save_paraview_point_binary_2d('PHIp_final.vts', Grid01, PHI, 'PHI')    
call get_fronto_final_phi( Phi_init, t_final  + T0 ) 
call save_paraview_point_binary_2d('PHIp_exact_final.vts', Grid01, Phi_init, 'PHIexact') 

 


 FpERR =  maxval (  abs(Phi_init%p( istart:iend,istart:iend) - &
                        PHI%p( istart:iend,istart:iend) ) )

  

if( isMM) then

 FyzERR =  maxval (  abs(Phi_init%yz( istart:iend,istart:iend) - &
                        PHI%yz( istart:iend,istart:iend) ) )

 FxzERR =  maxval (  abs(Phi_init%xz( istart:iend,istart:iend) - &
                        PHI%xz( istart:iend,istart:iend) ) )

 FxyzERR = maxval (  abs(Phi_init%xyz( istart:iend,istart:iend) - &
                         PHI%xyz( istart:iend,istart:iend) ) ) 
                         
 
call Write_MultiMoment_VTM( PHI%p, PHI%yz, PHI%xz, PHI%xyz, &
		                        Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'PHI_final' )
 
call Write_MultiMoment_VTM( Phi_init%p, Phi_init%yz, Phi_init%xz, Phi_init%xyz, &
		                        Grid01%Xp, Grid01%Yp, Grid01%Xc, Grid01%Yc, 'PHI_exact' )

endif

 
write(sCFL,'(F7.5)' ) CFL 


 write(filename, '(A,I0,A,I0,A,I0,A,I0,A,A,A)') &
    "REPORT_SCHM_", SptScheme, &
    "_NC_", METHNC, &
    "_TI_", TimeINTmethod, &
    "_NCELL_", NcellX, &
    "_CFL_", adjustl(trim(sCFL)), ".txt"
		 
		 

 


open(newunit=file_unit, file=filename, status="replace", action="write")

 write(6,'(A)') "****************************************************************************************************************"  
 write(6,'(A,4I6,F9.4)') "Scheme TimeINTmethod METHNC NcellX CFL ::", SptScheme ,TimeINTmethod, METHNC, NcellX, CFL
 write(6,*) "[1]           RK TIME  = ", TRK
 write(6,*) "[2]           BC TIME  = ", TBC
 write(6,*) "[3] SPATIAL (RHS) TIME = ", TRHS
 write(6,*) "[4]=[1]+[3]] Soltuion time except BC = ", TRK+ TRHS, "<<<<<<<<"
 write(6,'(A,3F12.4)') "     Check Measurement Error [1+2+3] vs TimeINT, DIFF= ", TtimeINT, TRK+TBC+ TRHS,  TtimeINT-( TRK+TBC+ TRHS)
 write(6,'(A)') "       SCHEME   NC     TI   NcellX           FpERR          FyzERR          FxzERR         FxyzERR        MaxERROR"
 write(6,'(A,I7,I5, I7, I9,5E16.6)') "FINREP", SptScheme , METHNC , TimeINTmethod, &
                          NcellX, FpERR, FyzERR, FxzERR, FxyzERR, max(  FpERR ,FyzERR,FxzERR ,FxyzERR )
 write(6,*) "///////////////////////////////^^^^^^^^^/////////////////////"
 write(6,'(A)') "****************************************************************************************************************"  
 

 write(file_unit,'(A)') &
            "****************************************************************************************************************"  
 write(file_unit,'(A,4I6)') "S1DT TimeINTmethod METHNC NcellX ", SptScheme ,TimeINTmethod, METHNC, NcellX
 write(file_unit,'(A,I6,4E18.5)') "NT  DX   DT   TF  CFL", n_steps, DX, DT,  t_final, CFL
 write(file_unit,*) "[1]           RK TIME  = ", TRK
 write(file_unit,*) "[2]           BC TIME  = ", TBC
 write(file_unit,*) "[3] SPATIAL (RHS) TIME = ", TRHS
 write(file_unit,*) "[4=[1]+[3]] Soltuion time except BC = ", TRK+ TRHS
 write(file_unit,'(A,3F12.4)') "     Check Measurement Error [1+2+3] vs TimeINT, DIFF= ", &
                                           TtimeINT, TRK+TBC+ TRHS,  TtimeINT-( TRK+TBC+ TRHS)
 write(file_unit,'(A)') " SCHEME   NC     TI   NcellX           FpERR          FyzERR          FxzERR         FxyzERR        MaxERROR"
 write(file_unit,'(A,I7,I5, I7, I9,5E16.6)') "FINREP", SptScheme , METHNC , TimeINTmethod, &
                          NcellX, FpERR, FyzERR, FxzERR, FxyzERR, max(  FpERR ,FyzERR,FxzERR ,FxyzERR )
 write(file_unit,'(A)') &
            "****************************************************************************************************************"  
 close(file_unit)
 

   
call COEF%free(); 
call  PHI%free();
call  Velocity%free();
call  TMP%free();
call  Grid01%free();




contains

   pure function uppercase(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: j, ia
      do j = 1, len(s)
         ia = iachar(s(j:j))
         if (ia >= iachar('a') .and. ia <= iachar('z')) then
            out(j:j) = achar(ia - 32)
         else
            out(j:j) = s(j:j)
         end if
      end do
   end function uppercase
   
 subroutine read_setup(fname, ng, ncell, cfl, u,K_diff, t_final, scase, method, scheme, METHNCs)

   implicit none

   character(len=*), intent(in)  :: fname
   integer,          intent(out) :: ng   
   integer,          intent(out) :: ncell
   real(wp),         intent(out) :: cfl
   real(wp),         intent(out) :: u, K_diff
   real(wp),         intent(out) :: t_final
   character(len=*), intent(out) :: scase
   character(len=*), intent(out) :: method
   character(len=*), intent(out) :: scheme
   character(len=*), intent(out) :: METHNCs   

   integer :: iu

   open(newunit=iu, file=fname, status='old', action='read')

   read(iu,*) ng
   read(iu,*) ncell
   read(iu,*) cfl
   read(iu,*) u, K_diff
   read(iu,*) t_final
   read(iu,'(A)') scase
   read(iu,'(A)') method
   read(iu,'(A)') scheme
   read(iu,'(A)') METHNCs   

   close(iu)

   scase  = adjustl(trim(scase))
   method = adjustl(trim(method))
   scheme = adjustl(trim(scheme))
   METHNCs = adjustl(trim(METHNCs))   
   
   call to_lower(scase) 
   call to_lower(method) 
   call to_lower(scheme)       
   call to_lower(METHNCs)     

end subroutine read_setup 

	subroutine to_lower(str)
	  character(len=*), intent(inout) :: str  ! intent(inout) allows modification
	  integer :: i, code

	  do i = 1, len(str)
	     code = iachar(str(i:i))
	     ! If the character is uppercase ('A' through 'Z')
	     if (code >= iachar('A') .and. code <= iachar('Z')) then
	        str(i:i) = achar(code + 32) ! Shift it down to lowercase
	     end if
	  end do
   end subroutine to_lower 


end program MMM_1Dconvect
