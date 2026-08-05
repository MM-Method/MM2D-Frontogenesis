module MMM_module_spatial
   use module_SchemeSetup, only:wp 
   use VTK_MultiMoment_Writer
   use module_physics
   implicit none
   private

   public :: get_erf_cellvalue, get_gauss_PV, GET_GAUSS_CV, get_fronto_initial_phi, get_fronto_initial_Vel, get_fronto_final_phi
   public :: save_paraview_cell_2d, set_MMMscheme, save_paraview_point_binary_2d_w_GC
   public :: MeshDomain, scalar_spatial_data,  velocity_spatial_data 
   public ::  save_paraview_point_binary_2d, Copy_SCALAR_AtoB
   public :: Coefficient_Differentiation ,  set_fronto_timedep_BC_simp
   public :: Compute_FirstDerivativeX_FD, Compute_FirstDerivativeY_FD
   public :: Compute_FirstDerivativeX_YZ,  Compute_FirstDerivativeY_XZ , ScaleField_scalar, update_scalar_ghosts_periodic


! ADD THESE TWO LINES:
   public :: Compute_SecondDerivativeX_FD, Compute_SecondDerivativeY_FD
   public :: Compute_SecondDerivativeX_YZ, Compute_SecondDerivativeY_XZ

    real(wp), parameter :: ep = 1.0e-14_wp
    real(wp) :: current_bc_time    
    
    ! Cache the constant factor 1 / ((2/9) * sqrt(3)) to save CPU cycles
 real(wp), parameter :: CONST_FACTOR = 1.0_wp / ((2.0_wp / 9.0_wp) * 1.732050807568877_wp)
   
   
!-------------------------------------------------------------------- type :: MeshDomain
   type :: MeshDomain
       integer iMeshDomain   ! Numerical ID of the domain
       integer NPx,NPy !Number of Points in x- and y-direction, respectively.
       integer NcellX, NcellY !Number of Cells in x- and y-direction, respectivel
       integer NPxtot,NPytot !Number of Points in x- and y-direction, respectively.
       integer NcellXtot, NcellYtot !Number of Cells in x- and y-direction, respectivel
       
       
       
       real (wp) ::    DX, DY !Grid Spacing in x- and y-direction. Currently, uniform  grid only
       integer NG                       ! Number of ghost cells
       integer IcSTART, IcEND ! first cell in X, last cell in X
       integer JcSTART, JcEND ! first cell in Y, last cell in Y
       
       integer IfSTART, IfEND ! first face (point) in X, last face (point) in X
       integer JfSTART, JfEND ! first face (point) in Y, last face (point) in Y       

       real(wp) :: X0, XL, Y0, YL
       real(wp) , dimension(:), allocatable :: Xp, Yp !location of nodal points (Xp,Yp) 
       real(wp) , dimension(:), allocatable :: Xc, Yc !location of cell-center (for visualization purpose)
 
    contains
        procedure :: init => allocate_Mesh
        procedure :: free => deallocate_Mesh        
   end type MeshDomain
 
 !--------------------------------------------------------------------- type :: scalar_spatial_data  
    type :: scalar_spatial_data
    
       type(MeshDomain), pointer :: Mesh => null()
       integer MeshID    ! MeshDomain this scalar data belong to  
       integer ProcID          ! ID of processor solving this grid
!       logical isMM            ! Is this a multi-moment data  
       integer SptScheme       !Spatialscheme
       
           
        real(wp), allocatable :: xyz(:,:) ! Cell values
        real(wp), allocatable :: xz(:,:)  ! XZ face values
        real(wp), allocatable :: yz(:,:)  ! YZ face values
        real(wp), allocatable :: p(:,:)   ! Point values
    contains
        procedure :: init => allocate_scalar
        procedure :: free => deallocate_scalar  
        procedure :: setZero => setZero_scalar          
        procedure :: ScaleField => ScaleField_scalar            
        procedure :: update_ghosts_periodic => update_scalar_ghosts_periodic    
    end type scalar_spatial_data
 !--------------------------------------------------------------------- type :: velocity_spatial_data  
    type :: velocity_spatial_data
    
       type(MeshDomain), pointer :: Mesh => null()
       integer MeshID    ! MeshDomain this scalar data belong to  
       integer ProcID          ! ID of processor solving this grid 
       logical isMM            ! Is this a multi-moment data  
       integer SptScheme       ! Spatialscheme 
       

       type(scalar_spatial_data) :: U,V

    contains
        procedure :: init => allocate_velocity
        procedure :: free => deallocate_velocity    
    end type velocity_spatial_data   

  !--------------------------------------------------------- type ::Coefficient_Differentiation 
	type Coefficient_Differentiation  
 
	    
	    real(wp), allocatable, dimension(:) :: COEF_FDX_B1_L, COEF_FDX_B1_R, COEF_FDX_B2_L, COEF_FDX_B2_R
	    real(wp), allocatable, dimension(:) :: COEF_SDX_B1_L, COEF_SDX_B1_R, COEF_SDX_B2_L, COEF_SDX_B2_R
	    real(wp), allocatable, dimension(:) :: COEF_FDX_B1all,COEF_FDX_B2all    
	    real(wp), allocatable, dimension(:) :: COEF_SDX_B1all,COEF_SDX_B2all    
	    real(wp), allocatable, dimension(:) :: COEF_FDY_B1_L, COEF_FDY_B1_R, COEF_FDY_B2_L, COEF_FDY_B2_R
	    real(wp), allocatable, dimension(:) :: COEF_SDY_B1_L, COEF_SDY_B1_R, COEF_SDY_B2_L, COEF_SDY_B2_R
	    real(wp), allocatable, dimension(:) :: COEF_FDY_B1all,COEF_FDY_B2all    
	    real(wp), allocatable, dimension(:) :: COEF_SDY_B1all,COEF_SDY_B2all      
	    
	    real(wp) COEF_SDX_B2_M, COEF_FDY_B1_M
	    real(wp) COEF_SDY_B2_M, COEF_FDY_B2_M

	    integer :: FDB1BW, FDB2BW, SDB1BW, SDB2BW !approximation bandwidth of B1&B2 of first derivative and second-derivative
	 
	    contains
		procedure :: init => allocate_coefdiff
		procedure :: free => deallocate_coefdiff  
		    
	end type Coefficient_Differentiation
  !--------------------------------------------------------------------- type :: Coefficient_Diffe
   
!--------------------------------------------------   
contains
!-------------------------------------------------



	SUBROUTINE allocate_coefdiff( G, Mesh )
		use module_SchemeSetup 
		class(Coefficient_Differentiation), intent(inout) :: G
		type(MeshDomain)  :: Mesh
		integer nx,ny,i,j
		real (wp) :: DX, DY, DX2, DY2
		
		nx = Mesh%NPxtot
		ny = Mesh%NPytot
		DX = Mesh%DX; DX2 = (Mesh%DX)**2.0;
		DY = Mesh%DY; DY2 = (Mesh%DY)**2.0;
	 
!		print *, "allocate coeff FDB1BW, FDB2B",    FDB1BW, FDB2BW    
		allocate( &
		          G%COEF_FDX_B2_L( FDB2BW ), G%COEF_FDX_B2_R( FDB2BW ), &
		          G%COEF_FDY_B2_L( FDB2BW ), G%COEF_FDY_B2_R( FDB2BW ),  &
		          G%COEF_SDX_B2_L( SDB2BW ), G%COEF_SDX_B2_R( SDB2BW ), &
		          G%COEF_SDY_B2_L( SDB2BW ), G%COEF_SDY_B2_R( SDB2BW ) &
		          )
		         
		 if(isMM) allocate( &
		          G%COEF_FDX_B1_L( FDB1BW ), G%COEF_FDX_B1_R( FDB1BW ), &
		          G%COEF_FDY_B1_L( FDB1BW ), G%COEF_FDY_B1_R( FDB1BW ), &
		          G%COEF_SDX_B1_L( SDB1BW ), G%COEF_SDX_B1_R( SDB1BW ), &
		          G%COEF_SDY_B1_L( SDB1BW ), G%COEF_SDY_B1_R( SDB1BW ) &
		          )  
		          
		CALL gen_FDIFF_coef_uni( G%COEF_FDX_B1_L, G%COEF_FDX_B1_R,  &
		                         G%COEF_FDX_B2_L, G%COEF_FDX_B2_R, DX )
		CALL gen_FDIFF_coef_uni( G%COEF_FDY_B1_L, G%COEF_FDY_B1_R, &
		                         G%COEF_FDY_B2_L, G%COEF_FDY_B2_R, DY ) 
		                         
		
	     contains
	      
     SUBROUTINE gen_FDIFF_coef_uni( COEF_FDX_B1_L, COEF_FDX_B1_R, COEF_FDX_B2_L, COEF_FDX_B2_R, DX  )
         use module_SchemeSetup	     
         IMPLICIT NONE
                      
		 real (wp) COEF_FDX_B1_L( : ), COEF_FDX_B1_R( : )
		 real (wp) COEF_FDX_B2_L( FDB2BW ), COEF_FDX_B2_R( FDB2BW )
		 real (wp) DX

          if(isMM) THEN 
                COEF_FDX_B1_L(:) = 0.0 ;  COEF_FDX_B1_R(:) = 0.0 
          ENDIF                
                COEF_FDX_B2_L(:) = 0.0 ;  COEF_FDX_B2_R(:) = 0.0 
		select case ( SptScheme ) 
		!-------------------------------MM Method
			case(3002) 
				 COEF_FDX_B1_L(1) = -1.0_wp/DX 
				 COEF_FDX_B1_R(1) =  1.0_wp/DX 
			case(3004) 
				 COEF_FDX_B2_L(1) =  0.5_wp/DX         	
				 COEF_FDX_B1_L(1) = -2.0_wp/DX 
				 COEF_FDX_B1_R(1) =  2.0_wp/DX  
				 COEF_FDX_B2_R(1) = -0.5_wp/DX 
!				 print *, "COEF MMEX4 : for DX", DX
!				 write(6,'(4E23.16)' ),    COEF_FDX_B2_L(1), COEF_FDX_B1_L(1),COEF_FDX_B1_R(1),   COEF_FDX_B2_R(1)
				 
			case(3006)
		 		 COEF_FDX_B1_L(2) = -1.0_wp/36.0_wp /DX       		  
				 COEF_FDX_B2_L(1) =  2.0_wp/ 3.0_wp /DX 
				 COEF_FDX_B1_L(1) = -9.0_wp/ 4.0_wp /DX 
				 COEF_FDX_B1_R(1) =  9.0_wp/ 4.0_wp /DX  
				 COEF_FDX_B2_R(1) = -2.0_wp/ 3.0_wp /DX
	     		 COEF_FDX_B1_R(2) =  1.0_wp/36.0_wp /DX  
	     		case(3008)	
	     			 COEF_FDX_B1_L(1) = -5.0_wp/ 2.0_wp /DX 
	     			 COEF_FDX_B1_L(2) = -7.0_wp/54.0_wp /DX 
	     			 COEF_FDX_B1_R(1:2) = -COEF_FDX_B1_L(1:2)
	     			 
	     			 COEF_FDX_B2_L(1) = 8.0_wp /  9.0_wp /DX
	     			 COEF_FDX_B2_L(2) = 1.0_wp / 36.0_wp /DX
	     			 COEF_FDX_B2_R   = -COEF_FDX_B2_L
	     		case(3010)      			 
	     			 COEF_FDX_B1_L(1) = -47.0_wp/ 18.0_wp /DX 
	     			 COEF_FDX_B1_L(2) =  -7.0_wp/ 36.0_wp /DX 
	     			 COEF_FDX_B1_L(3) = -1.0_wp/ 900.0_wp /DX 		 
	     			 COEF_FDX_B1_R(1:3) = -COEF_FDX_B1_L(1:3)
	     			 
	     			 COEF_FDX_B2_L(1) = 1.0_wp /DX
	     			 COEF_FDX_B2_L(2) = 1.0_wp / 20.0_wp /DX
	     			 COEF_FDX_B2_R   = -COEF_FDX_B2_L 
	      		case(3012)      			 
	     			 COEF_FDX_B1_L(1) =  -49.0_wp/    18.0_wp /DX 
	     			 COEF_FDX_B1_L(2) =  -41.0_wp/   144.0_wp /DX 
	     			 COEF_FDX_B1_L(3) = -157.0_wp/ 18000.0_wp /DX 		 
	     			 COEF_FDX_B1_R(1:3) = -COEF_FDX_B1_L(1:3)
	     			 
	     			 COEF_FDX_B2_L(1) = 9.0_wp /  8.0_wp/DX
	     			 COEF_FDX_B2_L(2) = 9.0_wp /100.0_wp /DX
	     			 COEF_FDX_B2_L(3) = 1.0_wp /600.0_wp /DX
	     			 COEF_FDX_B2_R(:)   = -COEF_FDX_B2_L(:)  
		!-------------------------------FD Method
                       case(7002) 
						 COEF_FDX_B2_L(1) = -0.5_wp/DX 
						 COEF_FDX_B2_R(1) =  0.5_wp/DX
                       case(7004) 
						 COEF_FDX_B2_L(2) =  1.0_wp/ (12.0_wp*DX) 
						 COEF_FDX_B2_L(1) = -8.0_wp/ (12.0_wp*DX)
						 COEF_FDX_B2_R(:) = -COEF_FDX_B2_L(:)	 
                         case(7006)
                           COEF_FDX_B2_L(3) = -1.0_wp/60.0_wp / DX 
                           COEF_FDX_B2_L(2) =  3.0_wp/20.0_wp / DX
                           COEF_FDX_B2_L(1) = -3.0_wp/ 4.0_wp / DX 
                           COEF_FDX_B2_R(:) = -COEF_FDX_B2_L(:) 
			 case(7008) ! 8th-order Finite Difference Scheme
			    COEF_FDX_B2_L(1) = -4.0_wp / 5.0_wp / DX
			    COEF_FDX_B2_L(2) =  1.0_wp / 5.0_wp / DX
			    COEF_FDX_B2_L(3) = -4.0_wp / 105.0_wp / DX
			    COEF_FDX_B2_L(4) =  1.0_wp / 280.0_wp / DX
			    COEF_FDX_B2_R(1:4) = -COEF_FDX_B2_L(1:4)
			case(7010) ! 10th-order Finite Difference Scheme
			    COEF_FDX_B2_L(1) = -5.0_wp / 6.0_wp / DX
			    COEF_FDX_B2_L(2) =  5.0_wp / 21.0_wp / DX
			    COEF_FDX_B2_L(3) = -5.0_wp / 84.0_wp / DX
			    COEF_FDX_B2_L(4) =  5.0_wp / 504.0_wp / DX
			    COEF_FDX_B2_L(5) = -1.0_wp / 1260.0_wp / DX
			    COEF_FDX_B2_R(1:5) = -COEF_FDX_B2_L(1:5)
			case(7012) ! 12th-order Finite Difference Scheme
			    COEF_FDX_B2_L(1) = -6.0_wp / 7.0_wp / DX
			    COEF_FDX_B2_L(2) =  15.0_wp / 56.0_wp / DX
			    COEF_FDX_B2_L(3) = -5.0_wp / 63.0_wp / DX
			    COEF_FDX_B2_L(4) =  1.0_wp / 56.0_wp / DX
			    COEF_FDX_B2_L(5) = -1.0_wp / 385.0_wp / DX
			    COEF_FDX_B2_L(6) =  1.0_wp / 5544.0_wp / DX
			    COEF_FDX_B2_R(1:6) = -COEF_FDX_B2_L(1:6)	 		 
			     			    			 
		case default  
		        error stop 'Unknown Spatial Scheme: gen_FDIFF_coef_uni'
		end select             
	    END SUBROUTINE gen_FDIFF_coef_uni 
 
	END SUBROUTINE allocate_coefdiff

 
SUBROUTINE gen_SDIFF_coef_uni( COEF_SDX_B1_L, COEF_SDX_B1_R, &
                               COEF_SDX_B2_L, COEF_SDX_B2_R, &
                               COEF_SDX_B2_M, COEF_SDX_B1_M, DX )
   use module_SchemeSetup
   IMPLICIT NONE

   real (wp) :: COEF_SDX_B1_L( : ), COEF_SDX_B1_R( : )
   real (wp) :: COEF_SDX_B2_L( SDB2BW ), COEF_SDX_B2_R( SDB2BW )
   real (wp) :: COEF_SDX_B2_M, COEF_SDX_B1_M
   real (wp) :: DX
   real (wp) :: iDX2

   iDX2 = 1.0_wp / (DX * DX)

   ! Zero out coefficients prior to assignment
   if (isMM) THEN 
      COEF_SDX_B1_L(:) = 0.0_wp; COEF_SDX_B1_R(:) = 0.0_wp
      COEF_SDX_B1_M    = 0.0_wp
   ENDIF               
   COEF_SDX_B2_L(:) = 0.0_wp; COEF_SDX_B2_R(:) = 0.0_wp
   COEF_SDX_B2_M    = 0.0_wp

   SELECT CASE ( SptScheme ) 

   !---------------------------------------------------------------------------
   ! Multi-Moment (MM) Second-Derivative Schemes (From Table 3)
   !---------------------------------------------------------------------------
   case(3002) ! 2nd-order MM
      COEF_SDX_B2_M    = -6.0_wp * iDX2
      COEF_SDX_B2_L(1) =  3.0_wp * iDX2
      COEF_SDX_B2_R(1) =  3.0_wp * iDX2

   case(3004) ! 4th-order MM
      COEF_SDX_B2_M    = -12.0_wp * iDX2
      
      ! Primary cell coefficients (b1)
      COEF_SDX_B2_L(1) = (15.0_wp / 2.0_wp) * iDX2
      COEF_SDX_B2_R(1) = COEF_SDX_B2_L(1)
      
      ! Face derivative coefficients (c1)
      COEF_SDX_B1_L(1) = (-3.0_wp / 2.0_wp) * iDX2
      COEF_SDX_B1_R(1) = COEF_SDX_B1_L(1)

   case(3006) ! 6th-order MM
      COEF_SDX_B2_M    = (-27.0_wp / 2.0_wp) * iDX2
      
      ! Primary cell coefficients (b1, b2)
      COEF_SDX_B2_L(1) = (209.0_wp / 24.0_wp) * iDX2
      COEF_SDX_B2_L(2) = (  1.0_wp / 24.0_wp) * iDX2
      COEF_SDX_B2_R(1:2) = COEF_SDX_B2_L(1:2)
      
      ! Face derivative coefficients (c1)
      COEF_SDX_B1_L(1) = -2.0_wp * iDX2
      COEF_SDX_B1_R(1) = COEF_SDX_B1_L(1)

   case(3008) ! 8th-order MM
      COEF_SDX_B2_M    = -15.0_wp * iDX2
      
      ! Primary cell coefficients (b1, b2)
      COEF_SDX_B2_L(1) = (1439.0_wp / 144.0_wp) * iDX2
      COEF_SDX_B2_L(2) = (  31.0_wp / 144.0_wp) * iDX2
      COEF_SDX_B2_R(1:2) = COEF_SDX_B2_L(1:2)
      
      ! Face derivative coefficients (c1, c2)
      COEF_SDX_B1_L(1) = (-8.0_wp /  3.0_wp) * iDX2
      COEF_SDX_B1_L(2) = (-1.0_wp / 24.0_wp) * iDX2
      COEF_SDX_B1_R(1:2) = COEF_SDX_B1_L(1:2)

   case(3010) ! 10th-order MM
      COEF_SDX_B2_M    = (-47.0_wp / 3.0_wp) * iDX2
      
      ! Primary cell coefficients (b1, b2, b3)
      COEF_SDX_B2_L(1) = (38083.0_wp / 3600.0_wp) * iDX2
      COEF_SDX_B2_L(2) = ( 1183.0_wp / 3600.0_wp) * iDX2
      COEF_SDX_B2_L(3) = (    1.0_wp /  900.0_wp) * iDX2
      COEF_SDX_B2_R(1:3) = COEF_SDX_B2_L(1:3)
      
      ! Face derivative coefficients (c1, c2)
      COEF_SDX_B1_L(1) = -3.0_wp * iDX2
      COEF_SDX_B1_L(2) = (-3.0_wp / 40.0_wp) * iDX2
      COEF_SDX_B1_R(1:2) = COEF_SDX_B1_L(1:2)

   case(3012) ! 12th-order MM
      COEF_SDX_B2_M    = (-49.0_wp / 3.0_wp) * iDX2
      
      ! Primary cell coefficients (b1, b2, b3)
      COEF_SDX_B2_L(1) = (201209.0_wp / 18000.0_wp) * iDX2
      COEF_SDX_B2_L(2) = (  4417.0_wp /  9000.0_wp) * iDX2
      COEF_SDX_B2_L(3) = (   167.0_wp / 18000.0_wp) * iDX2
      COEF_SDX_B2_R(1:3) = COEF_SDX_B2_L(1:3)
      
      ! Face derivative coefficients (c1, c2, c3)
      COEF_SDX_B1_L(1) = ( -27.0_wp /   8.0_wp) * iDX2
      COEF_SDX_B1_L(2) = ( -27.0_wp / 200.0_wp) * iDX2
      COEF_SDX_B1_L(3) = (  -1.0_wp / 600.0_wp) * iDX2
      COEF_SDX_B1_R(1:3) = COEF_SDX_B1_L(1:3)

   !---------------------------------------------------------------------------
   ! Standard Finite Difference (FD) Second-Derivative Schemes
   !---------------------------------------------------------------------------
   case(7002) ! 2nd-order FD
      COEF_SDX_B2_M    = -2.0_wp * iDX2
      COEF_SDX_B2_L(1) =  1.0_wp * iDX2
      COEF_SDX_B2_R(1) =  1.0_wp * iDX2

   case(7004) ! 4th-order FD
      COEF_SDX_B2_M    = (-30.0_wp / 12.0_wp) * iDX2
      COEF_SDX_B2_L(1) = ( 16.0_wp / 12.0_wp) * iDX2
      COEF_SDX_B2_L(2) = ( -1.0_wp / 12.0_wp) * iDX2
      COEF_SDX_B2_R(1:2) = COEF_SDX_B2_L(1:2)

   case(7006) ! 6th-order FD
      COEF_SDX_B2_M    = (-490.0_wp / 180.0_wp) * iDX2
      COEF_SDX_B2_L(1) = ( 270.0_wp / 180.0_wp) * iDX2
      COEF_SDX_B2_L(2) = ( -27.0_wp / 180.0_wp) * iDX2
      COEF_SDX_B2_L(3) = (   2.0_wp / 180.0_wp) * iDX2
      COEF_SDX_B2_R(1:3) = COEF_SDX_B2_L(1:3)

   case(7008) ! 8th-order FD
      COEF_SDX_B2_M    = (-14350.0_wp / 5040.0_wp) * iDX2
      COEF_SDX_B2_L(1) = (  8064.0_wp / 5040.0_wp) * iDX2
      COEF_SDX_B2_L(2) = ( -1008.0_wp / 5040.0_wp) * iDX2
      COEF_SDX_B2_L(3) = (   128.0_wp / 5040.0_wp) * iDX2
      COEF_SDX_B2_L(4) = (    -9.0_wp / 5040.0_wp) * iDX2
      COEF_SDX_B2_R(1:4) = COEF_SDX_B2_L(1:4)

   case(7010) ! 10th-order FD
      COEF_SDX_B2_M    = (-73766.0_wp / 25200.0_wp) * iDX2
      COEF_SDX_B2_L(1) = ( 42000.0_wp / 25200.0_wp) * iDX2
      COEF_SDX_B2_L(2) = ( -6000.0_wp / 25200.0_wp) * iDX2
      COEF_SDX_B2_L(3) = (  1000.0_wp / 25200.0_wp) * iDX2
      COEF_SDX_B2_L(4) = (  -125.0_wp / 25200.0_wp) * iDX2
      COEF_SDX_B2_L(5) = (     8.0_wp / 25200.0_wp) * iDX2
      COEF_SDX_B2_R(1:5) = COEF_SDX_B2_L(1:5)

   case(7012) ! 12th-order FD
      COEF_SDX_B2_M    = (-12012.0_wp / 3969.0_wp) * iDX2
      COEF_SDX_B2_L(1) = (  6912.0_wp / 3969.0_wp) * iDX2
      COEF_SDX_B2_L(2) = ( -1080.0_wp / 3969.0_wp) * iDX2
      COEF_SDX_B2_L(3) = (   224.0_wp / 3969.0_wp) * iDX2
      COEF_SDX_B2_L(4) = (   -36.0_wp / 3969.0_wp) * iDX2
      COEF_SDX_B2_L(5) = (    24.0_wp / 6615.0_wp) * iDX2
      COEF_SDX_B2_L(6) = (    -1.0_wp / 7938.0_wp) * iDX2
      COEF_SDX_B2_R(1:6) = COEF_SDX_B2_L(1:6)

   case default  
      error stop 'Unknown Spatial Scheme: gen_SDIFF_coef_uni'
   end select             

END SUBROUTINE gen_SDIFF_coef_uni
 


	SUBROUTINE deallocate_coefdiff( G )
                use module_SchemeSetup	 , only: isMM
		class(Coefficient_Differentiation), intent(inout) :: G
 

		IF( allocated( G%COEF_FDX_B2_L ) ) THEN
			deallocate( &
				  G%COEF_FDX_B2_L, G%COEF_FDX_B2_R, &
				  G%COEF_FDY_B2_L, G%COEF_FDY_B2_R,  &
				  G%COEF_SDX_B2_L, G%COEF_SDX_B2_R, &
				  G%COEF_SDY_B2_L, G%COEF_SDY_B2_R &
				  )
				 
			 if(isMM) deallocate( &
				  G%COEF_FDX_B1_L, G%COEF_FDX_B1_R, &
				  G%COEF_FDY_B1_L, G%COEF_FDY_B1_R, &
				  G%COEF_SDX_B1_L, G%COEF_SDX_B1_R, &
				  G%COEF_SDY_B1_L, G%COEF_SDY_B1_R &
				  )  
		 ENDIF


	END SUBROUTINE deallocate_coefdiff
	
	
 

!!!!    Velocity contain metadata indentical to the scalar. Thus data are duplicate
 
subroutine set_MMMscheme(   method, METHNCs, SptScheme, isMM , METHNC )
 
	character(len=16), intent(in) ::   method , METHNCs
	integer, intent(out) ::   SptScheme, METHNC
	logical, intent(out) ::   isMM
	


	
      isMM =.FALSE. !default	
      select case ( method )  
      	case ( 'mmex4' ) !Explicit 4th-order 
      		SptScheme = 3004
      		isMM = .TRUE. ;  
      	case ( 'mmex6' ) !Explicit 4th-order 
      		SptScheme = 3006
      		isMM = .TRUE. ;    		 
      	case ( 'mmex8' ) !standard RK5
      		SptScheme = 3008 
      		isMM = .TRUE. ;    
      	case ( 'mmex10' ) !standard RK5
      		SptScheme = 3010 
      		isMM = .TRUE. ;        		
      	case ( 'mmex12' ) !standard RK5
      		SptScheme = 3012       	
      		isMM = .TRUE. ;      			
      	case ( 'mmex14' ) !standard RK5
      		SptScheme = 3014      		
      	case ( 'mmex24' ) !standard RK5
      		SptScheme = 3024   
!------------------- FDM EXPLICIT ---------------
      	case ( 'fdex2' ) !Explicit 4th-order 
      		SptScheme = 7002
      	case ( 'fdex4' ) !Explicit 4th-order 
      		SptScheme = 7004
      	case ( 'fdex6' ) !Explicit 4th-order 
      		SptScheme = 7006 
      	case ( 'fdex8' ) !standard RK5
      		SptScheme = 7008 
      	case ( 'fdex10' ) !standard RK5
      		SptScheme = 7010 
!------------------- FDM Implicit ---------------
      	case ( 'fdcom4' ) !Compact 4th-order 
      		SptScheme = 8004
      	case ( 'fdcom6' ) !Compact 4th-order 
      		SptScheme = 8006 
      	case ( 'fdcom8' ) !Compact RK5
      		SptScheme = 8008       		   		
        case default
                print *, "INPUT Spatial Scheme METHOD IS:", method
                error stop 'Unknown Spatial Scheme: set_MMMscheme'
        end select     
        
        if( SptScheme .LT. 7000) isMM = .TRUE.
        
        
            select case ( METHNCs )  
            case ( 'nc0' ) !Compact 4th-order 
            METHNC = 0
            case ( 'nc4' ) !Compact 4th-order 
            METHNC = 4
            case ( 'nc6' ) !Compact 4th-order 
            METHNC = 6            
			case default
			        print *, "Only nc, nc4, and nc0 are accepted. This was given:", METHNCs
			        error stop 'Unknown Nonlinear correction: set_MMMscheme'
			end select                     	
                  		      		      			
end subroutine set_MMMscheme



	SUBROUTINE Copy_SCALAR_AtoB( A, B)
  	    use module_SchemeSetup, only: isMM
		type(scalar_spatial_data), intent(in) :: A
		type(scalar_spatial_data), intent(OUT) :: B  
		integer NY,NX,i
		
		NY = A%Mesh%NPytot !SIZE(  A%p, 1)
		NX = A%Mesh%NPxtot !SIZE(  A%p, 2)
		B%Mesh%NPxtot = NX
		B%Mesh%NPytot = NY
 
if(1==1) THEN		
		IF(  (NX ==  B%Mesh%NPxtot ) .AND. (NY ==  B%Mesh%NPytot) ) THEN

		   B%MeshID =A%MeshID    ! MeshDomain this scalar data belong to    
		   B%SptScheme =A%SptScheme       !Spatialscheme	
		   IF( isMM ) THEN 
				DO i = 1,NX - 1			
					 B%xz(:,i)=A%xz(:,i)  ! XZ face values			 		
				ENDDO
				
				DO i = 1,NX - 1			
					 B%xyz(:,i)=A%xyz(:,i) ! Cell values	 		
				ENDDO			

				DO i = 1,NX		              	
					 B%yz(:,i)=A%yz(:,i)  ! YZ face values			
				ENDDO	
			ENDIF	
			DO i = 1,NX		              
				 B%p (:,i)=A%p(:,i)   ! Point values			
			ENDDO						
			
		ELSE	
			print *, "ERROR: SIZE MISMATCHED"	
		ENDIF
ELSE

   B%p = A%p
   IF (isMM) THEN
      B%xyz = A%xyz
      B%xz  = A%xz
      B%yz  = A%yz
   END IF

ENDIF		

	END SUBROUTINE Copy_SCALAR_AtoB

 
   subroutine allocate_scalar(this, source_mesh)
        use module_SchemeSetup
        class(scalar_spatial_data), intent(inout) :: this
        type(MeshDomain),  intent(in), target :: source_mesh
        integer nx,ny
        
             
        this%Mesh => source_mesh  
        nx = source_mesh%NPx +  2*source_mesh%ng;
        ny = source_mesh%NPy +  2*source_mesh%ng;

	allocate(  this%p  ( ny    , nx    )  )    
	this%p = 0.0_wp;       
	
        IF( isMM ) THEN
        	allocate(this%xyz( ny - 1, nx - 1), &
		         this%xz ( ny    , nx - 1    ), &
		         this%yz ( ny - 1, nx    )   )       
		this%xz =  0.0_wp; 
		this%yz =  0.0_wp; 
		this%xyz  =  0.0_wp;		          
        
        ENDIF

    end subroutine allocate_scalar

    subroutine allocate_velocity(this, source_mesh) 
        class(velocity_spatial_data), intent(inout) :: this
        type(MeshDomain),  intent(in), target :: source_mesh
        integer  :: nx, ny

        this%Mesh => source_mesh          
        nx = source_mesh%NPx ;
        ny = source_mesh%NPy        
        
        call allocate_scalar(this%U, source_mesh )
        call allocate_scalar(this%V, source_mesh )
    end subroutine allocate_velocity    

!--------------- allocate_Mesh allocate_Mesh
subroutine allocate_Mesh(this, iMesh, ng, nx, ny, X0, XL, Y0, YL)
        class(MeshDomain), intent(inout) :: this
        integer, intent(in) :: iMesh, nx, ny, ng
        real (wp), intent(in) ::   X0, XL, Y0, YL
        
        real (wp) :: DX,DY
        integer i, NPXtot, NPYtot
        
         !print *, iMesh, ng,nx,ny
        
! store mesh information in the mesdata
         this%iMeshDomain = iMesh ;
         this%Ng = ng ; 
         this%NPx = nx;       
         this%NPy = ny;
         this%NcellX = nx - 1;
         this%NcellY = ny - 1; 
         this%NPxtot = nx +2*ng ;
         this%NPytot = ny +2*ng ;          
         this%NcellXtot = this%NcellX +2*ng 
         this%NcellYtot = this%NcellY +2*ng    
                    
         this%X0 = X0 ;       this%Y0 = Y0;
         this%XL = XL ;       this%YL = YL  ;

         this%IcStart = Ng+1           ;this%IfStart  = Ng+1 ;
         this%JcStart = Ng+1           ;this%JfStart  = Ng+1   ;
         this%IcEnd   = Ng+this%NcellX ;this%IfEnd   = Ng+this%NcellX+1 ;
         this%JcEnd   = Ng+this%NcellY ;this%JfEnd   = Ng+this%NcellY+1 ;         
         
         NPXtot =     this%NPx +2*ng ; NPYtot =     this%NPY +2*ng
         !print *, (XL-X0), this%NcellX, real(   this%NcellX ,kind=wp ) 
         DX = (XL-X0)/real(   this%NcellX ,kind=wp )
         DY = (YL-Y0)/real(   this%NcellY ,kind=wp )
         
         this%DX = DX
         this%DY = DY
         
        allocate( this%Xp( NPXtot  ), &  
                  this%Yp( NPYtot  ), & 
                  this%Xc( NPXtot - 1 ),  &
                  this%Yc( NPYtot - 1 )  )
                  
        this%Xp(:) = 0.0_wp
        this%Yp(:) = 0.0_wp
        this%Xc(:) = 0.0_wp
        this%Yc(:) = 0.0_wp
        
        this%Xp(1) = X0 -     real(   ng,kind=wp ) *DX
        DO i=2,NPXtot
        	this%Xp(i)   =   this%Xp(i-1) + DX
        	this%Xc(i-1) =  0.5*(this%Xp(i-1)  + this%Xp(i))
        ENDDO

        this%Yp(1) = Y0 -     real(   ng,kind=wp ) *DY    
        DO i=2,NPYtot
        	this%Yp(i) =   this%Yp(i-1) + DY
        	this%Yc(i-1) =  0.5*(this%Yp(i-1)  + this%Yp(i))        	
        ENDDO
        
        
      
            
end  subroutine allocate_Mesh


subroutine deallocate_Mesh(this)
        class(MeshDomain), intent(inout) :: this
        if (allocated(this%Xp)) deallocate(this%Xp)
        if (allocated(this%Yp)) deallocate(this%Yp)
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Yc)) deallocate(this%Yc)
end subroutine deallocate_Mesh

!--------------- end deallocate_Mesh deallocate_Mesh

 




    subroutine deallocate_scalar(this)
        use module_SchemeSetup
        class(scalar_spatial_data), intent(inout) :: this
        IF (allocated(this%p)) then 
        	deallocate( this%p )
       	        IF( isMM ) THEN
	        	deallocate(this%xyz, this%xz, this%yz )       	        
       	        ENDIF
       	ENDIF
        	
        
    end subroutine deallocate_scalar


    subroutine setZero_scalar(this)
    use module_SchemeSetup, only:wp, isMM 
    class(scalar_spatial_data), intent(inout) :: this 
    
      IF( isMM ) THEN  
		  this%xyz = 0.0_wp
		  this%xz = 0.0_wp
		  this%yz = 0.0_wp
      ENDIF
      this%p  = 0.0_wp      

    end subroutine setZero_scalar    
    
    subroutine ScaleField_scalar(this, fak )
    use module_SchemeSetup, only:wp, isMM 
    class(scalar_spatial_data), intent(inout) :: this 
    REAL(WP) fak
    
      IF( isMM ) THEN  
		  this%xyz(:,:) = this%xyz(:,:) * fak
		  this%xz(:,:) = this%xz(:,:) * fak
		  this%yz(:,:) = this%yz(:,:) * fak
      ENDIF
      this%p(:,:)  = this%p(:,:) * fak      

    end subroutine ScaleField_scalar     
    
    
    
    



    subroutine deallocate_velocity(this)
        class(velocity_spatial_data), intent(inout) :: this
        call deallocate_scalar(this%U)
        call deallocate_scalar(this%V)
    end subroutine deallocate_velocity






   subroutine update_scalar_ghosts_periodic(F)
         class(scalar_spatial_data), intent(inout) :: F  
   

      integer :: nx, ny, ng
      
      nx = F%Mesh%NcellX 
      ny = F%Mesh%NcellY
      ng = F%Mesh%Ng
     
     !FV direction : START
      F%xyz(  : ,  1:ng    ) = F%xyz( : , (nx+1):(nx+ng)    )
       F%xz(  : ,  1:ng    ) =  F%xz( : , (nx+1):(nx+ng)    )
      F%xyz(       1:ng, : ) = F%xyz(     (ny+1):(ny+ng), : )
       F%xz(       1:ng, : ) =  F%xz(     (ny+1):(ny+ng), : )

     !FD direction: START 
      F%xyz(       1:ng ,  :  ) = F%xyz(     (ny+1):(ny+ng), : )
       F%yz(       1:ng ,  :  ) =  F%yz(     (ny+1):(ny+ng), : )
      
      ! FV direction: END
      F%xyz(  : , (nx+ng+1):(nx+2*ng)    ) = F%xyz( : ,(ng+1):(2*ng)    )
       F%xz(  : , (nx+ng+1):(nx+2*ng)    ) =  F%xz( : ,(ng+1):(2*ng)    )
      F%xyz(      (ny+ng+1):(ny+2*ng), : ) = F%xyz(    (ng+1):(2*ng), : )
       F%yz(      (ny+ng+1):(ny+2*ng), : ) =  F%yz(    (ng+1):(2*ng), : )
       
       ! FD direction: END
       F%yz(  : , (nx+ng+1):(nx+2*ng+1)    )  =  F%yz(  : , (ng+1):(2*ng+1)    )
       F%xz(      (ny+ng+1):(ny+2*ng+1), : )  =  F%xz(      (ng+1):(2*ng+1), : )  
       
      
      

      ! FDM
      F%P(  : ,  1:ng       ) = F%P( : , (nx+1):(nx+ng)    )
      F%P(       1:ng ,  :  ) = F%P(     (ny+1):(ny+ng), : )
      F%P(  : , (nx+ng+1):(nx+2*ng+1)    )  =  F%P(  : , (ng+1):(2*ng+1)  )
      F%P(      (ny+ng+1):(ny+2*ng+1), : )  =  F%P(      (ng+1):(2*ng+1), :   )
      
 
   end subroutine update_scalar_ghosts_periodic

  


 elemental function get_gauss_PV(A, sigma0, K, x, t) result(value)
       real(wp) :: value
       real(wp), intent(in) :: A,sigma0, K,t,  x  
              
  value = A*sigma0*exp(-x**2/(4.0_wp*K*t+sigma0**2))/sqrt(4.0_wp*K*t+sigma0**2)
 end function get_gauss_PV 

 elemental function get_gauss_CV(A, sigma0, K, xa,xb, t) result(value)
       real(wp) :: value,dx
       real(wp), intent(in) :: A,sigma0, K,t,  xa,xb  
   real(wp), parameter :: pi = acos(-1.0_wp)       
       dx=xb-xa
              
  value = -(0.5_wp)*A*sigma0*sqrt(Pi)*(erf(xa/sqrt(4.0_wp*K*t+sigma0*sigma0))-erf(xb/sqrt(4.0_wp*K*t+sigma0*sigma0)))/dx
  
 end function get_gauss_CV  



   subroutine get_fronto_initial_phi( PHI )  
   use module_SchemeSetup, only: isMM, NG
   type(scalar_spatial_data), intent(inout) :: PHI
   real(wp) y1,y2
   integer i,j   
   ! AT INITIAL CONDITON PHI(x,y) = -tanh( y/2) 
   
!        print *, "/X/: ", PHI%Mesh%Xp
!        print *, "/Y/: ", PHI%Mesh%Yp        


!!   Do i=1,PHI%Mesh%NcellXtot   
	Do j=1,PHI%Mesh%NPytot   
  		PHI%P(j,:)=  -tanh( 0.5_wp * PHI%Mesh%Yp(j) )
	ENDDO	
!!   ENDDO
	
   IF( isMM ) THEN
!!	Do i=1,PHI%Mesh%NcellXtot   
		Do j=1,PHI%Mesh%NPytot   
			PHI%xz(j,:)=    PHI%P(j, NG+1 ) !!-tanh( 0.5_wp * PHI%Mesh%Yp(j) )
		ENDDO	
	!!ENDDO	

!!	Do i=1,PHI%Mesh%NcellXtot   
		Do j=1,PHI%Mesh%NcellYtot   
			y1 = PHI%Mesh%Yp(j)
			y2 = PHI%Mesh%Yp(j+1)	        
			PHI%yz(j,:) =  &
				(log(abs( tanh((0.5_wp)*y1)-1.0_wp))+&
				 log(abs( tanh((0.5_wp)*y1)+1.0_wp))-&
				 log(abs( tanh((0.5_wp)*y2)-1.0_wp))-&
				 log(abs( tanh((0.5_wp)*y2)+1.0_wp)))/(y1-y2)
			PHI%xyz(j,:) =  PHI%yz(j,NG+1)
		ENDDO	
!!	ENDDO	
 
	 		
   ENDIF

   end subroutine get_fronto_initial_phi   
   
   
   
   
   !--------------------------------------------------get_fronto_final_phi
   subroutine get_fronto_final_phi( PHI, t )  
   use module_SchemeSetup, only: isMM
   
   Implicit None
   type(scalar_spatial_data), intent(inout) :: PHI
   real(wp), intent(in)                     :: t
   
   real(wp) :: x1, x2, y1, y2, tol, DX, DY, DXDY
   integer  :: i, j
   
   tol = 1E-13   
   DX  = PHI%Mesh%DX
   DY  = PHI%Mesh%DY
   DXDY = DX * DY     

   ! 1. Calculate Point Values: PHI%P(j,i)
   Do i = 1, PHI%Mesh%NPXtot    
      x1 = PHI%Mesh%Xp(i)   
      Do j = 1, PHI%Mesh%NPYtot 
         y1 = PHI%Mesh%Yp(j)
         PHI%P(j,i) = wrapper_phi( x1, y1 )
      End Do
   End Do
      
   !----------------------------- MM MM MM ---------------------------------------------------
   IF( isMM ) THEN 
      
      ! 2. Calculate Cell Averages: PHI%xyz(j,i) (Integrated over dX and dY)
      Do i = 1, PHI%Mesh%NcellXtot 
         x1 = PHI%Mesh%Xp(i)
         x2 = PHI%Mesh%Xp(i+1)    
         Do j = 1, PHI%Mesh%NcellYtot   
            y1 = PHI%Mesh%Yp(j)
            y2 = PHI%Mesh%Yp(j+1)        
            PHI%xyz(j,i) = AdaptiveNewtonCote2D( wrapper_phi, x1, x2, y1, y2, tol ) / DXDY
         End Do   
      End Do

      ! 3. Calculate Line Averages along Y: PHI%yz(j,i) (Integrated over dY)
      Do i = 1, PHI%Mesh%NPXtot 
         x1 = PHI%Mesh%Xp(i)
         Do j = 1, PHI%Mesh%NcellYtot   
            y1 = PHI%Mesh%Yp(j)
            y2 = PHI%Mesh%Yp(j+1)           
            PHI%yz(j,i) = AdaptiveNewtonCote1D_Y( wrapper_phi, x1, y1, y2, tol ) / DY
         End Do   
      End Do   

      ! 4. Calculate Line Averages along X: PHI%xz(j,i) (Integrated over dX)
      Do i = 1, PHI%Mesh%NcellXtot 
         x1 = PHI%Mesh%Xp(i)
         x2 = PHI%Mesh%Xp(i+1)            
         Do j = 1, PHI%Mesh%NPYtot   
            y1 = PHI%Mesh%Yp(j)           
            PHI%xz(j,i) = AdaptiveNewtonCote1D_X( wrapper_phi, y1, x1, x2, tol ) / DX
         End Do   
      End Do   

   ENDIF    
!----------------------------- MM MM MM ---------------------------------------------------   

 CONTAINS

    ! Non-elemental wrapper to bypass Fortran's dummy argument restrictions on elemental procedures
    pure function wrapper_phi(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp)             :: val
        val = get_analytical_phi(x, y, t)
    end function wrapper_phi
   
   end subroutine get_fronto_final_phi
   
   
   
   
   
  


subroutine get_fronto_initial_Vel( VEL )  
!!   use froto_vortex_functions, only: get_myU, get_myV   
   use module_SchemeSetup, only: isMM
   
   Implicit None
   type(velocity_spatial_data), intent(inout) :: VEL
   real(wp) :: x1,x2,y1,y2, tol ,DX,DY,DXDY
   integer i,j
   
	Do i=1,VEL%Mesh%NPXtot 	
		x1 = VEL%Mesh%Xp(i)	
	 	Do j=1,VEL%Mesh%NPYtot 
		y1 = VEL%Mesh%Yp(j)
 	
	 	 VEL%U%P(j,i) = wrapper_myU( x1,y1 )
	 	 VEL%V%P(j,i) = wrapper_myV( x1,y1 )
 	
	 	ENDDO
 	ENDDO
	 		
!----------------------------- MM MM MM ---------------------------------------------------
   IF( isMM ) THEN 
	   ! AT INITIAL CONDITON PHI(x,y) = -tanh( y/2) 
        tol = 1E-13	   
	DX = VEL%Mesh%DX
	DY = VEL%Mesh%DY
	DXDY = DX*DY     
   !Fxyz 
	Do i=1,VEL%Mesh%NcellXtot 
		x1 = VEL%Mesh%Xp(i)
		x2 = VEL%Mesh%Xp(i+1)	  
		Do j=1,VEL%Mesh%NcellYtot   
			y1 = VEL%Mesh%Yp(j)
			y2 = VEL%Mesh%Yp(j+1)	 
!!			print '(2i5,4F12.6)', i,j, x1,x2,y1,y2       
			VEL%U%xyz(j,i) =  AdaptiveNewtonCote2D(  wrapper_myU, x1,x2,y1,y2,tol)/DXDY
			VEL%V%xyz(j,i) =  AdaptiveNewtonCote2D(  wrapper_myV, x1,x2,y1,y2,tol)/DXDY			
		ENDDO	
	ENDDO
   !Fyz 
	Do i=1,VEL%Mesh%NPXtot 
		x1 = VEL%Mesh%Xp(i)
		Do j=1,VEL%Mesh%NcellYtot   
			y1 = VEL%Mesh%Yp(j)
			y2 = VEL%Mesh%Yp(j+1)	        
			VEL%U%yz(j,i) =  AdaptiveNewtonCote1D_Y(  wrapper_myU, x1,y1,y2,tol)/DY
			VEL%V%yz(j,i) =  AdaptiveNewtonCote1D_Y(  wrapper_myV, x1,y1,y2,tol)/DY			
		ENDDO	
	ENDDO	
   !Fxz 
	Do i=1,VEL%Mesh%NcellXtot 
		x1 = VEL%Mesh%Xp(i)
		x2 = VEL%Mesh%Xp(i+1)			
		Do j=1,VEL%Mesh%NPYtot   
			y1 = VEL%Mesh%Yp(j)	        
			VEL%U%xz(j,i) =  AdaptiveNewtonCote1D_X(  wrapper_myU, y1,x1,x2,tol)/DX
			VEL%V%xz(j,i) =  AdaptiveNewtonCote1D_X(  wrapper_myV, y1,x1,x2,tol)/DX			
		ENDDO	
	ENDDO	
    ENDIF	
!----------------------------- MM MM MM ---------------------------------------------------	

 CONTAINS

    ! Non-elemental wrappers to bypass standard restrictions
    pure function wrapper_myU(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp)             :: val
        val = get_myU(x, y)
    end function wrapper_myU

    pure function wrapper_myV(x, y) result(val)
        real(wp), intent(in) :: x, y
        real(wp)             :: val
        val = get_myV(x, y)
    end function wrapper_myV
   
   end subroutine get_fronto_initial_Vel





SUBROUTINE set_fronto_timedep_BC_simp(PHI, time)

!!    USE froto_vortex_functions, ONLY: wp
!    USE bc_integration_wrapper, ONLY: current_bc_time, phi_time_wrapped, phi_time_wrapped2
    USE module_SchemeSetup,     ONLY: isMM, NG

    IMPLICIT NONE

  
    real (wp), allocatable, dimension(:,:) :: HorX, HorM, HorY, VertX, VertY, VertM !vector for storing vectorize coordinate

    type(scalar_spatial_data), intent(inout) :: PHI
    real(wp),                  intent(in)    :: time

    real(wp) :: x1, x2, y1, y2, t1,t2
    integer  :: i, j
    integer  :: nxp, nyp, ncx, ncy

    ! Synchronize the analytical boundary function with the RK stage time.
    current_bc_time = time

    nxp = PHI%Mesh%NPXtot
    nyp = PHI%Mesh%NPYtot
    ncx = PHI%Mesh%NcellXtot
    ncy = PHI%Mesh%NcellYtot
    
    allocate( HorX( NG+1, nxp), HorY( NG+1, nxp), HorM( NG+1, nxp)  )
    allocate( VertX( nyp ,NG+1), VertY( nyp, NG+1 ) , VertM( nyp, NG+1 ) )    

    !=======================================================================
    ! 1. Point values: exact nodal evaluation
    !=======================================================================
    
!    print *,  size( PHI%Mesh%Xp(:) )

    

    

!    call cpu_time(t1)
    ! Bottom boundary
        DO j = 1, NG + 1
    	HorX(j,:) = PHI%Mesh%Xp(:)
    	HorY(j,:) = PHI%Mesh%Yp(j)    	
       ENDDO
        PHI%p(1:(NG + 1), :) = phi_time_wrapped(   HorX, HorY    )   


    ! Top boundary
       DO j = 1,NG+1
    	HorY(j,:) = PHI%Mesh%Yp(nyp-NG-1+j)    	
       ENDDO
      PHI%p( (nyp - NG): nyp, :) = phi_time_wrapped(   HorX, HorY  )         

    ! Left boundary
    DO i = 1, NG + 1
    	VertX(:,i) = PHI%Mesh%Xp(i)
    	VertY(:,i) = PHI%Mesh%Yp(:)    	
    ENDDO    
    PHI%p ( :, 1:(NG+1)) = phi_time_wrapped(   VertX, VertY )

    ! Right boundary
    DO i = 1, NG + 1
    	VertX(:,i) = PHI%Mesh%Xp(nxp - NG -1 +i)  	
    ENDDO    
    PHI%p(:,(nxp - NG):nxp ) = phi_time_wrapped(   VertX, VertY )

    
!    call cpu_time(t2)
!    write(6,*) "FDM: ", t2-t1  

    IF (.NOT. isMM) RETURN

    !=======================================================================
    ! 2. YZ line averages: integration along y at fixed x
    !=======================================================================
!        call cpu_time(t1)




   ! Bottom boundary
    DO j = 1, NG + 1
	HorX(j,:) = PHI%Mesh%Xp(:)
	HorY(j,:) = PHI%Mesh%Yp(j)    	
    ENDDO

    PHI%yz(1:(NG), 1:nxp) =     (                                         &
              phi_time_wrapped( HorX(1:NG,:),        HorY(1:(NG ),1:nxp)  )            &
   + 4.0_wp * phi_time_wrapped( HorX(1:NG,:),  0.5*( HorY(1:(NG ),1:nxp)+HorY(2:(NG+1),1:nxp)) ) &
        +     phi_time_wrapped( HorX(1:NG,:),                            HorY(2:(NG+1),1:nxp)  )    &
    ) / 6.0_wp 

    PHI%xz(1:(NG + 1), 1:ncx) = simpson_average_x( HorY( 1:(NG + 1), 1:ncx ), HorX( 1:(NG+1), 1:(ncx  ) ), &
                                                                              HorX( 1:(NG+1), 2:(ncx+1)) )          




! Top boundary
   DO j = 1,NG+1
	HorY(j,:) = PHI%Mesh%Yp(nyp-NG-1+j)    	
   ENDDO
   
    
  PHI%yz( (ncy-NG+1):ncy, :) =    (                                         &
              phi_time_wrapped( HorX(1:NG,:),        HorY(1:(NG ),1:nxp)  )            &
   + 4.0_wp * phi_time_wrapped( HorX(1:NG,:),  0.5*( HorY(1:(NG ),1:nxp)+HorY(2:(NG+1),1:nxp)) ) &
        +     phi_time_wrapped( HorX(1:NG,:),                            HorY(2:(NG+1),1:nxp)  )    &
    ) / 6.0_wp         

    PHI%xz((nyp - NG): nyp, 1:ncx ) = simpson_average_x( HorY( 1:(NG + 1), 1:ncx ), &
                                                                          HorX( 1:(NG+1), 1:(ncx  ) ), &
                                                                          HorX( 1:(NG+1), 2:(ncx+1)) )    

    ! Left boundary
    DO i = 1, NG + 1
    	VertX(:,i) = PHI%Mesh%Xp(i)
    	VertY(:,i) = PHI%Mesh%Yp(:)    	
    ENDDO    

    PHI%yz( :, 1:(NG+1) ) = simpson_average_y( VertX(1:ncy, : ), VertY(1:ncy,:), VertY( 2:(ncy+1),:) )   
    PHI%xz( : ,1:NG )     = simpson_average_x( VertY(:, 1:NG ), VertX(:, 1:NG) , VertX(:, 2:(NG+1)))

    ! Right boundary

    DO i = 1, NG + 1
    	VertX(:,i) = PHI%Mesh%Xp(nxp - NG -1 +i)  	
    ENDDO    

    PHI%yz( :,  (nxp - NG): nxp )    = simpson_average_y( VertX(1:ncy,:), VertY(1:ncy,:), VertY(2:(ncy+1),:) )        
    PHI%xz( : , (ncx-ng+1):ncx )     = simpson_average_x( VertY(:, 1:NG ), VertX(:, 1:NG) , VertX(:, 2:(NG+1)))   






    !=======================================================================
    ! 4. XYZ cell averages: tensor-product Simpson rule
    !=======================================================================
  
        
         ! Bottom
        DO j = 1, NG + 1
    	HorX(j,:) = PHI%Mesh%Xp(:)
    	HorY(j,:) = PHI%Mesh%Yp(j)    	
       ENDDO        
        
        PHI%xyz(1:NG, :) = simpson_average_xy( HorX(1:NG,1:NCX), HorX( 1: NG   ,2:(NCX+1) ), &
                                               HorY(1:NG,1:NCX), HorY( 2:(NG+1),2:(NCX+1) )  )
                                               
       ! print *,  PHI%xyz(1:NG, :)     
        ! print *,"--------------------"                                  
    ! Top boundary
       DO j = 1,NG+1
    	HorY(j,:) = PHI%Mesh%Yp(nyp-NG-1+j)    	
       ENDDO        
        PHI%xyz( (NCY-NG+1):NCY, 1:NCX) = simpson_average_xy( HorX(1:NG,1:NCX), HorX( 1: NG   ,2:(NCX+1) ), &
                                                   HorY(1:NG,1:NCX), HorY( 2:(NG+1),1:NCX)    )      
 
    ! Left boundary
    DO i = 1, NG + 1
    	VertX(:,i) = PHI%Mesh%Xp(i)
    	VertY(:,i) = PHI%Mesh%Yp(:)    	
    ENDDO  
          PHI%xyz(1:ncy, 1:NG ) = simpson_average_xy( VertX(1:ncy, 1:NG) , VertX( 1:ncy    , 2:(NG+1) ), &
                                                      VertY(1:Ncy, 1:NG),  VertY( 2:(ncy+1), 1:NG  ) ) 
    ! Right boundary
    DO i = 1, NG + 1
    	VertX(:,i) = PHI%Mesh%Xp(nxp - NG -1 +i)  	
    ENDDO     
          PHI%xyz(1:ncy, 1:NG ) = simpson_average_xy(VertX(1:ncy, 1:NG) , VertX( 1:ncy, 2:(NG+1) ), &
                                                     VertY(1:Ncy, 1:NG),  VertY( 2:(ncy+1), 1:NG ) )        
    ! Bottom and top y boundaries
    !DO i = 1, ncx
    !    x1 = PHI%Mesh%Xp(i)
    !    x2 = PHI%Mesh%Xp(i + 1)
    ! 
    !    PHI%xyz(1:NG, i) = simpson_average_xy(x1, x2, PHI%Mesh%Yp(1:NG),  PHI%Mesh%Yp(2:(NG+1)) )
    !    PHI%xyz( (ncy - NG + 1):ncy, i) = simpson_average_xy(x1, x2, PHI%Mesh%Yp((ncy - NG + 1):ncy), &
    !!                                                                PHI%Mesh%Yp((ncy - NG + 2):(ncy+1) )) 
    !END DO

    ! Left and right x boundaries
 !   DO j = 1, ncy
 !       y1 = PHI%Mesh%Yp(j)
!         y2 = PHI%Mesh%Yp(j + 1)
! 
   !       PHI%xyz(j, 1:NG ) = simpson_average_xy(PHI%Mesh%Xp(1:NG), PHI%Mesh%Xp(2:(NG+1)), y1, y2)
   !       PHI%xyz(j,  (ncx - NG + 1):ncx) = simpson_average_xy(PHI%Mesh%Xp( (ncx - NG + 1):ncx), &
   !                                          PHI%Mesh%Xp( (ncx - NG + 2):(ncx+1) )  , y1, y2)
   ! END DO
    
    !call cpu_time(t2)
    !write(6,*) "FVM: ", t2-t1  
    !stop "CHECK XYZ"  
    
    deallocate( HorX, HorY, VertX, VertY )

CONTAINS

    !=======================================================================
    ! Simpson approximation to the average along x
    !=======================================================================
    elemental FUNCTION simpson_average_x(y, xa, xb) RESULT(f_average)

        real(wp), intent(in) :: y, xa, xb
        real(wp)             :: f_average
        real(wp)             :: xm

        xm = 0.5_wp * (xa + xb)

        f_average = (                                             &
              phi_time_wrapped(xa, y)                             &
            + 4.0_wp * phi_time_wrapped(xm, y)                    &
            + phi_time_wrapped(xb, y)                             &
        ) / 6.0_wp

    END FUNCTION simpson_average_x

    !=======================================================================
    ! Simpson approximation to the average along y
    !=======================================================================
    elemental FUNCTION simpson_average_y(x, ya, yb) RESULT(f_average)

        real(wp), intent(in) :: x, ya, yb
        real(wp)             :: f_average
        real(wp)             :: ym

        ym = 0.5_wp * (ya + yb)

        f_average = (                                             &
                       phi_time_wrapped(x, ya)                             &
            + 4.0_wp * phi_time_wrapped(x, ym)                    &
            +          phi_time_wrapped(x, yb)                             &
        ) / 6.0_wp

    END FUNCTION simpson_average_y

    !=======================================================================
    ! Tensor-product Simpson approximation to a cell average
    !=======================================================================
    elemental FUNCTION simpson_average_xy(xa, xb, ya, yb) RESULT(f_average)

        real(wp), intent(in) :: xa, xb, ya, yb
        real(wp)             :: f_average
        real(wp)             :: xm, ym

        xm = 0.5_wp * (xa + xb)
        ym = 0.5_wp * (ya + yb)

        f_average = (                                             &
              phi_time_wrapped(xa, ya)                            &
            + 4.0_wp * phi_time_wrapped(xm, ya)                   &
            + phi_time_wrapped(xb, ya)                            &
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            + 4.0_wp * phi_time_wrapped(xa, ym)                   &
            + 16.0_wp * phi_time_wrapped(xm, ym)                  &
            + 4.0_wp * phi_time_wrapped(xb, ym)                   &
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            +          phi_time_wrapped(xa, yb)                            &
            + 4.0_wp * phi_time_wrapped(xm, yb)                   &
            +          phi_time_wrapped(xb, yb)                            &
        ) / 36.0_wp


 
    END FUNCTION simpson_average_xy

END SUBROUTINE set_fronto_timedep_BC_simp





 



!!!!!!!!!!!!!!! FOR FDM  FDM  FDM FDM FDM FDM FDM FDM FDM FDM

  subroutine Compute_FirstDerivativeX_FD( DFDX, Fp,  iMesh, COEF, NX, NY )
  
  use module_SchemeSetup, only: isMM, SptScheme 
  type(Coefficient_Differentiation), intent(in) :: COEF
  type(MeshDomain), intent(in) :: iMesh
  integer , intent(in) :: NX, NY  ! Dimension of the result The source will have DIM:NX-1, but same NY
 
  real (wp), intent(in) ::      Fp( NY, NX ) 
  real (wp), intent(out) :: DFDX( NY, NX ) 

  integer istart, istop, jstart, jstop, i, j, if0, ifL, jf0, jfL
 
! Local index variables for optimized vector slicing (X and Y axes)


      INTEGER :: iF0M1, iFLM1, iF0M2, iFLM2, iF0M3, iFLM3, &
                 iF0M4, iFLM4, iF0M5, iFLM5, iF0M6, iFLM6
                 
      INTEGER :: iF0P1, iFLP1, iF0P2, iFLP2, iF0P3, iFLP3, &
                 iF0P4, iFLP4, iF0P5, iFLP5, iF0P6, iFLP6 
 
  
  ! This will over compute the cell by 1 row, but the code is much easier to maintain
  !  comparing to split the code, or add IF-clause
  if0=iMesh%IfSTART
  ifL=iMesh%IfEND
    
  jf0=iMesh%JfSTART
  jfL=iMesh%JfEND
  
!---------------------------------------------- FDM METHOD   IF:LINE 681  
 
 ! Pre-calculate array slice bounds to eliminate loop branching overhead
      iF0M1 = iF0 - 1;   iFLM1 = iFL - 1
      iF0M2 = iF0 - 2;   iFLM2 = iFL - 2
      iF0M3 = iF0 - 3;   iFLM3 = iFL - 3
      iF0M4 = iF0 - 4;   iFLM4 = iFL - 4
      iF0M5 = iF0 - 5;   iFLM5 = iFL - 5
      iF0M6 = iF0 - 6;   iFLM6 = iFL - 6

      iF0P1 = iF0 + 1;   iFLP1 = iFL + 1
      iF0P2 = iF0 + 2;   iFLP2 = iFL + 2
      iF0P3 = iF0 + 3;   iFLP3 = iFL + 3
      iF0P4 = iF0 + 4;   iFLP4 = iFL + 4
      iF0P5 = iF0 + 5;   iFLP5 = iFL + 5
      iF0P6 = iF0 + 6;   iFLP6 = iFL + 6

      SELECT CASE ( SptScheme ) 
         CASE( 7002 )
            DFDX(jf0:jfL, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) 
		 CASE( 7004 )
            DFDX(jf0:jfL, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2)

         CASE( 7006 )
            DFDX(jf0:jfL, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
               COEF%COEF_FDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3)   

         CASE( 7008 )
            DFDX(jf0:jfL, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(4) * Fp(jf0:jfL, iF0M4:iFLM4) + &
               COEF%COEF_FDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
               COEF%COEF_FDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3) + &
               COEF%COEF_FDX_B2_R(4) * Fp(jf0:jfL, iF0P4:iFLP4)   

         CASE( 7010 )
            DFDX(jf0:jfL, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(5) * Fp(jf0:jfL, iF0M5:iFLM5) + &
               COEF%COEF_FDX_B2_L(4) * Fp(jf0:jfL, iF0M4:iFLM4) + &
               COEF%COEF_FDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
               COEF%COEF_FDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3) + &
               COEF%COEF_FDX_B2_R(4) * Fp(jf0:jfL, iF0P4:iFLP4) + &
               COEF%COEF_FDX_B2_R(5) * Fp(jf0:jfL, iF0P5:iFLP5)   

         CASE( 7012 )
            DFDX(jf0:jfL, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(6) * Fp(jf0:jfL, iF0M6:iFLM6) + &
               COEF%COEF_FDX_B2_L(5) * Fp(jf0:jfL, iF0M5:iFLM5) + &
               COEF%COEF_FDX_B2_L(4) * Fp(jf0:jfL, iF0M4:iFLM4) + &
               COEF%COEF_FDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
               COEF%COEF_FDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3) + &
               COEF%COEF_FDX_B2_R(4) * Fp(jf0:jfL, iF0P4:iFLP4) + &
               COEF%COEF_FDX_B2_R(5) * Fp(jf0:jfL, iF0P5:iFLP5) + &
               COEF%COEF_FDX_B2_R(6) * Fp(jf0:jfL, iF0P6:iFLP6)   

         CASE DEFAULT
            ERROR STOP 'Unknown Spatial Scheme: Compute_FirstDerivativeX_FD'
      END SELECT		        	 				 	 		 	
	  
    	
!-------------------------- END FDM  FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM

  END   subroutine Compute_FirstDerivativeX_FD











  
subroutine Compute_FirstDerivativeY_FD( DFDY,  Fp, iMesh, COEF, NX, NY )
  use module_SchemeSetup, only: isMM, SptScheme 
  type(MeshDomain), intent(in) :: iMesh
  type(Coefficient_Differentiation), intent(in) :: COEF  
  integer , intent(in) :: NX,NY ! Dimension of the result The cell source will have DIM:NY-1, but same  NX

  real (wp),  intent(in) ::      Fp( NY  , NX )
  real (wp),  intent(inout) :: DFDY( NY  , NX )
  
  integer istart, istop, jstart, jstop, i, j, if0, ifL, jf0, jfL
! Local index variables for optimized vector slicing
      INTEGER :: jF0M1, jFLM1, jF0M2, jFLM2, jF0M3, jFLM3, &
                 jF0M4, jFLM4, jF0M5, jFLM5, jF0M6, jFLM6
                 
      INTEGER :: jF0P1, jFLP1, jF0P2, jFLP2, jF0P3, jFLP3, &
                 jF0P4, jFLP4, jF0P5, jFLP5, jF0P6, jFLP6

 
  if0=iMesh%IfSTART
  ifL=iMesh%IfEND

  jf0=iMesh%JfSTART
  jfL=iMesh%JfEND
 
! Define and collect slice index tracking variables
      jF0M1 = jF0 - 1;   jFLM1 = jFL - 1
      jF0M2 = jF0 - 2;   jFLM2 = jFL - 2
      jF0M3 = jF0 - 3;   jFLM3 = jFL - 3
      jF0M4 = jF0 - 4;   jFLM4 = jFL - 4
      jF0M5 = jF0 - 5;   jFLM5 = jFL - 5
      jF0M6 = jF0 - 6;   jFLM6 = jFL - 6

      jF0P1 = jF0 + 1;   jFLP1 = jFL + 1
      jF0P2 = jF0 + 2;   jFLP2 = jFL + 2
      jF0P3 = jF0 + 3;   jFLP3 = jFL + 3
      jF0P4 = jF0 + 4;   jFLP4 = jFL + 4
      jF0P5 = jF0 + 5;   jFLP5 = jFL + 5
      jF0P6 = jF0 + 6;   jFLP6 = jFL + 6

      SELECT CASE ( SptScheme ) 
         CASE( 7002 )
            DFDY(jF0:jFL, if0:ifL) = &
               COEF%COEF_FDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) 
         CASE( 7004 )
            DFDY(jF0:jFL, if0:ifL) = &
               COEF%COEF_FDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
               COEF%COEF_FDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL)

         CASE( 7006 )
            DFDY(jF0:jFL, if0:ifL) = &
               COEF%COEF_FDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
               COEF%COEF_FDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
               COEF%COEF_FDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
               COEF%COEF_FDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL)   

         CASE( 7008 )
            DFDY(jF0:jFL, if0:ifL) = &
               COEF%COEF_FDY_B2_L(4) * Fp(jF0M4:jFLM4, if0:ifL) + &
               COEF%COEF_FDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
               COEF%COEF_FDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
               COEF%COEF_FDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
               COEF%COEF_FDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL) + &
               COEF%COEF_FDY_B2_R(4) * Fp(jF0P4:jFLP4, if0:ifL)   

         CASE( 7010 )
            DFDY(jF0:jFL, if0:ifL) = &
               COEF%COEF_FDY_B2_L(5) * Fp(jF0M5:jFLM5, if0:ifL) + &     
               COEF%COEF_FDY_B2_L(4) * Fp(jF0M4:jFLM4, if0:ifL) + &
               COEF%COEF_FDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
               COEF%COEF_FDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
               COEF%COEF_FDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
               COEF%COEF_FDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL) + &
               COEF%COEF_FDY_B2_R(4) * Fp(jF0P4:jFLP4, if0:ifL) + &
               COEF%COEF_FDY_B2_R(5) * Fp(jF0P5:jFLP5, if0:ifL)   

         CASE( 7012 )
            DFDY(jF0:jFL, if0:ifL) = &
               COEF%COEF_FDY_B2_L(6) * Fp(jF0M6:jFLM6, if0:ifL) + &    
               COEF%COEF_FDY_B2_L(5) * Fp(jF0M5:jFLM5, if0:ifL) + &     
               COEF%COEF_FDY_B2_L(4) * Fp(jF0M4:jFLM4, if0:ifL) + &
               COEF%COEF_FDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
               COEF%COEF_FDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
               COEF%COEF_FDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
               COEF%COEF_FDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
               COEF%COEF_FDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL) + &
               COEF%COEF_FDY_B2_R(4) * Fp(jF0P4:jFLP4, if0:ifL) + &
               COEF%COEF_FDY_B2_R(5) * Fp(jF0P5:jFLP5, if0:ifL) + &
               COEF%COEF_FDY_B2_R(6) * Fp(jF0P6:jFLP6, if0:ifL)        

         CASE DEFAULT
            ERROR STOP 'Unknown Spatial Scheme: Compute_FirstDerivativeY_FD'
      END SELECT 
      
  END   subroutine Compute_FirstDerivativeY_FD   
  
  
  
  
  
  


 
subroutine Compute_SecondDerivativeX_FD( DFDX, Fp, iMesh, COEF, NX, NY )
  
  use module_SchemeSetup, only: isMM, SptScheme 
  type(Coefficient_Differentiation), intent(in) :: COEF
  type(MeshDomain), intent(in) :: iMesh
  integer , intent(in) :: NX, NY  ! Dimension of the result

  real (wp), intent(in)  :: Fp( NY, NX ) 
  real (wp), intent(out) :: DFDX( NY, NX ) 

  integer :: if0, ifL, jf0, jfL

  ! Local index variables for optimized vector slicing (X and Y axes)
  integer :: iF0M1, iFLM1, iF0M2, iFLM2, iF0M3, iFLM3, &
             iF0M4, iFLM4, iF0M5, iFLM5, iF0M6, iFLM6
             
  integer :: iF0P1, iFLP1, iF0P2, iFLP2, iF0P3, iFLP3, &
             iF0P4, iFLP4, iF0P5, iFLP5, iF0P6, iFLP6 

  if0 = iMesh%IfSTART
  ifL = iMesh%IfEND
    
  jf0 = iMesh%JfSTART
  jfL = iMesh%JfEND
  
  ! Pre-calculate array slice bounds to eliminate loop branching overhead
  iF0M1 = iF0 - 1;   iFLM1 = iFL - 1
  iF0M2 = iF0 - 2;   iFLM2 = iFL - 2
  iF0M3 = iF0 - 3;   iFLM3 = iFL - 3
  iF0M4 = iF0 - 4;   iFLM4 = iFL - 4
  iF0M5 = iF0 - 5;   iFLM5 = iFL - 5
  iF0M6 = iF0 - 6;   iFLM6 = iFL - 6

  iF0P1 = iF0 + 1;   iFLP1 = iFL + 1
  iF0P2 = iF0 + 2;   iFLP2 = iFL + 2
  iF0P3 = iF0 + 3;   iFLP3 = iFL + 3
  iF0P4 = iF0 + 4;   iFLP4 = iFL + 4
  iF0P5 = iF0 + 5;   iFLP5 = iFL + 5
  iF0P6 = iF0 + 6;   iFLP6 = iFL + 6

  SELECT CASE ( SptScheme ) 
     CASE( 7002 )
        DFDX(jf0:jfL, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Fp(jf0:jfL, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) 
	 CASE( 7004 )
        DFDX(jf0:jfL, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Fp(jf0:jfL, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2)

     CASE( 7006 )
        DFDX(jf0:jfL, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Fp(jf0:jfL, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
           COEF%COEF_SDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3)   

     CASE( 7008 )
        DFDX(jf0:jfL, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(4) * Fp(jf0:jfL, iF0M4:iFLM4) + &
           COEF%COEF_SDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Fp(jf0:jfL, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
           COEF%COEF_SDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3) + &
           COEF%COEF_SDX_B2_R(4) * Fp(jf0:jfL, iF0P4:iFLP4)   

     CASE( 7010 )
        DFDX(jf0:jfL, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(5) * Fp(jf0:jfL, iF0M5:iFLM5) + &
           COEF%COEF_SDX_B2_L(4) * Fp(jf0:jfL, iF0M4:iFLM4) + &
           COEF%COEF_SDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Fp(jf0:jfL, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
           COEF%COEF_SDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3) + &
           COEF%COEF_SDX_B2_R(4) * Fp(jf0:jfL, iF0P4:iFLP4) + &
           COEF%COEF_SDX_B2_R(5) * Fp(jf0:jfL, iF0P5:iFLP5)   

     CASE( 7012 )
        DFDX(jf0:jfL, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(6) * Fp(jf0:jfL, iF0M6:iFLM6) + &
           COEF%COEF_SDX_B2_L(5) * Fp(jf0:jfL, iF0M5:iFLM5) + &
           COEF%COEF_SDX_B2_L(4) * Fp(jf0:jfL, iF0M4:iFLM4) + &
           COEF%COEF_SDX_B2_L(3) * Fp(jf0:jfL, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B2_L(2) * Fp(jf0:jfL, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Fp(jf0:jfL, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Fp(jf0:jfL, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Fp(jf0:jfL, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Fp(jf0:jfL, iF0P2:iFLP2) + &
           COEF%COEF_SDX_B2_R(3) * Fp(jf0:jfL, iF0P3:iFLP3) + &
           COEF%COEF_SDX_B2_R(4) * Fp(jf0:jfL, iF0P4:iFLP4) + &
           COEF%COEF_SDX_B2_R(5) * Fp(jf0:jfL, iF0P5:iFLP5) + &
           COEF%COEF_SDX_B2_R(6) * Fp(jf0:jfL, iF0P6:iFLP6)   

     CASE DEFAULT
        ERROR STOP 'Unknown Spatial Scheme: Compute_SecondDerivativeX_FD'
  END SELECT                       

end subroutine Compute_SecondDerivativeX_FD

 
subroutine Compute_SecondDerivativeY_FD( DFDY, Fp, iMesh, COEF, NX, NY )
  use module_SchemeSetup, only: isMM, SptScheme 
  type(MeshDomain), intent(in) :: iMesh
  type(Coefficient_Differentiation), intent(in) :: COEF  
  integer , intent(in) :: NX, NY ! Dimension of the result

  real (wp),  intent(in)    :: Fp( NY, NX )
  real (wp),  intent(inout) :: DFDY( NY, NX )
  
  integer :: if0, ifL, jf0, jfL

  ! Local index variables for optimized vector slicing
  integer :: jF0M1, jFLM1, jF0M2, jFLM2, jF0M3, jFLM3, &
             jF0M4, jFLM4, jF0M5, jFLM5, jF0M6, jFLM6
             
  integer :: jF0P1, jFLP1, jF0P2, jFLP2, jF0P3, jFLP3, &
             jF0P4, jFLP4, jF0P5, jFLP5, jF0P6, jFLP6

  if0 = iMesh%IfSTART
  ifL = iMesh%IfEND

  jf0 = iMesh%JfSTART
  jfL = iMesh%JfEND
 
  ! Define and collect slice index tracking variables
  jF0M1 = jF0 - 1;   jFLM1 = jFL - 1
  jF0M2 = jF0 - 2;   jFLM2 = jFL - 2
  jF0M3 = jF0 - 3;   jFLM3 = jFL - 3
  jF0M4 = jF0 - 4;   jFLM4 = jFL - 4
  jF0M5 = jF0 - 5;   jFLM5 = jFL - 5
  jF0M6 = jF0 - 6;   jFLM6 = jFL - 6

  jF0P1 = jF0 + 1;   jFLP1 = jFL + 1
  jF0P2 = jF0 + 2;   jFLP2 = jFL + 2
  jF0P3 = jF0 + 3;   jFLP3 = jFL + 3
  jF0P4 = jF0 + 4;   jFLP4 = jFL + 4
  jF0P5 = jF0 + 5;   jFLP5 = jFL + 5
  jF0P6 = jF0 + 6;   jFLP6 = jFL + 6

  SELECT CASE ( SptScheme ) 
     CASE( 7002 )
        DFDY(jF0:jFL, if0:ifL) = &
           COEF%COEF_SDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
           COEF%COEF_SDY_B2_M    * Fp(jF0  :jFL  , if0:ifL) + &
           COEF%COEF_SDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) 
	 CASE( 7004 )
        DFDY(jF0:jFL, if0:ifL) = &
           COEF%COEF_SDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
           COEF%COEF_SDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
           COEF%COEF_SDY_B2_M    * Fp(jF0  :jFL  , if0:ifL) + &
           COEF%COEF_SDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
           COEF%COEF_SDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL)

     CASE( 7006 )
        DFDY(jF0:jFL, if0:ifL) = &
           COEF%COEF_SDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
           COEF%COEF_SDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
           COEF%COEF_SDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
           COEF%COEF_SDY_B2_M    * Fp(jF0  :jFL  , if0:ifL) + &
           COEF%COEF_SDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
           COEF%COEF_SDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
           COEF%COEF_SDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL)   

     CASE( 7008 )
        DFDY(jF0:jFL, if0:ifL) = &
           COEF%COEF_SDY_B2_L(4) * Fp(jF0M4:jFLM4, if0:ifL) + &
           COEF%COEF_SDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
           COEF%COEF_SDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
           COEF%COEF_SDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
           COEF%COEF_SDY_B2_M    * Fp(jF0  :jFL  , if0:ifL) + &
           COEF%COEF_SDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
           COEF%COEF_SDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
           COEF%COEF_SDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL) + &
           COEF%COEF_SDY_B2_R(4) * Fp(jF0P4:jFLP4, if0:ifL)   

     CASE( 7010 )
        DFDY(jF0:jFL, if0:ifL) = &
           COEF%COEF_SDY_B2_L(5) * Fp(jF0M5:jFLM5, if0:ifL) + &
           COEF%COEF_SDY_B2_L(4) * Fp(jF0M4:jFLM4, if0:ifL) + &
           COEF%COEF_SDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
           COEF%COEF_SDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
           COEF%COEF_SDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
           COEF%COEF_SDY_B2_M    * Fp(jF0  :jFL  , if0:ifL) + &
           COEF%COEF_SDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
           COEF%COEF_SDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
           COEF%COEF_SDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL) + &
           COEF%COEF_SDY_B2_R(4) * Fp(jF0P4:jFLP4, if0:ifL) + &
           COEF%COEF_SDY_B2_R(5) * Fp(jF0P5:jFLP5, if0:ifL)   

     CASE( 7012 )
        DFDY(jF0:jFL, if0:ifL) = &
           COEF%COEF_SDY_B2_L(6) * Fp(jF0M6:jFLM6, if0:ifL) + &
           COEF%COEF_SDY_B2_L(5) * Fp(jF0M5:jFLM5, if0:ifL) + &
           COEF%COEF_SDY_B2_L(4) * Fp(jF0M4:jFLM4, if0:ifL) + &
           COEF%COEF_SDY_B2_L(3) * Fp(jF0M3:jFLM3, if0:ifL) + &
           COEF%COEF_SDY_B2_L(2) * Fp(jF0M2:jFLM2, if0:ifL) + &
           COEF%COEF_SDY_B2_L(1) * Fp(jF0M1:jFLM1, if0:ifL) + &
           COEF%COEF_SDY_B2_M    * Fp(jF0  :jFL  , if0:ifL) + &
           COEF%COEF_SDY_B2_R(1) * Fp(jF0P1:jFLP1, if0:ifL) + &
           COEF%COEF_SDY_B2_R(2) * Fp(jF0P2:jFLP2, if0:ifL) + &
           COEF%COEF_SDY_B2_R(3) * Fp(jF0P3:jFLP3, if0:ifL) + &
           COEF%COEF_SDY_B2_R(4) * Fp(jF0P4:jFLP4, if0:ifL) + &
           COEF%COEF_SDY_B2_R(5) * Fp(jF0P5:jFLP5, if0:ifL) + &
           COEF%COEF_SDY_B2_R(6) * Fp(jF0P6:jFLP6, if0:ifL)   

     CASE DEFAULT
        ERROR STOP 'Unknown Spatial Scheme: Compute_SecondDerivativeY_FD'
  END SELECT 
  
END subroutine Compute_SecondDerivativeY_FD  
  
  
  
  
  
  
!-------------------------- END FDM  FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM
!-------------------------- END FDM  FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM FDM



 

   
!  MM X-Derivatives     MM X-Derivatives MM X-Derivatives MM X-Derivatives MM X-Derivatives 
!  MM X-Derivatives     MM X-Derivatives MM X-Derivatives MM X-Derivatives MM X-Derivatives 
!  MM X-Derivatives     MM X-Derivatives MM X-Derivatives MM X-Derivatives MM X-Derivatives 
  subroutine Compute_FirstDerivativeX_YZ( DFDX, Uxz, Up,  iMesh, COEF, NX, NY )

  use module_SchemeSetup, only: isMM, SptScheme , NG
  type(Coefficient_Differentiation), intent(in) :: COEF
  type(MeshDomain), intent(in) :: iMesh
  integer , intent(in) :: NX, NY  ! Dimension of the result The source will have DIM:NX-1, but same NY
  real (wp), intent(in) ::     Uxz( NY, NX-1 ) 
  real (wp), intent(in) ::      Up( NY, NX ) 
  real (wp), intent(inout) :: DFDX( NY, NX ) 

  integer istart, istop, jstart, jstop, i, j, if0, ifL, jf0, jfL, jC0, jCL, lpos
 
 ! Local index variables for optimized X-vector tracking slices
      INTEGER :: iF0M1, iFLM1, iF0M2, iFLM2, iF0M3, iFLM3, &
                 iF0M4, iFLM4, iF0M5, iFLM5, iF0M6, iFLM6
                 
      INTEGER :: iF0P1, iFLP1, iF0P2, iFLP2, iF0P3, iFLP3, &
                 iF0P4, iFLP4, iF0P5, iFLP5, iF0P6, iFLP6
  
  ! This will over compute the cell by 1 row, but the code is much easier to maintain
  !  comparing to split the code, or add IF-clause
  jC0=iMesh%JcSTART
  jCL=iMesh%JcEND
  
  if0=iMesh%IfSTART 
  ifL=iMesh%IfEND
  
 ! Since NY is varied by the input type, the last position will be varied by NY as well
 ! The last position is thus NY - NG
 lpos = NY - NG
  
 
 
! The index boundary assignments needed for this block:
       iF0M1 = iF0 - 1;   iFLM1 = iFL - 1
       iF0M2 = iF0 - 2;   iFLM2 = iFL - 2
       iF0M3 = iF0 - 3;   iFLM3 = iFL - 3
       iF0P1 = iF0 + 1;   iFLP1 = iFL + 1
       iF0P2 = iF0 + 2;   iFLP2 = iFL + 2
       iF0P3 = iF0 + 3;   iFLP3 = iFL + 3

 
 
      SELECT CASE ( SptScheme )   
         CASE( 3004 )
            DFDX(jC0:lpos, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
               ! 
               COEF%COEF_FDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B1_R(1) * Uxz(jC0:lpos, iF0:iFL)

 

         CASE( 3006 )
            DFDX(jC0:lpos, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
               ! 
               COEF%COEF_FDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &              
               COEF%COEF_FDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL    ) + &
               COEF%COEF_FDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1)

         CASE( 3008 )
            DFDX(jC0:lpos, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(2) * Up(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Up(jC0:lpos, iF0P2:iFLP2) + &
               ! 
               COEF%COEF_FDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL) + &
               COEF%COEF_FDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1)

         CASE( 3010 )
            DFDX(jC0:lpos, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(2) * Up(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Up(jC0:lpos, iF0P2:iFLP2) + &
               COEF%COEF_FDX_B1_L(3) * Uxz(jC0:lpos, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B1_R(1) * Uxz(jC0:lpos, iF0:iFL) + &
               COEF%COEF_FDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B1_R(3) * Uxz(jC0:lpos, iF0P2:iFLP2)

         CASE( 3012 )
            DFDX(jC0:lpos, iF0:iFL) = &
               COEF%COEF_FDX_B2_L(3) * Up(jC0:lpos, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B2_L(2) * Up(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
               COEF%COEF_FDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B2_R(2) * Up(jC0:lpos, iF0P2:iFLP2) + &
               COEF%COEF_FDX_B2_R(3) * Up(jC0:lpos, iF0P3:iFLP3) + &
               !
               COEF%COEF_FDX_B1_L(3) * Uxz(jC0:lpos, iF0M3:iFLM3) + &
               COEF%COEF_FDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
               COEF%COEF_FDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
               !
               COEF%COEF_FDX_B1_R(1) * Uxz(jC0:lpos, iF0:iFL) + &
               COEF%COEF_FDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1) + &
               COEF%COEF_FDX_B1_R(3) * Uxz(jC0:lpos, iF0P2:iFLP2)

         CASE DEFAULT
            ERROR STOP 'Unknown Spatial Scheme: Compute_FirstDerivativeX_YZ'
      END SELECT

 
 
!-------------------------- END Multi-moment method  

  END   subroutine Compute_FirstDerivativeX_YZ 
   
   
   
  
   
!MM-Y-Derivatives MM-Y-Derivatives    MM-Y-Derivatives    MM-Y-Derivatives    MM-Y-Derivatives       
!MM-Y-Derivatives MM-Y-Derivatives    MM-Y-Derivatives    MM-Y-Derivatives    MM-Y-Derivatives 
!MM-Y-Derivatives MM-Y-Derivatives    MM-Y-Derivatives    MM-Y-Derivatives    MM-Y-Derivatives 
  subroutine Compute_FirstDerivativeY_XZ( DFDY, Uyz, Up, iMesh, COEF, NX, NY )
  use module_SchemeSetup, only: isMM, SptScheme , NG
  
  IMPLICIT NONE
  type(MeshDomain), intent(in) :: iMesh
  type(Coefficient_Differentiation), intent(in) :: COEF  
 
  integer , intent(in) :: NX,NY ! Dimension of the result The cell source will have DIM:NY-1, but same  NX
  real (wp),  intent(in) ::     Uyz( NY-1, NX )
  real (wp),  intent(in) ::      Up( NY  , NX )
  real (wp),  intent(inout) :: DFDY( NY  , NX )
  
  integer istart, istop, jstart, jstop, i, j, if0, ifL, iC0, iCL, jf0, jfL, ilpos, jlpos, lpos, lposP1, iC0M1

 ! Local index variables for optimized Y-vector tracking slices
      INTEGER :: jF0M1, jFLM1, jF0M2, jFLM2, jF0M3, jFLM3, &
                 jF0M4, jFLM4, jF0M5, jFLM5, jF0M6, jFLM6
                 
      INTEGER :: jF0P1, jFLP1, jF0P2, jFLP2, jF0P3, jFLP3, &
                 jF0P4, jFLP4, jF0P5, jFLP5, jF0P6, jFLP6
 
 
  jf0=iMesh%JfSTART
  jfL=iMesh%JfEND 
  
  if0=iMesh%IfSTART
  ifL=iMesh%IfEND 
  
  iC0=iMesh%IcSTART
  iCL=iMesh%IcEND   
  iC0M1=IC0 - 1
 
 ! Since NX is varied by the input type, the last position will be varied by NX as well
 ! The last position is thus NX - NG

 lpos= NX - NG 
 lposP1=lpos + 1
 
 
! Assumes these index boundary assignments are calculated prior to this block:
       jF0M1 = jF0 - 1;   jFLM1 = jFL - 1
       jF0M2 = jF0 - 2;   jFLM2 = jFL - 2
       jF0M3 = jF0 - 3;   jFLM3 = jFL - 3
       jF0P1 = jF0 + 1;   jFLP1 = jFL + 1
       jF0P2 = jF0 + 2;   jFLP2 = jFL + 2
       jF0P3 = jF0 + 3;   jFLP3 = jFL + 3

     SELECT CASE ( SptScheme )   
         CASE( 3004 )
            DFDY(jF0:jFL, iC0M1:lposP1) = &
               COEF%COEF_FDY_B2_L(1) *  Up(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(1) *  Up(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(1) * Uyz(jF0:jFL,     iC0M1:lposP1)

         CASE( 3006 )
            DFDY(jF0:jFL, iC0M1:lposP1) = &
               COEF%COEF_FDY_B2_L(1) *  Up(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(1) *  Up(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(1) * Uyz(jF0:jFL,     iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1)

         CASE( 3008 )
            DFDY(jF0:jFL, iC0M1:lposP1) = &
               COEF%COEF_FDY_B2_L(2) *  Up(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_L(1) *  Up(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(1) *  Up(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(2) *  Up(jF0P2:jFLP2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(1) * Uyz(jF0:jFL,     iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1)

         CASE( 3010 )
            DFDY(jF0:jFL, iC0M1:lposP1) = &
               COEF%COEF_FDY_B2_L(2) *  Up(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_L(1) *  Up(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(1) *  Up(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(2) *  Up(jF0P2:jFLP2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(3) * Uyz(jF0M3:jFLM3, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(1) * Uyz(jF0:jFL,     iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(3) * Uyz(jF0P2:jFLP2, iC0M1:lposP1)

         CASE( 3012 )
            DFDY(jF0:jFL, iC0M1:lposP1) = &
               COEF%COEF_FDY_B2_L(3) *  Up(jF0M3:jFLM3, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_L(2) *  Up(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_L(1) *  Up(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(1) *  Up(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(2) *  Up(jF0P2:jFLP2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B2_R(3) *  Up(jF0P3:jFLP3, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(3) * Uyz(jF0M3:jFLM3, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(1) * Uyz(jF0:jFL,     iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1) + &
               COEF%COEF_FDY_B1_R(3) * Uyz(jF0P2:jFLP2, iC0M1:lposP1)

         CASE DEFAULT
            ERROR STOP 'Unknown Spatial Scheme: Compute_FirstDerivativeY_XZ'
      END SELECT 
 
  
 		
  END   subroutine Compute_FirstDerivativeY_XZ
 
   
   
   
   
   
!MM-X-SECOND Derivatives MM-X-SECOND Derivatives MM-X-SECOND Derivatives  MM-X-SECOND Derivatives   MM-X-SECOND Derivatives      
!MM-X-SECOND Derivatives MM-X-SECOND Derivatives MM-X-SECOND Derivatives  MM-X-SECOND Derivatives   MM-X-SECOND Derivatives      
!MM-X-SECOND Derivatives MM-X-SECOND Derivatives MM-X-SECOND Derivatives  MM-X-SECOND Derivatives   MM-X-SECOND Derivatives  

 
subroutine Compute_SecondDerivativeX_YZ( DFDX, Uxz, Up, iMesh, COEF, NX, NY )
  
  use module_SchemeSetup, only: isMM, SptScheme, NG
  type(Coefficient_Differentiation), intent(in) :: COEF
  type(MeshDomain), intent(in) :: iMesh
  integer , intent(in) :: NX, NY  ! Dimension of the result
  real (wp), intent(in)    :: Uxz( NY, NX-1 ) 
  real (wp), intent(in)    :: Up( NY, NX ) 
  real (wp), intent(inout) :: DFDX( NY, NX ) 

  integer :: if0, ifL, jC0, lpos
 
  ! Local index variables for optimized X-vector tracking slices
  integer :: iF0M1, iFLM1, iF0M2, iFLM2, iF0M3, iFLM3
  integer :: iF0P1, iFLP1, iF0P2, iFLP2, iF0P3, iFLP3
  
  jC0 = iMesh%JcSTART
  if0 = iMesh%IfSTART 
  ifL = iMesh%IfEND
  
  lpos = NY - NG
  
  ! Index boundary assignments
  iF0M1 = iF0 - 1;   iFLM1 = iFL - 1
  iF0M2 = iF0 - 2;   iFLM2 = iFL - 2
  iF0M3 = iF0 - 3;   iFLM3 = iFL - 3
  iF0P1 = iF0 + 1;   iFLP1 = iFL + 1
  iF0P2 = iF0 + 2;   iFLP2 = iFL + 2
  iF0P3 = iF0 + 3;   iFLP3 = iFL + 3

  SELECT CASE ( SptScheme )   
     CASE( 3004 )
        DFDX(jC0:lpos, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Up(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
           !
           COEF%COEF_SDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL  )

     CASE( 3006 )
        DFDX(jC0:lpos, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Up(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
           !
           COEF%COEF_SDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1)

     CASE( 3008 )
        DFDX(jC0:lpos, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(2) * Up(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Up(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Up(jC0:lpos, iF0P2:iFLP2) + &
           !
           COEF%COEF_SDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1)

     CASE( 3010 )
        DFDX(jC0:lpos, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(2) * Up(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Up(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Up(jC0:lpos, iF0P2:iFLP2) + &
           !
           COEF%COEF_SDX_B1_L(3) * Uxz(jC0:lpos, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B1_R(3) * Uxz(jC0:lpos, iF0P2:iFLP2)

     CASE( 3012 )
        DFDX(jC0:lpos, iF0:iFL) = &
           COEF%COEF_SDX_B2_L(3) * Up(jC0:lpos, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B2_L(2) * Up(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B2_L(1) * Up(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B2_M    * Up(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B2_R(1) * Up(jC0:lpos, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B2_R(2) * Up(jC0:lpos, iF0P2:iFLP2) + &
           COEF%COEF_SDX_B2_R(3) * Up(jC0:lpos, iF0P3:iFLP3) + &
           !
           COEF%COEF_SDX_B1_L(3) * Uxz(jC0:lpos, iF0M3:iFLM3) + &
           COEF%COEF_SDX_B1_L(2) * Uxz(jC0:lpos, iF0M2:iFLM2) + &
           COEF%COEF_SDX_B1_L(1) * Uxz(jC0:lpos, iF0M1:iFLM1) + &
           COEF%COEF_SDX_B1_R(1) * Uxz(jC0:lpos, iF0  :iFL  ) + &
           COEF%COEF_SDX_B1_R(2) * Uxz(jC0:lpos, iF0P1:iFLP1) + &
           COEF%COEF_SDX_B1_R(3) * Uxz(jC0:lpos, iF0P2:iFLP2)

     CASE DEFAULT
        ERROR STOP 'Unknown Spatial Scheme: Compute_SecondDerivativeX_YZ'
  END SELECT

END subroutine Compute_SecondDerivativeX_YZ

!MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives
!MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives
!MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives MM-Y-SECOND Derivatives
 
 
subroutine Compute_SecondDerivativeY_XZ( DFDY, Uyz, Up, iMesh, COEF, NX, NY )
  use module_SchemeSetup, only: isMM, SptScheme, NG
  
  implicit none
  type(MeshDomain), intent(in) :: iMesh
  type(Coefficient_Differentiation), intent(in) :: COEF  
 
  integer , intent(in) :: NX, NY ! Dimension of the result
  real (wp),  intent(in)    :: Uyz( NY-1, NX )
  real (wp),  intent(in)    :: Up( NY  , NX )
  real (wp),  intent(inout) :: DFDY( NY  , NX )
  
  integer :: jf0, jfL, iC0, lpos, lposP1, iC0M1

  ! Local index variables for optimized Y-vector tracking slices
  integer :: jF0M1, jFLM1, jF0M2, jFLM2, jF0M3, jFLM3
  integer :: jF0P1, jFLP1, jF0P2, jFLP2, jF0P3, jFLP3
 
  jf0 = iMesh%JfSTART
  jfL = iMesh%JfEND 
  
  iC0   = iMesh%IcSTART
  iC0M1 = iC0 - 1
 
  lpos   = NX - NG 
  lposP1 = lpos + 1

  ! Index boundary assignments
  jF0M1 = jF0 - 1;   jFLM1 = jFL - 1
  jF0M2 = jF0 - 2;   jFLM2 = jFL - 2
  jF0M3 = jF0 - 3;   jFLM3 = jFL - 3
  jF0P1 = jF0 + 1;   jFLP1 = jFL + 1
  jF0P2 = jF0 + 2;   jFLP2 = jFL + 2
  jF0P3 = jF0 + 3;   jFLP3 = jFL + 3

  SELECT CASE ( SptScheme )   
     CASE( 3004 )
        DFDY(jF0:jFL, iC0M1:lposP1) = &
           COEF%COEF_SDY_B2_L(1) * Up(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_M    * Up(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(1) * Up(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(1) * Uyz(jF0  :jFL  , iC0M1:lposP1)

     CASE( 3006 )
        DFDY(jF0:jFL, iC0M1:lposP1) = &
           COEF%COEF_SDY_B2_L(1) * Up(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_M    * Up(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(1) * Up(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(1) * Uyz(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1)

     CASE( 3008 )
        DFDY(jF0:jFL, iC0M1:lposP1) = &
           COEF%COEF_SDY_B2_L(2) * Up(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_L(1) * Up(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_M    * Up(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(1) * Up(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(2) * Up(jF0P2:jFLP2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(1) * Uyz(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1)

     CASE( 3010 )
        DFDY(jF0:jFL, iC0M1:lposP1) = &
           COEF%COEF_SDY_B2_L(2) * Up(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_L(1) * Up(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_M    * Up(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(1) * Up(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(2) * Up(jF0P2:jFLP2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(3) * Uyz(jF0M3:jFLM3, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(1) * Uyz(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(3) * Uyz(jF0P2:jFLP2, iC0M1:lposP1)

     CASE( 3012 )
        DFDY(jF0:jFL, iC0M1:lposP1) = &
           COEF%COEF_SDY_B2_L(3) * Up(jF0M3:jFLM3, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_L(2) * Up(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_L(1) * Up(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_M    * Up(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(1) * Up(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(2) * Up(jF0P2:jFLP2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B2_R(3) * Up(jF0P3:jFLP3, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(3) * Uyz(jF0M3:jFLM3, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(2) * Uyz(jF0M2:jFLM2, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_L(1) * Uyz(jF0M1:jFLM1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(1) * Uyz(jF0  :jFL  , iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(2) * Uyz(jF0P1:jFLP1, iC0M1:lposP1) + &
           COEF%COEF_SDY_B1_R(3) * Uyz(jF0P2:jFLP2, iC0M1:lposP1)

     CASE DEFAULT
        ERROR STOP 'Unknown Spatial Scheme: Compute_SecondDerivativeY_XZ'
  END SELECT 

END subroutine Compute_SecondDerivativeY_XZ

    
   
!MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  
!MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  
!MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS  MISC FUNCTIONS   
   


   elemental function get_erf_cellvalue(a, b, rho, c0, c1) result(value)
      real(wp), intent(in) :: a, b, rho, c0, c1
      real(wp) :: value
      real(wp) :: a_minus, a_plus, b_minus, b_plus
      real(wp) :: erf_group, exp_group, denominator
      real(wp), parameter :: pi = acos(-1.0_wp)

      a_minus = rho * (a - c0)
      a_plus  = rho * (a + c1)
      b_minus = rho * (b - c0)
      b_plus  = rho * (b + c1)

      erf_group = sqrt(pi) * ( &
           b_minus * erf(b_minus) &
         - b_plus  * erf(b_plus)  &
         - a_minus * erf(a_minus) &
         + a_plus  * erf(a_plus) )

      exp_group = &
         - exp(-a_minus*a_minus) &
         + exp(-b_minus*b_minus) &
         + exp(-a_plus*a_plus)   &
         - exp(-b_plus*b_plus)

      denominator = 2.0_wp* rho * sqrt(pi) * (a - b)
      value = (erf_group + exp_group) / denominator
   end function get_erf_cellvalue

 

subroutine save_paraview_2d(filename, Mesh, F)
    implicit none
    character(len=*), intent(in) :: filename
    type(MeshDomain), intent(in) :: Mesh
    type(scalar_spatial_data), intent(in) :: F

    integer :: iu, i, j, ios
    integer :: npx_active, npy_active, ncellx_active, ncelly_active
    
    ! We only want to export the PHYSICAL domain, excluding ghost cells.
    ! Active Points count
    npx_active = Mesh%IfEnd - Mesh%IfStart + 1
    npy_active = Mesh%JfEnd - Mesh%JfStart + 1
    
    ! Active Cells count
    ncellx_active = Mesh%IcEnd - Mesh%IcStart + 1
    ncelly_active = Mesh%JcEnd - Mesh%JcStart + 1

    open(newunit=iu, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file for ParaView export: ", trim(filename)
        return
    end if

    ! 1. Write VTK Header
    write(iu, '(A)') '# vtk DataFile Version 3.0'
    write(iu, '(A)') 'MMM 2D Spatial Data Export'
    write(iu, '(A)') 'ASCII'
    write(iu, '(A)') 'DATASET STRUCTURED_GRID'
    
    ! VTK expects dimensions as (X, Y, Z). Since it's 2D, Z dimension is 1.
    write(iu, '(A, 3I8)') 'DIMENSIONS ', npx_active, npy_active, 1

    ! 2. Write Grid Points (Nodal Coordinates)
    write(iu, '(A, I8, A)') 'POINTS ', npx_active * npy_active, ' double'
    do j = Mesh%JfStart, Mesh%JfEnd
        do i = Mesh%IfStart, Mesh%IfEnd
            ! VTK requires 3D coordinates, so we set Z component to 0.0
            write(iu, '(3E16.8)') Mesh%Xp(i), Mesh%Yp(j), 0.0_wp
        end do
    end do

    ! 3. Write POINT DATA (Node-centered values)
    write(iu, '(A, I8)') 'POINT_DATA ', npx_active * npy_active
    write(iu, '(A)') 'SCALARS Phi_Point double 1'
    write(iu, '(A)') 'LOOKUP_TABLE default'
    do j = Mesh%JfStart, Mesh%JfEnd
        do i = Mesh%IfStart, Mesh%IfEnd
            write(iu, '(E16.8)') F%p(j, i)
        end do
    end do

    ! 4. Write CELL DATA (Cell-centered values)
    write(iu, '(A, I8)') 'CELL_DATA ', ncellx_active * ncelly_active
    write(iu, '(A)') 'SCALARS Phi_Cell double 1'
    write(iu, '(A)') 'LOOKUP_TABLE default'
    do j = Mesh%JcStart, Mesh%JcEnd
        do i = Mesh%IcStart, Mesh%IcEnd
            write(iu, '(E16.8)') F%xyz(j, i)
        end do
    end do

    close(iu)
    print *, "Successfully saved ParaView file: ", trim(filename)
end subroutine save_paraview_2d


subroutine save_paraview_cell_2d(filename, Grid, cell_data_array, var_name)
    implicit none
    character(len=*), intent(in) :: filename
    type(MeshDomain), intent(in) :: Grid
    real(wp), dimension(:,:), intent(in) :: cell_data_array ! Dimensions should cover active cells
    character(len=*), intent(in) :: var_name

    integer :: i, j, out_unit, ios
    integer :: nx_points, ny_points
    real(wp) :: x_origin, y_origin

    ! VTK Dimensions describe the POINT grid topology.
    ! For N cells, we have N + 1 grid points bounding them.
    nx_points = Grid%NcellX + 1
    ny_points = Grid%NcellY + 1

    ! Calculate global physical origin from grid structure
    ! (Assuming IcSTART/JcSTART index references your active bounding domain)
    x_origin = Grid%Xp(Grid%IfSTART)
    y_origin = Grid%Yp(Grid%JfSTART)

    open(newunit=out_unit, file=trim(filename), status='unknown', action='write', iostat=ios)
    if (ios /= 0) then
        print *, "Error: Could not open file for ParaView output: ", trim(filename)
        return
    end if

    ! 1. Write the mandatory VTK standard header
    write(out_unit, '(A)') "# vtk DataFile Version 3.0"
    write(out_unit, '(A)') "2D Cell-Averaged Simulation Data Output"
    write(out_unit, '(A)') "ASCII"
    write(out_unit, '(A)') "DATASET STRUCTURED_POINTS"

    ! 2. Define geometry and layout topology
    ! For 2D datasets, Z dimension is always 1 point, thickness spacing 0.0
    write(out_unit, '(A, 3I8)') "DIMENSIONS ", nx_points, ny_points, 1
    write(out_unit, '(A, 3E16.7)') "ORIGIN     ", x_origin, y_origin, 0.0_wp
    write(out_unit, '(A, 3E16.7)') "SPACING    ", Grid%DX, Grid%DY, 0.0_wp

    ! 3. Write out the CELL_DATA segment
    ! Total cells = NcellX * NcellY
    write(out_unit, '(A, I10)') "CELL_DATA ", Grid%NcellX * Grid%NcellY
    write(out_unit, '(A, A)') "SCALARS ", trim(var_name) // " double 1"
    write(out_unit, '(A)') "LOOKUP_TABLE default"

    ! VTK orders cells nested loops: X varies fastest, then Y, then Z.
    do j = Grid%JcSTART, Grid%JcEND
        do i = Grid%IcSTART, Grid%IcEND
            write(out_unit, '(E17.9)') cell_data_array(i, j)
        end do
    end do

    close(out_unit)
    print *, "Successfully wrote cell data to: ", trim(filename)
end subroutine save_paraview_cell_2d



 !******************************************************
  




 !******************************************************
 !******************************************************
subroutine save_paraview_point_binary_2d(filename, Mesh, F, var_name)
    implicit none
    character(len=*), intent(in) :: filename
    type(MeshDomain), intent(in) :: Mesh
    type(scalar_spatial_data), intent(in) :: F
    character(len=*), intent(in), optional :: var_name   ! <-- Added optional argument
    
    character(len=256) :: char_buffer    
    character(len=64)  :: actual_var_name                 ! <-- Local variable to store the name
    integer(8) :: offset_points    
    integer(8) :: block_size    

    integer :: iu, i, j, ios
    integer :: npx, npy
    
    ! Determine the variable name to use
    if (present(var_name)) then
        actual_var_name = trim(adjustl(var_name))
    else
        actual_var_name = "Phi_Point"
    end if

    ! Extract active physical domain bounds for points
    npx = Mesh%IfEnd - Mesh%IfStart + 1
    npy = Mesh%JfEnd - Mesh%JfStart + 1

    ! Open file using standard stream access for continuous raw writing
    open(newunit=iu, file=trim(filename), status='replace', &
         access='stream', form='unformatted', iostat=ios)
         
    if (ios /= 0) then
        print *, "Error opening file for binary export: ", trim(filename)
        return
    end if

    ! 1. XML Header and Grid Structure Definition
    call write_str(iu, '<?xml version="1.0"?>'//new_line('a'))
    call write_str(iu, '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//new_line('a'))
    
    ! Extent defines array indices: MinX MaxX MinY MaxY MinZ MaxZ (0-based)
    write(char_buffer, '("<StructuredGrid WholeExtent=""0 ",I0," 0 ",I0," 0 0"">")') npx-1, npy-1
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    write(char_buffer, '("<Piece Extent=""0 ",I0," 0 ",I0," 0 0"">")') npx-1, npy-1
    call write_str(iu, trim(char_buffer)//new_line('a'))

    ! 2. Point Data (Nodal Fields) Descriptor using the dynamic name
    write(char_buffer, '("<PointData Scalars=""",A,""">")') trim(actual_var_name)
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    write(char_buffer, '("<DataArray type=""Float64"" Name=""",A,""" format=""appended"" offset=""0""/>")') trim(actual_var_name)
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    call write_str(iu, '</PointData>'//new_line('a'))

    ! 3. Grid Geometry Points Descriptor
    offset_points = 8 + (int(npx, 8) * int(npy, 8) * 8)
    
    call write_str(iu, '<Points>'//new_line('a'))
    write(char_buffer, '("<DataArray type=""Float64"" Name=""Points"" NumberOfComponents=""3"" format=""appended"" offset=""",I0,"""/>")') offset_points
    call write_str(iu, trim(char_buffer)//new_line('a'))
    call write_str(iu, '</Points>'//new_line('a'))

    ! Close out XML geometry structural blocks
    call write_str(iu, '</Piece>'//new_line('a'))
    call write_str(iu, '</StructuredGrid>'//new_line('a'))

    ! 4. APPENDED BINARY DATA SECTION
    call write_str(iu, '<AppendedData encoding="raw">'//new_line('a')//'_')

    ! Block A: Write Point Values
    block_size = int(npx, 8) * int(npy, 8) * 8_8
    write(iu) block_size
    do j = Mesh%JfStart, Mesh%JfEnd
        do i = Mesh%IfStart, Mesh%IfEnd
            write(iu) real(F%p(j, i), wp)
        end do
    end do

    ! Block B: Write Coordinates (Points X, Y, Z)
    block_size = int(npx, 8) * int(npy, 8) * 3_8 * 8_8
    write(iu) block_size
    do j = Mesh%JfStart, Mesh%JfEnd
        do i = Mesh%IfStart, Mesh%IfEnd
            write(iu) real(Mesh%Xp(i), wp)
            write(iu) real(Mesh%Yp(j), wp)
            write(iu) 0.0_wp
        end do
    end do

    ! Close out file tags
    call write_str(iu, new_line('a')//'</AppendedData>'//new_line('a'))
    call write_str(iu, '</VTKFile>'//new_line('a'))
    close(iu)

contains

    subroutine write_str(unit_idx, str)
        integer, intent(in) :: unit_idx
        character(len=*), intent(in) :: str
        integer :: k
        do k = 1, len(str)
            write(unit_idx) str(k:k)
        end do
    end subroutine write_str

end subroutine save_paraview_point_binary_2d





 !******************************************************
 !******************************************************
subroutine save_paraview_point_binary_2d_w_GC(filename, Mesh, F, var_name) !VTS with ghost cells
    implicit none
    character(len=*), intent(in) :: filename
    type(MeshDomain), intent(in) :: Mesh
    type(scalar_spatial_data), intent(in) :: F
    character(len=*), intent(in), optional :: var_name   ! <-- Added optional argument
    
    character(len=256) :: char_buffer    
    character(len=64)  :: actual_var_name                 ! <-- Local variable to store the name
    integer(8) :: offset_points    
    integer(8) :: block_size    

    integer :: iu, i, j, ios
    integer :: npx, npy
    
    ! Determine the variable name to use
    if (present(var_name)) then
        actual_var_name = trim(adjustl(var_name))
    else
        actual_var_name = "Phi_Point"
    end if

    ! Extract active physical domain bounds for points
    npx =   Mesh%NPXtot 
    npy =   Mesh%NPYtot 

    ! Open file using standard stream access for continuous raw writing
    open(newunit=iu, file=trim(filename), status='replace', &
         access='stream', form='unformatted', iostat=ios)
         
    if (ios /= 0) then
        print *, "Error opening file for binary export: ", trim(filename)
        return
    end if

    ! 1. XML Header and Grid Structure Definition
    call write_str(iu, '<?xml version="1.0"?>'//new_line('a'))
    call write_str(iu, '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//new_line('a'))
    
    ! Extent defines array indices: MinX MaxX MinY MaxY MinZ MaxZ (0-based)
    write(char_buffer, '("<StructuredGrid WholeExtent=""0 ",I0," 0 ",I0," 0 0"">")') npx-1, npy-1
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    write(char_buffer, '("<Piece Extent=""0 ",I0," 0 ",I0," 0 0"">")') npx-1, npy-1
    call write_str(iu, trim(char_buffer)//new_line('a'))

    ! 2. Point Data (Nodal Fields) Descriptor using the dynamic name
    write(char_buffer, '("<PointData Scalars=""",A,""">")') trim(actual_var_name)
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    write(char_buffer, '("<DataArray type=""Float64"" Name=""",A,""" format=""appended"" offset=""0""/>")') trim(actual_var_name)
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    call write_str(iu, '</PointData>'//new_line('a'))

    ! 3. Grid Geometry Points Descriptor
    offset_points = 8 + (int(npx, 8) * int(npy, 8) * 8)
    
    call write_str(iu, '<Points>'//new_line('a'))
    write(char_buffer, '("<DataArray type=""Float64"" Name=""Points"" NumberOfComponents=""3"" format=""appended"" offset=""",I0,"""/>")') offset_points
    call write_str(iu, trim(char_buffer)//new_line('a'))
    call write_str(iu, '</Points>'//new_line('a'))

    ! Close out XML geometry structural blocks
    call write_str(iu, '</Piece>'//new_line('a'))
    call write_str(iu, '</StructuredGrid>'//new_line('a'))

    ! 4. APPENDED BINARY DATA SECTION
    call write_str(iu, '<AppendedData encoding="raw">'//new_line('a')//'_')

    ! Block A: Write Point Values
    block_size = int(npx, 8) * int(npy, 8) * 8_8
    write(iu) block_size
    do j = 1, npy
        do i = 1, npx
            write(iu) real(F%p(j, i), wp)
        end do
    end do

    ! Block B: Write Coordinates (Points X, Y, Z)
    block_size = int(npx, 8) * int(npy, 8) * 3_8 * 8_8
    write(iu) block_size
    do j = 1, npy
        do i = 1, npx
            write(iu) real(Mesh%Xp(i), wp)
            write(iu) real(Mesh%Yp(j), wp)
            write(iu) 0.0_wp
        end do
    end do

    ! Close out file tags
    call write_str(iu, new_line('a')//'</AppendedData>'//new_line('a'))
    call write_str(iu, '</VTKFile>'//new_line('a'))
    close(iu)

contains

    subroutine write_str(unit_idx, str)
        integer, intent(in) :: unit_idx
        character(len=*), intent(in) :: str
        integer :: k
        do k = 1, len(str)
            write(unit_idx) str(k:k)
        end do
    end subroutine write_str

end subroutine save_paraview_point_binary_2d_w_GC






subroutine save_paraview_binary_2d_cell_point(filename, Mesh, F)
    implicit none
    character(len=*), intent(in) :: filename
   character(len=256) :: char_buffer    
    type(MeshDomain), intent(in) :: Mesh
    type(scalar_spatial_data), intent(in) :: F
    integer(8) :: offset_cell, offset_points    
    integer(8) :: block_size    

    integer :: iu, i, j, ios
    integer :: npx, npy, ncellx, ncelly
    
    ! Extract active physical domain bounds
    npx = Mesh%IfEnd - Mesh%IfStart + 1
    npy = Mesh%JfEnd - Mesh%JfStart + 1
    ncellx = Mesh%IcEnd - Mesh%IcStart + 1
    ncelly = Mesh%JcEnd - Mesh%JcStart + 1

    ! Open file using standard stream access for continuous raw writing
    open(newunit=iu, file=trim(filename), status='replace', &
         access='stream', form='unformatted', iostat=ios)
         
    if (ios /= 0) then
        print *, "Error opening file for binary export: ", trim(filename)
        return
    end if

    ! 1. XML XML Header and Grid Structure Definition
    ! Note: 'Appended' data allows us to write XML descriptors first, followed by raw binary blocks.
    call write_str(iu, '<?xml version="1.0"?>'//new_line('a'))
    call write_str(iu, '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//new_line('a'))
    
    ! Extent defines array indices: MinX MaxX MinY MaxY MinZ MaxZ (0-based)
    write(char_buffer, '("<StructuredGrid WholeExtent=""0 ",I0," 0 ",I0," 0 0"">")') npx-1, npy-1
    call write_str(iu, trim(char_buffer)//new_line('a'))
    
    write(char_buffer, '("<Piece Extent=""0 ",I0," 0 ",I0," 0 0"">")') npx-1, npy-1
    call write_str(iu, trim(char_buffer)//new_line('a'))

    ! 2. Point Data (Nodal Fields) Descriptor
    call write_str(iu, '<PointData Scalars="Phi_Point">'//new_line('a'))
    call write_str(iu, '<DataArray type="Float64" Name="Phi_Point" format="appended" offset="0"/>'//new_line('a'))
    call write_str(iu, '</PointData>'//new_line('a'))

    ! 3. Cell Data (Element Fields) Descriptor
    ! Offset tracks where the binary data block begins relative to the <AppendedData> section start tag.
    ! Size of Phi_Point block = 8 bytes (Header) + (npx * npy * 8 bytes)

    offset_cell = 8 + (int(npx, 8) * int(npy, 8) * 8)
    
    call write_str(iu, '<CellData Scalars="Phi_Cell">'//new_line('a'))
    write(char_buffer, '("<DataArray type=""Float64"" Name=""Phi_Cell"" format=""appended"" offset=""",I0,"""/>")') offset_cell
    call write_str(iu, trim(char_buffer)//new_line('a'))
    call write_str(iu, '</CellData>'//new_line('a'))

    ! 4. Grid Geometry Points Descriptor
    ! Size of Phi_Cell block = ncellx * ncelly * 8 bytes
    offset_points = offset_cell + 8 + (int(ncellx, 8) * int(ncelly, 8) * 8)
    
    call write_str(iu, '<Points>'//new_line('a'))
    write(char_buffer, '("<DataArray type=""Float64"" Name=""Points"" NumberOfComponents=""3"" format=""appended"" offset=""",I0,"""/>")') offset_points
    call write_str(iu, trim(char_buffer)//new_line('a'))
    call write_str(iu, '</Points>'//new_line('a'))

    ! Close out XML geometry structural blocks
    call write_str(iu, '</Piece>'//new_line('a'))
    call write_str(iu, '</StructuredGrid>'//new_line('a'))

    ! 5. APPENDED BINARY DATA SECTION
    ! Every separate block requires an initial 8-byte integer header denoting the data size in bytes.
    call write_str(iu, '<AppendedData encoding="raw">'//new_line('a')//'_')

    ! Block A: Write Point Values (Phi_Point)

    block_size = int(npx, 8) * int(npy, 8) * 8_8
    write(iu) block_size
    do j = Mesh%JfStart, Mesh%JfEnd
        do i = Mesh%IfStart, Mesh%IfEnd
            write(iu) real(F%p(j, i), wp)
        end do
    end do

    ! Block B: Write Cell Values (Phi_Cell)
    block_size = int(ncellx, 8) * int(ncelly, 8) * 8_8
    write(iu) block_size
    do j = Mesh%JcStart, Mesh%JcEnd
        do i = Mesh%IcStart, Mesh%IcEnd
            write(iu) real(F%xyz(j, i), wp)
        end do
    end do

    ! Block C: Write Coordinates (Points X, Y, Z)
    block_size = int(npx, 8) * int(npy, 8) * 3_8 * 8_8
    write(iu) block_size
    do j = Mesh%JfStart, Mesh%JfEnd
        do i = Mesh%IfStart, Mesh%IfEnd
            write(iu) real(Mesh%Xp(i), wp)
            write(iu) real(Mesh%Yp(j), wp)
            write(iu) 0.0_wp
        end do
    end do

    ! Close out file tags
    call write_str(iu, new_line('a')//'</AppendedData>'//new_line('a'))
    call write_str(iu, '</VTKFile>'//new_line('a'))
    close(iu)

    print *, "Successfully saved raw binary VTS file: ", trim(filename)

contains

    subroutine write_str(unit_idx, str)
        integer, intent(in) :: unit_idx
        character(len=*), intent(in) :: str
        integer :: k
        do k = 1, len(str)
            write(unit_idx) str(k:k)
        end do
    end subroutine write_str

   end subroutine save_paraview_binary_2d_cell_point

   ! =========================================================================
   ! INTERMEDIATE SPATIAL DRIVER
   ! Aggregates Advection and Diffusion calculations into dFpdx and NetFLUXdFpdx
   ! =========================================================================


! =========================================================================
! 1D ADAPTIVE QUADRATURE ALONG X (Y is a fixed scalar constant)
! =========================================================================
! =========================================================================
! 1D ADAPTIVE QUADRATURE ALONG X (Y is a fixed scalar constant)
! =========================================================================
recursive function AdaptiveNewtonCote1D_X(f, y_const, x1, x2, tol) result(integral)
    implicit none
    interface
         function f(x, y) result(res)
            import :: wp
            real(wp), intent(in) :: x, y
            real(wp) :: res
        end function f
    end interface
    real(wp), intent(in) :: y_const, x1, x2, tol
    real(wp) :: integral
    
    real(wp) :: h, mid, coarse_int, fine_int, left_int, right_int
    real(wp), dimension(7) :: x_pts, f_pts
    integer :: i
    
    mid = 0.5_wp * (x1 + x2)
    
    ! --- Coarse Mesh evaluation (1 panel, 7 points) ---
    h = (x2 - x1) / 6.0_wp
    do i = 1, 7
        x_pts(i) = x1 + real(i-1, wp) * h
        f_pts(i) = f(x_pts(i), y_const)
    end do
    coarse_int = (3.0_wp * h / 10.0_wp) * (f_pts(1) + 5.0_wp*f_pts(2) + f_pts(3) + &
                                           6.0_wp*f_pts(4) + f_pts(5) + 5.0_wp*f_pts(6) + f_pts(7))
    
    ! --- Fine Mesh refinement (split into 2 panels) ---
    left_int  = Weddle1D_X_Panel(x1, mid)
    right_int = Weddle1D_X_Panel(mid, x2)
    fine_int  = left_int + right_int
    
    ! Error checking and adaptive recursion
    if (abs(fine_int - coarse_int) < tol) then
        integral = fine_int
    else
        integral = AdaptiveNewtonCote1D_X(f, y_const, x1, mid, 0.5_wp*tol) + &
                   AdaptiveNewtonCote1D_X(f, y_const, mid, x2, 0.5_wp*tol)
    end if
contains
    function Weddle1D_X_Panel(xa, xb) result(res)
        real(wp), intent(in) :: xa, xb
        real(wp) :: res, hh
        hh = (xb - xa) / 6.0_wp
        ! Directly evaluates 'f' from parent scope
        res = f(xa, y_const) + 5.0_wp*f(xa+hh, y_const) + f(xa+2.0_wp*hh, y_const) + &
              6.0_wp*f(xa+3.0_wp*hh, y_const) + f(xa+4.0_wp*hh, y_const) + 5.0_wp*f(xa+5.0_wp*hh, y_const) + f(xb, y_const)
        res = res * (3.0_wp * hh / 10.0_wp)
    end function Weddle1D_X_Panel
end function AdaptiveNewtonCote1D_X

! =========================================================================
! 1D ADAPTIVE QUADRATURE ALONG Y (X is a fixed scalar constant)
! =========================================================================
recursive function AdaptiveNewtonCote1D_Y(f, x_const, y1, y2, tol) result(integral)
    implicit none
    interface
         function f(x, y) result(res)
!!!!            use, intrinsic :: iso_fortran_env, only: wp
            import :: wp
            real(wp), intent(in) :: x, y
            real(wp) :: res
        end function f
    end interface
    real(wp), intent(in) :: x_const, y1, y2, tol
    real(wp) :: integral
    
    real(wp) :: h, mid, coarse_int, fine_int, left_int, right_int
    real(wp), dimension(7) :: y_pts, f_pts
    integer :: i
    
    mid = 0.5_wp * (y1 + y2)
    
    h = (y2 - y1) / 6.0_wp
    do i = 1, 7
        y_pts(i) = y1 + real(i-1, wp) * h
        f_pts(i) = f(x_const, y_pts(i))
    end do
    coarse_int = (3.0_wp * h / 10.0_wp) * (f_pts(1) + 5.0_wp*f_pts(2) + f_pts(3) + &
                                           6.0_wp*f_pts(4) + f_pts(5) + 5.0_wp*f_pts(6) + f_pts(7))
    
    left_int  = Weddle1D_Y_Panel(y1, mid)
    right_int = Weddle1D_Y_Panel(mid, y2)
    fine_int  = left_int + right_int
    
    if (abs(fine_int - coarse_int) < tol) then
        integral = fine_int
    else
        integral = AdaptiveNewtonCote1D_Y(f, x_const, y1, mid, 0.5_wp*tol) + &
                   AdaptiveNewtonCote1D_Y(f, x_const, mid, y2, 0.5_wp*tol)
    end if
contains
    function Weddle1D_Y_Panel(ya, yb) result(res)
        real(wp), intent(in) :: ya, yb
        real(wp) :: res, hh
        hh = (yb - ya) / 6.0_wp
        ! Directly evaluates 'f' from parent scope
        res = f(x_const, ya) + 5.0_wp*f(x_const, ya+hh) + f(x_const, ya+2.0_wp*hh) + &
              6.0_wp*f(x_const, ya+3.0_wp*hh) + f(x_const, ya+4.0_wp*hh) + 5.0_wp*f(x_const, ya+5.0_wp*hh) + f(x_const, yb)
        res = res * (3.0_wp * hh / 10.0_wp)
    end function Weddle1D_Y_Panel
end function AdaptiveNewtonCote1D_Y

! =========================================================================
! 2D ADAPTIVE QUADRATURE (Performs nested adaptive 1D calls)
! =========================================================================
! =========================================================================
! TRUE 2D ADAPTIVE QUADRATURE (7x7 Point Patches via Quad-Tree Subdivision)
! =========================================================================


recursive function AdaptiveSimpson2D(f, x1, x2, y1, y2, tol) result(integral)
    implicit none
    interface
         function f(x, y) result(res)
!!!            use, intrinsic :: iso_fortran_env, only: wp
            import :: wp
            real(wp), intent(in) :: x, y
            real(wp) :: res
        end function f
    end interface
    real(wp), intent(in) :: x1, x2, y1, y2, tol
    real(wp) :: integral
    
    real(wp) :: mid_x, mid_y, coarse_int, fine_int
    real(wp) :: q1, q2, q3, q4
    
    mid_x = 0.5_wp * (x1 + x2)
    mid_y = 0.5_wp * (y1 + y2)
    
    ! --- Coarse Mesh Evaluation: 1 single 2D Patch (7x7 = 49 points) ---
    coarse_int = Weddle2D_Patch(x1, x2, y1, y2)
    integral = coarse_int
     

contains

    ! Helper function to evaluate a full 7x7 Newton-Cotes / Weddle grid over an area
    function Weddle2D_Patch(xa, xb, ya, yb) result(res)
        real(wp), intent(in) :: xa, xb, ya, yb
        real(wp) :: res, wsum, locx(3), locy(3), wij
        real(wp) :: hx, hy, x_val, y_val, f_val
        real(wp), dimension(3) :: w
        integer :: i, j
        
        ! Weddle 1D base weights
        w = [1.0_wp, 4.0_wp, 1.0_wp]
        locx = [xa, 0.5_wp*(xa+xb), xb ]
        locy = [ya, 0.5_wp*(ya+yb), yb ]        
        
        hx = (xb - xa)  
        hy = (yb - ya) 
        
        res = 0.0_wp
        wsum = 0.0_wp
        
        ! Perform the 7x7 double summation
        do j = 1, 3
            do i = 1, 3

                !36
                ! Evaluate the parent function 'f' via host association
                wij = (w(i) * w(j))
                f_val = f( locx(i), locy(j) )
                wsum = wsum + wij
                
                ! Add to accumulator: (w_i * w_j) * f(x,y)
                res = res + wij * f_val
            end do
        end do
        
        ! Apply the 2D scaling factor: (3*hx/10) * (3*hy/10)
        res = res * ( hx * hy / wsum )
    end function Weddle2D_Patch

end function AdaptiveSimpson2D












recursive function AdaptiveNewtonCote2D(f, x1, x2, y1, y2, tol) result(integral)
    implicit none
    interface
         function f(x, y) result(res)
!!!            use, intrinsic :: iso_fortran_env, only: wp
            import :: wp
            real(wp), intent(in) :: x, y
            real(wp) :: res
        end function f
    end interface
    real(wp), intent(in) :: x1, x2, y1, y2, tol
    real(wp) :: integral
    
    real(wp) :: mid_x, mid_y, coarse_int, fine_int
    real(wp) :: q1, q2, q3, q4
    
    mid_x = 0.5_wp * (x1 + x2)
    mid_y = 0.5_wp * (y1 + y2)
    
    ! --- Coarse Mesh Evaluation: 1 single 2D Patch (7x7 = 49 points) ---
    coarse_int = Weddle2D_Patch(x1, x2, y1, y2)
    
    ! --- Fine Mesh Evaluation: 4 Sub-Patches (Quadrants) ---
    q1 = Weddle2D_Patch(x1,    mid_x, y1,    mid_y) ! Bottom-Left
    q2 = Weddle2D_Patch(mid_x, x2,    y1,    mid_y) ! Bottom-Right
    q3 = Weddle2D_Patch(x1,    mid_x, mid_y, y2)    ! Top-Left
    q4 = Weddle2D_Patch(mid_x, x2,    mid_y, y2)    ! Top-Right
    fine_int = q1 + q2 + q3 + q4
    
    ! Error checking and 2D quad-tree adaptive recursion
    if (abs(fine_int - coarse_int) < tol) then
        integral = fine_int
    else
        ! Subdivide tolerance by 4 since we split into 4 pieces
        integral = AdaptiveNewtonCote2D(f, x1,    mid_x, y1,    mid_y, 0.25_wp*tol) + &
                   AdaptiveNewtonCote2D(f, mid_x, x2,    y1,    mid_y, 0.25_wp*tol) + &
                   AdaptiveNewtonCote2D(f, x1,    mid_x, mid_y, y2,    0.25_wp*tol) + &
                   AdaptiveNewtonCote2D(f, mid_x, x2,    mid_y, y2,    0.25_wp*tol)
    end if

contains

    ! Helper function to evaluate a full 7x7 Newton-Cotes / Weddle grid over an area
    function Weddle2D_Patch(xa, xb, ya, yb) result(res)
        real(wp), intent(in) :: xa, xb, ya, yb
        real(wp) :: res
        real(wp) :: hx, hy, x_val, y_val, f_val
        real(wp), dimension(7) :: w
        integer :: i, j
        
        ! Weddle 1D base weights
        w = [1.0_wp, 5.0_wp, 1.0_wp, 6.0_wp, 1.0_wp, 5.0_wp, 1.0_wp]
        
        hx = (xb - xa) / 6.0_wp
        hy = (yb - ya) / 6.0_wp
        
        res = 0.0_wp
        
        ! Perform the 7x7 double summation
        do j = 1, 7
            y_val = ya + real(j-1, wp) * hy
            do i = 1, 7
                x_val = xa + real(i-1, wp) * hx
                
                ! Evaluate the parent function 'f' via host association
                f_val = f(x_val, y_val)
                
                ! Add to accumulator: (w_i * w_j) * f(x,y)
                res = res + (w(i) * w(j)) * f_val
            end do
        end do
        
        ! Apply the 2D scaling factor: (3*hx/10) * (3*hy/10)
        res = res * (9.0_wp * hx * hy / 100.0_wp)
    end function Weddle2D_Patch

end function AdaptiveNewtonCote2D




    !-----------------------------------------------------------
    elemental function get_R(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        res = sqrt(x**2 + y**2) + ep
    end function get_R

    !-----------------------------------------------------------
    elemental function get_VT(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        real(wp) :: r_val
        
        r_val = get_R(x, y)
        ! cosh(x) is standard in Fortran, sech(x) = 1.0 / cosh(x)
        res = CONST_FACTOR * (1.0_wp / (cosh(r_val)+ep) )**2 * tanh(r_val) 
    end function get_VT

    !-----------------------------------------------------------
    elemental function get_Omeg(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        res = get_VT(x, y) / get_R(x, y)
    end function get_Omeg

    !-----------------------------------------------------------
    elemental function get_myU(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        res = -y*get_VT(x, y) / get_R(x, y)
    end function get_myU

    !-----------------------------------------------------------
    elemental function get_myV(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        res = x * get_VT(x, y)  / get_R(x, y)
    end function get_myV
    
    
!-----------------------------------------------------------
    ! New: Exact time-dependent analytical solution phi(x,y,t)
    !-----------------------------------------------------------
    elemental function get_analytical_phi(x, y, t) result(res)
        real(wp), intent(in) :: x, y, t
        real(wp) :: res
        real(wp) :: omega 
        
        omega = get_Omeg(x, y)
        res = -tanh( 0.5_wp * (y * cos(omega * t) - x * sin(omega * t)) )
    end function get_analytical_phi   
    
 
    
    
    ! New: Approximate time-dependent analytical boundary condition for phi(x,y,t)
    !-----------------------------------------------------------
    elemental function get_analytical_phi_BC(x, y, t) result(res)
        use tanh_spline_mod    , only: tanh_spline    
        real(wp), intent(in) :: x, y, t
        real(wp) :: res
        real(wp) :: omega, xi, d_xi
        
		real(wp), parameter :: C0 = -0.4621171572600098_wp
		real(wp), parameter :: C1 = -0.7864477329659274_wp
		real(wp), parameter :: C2 =  0.3634277259169493_wp
		real(wp), parameter :: C3 = -0.0941328404222044_wp
		real(wp), parameter :: C4 =  0.1306141312217036_wp         
        
        omega = get_Omeg(x, y)        
		! 1. Compute the cheap exact trigonometric terms
!		xi = 0.5_wp * ( cos(omega * y * t) - sin(omega * x * t) )
!		d_xi = xi - 0.5_wp

        res = -tanh_spline( 0.5_wp * (y * cos(omega * t) - x * sin(omega * t)) )

		! 2. Fast polynomial evaluation via Horner's scheme (Zero tanh calls!)
!		res = C0 + d_xi * (C1 + d_xi * (C2 + d_xi * (C3 + C4 * d_xi)))        


    end function get_analytical_phi_BC        
     
    elemental function tanhf_wp(x) result(res)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    
    ! Define 'wp' (Working Precision) as double precision
    integer, parameter :: wp = real64
    
    ! Input and output types
    real(kind=wp), intent(in) :: x
    real(kind=wp)             :: res
    
    ! Local variables for the Taylor series terms
    real(kind=wp) :: t1, t2, t3, t4, t5, t6
    
    ! Snap to 1.0 (with sign) when polynomial begins to deviate at the tails
    if (abs(x) > 1.365_wp) then
        res = sign(1.0_wp, x)
        return
    end if
    
    ! Polynomial approximation
    t1 = x
    t2 = x**3 * (1.0_wp / 3.0_wp)
    t3 = x**5 * (2.0_wp / 15.0_wp)
    t4 = x**7 * (17.0_wp / 315.0_wp)
    t5 = x**9 * (62.0_wp / 2835.0_wp)
    t6 = x**11 * (1382.0_wp / 155925.0_wp)
    
    res = t1 - t2 + t3 - t4 + t5 - t6
    
end function tanhf_wp

 
    elemental function phi_time_wrapped(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        res = get_analytical_phi(x, y, current_bc_time)
    end function phi_time_wrapped
    
    elemental function phi_time_wrapped2(x, y) result(res)
        real(wp), intent(in) :: x, y
        real(wp) :: res
        res = get_analytical_phi_BC(x, y, current_bc_time)
    end function phi_time_wrapped2    
   
end module MMM_module_spatial

