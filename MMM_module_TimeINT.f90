MODULE MMM_module_TimeINT
   USE, INTRINSIC :: iso_fortran_env, ONLY: real64
   USE MMM_module_spatial 
   USE module_SchemeSetup, ONLY: wp   
   use matlab_io_mod, only:save_bin_matlab   
   use module_physics
   
   !public ::set_fronto_timedep_BC
   IMPLICIT NONE
!   public :: TimeNOW
   REAL(wp) ::  TRK=0.0_wp,  TRHS=0.0_wp, TtimeINT =0.0_wp,  TBC =0.0_wp
      
      TYPE :: TimeIntArray
      TYPE(scalar_spatial_data), ALLOCATABLE :: K(:)
      TYPE(scalar_spatial_data) :: VAR0 ! For initial data (TimeN)
      TYPE(scalar_spatial_data) :: VAR1, VAR2   ! Temporary derivative pointers
      TYPE(scalar_spatial_data) :: VAR3, VAR4   ! Temporary second derivative pointers      
      


      INTEGER :: Nstage
      INTEGER :: Method
      REAL (WP) :: DT, TimeNOW      
      
      !----------------------Butcher Table
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: a
      REAL(wp), DIMENSION(:),   ALLOCATABLE :: b, c

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: F0xyz, F0xz, F0yz, F0p      
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: F1xyz, F1xz, F1yz, F1p 
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: F2xyz, F2xz, F2yz, F2p  
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: F3xyz, F3xz, F3yz, F3p  
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: F4xyz, F4xz, F4yz, F4p  
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: F5xyz, F5xz, F5yz, F5p   
            
   CONTAINS
      PROCEDURE :: init => ALLOCATE_TimeIntArray    
        
      PROCEDURE :: free => deALLOCATE_TimeIntArray      
   END TYPE TimeIntArray

CONTAINS

! =====================================================================
! CLEAN TIME INTEGRATION DRIVER
! =====================================================================
! =====================================================================
! OPTIMIZED TIME INTEGRATION DRIVER
! =====================================================================
SUBROUTINE TimeIntegration( iMesh, Velocity, DUDX, DUDY, DVDX, DVDY, F, TMP, COEF, TimeNOW )
   USE module_SchemeSetup, ONLY: isMM
   USE MMM_module_spatial, ONLY: Coefficient_DIFferentiation
   
   IMPLICIT NONE
   TYPE(MeshDomain), INTENT(in) :: iMesh
   TYPE(velocity_spatial_data), TARGET, INTENT(in) :: Velocity
   TYPE(scalar_spatial_data), INTENT(inout) :: F                      ! Scalar Phi (F^n -> F^{n+1})
   TYPE(scalar_spatial_data), INTENT(in) :: DUDX, DUDY, DVDX, DVDY    ! Velocity Gradients
   TYPE(TimeIntArray), TARGET, INTENT(inout) :: TMP                  ! Time Integration Arrays
   TYPE(Coefficient_DIFferentiation), INTENT(in) :: COEF 
    

   ! Local Pointers as Clean Aliases
   TYPE(scalar_spatial_data), POINTER :: TimeN => NULL()
   real(wp) ::  DT, TRK0, TRK1,  TRHS0, TRHS1,  TtimeINT0, TtimeINT1, TBC0, TBC1
   real(wp), intent(in) ::  TimeNOW
   INTEGER :: iST, s, i


   call CPU_TIME( TtimeINT0 )
   ! Set pointer to our initial step buffer
   TimeN => TMP%VAR0

   ! 1. Save the starting state of this timestep: TimeN = F^n
   !    Using whole-array assignment for optimal performance
   
   !CALL Copy_SCALAR_AtoB( F , TimeN )
   TimeN = F
 
   
   DT = TMP%DT
   
   call CPU_TIME( TtimeINT1 )   
   TRK = TRK  + (TtimeINT1 - TtimeINT0)

   ! 2. Loop Through Runge-Kutta substep ------------------------------------------

!CALL set_fronto_timedep_BC( F, TimeNOW  )
   DO iST = 1, TMP%Nstage
!!  print *,  "iST:Maxval Fp" , ist, maxval(  F%p )      
      !set time-dependent boundary coditions
      call CPU_TIME( TBC0 )
        CALL set_fronto_timedep_BC_simp( F, TimeNOW  + TMP%c(iST) *DT )
      call CPU_TIME( TBC1 )
      TBC = TBC + TBC1 - TBC0


      ! Compute Spatial Physics RHS on the current intermediate state F
      call CPU_TIME( TRHS0 )
        CALL Compute_RHS(iMesh, Velocity, DUDX, DUDY, DVDX, DVDY, F, TMP, COEF, TMP%K(iST), TimeNOW )
      call CPU_TIME( TRHS1 )
      
      TRHS = TRHS + TRHS1 - TRHS0
      
      
     call CPU_TIME( TRK0 )

	     CALL ScaleField_scalar( TMP%K(iST), DT )
	     	 
	      IF (iST < TMP%Nstage) THEN
		 CALL Compute_RK_NextSubstep( F, TimeN, TMP%a(1:iST, ist), TMP%K(1:iST), iST  )
	      END IF
      
      call CPU_TIME( TRK1 )
      TRK = TRK + TRK1 - TRK0  
 
   END DO
 
  ! 3. Final Step Update: Accumulate results directly to get F^{n+1}
   call CPU_TIME( TRK0 )
        CALL Accumulate_Final_Step( F, TimeN, TMP%b, TMP%K )
   call CPU_TIME( TRK1 )
   TRK = TRK + TRK1 - TRK0
   !    Accumulate final step results

!!  print *,  "iST:Maxval Fp" , ist, maxval(  F%p )    

   ! Safe cleanup of local pointers
   TimeN => NULL()
   
   call CPU_TIME( TtimeINT1 )
   TtimeINT = TtimeINT + TtimeINT1 - TtimeINT0

END SUBROUTINE TimeIntegration


! =====================================================================
! SPECIALIZED STAGE ACCUMULATOR 
! =====================================================================
SUBROUTINE Compute_RK_NextSubstep( F, TimeN, coeffs, K_stages, iST )
   USE module_SchemeSetup, ONLY: isMM, TimeINTmethod
 
   IMPLICIT NONE
   TYPE(scalar_spatial_data), INTENT(inout) :: F
   TYPE(scalar_spatial_data), INTENT(in)    :: TimeN
   REAL(wp), DIMENSION(:), INTENT(in)       :: coeffs    ! Vector of RK coefficients
   TYPE(scalar_spatial_data), DIMENSION(:), INTENT(in)  :: K_stages  ! Array of K stages
   CHARACTER(len=32)  :: fname 
 !!  REAL(wp), intent(in) :: DT
   
   INTEGER :: s, N_stages, iST

   N_stages = iST !!SIZE(coeffs)
 
   IF( TimeINTmethod == 41 ) THEN
 
		    F%p = TimeN%p + coeffs( iST ) * K_stages( iST )%p 

		    IF (isMM) THEN
		        F%xyz  = TimeN%xyz + coeffs( N_stages ) *  K_stages( N_stages )%xyz  
		         F%yz  =  TimeN%yz + coeffs( N_stages ) *  K_stages( N_stages )%yz  
		         F%xz  =  TimeN%xz + coeffs( N_stages ) *  K_stages( N_stages )%xz  
		    END IF 
		     

   ELSEIF( TimeINTmethod == 51 ) THEN

   
      SELECT CASE ( N_stages )
      CASE ( 1 )
         ! Row 1: coeffs = [ 1/4, 0, 0, 0, 0 ]
         F%p = F%p + (1.0_wp/4.0_wp) * K_stages(1)%p
         IF (isMM) THEN
            F%xyz = TimeN%xyz + (1.0_wp/4.0_wp) * K_stages(1)%xyz
            F%xz  = TimeN%xz  + (1.0_wp/4.0_wp) * K_stages(1)%xz
            F%yz  = TimeN%yz  + (1.0_wp/4.0_wp) * K_stages(1)%yz
         END IF

      CASE ( 2 )
         ! Row 2: coeffs = [ 1/8, 1/8, 0, 0, 0 ]
         F%p = F%p + (1.0_wp/8.0_wp) * K_stages(2)%p
         IF (isMM) THEN
            F%xyz = TimeN%xyz + (1.0_wp/8.0_wp) * K_stages(1)%xyz + (1.0_wp/8.0_wp) * K_stages(2)%xyz
            F%xz  = TimeN%xz  + (1.0_wp/8.0_wp) * K_stages(1)%xz  + (1.0_wp/8.0_wp) * K_stages(2)%xz
            F%yz  = TimeN%yz  + (1.0_wp/8.0_wp) * K_stages(1)%yz  + (1.0_wp/8.0_wp) * K_stages(2)%yz
         END IF

      CASE ( 3 )
         ! Row 3: coeffs = [ 0, -1/2, 1, 0, 0 ] -> K(1) is exactly skipped!
         F%p = F%p + 1.0_wp * K_stages(3)%p
         IF (isMM) THEN
            F%xyz = TimeN%xyz + (-0.5_wp * K_stages(2)%xyz + 1.0_wp * K_stages(3)%xyz)
            F%xz  = TimeN%xz  + (-0.5_wp * K_stages(2)%xz  + 1.0_wp * K_stages(3)%xz)
            F%yz  = TimeN%yz  + (-0.5_wp * K_stages(2)%yz  + 1.0_wp * K_stages(3)%yz)
         END IF

      CASE ( 4 )
         ! Row 4: coeffs = [ 3/16, 0, 0, 9/16, 0 ] -> K(2) & K(3) are exactly skipped!
         F%p = F%p + (9.0_wp/16.0_wp) * K_stages(4)%p
         IF (isMM) THEN
            F%xyz = TimeN%xyz + ((3.0_wp/16.0_wp) * K_stages(1)%xyz + (9.0_wp/16.0_wp) * K_stages(4)%xyz)
            F%xz  = TimeN%xz  + ((3.0_wp/16.0_wp) * K_stages(1)%xz  + (9.0_wp/16.0_wp) * K_stages(4)%xz)
            F%yz  = TimeN%yz  + ((3.0_wp/16.0_wp) * K_stages(1)%yz  + (9.0_wp/16.0_wp) * K_stages(4)%yz)
         END IF

      CASE ( 5 )
         ! Row 5: coeffs = [ -3/7, 2/7, 12/7, -12/7, 8/7 ] (All non-zero)
         F%p = F%p + (8.0_wp/7.0_wp) * K_stages(5)%p
         IF (isMM) THEN
            F%xyz = TimeN%xyz +  (&
                            - (3.0_wp/7.0_wp)  * K_stages(1)%xyz + (2.0_wp/7.0_wp) * K_stages(2)%xyz &
                            + (12.0_wp/7.0_wp) * K_stages(3)%xyz - (12.0_wp/7.0_wp) * K_stages(4)%xyz &
                            + (8.0_wp/7.0_wp)  * K_stages(5)%xyz)
            F%xz  = TimeN%xz   + (- (3.0_wp/7.0_wp) * K_stages(1)%xz  + (2.0_wp/7.0_wp) * K_stages(2)%xz  &
                            + (12.0_wp/7.0_wp) * K_stages(3)%xz  - (12.0_wp/7.0_wp) * K_stages(4)%xz  &
                            + (8.0_wp/7.0_wp) * K_stages(5)%xz)
            F%yz  = TimeN%yz   + (- (3.0_wp/7.0_wp) * K_stages(1)%yz  + (2.0_wp/7.0_wp) * K_stages(2)%yz  &
                            + (12.0_wp/7.0_wp) * K_stages(3)%yz  - (12.0_wp/7.0_wp) * K_stages(4)%yz  &
                            + (8.0_wp/7.0_wp) * K_stages(5)%yz)
         END IF

      END SELECT   
   

   ELSE 
   
   ! NEED TO USE F itself as the accumulator
		F%p = TimeN%p
		IF (isMM) THEN
			F%xyz = TimeN%xyz
			F%xz  = TimeN%xz
			F%yz  = TimeN%yz
		END IF   
   
	   DO s = 1, N_stages
		  IF (ABS(coeffs(s)) > 2.0_wp*epsilon(1.0_wp)) THEN
		     ! p accumulates incrementally
		     F%p = F%p +  coeffs(s) * K_stages(s)%p
		     
		     ! xyz, xz, yz build directly from the constant initial state (TimeN)
		     IF (isMM) THEN
		        F%xyz = F%xyz +   coeffs(s) * K_stages(s)%xyz
		        F%xz  = F%xz  +   coeffs(s) * K_stages(s)%xz
		        F%yz  = F%yz  +   coeffs(s) * K_stages(s)%yz
		     END IF
		  END IF
	   END DO
	ENDIF
	 
END SUBROUTINE Compute_RK_NextSubstep

! =====================================================================
! OPTIMIZED FINAL STEP ACCUMULATOR 
! ===================================================================== 

 
SUBROUTINE Accumulate_Final_Step( F, TimeN, coeffs, K_stages )
   USE module_SchemeSetup, ONLY: isMM, TimeINTmethod
   IMPLICIT NONE
   TYPE(scalar_spatial_data), INTENT(inout) :: F
   TYPE(scalar_spatial_data), INTENT(in)    :: TimeN     ! The initial state F0
   REAL(wp), DIMENSION(:), INTENT(in)       :: coeffs    ! Vector of RK final weights (b)
   TYPE(scalar_spatial_data), DIMENSION(:), INTENT(in)  :: K_stages  ! Array of K stages
   
   INTEGER :: s, N_stages

   N_stages = SIZE(coeffs)

   ! Select execution branch based on Time Integration Method
   IF ( TimeINTmethod == 41 ) THEN
      ! ------------------------------------------------------------------
      ! BRANCH 41: CLASSIC RK4 (Unrolled, completely loop-free, direct FMA)
      ! ------------------------------------------------------------------
      F%p = TimeN%p + coeffs(1) * K_stages(1)%p &
                    + coeffs(2) * K_stages(2)%p &
                    + coeffs(3) * K_stages(3)%p &
                    + coeffs(4) * K_stages(4)%p

      IF (isMM) THEN
         F%xyz = TimeN%xyz + coeffs(1) * K_stages(1)%xyz &
                           + coeffs(2) * K_stages(2)%xyz &
                           + coeffs(3) * K_stages(3)%xyz &
                           + coeffs(4) * K_stages(4)%xyz

         F%xz  = TimeN%xz  + coeffs(1) * K_stages(1)%xz &
                           + coeffs(2) * K_stages(2)%xz &
                           + coeffs(3) * K_stages(3)%xz &
                           + coeffs(4) * K_stages(4)%xz

         F%yz  = TimeN%yz  + coeffs(1) * K_stages(1)%yz &
                           + coeffs(2) * K_stages(2)%yz &
                           + coeffs(3) * K_stages(3)%yz &
                           + coeffs(4) * K_stages(4)%yz
      END IF

   ELSEIF ( TimeINTmethod == 51 ) THEN
      ! ------------------------------------------------------------------
      ! BRANCH 51: LAWSON RK5 (Unrolled, direct FMA, skips Stage 2)
      ! ------------------------------------------------------------------
      F%p = TimeN%p + coeffs(1) * K_stages(1)%p &
                    + coeffs(3) * K_stages(3)%p &
                    + coeffs(4) * K_stages(4)%p &
                    + coeffs(5) * K_stages(5)%p &
                    + coeffs(6) * K_stages(6)%p

      IF (isMM) THEN
         F%xyz = TimeN%xyz + coeffs(1) * K_stages(1)%xyz &
                           + coeffs(3) * K_stages(3)%xyz &
                           + coeffs(4) * K_stages(4)%xyz &
                           + coeffs(5) * K_stages(5)%xyz &
                           + coeffs(6) * K_stages(6)%xyz

         F%xz  = TimeN%xz  + coeffs(1) * K_stages(1)%xz &
                           + coeffs(3) * K_stages(3)%xz &
                           + coeffs(4) * K_stages(4)%xz &
                           + coeffs(5) * K_stages(5)%xz &
                           + coeffs(6) * K_stages(6)%xz

         F%yz  = TimeN%yz  + coeffs(1) * K_stages(1)%yz &
                           + coeffs(3) * K_stages(3)%yz &
                           + coeffs(4) * K_stages(4)%yz &
                           + coeffs(5) * K_stages(5)%yz &
                           + coeffs(6) * K_stages(6)%yz
      END IF

   ELSE
      ! ------------------------------------------------------------------
      ! GENERAL BRANCH: High-Order (RK7 - RK10) Sparse-Aware Loop
      ! ------------------------------------------------------------------
      ! 1. Initialize F back to the starting step state (F0)
      F%p = TimeN%p
      IF (isMM) THEN
         F%xyz = TimeN%xyz
         F%xz  = TimeN%xz
         F%yz  = TimeN%yz
      END IF
 
      ! 2. Loop over stages
      DO s = 1, N_stages
      IF (ABS(coeffs(s)) > 2.0_wp * EPSILON(1.0_wp)) THEN
      
            F%p = F%p + coeffs(s) * K_stages(s)%p
            IF (isMM) THEN
               F%xyz = F%xyz + coeffs(s) * K_stages(s)%xyz
               F%xz  = F%xz  + coeffs(s) * K_stages(s)%xz
               F%yz  = F%yz  + coeffs(s) * K_stages(s)%yz
            END IF
         END IF
      END DO

   END IF
END SUBROUTINE Accumulate_Final_Step 
 
 
 
 
 
 
 
 
 
 
! =====================================================================
! ENCAPSULATED SPATIAL RHS CALCULATOR (The heavy physics loops)
! =====================================================================
SUBROUTINE Compute_RHS(iMesh, Velocity, DUDX, DUDY, DVDX, DVDY, F, TMP, COEF, K_stage, TimeNow)
   USE module_SchemeSetup, ONLY: isMM, METHNC
   USE MMM_module_spatial, ONLY: Coefficient_DIFferentiation
 

   IMPLICIT NONE
   TYPE(MeshDomain), INTENT(in) :: iMesh
   TYPE(velocity_spatial_data), TARGET, INTENT(in) :: Velocity
   TYPE(scalar_spatial_data), INTENT(in) :: DUDX, DUDY, DVDX, DVDY
   TYPE(scalar_spatial_data), INTENT(in) :: F
   TYPE(TimeIntArray), TARGET, INTENT(inout) :: TMP
   TYPE(Coefficient_DIFferentiation), INTENT(in) :: COEF  
   TYPE(scalar_spatial_data), INTENT(inout) :: K_stage ! K stage to compute

   ! Local Pointers to avoid deep-scoping overhead
   TYPE(scalar_spatial_data), POINTER :: U => NULL()
   TYPE(scalar_spatial_data), POINTER :: V => NULL()
   TYPE(scalar_spatial_data), POINTER :: DFDX => NULL(), SFSX => NULL()
   TYPE(scalar_spatial_data), POINTER :: DFDY => NULL(), SFSY => NULL()
   
   real (wp) t0,t1,t2, TimeNow, t3,t4

   INTEGER :: NcellXtot, NcellYtot, NPxtot, NPytot
   INTEGER :: iF0, iFL, jF0, jFL, iC0, iCL, jC0, jCL
   INTEGER :: iC0P1, iCLP1, jC0P1, jCLP1
   INTEGER :: iF0e, iFLe, jF0e, jFLe, iC0e, iCLe, jC0e, jCLe, jF0P1, jFLP1
   INTEGER :: iF0eP1, iFLeP1, jF0eP1, jFLeP1, iC0eP1, iCLeP1, jC0eP1, jCLeP1 
   INTEGER :: iF0eM1, iFLeM1, jF0eM1, jFLeM1, iC0eM1, iCLeM1, jC0eM1, jCLeM1   
   INTEGER :: i,j
   REAL(wp) :: iDX, iDY, DT, iDXDT, iDYDT, tlin, tnon

   REAL(wp), PARAMETER :: FAC4O = 1.0_wp / 12.0_wp    
   REAL(wp), PARAMETER :: FAC6O_Q22 = 36.0_wp  / 720.0_wp    
   REAL(wp), PARAMETER :: FAC6O_Q13 = 3.0_wp / 720.0_wp    
   
	! --- YZ Direction Registers ---
	real(wp) :: u_p0, u_p1, f_p0_y, f_p1_y
	real(wp) :: u_yz_m1, u_yz_0, u_yz_p1
	real(wp) :: f_yz_m1, f_yz_0, f_yz_p1
	real(wp) :: u_diff, f_diff_y, u_dev, f_dev_y
	real(wp) :: corr_yz

	! --- XZ Direction Registers ---
	real(wp) :: v_p0, v_p1, f_p0_x, f_p1_x
	real(wp) :: v_xz_m1, v_xz_0, v_xz_p1
	real(wp) :: f_xz_m1, f_xz_0, f_xz_p1
	real(wp) :: v_diff, f_diff_x, v_dev, f_dev_x
	real(wp) :: corr_xz, UF_curr, VF_curr, DFDY_curr, DFDX_curr, Fp_curr, Up_curr, Vp_curr
 
	REAL(wp) :: Fp_jp1, Fp_ip1
	REAL(wp) :: Vp_jp1, Up_ip1	, KiDX, KiDY

   
   DT = TMP%DT
   iDX = 1.0_wp / iMesh%DX
   iDY = 1.0_wp / iMesh%DY
   iDXDT = iDX*DT
   iDYDT = iDY*DT
   
   KiDX = KDIFF*iDX
   KiDY = KDIFF*iDY   
   ! 1. Index Configurations

   NcellXtot = iMesh%NcellXtot; 
   NcellYtot = iMesh%NcellYtot  
   NPxtot = iMesh%NPxtot; 
   NPytot = iMesh%NPytot

   iF0 = iMesh%IfSTART; iFL = iMesh%IfEND 
   iC0 = iMesh%IcSTART; iCL = iMesh%IcEND
   jF0 = iMesh%JfSTART; jFL = iMesh%JfEND 
   jC0 = iMesh%JcSTART; jCL = iMesh%JcEND

   iC0P1 = iC0 + 1; iCLP1 = iCL + 1	
   jC0P1 = jC0 + 1; jCLP1 = jCL + 1

   iF0e = iF0 - 1; iFLe = iFL + 1; iC0e = iC0 - 1; iCLe = iCL + 1
   jF0e = jF0 - 1; jFLe = jFL + 1; jC0e = jC0 - 1; jCLe = jCL + 1    

   iF0eP1 = iF0e + 1; iFLeP1 = iFLe + 1; iC0eP1 = iC0e + 1; iCLeP1 = iCLe + 1
   iF0eM1 = iF0e - 1; iFLeM1 = iFLe - 1; iC0eM1 = iC0e - 1; iCLeM1 = iCLe - 1
   jF0eP1 = jF0e + 1; jFLeP1 = jFLe + 1; jC0eP1 = jC0e + 1; jCLeP1 = jCLe + 1
   jF0eM1 = jF0e - 1; jFLeM1 = jFLe - 1; jC0eM1 = jC0e - 1; jCLeM1 = jCLe - 1    
   
   jF0P1 = jF0 + 1
   jFLP1 = jFL + 1

   ! Link Local temporary pointers
   DFDX => TMP%VAR1; DFDY => TMP%VAR2  
   SFSX => TMP%VAR3; SFSY => TMP%VAR4 
   U => Velocity%U;  V => Velocity%V 

   ! 2. Branching Spatial Physics (FDM vs Multi-Moment)
   IF (.NOT. isMM) THEN
!        call cpu_time(t1)
      ! Pure FDM CODE   Pure FDM CODE   Pure FDM CODE   Pure FDM CODE   Pure FDM CODE   Pure FDM CODE   
      CALL Compute_FirstDerivativeX_FD(DFDX%p, F%p, iMesh, COEF, NPxtot, NPytot)   
      CALL Compute_FirstDerivativeY_FD(DFDY%p, F%p, iMesh, COEF, NPxtot, NPytot)     	

!!      IF( KDIFF .GT. 1E-15_WP)      THEN
        CALL Compute_SecondDerivativeX_FD(SFSX%p, F%p, iMesh, COEF, NPxtot, NPytot)   
        CALL Compute_SecondDerivativeY_FD(SFSY%p, F%p, iMesh, COEF, NPxtot, NPytot)     	      
!!      ENDIF   



IF(1==1) THEN
	DO i = iF0, iFL
	    DO j = jF0, jFL
	!! The multi-moment method offers multiple opportunities for data reuse.
	!! Explicitly storing frequently reused values in registers can improve performance.
	!! To eliminate any concern of implementation bias, this FDM version uses the same
	!! coding style as the finite-volume implementation.
	!! In practice, however, this implementation and the cleaner vector-slice version (see further below)
	!! exhibit nearly identical performance.

		! Current-cell values
		Fp_curr   = F%p(j,i)
		Up_curr   = U%p(j,i)
		Vp_curr   = V%p(j,i)
		DFDX_curr = DFDX%p(j,i)
		DFDY_curr = DFDY%p(j,i)

	 

		! Primary RHS
		K_stage%p(j,i) = -( &
		      Up_curr *  DFDX_curr &
		    + Fp_curr * (DUDX%p(j,i) + DVDY%p(j,i)) &
		    + Vp_curr *                DFDY_curr ) &
		    + KDIFF * (SFSX%p(j,i) + SFSY%p(j,i))

	 

	    END DO
	END DO
ELSE
!! NICE AND EASY TO READ CODE
      DO i = iF0, iFL    
         K_stage%p(jF0:jFL, i) =  (&
            - U%p(jF0:jFL, i) *   DFDX%p(jF0:jFL, i) &
            - F%p(jF0:jFL, i) * ( DUDX%p(jF0:jFL, i)+DVDY%p(jF0:jFL, i) ) &
            - V%p(jF0:jFL, i) *                      DFDY%p(jF0:jFL, i) &
            )  &             
            + KDIFF*(   SFSX%p(jF0:jFL, i) +  SFSY%p(jF0:jFL, i) ) 
     ENDDO          
ENDIF           


! ALTERNATIVELY, ONE CANE USE ARRAY-SLICE APPROACH !!!!!!!!!!!!!!!!!!
! THIS IS HIGHLY COMPILER DEPEDENT. A VERY FOOD COMPILER CAN FULLY VECTORIZE THE OPRATION
! BUT A BAD ONE WILL CREATE TEMP ARRAY AND DELIVER BAD PERFORMANCE
!         K_stage%p(jF0:jFL, iF0:iFL) =  (&
!            - U%p(jF0:jFL, iF0:iFL) * DFDX%p(jF0:jFL, iF0:iFL) &
!            - F%p(jF0:jFL, iF0:iFL) * DUDX%p(jF0:jFL, iF0:iFL) &
!            - V%p(jF0:jFL, iF0:iFL) * DFDY%p(jF0:jFL, iF0:iFL) &
!            - F%p(jF0:jFL, iF0:iFL) * DVDY%p(jF0:jFL, iF0:iFL) &
!            ) + KDIFF*(   SFSX%p(jF0:jFL,  iF0:iFL) +  SFSY%p(jF0:jFL,  iF0:iFL) ) 
 
   ELSE
   
 
! Multi-Moment CODE -------------------------------------------------------------------------
! Multi-Moment CODE -------------------------------------------------------------------------
      
! --- P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment 
! --- P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment 
! --- P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment P Moment 

!!!!!!!!!!!! FOR CONVECTIVE TERM      
      CALL Compute_FirstDerivativeX_YZ(DFDX%p, F%xz, F%p, iMesh, COEF, NPxtot, NPytot) 
      CALL Compute_FirstDerivativeY_XZ(DFDY%p, F%yz, F%p, iMesh, COEF, NPxtot, NPytot)	 

!!!!!!!!!!!! FOR DIFFUSION TERM
      CALL Compute_SecondDerivativeX_YZ(SFSX%p, F%xz, F%p, iMesh, COEF, NPxtot, NPytot) 
      CALL Compute_SecondDerivativeY_XZ(SFSY%p, F%yz, F%p, iMesh, COEF, NPxtot, NPytot)
 
      
	 
	 
	!! The multi-moment method offers multiple opportunities for data reuse.
	!! Explicitly storing frequently reused values in registers can improve performa
	DO i = iF0, iFL
	    DO j = jF0, jFL

		! Current-cell values
		Fp_curr   = F%p(j,i)
		Up_curr   = U%p(j,i)
		Vp_curr   = V%p(j,i)
		DFDX_curr = DFDX%p(j,i)
		DFDY_curr = DFDY%p(j,i)

		! Shared current-cell fluxes
		UF_curr = Up_curr * Fp_curr
		VF_curr = Vp_curr * Fp_curr

		! Neighbor values
		Fp_jp1 = F%p(j+1,i)
		Vp_jp1 = V%p(j+1,i)

		Fp_ip1 = F%p(j,i+1)
		Up_ip1 = U%p(j,i+1)

		! FDM OF Fp                             ***K_stage%p is initialized here ***
		K_stage%p(j,i) = -( &
		      Up_curr *  DFDX_curr &
		    + Fp_curr * (DUDX%p(j,i) + DVDY%p(j,i)) &
		    + Vp_curr *                DFDY_curr ) &
		    + KDIFF * (SFSX%p(j,i) + SFSY%p(j,i))

		! FVM-PART OF Fyz (Y-flux divergence)  ***K_stage%yz is initialized here ***
		K_stage%yz(j,i) = ( &
		    -(Vp_jp1 * Fp_jp1 - VF_curr) &
		    + KDIFF * (DFDY%p(j+1,i) - DFDY_curr) &
		    ) * iDY

		! FVM-PART OF Fxz (X-flux divergence)  ***K_stage%xz is initialized here ***
		K_stage%xz(j,i) = ( &
		    -(Up_ip1 * Fp_ip1 - UF_curr) &
		    + KDIFF * (DFDX%p(j,i+1) - DFDX_curr) &
		    ) * iDX

	    END DO
	END DO


  
! --- YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  
! --- YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment        
! --- YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment  YZ Moment 
   
      CALL Compute_FirstDerivativeX_YZ(  DFDX%yz, F%xyz, F%yz, iMesh, COEF, NPxtot, NPytot - 1) 
      CALL Compute_SecondDerivativeX_YZ( SFSX%yz, F%xyz, F%yz, iMesh, COEF, NPxtot, NPytot - 1)       

     
	!ACCUMMULATE THE RHS--------------------------------  
	DO i = iF0, iFL
	DO j = jC0, jCL

	DFDX_CURR = DFDX%yz(j, i)

	!FDM PART OF Fyz 
	K_stage%yz(j, i) = K_stage%yz(j, i) - ( &
	      F%yz(j, i) * DUDX%yz(j, i) &
	    + U%yz(j, i) * DFDX_CURR   ) &
	    + KDIFF * SFSX%yz(j, i)

	!DIFFUSION OF PART OF Fxyz     ***K_stage%xyz is initialized here***
	K_stage%xyz(j, i) = KiDX*(  DFDX%yz(j, i+1) - DFDX_CURR )

	END DO
	END DO

 
     
             
   

       !ACCUMMULATE THE RHS from NONLINEAR CORRECTION for Fyz-------------------------------
      IF (METHNC >= 4) THEN   

         CALL NONLIN_CORR_YZ_FDM_DIRECTION(K_stage%yz, TMP%F1yz, &
            F%p, F%yz, DFDX%p, DFDX%yz, &
            U%p, U%yz, DUDX%p, DUDX%yz )
      END IF



! --- XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  ---
! --- XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  ---
! --- XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  XZ Moment  ---
       

      CALL Compute_FirstDerivativeY_XZ (DFDY%xz, F%xyz, F%xz, iMesh, COEF, NPxtot - 1, NPytot)      
      CALL Compute_SecondDerivativeY_XZ(SFSY%xz, F%xyz, F%xz, iMesh, COEF, NPxtot - 1, NPytot)       
      

      

     !ACCUMMULATE THE RHS--------------------------------  
   DO i = iC0, iCL
    DO j = jF0, jFL
         DFDY_CURR =   DFDY%xz(j, i)
        
        !FDM PART OF Fxz 
        K_stage%xz(j, i) = K_stage%xz(j, i) - ( &
              F%xz(j, i) * DVDY%xz(j, i) &
            + V%xz(j, i) * DFDY_CURR   ) &
            + KDIFF * SFSY%xz(j, i)
        
        !DIFFUSION OF PART OF Fxyz     
        K_stage%xyz(j, i) = K_stage%xyz(j, i) + KiDY*( DFDY%xz(j+1, i) - DFDY_CURR ) 

    END DO
    END DO
 
    
 

  
      !ACCUMMULATE THE RHS from NONLINEAR Correction for Fxz -------------------------------             
      IF (METHNC >= 4) THEN   
 
         CALL NONLIN_CORR_XZ_FDM_DIRECTION(K_stage%xz, TMP%F1xz, &
            F%p, F%xz, DFDY%p, DFDY%xz, &
            V%p, V%xz, DVDY%p, DVDY%xz )
      ENDIF
 

! FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) 
! FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) 
! FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) FVM MOMENT (XYZ-MOMENT) 

   ! NO NONLINEAR CORRECTION
   IF (METHNC .LT.  4) THEN 
  
       do i = iF0, iFL
	       TMP%F0yz(jC0:jCL,i)  = U%yz(jC0:jCL,i)  * F%yz(jC0:jCL,i)      
	       TMP%F0xz(jF0:jFL,i)  = V%xz(jF0:jFL,i)  * F%xz(jF0:jFL,i) 
       end do
       
 
			do i = iC0, iCL
			   do j = jC0, jCL
			      K_stage%xyz(j, i) = K_stage%xyz(j, i) + &
				 - ( TMP%F0yz(j, i+1) - TMP%F0yz(j, i) ) * iDX &
				 - ( TMP%F0xz(j+1, i) - TMP%F0xz(j, i) ) * iDY
			end do
			end do       

   ELSEIF (METHNC == 4) THEN 
   
   
   
  
		 do i = iF0, iFL
		   do j = jF0, jFL

		      !=========================================================================
		      ! 1. Compute & Accumulate Nonlinear Flux Correction: YZ Direction
		      !=========================================================================
		      corr_yz = 0.0_wp

		     !! if (METHNC >= 4) then
			 u_p0   = U%p(j, i)
			 u_p1   = U%p(j+1, i)
			 f_p0_y = F%p(j, i)
			 f_p1_y = F%p(j+1, i)

			 u_diff   = u_p1 - u_p0
			 f_diff_y = f_p1_y - f_p0_y
		         corr_yz = FAC4O * ( f_diff_y * u_diff )
	 
			    u_yz_0  = U%yz(j, i)
			    f_yz_0  = F%yz(j, i)

			    ! 6th-order Q22 + Q31 correction terms
			    TMP%F0yz(j, i) =  f_yz_0*u_yz_0  + corr_yz

		      !=========================================================================
		      ! 2. Compute & Accumulate Nonlinear Flux Correction: XZ Direction
		      !=========================================================================
		      corr_xz = 0.0_wp


			 v_p0   = V%p(j, i)
			 v_p1   = V%p(j, i+1)
			 f_p0_x = F%p(j, i)
			 f_p1_x = F%p(j, i+1)

			 v_diff   = v_p1 - v_p0
			 f_diff_x = f_p1_x - f_p0_x


		         corr_xz = FAC4O * ( f_diff_x * v_diff )

                        TMP%F0xz(j, i)  =  V%xz(j, i) * F%xz(j, i) +  corr_xz
	 

		   end do
		end do

		!===============================================================================
		! 3. Final Divergence Calculation using Corrected Fluxes
		!===============================================================================
		do i = iC0, iCL
		   do j = jC0, jCL
		      K_stage%xyz(j, i) = K_stage%xyz(j, i) + &
			 - ( TMP%F0yz(j, i+1) - TMP%F0yz(j, i) ) * iDX &
			 - ( TMP%F0xz(j+1, i) - TMP%F0xz(j, i) ) * iDY
		   end do
		end do   

!------------NC6 NC6 NC6 NC6 NC6 NC6 NC6 NC6 
    ELSE 
!------------NC6 NC6 NC6 NC6 NC6 NC6 NC6 NC6 
	do i = iF0, iFL
	   do j = jF0, jFL

	      !=========================================================================
	      ! 1. Compute & Accumulate Nonlinear Flux Correction: YZ Direction
	      !=========================================================================
	      corr_yz = 0.0_wp

	     !! if (METHNC >= 4) then
		 u_p0   = U%p(j, i)
		 u_p1   = U%p(j+1, i)
		 f_p0_y = F%p(j, i)
		 f_p1_y = F%p(j+1, i)

		 u_diff   = u_p1 - u_p0
		 f_diff_y = f_p1_y - f_p0_y

		 ! 4th-order base correction term
		 corr_yz = FAC4O * ( f_diff_y * u_diff )

		!! if (METHNC >= 6) then
		    u_yz_m1 = U%yz(j-1, i)
		    u_yz_0  = U%yz(j, i)
		    u_yz_p1 = U%yz(j+1, i)

		    f_yz_m1 = F%yz(j-1, i)
		    f_yz_0  = F%yz(j, i)
		    f_yz_p1 = F%yz(j+1, i)

		    u_dev   = u_p1 + u_p0 - 2.0_wp * u_yz_0
		    f_dev_y = f_p1_y + f_p0_y - 2.0_wp * f_yz_0

		    ! 6th-order Q22 + Q31 correction terms
		    TMP%F0yz(j, i) =  f_yz_0*u_yz_0 +  corr_yz &
		       + FAC6O_Q22 * ( u_dev * f_dev_y ) &
		       - FAC6O_Q13 * ( &
		           ( u_yz_p1 - u_yz_m1 - 2.0_wp * u_diff ) * f_diff_y + &
		           ( f_yz_p1 - f_yz_m1 - 2.0_wp * f_diff_y ) * u_diff   &
		         )

	      !=========================================================================
	      ! 2. Compute & Accumulate Nonlinear Flux Correction: XZ Direction
	      !=========================================================================
	      corr_xz = 0.0_wp

!!	      if (METHNC >= 4) then
		 v_p0   = V%p(j, i)
		 v_p1   = V%p(j, i+1)
		 f_p0_x = F%p(j, i)
		 f_p1_x = F%p(j, i+1)

		 v_diff   = v_p1 - v_p0
		 f_diff_x = f_p1_x - f_p0_x

		 ! 4th-order base correction term
		 corr_xz = FAC4O * ( f_diff_x * v_diff )

!!		 if (METHNC >= 6) then
		    v_xz_m1 = V%xz(j, i-1)
		    v_xz_0  = V%xz(j, i)
		    v_xz_p1 = V%xz(j, i+1)

		    f_xz_m1 = F%xz(j, i-1)
		    f_xz_0  = F%xz(j, i)
		    f_xz_p1 = F%xz(j, i+1)

		    v_dev   = v_p1 + v_p0 - 2.0_wp * v_xz_0
		    f_dev_x = f_p1_x + f_p0_x - 2.0_wp * f_xz_0

		    ! 6th-order Q22 + Q31 correction terms
		    TMP%F0xz(j, i)  = v_xz_0*f_xz_0+corr_xz &
		       + FAC6O_Q22 * ( v_dev * f_dev_x ) &
		       - FAC6O_Q13 * ( &
		           ( v_xz_p1 - v_xz_m1 - 2.0_wp * v_diff ) * f_diff_x + &
		           ( f_xz_p1 - f_xz_m1 - 2.0_wp * f_diff_x ) * v_diff   &
		         )

	   end do
	end do

	!===============================================================================
	! 3. Final Divergence Calculation using Corrected Fluxes
	!===============================================================================
	do i = iC0, iCL
	   do j = jC0, jCL
	      K_stage%xyz(j, i) = K_stage%xyz(j, i) &
		 - ( TMP%F0yz(j, i+1) - TMP%F0yz(j, i)  ) * iDX &
		 - ( TMP%F0xz(j+1, i) - TMP%F0xz(j, i)  ) * iDY
	   end do
	end do

    ENDIF
! END END FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ 
! END END FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ FVM-XYZ 

!!!!!!!!!!!!!!!! END IF OF THE METHOD SELECTION!!!!!!!!!            
   END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
contains

 
 
 
 
  
SUBROUTINE NONLIN_CORR_YZ_FDM_DIRECTION( YZ_out, tmpYZ, &
   Fp, Fyz, DFDXp, DFDXyz,  &
   Up, Uyz, DUDXp, DUDXyz  )
   
   IMPLICIT NONE
   
 

   REAL(wp), DIMENSION( :, : ), INTENT(inout) :: YZ_out 
   REAL(wp), DIMENSION( :, : ), INTENT(in)    :: Fp, DFDXp, Up, DUDXp
   REAL(wp), DIMENSION( :, : ), INTENT(in)    :: DFDXyz, DUDXyz, Fyz, Uyz
   REAL(wp), DIMENSION( :, : ), INTENT(out)   :: tmpYZ        
   REAL(wp), PARAMETER :: FAC4O = 1.0_wp / 12.0_wp    
   REAL(wp), PARAMETER :: FAC6O_Q22 = 36.0_wp  / 720.0_wp    
   REAL(wp), PARAMETER :: FAC6O_Q13 = 3.0_wp / 720.0_wp     
   real(wp) u_diff, df_diff, f_diff, du_diff , t1,t2
   integer :: i, j
real(wp) :: u_diff0, df_diff0, f_diff0, du_diff0
real(wp) :: u_diff1, df_diff1, f_diff1, du_diff1

 

! Scalar registers for pre-fetched values
real(wp) :: u_p0, u_p1, df_p0, df_p1
real(wp) :: f_p0, f_p1, du_p0, du_p1

! Pre-computed differences and sums
real(wp) :: u_sum,  df_sum,  f_sum,  du_sum

 
 IF( METHNC .EQ. 4 ) THEN

      CALL DFpDY_x_DGpDY_YZ( Fp, DUDXp ,  YZ_out ) 
      CALL DFpDY_x_DGpDY_YZ( Up, DFDXp ,  YZ_out )  
   
ELSE

	do i = iF0e, iFLe
	   do j = jC0e, jCLe

	      ! 1. Fetch values ONCE into CPU registers
	      u_p0  = Up(j, i)
	      u_p1  = Up(j+1, i)
	      df_p0 = DFDXp(j, i)
	      df_p1 = DFDXp(j+1, i)

	      f_p0  = Fp(j, i)
	      f_p1  = Fp(j+1, i)
	      du_p0 = DUDXp(j, i)
	      du_p1 = DUDXp(j+1, i)

	      ! 2. Compute differences and 2*average terms in registers
	      u_diff  = u_p1 - u_p0
	      df_diff = df_p1 - df_p0
	      f_diff  = f_p1 - f_p0
	      du_diff = du_p1 - du_p0

	      u_sum  = u_p1 + u_p0 - 2.0_wp * Uyz(j, i)
	      df_sum = df_p1 + df_p0 - 2.0_wp * DFDXyz(j, i)
	      f_sum  = f_p1 + f_p0 - 2.0_wp * Fyz(j, i)
	      du_sum = du_p1 + du_p0 - 2.0_wp * DUDXyz(j, i)

	      ! 3. Refactored single-assignment evaluation
	      YZ_out(j, i) = YZ_out(j, i) - ( FAC4O * ( f_diff * du_diff + u_diff * df_diff ) &
		 ! Grouped Q22 terms
		 + FAC6O_Q22 * ( u_sum * df_sum + f_sum * du_sum ) &
		 ! Grouped Q31 terms
		 - FAC6O_Q13 * ( &
		     ( Uyz(j+1, i) - Uyz(j-1, i) - 2.0_wp * u_diff ) * df_diff + &
		     ( DFDXyz(j+1, i) - DFDXyz(j-1, i) - 2.0_wp * df_diff ) * u_diff + &
		     ( Fyz(j+1, i) - Fyz(j-1, i) - 2.0_wp * f_diff ) * du_diff + &
		     ( DUDXyz(j+1, i) - DUDXyz(j-1, i) - 2.0_wp * du_diff ) * f_diff &
		   )  )

	   end do
	end do
 

                       
   ENDIF

   
 
   

END SUBROUTINE NONLIN_CORR_YZ_FDM_DIRECTION

!*********************************************************************
! =====================================================================
! OVERLOADED SET 2: XZ VARIADIC FUNCTIONS
! =====================================================================
SUBROUTINE NONLIN_CORR_XZ_FDM_DIRECTION( XZ_out, tmpXZ, &
   Fp, Fxz, DFDYp, DFDYxz, &
   Vp, Vxz, DVDYp, DVDYxz  )
   IMPLICIT NONE   
   
   REAL(wp), DIMENSION( :, : ), INTENT(inout) :: XZ_out      
   REAL(wp), DIMENSION( :, : ), INTENT(out)   :: tmpXZ 
   REAL(wp), DIMENSION( :, : ),    INTENT(in)    :: Fp, DVDYp, Vp, DFDYp   
   REAL(wp), DIMENSION( :, : ), INTENT(in)  :: DFDYxz, DVDYxz, Fxz, Vxz
   REAL(wp), PARAMETER :: FAC4O = 1.0_wp / 12.0_wp    
   REAL(wp), PARAMETER :: FAC6O_Q22 = 36.0_wp  / 720.0_wp    
   REAL(wp), PARAMETER :: FAC6O_Q13 = 3.0_wp / 720.0_wp      
 
   REAL(wp)  t1,t2     , v_diff, dfdy_diff, f_diff, dvdy_diff
   integer :: i, j

real(wp) :: v_p0, v_p1, f_p0, f_p1

real(wp) :: v_sum,  dfdy_sum,  f_sum,  dvdy_sum

 
   IF( METHNC .EQ. 4 ) THEN
!!      call CPU_TIME( T1 )
      CALL DFpDX_x_DGpDX_XZ( Fp, DVDYp , XZ_out ) 
      CALL DFpDX_x_DGpDX_XZ( Vp, DFDYp , XZ_out )
 
 
 
   ELSE 

	! Pre-calculated scalar terms held in CPU registers


	do i = iC0e, iCLe
	   do j = jF0e, jFLe

	      ! 1. Fetch values once into registers
	      v_p0 = Vp(j, i)
	      v_p1 = Vp(j, i+1)
	      f_p0 = Fp(j, i)
	      f_p1 = Fp(j, i+1)

	      ! 2. Pre-compute reusable X-differences & sums
	      v_diff    = v_p1 - v_p0
	      f_diff    = f_p1 - f_p0
	      dfdy_diff = DFDYp(j, i+1) - DFDYp(j, i)
	      dvdy_diff = DVDYp(j, i+1) - DVDYp(j, i)

	      v_sum    = v_p1 + v_p0 - 2.0_wp * Vxz(j, i)
	      f_sum    = f_p1 + f_p0 - 2.0_wp * Fxz(j, i)
	      dfdy_sum = DFDYp(j, i+1) + DFDYp(j, i) - 2.0_wp * DFDYxz(j, i)
	      dvdy_sum = DVDYp(j, i+1) + DVDYp(j, i) - 2.0_wp * DVDYxz(j, i)

	      ! 3. Clean, register-bound arithmetic
	      XZ_out(j, i) =XZ_out(j, i) -( FAC4O * ( v_diff * dfdy_diff + f_diff * dvdy_diff ) &
		 ! Q22--1 & Q22--2
		 + FAC6O_Q22 * ( v_sum * dfdy_sum + f_sum * dvdy_sum ) &
		 ! Q31--1 (V & DFDY)
		 - FAC6O_Q13 * ( ( Vxz(j, i+1) - Vxz(j, i-1) - 2.0_wp * v_diff ) * dfdy_diff &
		               + ( DFDYxz(j, i+1) - DFDYxz(j, i-1) - 2.0_wp * dfdy_diff ) * v_diff ) &
		 ! Q31--2 (F & DVDY)
		 - FAC6O_Q13 * ( ( Fxz(j, i+1) - Fxz(j, i-1) - 2.0_wp * f_diff ) * dvdy_diff &
		               + ( DVDYxz(j, i+1) - DVDYxz(j, i-1) - 2.0_wp * dvdy_diff ) * f_diff ) )

	   end do
	end do

 
 
   
 
          
ENDIF
           
END SUBROUTINE NONLIN_CORR_XZ_FDM_DIRECTION





SUBROUTINE NONLIN_CORR_YZ_FVM_DIRECTION( YZ_out, tmpYZ, &
   Fp,      Fyz,   &
   Up,      Uyz  )
   
   REAL(wp), DIMENSION( :, : ), INTENT(inout) :: YZ_out 
   REAL(wp), DIMENSION( :, : ), INTENT(in)    :: Fp, Up
   REAL(wp), DIMENSION( :, : ), INTENT(in)    :: Fyz, Uyz
   REAL(wp), DIMENSION( :, : ), INTENT(out)   :: tmpYZ 
   REAL(wp), PARAMETER :: FAC4O = 1.0_wp / 12.0_wp          

   ! array for storing accumulate correction
 
   IF( METHNC .EQ. 4 ) THEN
      CALL DFpDY_x_DGpDY_YZ( Fp, Up ,   YZ_out ) 
   ENDIF

 
END SUBROUTINE NONLIN_CORR_YZ_FVM_DIRECTION




SUBROUTINE NONLIN_CORR_XZ_FVM_DIRECTION( XZ_out, tmpXZ, &
   Fp, Fxz, &
   Vp, Vxz)

   REAL(wp), DIMENSION(:, :), INTENT(inout) :: XZ_out
   REAL(wp), DIMENSION(:, :),    INTENT(in)    :: Fp, Vp
   REAL(wp), DIMENSION(:, :),    INTENT(in)    :: Fxz, Vxz
   REAL(wp), DIMENSION(:, :), INTENT(out)   :: tmpXZ 

 
   IF (METHNC == 4) THEN
      CALL DFpDX_x_DGpDX_XZ( Fp, Vp, XZ_out )
   END IF

 
END SUBROUTINE NONLIN_CORR_XZ_FVM_DIRECTION
 
 

SUBROUTINE Q22_YZ( Fp, Fyz, Gp, Gyz, NC6_YZ ) 
   REAL(wp), DIMENSION( NPytot,    NPxtot ), INTENT(in)    :: Fp, Gp      
   REAL(wp), DIMENSION( NcellYtot, NPxtot ), INTENT(in)    :: Fyz, Gyz 
   REAL(wp), DIMENSION( NcellYtot, NPxtot ), INTENT(inout) :: NC6_YZ	
   REAL(wp), PARAMETER :: FAC6O_Q22 = 36.0_wp  / 720.0_wp 	 
      
   NC6_YZ(jC0e:jCLe,iF0e:iFLe) = &
   NC6_YZ(jC0e:jCLe, iF0e:iFLe) + FAC6O_Q22 * &
              (         Fp(jC0eP1:jCLeP1, iF0e:iFLe) &
                       +Fp(jC0e  :jCLe  , iF0e:iFLe) &
               -2.0_wp*Fyz(jC0e  :jCLe  , iF0e:iFLe) ) * &
              (        +Gp(jC0eP1:jCLeP1, iF0e:iFLe) &
                       +Gp(jC0e  :jCLe  , iF0e:iFLe) &
               -2.0_wp*Gyz(jC0e  :jCLe  , iF0e:iFLe) )                         
END SUBROUTINE Q22_YZ  
      

SUBROUTINE Q22_XZ( Fp, Fxz, Gp, Gxz, NC6_XZ )
   REAL(wp), DIMENSION(NPytot, NPXtot),     INTENT(in)    :: Fp, Gp
   REAL(wp), DIMENSION(NPytot, NcellXtot),  INTENT(in)    :: Fxz, Gxz
   REAL(wp), DIMENSION(NPytot, NcellXtot),  INTENT(inout) :: NC6_XZ
   REAL(wp), PARAMETER :: FAC6O_Q22 = 36.0_wp / 720.0_wp

   NC6_XZ(jF0e:jFLe, iC0e:iCLe) = &
   NC6_XZ(jF0e:jFLe, iC0e:iCLe) + &
      ( &
                      Fp( jF0e:jFLe, iC0eP1:iCLeP1) &
                    + Fp( jF0e:jFLe, iC0e   :iCLe ) &
          - 2.0_wp * Fxz( jF0e:jFLe, iC0e:iCLe) &
      ) * &
      (   &
                      Gp(jF0e:jFLe, iC0eP1:iCLeP1) &
          +           Gp(jF0e:jFLe, iC0e   :iCLe ) &
          - 2.0_wp * Gxz(jF0e:jFLe, iC0e:iCLe) &
      ) * FAC6O_Q22 
END SUBROUTINE Q22_XZ         


SUBROUTINE Q31_YZ( Fp, Fyz,   Gp, Gyz,  NC6_YZ )          
   REAL(wp), DIMENSION( NPYtot,    NPxtot ), INTENT(in)    :: Fp, Gp      
   REAL(wp), DIMENSION( NcellYtot, NPxtot ), INTENT(in)    :: Fyz, Gyz
   REAL(wp), DIMENSION( NcellYtot, NPxtot ), INTENT(inout) :: NC6_YZ	  
   REAL(wp), PARAMETER :: FAC6O_Q13 = 3.0_wp / 720.0_wp 	 


!!-2* ( Qp[i+1] -Qp[i] )  +Qyz[i+1]+-Qyz[i-1]

   NC6_YZ(jC0e:jCLe, iF0e:iFLe) = &
   NC6_YZ(jC0e:jCLe, iF0e:iFLe) - FAC6O_Q13 * ( &
       ( &
             Fyz(jC0eP1:jCLeP1, iF0e:iFLe) &
           - Fyz(jC0eM1:jCLeM1, iF0e:iFLe) &
           - 2.0_wp * ( &
                 Fp(jC0eP1:jCLeP1, iF0e:iFLe) &
               - Fp(jC0e   :jCLe,   iF0e:iFLe) &
        ) &
        ) * (Gp(jC0eP1:jCLeP1, iF0e:iFLe) - Gp(jC0e:jCLe, iF0e:iFLe) )   &       
       + &
       ( &
             Gyz(jC0eP1:jCLeP1, iF0e:iFLe) &
           - Gyz(jC0eM1:jCLeM1, iF0e:iFLe) &
           - 2.0_wp * ( &
                 Gp(jC0eP1:jCLeP1, iF0e:iFLe) &
               - Gp(jC0e   :jCLe,   iF0e:iFLe) &
        ) &
        ) * (Fp(jC0eP1:jCLeP1, iF0e:iFLe) - Fp(jC0e:jCLe, iF0e:iFLe) )   &            
   )
END SUBROUTINE Q31_YZ                        


SUBROUTINE Q31_XZ( Fp, Fxz,   Gp, Gxz,  NC6_XZ )
   IMPLICIT NONE
   REAL(wp), DIMENSION(NPYtot, NPXtot),    INTENT(in)    :: Fp, Gp
   REAL(wp), DIMENSION(NPYtot, NcellXtot), INTENT(in)    :: Fxz, Gxz
   REAL(wp), DIMENSION(NPYtot, NcellXtot), INTENT(inout) :: NC6_XZ
   REAL(wp), PARAMETER :: FAC6O_Q13 = +3.0_wp / 720.0_wp

   NC6_XZ(jF0e:jFLe, iC0e:iCLe) = &
   NC6_XZ(jF0e:jFLe, iC0e:iCLe) - FAC6O_Q13 * ( &
     ( &
         Fxz(jF0e:jFLe,iC0eP1:iCLeP1) &
       - Fxz(jF0e:jFLe,iC0eM1:iCLeM1) &
       - 2.0_wp*(  Fp(jF0e:jFLe,iC0eP1:iCLeP1)  - Fp(jF0e:jFLe,iC0e:iCLe) ) &
     ) *  (Gp(jF0e:jFLe, iC0eP1:iCLeP1) - Gp(jF0e:jFLe, iC0e:iCLe)) &
   + &
     ( &
         Gxz(jF0e:jFLe, iC0eP1:iCLeP1) &
       - Gxz(jF0e:jFLe, iC0eM1:iCLeM1) &
       - 2.0_wp*( &
           Gp(jF0e:jFLe, iC0eP1:iCLeP1) &
         - Gp(jF0e:jFLe, iC0e:iCLe) ) &
     ) *  (Fp(jF0e:jFLe, iC0eP1:iCLeP1) - Fp(jF0e:jFLe, iC0e:iCLe)) &
   ) 
END SUBROUTINE Q31_XZ


SUBROUTINE DFpDY_x_DGpDY_YZ( Fp, Gp,   DFpDY_x_DGpDY ) 
   IMPLICIT NONE
   REAL(wp), DIMENSION( NPytot,    NPxtot ), INTENT(in)    :: Fp, Gp    
   REAL(wp), DIMENSION( NcellYtot, NPxtot ), INTENT(inout) :: DFpDY_x_DGpDY 
   REAL(wp), PARAMETER :: FAC4O = 1.0_wp / 12.0_wp         

 

      DFpDY_x_DGpDY( jC0e:jCLe, iF0e:iFLe ) = DFpDY_x_DGpDY( jC0e:jCLe, iF0e:iFLe ) - &
         FAC4O * ( Fp( jC0eP1:jCLeP1, iF0e:iFLe ) - Fp( jC0e:jCLe, iF0e:iFLe ) ) * &
                 ( Gp( jC0eP1:jCLeP1, iF0e:iFLe ) - Gp( jC0e:jCLe, iF0e:iFLe ) )
 
END SUBROUTINE DFpDY_x_DGpDY_YZ
 

SUBROUTINE DFpDX_x_DGpDX_XZ( Fp, Gp,  DFpDX_x_DGpDX )
   IMPLICIT NONE
   REAL(wp), DIMENSION(NPYtot, NPXtot),    INTENT(in)    :: Fp, Gp
   REAL(wp), DIMENSION(NPYtot, NcellXtot), INTENT(inout) :: DFpDX_x_DGpDX
   REAL(wp), PARAMETER :: FAC4O = 1.0_wp / 12.0_wp
 
       DFpDX_x_DGpDX(jF0e:jFLe, iC0e:iCLe) = DFpDX_x_DGpDX(jF0e:jFLe, iC0e:iCLe) - &
         FAC4O * ( Fp(jF0e:jFLe, iC0eP1:iCLeP1) - Fp(jF0e:jFLe, iC0e:iCLe) ) * &
		         ( Gp(jF0e:jFLe, iC0eP1:iCLeP1) - Gp(jF0e:jFLe, iC0e:iCLe) )
 
END SUBROUTINE DFpDX_x_DGpDX_XZ 
 
   
   

END SUBROUTINE Compute_RHS
!*********************************************************
!*********************************************************





! =====================================================================
! HELPERS AND MEMORY CONSTRUCTORS/DESTRUCTORS
! =====================================================================
SUBROUTINE set_TimeINTmethod( method, TimeINTmethod )
	CHARACTER(len=16), INTENT(in) :: method 
	INTEGER, INTENT(out) :: TimeINTmethod
	
   SELECT CASE ( method )  
      CASE ( 'rk3' )	
         TimeINTmethod = 31
      CASE ( 'rk4' )
         TimeINTmethod = 41 
      CASE ( 'rk5' )
         TimeINTmethod = 51 
      CASE ( 'rk7' )
         TimeINTmethod = 71 
      CASE default
         ERROR STOP 'Unknown Time integration method: set_TimeINTmethod'
   END SELECT       		
END SUBROUTINE set_TimeINTmethod






SUBROUTINE ALLOCATE_TimeIntArray(this, TimeINTmethod , source_mesh )
   CLASS(TimeIntArray), INTENT(inout) :: this
   TYPE(MeshDomain),  INTENT(in), TARGET :: source_mesh
   INTEGER, INTENT(in) :: TimeINTmethod
   INTEGER :: i
   
   this%Method = TimeINTmethod
   
   SELECT CASE ( TimeINTmethod )  
      CASE ( 31 )
         this%Nstage = 3
         ALLOCATE( this%a( this%Nstage-1,  this%Nstage -1), &
                   this%b( this%Nstage  ),              &
                   this%c( this%Nstage  ) )   
         this%a(:,1) = [ 0.5_wp, 0.0_wp ]
         this%a(:,2) = [ -1.0_wp, 2.0_wp ]
         this%b = [ 1.0_wp/6.0_wp, 2.0_wp/3.0_wp, 1.0_wp/6.0_wp ]
         this%c = [ 0.0_wp, 0.5_wp, 1.0_wp ]        
                             	        
      CASE ( 41 )
         this%Nstage = 4
         ALLOCATE( this%a( this%Nstage-1,  this%Nstage -1), &
                   this%b( this%Nstage  ),              &
                   this%c( this%Nstage  ) )
         this%a(:,1) = [ 0.5_wp, 0.0_wp, 0.0_wp ]
         this%a(:,2) = [ 0.0_wp, 0.5_wp, 0.0_wp ]
         this%a(:,3) = [ 0.0_wp, 0.0_wp, 1.0_wp ]
         this%b = [ 1.0_wp/6.0_wp, 1.0_wp/3.0_wp, 1.0_wp/3.0_wp, 1.0_wp/6.0_wp ]
         this%c = [ 0.0_wp,        0.5_wp,        0.5_wp,        1.0_wp        ]

      CASE ( 51 )
         this%Nstage = 6    
         ALLOCATE( this%a( this%Nstage-1,  this%Nstage-1 ), &
                   this%b( this%Nstage  ),              &
                   this%c( this%Nstage  ) )		
         this%a(:,1) = [ 1.0_wp/4.0_wp,   0.0_wp,         0.0_wp,          0.0_wp,          0.0_wp ]
         this%a(:,2) = [ 1.0_wp/8.0_wp,   1.0_wp/8.0_wp,  0.0_wp,          0.0_wp,          0.0_wp ]
         this%a(:,3) = [ 0.0_wp,         -0.5_wp,         1.0_wp,          0.0_wp,          0.0_wp ]
         this%a(:,4) = [ 3.0_wp/16.0_wp,  0.0_wp,         0.0_wp,          9.0_wp/16.0_wp,  0.0_wp ]
         this%a(:,5) = [ -3.0_wp/7.0_wp,  2.0_wp/7.0_wp,  12.0_wp/7.0_wp, -12.0_wp/7.0_wp,  8.0_wp/7.0_wp ]
         this%b = [ 7.0_wp, 0.0_wp, 32.0_wp, 12.0_wp, 32.0_wp, 7.0_wp ] / 90.0_wp
         this%c = [ 0.0_wp, 1.0_wp/4.0_wp, 1.0_wp/4.0_wp, 0.5_wp, 3.0_wp/4.0_wp, 1.0_wp ]		

      CASE ( 71 ) 	
         this%Nstage = 9
         ALLOCATE( this%a( this%Nstage-1,  this%Nstage -1), &
                   this%b( this%Nstage  ),              &
                   this%c( this%Nstage  ) )      		
      CASE default
         ERROR STOP 'Unknown integration method.'
   END SELECT      		 

   ALLOCATE(this%K( this%Nstage ) )
   DO i=1,this%Nstage
      CALL this%K(i)%init( source_mesh )
   END DO
   
   CALL this%VAR0%init( source_mesh ) 
   CALL this%VAR1%init( source_mesh ) 
   CALL this%VAR2%init( source_mesh )
   CALL this%VAR3%init( source_mesh )
   CALL this%VAR4%init( source_mesh )      

   ALLOCATE( this%F0yz( source_mesh%NcellYtot, source_mesh%NPxtot ) ) 
   ALLOCATE( this%F1yz( source_mesh%NcellYtot, source_mesh%NPxtot ) )
   ALLOCATE( this%F2yz( source_mesh%NcellYtot, source_mesh%NPxtot ) )
   ALLOCATE( this%F3yz( source_mesh%NcellYtot, source_mesh%NPxtot ) )
   ALLOCATE( this%F4yz( source_mesh%NcellYtot, source_mesh%NPxtot ) )
   ALLOCATE( this%F5yz( source_mesh%NcellYtot, source_mesh%NPxtot ) )  

   ALLOCATE( this%F0xz( source_mesh%NPytot, source_mesh%NcellXtot ) )     
   ALLOCATE( this%F1xz( source_mesh%NPytot, source_mesh%NcellXtot ) )
   ALLOCATE( this%F2xz( source_mesh%NPytot, source_mesh%NcellXtot ) )
   ALLOCATE( this%F3xz( source_mesh%NPytot, source_mesh%NcellXtot ) )
   ALLOCATE( this%F4xz( source_mesh%NPytot, source_mesh%NcellXtot ) )
   ALLOCATE( this%F5xz( source_mesh%NPytot, source_mesh%NcellXtot ) )    
END SUBROUTINE ALLOCATE_TimeIntArray


SUBROUTINE deALLOCATE_TimeIntArray(this)
   CLASS(TimeIntArray), INTENT(inout) :: this
  
   IF ( ALLOCATED(this%K) ) THEN
      DEALLOCATE( this%K )   
      CALL this%VAR0%free()
      CALL this%VAR1%free(); CALL this%VAR2%free()
      CALL this%VAR3%free(); CALL this%VAR4%free()      
		
      DEALLOCATE( this%F1yz, this%F2yz, this%F3yz, this%F4yz, this%F5yz )
      DEALLOCATE( this%F1xz, this%F2xz, this%F3xz, this%F4xz, this%F5xz )			
   ENDIF
END SUBROUTINE deALLOCATE_TimeIntArray

END MODULE MMM_module_TimeINT
