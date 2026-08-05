!===============================================================================
! MODULE: VTK_MultiMoment_Writer
!
! Writes the four multi-moment data locations (Fp, Fxz, Fyz, Fxyz) of a
! scalar field to a single .vtm (VTK MultiBlock) file. Each location is
! stored as a separate VTK RectilinearGrid (.vtr) piece, since each field
! lives on a different (staggered) subgrid.
!
! Mapping between this code's field names, the multi-moment picture, and
! the underlying coordinate lines (Xp/Yp = node/point locations,
! Xc/Yc = cell-center locations):
!
!   code name   picture name    located at (x, y)       array shape (j,i)
!   ---------   ------------    -----------------        -----------------
!   F%p         phi_{i,j}       (Xp(i),  Yp(j))           (NPYtot,    NPXtot)
!   F%yz        phi^y_{i,jc}    (Xp(i),  Yc(jc))          (NcellYtot, NPXtot)
!   F%xz        phi^x_{ic,j}    (Xc(ic), Yp(j))           (NPYtot,    NcellXtot)
!   F%xyz       phi^xy_{ic,jc}  (Xc(ic), Yc(jc))          (NcellYtot, NcellXtot)
!
! All arrays are assumed to follow this codebase's (j,i) = (y-index, x-index)
! convention, i.e. j (first index) is the y-related index, i (second index)
! is the x-related index. If your F%p etc. are actually stored (i,j)
! instead, transpose the arguments before calling write_vtr_2d, or adjust
! the loop order inside it.
! 
! Output format is ASCII for simplicity/portability. For very large grids
! (e.g. 1000x1000+, called every few time steps) ASCII I/O becomes a real
! cost (large files, slow write). If this output routine is called often,
! consider switching the DataArray "format" to "appended" with raw binary
! payloads, or writing fewer snapshots.
!===============================================================================
MODULE VTK_MultiMoment_Writer

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: Write_MultiMoment_VTM

   INTEGER, PARAMETER :: wp = KIND(1.0d0)   ! <-- match this to your module_Precision's wp

CONTAINS

 
   !---------------------------------------------------------------------------
   ! Top-level driver: writes Fp, Fyz, Fxz, Fxyz as 4 .vtr pieces + 1 .vtm
   !
   ! Example call:
   !   CALL Write_MultiMoment_VTM( F%p, F%yz, F%xz, F%xyz, &
   !                               iMesh%Xp, iMesh%Yp, iMesh%Xc, iMesh%Yc, &
   !                               'output/step_0001' )
   !---------------------------------------------------------------------------
   SUBROUTINE Write_MultiMoment_VTM( Fp, Fyz, Fxz, Fxyz, Xp, Yp, Xc, Yc, basename )
      REAL(wp), DIMENSION(:,:), INTENT(in) :: Fp     ! (NPYtot,    NPXtot)
      REAL(wp), DIMENSION(:,:), INTENT(in) :: Fyz    ! (NcellYtot, NPXtot)
      REAL(wp), DIMENSION(:,:), INTENT(in) :: Fxz    ! (NPYtot,    NcellXtot)
      REAL(wp), DIMENSION(:,:), INTENT(in) :: Fxyz   ! (NcellYtot, NcellXtot)
      REAL(wp), DIMENSION(:),   INTENT(in) :: Xp, Yp ! node coordinates
      REAL(wp), DIMENSION(:),   INTENT(in) :: Xc, Yc ! cell-center coordinates
      CHARACTER(*), INTENT(in) :: basename           ! e.g. "output/step_0001"

      CHARACTER(LEN=256) :: fname_p, fname_yz, fname_xz, fname_xyz, fname_vtm
      INTEGER :: unit_vtm

      fname_p   = TRIM(basename)//'_Fp.vtr'
      fname_yz  = TRIM(basename)//'_Fyz.vtr'
      fname_xz  = TRIM(basename)//'_Fxz.vtr'
      fname_xyz = TRIM(basename)//'_Fxyz.vtr'
      fname_vtm = TRIM(basename)//'.vtm'

      ! 1. Write each field on its own subgrid
      CALL write_vtr_2d( fname_p,   Xp, Yp, Fp,   'Fp'   )
      CALL write_vtr_2d( fname_yz,  Xp, Yc, Fyz,  'Fyz'  )
      CALL write_vtr_2d( fname_xz,  Xc, Yp, Fxz,  'Fxz'  )
      
            !print *, "just b4 ", SIZE(Fxyz,1), SIZE(Fxyz,2)    
      
      CALL write_vtr_2d( fname_xyz, Xc, Yc, Fxyz, 'Fxyz' )

      ! 2. Write the multiblock container that ties them together
      OPEN( NEWUNIT=unit_vtm, FILE=TRIM(fname_vtm), STATUS='replace', ACTION='write' )
      WRITE(unit_vtm,'(A)') '<?xml version="1.0"?>'
      WRITE(unit_vtm,'(A)') '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
      WRITE(unit_vtm,'(A)') '  <vtkMultiBlockDataSet>'
      WRITE(unit_vtm,'(A)') '    <DataSet index="0" name="Fp"   file="'//TRIM(get_basefile(fname_p))//'"/>'
      WRITE(unit_vtm,'(A)') '    <DataSet index="1" name="Fyz"  file="'//TRIM(get_basefile(fname_yz))//'"/>'
      WRITE(unit_vtm,'(A)') '    <DataSet index="2" name="Fxz"  file="'//TRIM(get_basefile(fname_xz))//'"/>'
      WRITE(unit_vtm,'(A)') '    <DataSet index="3" name="Fxyz" file="'//TRIM(get_basefile(fname_xyz))//'"/>'
      WRITE(unit_vtm,'(A)') '  </vtkMultiBlockDataSet>'
      WRITE(unit_vtm,'(A)') '</VTKFile>'
      CLOSE(unit_vtm)

   END SUBROUTINE Write_MultiMoment_VTM

   !---------------------------------------------------------------------------
   ! Strips any directory path so the .vtm references the .vtr files by name
   ! only (assumes .vtr files sit alongside the .vtm; adjust if you keep them
   ! in a different subdirectory than the .vtm itself).
   !---------------------------------------------------------------------------
   FUNCTION get_basefile( fullpath ) RESULT( base )
      CHARACTER(*), INTENT(in) :: fullpath
      CHARACTER(LEN=256) :: base
      INTEGER :: k
      k = INDEX( fullpath, '/', BACK=.true. )
      base = fullpath(k+1:)
   END FUNCTION get_basefile

   !---------------------------------------------------------------------------
   ! Writes one 2D scalar field as a VTK RectilinearGrid (.vtr), ASCII format.
   !
   ! Xcoord: x-coordinate line, size nx
   ! Ycoord: y-coordinate line, size ny
   ! data2d: field values, shape (ny, nx)  -- (j,i) = (y,x) convention
   !---------------------------------------------------------------------------
   
   SUBROUTINE write_vtr_2d( filename, Xcoord, Ycoord, data2d, fieldname )
      CHARACTER(*), INTENT(in) :: filename, fieldname
      REAL(wp), DIMENSION(:),   INTENT(in) :: Xcoord, Ycoord
      REAL(wp), DIMENSION(:,:), INTENT(in) :: data2d

      INTEGER :: nx, ny, i, j, unit_vtr

      nx = SIZE(Xcoord)
      ny = SIZE(Ycoord)
      
!      print *, SIZE(Xcoord), SIZE(Ycoord)
!      print *, SIZE(data2d,1), SIZE(data2d,2)      

!      IF ( SIZE(data2d,1) /= ny .OR. SIZE(data2d,2) /= nx ) THEN
!         WRITE(*,*) 'write_vtr_2d ERROR: shape mismatch for ', TRIM(fieldname), &
!                    ' expected (',ny,',',nx,') got (',SIZE(data2d,1),',',SIZE(data2d,2),')'
!         RETURN
!      END IF

      OPEN( NEWUNIT=unit_vtr, FILE=TRIM(filename), STATUS='replace', ACTION='write' )

      WRITE(unit_vtr,'(A)') '<?xml version="1.0"?>'
      WRITE(unit_vtr,'(A)') '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'
      WRITE(unit_vtr,'(A,I0,A,I0,A)') '  <RectilinearGrid WholeExtent="0 ', nx-1, ' 0 ', ny-1, ' 0 0">'
      WRITE(unit_vtr,'(A,I0,A,I0,A)') '    <Piece Extent="0 ', nx-1, ' 0 ', ny-1, ' 0 0">'

      ! --- Coordinates (X, Y, and a single Z=0 plane) ---
      WRITE(unit_vtr,'(A)') '      <Coordinates>'
      WRITE(unit_vtr,'(A)') '        <DataArray type="Float64" Name="X" format="ascii">'
      WRITE(unit_vtr,'(6(1PE16.8))') ( Xcoord(i), i=1,nx )
      WRITE(unit_vtr,'(A)') '        </DataArray>'
      WRITE(unit_vtr,'(A)') '        <DataArray type="Float64" Name="Y" format="ascii">'
      WRITE(unit_vtr,'(6(1PE16.8))') ( Ycoord(j), j=1,ny )
      WRITE(unit_vtr,'(A)') '        </DataArray>'
      WRITE(unit_vtr,'(A)') '        <DataArray type="Float64" Name="Z" format="ascii">'
      WRITE(unit_vtr,'(1PE16.8)') 0.0_wp
      WRITE(unit_vtr,'(A)') '        </DataArray>'
      WRITE(unit_vtr,'(A)') '      </Coordinates>'

      ! --- Point data (the scalar field itself) ---
      WRITE(unit_vtr,'(A)') '      <PointData Scalars="'//TRIM(fieldname)//'">'
      WRITE(unit_vtr,'(A)') '        <DataArray type="Float64" Name="'//TRIM(fieldname)//'" format="ascii">'
      ! VTK expects data ordered with X varying fastest, then Y, then Z
      DO j = 1, ny
         WRITE(unit_vtr,'(6(1PE16.8))') ( data2d(j,i), i=1,nx )
      END DO
      WRITE(unit_vtr,'(A)') '        </DataArray>'
      WRITE(unit_vtr,'(A)') '      </PointData>'

      WRITE(unit_vtr,'(A)') '    </Piece>'
      WRITE(unit_vtr,'(A)') '  </RectilinearGrid>'
      WRITE(unit_vtr,'(A)') '</VTKFile>'

      CLOSE(unit_vtr)

   END SUBROUTINE write_vtr_2d

END MODULE VTK_MultiMoment_Writer
