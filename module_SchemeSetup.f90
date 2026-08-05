module module_SchemeSetup
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
     
    public :: wp, TimeINTmethod, SptScheme, FD_scheme, SD_scheme, METHNC
    public :: BWFD, BWSD,  FDB1BW, FDB2BW, SDB1BW, SDB2BW, isMM, NG
    public :: setBandWidth_coeff ! Added so it can be called elsewhere

    integer :: TimeINTmethod, SptScheme, FD_scheme, SD_scheme, METHNC, BWFD, BWSD
    integer :: FDB1BW, FDB2BW, SDB1BW, SDB2BW , NG
    logical :: isMM
    integer, parameter :: wp = real64

contains

    SUBROUTINE setBandWidth_coeff( inputNG )
    integer, intent(in) :: inputNG
    
         NG = inputNG
        ! Default fallback reset
        FDB1BW = 0; FDB2BW = 0
        SDB1BW = 0; SDB2BW = 0

        select case ( SptScheme ) 
        case(3002) 
            FDB1BW = 1; FDB2BW = 0
            SDB1BW = 1; SDB2BW = 0
        case(3004) 
            FDB1BW = 1; FDB2BW = 1
            SDB1BW = 1; SDB2BW = 1                        
        case(3006)
            FDB1BW = 2; FDB2BW = 1
            SDB1BW = 2; SDB2BW = 1   
        case(3008)    
            FDB1BW = 2; FDB2BW = 2
            SDB1BW = 2; SDB2BW = 2   
        case(3010)                    
            FDB1BW = 3; FDB2BW = 2
            SDB1BW = 3; SDB2BW = 2   
        case(3012)                    
            FDB1BW = 3; FDB2BW = 3
            SDB1BW = 3; SDB2BW = 3 
        ! FDM ------------------------------------       
        case(7002) 
            FDB2BW = 1; SDB2BW = FDB2BW
        case(7004) 
            FDB2BW = 2; SDB2BW = FDB2BW                     
        case(7006)
            FDB2BW = 3; SDB2BW = FDB2BW   
        case(7008)    
            FDB2BW = 4; SDB2BW = FDB2BW   
        case(7010)                    
            FDB2BW = 5; SDB2BW = FDB2BW   
        case(7012)                    
            FDB2BW = 6; SDB2BW = FDB2BW                                        
        case default
            error stop 'Unknown Spatial Scheme.'
        end select           
    
    END SUBROUTINE setBandWidth_coeff

end module module_SchemeSetup



module matlab_io_mod
    implicit none
    private
    
    ! Define working precision (Double Precision)
    integer, parameter, public :: wp = kind(1.0d0)
    
    public :: save_bin_matlab

contains

    subroutine save_bin_matlab(filename, var)
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: var(:,:) ! Assumed-shape 2D array
        
        integer :: file_unit, iostat
        character(len=256) :: iomsg

        ! Open file in raw binary stream mode
        open(newunit=file_unit, &
             file=trim(filename), &
             access='stream', &
             status='replace', &
             iostat=iostat, &
             iomsg=iomsg)
        
        ! Error checking for file opening
        if (iostat /= 0) then
            print *, "Error opening file: ", trim(filename)
            print *, trim(iomsg)
            return
        end if
        
        ! Write the entire matrix natively in column-major layout
        write(file_unit, iostat=iostat, iomsg=iomsg) var
        
        ! Error checking for writing data
        if (iostat /= 0) then
            print *, "Error writing matrix to file: ", trim(filename)
            print *, trim(iomsg)
        end if
        
        close(file_unit)
    end subroutine save_bin_matlab

end module matlab_io_mod

module module_physics
    use module_SchemeSetup, only:wp
    implicit none
    public :: KDIFF

    REAL (WP) :: KDIFF
end module
