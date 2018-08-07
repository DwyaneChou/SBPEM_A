MODULE module_para
    
    implicit none
    
    real*8,  parameter     :: omg0 = 7.292d-5		     ! angular velocity of the earth rotation
	real*8,  parameter     :: a    = 6371220d0		     ! radius of the earth
    real*8,  parameter     :: g    = 9.80616d0		     ! gravity acceleration
                                                        
    integer                :: nx                         ! zonal grid nmuber 
    integer                :: nx1
                                                        
    integer                :: ny                         ! meridional grid number
    integer                :: ny1
                                                        
	integer, parameter     :: t_start     = 0		     ! initial time
	integer                :: t_end                      ! final time
                              
	integer, parameter     :: case_select = 0            ! 0: using RH waves as initial conditon (IC); 1: read IC from files
                                                    
	real*8                 :: pi					     ! pi = datan(1d0)*4.d0 = 3.141592653589793238462643384279
    
    integer,parameter      :: full_pass = 0
    integer,parameter      :: fast_pass = 1
    integer,parameter      :: slow_pass = 2

    !namelist variables
    real*8                 :: zonal_resolution
    real*8                 :: meridional_resolution
    
	integer                :: time_step                  ! time stepsize
    integer                :: run_days
    integer                :: run_hours
    integer                :: run_minutes
    integer                :: run_seconds
    integer                :: history_interval           ! The time for saving the intermediate result
    
    character*20           :: integral_scheme
    character*20           :: split_scheme
    integer                :: split_step
    real*8                 :: split_time_step
    
    namelist /mesh_settings/     zonal_resolution,meridional_resolution
    namelist /time_settings/     time_step,run_days,run_hours,run_minutes,run_seconds,history_interval
    namelist /integral_settings/ integral_scheme,split_scheme,split_step
    
    contains
    
    subroutine read_namelist
        implicit none
        character*200 namelist_file
        
        namelist_file = trim(adjustl('namelist.input'))
        
        open(1,file = namelist_file)
        read(1,nml=time_settings    )
        read(1,nml=mesh_settings    )
        read(1,nml=integral_settings)
        close(1)
        
    end subroutine read_namelist
    
    subroutine preprocess_parameter
        implicit none
        real*8 d2r
        
        t_end = run_days*24*3600 + run_hours*3600 + run_minutes*60 + run_seconds
        
        pi    = datan(1d0)*4.d0
              
        d2r   = pi/180.d0
              
        nx    = int(360.d0/zonal_resolution)
        ny    = int(180.d0/meridional_resolution)+1
              
        nx1   = nx+1
        ny1   = ny+1
        
        split_time_step = dble(time_step)/split_step
    
    end subroutine preprocess_parameter
    
END MODULE module_para
    
