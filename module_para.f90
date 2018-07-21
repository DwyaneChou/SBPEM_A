MODULE module_para
    
    implicit none
    
    real*8,  parameter       ::    omg0 = 7.292d-5		     ! angular velocity of the earth rotation
	real*8,  parameter       ::    a    = 6371220d0		     ! radius of the earth
    real*8,  parameter       ::    g    = 9.80616d0		     ! gravity acceleration
                                                             
    integer, parameter       ::    nx          = 360         ! zonal grid nmuber 
    integer, parameter       ::    nx1         = nx+1        
                                                             
    integer, parameter       ::    ny          = 181         ! meridional grid number
    integer, parameter       ::    ny1         = ny+1 		  	   
                                                             
	integer, parameter       ::    time_step   = 15  	     ! time stepsize
	integer, parameter       ::    t_start     = 0		     ! initial time
	integer, parameter       ::    t_end       = 66*24*3600  ! final time
    
	integer, parameter       ::    case_select = 0           ! 0: using RH waves as initial conditon (IC); 1: read IC from files
	integer, parameter       ::    thalf       = 3600        ! The time for saving the intermediate result
                                                             
	character*7, parameter   ::    fu='uui.dat' 		     ! initial field of zonal wind for self-defined experiment
	character*7, parameter   ::    fv='vvi.dat'	   	         ! initial field of meridional wind for self-defined experiment
	character*7, parameter   ::    fh='hhi.dat'	 	         ! initial field of geopotential height for self-defined experiment
                                                             
	real*8                   ::    pi					     ! pi = datan(1d0)*4.d0 = 3.141592653589793238462643384279

END MODULE module_para
    
