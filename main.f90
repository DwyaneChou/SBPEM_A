!*****************************************************************************!
!			       The Barotropic Model on Arakawa A-grid					  !
!                               by Bin Wang	                                  !
!	           using the constant-time-stepsize explicit scheme 			  !
!                       with energy conservation                              !
!*****************************************************************************!	
!
PROGRAM MAIN

    use module_para
    use module_array
    use mesh
    use output_netCDF
!
    implicit none
!
    real*8                       :: tener,tener0               ! total energy at tn and t_start, respectively
    real*8                       :: tmass,tmass0               ! total mass at tn and t_start, respectively
                                                               
	real*4, dimension(1:nx,1:ny) :: pu,pv,ph                   ! for output
                                                               
	real*8                       :: dt,dtn,t_now               ! working variables
	integer                      :: i,j                        ! working variables
                                                               
	real*8, external             :: inner                      ! a external function to calculate inner product
                                 
    character*4                  :: output_row_length
    logical                      :: do_output
    integer                      :: output_num,output_idx
    integer                      :: time_1,time_2

    call SYSTEM_CLOCK(time_1)
    
    t_now = t_start
    
    output_num = t_end/thalf+1
    output_idx = 1
    
    print *,'Initial time is',t_now
	print *,'Final time is',t_end
	print *,'Time stepsize is',time_step

	if (case_select.eq.0) then
	   print *,'This is a RH-wave experiment,'
	   print *,'i.e., the IC is RH-wave'
    else
	   print *,'This is self-defined experiment,'
	   print *,'i.e., the IC is read from files'
    end if
    
!   generate mesh
    call generate_mesh
    
!	initial condition
!
	if (case_select.eq.0) then
    !Using Rossby-Haurwitz waves as initial condition
        call haurwitz
    else
        
    end if
    
    !IAP tranformation
    h  = dsqrt(wh)
    wu = u*h
    wv = v*h
    
    pu = u (1:nx,1:ny)
    pv = v (1:nx,1:ny)
    ph = wh(1:nx,1:ny)
    
    wu(0  ,:) = wu(nx,:)
    wu(nx1,:) = wu(1 ,:)
    
    wv(0  ,:) = wv(nx,:)
    wv(nx1,:) = wv(1 ,:)
    
    wh(0  ,:) = wh(nx,:)
    wh(nx1,:) = wh(1 ,:)
    
    call write_netCDF(pu,pv,ph,output_num,output_idx)

    ! Compute total energy and total mass
	tener0=inner(wu,wv,wh,wu,wv,wh)

    tmass0 = 0
	do j=1,ny
        do i=1,nx
            tmass0=tmass0+wh(i,j)*cos_lat(i,j)
        end do
	end do

    print *,'The total energy is ',tener0
    print *,'The total mass is ',tmass0

	dt =time_step
    dtn=time_step
!
 	print *,'The main part of this program has started'
!
    open(unit=21,file='out.txt')
!
!------------------------------------------------------
!	The main loop for time integration
!------------------------------------------------------
!
	do while(t_now<t_end)
        t_now = t_now+dt

!------------------------------------------------
!       The time integration
!------------------------------------------------
        call vary(dt,dtn)
        
        h = dsqrt(wh)
        u = wu/h
        v = wv/j
        
        do_output=(t_now.ge.thalf).and.(mod(t_now,dble(thalf))==0)
        if (do_output) then
            
            print *,    '-------------------------------------------------------------------------------'
            write(21,*) '-------------------------------------------------------------------------------'
            
            print *,    '       (The integral time is      ',t_now,')'
            write(21,*) '       (The integral time is      ',t_now,')'
            
            print *,    '    The Energy                The Total-Mass            The timestep'
            write(21,*) '    The Energy                The Total-Mass            The timestep'
            
            tener=inner(wu,wv,wh,wu,wv,wh)
            
            tmass = 0
            do j=1,ny
                do i=1,nx
                    tmass=tmass+wh(i,j)*cos_lat(i,j)
                end do
            end do
            
            print*,     tener,tmass,dtn
            write(21,*) tener,tmass,dtn
            
            ! Output data into netCDF file
            pu = u (1:nx,1:ny)
            pv = v (1:nx,1:ny)
            ph = wh(1:nx,1:ny)
            
            print*,'Output time ',output_idx,'/',output_num
            output_idx = output_idx+1
            call write_netCDF(pu,pv,ph,output_num,output_idx)
        end if

    end do

	close(21)
!------------------------- Main part of the program end------------
 	print *,'The main part of this program has ended'
!-------------------------- Out put the U,V,H ---------------------
	print *,'Now,the working is to output result'

    call SYSTEM_CLOCK(time_2)
    print*,'It took ',dble(time_2-time_1)/1000.0,' seconds to run this program'
    
    END PROGRAM MAIN
    
    
FUNCTION inner(u1,v1,h1,u2,v2,h2)
    use mesh
    use module_para, only : nx,ny,nx1,ny1
    
    implicit none
    
    real*8, dimension(0:nx1,0:ny1) :: u1,v1,h1
    real*8, dimension(0:nx1,0:ny1) :: u2,v2,h2
    
    real*8                         :: inner
    integer                        :: i,j
    
    inner = 0.d0
    
    do j=1,ny
        do i=1,nx
            inner=inner+(u1(i,j)*u2(i,j)+v1(i,j)*v2(i,j)+h1(i,j)*h2(i,j))*cos_lat(i,j)
        end do
    end do
    
END	FUNCTION inner
