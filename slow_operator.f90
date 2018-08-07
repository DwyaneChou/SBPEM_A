SUBROUTINE slow_operator(wu,wv,wh,du,dv,dh)
!
! This subroutine is to calculate the tendency of geopotential height
!
    use module_para      ,only : nx,nx1,ny,ny1,a
    use mesh             ,only : cos_lat,f1,f2,c11,c12,c13,c14,dtheta
!
 	implicit none
!
	real*8,dimension(0:nx1,0:ny1)             :: wu		 ! wu = U = h*u
	real*8,dimension(0:nx1,0:ny1)             :: wv		 ! wv = V =	h*v 
	real*8,dimension(0:nx1,0:ny1)             :: wh		 ! Geopotential height
	real*8,dimension(0:nx1,0:ny1),intent(out) :: du		 ! Tendency of wu
	real*8,dimension(0:nx1,0:ny1),intent(out) :: dv		 ! Tendency of wv
	real*8,dimension(0:nx1,0:ny1),intent(out) :: dh		 ! Tendency of wh
	real*8,dimension(0:nx1,0:ny1)             :: u		 ! zonal wind
	real*8,dimension(0:nx1,0:ny1)             :: v		 ! meridional wind
	real*8,dimension(0:nx1,0:ny1)             :: h		 ! velocity of gravity wave: h  = sqrt(wh)
    real*8,dimension(0:nx1,0:ny1)             :: uU		 ! u*U
    real*8,dimension(0:nx1,0:ny1)             :: uV		 ! u*V
    real*8,dimension(0:nx1,0:ny1)             :: vUcos	 ! v*U*cos(theta)
    real*8,dimension(0:nx1,0:ny1)             :: vVcos	 ! v*V*cos(theta)
	real*8,dimension(0:nx1,0:ny1)             :: vcos	 ! vv = v*cos(theta)
    real*8,dimension(0:nx1,0:ny1)             :: adv_u   ! u momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_v   ! v momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_u_x ! u momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_v_x ! v momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_u_y ! u momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_v_y ! v momentum advection
    
	integer                                   :: i,j     ! working variables
    
    do j=1,ny
        do i=1,nx
            !inverse IAP transformation
            h(i,j)     = dsqrt(wh(i,j))
            u(i,j)     = wu(i,j)/h(i,j)
            v(i,j)     = wv(i,j)/h(i,j)
                       
            uU(i,j)    = u(i,j)*wu(i,j)
            uV(i,j)    = u(i,j)*wv(i,j)
                       
            vcos(i,j)  = v(i,j)*cos_lat(i,j)
                       
            vUcos(i,j) = vcos(i,j)*wu(i,j)
            vVcos(i,j) = vcos(i,j)*wv(i,j)
        enddo
    enddo
    
    ! zonal boundary condition
    do j = 1,ny
        wu(nx1,j) = wu(1 ,j)
        wu(0  ,j) = wu(nx,j)
        
        uU(nx1,j) = uU(1 ,j)
        uU(0  ,j) = uU(nx,j)
        
        wv(nx1,j) = wv(1 ,j)
        wv(0  ,j) = wv(nx,j)
        
        uV(nx1,j) = uV(1 ,j)
        uV(0  ,j) = uV(nx,j)
        
        wh(nx1,j) = wh(1 ,j)
        wh(0  ,j) = wh(nx,j)
    enddo
    
    do j=1,ny
        do i=1,nx
            ! Compute Advection
            adv_u_x(i,j) = (u   (i,j)*(wu(i+1,j)-wu(i-1,j)) + uU   (i+1,j) - uU   (i-1,j))*c13(i,j)
            adv_u_y(i,j) = (vcos(i,j)*(wu(i,j+1)-wu(i,j-1)) + vUcos(i,j+1) - vUcos(i,j-1))*c14(i,j)
            adv_u  (i,j) = adv_u_x(i,j) + adv_u_y(i,j)
            
            adv_v_x(i,j) = (u   (i,j)*(wv(i+1,j)-wv(i-1,j)) + uV   (i+1,j) - uV   (i-1,j))*c13(i,j)
            adv_v_y(i,j) = (vcos(i,j)*(wv(i,j+1)-wv(i,j-1)) + vVcos(i,j+1) - vVcos(i,j-1))*c14(i,j)
            adv_v  (i,j) = adv_v_x(i,j) + adv_v_y(i,j)
        enddo
    enddo
        
    do j=1,ny
        do i=1,nx
            du(i,j) = adv_u (i,j)
            dv(i,j) = adv_v (i,j)
            dh(i,j) = 0.d0
        enddo
    enddo
    
    
END SUBROUTINE slow_operator
