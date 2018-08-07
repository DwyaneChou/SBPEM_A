SUBROUTINE fast_operator(wu,wv,wh,du,dv,dh)
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
	real*8,dimension(0:nx1,0:ny1)             :: hcos	 ! hh = h*cos(theta)
	real*8,dimension(0:nx1,0:ny1)             :: f 	     ! ff = 2*omg*sin(theta)+u*tan(theta)/a
    real*8,dimension(0:nx1,0:ny1)             :: pgf_u   ! pressure gradient force on u
    real*8,dimension(0:nx1,0:ny1)             :: pgf_v   ! pressure gradient force on v
    real*8,dimension(0:nx1,0:ny1)             :: fv      ! f*v : colioris force and curvature term on u
    real*8,dimension(0:nx1,0:ny1)             :: fu      ! f*u : colioris force and curvature term on v
    real*8,dimension(0:nx1,0:ny1)             :: flux_x  ! zonal mass flux
    real*8,dimension(0:nx1,0:ny1)             :: flux_y  ! meridional mass flux
    real*8,dimension(0:nx1,0:ny1)             :: hU      ! h*U
    real*8,dimension(0:nx1,0:ny1)             :: hVcos   ! h*V*cos(theta)
    
	integer                                   :: i,j     ! working variables
    
    do j=1,ny
        do i=1,nx
            !inverse IAP transformation
            h(i,j)     = dsqrt(wh(i,j))
            u(i,j)     = wu(i,j)/h(i,j)
            v(i,j)     = wv(i,j)/h(i,j)
                       
            hcos(i,j)  = h(i,j)*cos_lat(i,j)
                       
            f (i,j)    = f1(i,j)+u(i,j)*f2(i,j)
            
            hU   (i,j) = h(i,j)   *wu(i,j)
            hVcos(i,j) = hcos(i,j)*wv(i,j)
        enddo
    enddo
    
    ! zonal boundary condition
    do j = 1,ny
        wu(nx1,j) = wu(1 ,j)
        wu(0  ,j) = wu(nx,j)
        
        wv(nx1,j) = wv(1 ,j)
        wv(0  ,j) = wv(nx,j)
        
        wh(nx1,j) = wh(1 ,j)
        wh(0  ,j) = wh(nx,j)
        
        hU(nx1,j) = hU(1 ,j)
        hU(0  ,j) = hU(nx,j)
    enddo
        
    ! meridional boundary conditions
    do i = 1,nx
        hVcos(i,1 ) = 0.d0
        hVcos(i,ny) = 0.d0
    enddo
    
    do j=1,ny
        do i=1,nx            
            ! Compute Pressure Gradient Force
            pgf_u (i,j) = h(i,j)*(wh(i+1,j)-wh(i-1,j))*c11(i,j)
            pgf_v (i,j) = h(i,j)*(wh(i,j+1)-wh(i,j-1))/(2.d0*a*dtheta)
            
            ! Compute Colioris and Curvature term
            fv    (i,j) = f(i,j)*wv(i,j)
            fu    (i,j) = f(i,j)*wu(i,j)
            
            ! Compute mass flux
            flux_x(i,j) = (hU   (i+1,j) - hU   (i-1,j))*c11(i,j)
            flux_y(i,j) = (hVcos(i,j+1) - hVcos(i,j-1))*c12(i,j)
        enddo
    enddo
        
    do j=1,ny
        do i=1,nx
            du(i,j) = pgf_u (i,j) - fv(i,j)
            dv(i,j) = pgf_v (i,j) + fu(i,j)
            dh(i,j) = flux_x(i,j) + flux_y(i,j)
        enddo
    enddo
        
    ! For Southern Pole
    j = 1
    du(:,j) = 0.d0
    dv(:,j) = 0.d0
    dh(:,j) =  sum(hVcos(1:nx,j+1))/nx*c12(:,j)
    
    ! For Northern Pole
    j = ny
    du(:,j) = 0.d0
    dv(:,j) = 0.d0
    dh(:,j) = -sum(hVcos(1:nx,j-1))/nx*c12(:,j)
        
END SUBROUTINE fast_operator
