SUBROUTINE L_operator(wu,wv,wh,du,dv,dh)
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
	real*8,dimension(0:nx1,0:ny1)             :: hcos	 ! hh = h*cos(theta)
	real*8,dimension(0:nx1,0:ny1)             :: f 	     ! ff = 2*omg*sin(theta)+u*tan(theta)/a
    real*8,dimension(0:nx1,0:ny1)             :: adv_u   ! u momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_v   ! v momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_u_x ! u momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_v_x ! v momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_u_y ! u momentum advection
    real*8,dimension(0:nx1,0:ny1)             :: adv_v_y ! v momentum advection
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
                       
            uU(i,j)    = u(i,j)*wu(i,j)
            uV(i,j)    = u(i,j)*wv(i,j)
                       
            hcos(i,j)  = h(i,j)*cos_lat(i,j)
            vcos(i,j)  = v(i,j)*cos_lat(i,j)
                       
            vUcos(i,j) = vcos(i,j)*wu(i,j)
            vVcos(i,j) = vcos(i,j)*wv(i,j)
                       
            f (i,j)    = f1(i,j)+u(i,j)*f2(i,j)
            
            hU   (i,j) = h(i,j)   *wu(i,j)
            hVcos(i,j) = hcos(i,j)*wv(i,j)
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
            ! Compute Advection
            adv_u_x(i,j) = (u   (i,j)*(wu(i+1,j)-wu(i-1,j)) + uU   (i+1,j) - uU   (i-1,j))*c13(i,j)
            adv_u_y(i,j) = (vcos(i,j)*(wu(i,j+1)-wu(i,j-1)) + vUcos(i,j+1) - vUcos(i,j-1))*c14(i,j)
            adv_u  (i,j) = adv_u_x(i,j) + adv_u_y(i,j)
            
            adv_v_x(i,j) = (u   (i,j)*(wv(i+1,j)-wv(i-1,j)) + uV   (i+1,j) - uV   (i-1,j))*c13(i,j)
            adv_v_y(i,j) = (vcos(i,j)*(wv(i,j+1)-wv(i,j-1)) + vVcos(i,j+1) - vVcos(i,j-1))*c14(i,j)
            adv_v  (i,j) = adv_v_x(i,j) + adv_v_y(i,j)
            
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
            du(i,j) = adv_u (i,j) + pgf_u (i,j) - fv(i,j)
            dv(i,j) = adv_v (i,j) + pgf_v (i,j) + fu(i,j)
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
    
    ! Antisymmetry Check
    !adv_u(:,1 ) = 0.d0
    !adv_v(:,1 ) = 0.d0
    !adv_u(:,ny) = 0.d0
    !adv_v(:,ny) = 0.d0
    !print*,'                                                      '
    !print*,'----------------- Antisymmetry Check -----------------'
    !print*,'(ADV_U,U)              = ',sum(adv_u(1:nx,1:ny)*wu(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    !print*,'(ADV_V,V)              = ',sum(adv_v(1:nx,1:ny)*wv(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    !print*,'(PGF_U,U) + (FLUX_X,Z) = ',sum(pgf_u(1:nx,1:ny)*wu(1:nx,1:ny)*cos_lat(1:nx,1:ny)) + sum(flux_x(1:nx,1:ny)*wh(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    !print*,'(PGF_V,V) + (FLUX_Y,Z) = ',sum(pgf_v(1:nx,1:ny)*wv(1:nx,1:ny)*cos_lat(1:nx,1:ny)) + sum(flux_y(1:nx,1:ny)*wh(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    !print*,'(fV   ,U) - (fU    ,V) = ',sum(fv   (1:nx,1:ny)*wu(1:nx,1:ny)*cos_lat(1:nx,1:ny)) - sum(fu    (1:nx,1:ny)*wv(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    !print*,'(PGF_F,F) + (FLUX_Y,F) = ',sum(pgf_u(1:nx,1:ny)*wu(1:nx,1:ny)*cos_lat(1:nx,1:ny)) + sum(pgf_v (1:nx,1:ny)*wv(1:nx,1:ny)*cos_lat(1:nx,1:ny)) + sum(dh   (1:nx,1:ny)*wh(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    !print*,'(LF   ,F)              = ',sum(du   (1:nx,1:ny)*wu(1:nx,1:ny)*cos_lat(1:nx,1:ny)) + sum(dv    (1:nx,1:ny)*wv(1:nx,1:ny)*cos_lat(1:nx,1:ny)) + sum(dh   (1:nx,1:ny)*wh(1:nx,1:ny)*cos_lat(1:nx,1:ny))
    
    !do j = 1,ny
    !    print*,sum(adv_u_x(1:nx,j)*wu(1:nx,j)*cos_lat(1:nx,j))
    !enddo
    
    !do i = 1,nx
    !    print*,sum(adv_u_y(i,1:ny)*wu(i,1:ny)*cos_lat(i,1:ny))
    !enddo
        
    !print*,'max du                 = ',maxval(du   (1:nx,1:ny))
    !print*,'min du                 = ',minval(du   (1:nx,1:ny))
    !print*,'max dv                 = ',maxval(dv   (1:nx,1:ny))
    !print*,'min dv                 = ',minval(dv   (1:nx,1:ny))
    !print*,'max dh                 = ',maxval(dh   (1:nx,1:ny))
    !print*,'min dh                 = ',minval(dh   (1:nx,1:ny))
    !print*,'max u                  = ',maxval(u    (1:nx,1:ny))
    !print*,'min u                  = ',minval(u    (1:nx,1:ny))
    !print*,'max v                  = ',maxval(v    (1:nx,1:ny))
    !print*,'min v                  = ',minval(v    (1:nx,1:ny))
    !print*,'max U                  = ',maxval(wu   (1:nx,1:ny))
    !print*,'min U                  = ',minval(wu   (1:nx,1:ny))
    !print*,'max V                  = ',maxval(wv   (1:nx,1:ny))
    !print*,'min V                  = ',minval(wv   (1:nx,1:ny))
    !print*,'max vcos               = ',maxval(vcos (1:nx,1:ny))
    !print*,'min vcos               = ',minval(vcos (1:nx,1:ny))
    !print*,'max uU                 = ',maxval(uU   (1:nx,1:ny))
    !print*,'min uU                 = ',minval(uU   (1:nx,1:ny))
    !print*,'max uV                 = ',maxval(uV   (1:nx,1:ny))
    !print*,'min uV                 = ',minval(uV   (1:nx,1:ny))
    !print*,'max vUcos              = ',maxval(vUcos(1:nx,1:ny))
    !print*,'min vUcos              = ',minval(vUcos(1:nx,1:ny))
    !print*,'max vVcos              = ',maxval(vVcos(1:nx,1:ny))
    !print*,'min vVcos              = ',minval(vVcos(1:nx,1:ny))
    !print*,'                                                      '
    !pause
    
END SUBROUTINE L_operator
