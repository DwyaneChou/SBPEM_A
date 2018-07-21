module mesh
    use module_para,only: nx,nx1,ny,ny1
    use module_para,only: a,omg0,pi
    implicit none
    
    real*8,dimension(0:nx1,0:ny1)   :: latitude
    real*8,dimension(0:nx1,0:ny1)   :: longitude
    real*8,dimension(0:nx1,0:ny1)   :: cos_lat
    real*8,dimension(0:nx1,0:ny1)   :: sin_lat
    real*8,dimension(0:nx1,0:ny1)   :: tan_lat
                                    
    real*8                          :: dlambda,dtheta
    
    real*8,  dimension(0:nx1,0:ny1) :: c11, c12  ! c11(j)=1/(2*a*cos(theta(j))*dlambda),c12(j)=1/(2*a*cos(theta(j))*dtheta)
    real*8,  dimension(0:nx1,0:ny1) :: c13, c14  ! c13(j)=c11(j)/2, c14(j)=c12(j)/2
    real*8,  dimension(0:nx1,0:ny1) :: f1 , f2   ! f1(j)=2*omg0*sin(theta(j)), f2(j)=tan(theta(j))/a
        
    contains
    subroutine generate_mesh
        implicit none
        
        integer i,j
        
        pi=datan(1d0)*4.d0
        
        dlambda = 2*pi/nx
        dtheta  = pi/(ny-1)
        
        !print*,dlambda,dtheta
        do j=1,ny
            do i=1,nx
                longitude(i,j) = dble(i-1) * dlambda
                latitude (i,j) = -0.5d0*pi + dble(j-1) * dtheta
                
                cos_lat  (i,j) = dcos(latitude(i,j))
                sin_lat  (i,j) = dsin(latitude(i,j))
                tan_lat  (i,j) = dtan(latitude(i,j))
            enddo
        enddo
        
        longitude(0  ,:) = longitude(nx,:)
        longitude(nx1,:) = longitude(1 ,:)
        latitude (:,0  ) = -pi/2.d0
        latitude (:,ny1) =  pi/2.d0
        
        sin_lat  (:,1  ) = -1.d0
        sin_lat  (:,ny ) =  1.d0
        
        ! Boundary on pole
        cos_lat  (:,1  ) = cos(-pi/2.d0+0.5*dtheta)/4.d0
        cos_lat  (:,ny ) = cos( pi/2.d0-0.5*dtheta)/4.d0
        
        ! Compute parameters
        c11 = 1.d0/(2.d0*a*cos_lat*dlambda)
        c12 = 1.d0/(2.d0*a*cos_lat*dtheta )
        c13 = 0.5d0*c11
        c14 = 0.5d0*c12
        f1  = 2.d0*omg0*sin_lat
        f2  = tan_lat/a
        
        print*,'max longitude = ',maxval(longitude(1:nx,1:ny))/pi*180.d0
        print*,'min longitude = ',minval(longitude(1:nx,1:ny))/pi*180.d0
        print*,'max latitude  = ',maxval(latitude (1:nx,1:ny))/pi*180.d0
        print*,'min latitude  = ',minval(latitude (1:nx,1:ny))/pi*180.d0
        
    end subroutine generate_mesh
    
end module mesh
    