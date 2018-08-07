! ---------------------------------------------------------------------------------------!
! This subroutine is to provide initial condition using four-wave Rossby-Haurwitz waves	 !
! ---------------------------------------------------------------------------------------!
!
SUBROUTINE haurwitz

    use module_para, only : nx,nx1,ny,ny1
    use module_para, only : a,g
    use module_para, only : omg0
    
    use module_array,only : u,v,wh
    
    use mesh
    
    implicit none
    
    real*8,parameter              :: omg  = 7.848d-6         ! angular velocity of RH wave
    real*8,parameter              :: R    = 4d0              ! wave number of RH wave
    real*8,parameter              :: h0   = 8000.d0          ! wave number of RH wave
    
    real*8,dimension(0:nx1,0:ny1) :: u1,u2,u3                ! working array
    real*8,dimension(0:nx1,0:ny1) :: AA1,Ac,A21,A22,A23,Ah   ! working array
    real*8,dimension(0:nx1,0:ny1) :: Bc,BB1,BB2,Bh           ! working array
    real*8,dimension(0:nx1,0:ny1) :: CC,CC1,CC2,Ch           ! working array
    real*8,dimension(0:nx1,0:ny1) :: coslat                  ! working array
    
    integer                       :: i,j                     ! working variable
    
    coslat       = cos_lat
    coslat(:,1 ) = 0.d0
    coslat(:,ny) = 0.d0
    
    do j=1,ny
        do i=1,nx
            u1(i,j) = coslat(i,j)
            u2(i,j) = R*coslat(i,j)**(R-1)*sin_lat(i,j)**2*dcos(R*longitude(i,j))
            u3(i,j) = coslat(i,j)**(R+1)*dcos(R*longitude(i,j))
            u (i,j) = a*omg*(u1(i,j)+u2(i,j)-u3(i,j))
            
            v (i,j) = -a*omg*R*coslat(i,j)**(R-1)*sin_lat(i,j)*dsin(R*longitude(i,j))
            
            AA1 (i,j) = omg*0.5d0*(2.d0*omg0+omg)*coslat(i,j)**2
            Ac  (i,j) = 0.25*omg**2
            A21 (i,j) = (R+1.d0)*coslat(i,j)**(2.d0*R+2.d0)
            A22 (i,j) = (2.d0*R**2-R-2.d0)*coslat(i,j)**(2.d0*R)
            A23 (i,j) = 2.d0*R**2*coslat(i,j)**(2.d0*R-2)
            Ah  (i,j) = AA1(i,j)+Ac(i,j)*(A21(i,j)+A22(i,j)-A23(i,j))
            
            Bc  (i,j) = 2.*(omg0+omg)*omg/((R+1)*(R+2))*coslat(i,j)**R
            BB1 (i,j) = R**2+2.d0*R+2.d0
            BB2 (i,j) = (R+1.d0)**2.*coslat(i,j)**2
            Bh  (i,j) = Bc(i,j)*(BB1(i,j)-BB2(i,j))
            
            CC  (i,j) = 0.25*omg**2*coslat(i,j)**(2.d0*R)
            CC1 (i,j) = (R+1.d0)*coslat(i,j)**2;
            CC2 (i,j) = R+2.d0
            Ch  (i,j) = CC(i,j)*(CC1(i,j)-CC2(i,j))
            
            wh (i,j) = g*h0+a**2*(Ah(i,j) + Bh(i,j)*dcos(R*longitude(i,j)) + Ch(i,j)*dcos(2.0*R*longitude(i,j)))
        enddo
    enddo
    
    u (0  ,:) = u (nx,:)
    v (0  ,:) = v (nx,:)
    wh(0  ,:) = wh(nx,:)
    u (nx1,:) = u (1 ,:)
    v (nx1,:) = v (1 ,:)
    wh(nx1,:) = wh(1 ,:)
    
    u (:,0  ) = 0
    u (:,ny1) = 0
    v (:,0  ) = 0
    v (:,ny1) = 0
    wh(:,0  ) = 0
    wh(:,ny1) = 0
    
END SUBROUTINE haurwitz
