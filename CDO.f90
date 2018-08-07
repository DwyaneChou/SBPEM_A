SUBROUTINE CDO(pass_type,dt,dtn)
!Consisitent Dissipation Operator Scheme
    use module_para, only : nx,nx1,ny,ny1
    use module_array,only : wu,wv,wh
!
	implicit none
    integer,intent(in )            :: pass_type
    real*8 ,intent(in )            :: dt        ! time stepsize
    real*8 ,intent(out)            :: dtn	    
	real*8 ,dimension(0:nx1,0:ny1) :: du,dv,dh 	! Tendencies of wu, wv and wh, respectively.
	real*8 ,dimension(0:nx1,0:ny1) :: bu,bv,bh 	! Dissipation of wu, wv and wh, repsectively.
	real*8                         :: eps		! dissipation coefficient
	real*8                         :: dt2		! work variable
	real*8                         :: k1,k2,k3	! work variables
	real*8,dimension(0:nx1,0:ny1)  :: tu,tv,th 	! work arrays
    
	real*8,external                :: inner
    
	integer                        :: i,j		! working variables
    
    eps=0.5d0
    
	do j=1,ny
	   do i=1,nx
          tu(i,j)=wu(i,j)
          tv(i,j)=wv(i,j)
          th(i,j)=wh(i,j)
	   end do
    end do

    ! consistent dissipation operator
    call B_operator(pass_type,tu,tv,th,du,dv,dh,bu,bv,bh,dt)
    
	k3  =  inner(bu,bv,bh,tu,tv,th)
          
	k1  = -inner(du,dv,dh,du,dv,dh)/k3
	k2  =  inner(bu,bv,bh,du,dv,dh)/k3
	k3  = -inner(bu,bv,bh,bu,bv,bh)/k3
    
    dtn = 2.d0-k1/eps
    dtn = dtn/(k2+dsqrt(k2*k2+eps*dtn*k3))
    dt2 = eps*dtn*dtn
    
	do j=1,ny
	    do i=1,nx
           wu(i,j)=wu(i,j)-dtn*du(i,j)+dt2*bu(i,j)
           wv(i,j)=wv(i,j)-dtn*dv(i,j)+dt2*bv(i,j)
           wh(i,j)=wh(i,j)-dtn*dh(i,j)+dt2*bh(i,j)
	    end do
    end do
    
END SUBROUTINE CDO