
    subroutine RK4(pass_type,dt,dtn)
        use module_para ,only : nx,nx1,ny,ny1
        use module_array,only : wu,wv,wh
        
        implicit none
        integer,intent(in )           :: pass_type
        real*8 ,intent(in )           :: dt
        real*8 ,intent(out)           :: dtn
        
        real*8,dimension(0:nx1,0:ny1) :: tmp_u,tmp_v,tmp_h
        real*8,dimension(0:nx1,0:ny1) :: phi_u,phi_v,phi_h
        real*8,dimension(0:nx1,0:ny1) :: RK1_u,RK1_v,RK1_h
        real*8,dimension(0:nx1,0:ny1) :: RK2_u,RK2_v,RK2_h
        real*8,dimension(0:nx1,0:ny1) :: RK3_u,RK3_v,RK3_h
        real*8,dimension(0:nx1,0:ny1) :: RK4_u,RK4_v,RK4_h
        
        real*8                        :: phi_norm2,R1R2,R2R3,R3R4
        real*8                        :: beta
        
        real*8,external               :: inner
        real*8                        :: dt2
        
        integer i,j
        
        dt2 = 0.5d0*dt
        
        do j=1,ny
            do i=1,nx
                tmp_u(i,j) = wu(i,j)
                tmp_v(i,j) = wv(i,j)
                tmp_h(i,j) = wh(i,j)
            enddo
        enddo
                
        ! Compute RK1
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,RK1_u,RK1_v,RK1_h)
        
        ! Compute RK2
        do j=1,ny
            do i=1,nx
                tmp_u(i,j) = wu(i,j) - dt2*RK1_u(i,j)
                tmp_v(i,j) = wv(i,j) - dt2*RK1_v(i,j)
                tmp_h(i,j) = wh(i,j) - dt2*RK1_h(i,j)
            enddo
        enddo
        
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,RK2_u,RK2_v,RK2_h)
        
        ! Compute RK3
        do j=1,ny
            do i=1,nx
                tmp_u(i,j) = wu(i,j) - dt2*RK2_u(i,j)
                tmp_v(i,j) = wv(i,j) - dt2*RK2_v(i,j)
                tmp_h(i,j) = wh(i,j) - dt2*RK2_h(i,j)
            enddo
        enddo
        
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,RK3_u,RK3_v,RK3_h)
        
        ! Compute RK4
        do j=1,ny
            do i=1,nx
                tmp_u(i,j) = wu(i,j) - dt*RK3_u(i,j)
                tmp_v(i,j) = wv(i,j) - dt*RK3_v(i,j)
                tmp_h(i,j) = wh(i,j) - dt*RK3_h(i,j)
            enddo
        enddo
        
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,RK4_u,RK4_v,RK4_h)
        
        ! Compute dtn
        phi_u = (RK1_u + 2.d0*RK2_u + 2.d0*RK3_u + RK4_u)/6.d0
        phi_v = (RK1_v + 2.d0*RK2_v + 2.d0*RK3_v + RK4_v)/6.d0
        phi_h = (RK1_h + 2.d0*RK2_h + 2.d0*RK3_h + RK4_h)/6.d0
        
        phi_norm2 = inner(phi_u,phi_v,phi_h,phi_u,phi_v,phi_h)
        R1R2      = inner(RK1_u,RK1_v,RK1_h,RK2_u,RK2_v,RK2_h)
        R2R3      = inner(RK2_u,RK2_v,RK2_h,RK3_u,RK3_v,RK3_h)
        R3R4      = inner(RK3_u,RK3_v,RK3_h,RK4_u,RK4_v,RK4_h)
        
        beta      = 1.d0/(3.d0*phi_norm2)*(R1R2+R2R3+R3R4)
        dtn       = beta*dt
        
        do j=1,ny
            do i=1,nx
                wu(i,j) = wu(i,j) - dtn*phi_u(i,j)
                wv(i,j) = wv(i,j) - dtn*phi_v(i,j)
                wh(i,j) = wh(i,j) - dtn*phi_h(i,j)
            enddo
        enddo
                
    end subroutine RK4
    