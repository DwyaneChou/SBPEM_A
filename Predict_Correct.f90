
    subroutine Predict_Correct(pass_type,dt,dtn)
        use module_para ,only : nx,nx1,ny,ny1
        use module_array,only : wu,wv,wh
        
        implicit none
        integer,intent(in )           :: pass_type
        real*8 ,intent(in )           :: dt
        real*8 ,intent(out)           :: dtn
        
        real*8,dimension(0:nx1,0:ny1) :: tmp_u,tmp_v,tmp_h
        real*8,dimension(0:nx1,0:ny1) :: LU0,LV0,LZ0,&
                                         LU1,LV1,LZ1,&
                                         LU2,LV2,LZ2
        
        real*8                        :: LF1LF2,LF2_norm2
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
                
        !Predict
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,LU0,LV0,LZ0)
        
        do j=1,ny
            do i=1,nx
                tmp_u(i,j) = wu(i,j) - dt2*LU0(i,j)
                tmp_v(i,j) = wv(i,j) - dt2*LV0(i,j)
                tmp_h(i,j) = wh(i,j) - dt2*LZ0(i,j)
            enddo
        enddo
        
        !Correct 1
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,LU1,LV1,LZ1)
        
        do j=1,ny
            do i=1,nx
                tmp_u(i,j) = wu(i,j) - dt2*LU1(i,j)
                tmp_v(i,j) = wv(i,j) - dt2*LV1(i,j)
                tmp_h(i,j) = wh(i,j) - dt2*LZ1(i,j)
            enddo
        enddo
        
        !Correct 2
        call spatial_discrete(pass_type,tmp_u,tmp_v,tmp_h,LU2,LV2,LZ2)
        
        LF1LF2    = inner(LU1,LV1,LZ1,LU2,LV2,LZ2)
        LF2_norm2 = inner(LU2,LV2,LZ2,LU2,LV2,LZ2)
        
        do j=1,ny
            do i=1,nx
                beta = LF1LF2/LF2_norm2
            enddo
        enddo
        dtn = dt*beta
        
        do j=1,ny
            do i=1,nx
                wu(i,j) = wu(i,j) - dtn*LU2(i,j)
                wv(i,j) = wv(i,j) - dtn*LV2(i,j)
                wh(i,j) = wh(i,j) - dtn*LZ2(i,j)
            enddo
        enddo
        
    end subroutine Predict_Correct
    