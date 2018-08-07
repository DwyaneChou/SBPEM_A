    subroutine CSP2(dt,dtn)
        ! 2nd order conservative split pattern
        use module_para ,only: fast_pass,slow_pass,&
                               split_step,split_time_step
        use module_array,only: wu,wv,wh
        implicit none
        real*8 ,intent(in ) :: dt
        real*8 ,intent(out) :: dtn
        integer             :: pass_type
        
        real*8              :: dt2
        
        integer             :: i
        
        dt2 = 0.5*dt
        
        call integrator(slow_pass,dt2,dtn)
        
        do i=1,split_step
            call integrator(fast_pass,split_time_step ,dtn)
        enddo
        
        call integrator(slow_pass,dt2,dtn)
        
        dtn = dtn*2.d0
        
    end subroutine CSP2
    