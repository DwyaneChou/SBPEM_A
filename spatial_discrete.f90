    subroutine spatial_discrete(pass_type,wu,wv,wh,du,dv,dh)
        ! 2nd order conservative split pattern
        use module_para ,only:nx1,ny1,full_pass,fast_pass,slow_pass
        implicit none
        integer                      ,intent(in ) :: pass_type
        real*8,dimension(0:nx1,0:ny1),intent(in ) :: wu
        real*8,dimension(0:nx1,0:ny1),intent(in ) :: wv
        real*8,dimension(0:nx1,0:ny1),intent(in ) :: wh
        real*8,dimension(0:nx1,0:ny1),intent(out) :: du
        real*8,dimension(0:nx1,0:ny1),intent(out) :: dv
        real*8,dimension(0:nx1,0:ny1),intent(out) :: dh

        
        if    (pass_type == full_pass)then
            call L_operator   (wu,wv,wh,du,dv,dh)
        elseif(pass_type == fast_pass)then
            call fast_operator(wu,wv,wh,du,dv,dh)
        elseif(pass_type == slow_pass)then
            call slow_operator(wu,wv,wh,du,dv,dh)
        endif
        
    end subroutine spatial_discrete
    