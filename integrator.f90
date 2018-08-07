    subroutine integrator(pass_type,dt,dtn)
        use module_para,only: integral_scheme
        implicit none
        integer,intent(in ) :: pass_type
        integer,intent(in ) :: dt
        integer,intent(out) :: dtn
        
        if    (trim(adjustl(integral_scheme))=='CDO')then
            call CDO(pass_type,dt,dtn)
        elseif(trim(adjustl(integral_scheme))=='RK4')then
            call RK4(pass_type,dt,dtn)
        elseif(trim(adjustl(integral_scheme))=='Predict_Correct')then
            call Predict_Correct(pass_type,dt,dtn)
        endif
        
    end subroutine integrator
    