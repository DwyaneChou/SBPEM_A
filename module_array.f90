MODULE module_array

    use module_para

    implicit none

    real*8,dimension(0:nx1,0:ny1) :: u               ! zonal wind
    real*8,dimension(0:nx1,0:ny1) :: v               ! meridional wind
    real*8,dimension(0:nx1,0:ny1) :: wh              ! geopotential height
                                  
    real*8,dimension(0:nx1,0:ny1) :: wu              ! u*sqrt(wh)
    real*8,dimension(0:nx1,0:ny1) :: wv              ! v*sqrt(wh)
    real*8,dimension(0:nx1,0:ny1) :: h               ! sqrt(wh)
END MODULE module_array