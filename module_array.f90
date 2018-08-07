MODULE module_array

    use module_para

    implicit none

    real*8,dimension(:,:),allocatable :: u               ! zonal wind
    real*8,dimension(:,:),allocatable :: v               ! meridional wind
    real*8,dimension(:,:),allocatable :: wh              ! geopotential height
    real*8,dimension(:,:),allocatable :: wu              ! u*sqrt(wh)
    real*8,dimension(:,:),allocatable :: wv              ! v*sqrt(wh)
    real*8,dimension(:,:),allocatable :: h               ! sqrt(wh)
    
    contains
    
    subroutine init_array
        implicit none
        allocate(u (0:nx1,0:ny1),&
                 v (0:nx1,0:ny1),&
                 wh(0:nx1,0:ny1),&
                 wu(0:nx1,0:ny1),&
                 wv(0:nx1,0:ny1),&
                 h (0:nx1,0:ny1))
        
    end subroutine init_array
    
END MODULE module_array