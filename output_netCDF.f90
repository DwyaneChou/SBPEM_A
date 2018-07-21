MODULE output_netCDF
    use netcdf
    use module_para
    use mesh
    implicit none
    include 'netcdf.inc'
    
    contains
    
    subroutine write_netCDF(u,v,z,output_num,output_idx)
        implicit none
        integer                 ,intent(in):: output_num,output_idx
        real*4 ,dimension(nx,ny),intent(in):: u
        real*4 ,dimension(nx,ny),intent(in):: v
        real*4 ,dimension(nx,ny),intent(in):: z
                
        real          lon(nx),lat(ny)
        integer       ncid,status
        integer       longitude_dimid,latitude_dimid,output_num_dimid
        integer       u_id,v_id,z_id,longitude_id,latitude_id
        
        character*9:: nc_file='output.nc'
        if(output_idx==1)then
            
            !print*,'nf90_create'
            status = nf90_create(path=nc_file,cmode=NF90_CLOBBER,ncid=ncid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_def_dim'
            status = nf90_def_dim(ncid,'longitude',nx             ,longitude_dimid)
            status = nf90_def_dim(ncid,'latitude' ,ny             ,latitude_dimid)
            status = nf90_def_dim(ncid,'nt'       ,NF90_UNLIMITED ,output_num_dimid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_def_var'
            status = nf90_def_var(ncid,'longitude',NF90_FLOAT,(/longitude_dimid,latitude_dimid                 /),longitude_id)
            status = nf90_def_var(ncid,'latitude' ,NF90_FLOAT,(/longitude_dimid,latitude_dimid                 /),latitude_id )
            status = nf90_def_var(ncid,'u'        ,NF90_FLOAT,(/longitude_dimid,latitude_dimid,output_num_dimid/),u_id        )
            status = nf90_def_var(ncid,'v'        ,NF90_FLOAT,(/longitude_dimid,latitude_dimid,output_num_dimid/),v_id        )
            status = nf90_def_var(ncid,'z'        ,NF90_FLOAT,(/longitude_dimid,latitude_dimid,output_num_dimid/),z_id        )
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_att'
            !status = nf90_put_att()
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_enddef'
            status = nf90_enddef(ncid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_var'
            status = nf90_put_var(ncid,longitude_id,longitude                                  )
            status = nf90_put_var(ncid,latitude_id ,latitude                                   )
            status = nf90_put_var(ncid,u_id        ,u        ,start=(/1,1,1/),count=(/nx,ny,1/))
            status = nf90_put_var(ncid,v_id        ,v        ,start=(/1,1,1/),count=(/nx,ny,1/))
            status = nf90_put_var(ncid,z_id        ,z        ,start=(/1,1,1/),count=(/nx,ny,1/))
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_close'
            status = nf90_close(ncid)
            !if(status/=nf90_noerr) call handle_err(status)
        else
            !print*,'nf90_open'
            status = nf90_open(nc_file,NF90_WRITE,ncid)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_inq_varid'
            status = nf90_inq_varid(ncid,'u',u_id)
            status = nf90_inq_varid(ncid,'v',v_id)
            status = nf90_inq_varid(ncid,'z',z_id)
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_var'
            status = nf90_put_var(ncid,u_id,u,start=(/1,1,output_idx/),count=(/nx,ny,1/))
            status = nf90_put_var(ncid,v_id,v,start=(/1,1,output_idx/),count=(/nx,ny,1/))
            status = nf90_put_var(ncid,z_id,z,start=(/1,1,output_idx/),count=(/nx,ny,1/))
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_put_att'
            !status = nf90_put_att
            !if(status/=nf90_noerr) call handle_err(status)
            
            !print*,'nf90_close'
            status = nf90_close(ncid)
            !if(status/=nf90_noerr) call handle_err(status)
        endif
        
    end subroutine write_netCDF
        
    subroutine handle_err(status)
        implicit none
        integer,intent(in)::status
        
        if(status/=nf90_noerr)then
            print*, trim(nf90_strerror(status))
            stop "Stopped by netCDF"
        endif
                    
    endsubroutine handle_err
    
END MODULE output_netCDF