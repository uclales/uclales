!
! Merge raw particle data from UCLA-LES
! Command line input: 'prefix var nxprocs nyprocs'
! e.g: './merge_particles dir/rico x 8 8'
! Only one variable can be provided, allows for poor man's 
! parallelization using a number of serial batch jobs
!

program mergeparticles
  use netcdf
  implicit none
  
  character(100)                    :: pref,cnx,cny,path,var
  integer                           :: nx,ny
  integer                           :: nargs
  character(100)                    :: varname,path_out
  integer                           :: ncid_in,varid_in
  integer                           :: ncid_out,pdim_out,tdim_out,pid_out,tid_out,varid_out,start_in(2),start_out(2),count(2)
  integer                           :: npin=0,nptot=0,ntime=0,ntin=0
  integer                           :: i,j,t,loc
  character(8)                      :: procxy
  real, allocatable, dimension(:)   :: var_read
  
  real                              :: etime
  real*4                            :: elapsed(2)
  real                              :: total
  
  !-- MAIN CALLS
  ! Get input user 
  nargs = command_argument_count()
  if(nargs .lt. 4) then 
    print*,'provide -prefix varflg nxprox nyprocs- as arguments (see source code), STOPPING'
    stop
  else
    call get_command_argument(1,pref)
    call get_command_argument(2,var)
    call get_command_argument(3,cnx)
    call get_command_argument(4,cny)
    read (cnx,*) nx
    read (cny,*) ny
  end if
  
  ! Loop through files to read number of times/particles 
  do i=1,nx
    do j=1,ny
      write(procxy,'(i4.4,i4.4)') i-1,j-1
      path = trim(pref)//".particles."//procxy//".nc"
      call check(nf90_open(path,nf90_nowrite,ncid_in))
      call check(nf90_inq_dimid(ncid_in, "particles", varid_in))
      call check(nf90_inquire_dimension(ncid_in, varid_in, varname, npin))
      call check(nf90_inq_dimid(ncid_in, "time", varid_in))
      call check(nf90_inquire_dimension(ncid_in, varid_in, varname, ntin))
      call check(nf90_close(ncid_in))
      if(ntin > ntime) ntime = ntin
      nptot = nptot + npin 
    end do
  end do
  print*,'found ', nptot, 'particles, ', ntime, 'time steps'
 
  ! Create new NetCDF file
  print*,'creating new NetCDF file' 
  path_out = trim(pref)//".particles."//trim(var)//".nc"
  call check(nf90_create(path_out,nf90_clobber,ncid_out))
  ! Create dimensions
  call check(nf90_def_dim(ncid_out,'particles',nptot,pdim_out))
  call check(nf90_def_dim(ncid_out,'time',ntime,tdim_out))
  ! Create variables
  !call check(nf90_def_var(ncid_out,'particles',nf90_double,(/pdim_out/),pid_out))
  call check(nf90_def_var(ncid_out,'time',nf90_double,(/tdim_out/),tid_out))
  call check(nf90_def_var(ncid_out,trim(var),nf90_float,(/pdim_out,tdim_out/),varid_out))
  ! Leave define mode
  call check(nf90_enddef(ncid_out))

  ! Copy data
  print *,'merging files' 
  loc = 1
  do i=1,nx
    do j=1,ny
      write(procxy,'(i4.4,i4.4)') i-1,j-1
      path = trim(pref)//".particles."//procxy//".nc"
      call check(nf90_open(path,nf90_nowrite,ncid_in))
      call check(nf90_inq_dimid(ncid_in, "particles", varid_in))
      call check(nf90_inquire_dimension(ncid_in, varid_in, varname, npin))
  
      allocate(var_read(npin))
  
      call check(nf90_inq_varid(ncid_in, trim(var), varid_in))
      count  = (/npin,1/)
      do t=1,ntin
        start_in  = (/1,t/)
        start_out = (/loc,t/) 
        call check(nf90_get_var(ncid_in,  varid_in,  var_read, start=start_in, count=count))
        call check(nf90_put_var(ncid_out, varid_out, var_read, start=start_out, count=count)) 
      end do
  
      loc = loc + npin
      deallocate(var_read)

      print*,i,j

      if(i .eq. nx .and. j .eq. ny) then
        allocate(var_read(nptot))
        call check(nf90_inq_varid(ncid_in, 'time', varid_in))
        call check(nf90_get_var(ncid_in,  varid_in, var_read, start=(/1/), count=(/ntime/)))
        call check(nf90_put_var(ncid_out, tid_out, var_read, start=(/1/), count=(/ntime/))) 
        deallocate(var_read)
      end if

      call check(nf90_close(ncid_in))

      print*,'Processed ',trim(path),(i-1)*nx+j,'/',nx*ny 
    end do
  end do
  
  call check(nf90_close(ncid_out))

contains

  subroutine check(status)
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program mergeparticles
