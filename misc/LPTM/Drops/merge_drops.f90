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
  
  character(100)                       :: pref,cnx,cny,path,var
  integer                              :: nx,ny
  integer                              :: nargs
  character(100)                       :: varname,path_out
  character(100)                       :: lname,uname
  integer                              :: ncid_in,varid_in
  integer                              :: ncid_out,pdim_out,tdim_out,varid_out,start_in(2),start_out(2)
  integer                              :: count_in(2),count_out(2)
  integer                              :: tid_in,pid_in,pid_out,tid_out
  integer                              :: nptot=0,ntime=0,ntin=0,npin=0,npp=100000
  integer                              :: i,j,p,loc,k,npre
  character(8)                         :: procxy
  real, allocatable, dimension(:,:)    :: var_read, var_read_in, var_read_out, var_red
  logical, allocatable, dimension(:)   :: mask1d
  logical, allocatable, dimension(:,:) :: mask2d
    
  real                                 :: etime,fvalr
  real*4                               :: elapsed(2),fval
  real                                 :: total
 

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
  call check(nf90_def_dim(ncid_out,'time',ntime,tdim_out))
  call check(nf90_def_dim(ncid_out,'particles',nf90_unlimited,pdim_out))
  ! Create variables
  call check(nf90_def_var(ncid_out,'time',nf90_double,(/tdim_out/),tid_out))
  call check(nf90_def_var(ncid_out,trim(var),nf90_float,(/tdim_out,pdim_out/),varid_out))
  !! Leave define mode
  !call check(nf90_enddef(ncid_out))

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
  
  
      call check(nf90_inq_varid(ncid_in, trim(var), varid_in))
    
      if(i .eq. 1 .and. j .eq. 1) then
        ! Put attributes
	! dimensions
	call check(nf90_inq_dimid(ncid_in, "time", tid_in))
	call check(nf90_get_att(ncid_in,tid_in,'longname',lname))
	call check(nf90_put_att(ncid_out,tid_out,'longname',trim(lname)))
	call check(nf90_get_att(ncid_in,tid_in,'units',uname))
        call check(nf90_put_att(ncid_out,tid_out,'units',trim(uname)))
	call check(nf90_get_att(ncid_in,tid_in,'_FillValue',fvalr))
	call check(nf90_put_att(ncid_out,tid_out,'_FillValue',fvalr))
	!call check(nf90_inq_dimid(ncid_in, "particles", pid_in))
	!call check(nf90_get_att(ncid_in,pid_in,'longname',lname))
	!call check(nf90_put_att(ncid_out,pid_out,'longname',trim(lname)))
	!call check(nf90_get_att(ncid_in,pid_in,'units',uname))
        !call check(nf90_put_att(ncid_out,pid_out,'units',trim(uname)))
	!call check(nf90_get_att(ncid_in,pid_in,'_FillValue',fvalr))
	!call check(nf90_put_att(ncid_out,pid_out,'_FillValue',fvalr))
        ! variable
	call check(nf90_get_att(ncid_in,varid_in,'longname',lname))
	call check(nf90_put_att(ncid_out,varid_out,'longname',trim(lname)))
	call check(nf90_get_att(ncid_in,varid_in,'units',uname))
        call check(nf90_put_att(ncid_out,varid_out,'units',trim(uname)))
	call check(nf90_get_att(ncid_in,varid_in,'_FillValue',fval))
	call check(nf90_put_att(ncid_out,varid_out,'_FillValue',fval))
        ! Leave define mode
        call check(nf90_enddef(ncid_out))
      end if
     
      do p=1,npin/npp+1
	!prepare to read data
        start_in  = (/(p-1)*npp+1,1/)
        if (p*npp<npin) then
          count_in  = (/npp,ntin/)
	  allocate(var_read_in(npp,ntin))
	  allocate(mask1d(npp))
	  allocate(mask2d(npp,ntin))
 	else
          count_in  = (/npin-((p-1)*npp),ntin/)
          allocate(var_read_in(npin-((p-1)*npp),ntin))
          allocate(mask1d(npin-((p-1)*npp)))
          allocate(mask2d(npin-((p-1)*npp),ntin))
	end if
        
	!read data
	call check(nf90_get_var(ncid_in,  varid_in,  var_read_in,  start=start_in,  count=count_in))
	
	!reduce dims of var_read_in by deleting the empty particles
        mask1d = .not.all(var_read_in.eq.fval,2)
	npre = count(mask1d)
	print*,'number of active drops',npre
	mask2d = spread(mask1d,2,ntin)
	deallocate(mask1d)
	allocate(var_red(npre,ntin))
	allocate(var_read_out(ntin,npre))
	var_red = reshape(pack(var_read_in,mask2d),(/npre,ntin/))
	deallocate(mask2d)
	deallocate(var_read_in)
	var_read_out = transpose(var_red)
	deallocate(var_red)
	
	
	!prepare to write data
        start_out = (/1,loc/) 	
        count_out = (/ntin,npre/)
	
        !write data
        call check(nf90_put_var(ncid_out, varid_out, var_read_out, start=start_out, count=count_out)) 
	
	deallocate(var_read_out)
	loc = loc + npre
      end do
  
      
      

      print*,i,j

      !get and put time values
      if(i .eq. nx .and. j .eq. ny) then
        allocate(var_read(nptot,1))
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
