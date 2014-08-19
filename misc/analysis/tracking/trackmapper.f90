!
! compile with:
! fort -o trackmapperd trackmapper.f90 `/sw/sles10-x64/netcdf-latest-static-intel12/bin/nc-config --fflags --flibs`
!
! command line call 
! trackmapper rico_140.out.xy.track.nc rico_140.out.xy.w_cld.nc nrcloud w_cld smcloudlwp 2
!
! argument list is:
!   1. file name of tracking file
!   2. file name of file containing variable to be mapped
!   3. name of cloud tracking variable
!   4. name of variable to be mapped
!   5. some variable contained in file 1 with the same target array structure as variable 6
!   6. number of coarse graining grid points
!
!



module modnetcdf
  use netcdf
  implicit none
  real(kind=4),    parameter :: fillvalue_r = -1e9
  integer, parameter :: fillvalue_i = -1000000000
  integer(kind=8), parameter :: fillvalue_i64 = -1000000000
  integer(kind=2), parameter :: fillvalue_i16 = -huge(1_2)
  integer :: deflate_level = 1
  type netcdfvar
     integer       :: id
     integer, dimension(:), allocatable :: dim, dimids
     character(80) :: name,longname,units
  end type netcdfvar

  type(netcdfvar) :: xnc, ync, tnc

  interface write_ncvar
    module procedure write_ncvar_1D_r
    module procedure write_ncvar_2D_r
    module procedure write_ncvar_3D_r
    module procedure write_ncvar_1D_i
    module procedure write_ncvar_2D_i
    module procedure write_ncvar_3D_i
    module procedure write_ncvar_1D_i64
    module procedure write_ncvar_2D_i64
    module procedure write_ncvar_3D_i64
  end interface write_ncvar

  interface read_ncvar
    module procedure read_ncvar_1D
    module procedure read_ncvar_2D
    module procedure read_ncvar_3D
  end interface read_ncvar
contains

  subroutine inquire_ncvar(fid, ncvar)
    type(netcdfvar) :: ncvar
    integer :: fid, n, ndims , iret
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_inquire_variable(fid,ncvar%id, ndims=ndims  ))
    if (allocated(ncvar%dim   )) deallocate(ncvar%dim)
    if (allocated(ncvar%dimids)) deallocate(ncvar%dimids)
    allocate(ncvar%dim(ndims))
    allocate(ncvar%dimids(ndims))
    call check ( nf90_inquire_variable(fid,ncvar%id, dimids = ncvar%dimids))
    do n = 1, ndims
      call check (nf90_inquire_dimension(fid, ncvar%dimids(n), len = ncvar%dim(n)))
    end do
    iret = nf90_get_att(fid, ncvar%id, 'longname',ncvar%longname)
    if (iret /= nf90_noerr) then
      ncvar%longname = ''
    end if
    iret = nf90_get_att(fid, ncvar%id, 'units',ncvar%units)
    if (iret /= nf90_noerr) then
      ncvar%units = ''
    end if

  end subroutine inquire_ncvar

  subroutine read_ncvar_1D(fid,ncvar,array,start,count)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer,dimension(1), optional, intent(in) :: start,count
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(start)) then
      call check ( nf90_get_var(fid,ncvar%id,array,start = start, count = count) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_1D

  subroutine read_ncvar_2D(fid,ncvar,array,start,count)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer,dimension(2), optional, intent(in) :: start,count
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(start)) then
      call check ( nf90_get_var(fid,ncvar%id,array,start = start, count = count) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_2D

  subroutine read_ncvar_3D(fid,ncvar,array,start,count)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer,dimension(3), optional, intent(in) :: start,count
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(start)) then
      call check ( nf90_get_var(fid,ncvar%id,array,start = start, count = count) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if
  end subroutine read_ncvar_3D

  subroutine write_ncvar_1D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_1D_r

  subroutine write_ncvar_2D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,real(array)) )
  end subroutine write_ncvar_2D_r

  subroutine write_ncvar_3D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_3D_r

  subroutine write_ncvar_1D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_1D_i

  subroutine write_ncvar_2D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_2D_i

  subroutine write_ncvar_3D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_3D_i

  subroutine write_ncvar_1D_i64(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer(kind=8), dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_1D_i64

  subroutine write_ncvar_2D_i64(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer(kind=8), dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_2D_i64

  subroutine write_ncvar_3D_i64(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer(kind=8), dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_3D_i64

  subroutine define_ncdim(fid,ncvar, dimtype)
    type(netcdfvar)   :: ncvar
    integer           :: fid, iret, xtype
    integer, optional :: dimtype

    if (present(dimtype)) then
      xtype = dimtype
    else
      xtype = nf90_float
    end if
    !..open file and read data
    call check ( nf90_redef(fid))
    iret = nf90_def_dim(fid, trim(ncvar%name), ncvar%dim(1), ncvar%dimids(1))
    if (iret == nf90_enameinuse) then
      call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%dimids(1)) )
    else
      call check(iret)
    end if
    call check ( nf90_enddef(fid))
    call define_ncvar(fid, ncvar, xtype)

  end subroutine define_ncdim

  subroutine define_ncvar(fid,ncvar, xtype)
         INCLUDE 'netcdf.inc'
    type(netcdfvar) :: ncvar
    integer         :: fid, iret, xtype
    !..open file and read data
    call check ( nf90_redef(fid))
    iret =  nf90_def_var(fid, name = trim(ncvar%name),xtype=xtype, dimids = ncvar%dimids, &
                         varid = ncvar%id, deflate_level = deflate_level)
    if (iret == nf90_enameinuse) then
      call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    else
      call check(iret)
      call check ( nf90_put_att(fid, ncvar%id, 'longname', ncvar%longname))
      call check ( nf90_put_att(fid, ncvar%id, 'units', ncvar%units))
      select case(xtype)
      case(nf90_int)
        call check ( nf90_put_att(fid, ncvar%id, '_FillValue', fillvalue_i))
      case(nf90_int64)
        call check (nf_put_att_int(fid, ncvar%id, '_FillValue', nf90_int64, 1, (/fillvalue_i/)))
      case(nf90_float)
        call check ( nf90_put_att(fid, ncvar%id, '_FillValue', fillvalue_r))
      end select
    end if
    call check ( nf90_enddef(fid))

  end subroutine define_ncvar

  subroutine check(status, allowed)
    integer, intent (in) :: status
    integer, dimension(:), intent (in), optional :: allowed
    if(status /= nf90_noerr) then
      write(0,*)  status, trim(nf90_strerror(status))
      if (present(allowed)) then
        if (any(allowed == status)) return
      end if
      stop 2
    end if
  end subroutine check

end module modnetcdf

program trackmapper
  use modnetcdf
  implicit none
  character(100) :: ifile1, ifile2
  integer :: ifid1, ifid2, i, ii, j, jj, t, nx, nxr, ny, nyr, nt, nr, tmax, nrmax, ncoarsegrain
  real, allocatable, dimension(:,:,:) :: var1
  integer, allocatable, dimension(:) :: tt
  real, allocatable, dimension(:) :: time
  character(100) :: ctmp
  real, allocatable, dimension(:,:,:) :: var2
  real, allocatable, dimension(:,:) :: ovar, npts
  type(netcdfvar) :: ncvar1, ncvar2, ncvar3, ncovar, nctime


  !Read command line options
  call get_command_argument(1,ifile1)
  call get_command_argument(2,ifile2)
  call get_command_argument(3,ncvar1%name)
  call get_command_argument(4,ncvar2%name)
  call get_command_argument(5,ncvar3%name) !Something shapelike identical to the output variable
  call get_command_argument(6,ctmp)
  read(ctmp,'(i4)') ncoarsegrain

  !Open files
  write (*,*) ' ifile1 = ',trim(ifile1)
  write (*,*) ' ifile2 = ',trim(ifile2)

  call check ( nf90_open(trim(adjustl(ifile1)), NF90_WRITE, ifid1))
  call check ( nf90_open(trim(adjustl(ifile2)), NF90_NOWRITE, ifid2))

  !Find dimensions/sizes
  write (*,*) ' inquire ncvar1 = ',trim(ncvar1%name)
  call inquire_ncvar(ifid1, ncvar1)
  write (*,*) '   nx = ',ncvar1%dim(1)
  write (*,*) '   ny = ',ncvar1%dim(2)
  write (*,*) '   nt = ',ncvar1%dim(3)

  write (*,*) ' inquire ncvar3 = ',trim(ncvar3%name)
  call inquire_ncvar(ifid1, ncvar3)
  write (*,*) '   tmax  = ',ncvar3%dim(1)
  write (*,*) '   nrmax = ',ncvar3%dim(2)

  write (*,*) ' inquire ncvar2 = ',trim(ncvar2%name)
  call inquire_ncvar(ifid2, ncvar2)
  write (*,*) '   nx = ',ncvar2%dim(1)
  write (*,*) '   ny = ',ncvar2%dim(2)
  write (*,*) '   nt = ',ncvar2%dim(3)

! size of cloud field variable
  nx = ncvar1%dim(1)
  ny = ncvar1%dim(2)
  nt = ncvar1%dim(3)

! size of original data can be larger
  nxr = ncvar2%dim(1)
  nyr = ncvar2%dim(2)

  nctime%name = "time"
  write (*,*) ' inquire time coordinate'
  call inquire_ncvar(ifid1, nctime)
  write (*,*) '   ntime = ',nctime%dim(1)
  allocate(time(nctime%dim(1)))
  write (*,*) ' reading ',trim(nctime%name)
  call read_ncvar(ifid1, nctime, time)
  
  tmax  = ncvar3%dim(1)
  nrmax = ncvar3%dim(2)
  !Allocate and initializecldnrs
  allocate(tt(nrmax))
  tt = fillvalue_i
  allocate(ovar(tmax, nrmax))
  ovar = 0.
  allocate(npts(tmax, nrmax))
  npts = 0
  allocate(var1(nx,ny,nt))
  allocate(var2(nxr,nyr,nt))

  write (*,*) ' reading ',trim(ncvar1%name)
  call read_ncvar(ifid1, ncvar1, var1)
  write (*,*) ' reading ',trim(ncvar2%name)
  write (*,*) ' tstart ',nint(time(1)/(time(2)-time(1)))
  write (*,*) ' time(1) ',time(1)
  !call read_ncvar(ifid2, ncvar2, var2,(/1,1,nint(time(1)/(time(2)-time(1)))/),(/nxr,nyr,nt/))
  call read_ncvar(ifid2, ncvar2, var2,(/1,1,1/),(/nxr,nyr,nt/))
  !call read_ncvar(ifid2, ncvar2, var2)
  write (*,*) ' Done.'  

! Time loop
  do t = 1,nt
    if (mod(t,50).eq.0) write (*,'(a,i5)') '   Time = ',t
    do j=1,ny
      jj = (j-1) * ncoarsegrain + 1
      do i=1,nx
        ii = (i-1) * ncoarsegrain + 1
        if (var1(i,j,t) >= 1) then
          nr = var1(i,j,t)
          if (tt(nr)< 0) then ! A new cell
            tt(nr) = 1
            ovar(tt(nr),nr) = 0
          end if
          ovar(tt(nr),nr) = ovar(tt(nr),nr) + 1./(ncoarsegrain**2)*sum(var2(ii:ii+ncoarsegrain-1,jj:jj+ncoarsegrain-1,t))
          npts(tt(nr),nr) = npts(tt(nr),nr) + 1
        end if
      end do
    end do
    tt = tt + 1
  end do
  write (*,*) 'Finished Time Loop'

  write (*,*) 'Deallocate var1, var2'
  deallocate(var1)
  deallocate(var2)
  do t=1,tmax
     do nr=1,nrmax
        if (npts(t,nr).gt.0) then
           ovar(t,nr) = ovar(t,nr)/npts(t,nr)
        else
           ovar(t,nr) = -1e9
        end if
     end do
  end do
  write (*,*) '  max ovar = ',maxval(ovar)
  write (*,*) '  min ovar = ',minval(ovar)

!  where (npts > 0)
!    ovar = ovar/npts
!  elsewhere
!    ovar = fillvalue_r
!  end where

  allocate(ncovar%dim(2))
  allocate(ncovar%dimids(2))
  ncovar = ncvar3

  !Write to netcdf file: Cell base
  ncovar%name     = ncvar1%name(3:len_trim(ncvar1%name))//trim(ncvar2%name)
  ncovar%longname = ncvar1%longname(1:len_trim(ncvar1%longname)-6)//trim(ncvar2%longname)
  ncovar%units    = ncvar2%units
!  ncovar%dim      = ncvar3%dim
!  ncovar%dimids   = ncvar3%dimids

  write (*,*) '    ncovar%name = ',trim(ncovar%name)
  write (*,*) '    ncovar%longname = ',trim(ncovar%longname)
  write (*,*) '    ncovar%dim(1) = ',ncovar%dim(1)
  write (*,*) '    ncovar%dim(2) = ',ncovar%dim(2)
  write (*,*) '    ncovar%dimids(1) = ',ncovar%dimids(1)
  write (*,*) '    ncovar%dimids(2) = ',ncovar%dimids(2)
  write (*,*) '    ncvar3%dimids(1) = ',ncvar3%dimids(1)
  write (*,*) '    ncvar3%dimids(2) = ',ncvar3%dimids(2)

  write (*,*) ' Define ncovar in netcdf file'
  call define_ncvar(ifid1, ncovar, nf90_float)

  write (*,*) 'Write to netcdf file: '
  call write_ncvar(ifid1, ncovar, ovar)


end program trackmapper
