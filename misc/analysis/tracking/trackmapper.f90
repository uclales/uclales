
module modnetcdf
  use netcdf
  implicit none
  real(kind=4),    parameter :: fillvalue_r = -huge(1)
  integer, parameter :: fillvalue_i = -huge(1)
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
  end interface write_ncvar

  interface read_ncvar
    module procedure read_ncvar_1D_r
    module procedure read_ncvar_2D_r
    module procedure read_ncvar_3D_r
    module procedure read_ncvar_4D_r
    module procedure read_ncvar_1D_i
    module procedure read_ncvar_2D_i
    module procedure read_ncvar_3D_i
    module procedure read_ncvar_4D_i
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

  subroutine read_ncvar_1D_r(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer         :: fid
    integer, optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=(/dimstart/), count=(/dimrange/)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_1D_r

  subroutine read_ncvar_2D_r(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer         :: fid
    integer, dimension(:), optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=dimstart(1:2), count=dimrange(1:2)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_2D_r

  subroutine read_ncvar_3D_r(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=dimstart(1:3), count=dimrange(1:3)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_3D_r

  subroutine read_ncvar_4D_r(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=dimstart(1:4), count=dimrange(1:4)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_4D_r

  subroutine read_ncvar_1D_i(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    integer, dimension(:) :: array
    integer         :: fid
    integer, optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=(/dimstart/), count=(/dimrange/)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_1D_i

  subroutine read_ncvar_2D_i(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:) :: array
    integer         :: fid
    integer, dimension(:), optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=dimstart(1:2), count=dimrange(1:2)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_2D_i

  subroutine read_ncvar_3D_i(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=dimstart(1:3), count=dimrange(1:3)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_3D_i

  subroutine read_ncvar_4D_i(fid,ncvar,array, dimstart, dimrange)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:,:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart, dimrange

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_get_var(fid,ncvar%id,array, start=dimstart(1:4), count=dimrange(1:4)) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_4D_i

  subroutine write_ncvar_1D_r(fid,ncvar,array, dimstart)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart
    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_put_var(fid,ncvar%id,array, start=dimstart(1:1)) )
    else
      call check ( nf90_put_var(fid,ncvar%id,array) )
    end if
  end subroutine write_ncvar_1D_r

  subroutine write_ncvar_2D_r(fid,ncvar,array, dimstart)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_put_var(fid,ncvar%id,array, start=dimstart(1:2)) )
    else
      call check ( nf90_put_var(fid,ncvar%id,array) )
    end if
  end subroutine write_ncvar_2D_r

  subroutine write_ncvar_3D_r(fid,ncvar,array, dimstart)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_put_var(fid,ncvar%id,array, start=dimstart(1:3)) )
    else
      call check ( nf90_put_var(fid,ncvar%id,array) )
    end if
  end subroutine write_ncvar_3D_r

  subroutine write_ncvar_1D_i(fid,ncvar,array, dimstart)
    type(netcdfvar) :: ncvar
    integer, dimension(:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_put_var(fid,ncvar%id,array, start=dimstart(1:1)) )
    else
      call check ( nf90_put_var(fid,ncvar%id,array) )
    end if
  end subroutine write_ncvar_1D_i

  subroutine write_ncvar_2D_i(fid,ncvar,array, dimstart)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_put_var(fid,ncvar%id,array, start=dimstart(1:2)) )
    else
      call check ( nf90_put_var(fid,ncvar%id,array) )
    end if
  end subroutine write_ncvar_2D_i

  subroutine write_ncvar_3D_i(fid,ncvar,array, dimstart)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:,:) :: array
    integer         :: fid
    integer, dimension(:),optional :: dimstart

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(dimstart)) then
      call check ( nf90_put_var(fid,ncvar%id,array, start=dimstart(1:3)) )
    else
      call check ( nf90_put_var(fid,ncvar%id,array) )
    end if
  end subroutine write_ncvar_3D_i

  subroutine define_ncdim(fid,ncvar)
    type(netcdfvar) :: ncvar
    integer         :: fid, iret

    !..open file and read data
    iret = nf90_redef(fid)
    iret = nf90_def_dim(fid, trim(ncvar%name), ncvar%dim(1), ncvar%dimids(1))
    if (iret == nf90_enameinuse) then
      call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%dimids(1)) )
    else
      call check(iret)
    end if
    call check ( nf90_enddef(fid))
    call define_ncvar(fid, ncvar, nf90_float)

  end subroutine define_ncdim

  subroutine define_ncvar(fid,ncvar, xtype)
    type(netcdfvar) :: ncvar
    integer         :: fid, iret, xtype
    !..open file and read data
    call check ( nf90_redef(fid))
    iret =  nf90_def_var(fid, name = trim(ncvar%name), &
            xtype=xtype, dimids = ncvar%dimids, varid = ncvar%id)
    if (iret == nf90_enameinuse) then
      call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    else
      call check(iret)
      call check ( nf90_put_att(fid, ncvar%id, 'longname', ncvar%longname))
      call check ( nf90_put_att(fid, ncvar%id, 'units', ncvar%units))
      if (xtype == nf90_int) then
        call check ( nf90_put_att(fid, ncvar%id, '_FillValue', fillvalue_i))
      else
        call check ( nf90_put_att(fid, ncvar%id, '_FillValue', fillvalue_r))
      end if
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
  integer :: ifid1, ifid2, i, j, t, nx, ny, nt, nr, tmax, nrmax
  real, allocatable, dimension(:,:,:) :: var1
  integer, allocatable, dimension(:) :: tt

  real, allocatable, dimension(:,:,:) :: var2
  real, allocatable, dimension(:,:) :: ovar, npts
  type(netcdfvar) :: ncvar1, ncvar2, ncvar3, ncovar


  !Read command line options
  call get_command_argument(1,ifile1)
  call get_command_argument(2,ifile2)
  call get_command_argument(3,ncvar1%name)
  call get_command_argument(4,ncvar2%name)
  call get_command_argument(5,ncvar3%name) !Something shapelike identical to the output variable

  !Open files

  call check ( nf90_open(trim(adjustl(ifile1)), NF90_WRITE, ifid1))
  call check ( nf90_open(trim(adjustl(ifile2)), NF90_NOWRITE, ifid2))

  !Find dimensions/sizes
print *, ncvar1%name
  call inquire_ncvar(ifid1, ncvar1)
print *, ncvar3%name
  call inquire_ncvar(ifid1, ncvar3)
print *, ncvar2%name
  call inquire_ncvar(ifid2, ncvar2)
  nx = ncvar1%dim(1)
  ny = ncvar1%dim(2)
  nt = ncvar1%dim(3)

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
  allocate(var2(nx,ny,nt))
  call read_ncvar(ifid1, ncvar1, var1)
  call read_ncvar(ifid2, ncvar2, var2)
!   !Time loop

  do t = 1,nt
    write (*,*) 'Time = ',t
    do i=1,nx
      do j=1,ny
        if (var1(i,j,t) >= 1) then
          nr = var1(i,j,t)
          if (tt(nr)< 0) then ! A new cell
            tt(nr) = 1
            ovar(tt(nr),nr) = 0
          end if
          ovar(tt(nr),nr) = ovar(tt(nr),nr) + var2(i,j,t)
          npts(tt(nr),nr) = npts(tt(nr),nr) + 1
        end if
      end do
    end do
    tt = tt + 1
  end do
  where (npts > 0)
    ovar = ovar/npts
  elsewhere
    ovar = fillvalue_r
  end where
  allocate(ncovar%dim(2))
  allocate(ncovar%dimids(2))
  ncovar = ncvar3

  !Write to netcdf file: Cell base
  ncovar%name     = ncvar1%name(3:len_trim(ncvar1%name))//trim(ncvar2%name)
  ncovar%longname = ncvar1%longname(1:len_trim(ncvar1%longname)-6)//trim(ncvar2%longname)
  ncovar%units    = ncvar2%units
  call define_ncvar(ifid1, ncovar, nf90_float)
  call write_ncvar(ifid1, ncovar, ovar)


end program trackmapper
