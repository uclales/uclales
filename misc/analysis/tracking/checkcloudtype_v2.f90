
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
    module procedure read_ncvar_1D_r
    module procedure read_ncvar_2D_r
    module procedure read_ncvar_3D_r
    module procedure read_ncvar_1D_i
    module procedure read_ncvar_2D_i
    module procedure read_ncvar_3D_i
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

  subroutine read_ncvar_1D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_get_var(fid,ncvar%id,array) )

  end subroutine read_ncvar_1D_r

  subroutine read_ncvar_2D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_get_var(fid,ncvar%id,array) )

  end subroutine read_ncvar_2D_r

  subroutine read_ncvar_3D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_get_var(fid,ncvar%id,array) )

  end subroutine read_ncvar_3D_r

  subroutine read_ncvar_1D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_get_var(fid,ncvar%id,array) )

  end subroutine read_ncvar_1D_i

  subroutine read_ncvar_2D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_get_var(fid,ncvar%id,array) )

  end subroutine read_ncvar_2D_i

  subroutine read_ncvar_3D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:,:) :: array
    integer, dimension(:,:,:),allocatable :: arr
    integer         :: fid, i, nstep, nsize
    integer, parameter :: maxsize = 1e8
    !..open file and read data
    call inquire_ncvar(fid, ncvar)
    nstep = 100
!     do while(ncvar%dim(1) * ncvar%dim(2) * nstep > maxsize)
!       nstep  = nstep /2
!     end do

    allocate(arr(ncvar%dim(1),ncvar%dim(2),nstep))

    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    do i = 0, ceiling(1.*ncvar%dim(3)/nstep)-1
      write (*,*) 'i =',i
      nsize = min(nstep,ncvar%dim(3)-i*nstep)
      call check ( nf90_get_var(fid,ncvar%id,arr,start=(/1,1,1+i*nstep/), &
                   count=(/ncvar%dim(1),ncvar%dim(2),nsize/) ))
      array(:,:,i*nstep+1:i*nstep+nsize) = arr(:,:,1:nsize)
    end do

  end subroutine read_ncvar_3D_i

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

program checkcloudtype
  use modnetcdf
  implicit none
  !Declare variables
  integer :: i,j,t,n, nt,nx,ny, finput,  nsib, nsplit
  character(80) :: filename
  integer, allocatable, dimension(:,:,:) :: nrcloud, cloudtype
  integer, allocatable, dimension(:)    :: itype, sib, split
  type(netcdfvar) :: sibnc, splitnc, nrcloudnc, cloudtypenc

  !Open Netcdf
  filename = 'rico.out.xy.track.nc'
    !Read input
  write(*,*) 'Reading dimensions'
  call check ( nf90_open (trim(filename), NF90_WRITE, finput) )
  tnc%name = 'time'
  call inquire_ncvar(finput, tnc)
  nt = tnc%dim(1)

  xnc%name = 'xt'
  call inquire_ncvar(finput, xnc)
  nx = xnc%dim(1)

  ync%name = 'yt'
  call inquire_ncvar(finput, ync)
  ny = ync%dim(1)

  !Load nrsib, nrsplit
  splitnc%name = 'smcloudtype'
  call inquire_ncvar(finput, splitnc)
  nsplit = splitnc%dim(1)
  allocate(itype(nsplit))
  call read_ncvar(finput, splitnc, itype)
!
!   sibnc%name = 'smcloudsystemnr'
!   call inquire_ncvar(finput, sibnc)
!   nsib = sibnc%dim(1)
!   allocate(sib(nsib))
!   call read_ncvar(finput, sibnc, sib)
!   allocate(itype(nsib))
  !Load nrcloud
  nrcloudnc%name = 'nrcloud'
  call inquire_ncvar(finput, nrcloudnc)
  allocate(nrcloud(nx,ny,nt))
  call read_ncvar(finput, nrcloudnc, nrcloud)

  allocate(cloudtype(nx,ny,nt))
  cloudtype = nrcloud
!
!   do n = 1, nsplit
!     if (split(n)==0 .and. sib(n)==0) then
!       itype(n) = 1 !Passive
!     elseif (split(n)==1 .and. sib(n)==0) then
!       itype(n) = 2 !Single pulse
!     elseif (split(n)==0 .and. sib(n)>0) then
!       itype(n) = 3 !Outflow
!     elseif (split(n)==1 .and. sib(n)>0) then
!       itype(n) = 4 !Multi-pulse
!     end if
!   end do
  do t = 1,nt
    print *, 't=',t
    do j = 1,ny
      do i =1,nx
	if (cloudtype(i,j,t)>0) then
	  cloudtype(i,j,t) = itype(cloudtype(i,j,t))
	end if
      end do
    end do
  end do
  cloudtypenc = nrcloudnc
  cloudtypenc%name = 'cloudtype'
  cloudtypenc%longname = 'Cloud type. 1 for passive, 2 for single pulse, 3 for outflow, 4 for active'

  call define_ncvar(finput, cloudtypenc, nf90_int)
  call write_ncvar(finput, cloudtypenc, cloudtype)
  call check ( nf90_close(finput))

end program checkcloudtype