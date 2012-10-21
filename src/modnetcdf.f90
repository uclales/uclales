!  This file is part of MicroHH.
!
! MicroHH is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! MicroHH is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2010-2011 Chiel van Heerwaarden and Thijs Heus
!> Background routines to read and write NetCDF output
!! All calls to the netcdf library should be directed through here.
!! The module opens (with open_nc), closes (with close_nc),
!! reads (with readvar_nc) and writes (with writevar_nc) anything between 0D (e.g. timeseries)
!! and 4D (e.g. z,x,y,t) fields to file.
!! Inspired on the UCLA-LES routine by Bjorn Stevens.
!! \todo Parallel NETCDF
!! \todo Documentation
module modnetcdf
  use defs, only : single, short, double, int32
  implicit none

  real(kind=single), parameter    :: fillvalue_single = -32678. !< Fill value
  real(kind=double), parameter    :: fillvalue_double = -32678. !< Fill value
  integer(kind=short), parameter  :: fillvalue_short  = -32678  !< Fill value
  integer, parameter :: icomp_float = 0 !< Write back floats (default)
  integer, parameter :: icomp_dbl   = 1 !< Write back doubles instead of floats
  integer, parameter :: icomp_int   = 2 !< Write back integers instead of floats
  logical            :: lsync = .false. !< Synchronization switch
  integer            :: deflate_level = 0 !< Compression level

!> Interface to write into a netcdf file
  interface writevar_nc
    module procedure writevar0D_nc
    module procedure writevar1D_nc
    module procedure writevar2D_nc
    module procedure writevar3D_nc
  end interface writevar_nc

!> Interface to read a netcdf field
  interface readvar_nc
    module procedure readvar0D_nc
    module procedure readvar1D_nc
    module procedure readvar2D_nc
    module procedure readvar3D_nc
    module procedure readvar4D_nc
  end interface readvar_nc

!> Interface to put attributes in a netcdf file
  interface putatt_nc
    module procedure putatt_single_nc
    module procedure putatt_double_nc
    module procedure putatt_int_nc
    module procedure putatt_short_nc
    module procedure putatt_str_nc
  end interface putatt_nc

!> Interface to put attributes in a netcdf file
  interface getatt_nc
    module procedure getatt_single_nc
    module procedure getatt_double_nc
    module procedure getatt_int_nc
    module procedure getatt_short_nc
    module procedure getatt_str_nc
    module procedure getatt_str_fname_nc
  end interface getatt_nc

contains
  
!> Subroutine Open_NC: Opens a NetCDF File for writing
!! If the file already exists, the record number is being set to the current
!! time in the simulation
  subroutine open_nc(fname, ncid,nrec, rtimee, ldelete)
    use grid,   only : tname
    use netcdf, only : nf90_create, nf90_open, nf90_inquire, nf90_inquire_dimension, &
                       nf90_inq_varid, nf90_get_var, nf90_put_var, nf90_sync, nf90_share, &
                       nf90_write, nf90_noerr, nf90_inquire_variable, nf90_short
    integer, intent (out)          :: ncid   !< NetCDF file number
    integer, intent (out),optional :: nrec   !< The record number that corresponds to the current simulation time
    real, intent(in), optional     :: rtimee !< Simulation time
    logical, intent(in), optional  :: ldelete !< Whether to delete an old file when present
    character (len=*), intent (in) :: fname  !< File name
    integer :: ndims, nvars, n, nn, xtype, dimsize(7)
    integer, allocatable, dimension(:) :: start, dimids

    character (len=12) :: date, time
    integer            :: iret,  ncall, RecordDimID,timeid, totdimsize
    real, allocatable  :: xtimes(:)
    logical            :: exists, ldef

    if (present(rtimee)) then ! Check whether file exists
      inquire(file=trim(fname),exist=exists)
      if (exists .and. present(ldelete)) then
        if (ldelete) then
          exists = .false.
          open(1, file=trim(fname), status='old')
          close(1,status='delete')
        end if
      end if
    else
      exists = .false.
    end if
    ncall = 0
    if (.not.exists) then ! Write the header of the file, and ensure we're in data mode
      call date_and_time(date,time)
      iret = nf90_create(fname,NF90_SHARE,ncid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncid, 'title', trim(fname))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncid,'history','Created on '//date(1:8)//' at '//trim(time))
      ldef = ensuredata_nc(ncid)
    else ! If the file exists, only the record number needs to be set
      iret  = nf90_open (trim(fname), NF90_WRITE, ncid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_inquire_dimension(ncid, RecordDimID, len=nrec)
      if (iret==0) then
        if (nrec > 0) then
          iret  = nf90_inq_varid(ncid,tname,timeID)
          if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
          allocate (xtimes(0:nrec))
          iret  = nf90_get_var(ncid, timeId, xtimes(0:nrec-1))
          if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
          do while(ncall < nrec .and. &
                   xtimes(ncall) /= fillvalue_double .and. &
                   xtimes(ncall) <= rtimee - spacing(1.)) ! Step through the time dimension; stop when one is bigger
              ncall=ncall+1
          end do
          ldef = ensuredata_nc(ncid)
          xtimes = fillvalue_double
          iret = nf90_inquire(ncid,nvariables = nvars)
          do n = 1, nvars
            iret = nf90_inquire_variable(ncid, n, xtype = xtype, ndims = ndims)
            allocate(dimids(ndims))
            iret = nf90_inquire_variable(ncid, n, dimids = dimids)
            if (any(dimids == recorddimid)) then
              allocate(start(ndims))
              dimsize = 1
              start = 1
              do nn = 1, ndims
                if (dimids(nn) == recorddimid) then
                  start(nn) = ncall + 1
                  dimsize(nn) = nrec - ncall
                else
                  iret = nf90_inquire_dimension(ncid, dimids(nn), len = dimsize(nn))
                end if
              end do
              totdimsize = product(dimsize)
              select case (xtype)
              case(NF90_SHORT)
                iret = nf90_put_var(ncid, n, reshape((/(fillvalue_short, n = 1, totdimsize)/),dimsize), start = start)
              case default
                iret = nf90_put_var(ncid, n, reshape((/(fillvalue_double, n = 1, totdimsize)/),dimsize), start = start)
              end select
              deallocate(start)
            end if
            deallocate(dimids)
          end do
          if (ldef) ldef = ensuredefine_nc(ncid)
           
          deallocate(xtimes)
        end if
      end if
    end if
    if (present(nrec)) nrec = ncall
    iret = nf90_sync(ncid)
  end subroutine open_nc

!> Switch NetCDF file to define mode. Returns true if the dataset already was in define mode, and false if it was in data mode
  logical function ensuredefine_nc(ncid)
    use netcdf, only : nf90_redef, nf90_eindefine, nf90_noerr
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    ensuredefine_nc = .false.
    iret = nf90_redef(ncid)
    select case(iret)
    case(nf90_noerr)
      ensuredefine_nc = .false.
    case(nf90_eindefine)
      ensuredefine_nc = .true.
    case default
      call nchandle_error(ncid, iret)
    end select    
 end function ensuredefine_nc

!> Switch NetCDF file to data mode. Returns true if the dataset was in define mode, and false if it was already in data mode
  logical function ensuredata_nc(ncid)
    use netcdf, only : nf90_enddef, nf90_enotindefine, nf90_noerr
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    ensuredata_nc = .false.
    iret = nf90_enddef(ncid)
    select case(iret)
    case(nf90_noerr)
      ensuredata_nc = .true.
    case(nf90_enotindefine)
      ensuredata_nc = .false.
    case default
      call nchandle_error(ncid, iret)
    end select    
 end function ensuredata_nc

!> Switch NetCDF file to data mode
  subroutine enddefine_nc(ncid)
    use netcdf, only : nf90_enddef, nf90_noerr
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    iret = nf90_enddef(ncid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
  end subroutine enddefine_nc

!> Close a NetCDF file
  subroutine close_nc(ncid)
    use netcdf, only : nf90_close, nf90_inquire, nf90_noerr
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    iret = nf90_inquire(ncid)
    if (iret /= nf90_noerr) return
    iret = nf90_close(ncid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
  end subroutine close_nc

!> Add a variable to a NetCDF file. The procedure checks whether the variable(name) is already present,
!! and if not, adds it to the file. The same holds for the dimension(names) of the variable.
  subroutine addvar_nc(ncID, name, lname, unit, dimname, dimlongname, dimunit, dimsize, dimvalues, icompress, datamin, datamax)
    use netcdf, only : nf90_def_var, nf90_inq_varid, nf90_float, nf90_short, nf90_double, nf90_noerr, nf90_def_var_deflate
    integer, intent (in)                              :: ncID        !< NetCDF file number
    character (*), intent (in)                        :: name        !< Netcdf name of the variable
    character (*), intent (in)                        :: lname       !< Longname of the variable
    character (*), intent (in)                        :: unit        !< Unit of the variable
    character (*), dimension(:), intent(in), optional :: dimname     !< NetCDF names of the dimensions of the variable (in array form)
    character (*), dimension(:), intent(in), optional :: dimlongname !< Longnames of the dimensions of the variable (in array form)
    character (*), dimension(:), intent(in), optional :: dimunit     !< Units of the dimensions of the variable (in array form)
    integer, dimension(:), intent(in), optional       :: dimsize     !< List of dimension sizes; 0 for unlimited
    real, dimension(:,:), intent(in), optional        :: dimvalues   !< List of values of the dimension
    integer, optional                                 :: icompress   !< Choice of compression: 0 for 4 bytes float (default), 1 for 8 bytes double, 2 for 2 bytes integer
    real, optional                                    :: datamin     !< Minimum value in case of compression
    real, optional                                    :: datamax     !< Maximum value in case of compression
  
    real                                        :: scale_factor, add_offset
    integer                                     :: iret, n, nrdim,VarID, icomp, datatype
    integer, allocatable, dimension(:)          :: ncdim
    logical                                     :: ldef

    ldef = ensuredefine_nc(ncid) ! Need to be in define mode
    if (present(dimname)) then
      nrdim = size(dimname)
    else
      nrdim = 0
    end if
    if (present(icompress)) then
      icomp= icompress
    else
      icomp = icomp_float
    end if
    allocate (ncdim(nrdim))
    do n=1,nrdim ! For every dimension, check whether it exists in the file. If not, inqdim_nc will write it.
      if (present(dimvalues)) then
        ncdim(n) = inqdimid_nc(ncid ,dimname(n) ,dimlongname(n) ,dimunit(n) ,dimvalues(1:dimsize(n),n))
      else
        ncdim(n) = inqdimid_nc(ncid ,dimname(n) ,dimlongname(n) ,dimunit(n), (/0./))
      end if
    end do
    iret=nf90_inq_varid(ncid,trim(validate(name)),varid)
    if (iret /= nf90_noerr) then ! If the inquiry fails, the variable needs to be created
      select case (icomp)
      case (icomp_int)
        datatype = NF90_SHORT
      case (icomp_dbl)
        datatype = NF90_DOUBLE
      case default
        datatype = NF90_FLOAT
      end select
      iret=nf90_def_var(ncID, trim(validate(name)), datatype, ncdim,VarID, deflate_level = deflate_level)
      if (iret/=0) then
        write (*,*) 'Variable ',trim(name), ncdim
        call nchandle_error(ncid, iret)
      end if
      if (nrdim > 2 .and. deflate_level > 0) then
        iret=nf90_def_var_deflate(ncid, varid, 0, 1, deflate_level = deflate_level)
      end if
      iret = putatt_nc(ncID, 'longname', lname, trim(name))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, 'units', unit, trim(name))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      select case (icomp)
      case (icomp_int)
        scale_factor = (datamax - datamin) / (2**16 - 1)
        add_offset   = datamin + 2**15 * scale_factor
        iret = putatt_nc(ncID, 'scale_factor', scale_factor, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
        iret = putatt_nc(ncID, 'add_offset', add_offset, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
        iret = putatt_nc(ncID, '_FillValue', fillvalue_short, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      case (icomp_dbl)
        iret = putatt_nc(ncID, '_FillValue', fillvalue_double, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      case default
        iret = putatt_nc(ncID, '_FillValue', fillvalue_single, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      end select
    end if
    if (ldef .eqv. .false.) then
      call enddefine_nc(ncid) ! Leave define mode again
    end if
    deallocate (ncdim)
  end subroutine addvar_nc

!> Check whether in the give file a dimension with the given dimension name already exists.
!! If so, check whether the extent and values are correct.
!! If the dimension does not exist, create it and fill the related variable.
!! Inqdimid_nc sets the netcdf number of the resulting dimension as a return value.
  function inqdimid_nc(ncid,dimname,dimlongname,dimunit,dimvalues)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_put_var, &
                       nf90_def_dim, nf90_def_var, nf90_unlimited, nf90_double, nf90_noerr
    integer                         :: inqdimid_nc !< The netcdf dimension number (return value)
    integer, intent(in)             :: ncid        !< Netcdf file number
    character(*), intent(in)        :: dimname     !< Name of the dimension
    character(*), intent(in)        :: dimlongname !< Netcdf long name of the dimension
    character(*), intent(in)        :: dimunit     !< Unit of the dimension
    real, dimension(:), intent(in)  :: dimvalues   !< Values of the dimension. 0 for the unlimited dimension
    logical :: ltime, ldef
    integer :: iret,varid, dimsize
    real, dimension(:), allocatable :: dimvar
    character(80) :: fname

    if (all(dimvalues == 0)) then ! If all values of this dimension are zero, it is the unlimited (time) dimension
      ltime = .true.
    else
      ltime = .false.
    end if
    iret = 0
    iret = nf90_inq_dimid(ncid, trim(validate(dimname)), inqdimid_nc)
    if (iret == 0) then !If the dimension already exists....
      if (.not. ltime) then
        iret = nf90_inquire_dimension(ncID, inqdimid_nc, len = dimsize)
        allocate(dimvar(dimsize))
        if (dimsize == size(dimvalues)) then ! Check whether the number of points along this dimension is correct
          ldef = ensuredata_nc(ncid)
          iret =  nf90_inq_varid(ncid, trim(validate(dimname)), varid)
          iret  = nf90_get_var(ncid,varid,dimvar)
          ldef = ensuredefine_nc(ncid)
          if (any((abs((dimvar - dimvalues)/(dimvar+epsilon(1.)))) > 1e-4)) then ! Check whether dimension in file matches with the desired values
            iret = getatt_nc(ncid,'title', fname)
            if (iret /= nf90_noerr) fname = ''
            print *, 'NetCDF error in file ' // trim(fname)
            print *, 'For dimension ', trim(dimname)
            print *, "Inqdimid_nc: Dimensions don't match with what's already in the file"
            call appl_finalize (-1)
          end if
        else
          iret = getatt_nc(ncid, 'title', fname)
          if (iret /= nf90_noerr)  fname = ''
          print *, 'NetCDF error in file ' // trim(fname)
          print *, 'For dimension ', trim(dimname)
          print *, "Number of points doesn't fit"
          call appl_finalize(-1)
        end if
        deallocate(dimvar)
      end if
    else !If the dimension does not exist yet, we need to create it.
      if (ltime) then
        iret = nf90_def_dim(ncID, trim(validate(dimname)), NF90_UNLIMITED, inqdimid_nc)
      else
        iret = nf90_def_dim(ncID, trim(validate(dimname)), size(dimvalues), inqdimid_nc)
      endif
      if (iret/=0) then
        print *,'For dimension ',trim(dimname)
        call nchandle_error(ncid, iret)
      end if
      iret = nf90_def_var(ncID,trim(validate(dimname)), nf90_double, inqdimid_nc,VarID)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, 'longname', dimlongname, dimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, 'units', dimunit, dimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, '_FillValue', fillvalue_double, dimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
     if (.not. ltime) then !Fill the dimension-values for any dimension but time
        ldef = ensuredata_nc(ncid)
        iret = nf90_put_var(ncid, varID, dimvalues,(/1/))
        if (iret/=0) then
          print *,'For dimension ',trim(dimname)
          call nchandle_error(ncid, iret)
        end if
        ldef = ensuredefine_nc(ncid)
      end if
    end if
  end function inqdimid_nc

!> Write down a number of variables that depend on (possibly) time and 1 other dimension.
  subroutine writevar0D_nc(ncid,ncname,var,nrec)
    use netcdf, only : nf90_inquire, nf90_inquire_dimension, nf90_inq_varid, &
                       nf90_put_var, nf90_noerr, nf90_erange
    integer, intent(in)              :: ncid   !< Netcdf file number
    character(len=*),intent(in)      :: ncname !< The variable names
    real,intent(in)                  :: var    !< The variables to be written to file
    integer, intent(inout), optional :: nrec   !< Netcdf Record number

    integer :: iret,varid, udimid, loc(1)
    real    :: scale_factor, add_offset
    integer(kind=single)  :: var_i
    character(len=20) :: udimname
    iret = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (present(nrec)) then
      iret = nf90_inquire(ncid,unlimitedDimId=udimid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = nf90_inquire_dimension(ncid,udimid,name=udimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      if(trim(ncname)==trim(udimname)) nrec = nrec+1
      loc = nrec
    else
      loc = 1
    end if
    iret = getatt_nc(ncid, 'scale_factor', scale_factor, ncname)
    if (iret == 0) then !We're doing compression
      iret = getatt_nc(ncid, 'add_offset', add_offset, ncname)
      if (iret /=0) add_offset = 0.
      var_i = nint((var - add_offset) / scale_factor)
      iret = nf90_put_var(ncid, VarID, var_i,loc)
    else
      iret = nf90_put_var(ncid, VarID, var,loc)
    end if      
    if (iret /= nf90_noerr .and. iret /= nf90_erange) call nchandle_error(ncid, iret)
    iret = sync_nc(ncid)
  end subroutine writevar0D_nc

!> Write down a number of variables that depend on (possibly) time and 1 other dimension.
  subroutine writevar1D_nc(ncid,ncname,var,nrec)
    use netcdf, only : nf90_inquire_dimension, nf90_inquire_variable, nf90_inq_varid, nf90_put_var, nf90_noerr, nf90_erange
    integer, intent(in)           :: ncid   !< Netcdf file number
    character(len=*),intent(in)   :: ncname !< The variable names
    real,dimension(:),intent(in)  :: var    !< The variables to be written to file
    integer, intent(in), optional :: nrec   !< Netcdf Record number

    integer :: iret,varid
    integer :: dimids(2), dimsize(2)
    integer, dimension(:), allocatable :: loc
    real    :: scale_factor, add_offset
    integer(kind=single), allocatable, dimension(:) :: var_i
    dimsize(2) = 1 !The time dimension has size 1
    iret = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_variable(ncid,varid,dimids=dimids)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_dimension(ncid,dimids(1),len=dimsize(1))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (present(nrec)) then
      allocate(loc(2))
      loc = (/1, nrec/)
    else
      allocate(loc(1))
      loc = (/1/)
    end if
    iret = getatt_nc(ncid, 'scale_factor', scale_factor, ncname)
    if (iret == 0) then !We're doing compression
      iret = getatt_nc(ncid, 'add_offset', add_offset, ncname)
      if (iret /=0) add_offset = 0.
      if (.not. allocated(var_i)) allocate(var_i(1:dimsize(1)))
      var_i = nint((var(1:dimsize(1)) - add_offset) / scale_factor)
      iret = nf90_put_var(ncid, VarID, var_i(1:dimsize(1)),loc,dimsize)
    else
      iret = nf90_put_var(ncid, VarID, var(1:dimsize(1)),loc,dimsize)
    end if      
    if (iret /= nf90_noerr .and. iret /= nf90_erange) call nchandle_error(ncid, iret)
    iret = sync_nc(ncid)
    if (allocated(var_i)) deallocate(var_i)
  end subroutine writevar1D_nc

!> Write down a number of variables that depend on (possibly) time and 2 other dimension.
  subroutine writevar2D_nc(ncid,ncname,var,nrec)
    use netcdf, only : nf90_inquire_dimension, nf90_inquire_variable, nf90_inq_varid, nf90_put_var, nf90_noerr, nf90_erange
    integer, intent(in)             :: ncid   !< Netcdf file number
    character(len=*),intent(in)     :: ncname !< The variable names
    real,dimension(:,:),intent(in)  :: var    !< The variables to be written to file
    integer, intent(in), optional   :: nrec   !< Netcdf Record number

    integer :: iret,varid
    integer :: dimids(3), dimsize(3)
    integer, dimension(:), allocatable :: loc
    real    :: scale_factor, add_offset
    integer(kind=single), allocatable, dimension(:,:) :: var_i
    dimsize(3) = 1 !The time dimension has size 1
    iret = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_variable(ncid,varid,dimids=dimids)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_dimension(ncid,dimids(1),len=dimsize(1))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_dimension(ncid,dimids(2),len=dimsize(2))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (present(nrec)) then
      allocate(loc(3))
      loc = (/1, 1, nrec/)
    else
      allocate(loc(2))
      loc = (/1, 1/)
    end if
    iret = getatt_nc(ncid, 'scale_factor', scale_factor, ncname)
    if (iret == 0) then !We're doing compression
      iret = getatt_nc(ncid, 'add_offset', add_offset, ncname)
      if (iret /=0) add_offset = 0.
      allocate(var_i(1:dimsize(1), 1:dimsize(2)))
      var_i = nint((var(1:dimsize(1), 1:dimsize(2)) - add_offset) / scale_factor)
      iret = nf90_put_var(ncid, VarID, var_i(1:dimsize(1), 1:dimsize(2)),loc,dimsize)
      deallocate(var_i)
    else
      iret = nf90_put_var(ncid, VarID, var(1:dimsize(1),1:dimsize(2)),loc,dimsize)
    end if      
    if (iret /= nf90_noerr .and. iret /= nf90_erange) call nchandle_error(ncid, iret)
    iret = sync_nc(ncid)
    
  end subroutine writevar2D_nc

!> Write down a number of variables that depend on (possibly) time and 3 other dimension.
  subroutine writevar3D_nc(ncid,ncname,var,nrec)
    use netcdf, only : nf90_inquire_dimension, nf90_inquire_variable, nf90_inq_varid, nf90_put_var, nf90_noerr, nf90_erange
    integer, intent(in)              :: ncid   !< Netcdf file number
    character(len=*),intent(in)      :: ncname !< The variable names
    real,dimension(:,:,:),intent(in) :: var    !< The variables to be written to file
    integer, intent(in), optional    :: nrec   !< Netcdf Record number

    integer :: iret,varid
    integer :: dimids(4), dimsize(4)
    integer, dimension(:), allocatable :: loc
    real    :: scale_factor, add_offset
    integer(kind=single), allocatable, dimension(:,:,:) :: var_i
    dimsize(4) = 1 !The time dimension has size 1
    iret = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_variable(ncid,varid,dimids=dimids)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_dimension(ncid,dimids(1),len=dimsize(1))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_dimension(ncid,dimids(2),len=dimsize(2))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret = nf90_inquire_dimension(ncid,dimids(3),len=dimsize(3))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (present(nrec)) then
      allocate(loc(4))
      loc = (/1, 1, 1, nrec/)
    else
      allocate(loc(3))
      loc = (/1, 1, 1/)
    end if
    iret = getatt_nc(ncid, 'scale_factor', scale_factor, ncname)
    if (iret == 0) then !We're doing compression
      iret = getatt_nc(ncid, 'add_offset', add_offset, ncname)
      if (iret /=0) add_offset = 0.
      allocate(var_i(1:dimsize(1),1:dimsize(2),1:dimsize(3)))
      var_i = nint((var(1:dimsize(1),1:dimsize(2),1:dimsize(3)) - add_offset) / scale_factor)
      iret = nf90_put_var(ncid, VarID, var_i(1:dimsize(1),1:dimsize(2),1:dimsize(3)),loc,dimsize)
      deallocate(var_i)
    else
      iret = nf90_put_var(ncid, VarID, var(1:dimsize(1),1:dimsize(2),1:dimsize(3)),loc,dimsize)
    end if      
    if (iret /= nf90_noerr .and. iret /= nf90_erange) call nchandle_error(ncid, iret)

    iret = sync_nc(ncid)
    
  end subroutine writevar3D_nc

!> Reads out a variable that has 0 dimensions.
!! If the variable cannot be read, the netcdf errorcode will be passed on to the calling routine.
  function readvar0D_nc(fname,ncname,var)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_open, nf90_close, nf90_get_var, nf90_inquire_variable, nf90_inq_varid, nf90_nowrite, nf90_noerr
    integer                 :: readvar0D_nc !< Return value. Negative if variable cannot be read.
    character(*),intent(in) :: fname        !< File name
    character(*),intent(in) :: ncname       !< Variable name
    real,intent(out)        :: var          !< Variable value
    integer, parameter      :: ndim_exp = 0
    integer                 :: ncid, varid,  ndims, iret
    logical :: exists
    
    readvar0D_nc = 0
!Open file
    inquire(file=trim(fname),exist=exists)
    if (.not. exists) then
      readvar0D_nc = -1
      return
    end if
    iret  = nf90_open (trim(fname), NF90_NOWRITE, ncid)
    if (iret /= nf90_noerr) then
      readvar0D_nc = -1
      return
    end if
!Check whether variable is available
    iret  = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) then
      readvar0D_nc = iret
      iret         = nf90_close(ncid)
      return
    end if
!Check the number of dimensions in the file
    iret = nf90_inquire_variable (ncid,varid,ndims=ndims)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (ndims > ndim_exp) then
      print *, 'NetCDF Read error in file ' //trim(fname)
      print *, 'Too many dimensions (',ndims ,' instead of ',ndim_exp,')'
      print *, 'For variable ' // trim(ncname)
      call appl_finalize(-1)
    end if

!Read the variable
    iret = nf90_get_var(ncid,varid,var)
    if (iret /= nf90_noerr) then
      print *, 'Reading variable ', trim(ncname)
      call nchandle_error(ncid, iret)
    end if
    iret = nf90_close(ncid)
  end function readvar0D_nc

!> Reads out a variable that has up to 1 dimension.
!! If the variable cannot be read, the netcdf errorcode will be passed on to the calling routine.
!! If the dimensionality is lower than 1, the return value will be equal to the number of missing dimensions.
  function readvar1D_nc(fname, ncname, dimname, dimvar1, var, lslice)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_open, nf90_close, nf90_get_var, nf90_inquire_variable, nf90_inquire_dimension, &
                       nf90_inq_varid, nf90_nowrite, nf90_noerr
    integer                                          :: readvar1D_nc !< Return value. Negative if variable cannot be read, positive if less than 1 dimensions available.
    character(*),intent(in)                          :: fname        !< File name
    character(*),intent(in)                          :: ncname       !< Variable name
    character(*),dimension(1),intent(inout)          :: dimname      !< Name of the dimension (in array form)
    real,allocatable,dimension(:),intent(inout)      :: dimvar1      !< Values of the 1st dimension
    real,allocatable,dimension(:),intent(out)        :: var          !< Variable value
    logical, intent(in), optional                    :: lslice
    integer, parameter                     :: ndim_exp = 1
    integer                                :: ncid, varid, ndims, dimid, i, n, nn, iret, locmin, locmax
    integer, dimension(ndim_exp)           :: dimsize, dimids, start, nslice
    character(len=80), dimension(ndim_exp) :: dimname_in
    real, dimension(ndim_exp)              :: dimmax, dimmin
    real, allocatable, dimension(:,:)      :: dimvars
    logical :: exists
    
    readvar1D_nc = 0
    
    start         = 1
    dimsize       = 1
    dimmin        = -1
    dimmax        = dimmin - 1
    nslice        = 0
    dimname_in(:) = '########'
    if (present(lslice)) then
      if (lslice) then
        if (allocated(dimvar1)) then
          dimmax(1) = maxval(dimvar1)
          dimmin(1) = minval(dimvar1)
          dimname_in(1) = dimname(1)
        end if
      end if
    end if
    dimname = dimname_in
!Open file
    inquire(file=trim(validate(fname)),exist=exists)
    if (.not. exists) then
      readvar1D_nc = -1
      return
    end if
    iret  = nf90_open (trim(validate(fname)), NF90_NOWRITE, ncid)
    if (iret /= nf90_noerr) then
      readvar1D_nc = -1
      return
    end if
!Check whether variable is available
    iret  = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) then
      readvar1D_nc = iret
      iret         = nf90_close(ncid)
      return
    end if

!Check the dimensions of the variable
    iret = nf90_inquire_variable (ncid,varid,ndims=ndims)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (ndims < ndim_exp) then
      readvar1D_nc = ndim_exp - ndims
    else if (ndims > ndim_exp) then
      print *, 'NetCDF Read error in file ' //trim(fname)
      print *, 'Too many dimensions (',ndims ,' instead of ',ndim_exp,')'
      print *, 'For variable ' // trim(ncname)
      call appl_finalize(-1)
    end if
    if (ndims > 0 ) then
      iret = nf90_inquire_variable (ncid,varid,dimids=dimids)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      do n = 1, ndims
        iret = nf90_inquire_dimension(ncid,dimids(n),name=dimname(n),len=dimsize(n))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
        do nn = 1, ndim_exp 
          if (dimname(n)(1:len_trim(dimname_in(nn))) == trim(dimname_in(nn))) then
            nslice(n) = nn
            exit 
          end if
        end do
      end do
    end if
    allocate(dimvars(maxval(dimsize),ndim_exp))
    dimvars = 0.
    do n=1,ndims
      iret  = nf90_inq_varid(ncid, validate(dimname(n)), dimid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_get_var(ncid,dimid, dimvars(1:dimsize(n), n))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      if (nslice(n) > 0) then
        locmax = dimsize(n)
        locmin = 1
        do i = 1, dimsize(n)
          if (dimvars(i, n) >= dimmax(nslice(n))) locmax = max(i, locmax)
          if (dimvars(i, n) <= dimmin(nslice(n))) locmin = min(i, locmin)              
        end do
        dimsize(n) = locmax - locmin + 1
        start(n)   = locmin
      end if
    end do
    
    if (allocated(dimvar1)) deallocate(dimvar1)
    allocate(dimvar1(dimsize(1)))
    dimvar1 = dimvars(start(1):dimsize(1),1)
    
!Read the variable
    if (allocated(var    )) deallocate(var   )
    allocate(var(dimsize(1)))
    iret  = nf90_get_var(ncid,varid,var, start = start, count = dimsize)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret  = nf90_close(ncid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
  end function readvar1D_nc

!> Reads out a variable that has up to 2 dimensions.
!! If the variable cannot be read, the netcdf errorcode will be passed on to the calling routine.
!! If the dimensionality is lower than 2, the return value will be equal to the number of missing dimensions.
  function readvar2D_nc(fname, ncname, dimname, dimvar1, dimvar2, var, lslice)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_open, nf90_close, nf90_get_var, nf90_inquire_variable, nf90_inquire_dimension, &
                       nf90_inq_varid, nf90_nowrite, nf90_noerr
    integer                                          :: readvar2D_nc !< Return value. Negative if variable cannot be read, positive if less than 2 dimensions available.
    character(*),intent(in)                          :: fname        !< File name
    character(*),intent(in)                          :: ncname       !< Variable name
    character(*),dimension(2),intent(inout)          :: dimname      !< Name of the dimension (in array form)
    real,allocatable,dimension(:),intent(inout)      :: dimvar1      !< Values of the 1st dimension
    real,allocatable,dimension(:),intent(inout)      :: dimvar2      !< Values of the 2nd dimension
    real,allocatable,dimension(:,:),intent(out)      :: var          !< Variable value
    logical, intent(in), optional                    :: lslice
    integer, parameter                     :: ndim_exp = 2
    integer                                :: ncid, varid, ndims, dimid, i, n, nn, iret, locmin, locmax
    integer, dimension(ndim_exp)           :: dimsize, dimids, start, nslice
    character(len=80), dimension(ndim_exp) :: dimname_in
    real, dimension(ndim_exp)              :: dimmax, dimmin
    real, allocatable, dimension(:,:)      :: dimvars
    logical :: exists
    
    readvar2D_nc = 0
    
    start         = 1
    dimsize       = 1
    dimmin        = -1
    dimmax        = dimmin - 1
    nslice        = 0
    dimname_in(:) = '########'
    if (present(lslice)) then
      if (lslice) then
        if (allocated(dimvar1)) then
          dimmax(1) = maxval(dimvar1)
          dimmin(1) = minval(dimvar1)
          dimname_in(1) = dimname(1)
        end if
        if (allocated(dimvar2)) then
          dimmax(2) = maxval(dimvar2)
          dimmin(2) = minval(dimvar2)
          dimname_in(2) = dimname(2)
        end if
      end if
    end if
    dimname = dimname_in
!Open file
    inquire(file=trim(fname),exist=exists)
    if (.not. exists) then
      readvar2D_nc = -1
      return
    end if
    iret  = nf90_open (trim(fname), NF90_NOWRITE, ncid)
    if (iret /= nf90_noerr) then
      readvar2D_nc = -1
      return
    end if
!Check whether variable is available
    iret  = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) then
      readvar2D_nc = iret
      iret         = nf90_close(ncid)
      return
    end if

!Check the dimensions of the variable
    iret = nf90_inquire_variable (ncid,varid,ndims=ndims)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (ndims < ndim_exp) then
      readvar2D_nc = ndim_exp - ndims
    else if (ndims > ndim_exp) then
      print *, 'NetCDF Read error in file ' //trim(fname)
      print *, 'Too many dimensions (',ndims ,' instead of ',ndim_exp,')'
      print *, 'For variable ' // trim(ncname)
      call appl_finalize(-1)
    end if
    if (ndims > 0 ) then
      iret = nf90_inquire_variable (ncid,varid,dimids=dimids)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      do n = 1, ndims
        iret = nf90_inquire_dimension(ncid,dimids(n),name=dimname(n),len=dimsize(n))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
        do nn = 1, ndim_exp 
          if (dimname(n)(1:len_trim(dimname_in(nn))) == trim(dimname_in(nn))) then
            nslice(n) = nn
            exit 
          end if
        end do
      end do
    end if
    allocate(dimvars(maxval(dimsize),ndim_exp))
    dimvars = 0.
    do n=1,ndims
      iret  = nf90_inq_varid(ncid, validate(dimname(n)), dimid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_get_var(ncid,dimid, dimvars(1:dimsize(n), n))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      if (nslice(n) > 0) then
        locmax = dimsize(n)
        locmin = 1
        do i = 1, dimsize(n)
          if (dimvars(i, n) >= dimmax(nslice(n))) locmax = max(i, locmax)
          if (dimvars(i, n) <= dimmin(nslice(n))) locmin = min(i, locmin)              
        end do
        dimsize(n) = locmax - locmin + 1
        start(n)   = locmin
      end if
    end do
    
    if (allocated(dimvar1)) deallocate(dimvar1)
    allocate(dimvar1(dimsize(1)))
    dimvar1 = dimvars(start(1):dimsize(1),1)
    if (allocated(dimvar2)) deallocate(dimvar2)
    allocate(dimvar2(dimsize(2)))
    dimvar2 = dimvars(start(2):dimsize(2),2)
   
!Read the variable
    if (allocated(var    )) deallocate(var   )
    allocate(var(dimsize(1),dimsize(2)))
    iret  = nf90_get_var(ncid,varid,var, start = start, count = dimsize)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    iret  = nf90_close(ncid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
  end function readvar2D_nc

!> Reads out a variable that has up to 3 dimensions.
!! If the variable cannot be read, the netcdf errorcode will be passed on to the calling routine.
!! If the dimensionality is lower than 3, the return value will be equal to the number of missing dimensions.
  function readvar3D_nc(fname, ncname, dimname, dimvar1, dimvar2, dimvar3, var, lslice)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_open, nf90_close, nf90_get_var, nf90_inquire_variable, nf90_inquire_dimension, &
                       nf90_inq_varid, nf90_nowrite, nf90_noerr
    integer                                          :: readvar3D_nc !< Return value. Negative if variable cannot be read, positive if less than 3 dimensions available.
    character(*),intent(in)                          :: fname        !< File name
    character(*),intent(in)                          :: ncname       !< Variable name
    character(*),dimension(3),intent(inout)          :: dimname      !< Name of the dimension (in array form)
    real,allocatable,dimension(:),intent(inout)      :: dimvar1      !< Values of the 1st dimension
    real,allocatable,dimension(:),intent(inout)      :: dimvar2      !< Values of the 2nd dimension
    real,allocatable,dimension(:),intent(inout)      :: dimvar3      !< Values of the 3th dimension
    real,allocatable,dimension(:,:,:),intent(out)    :: var          !< Variable value
    logical, intent(in), optional                    :: lslice
    integer, parameter                     :: ndim_exp = 3
    integer                                :: ncid, varid, ndims, dimid, i, n, nn, iret, locmin, locmax
    integer, dimension(ndim_exp)           :: dimsize, dimids, start, nslice
    character(len=80), dimension(ndim_exp) :: dimname_in
    real, dimension(ndim_exp)              :: dimmax, dimmin
    real, allocatable, dimension(:,:)      :: dimvars
    logical                                :: exists

    readvar3D_nc = 0
    start         = 1
    dimsize       = 1
    dimmin        = -1
    dimmax        = dimmin - 1
    nslice        = 0
    dimname_in(:) = '########'
    if (present(lslice)) then
      if (lslice) then
        if (allocated(dimvar1)) then
          dimmax(1) = maxval(dimvar1)
          dimmin(1) = minval(dimvar1)
          dimname_in(1) = dimname(1)
        end if
        if (allocated(dimvar2)) then
          dimmax(2) = maxval(dimvar2)
          dimmin(2) = minval(dimvar2)
          dimname_in(2) = dimname(2)
        end if
        if (allocated(dimvar3)) then
          dimmax(3) = maxval(dimvar3)
          dimmin(3) = minval(dimvar3)
          dimname_in(3) = dimname(3)
        end if
      end if
    end if
    dimname = dimname_in
!Open file
    inquire(file=trim(fname),exist=exists)
    if (.not. exists) then
      readvar3D_nc = -1
      return
    end if
    iret  = nf90_open (trim(fname), NF90_NOWRITE, ncid)
    if (iret /= nf90_noerr) then
      readvar3D_nc = -1
      return
    end if
!Check whether variable is available
    iret  = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) then
      readvar3D_nc = iret
      iret         = nf90_close(ncid)
      return
    end if

!Check the dimensions of the variable
    iret = nf90_inquire_variable (ncid,varid,ndims=ndims)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
    if (ndims < ndim_exp) then
      readvar3D_nc = ndim_exp - ndims
    else if (ndims > ndim_exp) then
      print *, 'NetCDF Read error in file ' //trim(fname)
      print *, 'Too many dimensions (',ndims ,' instead of ',ndim_exp,')'
      print *, 'For variable ' // trim(ncname)
      call appl_finalize(-1)
    end if
    if (ndims > 0 ) then
     iret = nf90_inquire_variable (ncid,varid,dimids=dimids)
     if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
     do n = 1, ndims
        iret = nf90_inquire_dimension(ncid,dimids(n),name=dimname(n),len=dimsize(n))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
        do nn = 1, ndim_exp 
          if (dimname(n)(1:len_trim(dimname_in(nn))) == trim(dimname_in(nn))) then
            nslice(n) = nn
            exit 
          end if
        end do
      end do
    end if
    allocate(dimvars(maxval(dimsize),ndim_exp))
    dimvars = 0.
    do n=1,ndims
      iret  = nf90_inq_varid(ncid, validate(dimname(n)), dimid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_get_var(ncid,dimid, dimvars(1:dimsize(n), n))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      if (nslice(n) > 0) then
        locmax = dimsize(n)
        locmin = 1
        do i = 1, dimsize(n)
          if (dimvars(i, n) >= dimmax(nslice(n))) locmax = max(i, locmax)
          if (dimvars(i, n) <= dimmin(nslice(n))) locmin = min(i, locmin)              
        end do
        dimsize(n) = locmax - locmin + 1
        start(n)   = locmin
      end if
    end do
!Read the variable
    if (allocated(dimvar1)) deallocate(dimvar1)
    allocate(dimvar1(dimsize(1)))
    dimvar1 = dimvars(start(1):dimsize(1),1)
    if (allocated(dimvar2)) deallocate(dimvar2)
    allocate(dimvar2(dimsize(2)))
    dimvar2 = dimvars(start(2):dimsize(2),2)
    if (allocated(dimvar3)) deallocate(dimvar3)
    allocate(dimvar3(dimsize(3)))
    dimvar3 = dimvars(start(3):dimsize(3),3)
    if (allocated(var    )) deallocate(var   )
    allocate(var(dimsize(1),dimsize(2),dimsize(3)))
    iret  = nf90_get_var(ncid,varid,var, start = start, count = dimsize)
    iret  = nf90_close(ncid)
  end function readvar3D_nc

!> Reads out a variable that has up to 4 dimensions.
!! If the variable cannot be read, the netcdf errorcode will be passed on to the calling routine.
!! If the dimensionality is lower than 4, the return value will be equal to the number of missing dimensions.
  function readvar4D_nc(fname, ncname, dimname, dimvar1, dimvar2, dimvar3, dimvar4, var, lslice)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_open, nf90_close, nf90_get_var, nf90_inquire_variable, nf90_inquire_dimension, &
                       nf90_inq_varid, nf90_nowrite, nf90_noerr
    integer                                          :: readvar4D_nc !< Return value. Negative if variable cannot be read, positive if less than 4 dimensions available.
    character(*),intent(in)                          :: fname        !< File name
    character(*),intent(in)                          :: ncname       !< Variable name
    character(*),dimension(4),intent(inout)          :: dimname      !< Name of the dimension (in array form)
    real,allocatable,dimension(:),intent(inout)      :: dimvar1      !< Values of the 1st dimension
    real,allocatable,dimension(:),intent(inout)      :: dimvar2      !< Values of the 2nd dimension
    real,allocatable,dimension(:),intent(inout)      :: dimvar3      !< Values of the 3rd dimension
    real,allocatable,dimension(:),intent(inout)      :: dimvar4      !< Values of the 4th dimension
    real,allocatable,dimension(:,:,:,:),intent(out)  :: var          !< Variable value
    logical, intent(in), optional                    :: lslice
    integer, parameter                     :: ndim_exp = 4
    integer                                :: ncid, varid, ndims, dimid, i, n, nn, iret, locmin, locmax
    integer, dimension(ndim_exp)           :: dimsize, dimids, start, nslice
    character(len=80), dimension(ndim_exp) :: dimname_in
    real, dimension(ndim_exp)              :: dimmax, dimmin
    real, allocatable, dimension(:,:)      :: dimvars
    logical :: exists
    
    readvar4D_nc = 0
    
    start         = 1
    dimsize       = 1
    dimmin        = -1
    dimmax        = dimmin - 1
    nslice        = 0
    dimname_in(:) = '########'
    if (present(lslice)) then
      if (lslice) then
        if (allocated(dimvar1)) then
          dimmax(1) = maxval(dimvar1)
          dimmin(1) = minval(dimvar1)
          dimname_in(1) = dimname(1)
        end if
        if (allocated(dimvar2)) then
          dimmax(2) = maxval(dimvar2)
          dimmin(2) = minval(dimvar2)
          dimname_in(2) = dimname(2)
        end if
        if (allocated(dimvar3)) then
          dimmax(3) = maxval(dimvar3)
          dimmin(3) = minval(dimvar3)
          dimname_in(3) = dimname(3)
        end if
        if (allocated(dimvar4)) then
          dimmax(4) = maxval(dimvar4)
          dimmin(4) = minval(dimvar4)
          dimname_in(4) = dimname(4)
        end if
      end if
    end if
    dimname = dimname_in
!Open file
    inquire(file=trim(fname),exist=exists)
    if (.not. exists) then
      readvar4D_nc = -1
      return
    end if
    iret  = nf90_open (trim(fname), NF90_NOWRITE, ncid)
    if (iret /= nf90_noerr) then
      readvar4D_nc = -1
      return
    end if
!Check whether variable is available
    iret  = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) then
      readvar4D_nc = iret
      iret         = nf90_close(ncid)
      return
    end if

!Check the dimensions of the variable
    iret = nf90_inquire_variable (ncid,varid,ndims=ndims)
    if (ndims < ndim_exp) then
      readvar4D_nc = ndim_exp - ndims
    else if (ndims > ndim_exp) then
      print *, 'NetCDF Read error in file ' //trim(fname)
      print *, 'Too many dimensions (',ndims ,' instead of ',ndim_exp,')'
      print *, 'For variable ' // trim(ncname)
      call appl_finalize(-1)
    end if
    if (ndims > 0 ) then
      iret = nf90_inquire_variable (ncid,varid,dimids=dimids)
      do n = 1, ndims
        iret = nf90_inquire_dimension(ncid,dimids(n),name=dimname(n),len=dimsize(n))
        do nn = 1, ndim_exp 
          if (dimname(n)(1:len_trim(dimname_in(nn))) == trim(dimname_in(nn))) then
            nslice(n) = nn
            exit 
          end if
        end do
      end do
    end if
    allocate(dimvars(maxval(dimsize),ndim_exp))
    dimvars = 0.
    do n=1,ndims
      iret  = nf90_inq_varid(ncid, validate(dimname(n)), dimid)
      iret  = nf90_get_var(ncid,dimid, dimvars(1:dimsize(n), n))
      if (nslice(n) > 0) then
        locmax = dimsize(n)
        locmin = 1
        do i = 1, dimsize(n)
          if (dimvars(i, n) >= dimmax(nslice(n))) locmax = max(i, locmax)
          if (dimvars(i, n) <= dimmin(nslice(n))) locmin = min(i, locmin)              
        end do
        dimsize(n) = locmax - locmin + 1
        start(n)   = locmin
      end if
    end do
    
    if (allocated(dimvar1)) deallocate(dimvar1)
    allocate(dimvar1(dimsize(1)))
    dimvar1 = dimvars(start(1):dimsize(1),1)
    if (allocated(dimvar2)) deallocate(dimvar2)
    allocate(dimvar2(dimsize(2)))
    dimvar2 = dimvars(start(2):dimsize(2),2)
    if (allocated(dimvar3)) deallocate(dimvar3)
    allocate(dimvar3(dimsize(3)))
    dimvar3 = dimvars(start(3):dimsize(3),3)
    if (allocated(dimvar4)) deallocate(dimvar4)
    allocate(dimvar4(dimsize(4)))
    dimvar4 = dimvars(start(4):dimsize(4),4)
    
!Read the variable
    if (allocated(var    )) deallocate(var   )
    allocate(var(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
    iret  = nf90_get_var(ncid,varid,var, start = start, count = dimsize)
    iret  = nf90_close(ncid)
  end function readvar4D_nc
  
  integer function putatt_str_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_put_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    character(len=*), intent(in)           :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef
    
    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_str_nc = nf90_put_att(ncid, varid, validate(attrname), trim(attrval))
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_str_nc
  
  integer function putatt_single_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_put_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(kind=single), intent(in)          :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef
    
    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_single_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_single_nc
  
  integer function putatt_double_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_put_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(kind=double), intent(in)          :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef
    
    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_double_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_double_nc
  
  integer function putatt_int_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_put_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(kind=int32), intent(in)        :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef
    
    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_int_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_int_nc
  
  integer function putatt_short_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_put_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(kind=short), intent(in)        :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef
    
    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_short_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_short_nc
  
  integer function getatt_str_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_get_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    character(len=*), intent(out)          :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_str_nc = nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_str_nc
  
  integer function getatt_single_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_get_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(kind=single), intent(out)         :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_single_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_single_nc
  
  integer function getatt_double_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_get_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(kind=double), intent(out)         :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_double_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_double_nc
  
  integer function getatt_int_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_get_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(kind=int32), intent(out)       :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_int_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_int_nc
  
  integer function getatt_short_nc(ncid, attrname, attrval, varname)
    use netcdf, only : nf90_get_att, nf90_inq_varid, nf90_global
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(kind=short), intent(out)       :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_short_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_short_nc
  
  integer function getatt_str_fname_nc(fname, attrname, attrval, varname)
    use netcdf, only : nf90_open, nf90_get_att, nf90_inq_varid, nf90_global, nf90_nowrite
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: attrname
    character(len=*), intent(out)          :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid, ncid
    
    iret = nf90_open(trim(fname), NF90_NOWRITE, ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_str_fname_nc = nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_str_fname_nc
  
  function validate(input) !\todo Make this allocatable as soon as all common compilers allow it
    character(len=*), intent(in) :: input
    character(len=len_trim(input)):: validate
    character(len=2) :: cforbidden
    integer          :: n, nn
    cforbidden = '()'
    nn = 0
    do n=1,len_trim(input)
      validate(n:n) = ' '
      if (scan(input(n:n),cforbidden) == 0) then
        nn = nn + 1
        validate(nn:nn) = input(n:n)
      end if
    end do
  end function validate

!> Synchronize the netcdf file on disk with the buffer.
  integer function sync_nc(ncid)
    use netcdf, only : nf90_sync
    integer, intent(in) :: ncid            !< Netcdf file number
    sync_nc = 1
    if (lsync) sync_nc = nf90_sync(ncid)
  end function sync_nc

!> Write netcdf error to screen and stop the program
  subroutine nchandle_error(ncid, status)
    use mpi_interface, only : appl_finalize
    use netcdf, only : nf90_strerror, nf90_noerr
    integer, intent(in) :: ncid
    integer, intent(in) :: status !< NetCDF error code
    integer       :: iret
    character(80) :: fname
    if(status /= nf90_noerr) then
      iret = getatt_nc(ncid, 'title', fname)
      if (iret /= nf90_noerr) fname = ''
      print *, 'NetCDF error in file ' // trim(fname)
      print *, status, trim(nf90_strerror(status))
      call appl_finalize (-1)
    end if
  end subroutine nchandle_error

end module modnetcdf
