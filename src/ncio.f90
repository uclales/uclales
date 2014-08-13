 module ncio

  use netcdf
  use grid
  use mpi_interface, only : appl_abort, myid, pecount, wrxid, wryid
  use mcrp, only : cldw,rain,ice,snow,graupel,hail

  implicit none
  private

  public :: open_nc, define_nc, init_anal, close_anal, write_anal, deflate_level

  integer, private, save  :: nrec0, nvar0, nbase=15
  integer, save           :: ncid0,ncid_s, ncid_cross, deflate_level
  integer, save           :: crossx, crossy, crossz
  character (len=7), dimension(30) :: crossnames
  character (len=7),  private :: v_snm='sxx    '
  character (len=80), private :: fname

contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Open_NC: Opens a NetCDF File and identifies starting record
  !
  subroutine open_nc (fname, ename, time, npts, ncid, nrec, singleio)

    integer, intent(in)             :: npts
    integer, intent(out)            :: ncid
    integer, intent(out)            :: nrec
    real, intent (in)               :: time
    character (len=80), intent (in) :: fname, ename
    logical, intent(in), optional   :: singleio

    real, allocatable :: xtimes(:)

    character (len=8)  :: date
    character (len=88) :: lfname
    integer :: iret, ncall, VarID, RecordDimID
    logical :: exans
 
    if((present(singleio) .and. singleio .eqv. .true.) .or. pecount == 1) then
       write(lfname,'(a,a3)') trim(fname),'.nc'
    else
       write(lfname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',wrxid,wryid,'.nc'
    end if

    inquire(file=trim(lfname),exist=exans)

    ncall = 0
    if (.not.exans) then
       call date_and_time(date)
       iret = nf90_create(lfname,NF90_SHARE,ncid)

       iret = nf90_put_att(ncid,NF90_GLOBAL,'title',ename)
       iret = nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//date)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'Source','UCLA-LES Version 2.0')
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'Author','Bjorn Stevens')
      !iret = nf90_put_att(ncid, NF90_GLOBAL, '_FillValue',-999.)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'NPTS',npts)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'NPROCS',pecount)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'PROCID',myid)
    else
       iret = nf90_open (trim(lfname), NF90_WRITE, ncid)
       iret = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
       iret = nf90_inquire_dimension(ncid, RecordDimID, len=nrec)
       ncall=1
       iret = nf90_inq_varid(ncid,'time',VarID)
       allocate (xtimes(nrec+1))
       iret = nf90_get_var(ncid, VarId, xtimes(1:nrec))
       ncall = 1
       do while(ncall <= nrec .and. xtimes(ncall) < time - spacing(1.))
          ncall=ncall+1
       end do
       deallocate(xtimes)
    end if
    nrec = ncall
    iret = nf90_sync(ncid)

  end subroutine open_nc
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Define_NC: Defines the structure of the nc file (if not
  ! already open)
  !
  subroutine define_nc(ncID, nRec, nVar, sx, n1, n2, n3)

    integer, intent (in)           :: nVar, ncID
    integer, optional, intent (in) :: n1, n2, n3
    integer, intent (inout)        :: nRec
    character (len=7), intent (in) :: sx(nVar)   ! table with var names

    integer, save :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0,&
         zsoilID=0, dim_mttt(4)=0, dim_tmtt(4) = 0, dim_ttmt(4) = 0,          &
         dim_tttt(4) = 0, dim_tt(2)  = 0, dim_mt(2)  = 0, dim_mmt(3) = 0, &
         dim_soilmmt(4) = 0

    character (len=7) :: xnm
    integer :: iret, n, VarID
    !
    if (nRec == 0) then
       iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
       if (present(n1)) then
          iret = nf90_def_dim(ncID, 'zt', n1, ztID)
          iret = nf90_def_dim(ncID, 'zm', n1, zmID)
          iret = nf90_def_dim(ncID, 'zsoil', 4, zsoilID)
       end if
       if (present(n2)) then
          iret = nf90_def_dim(ncID, 'xt', n2, xtID)
          iret = nf90_def_dim(ncID, 'xm', n2, xmID)
       end if
       if (present(n3)) then
          iret = nf90_def_dim(ncID, 'yt', n3, ytID)
          iret = nf90_def_dim(ncID, 'ym', n3, ymID)
       end if
       dim_tt = (/ztId,timeId/)
       dim_mt = (/zmId,timeId/)
       dim_tttt= (/ztID,xtID,ytID,timeId/)       ! thermo point
       dim_mttt= (/zmID,xtID,ytID,timeId/)       ! wpoint
       dim_tmtt= (/ztID,xmID,ytID,timeId/)       ! upoint
       dim_ttmt= (/ztId,xtID,ymID,timeId/)       ! vpoint
       dim_mmt= (/xmID,ymID,timeId/)             ! srfc point
       dim_soilmmt= (/zsoilID,xmID,ymID,timeId/) ! soil point

       do n=1,nVar
          select case(trim(ncinfo(2,sx(n))))
          case ('time')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,timeId  ,VarID)
          case ('zt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,ztID    ,VarID)
          case ('zm')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,zmID    ,VarID)
          case ('xt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,xtID    ,VarID)
          case ('xm')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,xmID    ,VarID)
          case ('yt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,ytID    ,VarID)
          case ('ym')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,ymID    ,VarID)
          case ('tttt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttt,VarID)!, deflate_level = deflate_level)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
             if (sx(n).eq."l") then
                iret=nf90_put_att(ncID,VarID,'nu',cldw%nu)
                iret=nf90_put_att(ncID,VarID,'mu',cldw%mu)
                iret=nf90_put_att(ncID,VarID,'a_geo',cldw%a_geo)
                iret=nf90_put_att(ncID,VarID,'b_geo',cldw%b_geo)
                iret=nf90_put_att(ncID,VarID,'a_vel',cldw%a_vel)
                iret=nf90_put_att(ncID,VarID,'b_vel',cldw%b_vel)
                iret=nf90_put_att(ncID,VarID,'x_min',cldw%x_min)
                iret=nf90_put_att(ncID,VarID,'x_max',cldw%x_max)
             end if
             if (sx(n).eq."r") then
                iret=nf90_put_att(ncID,VarID,'nu',rain%nu)
                iret=nf90_put_att(ncID,VarID,'mu',rain%mu)
                iret=nf90_put_att(ncID,VarID,'a_geo',rain%a_geo)
                iret=nf90_put_att(ncID,VarID,'b_geo',rain%b_geo)
                iret=nf90_put_att(ncID,VarID,'a_vel',rain%a_vel)
                iret=nf90_put_att(ncID,VarID,'b_vel',rain%b_vel)
                iret=nf90_put_att(ncID,VarID,'x_min',rain%x_min)
                iret=nf90_put_att(ncID,VarID,'x_max',rain%x_max)
             end if
             if (sx(n).eq."rice") then
                iret=nf90_put_att(ncID,VarID,'nu',ice%nu)
                iret=nf90_put_att(ncID,VarID,'mu',ice%mu)
                iret=nf90_put_att(ncID,VarID,'a_geo',ice%a_geo)
                iret=nf90_put_att(ncID,VarID,'b_geo',ice%b_geo)
                iret=nf90_put_att(ncID,VarID,'a_vel',ice%a_vel)
                iret=nf90_put_att(ncID,VarID,'b_vel',ice%b_vel)
                iret=nf90_put_att(ncID,VarID,'x_min',ice%x_min)
                iret=nf90_put_att(ncID,VarID,'x_max',ice%x_max)
             end if
             if (sx(n).eq."rsnow") then
                iret=nf90_put_att(ncID,VarID,'nu',snow%nu)
                iret=nf90_put_att(ncID,VarID,'mu',snow%mu)
                iret=nf90_put_att(ncID,VarID,'a_geo',snow%a_geo)
                iret=nf90_put_att(ncID,VarID,'b_geo',snow%b_geo)
                iret=nf90_put_att(ncID,VarID,'a_vel',snow%a_vel)
                iret=nf90_put_att(ncID,VarID,'b_vel',snow%b_vel)
                iret=nf90_put_att(ncID,VarID,'x_min',snow%x_min)
                iret=nf90_put_att(ncID,VarID,'x_max',snow%x_max)
             end if
             if (sx(n).eq."rgrp") then
                iret=nf90_put_att(ncID,VarID,'nu',graupel%nu)
                iret=nf90_put_att(ncID,VarID,'mu',graupel%mu)
                iret=nf90_put_att(ncID,VarID,'a_geo',graupel%a_geo)
                iret=nf90_put_att(ncID,VarID,'b_geo',graupel%b_geo)
                iret=nf90_put_att(ncID,VarID,'a_vel',graupel%a_vel)
                iret=nf90_put_att(ncID,VarID,'b_vel',graupel%b_vel)
                iret=nf90_put_att(ncID,VarID,'x_min',graupel%x_min)
                iret=nf90_put_att(ncID,VarID,'x_max',graupel%x_max)
             end if
             if (sx(n).eq."rhail") then
                iret=nf90_put_att(ncID,VarID,'nu',hail%nu)
                iret=nf90_put_att(ncID,VarID,'mu',hail%mu)
                iret=nf90_put_att(ncID,VarID,'a_geo',hail%a_geo)
                iret=nf90_put_att(ncID,VarID,'b_geo',hail%b_geo)
                iret=nf90_put_att(ncID,VarID,'a_vel',hail%a_vel)
                iret=nf90_put_att(ncID,VarID,'b_vel',hail%b_vel)
                iret=nf90_put_att(ncID,VarID,'x_min',hail%x_min)
                iret=nf90_put_att(ncID,VarID,'x_max',hail%x_max)
             end if
          case ('mttt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mttt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('tmtt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tmtt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('ttmt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttmt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mt,VarID)
             end if
          case ('mmt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mmt,VarID)
          case ('soilmmt') 
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_soilmmt,VarID)
          case default
             if (myid == 0) print *, '  ABORTING: Bad dimensional information'
             call appl_abort(0)
          end select
          iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,sx(n)))
          iret=nf90_put_att(ncID,VarID,'_FillValue',-999._4)
          iret=nf90_put_att(ncID,VarID,'units'   ,ncinfo(1,sx(n)))
       end do
       iret  = nf90_enddef(ncID)

       iret  = nf90_sync(ncID)
       nRec = 1
    else
       iret = nf90_inquire(ncID, nVariables=n)
       if (n /= nVar) then
          iret = nf90_close(ncID)
          if (myid == 0) print *, '  ABORTING: Incompatible Netcdf File',n,nVar
          call appl_abort(0)
       else
          do n=1,nVar
             xnm=sx(n)
             iret = nf90_inquire_variable(ncID, n, name=xnm)
          end do
          iret = nf90_sync(ncID)
       end if
    end if

  end subroutine define_nc
  !
  ! ----------------------------------------------------------------------
  ! subroutine init_anal:  Defines the netcdf Analysis file
  !
  subroutine init_anal(time)

    use mpi_interface, only :myid

    integer, parameter :: nnames = 37
    character (len=7), save :: sbase(nnames) =  (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     '   ,& !1
         'ym     ','u0     ','v0     ','dn0    ','u      ','v      '   ,& !7 
         'w      ','t      ','p      ','q      ','l      ','r      '   ,& !13
         'n      ','rice   ','nice   ','rsnow  ','rgrp   ','nsnow  '   ,& !19
         'ngrp   ','rhail  ','nhail  ','rflx   ','lflxu  ','lflxd  '   ,& !25
         'mp_tlt ','mp_qt  ','mp_qr  ','mp_qi  ','mp_qs  ','mp_qg  '   ,& !31
         'mp_qh  '/)  !37

    real, intent (in) :: time
    integer           :: nbeg, nend

    nvar0 = nbase + naddsc
    if (level  >= 1)                   nvar0 = nvar0+1
    if (level  >= 2)                   nvar0 = nvar0+1
    if (level  >= 3)                   nvar0 = nvar0+2
    if (level  >= 4)                   nvar0 = nvar0+4
    if (level  >= 5)                   nvar0 = nvar0+4
    if (iradtyp > 1 .and. iradtyp < 5) nvar0 = nvar0+3
    if (lmptend)                       nvar0 = nvar0+7

    allocate (sanal(nvar0))
    sanal(1:nbase) = sbase(1:nbase)

    nvar0 = nbase
    !
    ! add additional scalars, in the order in which they appear in scalar
    ! table
    !
    if (level >= 1) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(nbase+1)
    end if
    !
    ! add liquid water, which is a diagnostic variable, first
    !
    if (level >= 2) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(nbase+2)
    end if

    if (level >= 3) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(nbase+3)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(nbase+4)
    end if

    if (level >= 4) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(20)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(21)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(22)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(23)
    end if

    if (level >= 5) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(24)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(25)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(26)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(27)
    end if

    if (iradtyp > 1 .and. iradtyp < 5) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(28)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(29)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(30)
    end if

    if (lmptend) then
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(31)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(32)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(33)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(34)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(35)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(36)
       nvar0 = nvar0+1
       sanal(nvar0) = sbase(37)
    end if

    nbeg = nvar0+1
    nend = nvar0+naddsc
    do nvar0 = nbeg, nend
       write(v_snm(2:3),'(i2.2)') nvar0-nbeg
       sanal(nvar0) = v_snm
    end do
    nvar0=nend
    fname =  trim(filprf)
    if(myid == 0) print                                                  &
            "(//' ',49('-')/,' ',/,'   Initializing: ',A20)",trim(fname)
    call open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid0, nrec0)
    call define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4)
    if (myid == 0) print *,'   ...starting record: ', nrec0

  end subroutine init_anal
  !
  ! ----------------------------------------------------------------------
  ! subroutine close_anal:  Closes netcdf anal file
  !
  integer function close_anal()

    use netcdf

    close_anal = nf90_close(ncid0)

  end function close_anal
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Write_anal:  Writes the netcdf Analysis file
  !
  subroutine write_anal(time)

    use netcdf
    use mpi_interface, only : myid, appl_abort
    use defs, only : cp, alvl

    real, intent (in) :: time

    integer :: iret, VarID, nn, n
    integer :: ibeg(4), icnt(4), i1, i2, j1, j2
    integer :: icntsfc(3),icntsoil(4),ibegsfc(3)

    !return
    icnt = (/nzp,nxp-4,nyp-4,1/)
    icntsfc = (/nxp-4,nyp-4,1/)
    icntsoil = (/4,nxp-4,nyp-4,1/)
    ibeg = (/1  ,1  ,1  ,nrec0/)
    ibegsfc = (/1,1,nrec0/)
    i1 = 3
    i2 = nxp-2
    j1 = 3
    j2 = nyp-2
    iret = nf90_inq_Varid(ncid0, sanal(1), VarID)
    iret = nf90_put_var(ncid0, VarID, time, start=(/nrec0/))
    if (nrec0 == 1) then
       iret = nf90_inq_varid(ncid0, sanal(2), VarID)
       iret = nf90_put_var(ncid0, VarID, zt, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(3), VarID)
       iret = nf90_put_var(ncid0, VarID, zm, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(4), VarID)
       iret = nf90_put_var(ncid0, VarID, xt(i1:i2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(5), VarID)
       iret = nf90_put_var(ncid0, VarID, xm(i1:i2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(6), VarID)
       iret = nf90_put_var(ncid0, VarID, yt(j1:j2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(7), VarID)
       iret = nf90_put_var(ncid0, VarID, ym(j1:j2), start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(8), VarID)
       iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(9), VarID)
       iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
       iret = nf90_inq_varid(ncid0, sanal(10), VarID)
       iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))
    end if
    iret = nf90_inq_varid(ncid0, sanal(11), VarID)
    iret = nf90_put_var(ncid0, VarID, a_up(:,i1:i2,j1:j2), start=ibeg,    &
         count=icnt)
    iret = nf90_inq_varid(ncid0, sanal(12), VarID)
    iret = nf90_put_var(ncid0, VarID, a_vp(:,i1:i2,j1:j2), start=ibeg,    &
         count=icnt)
    iret = nf90_inq_varid(ncid0, sanal(13), VarID)
    iret = nf90_put_var(ncid0, VarID, a_wp(:,i1:i2,j1:j2), start=ibeg,    &
         count=icnt)
    iret = nf90_inq_varid(ncid0, sanal(14), VarID)
    iret = nf90_put_var(ncid0, VarID, a_tp(:,i1:i2,j1:j2)+th00, start=ibeg, &
         count=icnt)
    iret = nf90_inq_varid(ncid0, sanal(15), VarID)
    iret = nf90_put_var(ncid0, VarID, press(:,i1:i2,j1:j2), start=ibeg, &
         count=icnt)

    if (level > 0) then
       iret = nf90_inq_varid(ncid0, sanal(16), VarID)
       iret = nf90_put_var(ncid0, VarID, a_rp(:,i1:i2,j1:j2), start=ibeg, &
           count=icnt)
    end if 

    if (level >= 2)  then
       iret = nf90_inq_varid(ncid0, sanal(17), VarID)
       iret = nf90_put_var(ncid0, VarID, liquid(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)
    end if
    nn = nbase+2

    if (level >=3) then
      do n = nbase+2, 18
        nn = nn+1
        call newvar(nn-12)
        iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
        iret = nf90_put_var(ncid0,VarID,a_sp(:,i1:i2,j1:j2), start=ibeg,   &
               count=icnt)
         if (myid==0) print*,"sanal(nn):",sanal(nn),nn
       end do
    endif  

    if (level >=4) then
      do n = 20, 23
        nn = nn+1
        call newvar(nn-12)
        iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
        iret = nf90_put_var(ncid0,VarID,a_sp(:,i1:i2,j1:j2), start=ibeg,   &
               count=icnt)
        if (myid==0) print*,"sanal(nn):",sanal(nn),nn
      end do
    endif  

    if (level >=5) then
      do n = 24, 27
        nn = nn+1
        call newvar(nn-12)
        iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
        iret = nf90_put_var(ncid0,VarID,a_sp(:,i1:i2,j1:j2), start=ibeg,   &
             count=icnt)
        if (myid==0) print*,"sanal(nn):",sanal(nn),nn
      end do
    endif  

    if (iradtyp > 1 .and. iradtyp < 5)  then
       nn = nn+1
       iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
       iret = nf90_put_var(ncid0, VarID, a_rflx(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)  
       nn = nn+1
       iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
       iret = nf90_put_var(ncid0, VarID, a_lflxu(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)
       nn = nn+1
       iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
       iret = nf90_put_var(ncid0, VarID, a_lflxd(:,i1:i2,j1:j2), start=ibeg, &
            count=icnt)
    end if

    if (lmptend)  then
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_tlt(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_qt(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_qr(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_qi(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_qs(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_qg(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
         nn = nn+1
         iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
         iret = nf90_put_var(ncid0, VarID, mp_qh(:,i1:i2,j1:j2), start=ibeg, &
              count=icnt)
      end if

!     if (nn /= nvar0) then
!        if (myid == 0) print *, 'ABORTING:  Anal write error'
!        call appl_abort(0)
!     end if

    if (myid==0) print "(//' ',12('-'),'   Record ',I3,' to: ',A60)",    &
         nrec0,fname

    iret  = nf90_sync(ncid0)
    nrec0 = nrec0+1

  end subroutine write_anal

  !
  ! ----------------------------------------------------------------------
  ! Subroutine nc_info: Gets long_name, units and dimension info given a
  ! short name.
  !
  character (len=80) function ncinfo(itype,short_name)

    character (len=40) :: v_lnm ='scalar xx mixing ratio                  '

    integer, intent (in) :: itype
    character (len=*), intent (in) :: short_name

    integer :: scalar_number

    select case (trim(short_name))
    case ('sxx')
       read (short_name(2:3),'(i2.2)') scalar_number
       write(v_lnm(8:9),'(i2.2)') scalar_number
       if (itype==0) ncinfo = v_lnm
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('time')
       if (itype==0) ncinfo = 'Time'
       if (itype==1) ncinfo = 'seconds since 2000-00-00 0000'
       if (itype==2) ncinfo = 'time'
    case('zt')
       if (itype==0) ncinfo = 'Vertical displacement of cell centers'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'zt'
    case('zm')
       if (itype==0) ncinfo = 'Vertical displacement of cell edges'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'zm'
    case('xt')
       if (itype==0) ncinfo = 'East-west displacement of cell centers'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'xt'
    case('xm')
       if (itype==0) ncinfo = 'East-west displacement of cell edges'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'xm'
    case('yt')
       if (itype==0) ncinfo = 'North-south displacement of cell centers'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'yt'
    case('ym')
       if (itype==0) ncinfo = 'North-south displacement of cell edges'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ym'
    case('u0')
       if (itype==0) ncinfo = 'Geostrophic zonal wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'zt'
    case('v0')
       if (itype==0) ncinfo = 'Geostrophic meridional wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'zt'
    case('dn0')
       if (itype==0) ncinfo = 'Base-state density'
       if (itype==1) ncinfo = 'kg/m^3'
       if (itype==2) ncinfo = 'zt'
    case('u')
       if (itype==0) ncinfo = 'Zonal wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'tmtt'
    case('v')
       if (itype==0) ncinfo = 'Meridional wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('w')
       if (itype==0) ncinfo = 'Vertical velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'mttt'
    case('t')
       if (itype==0) ncinfo = 'Liquid Water Potential temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('p')
       if (itype==0) ncinfo = 'Pressure'
       if (itype==1) ncinfo = 'Pa'
       if (itype==2) ncinfo = 'tttt'
    case('q')
       if (itype==0) ncinfo = 'Total water mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('l')
       if (itype==0) ncinfo = 'Liquid water mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('r')
       if (itype==0) ncinfo = 'Rain-water mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('RH')
       if (itype==0) ncinfo = 'Relative Humidity'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('inuc')
       if (itype==0) ncinfo = 'Number of ice nuclei'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rice')
       if (itype==0) ncinfo = 'Cloud Ice mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('nice')
       if (itype==0) ncinfo = 'Number of ice particles'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rsnow')
       if (itype==0) ncinfo = 'Snow mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('nsnow')
       if (itype==0) ncinfo = 'Number of snow particles'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rgrp')
       if (itype==0) ncinfo = 'Graupel mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('ngrp')
       if (itype==0) ncinfo = 'Number of graupel particles'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rhail')
       if (itype==0) ncinfo = 'Hail mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('nhail')
       if (itype==0) ncinfo = 'Number of hail particles'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rsup')
       if (itype==0) ncinfo = 'Supersaturation wrt ice'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('n')
       if (itype==0) ncinfo = 'Rain-drop number mixing ratio'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('cfl')
       if (itype==0) ncinfo = 'Courant number'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('maxdiv')
       if (itype==0) ncinfo = 'Maximum divergence'
       if (itype==1) ncinfo = '1/s'
       if (itype==2) ncinfo = 'time'
    case('zi1_bar')
       if (itype==0) ncinfo = 'Height of maximum theta gradient'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zi2_bar')
       if (itype==0) ncinfo = 'Height of maximum theta variance'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zi3_bar')
       if (itype==0) ncinfo = 'Height of minimum buoyancy flux'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('vtke')
       if (itype==0) ncinfo = 'Vertical integral of total TKE'
       if (itype==1) ncinfo = 'kg/s'
       if (itype==2) ncinfo = 'time'
       !irina
    case('tkeint')
       if (itype==0) ncinfo = 'Vertical integral of total TKE non-weighted'
       if (itype==1) ncinfo = 'm3/s2'
       if (itype==2) ncinfo = 'time'
    case('sfcbflx')
       if (itype==0) ncinfo = 'Surface Buoyancy Flux'
       if (itype==1) ncinfo = 'm/s^2'
       if (itype==2) ncinfo = 'time'
    case('wmax')
       if (itype==0) ncinfo = 'Maximum vertical velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'time'
    case('tsrf')
       if (itype==0) ncinfo = 'Surface temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('tsair')
       if (itype==0) ncinfo = 'Surface air temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('ustar')
       if (itype==0) ncinfo = 'Surface friction velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'time'
    case('shf_bar')
       if (itype==0) ncinfo = 'Sensible heat flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lhf_bar')
       if (itype==0) ncinfo = 'Latent heat flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('zi_bar')
       if (itype==0) ncinfo = 'Height of maximum scalar gradient'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('lwp_bar')
       if (itype==0) ncinfo = 'Liquid-water path'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('lwp_var')
       if (itype==0) ncinfo = 'Liquid-water path variance'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('zc')
       if (itype==0) ncinfo = 'Cloud-top height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zb')
       if (itype==0) ncinfo = 'Cloud-base height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zcmn')
       if (itype==0) ncinfo = 'Mean Cloud-top height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zbmn')
       if (itype==0) ncinfo = 'Mean Cloud-base height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('cfrac')
       if (itype==0) ncinfo = 'Cloud fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('lmax')
       if (itype==0) ncinfo = 'Maximum liquid water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('albedo')
       if (itype==0) ncinfo = 'Reflected (TOA) shortwave radiation'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('rwp_bar')
       if (itype==0) ncinfo = 'Rain-water path'
       !if (itype==1) ncinfo = 'kg/m^2'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       if (itype==2) ncinfo = 'time'
    case('prcp')
       if (itype==0) ncinfo = 'Surface precipitation rate'
       !irina
       if (itype==1) ncinfo = 'kg/kg m/s'
       !if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('pfrac')
       if (itype==0) ncinfo = 'Precipitation fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('CCN')
       if (itype==0) ncinfo = 'Cloud condensation nuclei'
       if (itype==1) ncinfo = '# per cc'
       if (itype==2) ncinfo = 'time'
    case('nrain')
       if (itype==0) ncinfo = 'Conditionally sampled rain number mixing ratio'
       if (itype==1) ncinfo = '# per liter'
       if (itype==2) ncinfo = 'time'
    case('nrcnt')
       if (itype==0) ncinfo = 'Rain cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('fsttm')
       if (itype==0) ncinfo = 'First sample time'
       if (itype==1) ncinfo = 's'
       if (itype==2) ncinfo = 'time'
    case('lsttm')
       if (itype==0) ncinfo = 'Last sample time'
       if (itype==1) ncinfo = 's'
       if (itype==2) ncinfo = 'time'
    case('nsmp')
       if (itype==0) ncinfo = 'Sample time counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('u_2')
       if (itype==0) ncinfo = 'Variance of u wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'tttt'
    case('v_2')
       if (itype==0) ncinfo = 'Variance of v wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'tttt'
    case('w_2')
       if (itype==0) ncinfo = 'Variance of w wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('t_2')
       if (itype==0) ncinfo = 'Variance of theta'
       if (itype==1) ncinfo = 'K^2'
       if (itype==2) ncinfo = 'tttt'
    case('w_3')
       if (itype==0) ncinfo = 'Third moment of w wind'
       if (itype==1) ncinfo = 'm^3/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('t_3')
       if (itype==0) ncinfo = 'Third moment of theta'
       if (itype==1) ncinfo = 'K^3'
       if (itype==2) ncinfo = 'tttt'
    case('tot_tw')
       if (itype==0) ncinfo = 'Total vertical flux of theta'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_tw')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of theta'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_tvw')
       if (itype==0) ncinfo = 'Total buoyancy flux'
       if (itype==1) ncinfo = 'Km/s'
       if (itype==2) ncinfo = 'ttmt'
    case('sgs_tvw')
       if (itype==0) ncinfo = 'Sub-filter scale buoyancy flux'
       if (itype==1) ncinfo = 'Km/s'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_uw')
       if (itype==0) ncinfo = 'Total vertical flux of u-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_uw')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of u-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_vw')
       if (itype==0) ncinfo = 'Total vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_vw')
       if (itype==0) ncinfo = 'SGS vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_ww')
       if (itype==0) ncinfo = 'Total vertical flux of w-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_ww')
       if (itype==0) ncinfo = 'SGS vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('km')
       if (itype==0) ncinfo = 'Eddy viscosity'
       if (itype==1) ncinfo = 'm^2/s'
       if (itype==2) ncinfo = 'ttmt'
    case('kh')
       if (itype==0) ncinfo = 'Eddy diffusivity'
       if (itype==1) ncinfo = 'm^2/s'
       if (itype==2) ncinfo = 'ttmt'
    case('lmbd')
       if (itype==0) ncinfo = 'Mixing lengthscale'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttmt'
    case('lmbde')
       if (itype==0) ncinfo = 'Dissipation lengthscale'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttmt'
    case('csmago')
       if (itype==0) ncinfo = 'Smagorinsky constant'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_tke')
       if (itype==0) ncinfo = 'Sub-filter scale TKE'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_boy')
       if (itype==0) ncinfo = 'Subfilter Buoyancy production of TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_shr')
       if (itype==0) ncinfo = 'Shear production of SGS TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('boy_prd')
       if (itype==0) ncinfo = 'Buoyancy production of resolved TKE '
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('shr_prd')
       if (itype==0) ncinfo = 'Shear production of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('trans')
       if (itype==0) ncinfo = 'Net transport of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('diss')
       if (itype==0) ncinfo = 'Dissipation rate of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('dff_u')
       if (itype==0) ncinfo = 'u(du/dt) from diffusion'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('dff_v')
       if (itype==0) ncinfo = 'v(dv/dt) from diffusion'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('dff_w')
       if (itype==0) ncinfo = 'w(dw/dt) from diffusion'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('adv_u')
       if (itype==0) ncinfo = 'u(du/dt) from advection'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('adv_v')
       if (itype==0) ncinfo = 'v(dv/dt) from advection'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('adv_w')
       if (itype==0) ncinfo = 'w(dw/dt) from advection'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('prs_u')
       if (itype==0) ncinfo = 'u(du/dt) from pressure'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('prs_v')
       if (itype==0) ncinfo = 'v(dv/dt) from pressure'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('prs_w')
       if (itype==0) ncinfo = 'w(dw/dt) from pressure'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('prd_uw')
       if (itype==0) ncinfo = 'uw shear production'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('storage')
       if (itype==0) ncinfo = 'Rate of increase of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('q_2')
       if (itype==0) ncinfo = 'Variance of total water'
       if (itype==1) ncinfo = 'kg^2/kg^2'
       if (itype==2) ncinfo = 'tttt'
    case('q_3')
       if (itype==0) ncinfo = 'Third moment of total water'
       if (itype==1) ncinfo = 'kg^3/kg^3'
       if (itype==2) ncinfo = 'tttt'
    case('tot_qw')
       if (itype==0) ncinfo = 'Total vertical flux of q'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_qw')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of q'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('rflx')
       if (itype==0) ncinfo = 'Total Radiative flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
       !irina
    case('lflxu')
       if (itype==0) ncinfo = 'Longwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('lflxd')
       if (itype==0) ncinfo = 'Longwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('mp_tlt')
       if (itype==0) ncinfo = 'T tendency from microphysics'
       if (itype==1) ncinfo = 'K/s'
       if (itype==2) ncinfo = 'tttt'
    case('mp_qt')
       if (itype==0) ncinfo = 'Qt tendency from microphysics'
       if (itype==1) ncinfo = 'kg/kg/s'
       if (itype==2) ncinfo = 'tttt'
    case('mp_qr')
       if (itype==0) ncinfo = 'Qr tendency from microphysics'
       if (itype==1) ncinfo = 'kg/kg/s'
       if (itype==2) ncinfo = 'tttt'
    case('mp_qi')
       if (itype==0) ncinfo = 'Qi tendency from microphysics'
       if (itype==1) ncinfo = 'kg/kg/s'
       if (itype==2) ncinfo = 'tttt'
    case('mp_qs')
       if (itype==0) ncinfo = 'Qs tendency from microphysics'
       if (itype==1) ncinfo = 'kg/kg/s'
       if (itype==2) ncinfo = 'tttt'
    case('mp_qg')
       if (itype==0) ncinfo = 'Qg tendency from microphysics'
       if (itype==1) ncinfo = 'kg/kg/s'
       if (itype==2) ncinfo = 'tttt'
    case('mp_qh')
       if (itype==0) ncinfo = 'Qh tendency from microphysics'
       if (itype==1) ncinfo = 'kg/kg/s'
       if (itype==2) ncinfo = 'tttt'
    case('lwuca')
       if (itype==0) ncinfo = 'Clear Air Longwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('lwdca')
       if (itype==0) ncinfo = 'Clear Air Longwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('lflxdt')
       if (itype==0) ncinfo =  'Top of Atmosphere Longwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lflxut')
       if (itype==0) ncinfo = 'Top of Atmosphere Longwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lflxutc')
       if (itype==0) ncinfo =  'Clear AirTop of Atmosphere Longwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lflxds')
       if (itype==0) ncinfo =  'Surface Longwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lflxus')
       if (itype==0) ncinfo =  'Surface Longwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lflxdsc')
       if (itype==0) ncinfo =  'Clear AirSurface Longwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lflxusc')
       if (itype==0) ncinfo =  'Clear AirSurface Longwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('rflx2')
       if (itype==0) ncinfo = 'Variance of total radiative flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('sflx')
       if (itype==0) ncinfo = 'Shortwave radiative flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
       !irina
    case('sflxu')
       if (itype==0) ncinfo = 'Shortwave radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('sflxd')
       if (itype==0) ncinfo = 'Shortwave radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('swuca')
       if (itype==0) ncinfo = 'Clear Air Shortwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('swdca')
       if (itype==0) ncinfo = 'Clear Air Shortwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'
    case('sflxut')
       if (itype==0) ncinfo = 'Top of Atmosphere Shortwave radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflxdt')
       if (itype==0) ncinfo = 'Top of Atmosphere Shortwave radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflxutc')
       if (itype==0) ncinfo = 'Clear Air Top of Atmosphere Shortwave radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflxds')
       if (itype==0) ncinfo =  'Surface Shortwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflxus')
       if (itype==0) ncinfo =  'Surface Shortwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflxdsc')
       if (itype==0) ncinfo =  'Clear AirSurface Shortwave Radiative flux DW'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflxusc')
       if (itype==0) ncinfo =  'Clear AirSurface Shortwave Radiative flux UP'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sflx2')
       if (itype==0) ncinfo = 'Variance of shortwave radiative flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'tttt'

    case('l_2')
       if (itype==0) ncinfo = 'Variance of liquid'
       if (itype==1) ncinfo = 'kg^2/kg^2'
       if (itype==2) ncinfo = 'tttt'
    case('l_3')
       if (itype==0) ncinfo = 'Third moment of liquid'
       if (itype==1) ncinfo = 'kg^3/kg^3'
       if (itype==2) ncinfo = 'tttt'
    case('tot_lw')
       if (itype==0) ncinfo = 'Resolved turbulent flux of liquid'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sed_lw')
       if (itype==0) ncinfo = 'Sedimentation flux of r_l'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('cs1')
       if (itype==0) ncinfo = 'Conditionally sampled fraction of flow'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('cnt_cs1')
       if (itype==0) ncinfo = 'Sum of I_cs1'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'tttt'
    case('w_cs1')
       if (itype==0) ncinfo = 'Conditional average of w over cs1'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttmt'
    case('tl_cs1')
       if (itype==0) ncinfo = 'Conditional average of theta_l over cs1'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'ttmt'
    case('tv_cs1')
       if (itype==0) ncinfo = 'Conditional average of theta_v over cs1'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('rt_cs1')
       if (itype==0) ncinfo = 'Conditional average of rt over cs1'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rl_cs1')
       if (itype==0) ncinfo = 'Conditional average of rl over cs1'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('wt_cs1')
       if (itype==0) ncinfo = 'Covariance of wtheta_l flux and cs1'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('wv_cs1')
       if (itype==0) ncinfo = 'Covariance of wtheta_v flux and cs1'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('wr_cs1')
       if (itype==0) ncinfo = 'Covariance of wr_t flux and cs1'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('cs2')
       if (itype==0) ncinfo = 'Conditionally sampled fraction of flow'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('cnt_cs2')
       if (itype==0) ncinfo = 'Sum of I_cs2'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'tttt'
    case('w_cs2')
       if (itype==0) ncinfo = 'Conditional average of w over cs2'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttmt'
    case('tl_cs2')
       if (itype==0) ncinfo = 'Conditional average of theta_l over cs2'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('tv_cs2')
       if (itype==0) ncinfo = 'Conditional average of theta_v over cs2'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('rt_cs2')
       if (itype==0) ncinfo = 'Conditional average of rt over cs2'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rl_cs2')
       if (itype==0) ncinfo = 'Conditional average of rl over cs2'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('wt_cs2')
       if (itype==0) ncinfo = 'Covariance of wtheta_l flux and cs2'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('wv_cs2')
       if (itype==0) ncinfo = 'Covariance of wtheta_v flux and cs2'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('wr_cs2')
       if (itype==0) ncinfo = 'Covariance of wr_t flux and cs2'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('Nc')
       if (itype==0) ncinfo = 'Cloud Number Concentration'
       if (itype==1) ncinfo = 'dm^-3'
       if (itype==2) ncinfo = 'tttt'
    case('Nr')
       if (itype==0) ncinfo = 'Rain Number Concentration'
       if (itype==1) ncinfo = 'dm^-3'
       if (itype==2) ncinfo = 'tttt'
    case('rr')
       if (itype==0) ncinfo = 'Rain water'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('i_nuc')
       if (itype==0) ncinfo = 'Ice nuclei Concentration'
       if (itype==1) ncinfo = 'dm^-3'
       if (itype==2) ncinfo = 'tttt'
    case('n_ice')
       if (itype==0) ncinfo = 'Ice Number Concentration'
       if (itype==1) ncinfo = 'dm^-3'
       if (itype==2) ncinfo = 'tttt'
    case('ice')
       if (itype==0) ncinfo = 'Cloud ice'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('snow')
       if (itype==0) ncinfo = 'Snow'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('graupel')
       if (itype==0) ncinfo = 'Graupel'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('hail')
       if (itype==0) ncinfo = 'hail mixing ratio'
       if (itype==1) ncinfo = 'g/kg'
       if (itype==2) ncinfo = 'tttt'
    case('prc_c')
       if (itype==0) ncinfo = 'Cloud water Precipitation Flux (positive downward)'
       !irina
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_r')
       if (itype==0) ncinfo = 'Rain Precipitation Flux (positive downward)'
       !irina
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_i')
       if (itype==0) ncinfo = 'Ice Precipitation Flux (positive downward)'
       !irina
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_s')
       if (itype==0) ncinfo = 'Snow Precipitation Flux (positive downward)'
       !irina
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_g')
       if (itype==0) ncinfo = 'Graupel Precipitation Flux (positive downward)'
       !irina
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_h')
       if (itype==0) ncinfo = 'Hail Precipitation Flux (positive downward)'
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('evap')
       if (itype==0) ncinfo = 'Net evap  of rain-water'
       if (itype==1) ncinfo = 's^-1'
       if (itype==2) ncinfo = 'tttt'
    case('frc_prc')
       if (itype==0) ncinfo = 'Conditionally sampled rain fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_prc')
       if (itype==0) ncinfo = 'Conditionally sampled rain rate'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('frc_ran')
       if (itype==0) ncinfo = 'Rain water fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('hst_srf')
       if (itype==0) ncinfo = 'Histogram of surface rain rates'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
       !irina
    case('cdsed')
       if (itype==0) ncinfo = 'Sedimentation Flux (positive downward)'
       if (itype==1) ncinfo = 'kg/kg m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wvp_bar')
       if (itype==0) ncinfo = 'Water vapor path'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('wvp_var')
       if (itype==0) ncinfo = 'Water vapor path variance'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('iwp_bar')
       if (itype==0) ncinfo = 'Cloud ice path'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('iwp_var')
       if (itype==0) ncinfo = 'Cloud ice path variance'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('swp_bar')
       if (itype==0) ncinfo = 'Snow path'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('swp_var')
       if (itype==0) ncinfo = 'Snow path variance'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('gwp_bar')
       if (itype==0) ncinfo = 'Graupel path'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('gwp_var')
       if (itype==0) ncinfo = 'Graupel path variance'
       !irina
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('hwp_bar')
       if (itype==0) ncinfo = 'Hail path'
       ! units may be wrong here
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('hwp_var')
       if (itype==0) ncinfo = 'Hail path variance'
       ! units may be wrong here
       if (itype==1) ncinfo = 'g/m^2'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('thl_int')
       if (itype==0) ncinfo = 'Integrated theta_l'
       !irina
       if (itype==1) ncinfo = 'Km'
       !if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('qt_th')    
       if (itype==0) ncinfo = 'Covariance of total water and liquid water potential temperature'
       if (itype==1) ncinfo = 'kg/kg K'
       if (itype==2) ncinfo = 'tttt'
    case('s_1')    
       if (itype==0) ncinfo = 'Mean of s (extended liquid water specific humidity)'
       if (itype==1) ncinfo = ''
       if (itype==2) ncinfo = 'tttt'
    case('s_2')    
       if (itype==0) ncinfo = 'Variance of s (extended liquid water specific humidity)'
       if (itype==1) ncinfo = ''
       if (itype==2) ncinfo = 'tttt'
    case('s_3')    
       if (itype==0) ncinfo = 'Third moment of s (extended liquid water specific humidity)'
       if (itype==1) ncinfo = ''
       if (itype==2) ncinfo = 'tttt'
    case('Qnet')    
       if (itype==0) ncinfo = 'Qnet'
       if (itype==1) ncinfo = 'w/m^2'
       if (itype==2) ncinfo = 'time'
    case('G0')    
       if (itype==0) ncinfo = 'Ground heat flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('tndskin')    
       if (itype==0) ncinfo = 'Skin heat tendency'
       if (itype==1) ncinfo = 'w/m^2'
       if (itype==2) ncinfo = 'time'
    case('ra')    
       if (itype==0) ncinfo = 'Aerodynamic resistance'
       if (itype==1) ncinfo = 's/m'
       if (itype==2) ncinfo = 'time'
    case('rsurf')    
       if (itype==0) ncinfo = 'Surface resistance'
       if (itype==1) ncinfo = 's/m'
       if (itype==2) ncinfo = 'time'
    case('rsveg')    
       if (itype==0) ncinfo = 'Vegetation resistance'
       if (itype==1) ncinfo = 's/m'
       if (itype==2) ncinfo = 'time'
    case('rssoil')    
       if (itype==0) ncinfo = 'Soil resistance'
       if (itype==1) ncinfo = 's/m'
       if (itype==2) ncinfo = 'time'
    case('tskinav')    
       if (itype==0) ncinfo = 'Average skin potential temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('qskinav')    
       if (itype==0) ncinfo = 'Average skin specific humidity'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('obl')    
       if (itype==0) ncinfo = 'Monin Obukhov Length'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('cliq')    
       if (itype==0) ncinfo = 'Fraction of vegetated surface covered with liquid water'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('a_Wl')    
       if (itype==0) ncinfo = 'Liquid water reservoir'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('shf')    
       if (itype==0) ncinfo = 'Surface Sensible Heat Flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'mmt'
    case('lhf')    
       if (itype==0) ncinfo = 'Surface Latent Heat Flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'mmt'
    case('ustars')    
       if (itype==0) ncinfo = 'Surface Friction Velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'mmt'
    case('a_tskin')    
       if (itype==0) ncinfo = 'Surface Skin Temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'mmt'
    case('a_qskin')    
       if (itype==0) ncinfo = 'Surface Skin Specific Humidity'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'mmt'
    case('tsoil')    
       if (itype==0) ncinfo = 'Soil Temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'soilmmt'
    case('phiw')    
       if (itype==0) ncinfo = 'Soil Moisture'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'soilmmt'
    case('a_Qnet')    
       if (itype==0) ncinfo = 'Surface Net Radiation'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'mmt'
    case('a_G0')    
       if (itype==0) ncinfo = 'Surface Ground Heat Flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'mmt'
    case('Tsoil1')    
       if (itype==0) ncinfo = 'Bulk soil temperature (11)'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('Tsoil2')    
       if (itype==0) ncinfo = 'Bulk soil temperature (2)'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('Tsoil3')    
       if (itype==0) ncinfo = 'Bulk soil temperature (3)'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('Tsoil4')    
       if (itype==0) ncinfo = 'Bulk soil temperature (4)'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'

    case default
       if (myid==0) print *, 'ABORTING: variable not found in ncinfo, ',trim(short_name)
       call appl_abort(0)
    end select

  end function ncinfo

end module ncio
