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
!> All kinds of auxilarily functions and variables for
module modcross

implicit none

  logical            :: lcross = .false., ldocross, lxy = .false., lxz = .false., lyz = .false.
  real               :: xcross = 0., ycross = 0., zcross(10) = 0.
  integer            :: icross,jcross,kcross(10) = 0, nkcross
  real               :: threstracer = 2
  real               :: thres_rw    = 1.e-8
  integer :: ncross = 0
  character(len=7), allocatable, dimension(:) :: crossname
  character(len=7),  dimension(10) :: hname
  character(len=80), dimension(10) :: hlname
  character(len=7) :: hname_prc
  integer, parameter :: nvar_all = 49
  character (len=7), dimension(nvar_all)  :: crossvars =  (/ &
         'u      ','v      ','w      ','t      ','r      ', & !1-5
         'l      ','rp     ','np     ','tv     ','ricep  ', & !6-10
         'nicep  ','rsnowp ','rgrpp  ','nsnowp ','ngrpp  ', & !11-15
         'rhailp ','nhailp ','lwp    ','rwp    ','iwp    ', & !16-20
         'swp    ','gwp    ','hwp    ','prc_acc','cnd_acc', & !21-25
         'cev_acc','rev_acc','cldbase','cldtop ','rwpbase', & !26-30
         'rwptop ','tracer ','trcpath','trcbase','trctop ', & !31-35
         'wdev_cl','wdev_sc','w_cld  ','tdev_cl','tdev_sc', & !36-40
         't_cld  ','qdev_cl','qdev_sc','q_cld  ','tv_cl  ', & !41-45
         'tv_sc  ','tv_cld ','core   ','th_e   '/)            !46-49
  integer :: nccrossxzid,nccrossyzid,nccrossxyid, nccrossrec, nvar

  interface writecross
    module procedure writecross_2D
    module procedure writecross_3D
  end interface writecross


contains

!-------------------
!CROSSSECTIONS
!-------------------
  subroutine initcross(rtimee, expname)
    use mpi_interface,   only : wrxid, wryid
    use modnetcdf,       only : open_nc, addvar_nc
    use grid,            only : nzp, nxp, nyp, zt, xt, yt, xm, ym, &
                                tname, tlongname, tunit

    real, intent(in) :: rtimee
    integer :: n
    character(len=*), intent(in) :: expname
    integer :: i, j, k
    character(len=4) :: cmpicoordx, cmpicoordy

!     return
    write(cmpicoordx,'(i4.4)') wrxid
    write(cmpicoordy,'(i4.4)') wryid
    if (lxy) then
      nkcross = 10
      docross: do n = 1, 10
        select case (nint(zcross(n)))
        case(0)
          nkcross = n - 1
          exit docross
        case(-1)
          hname(n)  = 'zi'
          hlname(n) = 'Maximum buoyancy variance'
        case(-2)
          hname(n)  = 'cb'
          hlname(n) = 'Cloud base'
        case(-3)
          hname(n)  = 'lcl'
          hlname(n) = 'at lifting condensation level'
        case default
          write(hname(n),'(i4.4)') nint(zcross(n))
          hlname(n) = ' at '//trim(hname(n))//' m'
          if (nint(zcross(n)) == 0) then
          elseif (zcross(n) < 0) then
          else
            do k=2,nzp-1
              if (zt(k)>zcross(n)) exit
            end do
            kcross(n) = k
          end if
        end select
      end do docross
    end if
    if (lxz) then
      if (ycross < ym(2) .or. ycross >= ym(nyp - 2)) then
        lxz = .false.
      else
        do j=3,nyp-2
          if (yt(j)>ycross) exit
        end do
        jcross = j
        call open_nc(trim(expname)//'.out.xz.'//cmpicoordx//'.'//cmpicoordy//'.nc', nccrossxzid, nccrossrec, rtimee)
      end if

    end if
    if (lyz) then
      if (xcross < xm(2) .or. xcross >= xm(nxp - 2)) then
        lyz = .false.
      else
        do i=3,nxp-2
          if (xt(i)>xcross) exit
        end do
        icross = i
        call open_nc(trim(expname)//'.out.yz.'//cmpicoordx//'.'//cmpicoordy//'.nc', nccrossyzid, nccrossrec, rtimee)
      end if
    end if

!     if (.not.(lxy .or. lxz .or. lyz)) lcross = .false.
    if (lcross) call open_nc(trim(expname)//'.out.xy.'//cmpicoordx//'.'//cmpicoordy//'.nc', nccrossxyid, nccrossrec, rtimee)
    do n=1,nkcross
     if (zcross(n) < 0) then
       call addvar_nc(nccrossxyid, trim(hname(n)),trim(hlname(n)), 'm', (/tname/), (/tlongname/), (/tunit/))
      end if
    end do
    do n=1,nvar_all
      call addcross(crossvars(n))
    end do
  end subroutine initcross

  subroutine addcross(name)
    use modnetcdf,     only : addvar_nc
    use mpi_interface, only : appl_abort
    use grid,          only : level, ictr, ihlf, nzp, nxp, nyp, zt, xt, yt, zm, xm, ym, &
                            zname, zlongname, zunit, &
                            xname, xlongname, xunit, &
                            yname, ylongname, yunit, &
                            tname, tlongname, tunit, &
                            lwaterbudget, lcouvreux, prc_lev

    character (*), intent(in)     :: name
    character (40), dimension(3) :: dimname, dimlongname, dimunit
    character (len=80), dimension(size(crossname)) :: ctmp
    character (len=80) :: longname, unit
    real, allocatable, dimension(:,:) :: dimvalues
    integer, dimension(3)             :: loc, dimsize
    integer :: n, k
    logical :: iscross

    if (lcross) then
      loc = (/ictr, ictr, ictr/)
      unit = 'kg/kg'
      iscross = .false.
      select case (trim(name))
      case('u')
        loc = (/ictr, ihlf, ictr/)
        longname =  'Zonal wind'
        unit = 'm/s'
        iscross = .true.
     case('v')
        loc = (/ictr, ictr, ihlf/)
        longname =  'Meridional wind'
        unit = 'm/s'
        iscross = .true.
      case('w')
        loc = (/ihlf, ictr, ictr/)
        longname =  'Meridional wind'
        unit = 'm/s'
        iscross = .true.
      case('t')
        longname =  'Liquid water potential temperature'
        unit = 'K'
        iscross = .true.
      case('tv')
        longname =  'Virtual potential temperature'
        unit = 'K'
        iscross = .true.
      case('r')
        if (level < 1) return
        longname =  'Total water content'
        iscross = .true.
      case('l')
        if (level < 2) return
        longname =  'Liquid water content'
        iscross = .true.
      case('cvrx')
        if (.not. lcouvreux) return
        longname =  'Couvreux Scalar'
        iscross = .true.
      case('rp')
        if (level < 3) return
        longname =  'Rain water content'
        iscross = .true.
      case('np')
        if (level < 3) return
        longname =  'Rain water number density'
        iscross = .true.
      case('ricep')
        if (level < 4) return
        longname =  'Cloud ice content'
        iscross = .true.
      case('nicep')
        if (level < 4) return
        longname =  'Cloud ice number density'
        iscross = .true.
      case('rsnowp')
        if (level < 4) return
        longname =  'Snow content'
        iscross = .true.
      case('nsnowp')
        if (level < 5) return
        longname =  'Snow number density'
        iscross = .true.
      case('rgrpp')
        if (level < 4) return
        longname =  'Graupel content'
        iscross = .true.
      case('ngrpp')
        if (level < 5) return
        longname =  'Graupel number density'
        iscross = .true.
      case('rhailp')
        if (level < 5) return
        longname =  'Hail content'
        iscross = .true.
      case('tracer')
        if (.not. lcouvreux) return
        longname =  'Tracer'
        iscross = .true.
        unit = '-'
      case('nhailp')
        if (level < 5) return
        longname =  'Hail number density'
        iscross = .true.
      case('lwp')
        if (level < 2) return
        longname = 'Liquid water path'
        unit = 'kg/m2'
      case('rwp')
        if (level < 3) return
        longname = 'Rain water path'
        unit = 'kg/m2'
      case('iwp')
        if (level < 4) return
        longname = 'Ice water path'
        unit = 'kg/m2'
      case('swp')
        if (level < 4) return
        longname = 'Snow water path'
        unit = 'kg/m2'
      case('gwp')
        if (level < 4) return
        longname = 'Graupel water path'
        unit = 'kg/m2'
      case('hwp')
        if (level < 5) return
        longname = 'Hail water path'
        unit = 'kg/m2'
      case ('prc_acc')
        if (level < 3) return
        longname = 'acc. precip'
        unit = 'kg/m2'
      case ('cnd_acc')
        if (.not.lwaterbudget) return
        longname = 'acc. condensation'
        unit = 'kg/m2'
      case ('cev_acc')
        if (.not.lwaterbudget) return
        longname = 'acc. evaporation of cloud water'
        unit = 'kg/m2'
      case ('rev_acc')
        if (level < 3) return
        longname = 'acc. evaporation of rain water'
        unit = 'kg/m2'
      case ('cldbase')
        if (level < 2) return
        longname = 'Cloud base height'
        unit = 'm'
      case ('w_cld')
        if (level < 2) return
        longname = 'Av. In cloud vertical velocity'
        unit = 'm/s'
      case ('wdev_cl')
        longname = 'Av. Cloud layer vertical velocity'
        unit = 'm/s'
      case ('wdev_sc')
        if (level < 2) return
        longname = 'SubCloud layer vertical velocity'
        unit = 'm/s'
      case ('t_cld')
        if (level < 2) return
        longname = 'Av. In cloud temperature deviation'
        unit = 'K'
      case ('tdev_cl')
        if (level < 2) return
        longname = 'Av. Cloud layer temperature deviation'
        unit = 'K'
      case ('tdev_sc')
        if (level < 2) return
        longname = 'SubCloud layer temperature deviation'
        unit = 'K'
      case ('tv_cld')
        if (level < 2) return
        longname = 'Av. In cloud virt. pot. temperature deviation'
        unit = 'K'
      case ('tvdev_cl')
        if (level < 2) return
        longname = 'Av. Cloud layer virt. temperature deviation'
        unit = 'K'
      case ('tvdev_sc')
        if (level < 2) return
        longname = 'SubCloud layer virt. temperature deviation'
        unit = 'K'
      case ('q_cld')
        if (level < 2) return
        longname = 'Av. In cloud humidity deviation'
        unit = 'kg/kg'
      case ('qdev_cl')
        if (level < 2) return
        longname = 'Av. Cloud layer humidity deviation'
        unit = 'kg/kg'
      case ('qdev_sc')
        if (level < 2) return
        longname = 'SubCloud layer humidity deviation'
        unit = 'kg/kg'
      case ('cldtop')
        if (level < 2) return
        longname = 'Cloud top height'
        unit = 'm'
      case ('rwpbase')
        if (level < 3) return
        longname = 'Rain base height'
        unit = 'm'
      case ('rwptop')
        if (level < 3) return
        longname = 'Rain top height'
        unit = 'm'
      case ('trcpath')
        if (.not.lcouvreux) return
        longname = 'Tracer path'
        unit = 'm^-2'
      case ('trcbase')
        if (.not.lcouvreux) return
        longname = 'Tracer base height'
        unit = 'm'
      case ('trctop')
        if (.not.lcouvreux) return
        longname = 'Tracer top height'
        unit = 'm'
      case ('core')
        if (level < 2) return
        longname = 'Max. in cloud buoyancy'
        unit = 'K'
      case ('th_e')
        if (level < 2) return
        longname = 'equivalent potential temperature'
        unit = 'K'
        iscross = .true.
      case default
        return
      end select
      if(ncross /= 0) then
        ctmp = crossname
        deallocate(crossname)
      end if
      ncross   = ncross + 1


      allocate(crossname(ncross))
      if (ncross > 1) crossname(1:ncross-1) = ctmp
      crossname(ncross)     = name

      allocate(dimvalues(max(nzp-2, nxp-4, nyp-4),3))
      dimname(3)            = tname
      dimlongname(3)        = tlongname
      dimunit(3)            = tunit
      dimsize    = 0
      dimvalues  = 0
      if (iscross) then
        if (lxz) then
          dimunit(1) = zunit
          dimunit(2) = xunit
          dimsize(1) = nzp - 2
          dimsize(2) = nxp - 4
          dimvalues(1:nzp-2,1)  = zt(2:nzp-1)
          dimname(1)            = zname
          dimlongname(1)        = zlongname

          dimvalues(1:nxp-4,2)  = xt(3:nxp-2)
          dimname(2)            = xname
          dimlongname(2)        = xlongname
          call addvar_nc(nccrossxzid, trim(name), 'xz crosssection of '//trim(longname), &
          unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
          crossname(ncross) = name
        end if
        if (lyz) then
          dimunit(1) = zunit
          dimunit(2) = yunit
          dimsize(1) = nzp - 2
          dimsize(2) = nyp - 4
          dimvalues(1:nzp-2,1)  = zt(2:nzp-1)
          dimname(1)            = zname
          dimlongname(1)        = zlongname
          dimvalues(1:nyp-4,2)  = yt(3:nyp-2)
          dimname(2)            = yname
          dimlongname(2)        = ylongname
          call addvar_nc(nccrossyzid, trim(name), 'yz crosssection of '//trim(longname), &
          unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
          crossname(ncross) = name
        end if

        if (lxy) then
          dimunit(1) = xunit
          dimunit(2) = yunit
          dimsize(1) = nxp - 4
          dimsize(2) = nyp - 4
          dimvalues(1:nxp-4,1)  = xt(3:nxp-2)
          dimname(1)            = xname
          dimlongname(1)        = xlongname
          dimvalues(1:nyp-4,2)  = yt(3:nyp-2)
          dimname(2)            = yname
          dimlongname(2)        = ylongname
          do n = 1, nkcross
            call addvar_nc(nccrossxyid, trim(name)//trim(hname(n)), 'xy crosssection of ' &
                           //trim(longname)//' at '//trim(hlname(n)), &
            unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
          end do
          crossname(ncross) = name
        end if
      else
        dimunit(1) = xunit
        dimunit(2) = yunit
        dimsize(1) = nxp - 4
        dimsize(2) = nyp - 4
        dimvalues(1:nxp-4,1)  = xt(3:nxp-2)
        dimname(1)            = xname
        dimlongname(1)        = xlongname
        dimvalues(1:nyp-4,2)  = yt(3:nyp-2)
        dimname(2)            = yname
        dimlongname(2)        = ylongname
        if (trim(name) == 'prc_acc') then
          do n = 1, count(prc_lev>=0)
            do k = 2,nzp
              if (zt(k) >= prc_lev(n)) then
                prc_lev(n) = k
                exit
              end if
            end do
            write(hname_prc,'(i4.4)') nint(zt(prc_lev(n)))
            call addvar_nc(nccrossxyid, trim(name)//'_'//trim(hname_prc), trim(longname), &
            unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
          end do
        else
          call addvar_nc(nccrossxyid, trim(name), trim(longname), &
          unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
        end if
        crossname(ncross) = name
      end if
    end if

  end subroutine addcross

  subroutine triggercross(rtimee)
    use grid,      only : level,nxp, nyp, nzp, tname, zt, zm, dzi_m, dzi_t, a_up, a_vp, a_wp, a_tp, a_rp, liquid, a_rpp, a_npp, &
       a_ricep, a_nicep, a_rsnowp, a_nsnowp, a_rgrp, a_ngrp, a_rhailp, a_nhailp, &
       prc_acc, cnd_acc, cev_acc, rev_acc, a_cvrxp, lcouvreux, a_theta, pi0, pi1, a_pexnr, prc_lev, umean, vmean, th00
    use modnetcdf, only : writevar_nc, fillvalue_double
    use util,      only : get_avg3, get_var3, calclevel
    use defs,      only : ep2,cp,cpr, p00
    use thrm,      only: rslf
    real, intent(in) :: rtimee
    real             :: tstar, exner, tk
    real, dimension(3:nxp-2,3:nyp-2) :: tmp
    real, dimension(nzp,nxp,nyp) :: tracer, tv, interp, th_e,p
    real, dimension(nzp)         :: c1, thvar, tvbar, tvenv, tvcld
    integer :: n, nn, i, j, k, ct, cb, zi, lcl

    if (.not. lcross) return

    call writevar_nc(nccrossxyid, tname, rtimee, nccrossrec)
    if (lxz) then
      nccrossrec = nccrossrec - 1
      call writevar_nc(nccrossxzid, tname, rtimee, nccrossrec)
    end if
    if (lyz) then
      nccrossrec = nccrossrec - 1
      call writevar_nc(nccrossyzid, tname, rtimee, nccrossrec)
    end if

    if (level == 0) then
      do j=3,nyp-2
        do i=3,nxp-2
            do k=1,nzp
              tv(k,i,j) = a_theta(k,i,j)
            end do
        end do
      end do
    else if (level == 1) then
      do j=3,nyp-2
        do i=3,nxp-2
            do k=1,nzp
              tv(k,i,j) = a_theta(k,i,j)*(1.+ep2*a_rp(k,i,j))
            end do
        end do
      end do
    else
      do j=3,nyp-2
        do i=3,nxp-2
            do k=1,nzp
              tv(k,i,j) = a_theta(k,i,j)*(1.+ep2*a_rp(k,i,j) - liquid(k,i,j))
            end do
        end do
      end do
     do j=3,nyp-2
       do i=3,nxp-2
          do k=1,nzp
             exner = (pi0(k)+pi1(k)+a_pexnr(k,i,j))/cp
             p(k,i,j) = p00 * (exner)**cpr
             tk=a_theta(k,i,j)*exner
             tstar = 1./(1./(tk-55.)-log(a_rp(k,i,j)/rslf(p(k,i,j),tk))/2840.)+55.
             th_e(k,i,j) = tk*(p00/p(k,i,j))**(0.2854*(1.-0.28*a_rp(k,i,j)))* &
                           exp(a_rp(k,i,j)*(1.+0.81*a_rp(k,i,j))*(3376./tstar-2.54))
          end do
       end do
      end do

    end if


    call get_avg3(nzp,nxp,nyp,tv,tvbar)
    do j=3,nyp-2
       do i=3,nxp-2
          do k=1,nzp
             tv(k,i,j) = tv(k,i,j) - tvbar(k)
          end do
       end do
    end do


    call get_avg3(nzp,nxp,nyp, tv, c1)
    call get_var3(nzp,nxp,nyp, tv, c1, thvar)
    zi = maxloc(thvar,1)
    if (level > 1) then
      call calclevel(liquid, cb, 'base')
      call calclevel(liquid, ct, 'top')

      if (cb >= nzp-1) then
        cb = zi
      end if

  !-------- calc lcl -------------------
      lcl=cb+nint(200.*dzi_m(cb))
  !-------------------------------------
      if (lcl >= nzp-1) then
        lcl = 0
      end if
    end if

    do n = 1, nkcross
      select case (nint(zcross(n)))
      case(-1)
        kcross(n) = zi
        call writevar_nc(nccrossxyid, trim(hname(n)), zt(zi), nccrossrec)

      case(-2)
        call writevar_nc(nccrossxyid, trim(hname(n)), zt(cb), nccrossrec)
        kcross(n) = cb
      case(-3)
        call writevar_nc(nccrossxyid, trim(hname(n)), zt(lcl), nccrossrec)
        kcross(n) = lcl
      case default
      end select
    end do
    if (lcouvreux) then
      call scalexcess(a_cvrxp, tracer)
    end if

    do n = 1, ncross
      select case(trim(crossname(n)))
      case('u')
        do j=3,nyp-2
          do i=3,nxp-2
              do k=1,nzp
                interp(k,i,j) = 0.5*(a_up(k,i-1,j) + a_up(k,i,j)) + umean
              end do
          end do
        end do
        call writecross(crossname(n), interp)
      case('v')
        do j=3,nyp-2
          do i=3,nxp-2
              do k=1,nzp
                interp(k,i,j) = 0.5*(a_vp(k,i,j-1) + a_vp(k,i,j)) + vmean
              end do
          end do
        end do
        call writecross(crossname(n), interp)
      case('w')
        do j=3,nyp-2
          do i=3,nxp-2
              do k=2,nzp
                interp(k,i,j) = 0.5*dzi_t(k) * (a_wp(k-1,i,j) / dzi_m(k) + a_wp(k,i,j) / dzi_m(k-1))
              end do
          end do
        end do
        call writecross(crossname(n), interp)
      case('t')
        call writecross(crossname(n), a_tp+th00)
      case('r')
        call writecross(crossname(n), a_rp)
      case('l')
        call writecross(crossname(n), liquid)
      case('tv')
        call writecross(crossname(n), tv)
      case('rp')
        call writecross(crossname(n), a_rpp)
      case('np')
        call writecross(crossname(n), a_npp)
      case('ricep')
        call writecross(crossname(n), a_ricep)
      case('nicep')
        call writecross(crossname(n), a_nicep)
      case('rsnowp')
        call writecross(crossname(n), a_rsnowp)
      case('nsnowp')
        call writecross(crossname(n), a_nsnowp)
      case('rgrpp')
        call writecross(crossname(n), a_rgrp)
      case('ngrpp')
        call writecross(crossname(n), a_ngrp)
      case('rhailp')
        call writecross(crossname(n), a_rhailp)
      case('nhailp')
        call writecross(crossname(n), a_nhailp)
      case ('th_e')
        call writecross(crossname(n), th_e)
      case('tracer')
        call writecross(crossname(n), tracer)
      case('lwp')
        call calcintpath(liquid, tmp)
        call writecross(crossname(n), tmp)
      case('rwp')
        call calcintpath(a_rpp, tmp, thres_rw)
        call writecross(crossname(n), tmp)
      case('iwp')
        call calcintpath(a_ricep, tmp)
        call writecross(crossname(n), tmp)
      case('swp')
        call calcintpath(a_rsnowp, tmp)
        call writecross(crossname(n), tmp)
      case('gwp')
        call calcintpath(a_rgrp, tmp)
        call writecross(crossname(n), tmp)
      case('hwp')
        call calcintpath(a_rhailp, tmp)
        call writecross(crossname(n), tmp)
      case('prc_acc')
        do nn = 1,count(prc_lev>0)
          write(hname_prc,'(i4.4)') nint(zt(prc_lev(nn)))
          tmp = prc_acc(3:nxp-2, 3:nyp-2,nn)
          call writecross(trim(crossname(n))//'_'//trim(hname_prc), tmp)
          prc_acc(:,:,nn) = 0.
        end do
      case('cnd_acc')
        tmp = cnd_acc(3:nxp-2, 3:nyp-2)
        call writecross(crossname(n), tmp)
        cnd_acc = 0.
      case('cev_acc')
        tmp = cev_acc(3:nxp-2, 3:nyp-2)
        call writecross(crossname(n), tmp)
        cev_acc = 0.
      case('rev_acc')
        tmp = rev_acc(3:nxp-2, 3:nyp-2)
        call writecross(crossname(n), tmp)
        rev_acc = 0.
      case ('cldbase')
        call calcbase(liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('cldtop')
        call calctop(liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('clddept')
        call calcdepth(liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('rwpbase')
        call calcbase(a_rpp, tmp, thres_rw)
        call writecross(crossname(n), tmp)
      case ('rwptop')
        call calctop(a_rpp, tmp, thres_rw)
        call writecross(crossname(n), tmp)
      case ('rwpdept')
        call calcdepth(a_rpp, tmp, thres_rw)
        call writecross(crossname(n), tmp)
      case('trcpath')
        call calcintpath(tracer, tmp)
        call writecross(crossname(n), tmp)
      case ('trcbase')
        call calcbase(tracer, tmp)
        call writecross(crossname(n), tmp)
      case ('trctop')
        call calctop(tracer, tmp)
        call writecross(crossname(n), tmp)
      case ('trcdept')
        call calcdepth(tracer, tmp)
        call writecross(crossname(n), tmp)
      case ('w_cld')
        call calcdevmask(a_wp, liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('wdev_cl')
        call calcdev(a_wp, cb, ct, tmp)
        call writecross(crossname(n), tmp)
      case ('wdev_sc')
        call calcdev(a_wp, 2, cb-1, tmp)
        call writecross(crossname(n), tmp)
      case ('t_cld')
        call calcdevmask(a_tp, liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('tdev_cl')
        call calcdev(a_tp, cb, ct, tmp)
        call writecross(crossname(n), tmp)
      case ('tdev_sc')
        call calcdev(a_tp, 2, cb-1, tmp)
        call writecross(crossname(n), tmp)
      case ('tv_cld')
        call calcdevmask(tv, liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('tv_cl')
        call calcdev(tv, cb, ct, tmp)
        call writecross(crossname(n), tmp)
      case ('tv_sc')
        call calcdev(a_tp, 2, cb-1, tmp)
        call writecross(crossname(n), tmp)
      case ('q_cld')
        call calcdevmask(a_rp, liquid, tmp)
        call writecross(crossname(n), tmp)
      case ('qdev_cl')
        call calcdev(a_rp, cb, ct, tmp)
        call writecross(crossname(n), tmp)
      case ('qdev_sc')
        call calcdev(a_rp, 2, cb-1, tmp)
        call writecross(crossname(n), tmp)
      case ('core')
        call calcmax(tv, liquid, tmp)
        call writecross(crossname(n), tmp)
      end select
    end do
  end subroutine triggercross

  subroutine writecross_3D(crossname, am)
    use grid,     only : nxp, nyp, nzp
    use modnetcdf,       only : writevar_nc
    use mpi_interface,   only : wrxid, wryid
    character(*), intent(in)                :: crossname
    real, dimension(1:nzp, 1:nxp, 1:nyp), intent(in) :: am
    real, dimension(:,:), allocatable  :: cross
    integer :: n
  ! XZ crosssection

    if (lxz) then
      allocate(cross(2:nzp-1, 3:nxp-2))
      cross = am(2:nzp-1, 3:nxp-2, jcross)
      call writevar_nc(nccrossxzid, trim(crossname), cross, nccrossrec)
      deallocate(cross)
    end if

  ! YZ crosssection
    if (lyz) then
      allocate(cross(2:nzp-1, 3:nyp-2))
      cross = am(2:nzp-1, icross, 3:nyp-2)
      call writevar_nc(nccrossyzid, trim(crossname), cross, nccrossrec)
      deallocate(cross)
    end if

  ! XY crosssection
    if (lxy) then
      allocate(cross(3:nxp-2, 3:nyp-2))
      do n=1,nkcross
        cross = am(kcross(n), 3:nxp-2, 3:nyp-2)
        call writevar_nc(nccrossxyid, trim(crossname)//trim(hname(n)), cross, nccrossrec)
      end do
      deallocate(cross)
    end if

  end subroutine writecross_3D


  subroutine writecross_2D(crossname, am)
    use grid,     only : nxp, nyp, nzp
    use modnetcdf,       only : writevar_nc

    character(*), intent(in)                :: crossname
    real, dimension(:,:), intent(in) :: am
    call writevar_nc(nccrossxyid, trim(crossname), am, nccrossrec)

  end subroutine writecross_2D

  subroutine exitcross
    use modnetcdf, only : close_nc
    if (lcross) call close_nc(nccrossxyid)
    if (lxz)    call close_nc(nccrossxzid)
    if (lyz)    call close_nc(nccrossyzid)
  end subroutine exitcross

  subroutine calcintpath(varin, varout, threshold)
    use modnetcdf, only : fillvalue_double
    use grid, only : nzp, nxp, nyp, dn0, zm
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(out), dimension(3:,3:)  :: varout
    real, intent(in), optional :: threshold
    integer :: i, j, k, km1
    real :: thres
    varout = 0.
    if (present(threshold)) then
      thres = threshold
    else
      thres = 0.0
    end if
    do j=3,nyp-2
      do i=3,nxp-2
        do k=2,nzp-1
          km1=max(1,k-1)
          if (varin(k,i,j) > thres) then
            varout(i,j) = varout(i,j)+varin(k,i,j)*(zm(k)-zm(km1))*dn0(k)
          end if
        enddo
        if (varout(i,j) > 1e10) varout(i,j) = fillvalue_double
      end do
    end do
  end subroutine calcintpath

  subroutine calcmax(varin, mask, varout)
    use modnetcdf, only : fillvalue_double
    use grid, only : nzp, nxp, nyp
    real, intent(in), dimension(:,:,:) :: varin, mask
    real, intent(out), dimension(3:,3:)  :: varout
    integer :: i, j, k, n

    varout = fillvalue_double
    do j=3,nyp-2
      do i=3,nxp-2
        do k=2,nzp-1
          if (mask(k,i,j) > 0.) then
            varout(i,j) = max(varout(i,j),varin(k,i,j))
          end if
        enddo
      end do
    end do
  end subroutine calcmax

  subroutine calcavg(varin, mask, varout)
    use modnetcdf, only : fillvalue_double
    use grid, only : nzp, nxp, nyp
    use mpi_interface, only : real_array_par_sum
    real, intent(in), dimension(:,:,:) :: varin, mask
    real, intent(out), dimension(:)  :: varout
    real, dimension(SIZE(varout))  :: varoutsum, lvaroutsum, gvaroutsum
    integer :: i, j, k, n
    real, dimension(SIZE(varout)) :: ncl, lncl, gncl

    varout = fillvalue_double
    do k=1,nzp
      ncl(k)=0.
      varoutsum(k)=0.
        do j=3,nyp-2
        do i=3,nxp-2
          if (mask(k,i,j) == 0.) then
            varoutsum(k) = varoutsum(k)+varin(k,i,j)
            ncl(k) = ncl(k)+1
          else
            varoutsum(k) = varoutsum(k)
            ncl(k)= ncl(k)
          end if
        end do
        end do
      end do

    n=size(varout)
    lvaroutsum=varoutsum
    call real_array_par_sum(lvaroutsum,gvaroutsum,n)
    lncl=ncl
    call real_array_par_sum(lncl,gncl,n)

    do k=1,nzp
       if (gncl(k) == 0. ) then
         varout(k) = 0.
       else
         varout(k) = gvaroutsum(k)/gncl(k)
       end if
    end do
  end subroutine calcavg

  subroutine calcavgcld(varin, mask, varout)
    use modnetcdf, only : fillvalue_double
    use grid, only : nzp, nxp, nyp
    use mpi_interface, only : real_array_par_sum
    real, intent(in), dimension(:,:,:) :: varin, mask
    real, intent(out), dimension(:)  :: varout
    real, dimension(SIZE(varout))  :: varoutsum, lvaroutsum, gvaroutsum
    integer :: i, j, k, n
    real, dimension(SIZE(varout)) :: ncl, lncl, gncl

    do k=1,nzp
      ncl(k)=0.
      varoutsum(k)=0.
        do j=3,nyp-2
        do i=3,nxp-2
          if (mask(k,i,j) > 0. ) then
            varoutsum(k) = varoutsum(k)+varin(k,i,j)
            ncl(k)=ncl(k)+1
          else
            varoutsum(k) = varoutsum(k)
            ncl(k)= ncl(k)
          end if
        end do
        end do
      end do

    n=size(varout)
    lvaroutsum=varoutsum
    call real_array_par_sum(lvaroutsum,gvaroutsum,n)
    lncl=ncl
    call real_array_par_sum(lncl,gncl,n)

    do k=1,nzp
       if (gncl(k) == 0. ) then
         varout(k) = 0
       else
         varout(k) = gvaroutsum(k)/gncl(k)
       end if
    end do
  end subroutine calcavgcld


  subroutine calcbase(varin, varout, threshold)
    use grid, only : nzp, nxp, nyp, zt
    use modnetcdf, only : fillvalue_double
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(out), dimension(3:,3:)  :: varout
    real, intent(in), optional :: threshold
    integer :: i, j, k, km1
    real :: thres
    varout = fillvalue_double
    if (present(threshold)) then
      thres = threshold
    else
      thres = 0.0
    end if
    do j = 3, nyp - 2
      do i = 3, nxp - 2
        base:do k = 2, nzp - 1
          if (varin(k,i,j) > thres) then
            varout(i,j) = zt(k)
            exit base
          end if
        end do base
      end do
    end do
  end subroutine calcbase

  subroutine calctop(varin, varout, threshold)
    use grid, only : nzp, nxp, nyp, zt
    use modnetcdf, only : fillvalue_double
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(out), dimension(3:,3:)  :: varout
    real, intent(in), optional :: threshold
    integer :: i, j, k, km1
    real :: thres
    varout = fillvalue_double
    if (present(threshold)) then
      thres = threshold
    else
      thres = 0.0
    end if
    do j = 3, nyp - 2
      do i = 3, nxp - 2
        top:do k = nzp - 1, 2, -1
          if (varin(k,i,j) > thres) then
            varout(i,j) = zt(k)
            exit top
          end if
        end do top
      end do
    end do
  end subroutine calctop

  subroutine calcdepth(varin, varout, threshold)
    use grid, only : nzp, nxp, nyp, zm
    use modnetcdf, only : fillvalue_double
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(out), dimension(3:,3:)  :: varout
    real, intent(in), optional :: threshold
    integer :: i, j, k, km1
    real :: thres
    varout = 0.
    if (present(threshold)) then
      thres = threshold
    else
      thres = 0.0
    end if
    do j = 3, nyp - 2
      do i = 3, nxp - 2
        do k = 2, nzp - 1
          km1=max(1,k-1)
          if (varin(k,i,j) > thres) then
            varout(i,j) = varout(i,j) + zm(k)-zm(k-1)
          end if
        end do
        if (varout(i,j) == 0.) varout(i,j) = fillvalue_double
      end do
    end do
  end subroutine calcdepth

  subroutine calcdev(varin, base, top, varout)
    use grid, only : nzp, nxp, nyp, zm, zt, a_wp, dzi_t
    use util, only : get_avg3
    use modnetcdf, only : fillvalue_double
    real, intent(in), dimension(:,:,:) :: varin
    integer, intent(in) :: base, top
    real, intent(out), dimension(3:,3:)  :: varout
    real, dimension(nzp) :: mean, div, divmin
    integer :: i, j, k, km1
    varout = 0.
    divmin = 0.
    if (base >= nzp-2 .or. base <= 1 .or. top >= nzp-2 .or. top <= 1 ) then
      varout = fillvalue_double
      return
    end if
    call get_avg3(nzp, nyp, nxp,varin,mean)
      do j=3,nyp-2
        do i=3,nxp-2
          do k=base,top
            varout(i,j) = varout(i,j) + (varin(k,i,j) - mean(k))/(dzi_t(k) * (zm(top) - zm(base-1)))
          end do
        end do
      end do
  end subroutine calcdev

  subroutine calcdevmask(varin, mask, varout)
    use grid, only : nzp, nxp, nyp, zm, zt, a_wp, dzi_t
    use util, only : get_avg3
    use modnetcdf, only : fillvalue_double
    real, intent(in), dimension(:,:,:) :: varin, mask
    real, intent(out), dimension(3:,3:)  :: varout
    real, dimension(nzp) :: mean, div, divmin
    integer :: i, j, k, nr
    varout = 0.
    divmin = 0.
    call get_avg3(nzp, nyp, nxp,varin,mean)
      do j=3,nyp-2
        do i=3,nxp-2
          nr = count(mask(2:nzp-1,i,j)>0)
          if (nr == 0) then
            varout(i,j) = fillvalue_double
          else
            do k=2,nzp-1
              if (mask(k,i,j) > 0) then
                varout(i,j) = varout(i,j) + varin(k,i,j) - mean(k) 
              end if
            end do
            varout(i,j) = varout(i,j)/real(nr)
          end if
        end do
      end do
  end subroutine calcdevmask
!The output is the number of std. deviations over 1 that the local value of the local value
!is larger than the slab average. Only for points with an upward positive velocity.
  subroutine scalexcess(varin, varout)
    use grid, only : nzp, nxp, nyp, zm, zt, a_wp
    use util, only : get_avg3, get_var3
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(out), dimension(:,:,:)  :: varout
    real, dimension(nzp) :: mean, div, divmin
    integer :: i, j, k, km1
    varout = 0.
    divmin = 0.
    call get_avg3(nzp, nyp, nxp,varin,mean)
    call get_var3(nzp, nyp, nxp,varin,mean, div)
    div = sqrt(div)

    do k = 2, nzp -1
      divmin(k:nzp-1) = divmin(k:nzp-1) + div(k) * 0.05 * (zm(k)-zm(k-1))/zt(k)
      div(k) = 1./max(1e-10, max(divmin(k), div(k)))
    end do
    do j=3,nyp-2
      do i=3,nxp-2
        do k=2,nzp-1
!           if (0.5 * (a_wp(k,i,j) + a_wp(k-1,i,j))> 0.) then
            varout(k,i,j) = (varin(k,i,j) - mean(k)) * div(k) - threstracer
            varout(k,i,j) = max(0.,varout(k,i,j))
!           end if
        enddo
      end do
    end do
  end subroutine scalexcess


end module modcross

