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
  real               :: dtcross = 60, xcross = 0., ycross = 0., zcross = 0.
  integer            :: icross,jcross,kcross

  integer :: ncross = 0
  character(len=7), allocatable, dimension(:) :: crossname
  integer, parameter :: nvar_all = 26
  character (len=7), dimension(nvar_all)  :: crossvars =  (/ &
         'u      ','v      ','w      ','t      ','r      ', & !1-5
         'l      ','rp     ','np     ','ricep  ','nicep  ', & !6-10
         'rsnowp ','rgrpp  ','nsnowp ','ngrpp  ','rhailp ', & !11-15
         'nhailp ','lwp    ','rwp    ','iwp    ','swp    ', & !16-20
         'gwp    ','hwp    ','prc_acc','cnd_acc','cev_acc', & !21-25
         'rev_acc'/)                                          !26
  integer :: nccrossid, nccrossrec, nvar
  
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
    use modnetcdf,       only : open_nc
    use grid,            only : nzp, nxp, nyp, zt, xt, yt, xm, ym
    real, intent(in) :: rtimee
    integer :: n
    character(len=*), intent(in) :: expname
    integer :: i, j, k
    character(len=4) :: cmpicoordx, cmpicoordy
    
!     return
    write(cmpicoordx,'(i4.4)') wrxid
    write(cmpicoordy,'(i4.4)') wryid
    if (lxy) then
      do k=2,nzp-1
        if (zt(k)>zcross) exit
      end do
      kcross = k
    end if
    if (lxz) then
      if (ycross < ym(2) .or. ycross >= ym(nyp - 2)) then
        lxz = .false.
      else
        do j=3,nyp-2
          if (yt(j)>ycross) exit
        end do
        jcross = j
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
      end if
    end if
!     if (.not.(lxy .or. lxz .or. lyz)) lcross = .false.
    if (lcross) call open_nc(trim(expname)//'.out.cross.'//cmpicoordx//'.'//cmpicoordy//'.nc', nccrossid, nccrossrec, rtimee)
  
    do n=1,nvar_all
      call addcross(crossvars(n))
    end do
  end subroutine initcross

  subroutine addcross(name)
    use modnetcdf,     only : addvar_nc
    use mpi_interface, only : appl_abort
    use grid,          only : level, ictr, ihlf, nzp, nxp, nyp, zt, xt, yt, zm, xm, ym, &
                            zname, zlongname, zunit, zhname, zhlongname, &
                            xname, xlongname, xunit, xhname, xhlongname, &
                            yname, ylongname, yunit, yhname, yhlongname, &
                            tname, tlongname, tunit, &
                            lwaterbudget

    character (*), intent(in)     :: name
    character (40), dimension(3) :: dimname, dimlongname, dimunit
    character (len=80), dimension(size(crossname)) :: ctmp
    character (len=80) :: longname, unit
    real, allocatable, dimension(:,:) :: dimvalues
    integer, dimension(3)             :: loc, dimsize
    integer :: n
    logical :: ldoxy, ldoxz, ldoyz
    if (lcross) then
      loc = (/ictr, ictr, ictr/)
      unit = 'kg/kg'
      do n = 1, nvar_all
        select case (trim(name))
        case('u')
          loc = (/ictr, ihlf, ictr/)
          longname =  'Zonal wind'
          unit = 'm/s'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('v')
          loc = (/ictr, ictr, ihlf/)
          longname =  'Meridional wind'
          unit = 'm/s'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('w')
          loc = (/ihlf, ictr, ictr/)
          longname =  'Meridional wind'
          unit = 'm/s'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('t')
          longname =  'Liquid water potential temperature'
          unit = 'K'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('r')
          if (level < 1) return
          longname =  'Total water content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('l')
          if (level < 2) return
          longname =  'Liquid water content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('rp')
          if (level < 3) return
          longname =  'Rain water content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('np')
          if (level < 3) return
          longname =  'Rain water number density'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('ricep')
          if (level < 4) return
          longname =  'Cloud ice content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('nicep')
          if (level < 4) return
          longname =  'Cloud ice number density'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('rsnowp')
          if (level < 4) return
          longname =  'Snow content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('nsnowp')
          if (level < 5) return
          longname =  'Snow number density'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('rgrpp')
          if (level < 4) return
          longname =  'Graupel content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('ngrpp')
          if (level < 5) return
          longname =  'Graupel number density'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('rhailp')
          if (level < 5) return
          longname =  'Hail content'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('nhailp')
          if (level < 5) return
          longname =  'Hail number density'
          ldoxz = lxz
          ldoyz = lyz
          ldoxy = lxy
        case('lwp')
          if (level < 2) return
          longname = 'Liquid water path'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case('rwp')
          if (level < 3) return
          longname = 'Rain water path'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case('iwp')
          if (level < 4) return
          longname = 'Ice water path'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case('swp')
          if (level < 4) return
          longname = 'Snow water path'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case('gwp')
          if (level < 4) return
          longname = 'Graupel water path'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case('hwp')
          if (level < 5) return
          longname = 'Hail water path'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case ('prc_acc')
          if (.not.lwaterbudget) return
          longname = 'acc. precip'
          unit = 'kg/m2'
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case ('cnd_acc')
          if (.not.lwaterbudget) return
          longname = 'acc. condensation'
          unit = 'kg/m2'          
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case ('cev_acc')
          if (.not.lwaterbudget) return
          longname = 'acc. evaporation of cloud water'
          unit = 'kg/m2'          
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case ('rev_acc')
          if (.not.lwaterbudget) return
          longname = 'acc. evaporation of rain water'
          unit = 'kg/m2'          
          ldoxz = .false.
          ldoyz = .false.
          ldoxy = lxy
          if (.not. ldoxy) return
        case default
          return
        end select
      end do
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

      if (ldoxz) then
        dimunit(1) = zunit
        dimunit(2) = xunit
        dimsize(1) = nzp - 2
        dimsize(2) = nxp - 4
        if(loc(1)==ictr) then
          dimvalues(1:nzp-2,1)  = zt(2:nzp-1)
          dimname(1)            = zname
          dimlongname(1)        = zlongname
        else
          dimvalues(1:nzp-2,1)  = zm(2:nzp-1)
          dimname(1)            = zhname
          dimlongname(1)        = zhlongname
        end if
        if(loc(2)==ictr) then
          dimvalues(1:nxp-4,2)  = xt(3:nxp-2)
          dimname(2)            = xname
          dimlongname(2)        = xlongname
        else
          dimvalues(1:nxp-4,2)  = xm(3:nxp-2)
          dimname(2)            = xhname
          dimlongname(2)        = xhlongname
        end if
        call addvar_nc(nccrossid, trim(name)//'xz', 'xz crosssection of '//trim(longname), &
        unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
        crossname(ncross) = name
      end if
      if (ldoyz) then
        dimunit(1) = zunit
        dimunit(2) = yunit
        dimsize(1) = nzp - 2
        dimsize(2) = nyp - 4
        if(loc(1)==ictr) then
          dimvalues(1:nzp-2,1)  = zt(2:nzp-1)
          dimname(1)            = zname
          dimlongname(1)        = zlongname
        else
          dimvalues(1:nzp-2,1)  = zm(2:nzp-1)
          dimname(1)            = zhname
          dimlongname(1)        = zhlongname
        end if
        if(loc(2)==ictr) then
          dimvalues(1:nyp-4,2)  = yt(3:nyp-2)
          dimname(2)            = yname
          dimlongname(2)        = ylongname
        else
          dimvalues(1:nyp-4,2)  = ym(3:nyp-2)
          dimname(2)            = yhname
          dimlongname(2)        = yhlongname
        end if
        call addvar_nc(nccrossid, trim(name)//'yz', 'yz crosssection of '//trim(longname), &
        unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
        crossname(ncross) = name
      end if
      
      if (ldoxy) then
        dimunit(1) = xunit
        dimunit(2) = yunit
        dimsize(1) = nxp - 4
        dimsize(2) = nyp - 4
        if(loc(1)==ictr) then
          dimvalues(1:nxp-4,1)  = xt(3:nxp-2)
          dimname(1)            = xname
          dimlongname(1)        = xlongname
        else
          dimvalues(1:nxp-4,1)  = xm(3:nxp-2)
          dimname(1)            = xhname
          dimlongname(1)        = xhlongname
        end if
        if(loc(2)==ictr) then
          dimvalues(1:nyp-4,2)  = yt(3:nyp-2)
          dimname(2)            = yname
          dimlongname(2)        = ylongname
        else
          dimvalues(1:nyp-4,2)  = ym(3:nyp-2)
          dimname(2)            = yhname
          dimlongname(2)        = yhlongname
        end if
        call addvar_nc(nccrossid, trim(name)//'xy', 'xy crosssection of '//trim(longname), &
        unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
        crossname(ncross) = name
      end if
    end if
    
  end subroutine addcross

  subroutine triggercross(rtimee)
    use grid,      only : nxp, nyp, tname, a_up, a_vp, a_wp, a_tp, a_rp, liquid, a_rpp, a_npp, &
       a_ricep, a_nicep, a_rsnowp, a_nsnowp, a_rgrp, a_ngrp, a_rhailp, a_nhailp, &
       prc_acc, cnd_acc, cev_acc, rev_acc
    use modnetcdf, only : writevar_nc
    real, intent(in) :: rtimee
    real, dimension(3:nxp-2,3:nyp-2) :: tmp
    integer :: n
    if (.not. lcross) return
    call writevar_nc(nccrossid, tname, rtimee, nccrossrec)
    do n = 1, ncross
      select case(trim(crossname(n)))
      case('u')
        call writecross(crossname(n), a_up)
      case('v')
        call writecross(crossname(n), a_vp)
      case('w')
        call writecross(crossname(n), a_wp)
      case('t')
        call writecross(crossname(n), a_tp)
      case('r')
        call writecross(crossname(n), a_rp)
      case('l')
        call writecross(crossname(n), liquid)
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
      case('lwp')
        call calcintpath(liquid, tmp)
        call writecross(crossname(n), tmp)
      case('rwp')
        call calcintpath(a_rpp, tmp)
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
        tmp = prc_acc(3:nxp-2, 3:nyp-2)
        call writecross(crossname(n), tmp)
        prc_acc = 0.
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
      end select
    end do

  end subroutine triggercross

  subroutine writecross_3D(crossname, am)
    use grid,     only : nxp, nyp, nzp
    use modnetcdf,       only : writevar_nc

    character(*), intent(in)                :: crossname
    real, dimension(1:nzp, 1:nxp, 1:nyp), intent(in) :: am
    real, dimension(:,:), allocatable  :: cross

  ! XZ crosssection
    if (lxz) then
      allocate(cross(2:nzp-1, 3:nxp-2))
      cross = am(2:nzp-1, 3:nxp-2, jcross)
      call writevar_nc(nccrossid, trim(crossname)//'xz', cross, nccrossrec)
      deallocate(cross)
    end if

  ! YZ crosssection
    if (lyz) then
      allocate(cross(2:nzp-1, 3:nyp-2))
      cross = am(2:nzp-1, icross, 3:nyp-2)
      call writevar_nc(nccrossid, trim(crossname)//'yz', cross, nccrossrec)
      deallocate(cross)
    end if

  ! XY crosssection
    if (lxy) then
      allocate(cross(3:nxp-2, 3:nyp-2))
      cross = am(kcross, 3:nxp-2, 3:nyp-2)
      call writevar_nc(nccrossid, trim(crossname)//'xy', cross, nccrossrec)
      deallocate(cross)
    end if

  end subroutine writecross_3D
  

  subroutine writecross_2D(crossname, am)
    use grid,     only : nxp, nyp, nzp
    use modnetcdf,       only : writevar_nc

    character(*), intent(in)                :: crossname
    real, dimension(:,:), intent(in) :: am
    print *, crossname
    call writevar_nc(nccrossid, trim(crossname)//'xy', am, nccrossrec)

  end subroutine writecross_2D
  
  subroutine exitcross
    use modnetcdf, only : close_nc
    call close_nc(nccrossid)
  end subroutine exitcross

  subroutine calcintpath(varin, varout)
    use grid, only : nzp, nxp, nyp, dn0, zm
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(out), dimension(3:,3:)  :: varout
    integer :: i, j, k, km1
    varout = 0.
    do j=3,nyp-2
      do i=3,nxp-2
        do k=2,nzp-1
          km1=max(1,k-1)
          varout(i,j) = varout(i,j)+varin(k,i,j)*(zm(k)-zm(km1))*dn0(k)
        enddo
      end do
    end do
  end subroutine calcintpath
end module modcross

