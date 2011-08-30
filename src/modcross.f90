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
use defs, only : long

implicit none

  logical            :: lcross = .false., ldocross, lxy = .false., lxz = .false., lyz = .false.
  real               :: dtcross = 60, xcross = 0., ycross = 0., zcross = 0.
  integer            :: icross,jcross,kcross
  integer(kind=long) :: idtcross, itimenextcross

  integer :: ncross = 0
  character(len=80), allocatable, dimension(:) :: crossname

  integer :: nccrossid, nccrossrec

contains

!-------------------
!CROSSSECTIONS
!-------------------
  subroutine initcross(rtimee, expname)
    use mpi_interface,   only : wrxid, wryid
    use modnetcdf,       only : open_nc
    use grid,            only : nzp, nxp, nyp, zt, xt, yt, xm, ym
    real, intent(in) :: rtimee
    character(len=*), intent(in) :: expname


    integer :: i, j, k
    character(len=4) :: cmpicoordx, cmpicoordy
    write(cmpicoordx,'(i4.4)') wrxid
    write(cmpicoordy,'(i4.4)') wryid
    if (lxy) then
      do k=2,nzp-1
        if (zt(k)>zcross) exit
      end do
      kcross = k
    end if
    if (lxz) then
      if (ycross < ym(3) .or. ycross >= ym(nyp - 1)) then
        lxz = .false.
      else
        do j=3,nyp-2
          if (yt(j)>ycross) exit
        end do
        jcross = j
      end if
    end if
    if (lyz) then
      if (xcross < xm(3) .or. xcross >= xm(nxp - 1)) then
        lyz = .false.
      else
        do i=3,nxp-2
          if (xt(i)>xcross) exit
        end do
        icross = i
      end if
    end if

    if (.not.(lxy .or. lxz .or. lyz)) lcross = .false.
    if (lcross) call open_nc(trim(expname)//'.out.cross.'//cmpicoordx//'.'//cmpicoordy//'.nc', nccrossid, nccrossrec, rtimee)

  end subroutine initcross

  integer function addcross(name, longname, unit, loc)
    use modnetcdf,   only : addvar_nc
    use grid,        only : ictr, ihlf, nzp, nxp, nyp, zt, xt, yt, zm, xm, ym, &
                            zname, zlongname, zunit, zhname, zhlongname, &
                            xname, xlongname, xunit, xhname, xhlongname, &
                            yname, ylongname, yunit, yhname, yhlongname, &
                            tname, tlongname, tunit

    character (*), intent(in)     :: name, longname, unit
    integer, intent(in)           :: loc(3)
    character (40), dimension(3) :: dimname, dimlongname, dimunit
    character (len=80), dimension(size(crossname)) :: ctmp
    real, allocatable, dimension(:,:) :: dimvalues
    integer, dimension(3)             :: dimsize

    if (lcross) then
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

      if (lxz) then
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
      if (lyz) then
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
          dimvalues(1:nyp-4,2)  = ym(3:nxp-2)
          dimname(2)            = yhname
          dimlongname(2)        = yhlongname
        end if
        call addvar_nc(nccrossid, trim(name)//'yz', 'yz crosssection of '//trim(longname), &
        unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
        crossname(ncross) = name
      end if
      if (lxz) then
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
      
      if (lxy) then
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
          dimvalues(1:nyp-4,2)  = ym(3:nxp-2)
          dimname(2)            = yhname
          dimlongname(2)        = yhlongname
        end if
        call addvar_nc(nccrossid, trim(name)//'yz', 'yz crosssection of '//trim(longname), &
        unit, dimname, dimlongname, dimunit, dimsize, dimvalues)
        crossname(ncross) = name
      end if
    end if
    addcross = ncross
  end function addcross

  subroutine triggercross(rtimee)
    use grid,      only : tname
    use modnetcdf, only : writevar_nc
    real, intent(in) :: rtimee
    
    call writevar_nc(nccrossid, tname, rtimee, nccrossrec)

  end subroutine triggercross

  subroutine writecross(nr, am)
    use grid,     only : nxp, nyp, nzp
    use modnetcdf,       only : writevar_nc

    integer, intent(in)                :: nr
    real, dimension(1:nzp, 1:nxp, 1:nyp), intent(in) :: am
    real, dimension(:,:), allocatable  :: cross

  ! XZ crosssection
    if (lxz) then
      allocate(cross(2:nzp-1, 3:nxp-2))
      cross(2:nzp,2:nxp) = am(2:nzp-1, 3:nxp-2, jcross)
      call writevar_nc(nccrossid, trim(crossname(nr))//'xz', cross, nccrossrec)
      deallocate(cross)
    end if

  ! YZ crosssection
    if (lyz) then
      allocate(cross(2:nzp-1, 3:nyp-2))
      cross(2:nzp,2:nxp) = am(2:nzp-1, icross, 3:nyp-2)
      call writevar_nc(nccrossid, trim(crossname(nr))//'yz', cross, nccrossrec)
      deallocate(cross)
    end if

  ! XY crosssection
    if (lxy) then
      allocate(cross(3:nxp-2, 3:nyp-2))
      cross(2:nzp,2:nxp) = am(kcross, 3:nzp-2, 3:nyp-2)
      call writevar_nc(nccrossid, trim(crossname(nr))//'xy', cross, nccrossrec)
      deallocate(cross)
    end if

  end subroutine writecross

end module modcross

