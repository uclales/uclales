!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module cldwtr

  use defs, only : nv, mb
  implicit none
  integer, save :: nsizes
  logical, save :: Initialized = .False.
  logical, save :: iceInitialized = .False.
  logical, save :: grpInitialized = .False.
  integer, save :: mbs,mbir

  real, allocatable    :: re(:), fl(:), bz(:,:), wz(:,:), gz(:,:)
  real, allocatable    :: ap(:,:), bp(:,:), cps(:,:,:), dps(:,:), cpir(:,:)
  real, allocatable    :: bg(:), wgf(:), gg(:)
  real :: gwc

contains
  !
  !---------------------------------------------------------------------------
  ! Subroutine cloud_init initialize data arrays for the cloud model, 
  ! checking for consistency between band structure of cloud model and CKD
  !
  subroutine init_cldwtr

    use ckd, only : band, center
    integer, parameter  :: nrec = 21600

    real, dimension(mb) :: cntrs

    integer             :: ib, i, nbands
    character (len=12)  :: frmt

    open ( unit = 71, file = 'datafiles/cldwtr.dat', status = 'old', recl=nrec)
    read (71,'(2I3)') nsizes, nbands
    if (nbands /= mb .or. nsizes*nbands*15 > nrec) &
         stop 'TERMINATING: incompatible cldwtr.dat file'

    allocate (re(nsizes),fl(nsizes),bz(nsizes,mb),wz(nsizes,mb),gz(nsizes,mb))
    write(frmt,'(A1,I2.2,A8)') '(',mb,'E15.7)    '
    read (71,frmt) (cntrs(i), i=1,mb)
    do i=1,mb
       if (spacing(1.) < abs(cntrs(i)- center(band(i))) ) &
            stop 'TERMINATING: cloud properties not matched to band structure'
    end do

    write(frmt,'(A1,I2.2,A9)') '(',nsizes,'E15.7)   '
    read (71,frmt) (re(i), i=1,nsizes)
    read (71,frmt) (fl(i), i=1,nsizes)

    write(frmt,'(A1,I4.4,A7)') '(',nsizes*mb,'E15.7) '
    read (71,frmt) ((bz(i,ib), i=1,nsizes), ib=1,mb)
    read (71,frmt) ((wz(i,ib), i=1,nsizes), ib=1,mb)
    read (71,frmt) ((gz(i,ib), i=1,nsizes), ib=1,mb)
    close (71)

    if (minval((/bz,wz,gz/)) < 0.) &
         stop 'TERMINATING: cloud properties out of bounds'

    Initialized = .True.

  end subroutine init_cldwtr

  !
  !---------------------------------------------------------------------------
  ! Subroutine cloud_init initialize data arrays for the cloud model, 
  ! checking for consistency between band structure of cloud model and CKD
  !
  subroutine init_cldice

    use ckd, only : center
    integer, parameter  :: nrec = 21600


    integer             :: i, j

    mbs=6
    mbir=12

    allocate (ap(3,mb),bp(4,mb),cps(4,4,mbs),dps(4,mbs),cpir(4,mbir))
    open ( unit = 71, file = 'datafiles/cldice.dat', status = 'old', recl=nrec)
    do i=1,mb
       read (71,'(3E10.3)') ap(1,i), ap(2,i), ap(3,i)
    enddo
    do i=1,mb
       read (71,'(4E12.5)') bp(1,i), bp(2,i), bp(3,i), bp(4,i)
    enddo
    do i=1,mbs
       do j=1,4
          read (71,'(4E12.5)') cps(1,j,i), cps(2,j,i), cps(3,j,i), cps(4,j,i)
       enddo
    enddo
    do i=1,mbs
       read (71,'(4E12.5)') dps(1,i), dps(2,i), dps(3,i), dps(4,i)
    enddo
    do i=1,mbir
       read (71,'(4E13.4)') cpir(1,i), cpir(2,i), cpir(3,i), cpir(4,i)
    enddo

    close (71)

    iceInitialized = .True.

  end subroutine init_cldice

  !
  !---------------------------------------------------------------------------
  ! Subroutine cloud_init initialize data arrays for the cloud model, 
  ! checking for consistency between band structure of cloud model and CKD
  !
  subroutine init_cldgrp

    integer, parameter  :: nrec = 21600
    integer             :: i
    character (len=12)  :: frmt
    
    gwc  = 1.5e10

    allocate (bg(mb),wgf(mb),gg(mb))
    open ( unit = 71, file = 'datafiles/cldgrp.dat', status = 'old', recl=nrec)
    write(frmt,'(A1,I2.2,A8)') '(',mb,'E10.3)    '
    read (71,frmt) (bg(i), i=1,mb)
    write(frmt,'(A1,I2.2,A8)') '(',mb,'F7.4)    '
    read (71,frmt) (wgf(i), i=1,mb)
    read (71,frmt) (gg(i), i=1,mb)

    close (71)

    grpInitialized = .True.

  end subroutine init_cldgrp

  ! -----------------------------------------------------------------------
  ! Subroutine cloud_water:  calculates the optical depth (tw), single 
  ! scattering albedo (ww), and phase function (www(4)) given the cloud 
  ! water [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  subroutine cloud_water ( ib, pre, pcw, dz, tw, ww, www )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pre, pcw, dz
    real, intent (out) :: tw(nv), ww(nv), www(nv,4)

    integer :: k, j, j0, j1
    real    :: gg, wght, cwmks

    if (.not.Initialized) stop 'TERMINATING: Cloud not Initialized'

    do k = 1, nv
       cwmks = pcw(k)*1.e-3!convert to m to km
       if ( cwmks .ge. 1.e-8) then
          j = 0
          do while (j<nsizes .and. pre(k) > re(j+1))
             j = j + 1
          end do
          if (j >= 1 .and. j < nsizes) then
             j1 = j+1
             wght = (pre(k)-re(j))/(re(j1)-re(j))
             tw(k) = dz(k) * cwmks * ( bz(j,ib) / fl(j) +   &
                  ( bz(j1,ib) / fl(j1) - bz(j,ib) / fl(j) ) /    &
                  ( 1.0 / re(j1) - 1.0 / re(j) ) * ( 1.0 / pre(k) &
                  - 1.0 / re(j) ) )
             ww(k) = wz(j,ib) + (wz(j1,ib) - wz(j,ib) ) * wght
             gg    = gz(j,ib) + (gz(j1,ib) - gz(j,ib) ) * wght
          else
             j0 = max(j,1)
             tw(k) = dz(k) * cwmks * (bz(j0,ib)/fl(j0))
             ww(k) = wz(j0,ib) 
             gg    = gz(j0,ib)          
          end if
          www(k,1) = 3.0 * gg
          do j=2,4
             wght = real(2*j+1)/real(2*j-1)
             www(k,j) = www(k,j-1) * gg * wght
          end do
       else
          www(k,:) = 0.0
          tw(k) = 0.0
          ww(k) = 0.0
          gg    = 0.
       end if
    end do

    return
  end subroutine cloud_water

  ! -----------------------------------------------------------------------
  ! Subroutine cloud_ice:  calculates the optical depth (ti), single 
  ! scattering albedo (wi), and phase function (wwi(4)) given the cloud 
  ! ice [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  subroutine cloud_ice ( ib, pde, pci, dz, ti, wi, wwi )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pde, pci, dz
    real, intent (out) :: ti(nv), wi(nv), wwi(nv,4)

    integer :: k
    real    :: gg, wght, cwmks
    real    :: fw1, fw2, fw3, wf1, wf2, wf3, wf4, x1, x2, x3, x4, ibr, fd

    if (.not.iceInitialized) stop 'TERMINATING: Ice not Initialized'

    do k = 1, nv
       cwmks = pci(k)!here we don't need the factor 1000
       if ( (cwmks .ge. 1.e-5).and.(pde(k).gt.0.)) then
	     fw1 = pde(k)
	     fw2 = fw1 * pde(k)
	     fw3 = fw2 * pde(k)
             ti(k) = dz(k) * cwmks * ( ap(1,ib) + &
      	     ap(2,ib) / fw1 + ap(3,ib) / fw2 )
             wi(k) = 1.0 - ( bp(1,ib) + bp(2,ib) * fw1 + &
      	     bp(3,ib) * fw2 + bp(4,ib) * fw3 )
             if (wi(k).lt.0.) print*,'bad wi, ',wi(k),ib,k,bp(1,ib),bp(2,ib),bp(3,ib),bp(4,ib),fw1,fw2,fw3
             if (ti(k).lt.0.) print*,'bad ti, ',ti(k),ib,k,cwmks,dz(k),ap(1,ib),ap(2,ib),ap(3,ib),fw1,fw2
	     if ( ib .le. mbs ) then ! shortwave
	       fd = dps(1,ib) + dps(2,ib) * fw1 + &
               dps(3,ib) * fw2 + dps(4,ib) * fw3
               wf1 = cps(1,1,ib) + cps(2,1,ib) * fw1 + &
               cps(3,1,ib) * fw2 + cps(4,1,ib) * fw3
               wwi(k,1) = ( 1.0 - fd ) * wf1 + 3.0 * fd
	       wf2 = cps(1,2,ib) + cps(2,2,ib) * fw1 + &
               cps(3,2,ib) * fw2 + cps(4,2,ib) * fw3
               wwi(k,2) = ( 1.0 - fd ) * wf2 + 5.0 * fd
       	       wf3 = cps(1,3,ib) + cps(2,3,ib) * fw1 + &
               cps(3,3,ib) * fw2 + cps(4,3,ib) * fw3
               wwi(k,3) = ( 1.0 - fd ) * wf3 + 7.0 * fd
               wf4 = cps(1,4,ib) + cps(2,4,ib) * fw1 + &
               cps(3,4,ib) * fw2 + cps(4,4,ib) * fw3
               wwi(k,4) = ( 1.0 - fd ) * wf4 + 9.0 * fd
             else ! longwave
               ibr = ib - mbs
               gg = cpir(1,ibr) + cpir(2,ibr) * fw1 + &
               cpir(3,ibr) * fw2 + cpir(4,ibr) * fw3
	       x1 = gg
               x2 = x1 * gg
               x3 = x2 * gg
               x4 = x3 * gg
               wwi(k,1) = 3.0 * x1
	       wwi(k,2) = 5.0 * x2
               wwi(k,3) = 7.0 * x3
               wwi(k,4) = 9.0 * x4
	     endif
       else
          wwi(k,:) = 0.0
          ti(k) = 0.0
          wi(k) = 0.0
          gg    = 0.
       end if
    end do

    return
  end subroutine cloud_ice


  ! -----------------------------------------------------------------------
  ! Subroutine cloud_grp:  calculates the optical depth (ti), single 
  ! scattering albedo (wi), and phase function (wwi(4)) given the cloud 
  ! ice [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  ! tgr, wgr, and wwgr are the optical depth, single scattering albedo,
  ! and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
  ! due to the Mie scattering of graupel for a given layer. 
  !                        Jan. 19, 1993
  ! -----------------------------------------------------------------------

  subroutine cloud_grp ( ib, pcg, dz, tgr, wgr, wwgr )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pcg, dz
    real, intent (out) :: tgr(nv), wgr(nv), wwgr(nv,4)

    integer :: i
    real    :: cwmks
    real    :: y1, y2, y3, y4, x1, x2, x3, x4, ibr, fd
    if (.not.grpInitialized) stop 'TERMINATING: Ice not Initialized'

    x1 = gg(ib)
    x2 = x1 * gg(ib)
    x3 = x2 * gg(ib)
    x4 = x3 * gg(ib)
    y1 = 3.0 * x1
    y2 = 5.0 * x2
    y3 = 7.0 * x3
    y4 = 9.0 * x4
    do i = 1, nv
       cwmks = pcg(i)*1.e-3! convert to km
       if ( cwmks .lt. 1.0e-8 ) then
          tgr(i) = 0.0
          wgr(i) = 0.0
          wwgr(i,1) = 0.0
          wwgr(i,2) = 0.0
          wwgr(i,3) = 0.0
          wwgr(i,4) = 0.0
       else
          tgr(i) = dz(i) * cwmks * bg(ib) / gwc
          wgr(i) = wgf(ib)
          wwgr(i,1) = y1
          wwgr(i,2) = y2
          wwgr(i,3) = y3
          wwgr(i,4) = y4
       endif       
    end do

  end subroutine cloud_grp

  ! ---------------------------------------------------------------------------
  ! linear interpolation between two points, returns indicies of the 
  ! interpolation points and weights
  !
  subroutine interpolate(x,ny,y,i1,i2,alpha)

    integer, intent (in) :: ny
    real, intent (in)    :: x, y(ny)

    integer, intent (out) :: i1, i2
    real, intent (out)    :: alpha

    if (y(1) < y(2)) stop 'TERMINATING: band centers increasing'

    i2 = 1
    do while (x < y(i2) .and. i2 < ny)
       i2 = i2+1
    end do
    i1 = max(1,i2-1)
    alpha = 1.

    if(i2.ne.i1) alpha = (x-y(i1))/(y(i2)-y(i1))
    if (alpha <0 .or. alpha >1) print 600, x, y(1), y(ny), alpha

    return

600 format(/'CLOUD_INIT WARNING:  Extrapolating because data out of range', &
         /1x,'x = ',F8.1,', ymax = ',F7.0,', ymin =',F7.0,', alpha = ',F6.3)
  end subroutine interpolate

end module cldwtr
