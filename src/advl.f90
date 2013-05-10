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
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module advl

  implicit none

  integer :: advm = 4   !< advection scheme; 2=2nd(O) centered, all else defaults to 4th(O) centered

contains
  !
  ! ----------------------------------------------------------------------
  ! LADVECT:  This does the leap-frog advection using a fourth
  ! order centered in space algorithm.  Boundary point fluxes are not 
  ! computed... hence the routine is hardwired to a cyclic domain.  
  ! Vertical advection is density weighted consistent with the anelastic
  ! approximation
!
  subroutine ladvect

    use grid, only : a_ut, a_vt, a_wt, a_scr1, a_scr2, a_up,a_vp,a_wp,      &
         nxp, nyp, nzp, dzi_t, dzi_m, dxi, dyi, dn0
    use stat, only : sflg, updtst, acc_tend
    use util, only : get_avg3

    real, allocatable ::  dzmri(:), dztri(:)
    real, allocatable ::  v1(:), v2(:), v3(:), v4(:)
    integer :: ret

    allocate(dzmri(nzp), dztri(nzp), v1(nzp), v2(nzp), v3(nzp), v4(nzp), stat=ret)

    if(ret.ne.0) then
        print *, 'ladvect: auxiliary-vectors allocation failed'
        return
    endif

    if (sflg) call acc_tend(nzp,nxp,nyp,a_up,a_vp,a_wp,a_ut,a_vt,a_wt,        &
         v1,v2,v3,1,'adv')

    !
    ! prepare density weights for use vertical advection
    !
    call advl_prep(nzp,nxp,nyp,a_wp,a_scr1,dn0,dzi_t,dzi_m,dztri,dzmri)
    
    !
    ! advection of u by (u,v,w) all at current timelevel.  also when flag
    ! is set updated statistical array with uw flux derived from ladvzu
    !
    if(advm==2) then
      call ladvxu2nd(nzp,nxp,nyp,a_up,a_ut,a_scr2,dxi)
      call ladvyu2nd(nzp,nxp,nyp,a_up,a_ut,a_vp,a_scr2,dyi)
      call ladvzu2nd(nzp,nxp,nyp,a_up,a_ut,a_scr1,a_scr2,dztri)
    else
      call ladvxu(nzp,nxp,nyp,a_up,a_ut,a_scr2,dxi)
      call ladvyu(nzp,nxp,nyp,a_up,a_ut,a_vp,a_scr2,dyi)
      call ladvzu(nzp,nxp,nyp,a_up,a_ut,a_scr1,a_scr2,dztri)
    end if 

    if (sflg) then
       call get_avg3(nzp,nxp,nyp,a_scr2,v4)
       call updtst(nzp,'adv',-1,v4,1)
    end if
    
    !
    ! advection of v by (u,v,w) all at current timelevel.  also when flag
    ! is set updated statistical array with uw flux derived from ladvzu
    !    
    if(advm==2) then
      call ladvxv2nd(nzp,nxp,nyp,a_up,a_vp,a_vt,a_scr2,dxi)
      call ladvyv2nd(nzp,nxp,nyp,a_vp,a_vt,a_scr2,dyi)
      call ladvzv2nd(nzp,nxp,nyp,a_vp,a_vt,a_scr1,a_scr2,dztri)  
    else
      call ladvxv(nzp,nxp,nyp,a_up,a_vp,a_vt,a_scr2,dxi)
      call ladvyv(nzp,nxp,nyp,a_vp,a_vt,a_scr2,dyi)
      call ladvzv(nzp,nxp,nyp,a_vp,a_vt,a_scr1,a_scr2,dztri)  
    end if   
 
    if (sflg) then
       call get_avg3(nzp,nxp,nyp,a_scr2,v4)
       call updtst(nzp,'adv',-2,v4,1)
    end if
    
    !
    ! advection of w by (u,v,w) all at current timelevel.  also when flag
    ! is set updated statistical array with uw flux derived from ladvzu
    !
    if(advm==2) then
      call ladvxw2nd(nzp,nxp,nyp,a_up,a_wp,a_wt,a_scr2,dxi)
      call ladvyw2nd(nzp,nxp,nyp,a_vp,a_wp,a_wt,a_scr2,dyi)
      call ladvzw2nd(nzp,nxp,nyp,a_wp,a_wt,a_scr1,a_scr2,dzmri)
    else   
      call ladvxw(nzp,nxp,nyp,a_up,a_wp,a_wt,a_scr2,dxi)
      call ladvyw(nzp,nxp,nyp,a_vp,a_wp,a_wt,a_scr2,dyi)
      call ladvzw(nzp,nxp,nyp,a_wp,a_wt,a_scr1,a_scr2,dzmri)
    end if
 
    if (sflg) then
       call get_avg3(nzp,nxp,nyp,a_scr2,v4)
       call updtst(nzp,'adv',-3,v4,1)
    end if

    if (sflg) call acc_tend(nzp,nxp,nyp,a_up,a_vp,a_wp,a_ut,a_vt,a_wt,        &
         v1,v2,v3,2,'adv')
         
    deallocate(dzmri, dztri, v1, v2, v3, v4)

  end subroutine ladvect

  ! ###########################
  ! 4th order centered
  ! ###########################
  !
  ! ----------------------------------------------------------------------
  ! Ladvxu: Advection of U, by U, div is a scratch array, 
  ! tendencies are accumulated if ut variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxu(n1,n2,n3,u,ut,flx,dxi)

    integer, intent (in) ::  n1,n2,n3
    real, intent(in)     :: u(n1,n2,n3),dxi
    real, intent(inout)  :: ut(n1,n2,n3)
    real, intent(out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=3,n2-1
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(u(k,i,j)+u(k,i-1,j)) - 0.083333          &
                  *(u(k,i+1,j)+u(k,i-2,j)))*.5*(u(k,i,j)+u(k,i-1,j))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i+1,j)-flx(k,i,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxu
  !
  ! ----------------------------------------------------------------------
  ! LADVYU: Advection of U, by V, div is a scratch array, 
  ! tendencies are accumulated if fu variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyu(n1,n2,n3,u,ut,v,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),v(n1,n2,n3),dyi
    real, intent (inout)  :: ut(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=2,n3-2
       do i=1,n2-1
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(u(k,i,j)+u(k,i,j+1)) - 0.083333          &
                  *(u(k,i,j-1)+u(k,i,j+2)))*.5*(v(k,i,j)+v(k,i+1,j))
          end do
       end do
    end do

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i,j)-flx(k,i,j-1))*dyi
          end do
       end do
    end do

  end subroutine ladvyu
  !
  ! ----------------------------------------------------------------------
  ! Subroutine ladvzu: Advection of U, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if ut variable, v1 is the inverse density
  ! array valid at thermo points
  !
  subroutine ladvzu(n1,n2,n3,u,ut,wm,flx,v1)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),wm(n1,n2,n3),v1(n1)
    real, intent (inout)  :: ut(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=1,n2-1
          flx(1,i,j)=0.
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(u(k,i,j)+u(k+1,i,j)) - 0.083333          &
                  *(u(k-1,i,j)+u(k+2,i,j)))*.5*(wm(k,i,j)+wm(k,i+1,j))
          end do
          flx(n1-1,i,j)=0.
          flx(n1,i,j)=0.
       end do

       do i=1,n2-1
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*v1(k)
          end do
       end do
    end do

  end subroutine ladvzu
  !
  ! ----------------------------------------------------------------------
  ! LADVXV: Advection of V, by U, div is a scratch array, 
  ! tendencies are accumulated if vt variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxv(n1,n2,n3,u,v,vt,flx,dxi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),v(n1,n2,n3),dxi
    real, intent (inout)  :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3-1
       do i=2,n2-2
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(v(k,i,j)+v(k,i+1,j)) - 0.083333          &
                  *(v(k,i-1,j)+v(k,i+2,j)))*.5*(u(k,i,j)+u(k,i,j+1))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j)-flx(k,i-1,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxv
  !
  ! ----------------------------------------------------------------------
  ! LADVYV: Advection of V, by V, div is a scratch array, 
  ! tendencies are accumulated if vt variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyv(n1,n2,n3,v,vt,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: v(n1,n2,n3),dyi
    real, intent (inout)  :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=3,n3-1
       do i=1,n2
          do k=2,n1-1
             flx(k,i,j)=(0.583333*(v(k,i,j)+v(k,i,j-1))-0.083333           &
                  *(v(k,i,j+1)+v(k,i,j-2)))*.5*(v(k,i,j)+v(k,i,j-1))
          end do
       end do
    end do

    do j=3,n3-2
       do i=1,n2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j+1)-flx(k,i,j))*dyi
          end do
       end do
    end do

  end subroutine ladvyv
  !
  ! ----------------------------------------------------------------------
  ! LADVZV: Advection of V, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if vt variable, v1 is the inverse density
  ! array valid at thermo points
  !
  subroutine ladvzv(n1,n2,n3,v,vt,wm,flx,v1)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: v(n1,n2,n3),wm(n1,n2,n3),v1(n1)
    real, intent (inout)  :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3-1
       do i=1,n2
          flx(1,i,j)=0.
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(v(k,i,j)+v(k+1,i,j))-0.083333            &
                  *(v(k-1,i,j)+v(k+2,i,j)))*.5*(wm(k,i,j)+wm(k,i,j+1))
          end do
          flx(n1-1,i,j)=0.
          flx(n1,i,j)=0.
       end do

       do i=1,n2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*v1(k)
          end do
       end do
    end do

  end subroutine ladvzv
  !
  ! ----------------------------------------------------------------------
  ! LADVXW: Advection of W, by U, div is a scratch array, 
  ! tendencies are accumulated if wt variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxw(n1,n2,n3,u,w,wt,flx,dxi)

    integer, intent (in)  :: n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),w(n1,n2,n3),dxi
    real, intent (inout)  :: wt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=2,n2-2
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(w(k,i,j)+w(k,i+1,j))-0.083333*          &
                  (w(k,i-1,j)+w(k,i+2,j)))*.5*(u(k,i,j)+u(k+1,i,j))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k,i,j)-flx(k,i-1,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxw
  !
  ! ----------------------------------------------------------------------
  ! LADVYW: Advection of W, by V, div is a scratch array, 
  ! tendencies are accumulated if wt variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyw(n1,n2,n3,vm,w,wt,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: vm(n1,n2,n3),w(n1,n2,n3),dyi
    real, intent (inout)  :: wt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=2,n3-2
       do i=1,n2
          do k=2,n1-2
             flx(k,i,j)=(0.583333*(w(k,i,j)+w(k,i,j+1))-0.083333            &
                  *(w(k,i,j-1)+w(k,i,j+2)))*.5*(vm(k,i,j)+vm(k+1,i,j))
          end do
       end do
    end do

    do j=3,n3-2
       do i=1,n2
          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k,i,j)-flx(k,i,j-1))*dyi
          end do
       end do
    end do

  end subroutine ladvyw
  !
  ! ----------------------------------------------------------------------
  ! LADVZW: Advection of W, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if wt variable, v2 is the inverse density
  ! array valid at w-points
  !
  subroutine ladvzw(n1,n2,n3,w,wt,wm,flx,v2)

    integer, intent (in) ::  n1,n2,n3
    real, intent (in)    :: wm(n1,n2,n3),w(n1,n2,n3),v2(n1)
    real, intent (inout) :: wt(n1,n2,n3)
    real, intent (out)   :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=1,n2
          flx(1,i,j)=0.
          flx(2,i,j)=wm(2,i,j)*.5*(0.583333*(w(2,i,j))-0.083333*(w(3,i,j)))
          do k=3,n1-1
             flx(k,i,j)=(0.583333*(w(k,i,j)+w(k-1,i,j))-0.083333            &
                  *(w(k+1,i,j)+w(k-2,i,j)))*.5*(wm(k,i,j)+wm(k-1,i,j))
          end do
          flx(n1,i,j)=0.

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k+1,i,j)-flx(k,i,j))*v2(k)
          end do
       end do
    end do

  end subroutine ladvzw

  ! ###########################
  ! 2nd order centered
  ! ###########################
  !
  ! ----------------------------------------------------------------------
  ! Ladvxu: Advection of U, by U, div is a scratch array, 
  ! tendencies are accumulated if ut variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxu2nd(n1,n2,n3,u,ut,flx,dxi)

    integer, intent (in) ::  n1,n2,n3
    real, intent(in)     :: u(n1,n2,n3),dxi
    real, intent(inout)  :: ut(n1,n2,n3)
    real, intent(out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=3,n2-1
          do k=2,n1-1
             flx(k,i,j) = (0.5*(u(k,i,j)+u(k,i-1,j)))*.5*(u(k,i,j)+u(k,i-1,j))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i+1,j)-flx(k,i,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxu2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVYU: Advection of U, by V, div is a scratch array, 
  ! tendencies are accumulated if fu variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyu2nd(n1,n2,n3,u,ut,v,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),v(n1,n2,n3),dyi
    real, intent (inout)  :: ut(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=2,n3-2
       do i=1,n2-1
          do k=2,n1-1
             flx(k,i,j)=(0.5*(u(k,i,j)+u(k,i,j+1)))*.5*(v(k,i,j)+v(k,i+1,j))
          end do
       end do
    end do

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i,j)-flx(k,i,j-1))*dyi
          end do
       end do
    end do

  end subroutine ladvyu2nd
  !
  ! ----------------------------------------------------------------------
  ! Subroutine ladvzu: Advection of U, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if ut variable, v1 is the inverse density
  ! array valid at thermo points
  !
  subroutine ladvzu2nd(n1,n2,n3,u,ut,wm,flx,v1)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),wm(n1,n2,n3),v1(n1)
    real, intent (inout)  :: ut(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=1,n2-1
          flx(1,i,j)=0.
          do k=2,n1-2
             flx(k,i,j)=(0.5*(u(k,i,j)+u(k+1,i,j)))*.5*(wm(k,i,j)+wm(k,i+1,j))
          end do
          flx(n1-1,i,j)=0.
          flx(n1,i,j)=0.
       end do

       do i=1,n2-1
          do k=2,n1-1
             ut(k,i,j)=ut(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*v1(k)
          end do
       end do
    end do

  end subroutine ladvzu2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVXV: Advection of V, by U, div is a scratch array, 
  ! tendencies are accumulated if vt variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxv2nd(n1,n2,n3,u,v,vt,flx,dxi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),v(n1,n2,n3),dxi
    real, intent (inout)  :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3-1
       do i=2,n2-2
          do k=2,n1-1
             flx(k,i,j)=(0.5*(v(k,i,j)+v(k,i+1,j)))*.5*(u(k,i,j)+u(k,i,j+1))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j)-flx(k,i-1,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxv2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVYV: Advection of V, by V, div is a scratch array, 
  ! tendencies are accumulated if vt variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyv2nd(n1,n2,n3,v,vt,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: v(n1,n2,n3),dyi
    real, intent (inout)  :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=3,n3-1
       do i=1,n2
          do k=2,n1-1
             flx(k,i,j)=(0.5*(v(k,i,j)+v(k,i,j-1)))*.5*(v(k,i,j)+v(k,i,j-1))
          end do
       end do
    end do

    do j=3,n3-2
       do i=1,n2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j+1)-flx(k,i,j))*dyi
          end do
       end do
    end do

  end subroutine ladvyv2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVZV: Advection of V, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if vt variable, v1 is the inverse density
  ! array valid at thermo points
  !
  subroutine ladvzv2nd(n1,n2,n3,v,vt,wm,flx,v1)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: v(n1,n2,n3),wm(n1,n2,n3),v1(n1)
    real, intent (inout)  :: vt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3-1
       do i=1,n2
          flx(1,i,j)=0.
          do k=2,n1-2
             flx(k,i,j)=(0.5*(v(k,i,j)+v(k+1,i,j)))*.5*(wm(k,i,j)+wm(k,i,j+1))
          end do
          flx(n1-1,i,j)=0.
          flx(n1,i,j)=0.
       end do

       do i=1,n2
          do k=2,n1-1
             vt(k,i,j)=vt(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*v1(k)
          end do
       end do
    end do

  end subroutine ladvzv2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVXW: Advection of W, by U, div is a scratch array, 
  ! tendencies are accumulated if wt variable, dxi is the inverse of 
  ! delta-x the grid spacing in the horizontal direction.
  !
  subroutine ladvxw2nd(n1,n2,n3,u,w,wt,flx,dxi)

    integer, intent (in)  :: n1,n2,n3
    real, intent (in)     :: u(n1,n2,n3),w(n1,n2,n3),dxi
    real, intent (inout)  :: wt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=2,n2-2
          do k=2,n1-2
             flx(k,i,j)=(0.5*(w(k,i,j)+w(k,i+1,j)))*.5*(u(k,i,j)+u(k+1,i,j))
          end do
       end do

       do i=3,n2-2
          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k,i,j)-flx(k,i-1,j))*dxi
          end do
       end do
    end do

  end subroutine ladvxw2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVYW: Advection of W, by V, div is a scratch array, 
  ! tendencies are accumulated if wt variable, dyi is the inverse of 
  ! delta-y.
  !
  subroutine ladvyw2nd(n1,n2,n3,vm,w,wt,flx,dyi)

    integer, intent (in)  ::  n1,n2,n3
    real, intent (in)     :: vm(n1,n2,n3),w(n1,n2,n3),dyi
    real, intent (inout)  :: wt(n1,n2,n3)
    real, intent (out)    :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=2,n3-2
       do i=1,n2
          do k=2,n1-2
             flx(k,i,j)=(0.5*(w(k,i,j)+w(k,i,j+1)))*.5*(vm(k,i,j)+vm(k+1,i,j))
          end do
       end do
    end do

    do j=3,n3-2
       do i=1,n2
          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k,i,j)-flx(k,i,j-1))*dyi
          end do
       end do
    end do

  end subroutine ladvyw2nd
  !
  ! ----------------------------------------------------------------------
  ! LADVZW: Advection of W, by W*rho, div is a scratch array, 
  ! tendencies are accumulated if wt variable, v2 is the inverse density
  ! array valid at w-points
  !
  subroutine ladvzw2nd(n1,n2,n3,w,wt,wm,flx,v2)

    integer, intent (in) ::  n1,n2,n3
    real, intent (in)    :: wm(n1,n2,n3),w(n1,n2,n3),v2(n1)
    real, intent (inout) :: wt(n1,n2,n3)
    real, intent (out)   :: flx(n1,n2,n3)

    integer :: i,j,k

    do j=1,n3
       do i=1,n2
          flx(1,i,j)=0.
          flx(2,i,j)=wm(2,i,j)*.5*(0.5*(w(2,i,j)))
          do k=3,n1-1
             flx(k,i,j)=(0.5*(w(k,i,j)+w(k-1,i,j)))*.5*(wm(k,i,j)+wm(k-1,i,j))
          end do
          flx(n1,i,j)=0.

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)-(flx(k+1,i,j)-flx(k,i,j))*v2(k)
          end do
       end do
    end do

  end subroutine ladvzw2nd

  !
  ! ----------------------------------------------------------------------
  ! ADVL_PREP: prepares two scratch arrays with the inverse
  ! densities as they locate on thermo levels (v1) and w-levels (v2)
  !
  subroutine advl_prep(n1,n2,n3,w,wm,dn0,dzi_t,dzi_m,v1,v2)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    ::  w(n1,n2,n3),dn0(n1),dzi_t(n1),dzi_m(n1)
    real, intent (out)   ::  wm(n1,n2,n3),v1(n1),v2(n1)

    integer :: k

    do k=1,n1-1
       wm(k,:,:)=w(k,:,:)*(dn0(k)+dn0(k+1))*.5
       v1(k)=dzi_t(k)/dn0(k)
       v2(k)=2.*dzi_m(k)/(dn0(k)+dn0(k+1))
    end do

  end subroutine advl_prep

end module advl
