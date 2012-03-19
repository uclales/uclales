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
module advf

  implicit none

  integer :: lmtr = 3 

contains
  !
  !----------------------------------------------------------------------
  ! subroutine fadvect: This is the driver for scalar advection.  It 
  ! advects using the average of the velocities at the current and past
  ! times.
  !
  subroutine fadvect

    use grid, only : a_up, a_vp, a_wp, a_sp, a_st, liquid, a_scr1, a_scr2,    &
         dn0 , nxp, nyp, nzp, nxyzp, dt, dzi_t, dzi_m, zt, dxi, dyi, level, nscl, &
         newvar, nstep

    use stat, only      : sflg, updtst
    use util, only      : atob, get_avg3

    real    :: v1da(nzp)
    integer :: n
    !
    ! diagnose liquid water flux
    !
    if (sflg .and. level > 1) then
       call atob(nxyzp,liquid,a_scr1)
       call atob(nxyzp,a_wp,a_scr2)
       call mamaos(nzp,nxp,nyp,a_scr2,liquid,a_scr1,zt,dzi_m,dn0,dt,.false.)
       call get_avg3(nzp,nxp,nyp,a_scr2,v1da)
       call updtst(nzp,'adv',0,v1da,1)
    end if
    !
    ! loop through the scalar table, setting iscp and isct to the 
    ! appropriate scalar pointer and do the advection, also add large
    ! scale subsidence.  Don't advect TKE here since it resides at a
    ! w-point
    !
    do n=4,nscl
       call newvar(n,istep=nstep)
       call atob(nxyzp,a_sp,a_scr1)

       call mamaos_y(nzp,nxp,nyp,a_vp,a_sp,a_scr1,dyi,dt)
       call mamaos_x(nzp,nxp,nyp,a_up,a_sp,a_scr1,dxi,dt)
       call atob(nxyzp,a_wp,a_scr2)
       call mamaos(nzp,nxp,nyp,a_scr2,a_sp,a_scr1,dzi_t,dzi_m,dn0,dt,.false.)

       if (sflg) then
          call get_avg3(nzp,nxp,nyp,a_scr2,v1da)
          call updtst(nzp,'adv',n-3,v1da,1)
       end if

       call advtnd(nzp,nxp,nyp,a_sp,a_scr1,a_st,dt)
    end do

  end subroutine fadvect
  !
  ! ----------------------------------------------------------------------
  ! Subroutine advtnd: Backs out the advective tendencies
  !
  subroutine advtnd(n1,n2,n3,varo,varn,tnd,dt)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)   :: varo(n1,n2,n3),varn(n1,n2,n3),dt
    real, intent(inout)  :: tnd(n1,n2,n3)

    real :: dti
    integer :: i,j,k

    dti=1./dt

    do j=3,n3-2
       do i=3,n2-2
          tnd(1,i,j)  = 0.
          do k=2,n1-1
             tnd(k,i,j)=tnd(k,i,j)+(varn(k,i,j)-varo(k,i,j))*dti
          end do
          tnd(n1,i,j) = 0.
       enddo
    enddo

  end subroutine advtnd
  !
  !---------------------------------------------------------------------- 
  ! Subroutine mamaos: An alternative second order flux limited scheme 
  ! written by Verica and Christiana as part of the MAMAOS program.  
  ! 
  ! July 21, 2003
  !
  subroutine mamaos(n1,n2,n3,w,scp0,scp,dzi_t,dzi_m,dn0,dt,lwpt)

    use mpi_interface, only : myid, appl_abort

    integer, intent (in)    :: n1,n2,n3
    real, intent (in)       :: scp0(n1,n2,n3)
    real, intent (in)       :: dn0(n1),dzi_t(n1),dzi_m(n1)
    real, intent (in)       :: dt
    logical, intent (in)    :: lwpt
    real, intent (inout)    :: w(n1,n2,n3),scp(n1,n2,n3)

    real    :: density(n1)   ! averaged density
    real    :: dzi_t_local(n1) ! grid spacing for scalars
    real    :: dzi_m_local(n1) ! grid spacing for velocity
    real    :: cfl(n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n1)         ! limiter    
    real    :: r(n1)         ! slope ratio
    real    :: wpdn(n1)      ! momentum: wp*density
    integer :: i, j, k, kp1, k1, k2
    integer :: gamma
    !
    ! initialize fields for use later
    !
    do k = 1, n1
       kp1 = min(k+1,n1)
       density(k) = 0.5 * (dn0(k) + dn0(kp1))
       if (lwpt) then
          dzi_t_local(k) = dzi_m(k)
          dzi_m_local(k) = dzi_t(kp1)
       else
          dzi_t_local(k) = dzi_t(k)
          dzi_m_local(k) = dzi_m(k)
       endif
    enddo

    do j = 3, n3-2
       do i = 3, n2-2
          !
          ! compute CFL and momentum
          !
          do k = 1, n1-1
             cfl(k)  = w(k,i,j) * dt * dzi_m_local(k)
             wpdn(k) = w(k,i,j) * density(k)
             if (abs(cfl(k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_z'
                call appl_abort (0)
             end if
          enddo
          !
          ! calculate the ratio of slopes
          !
          do k = 1, n1-1
             gamma = -sign(1.,cfl(k))
             if (abs(scp0(k+1,i,j)-scp0(k,i,j)) > spacing(scp0(k,i,j))) then
                k2 = max(1,k+gamma)
                k1 = min(n1,k+gamma+1)
                r(k) = (scp0(k1,i,j) - scp0(k2,i,j)) / &
                     (scp0(k+1,i,j) - scp0(k,i,j))
             else
                r(k) = 0.
             endif
          enddo
          !
          ! calculate the flux limiters
          !
          select case (lmtr)
          case (1) ! minmod
             do k = 1, n1-2
                C(k) = max(0., min(1., r(k)))
             enddo
          case(2)  ! superbee
             do k = 1, n1-2
                C(k) = max(0., min(1., 2.*r(k)), min(2., r(k)))
             enddo
          case(3)  ! mc
             do k = 1, n1-2
                C(k) = max(0., min(2.*r(k),(1.+r(k))/2., 2.))
             enddo
          case(4)  ! van Leer
             do k = 1, n1-2
                C(k) = (r(k) + abs(r(k)))/(1. + abs(r(k)))
             enddo
          case default ! no limiter
             do k = 1, n1-2
                C(k) = 1.0
             enddo
          end select

          w(1,i,j) = 0.
          w(n1-1,i,j) = 0.
          do k = 2, n1-2
             w(k,i,j) = 0.5 * wpdn(k) * (scp0(k+1,i,j)+scp0(k,i,j)) -  &
                  0.5 * (scp0(k+1,i,j)-scp0(k,i,j)) *                  &
                  ((1.-C(k))*abs(wpdn(k)) + wpdn(k)*cfl(k)*C(k))
          end do
          do k = 2,n1-1
             scp(k,i,j) = scp(k,i,j) - ((w(k,i,j)-w(k-1,i,j)) -        & 
                  scp0(k,i,j)*(wpdn(k)-wpdn(k-1))) *                   &
                  dt*dzi_t_local(k)/dn0(k)
          enddo

       enddo
    enddo

  end subroutine mamaos
  !
  !---------------------------------------------------------------------- 
  ! Subroutine mamaos_x: An alternative second order flux limited scheme 
  ! for advection in the x direction.  (adapted from mamaos)
  ! 
  ! September 3, 2003
  !
  subroutine mamaos_x(n1,n2,n3,u,scp0,scp,dxi,dt)

    use mpi_interface, only : myid, appl_abort

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dxi,dt
    real, intent (in)    :: scp0(n1,n2,n3),u(n1,n2,n3)
    real, intent (inout) :: scp(n1,n2,n3)

    real    :: cfl(n2,n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n2,n1)         ! limiter    
    real    :: r(n2,n1)         ! slope ratio
    real    :: scr(n2,n1)       ! flux scratch array
    integer :: i, j, k, i1, i2
    integer :: gamma
    !

    do j = 3,n3-2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       do k = 2, n1-1
          do i = 1,n2-1
             cfl(i,k)  = u(k,i,j) * dt * dxi
             scr(i,k)  = scp0(k,i+1,j)
             if (abs(cfl(i,k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_x'
                call appl_abort(0)
             end if
          end do
       end do     
          !
          ! calculate the ratio of slopes
          !
       do k = 2, n1-1 
          do i = 2,n2-2
             gamma = int(-sign(1.,cfl(i,k)))
             if (abs(scr(i,k) - scp0(k,i,j)) > spacing(scr(i,k))) then
                i2 = i+gamma
                i1 = i+gamma+1
                r(i,k) = (scp0(k,i1,j)-scp0(k,i2,j))/(scr(i,k)-scp0(k,i,j))
             else
                r(i,k) = 0.
             endif

             select case (lmtr)
             case (1) ! minmod
                C(i,k) = max(0., min(1., r(i,k)))
             case(2)  ! superbee
                C(i,k) = max(0., min(1., 2.*r(i,k)), min(2., r(i,k)))
             case(3)  ! mc
                C(i,k) = max(0., min(2.*r(i,k),(1.+r(i,k))/2., 2.))
             case(4)  ! van Leer
                C(i,k) = (r(i,k) + abs(r(i,k)))/(1. + abs(r(i,k)))
             case default ! no limiter
                C(i,k) = 1.0
             end select

             scr(i,k) = 0.5 * u(k,i,j) * (scr(i,k)+scp0(k,i,j)) -      &
                  0.5 * (scr(i,k)-scp0(k,i,j)) *                        &
                  ((1.-C(i,k))*abs(u(k,i,j)) + u(k,i,j)*cfl(i,k)*C(i,k))
          end do

          do i = 3,n2-2
             scp(k,i,j) = scp(k,i,j) - ((scr(i,k)-scr(i-1,k)) -         &
                  scp0(k,i,j)*(u(k,i,j)-u(k,i-1,j)))*dt*dxi
          enddo
       enddo

    enddo

  end subroutine mamaos_x
  !
  !---------------------------------------------------------------------- 
  ! Subroutine mamaos_y: An alternative second order flux limited scheme 
  ! for advection in the y direction.  (adapted from mamaos)
  ! 
  ! September 3, 2003
  !
  subroutine mamaos_y(n1,n2,n3,v,scp0,scp,dyi,dt)

    use mpi_interface, only : myid, appl_abort

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dyi,dt
    real, intent (in)    :: scp0(n1,n2,n3),v(n1,n2,n3)
    real, intent (inout) :: scp(n1,n2,n3)

    real    :: cfl(n3,n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n3,n1)         ! limiter    
    real    :: r(n3,n1)         ! slope ratio
    real    :: scr(n3,n1)       ! flux scratch array
    integer :: i, j, k, j1, j2
    integer :: gamma
    !

    do i = 1, n2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       do k = 2, n1-1
          do j = 1,n3-1
             cfl(j,k)  = v(k,i,j) * dt * dyi
             scr(j,k)  = scp0(k,i,j+1)
             if (abs(cfl(j,k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_y'
                call appl_abort(0)
             end if
          end do
       end do
          !
          ! calculate the ratio of slopes
          !
       do k = 2, n1-1 
          do j = 2,n3-2
             gamma = int(-sign(1.,cfl(j,k)))
             if (abs(scr(j,k) - scp0(k,i,j)) > spacing(scr(j,k))) then
                j2 = j+gamma
                j1 = j+gamma+1
                r(j,k) = (scp0(k,i,j1)-scp0(k,i,j2))/(scr(j,k)-scp0(k,i,j))
             else
                r(j,k) = 0.
             endif

             select case (lmtr)
             case (1) ! minmod
                C(j,k) = max(0., min(1., r(j,k)))
             case(2)  ! superbee
                C(j,k) = max(0., min(1., 2.*r(j,k)), min(2., r(j,k)))
             case(3)  ! mc
                C(j,k) = max(0., min(2.*r(j,k),(1.+r(j,k))/2., 2.))
             case(4)  ! van Leer
                C(j,k) = (r(j,k) + abs(r(j,k)))/(1. + abs(r(j,k)))
             case default ! no limiter
                C(j,k) = 1.0
             end select

             scr(j,k) = 0.5 * v(k,i,j) * (scr(j,k)+scp0(k,i,j)) -      &
                  0.5 * (scr(j,k)-scp0(k,i,j)) *                        &
                  ((1.-C(j,k))*abs(v(k,i,j)) + v(k,i,j)*cfl(j,k)*C(j,k))
          end do

          do j = 3,n3-2
             scp(k,i,j) = scp(k,i,j) - ((scr(j,k)-scr(j-1,k)) -         &
                  scp0(k,i,j)*(v(k,i,j)-v(k,i,j-1)))*dt*dyi
          enddo
       enddo

    enddo

  end subroutine mamaos_y

end module advf
