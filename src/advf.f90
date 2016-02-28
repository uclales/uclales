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

    use grid, only : a_up, a_vp, a_wp, a_sp, a_st, liquid, a_scr1, a_scr2, a_scr3, &
         dn0 , nxp, nyp, nzp, nxyzp, dt, dzi_t, dzi_m, zt, dxi, dyi, level, nscl, &
         newvar, newvar_ls, nstep, levset, a_beta,a_betabar,ls_cellkind, &
         a_xp, ls_q0, ls_q1

    use stat, only      : sflg, updtst
    use util, only      : atob, azero, get_avg3

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
    if (levset >= 1) then
       do n=4,nscl
          ! setting tendency pointers
          call newvar(n,istep=nstep)

          ! saving old state
          call atob(nxyzp,a_xp(:,:,:,n),a_scr1)

          call azero(nxyzp, a_scr3)

          ! this is lstep 2
          call mamaos_y_ls(nzp,nxp,nyp,a_vp,ls_q0(:,:,:,n),a_scr1,a_betabar,2,dyi,dt)
          call mamaos_x_ls(nzp,nxp,nyp,a_up,ls_q0(:,:,:,n),a_scr1,a_betabar,2,dxi,dt)
          call atob(nxyzp,a_wp,a_scr2)
          call mamaos_ls(nzp,nxp,nyp,a_scr2,ls_q0(:,:,:,n),a_scr1,a_betabar,2,dzi_t,dzi_m,dn0,dt,.false.)
          a_scr3(:,:,:) = a_scr3(:,:,:) + a_scr2(:,:,:)

          ! this is lstep 1
          call mamaos_y_ls(nzp,nxp,nyp,a_vp,ls_q1(:,:,:,n),a_scr1,a_betabar,1,dyi,dt)
          call mamaos_x_ls(nzp,nxp,nyp,a_up,ls_q1(:,:,:,n),a_scr1,a_betabar,1,dxi,dt)
          call atob(nxyzp,a_wp,a_scr2)
          call mamaos_ls(nzp,nxp,nyp,a_scr2,ls_q1(:,:,:,n),a_scr1,a_betabar,1,dzi_t,dzi_m,dn0,dt,.false.)
          a_scr3(:,:,:) = a_scr3(:,:,:) + a_scr2(:,:,:)

          
          call advection_fix(nzp, nxp, nyp, a_scr1, ls_q0(:,:,:,n), ls_q1(:,:,:,n), n)
           
          if (sflg) then
             call get_avg3(nzp,nxp,nyp,a_scr3,v1da)
             call updtst(nzp,'adv',n-3,v1da,1)
          end if
          
          ! compute tendencies based on new (a_xp) and old state (a_scr1).
          call advtnd_ls(nzp,nxp,nyp,a_xp(:,:,:,n),a_scr1,a_st,dt,ls_cellkind)
       end do
    else
       do n=4,nscl
         call newvar(n,istep=nstep)
         call atob(nxyzp,a_sp,a_scr1)

         call mamaos_y(nzp,nxp,nyp,a_vp,a_sp,a_scr1,dyi,dt)
         call mamaos_x(nzp,nxp,nyp,a_up,a_sp,a_scr1,dxi,dt)
         call atob(nxyzp,a_wp,a_scr2)
         call mamaos(nzp,nxp,nyp,a_scr2,a_sp,a_scr1,dzi_t,dzi_m,dn0,dt,.false.)

         ! update advective fluxes of th_l (n=4) and q (n=5)
         if (sflg) then
            call get_avg3(nzp,nxp,nyp,a_scr2,v1da)
            call updtst(nzp,'adv',n-3,v1da,1)
         end if

         call advtnd(nzp,nxp,nyp,a_sp,a_scr1,a_st,dt)
      end do
    end if

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


  subroutine advtnd_ls(n1,n2,n3,varo,varn,tnd,dt,ls_cellkind)

    integer, intent(in) :: n1,n2,n3,ls_cellkind(n1,n2,n3)
    real, intent(in)    :: varo(n1,n2,n3),varn(n1,n2,n3),dt
    real, intent(inout) :: tnd(n1,n2,n3)

    real :: dti
    integer :: i,j,k

    dti=1./dt

    do j=3,n3-2
        do i=3,n2-2
          tnd(1,i,j) = 0.
          do k=2,n1-1
            tnd(k,i,j)=tnd(k,i,j)+(varn(k,i,j)-varo(k,i,j))*dti
          end do
          tnd(n1,i,j) = 0.
        enddo
    enddo

  end subroutine advtnd_ls
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

          do k = 1, n1-2
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
  ! Subroutine mamaos_ls: Second-order flux limited scheme for x advection
  !   based on LeVeque's piecewise linear reconstruction (Finite Volume Methods
  !   for Hyperbolic Problems. p. 113 ff)
  !
  subroutine mamaos_ls(n1,n2,n3,w,scp0,scp,beta,lstep,dzi_t,dzi_m,dn0,dt,lwpt)

    use mpi_interface, only : myid, appl_abort

    integer, intent (in)    :: n1,n2,n3
    real, intent (in)       :: scp0(n1,n2,n3),beta(n1,n2,n3,3)
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
    integer :: i, j, k, kp1, k1, k2, lstep
    integer :: gamma
    real    :: w1, w2              ! these correspond to beta or 1-beta, depending on lstep at i and i-1, resp.
                                   ! w1 = weight at k, w2 = weight at (k-1)
    real, parameter :: a(2) = (/ 0.0,  1.0 /),&
                       b(2) = (/ 1.0, -1.0 /)
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

          !vertical fluxes weighted by vertical face farctions
          w(1,i,j) = 0.
          w(n1-1,i,j) = 0.
          do k = 2, n1-2
              w(k,i,j) = ( a(lstep) + b(lstep)*beta(k,i,j,1) )    &
                       * ( min(0.0,wpdn(k)) * scp0(k+1,i,j)       &
                         + max(0.0,wpdn(k)) * scp0(k,i,j)         &
                         + 0.5*abs(wpdn(k))*(1.-abs(cfl(k)))*C(k) &
                              *(scp0(k+1,i,j) - scp0(k,i,j))      &
                         )
          end do
          do k = 2, n1-1 
            scp(k,i,j) = scp(k,i,j) - (w(k,i,j)-w(k-1,i,j))*dt*dzi_t_local(k)/dn0(k)
          enddo
       enddo
    enddo

  end subroutine mamaos_ls
  !
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
                if (myid == 0) print *, '  ABORTING: mamaos_x',dt,u(k,i,j)
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
  ! Subroutine mamaos_x_ls: Second-order flux limited scheme for x advection
  !   based on LeVeque's piecewise linear reconstruction (Finite Volume Methods
  !   for Hyperbolic Problems. p. 113 ff)
  !
  subroutine mamaos_x_ls(n1,n2,n3,u,scp0,scp,beta,lstep,dxi,dt)
    use mpi_interface, only : myid, appl_abort
    use grid, only : a_alpha,a_beta

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dxi,dt
    real, intent (in)    :: scp0(n1,n2,n3),u(n1,n2,n3),beta(n1,n2,n3,3)
    real, intent (inout) :: scp(n1,n2,n3)

    real    :: cfl(n2,n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n2,n1)         ! limiter    
    real    :: r(n2,n1)         ! slope ratio
    real    :: scr(n2,n1)       ! flux scratch array
    integer :: i, j, k, i1, i2, lstep
    integer :: gamma
    real    :: w1, w2              ! correspond to beta or 1-beta, depending on lstep at i and i-1, reps..
    real, parameter :: a(2) = (/ 0.0,  1.0 /),&
                       b(2) = (/ 1.0, -1.0 /)
   
    do j = 3,n3-2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       do k = 2, n1-1
          do i = 1,n2-1
             cfl(i,k)  = u(k,i,j) * dt * dxi
             scr(i,k)  = scp0(k,i+1,j)
             if (abs(cfl(i,k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_x_ls'
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

             scr(i,k) = min(0.,u(k,i,j)) * scp0(k,i+1,j)                        &
                      + max(0.,u(k,i,j)) * scp0(k,i,  j)                        &
                      + 0.5 * abs(u(k,i,j)) * (1. - abs(cfl(i,k))) * C(i,k)     &
                      * (scp0(k,i+1,j) - scp0(k,i,j))
          end do

          do i = 3,n2-2
            ! w1 and w2 are the appropriate cell face fractions 'beta' or
            !   '1-beta' depending on the current fluid phase (lstep) for the
            !   right (w1), and left (w2) cell face.
            w1 = a(lstep) + b(lstep)*beta(k,i  ,j,2)
            w2 = a(lstep) + b(lstep)*beta(k,i-1,j,2)
            scp(k,i,j) = scp(k,i,j) - ((w1*scr(i,k) - w2*scr(i-1,k)))*dt*dxi
          enddo
       enddo

    enddo
  end subroutine mamaos_x_ls
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
  !
  !---------------------------------------------------------------------- 
  ! Subroutine mamaos_y_ls: Second-order flux limited scheme for x advection
  !   based on LeVeque's piecewise linear reconstruction (Finite Volume Methods
  !   for Hyperbolic Problems. p. 113 ff)
  !
  subroutine mamaos_y_ls(n1,n2,n3,v,scp0,scp,beta,lstep,dyi,dt)
    use mpi_interface, only : myid, appl_abort

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: dyi,dt
    real, intent (in)    :: scp0(n1,n2,n3),v(n1,n2,n3),beta(n1,n2,n3,3)
    real, intent (inout) :: scp(n1,n2,n3)

    real    :: cfl(n3,n1)       ! cfl numbers at the interface (staggered)
    real    :: C(n3,n1)         ! limiter    
    real    :: r(n3,n1)         ! slope ratio
    real    :: scr(n3,n1)       ! flux scratch array
    integer :: i, j, k, j1, j2, lstep
    integer :: gamma
    real    :: w1, w2           ! correspond to beta or 1-beta, depending on lstep at i and i-1, reps..
    real, parameter :: a(2) = (/ 0.0,  1.0 /),&
                       b(2) = (/ 1.0, -1.0 /)

    do i = 3, n2-2
       !
       ! compute CFL and scr array for down-grid value of scalar
       !
       do k = 2, n1-1
          do j = 1,n3-1
             cfl(j,k)  = v(k,i,j) * dt * dyi
             scr(j,k)  = scp0(k,i,j+1)
             if (abs(cfl(j,k)) > 1.0) then
                if (myid == 0) print *, '  ABORTING: mamaos_y_ls'
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
             scr(j,k) = min(0.,v(k,i,j)) * scp0(k,i,j+1)                        &
                      + max(0.,v(k,i,j)) * scp0(k,i,j)                          &
                      + 0.5*abs(v(k,i,j)) * (1. - abs(cfl(j,k))) * C(j,k)       &
                      * (scp0(k,i,j+1) - scp0(k,i,j))
          end do
          do j = 3,n3-2
            ! w1 and w2 are the appropriate cell face fractions 'beta' or
            !   '1-beta' depending on the current fluid phase (lstep) for the
            !   right (w1), and left (w2) cell face.
            w1 = a(lstep) + b(lstep)*beta(k,i,j,  3)
            w2 = a(lstep) + b(lstep)*beta(k,i,j-1,3)
            scp(k,i,j) = scp(k,i,j) - (w1*scr(j,k) - w2*scr(j-1,k))*dt*dyi
          enddo
       enddo

    enddo

  end subroutine mamaos_y_ls

  subroutine advection_fix(nzp, nxp, nyp, q, q0, q1, n)
  
    use grid, only : ls_G, ls_nG, levset, a_alphatilde, deltax, deltay, dzi_t, &
        dn0, th_change, qt_change, filprf
    
    real, intent(in), dimension(nzp,nxp,nyp) :: q0, q1
    integer, intent(in)                      :: nzp, nxp, nyp, n
    real, intent(inout)                      :: q(nzp,nxp,nyp) 
    integer :: k, i, j, l, typ
    real :: temp, dxy, dV_k
    real :: q_bottom = 0.0
    real :: q_top = 1.0
    
    dxy = deltax * deltay
    if (trim(filprf)=='wobble') then
        q_bottom = 0.0
        q_top = 1.0
    else
        q_bottom = 1e-3
        q_top = 0.0
    end if
    
    select case(n)
    
    case default
    do l = 1, ls_nG
        k = ls_G(l, 1) 
        i = ls_G(l, 2) 
        j = ls_G(l, 3)
        dV_k     = dxy / dzi_t(k)
        q(k,i,j) = a_alphatilde(k,i,j) * q1(k,i,j) + (1.0 - a_alphatilde(k,i,j)) * q0(k,i,j)
    end do
    
    case(4)
    do l = 1, ls_nG
        k = ls_G(l, 1) 
        i = ls_G(l, 2) 
        j = ls_G(l, 3)
        dV_k      = dxy / dzi_t(k)
        temp      = a_alphatilde(k,i,j) * q1(k,i,j) + (1.0 - a_alphatilde(k,i,j)) * q0(k,i,j)
        th_change = th_change + dn0(k) * (temp - q(k,i,j)) * dV_k
        q(k,i,j)  = temp
    end do
    
    case(5)
    do l = 1, ls_nG
        k = ls_G(l, 1) 
        i = ls_G(l, 2) 
        j = ls_G(l, 3)
        dV_k      = dxy / dzi_t(k)
        temp      = a_alphatilde(k,i,j) * q1(k,i,j) + (1.0 - a_alphatilde(k,i,j)) * q0(k,i,j)
        qt_change = qt_change + dn0(k) * (temp - q(k,i,j)) * dV_k
        q(k,i,j)  = temp
    end do
    end select
    
  end subroutine advection_fix

end module advf
