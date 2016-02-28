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
module lset

  use grid, only : bandw,gkmin,gkmax,fdmthd,rimthd
  use reinit, only : reinitialize
  implicit none

  integer          :: something = 42
  real, parameter  :: eps       = 1e-15
  real, parameter  :: maxiter   = 20   ! max number of spreading iterations

  !Run-mode for the computation of level set velocities
  ! mode = 0: normal mode (= interpolate velocities)
  ! mode = 1: test mode   (= generate a velocity field for test cases)
  ! mode = 2: test mode   (= oscillating vertical velocity field)
  integer :: mode = 0

  !FD method to use to discretize level set equation (\phi_t + V.\nabla \phi = 0)
  ! fdmthd = 0: first-order upwind
  ! fdmthd = 1: 2nd-order polynomial upwind
  ! fdmthd = 2: second-order ENO

  !Reinitialization method
  ! rimthd = 0: no reinitialization
  ! rimthd = 1: column-wise vertical direct reinitialization
  ! rimthd = 2: Sethian's equation
 
  !New VolumeFraction Routines? 1=yes, 0=no.
  integer :: vfmode = 1
contains
  !
  ! ls_ghost: generates ghost fluid
  !
  subroutine ls_ghost(conservative)
    use grid, only : nscl,nxyzp,levset,ls_q0,ls_q1,ls_G_currently, ls_Gm, ls_Gp, ls_Gpc, ls_UCN, &
       a_xp,a_sp,ls_nG_currently,ls_nGm,ls_nGp,ls_nGpc,ls_nUCN

    logical, intent(in) :: conservative
    integer :: n

    if (levset>=1) then
        ls_q1 = a_xp
        ls_q0 = a_xp
        if (conservative) then
            do n = 4, nscl
                call ls_spread2(ls_q1(:,:,:,n),ls_Gm,ls_nGm,type='thrm',sig=-1)
                call ls_reconstruct(ls_q0(:,:,:,n), ls_q1(:,:,:,n), ls_G_currently,ls_nG_currently)
                call ls_spread2(ls_q0(:,:,:,n),ls_Gpc,ls_nGpc,type='thrm',sig=+1)
            end do
        else
            do n = 4, nscl
                call ls_spread2(ls_q1(:,:,:,n),ls_Gm,ls_nGm,type='thrm',sig=-1)
                call ls_spread2(ls_q0(:,:,:,n),ls_Gp,ls_nGp,type='thrm',sig=+1)
            end do
        end if
        !ls_q1(:,:,:,:) = 1.0
        !ls_q0(:,:,:,:) = 0.0
    end if
  end subroutine ls_ghost
  !
  ! ls_reconstruct: compute ghost fluid state 0 (bottom) based on ghost fluid 1 and
  !     volume fractions imersed in the respective fluids. On entering this routine,
  !     S0 contains a copy of the original cell averages. On leaving the routine,
  !     S0 will carry the reconstructed value in the mixed cells.
  !
  subroutine ls_reconstruct(S0, S1, G_currently, nG_currently)

    use grid, only : nxp,nyp,nzp,nxyzp,a_alpha

    real,    dimension(nzp,nxp,nyp), intent(in)    :: S1
    real,    dimension(nzp,nxp,nyp), intent(inout) :: S0
    integer, dimension(nxyzp,3),     intent(in)    :: G_currently
    integer, intent(in) :: nG_currently
    integer             :: k,i,j,l
    real :: tmp

    do l = 1, nG_currently
        call unpackIndices(G_currently, l, k, i, j)
        tmp       = S0(k,i,j)
        S0(k,i,j) = (S0(k,i,j) - S1(k,i,j)*a_alpha(k,i,j))/(1.-a_alpha(k,i,j)+tiny(1.0))
    end do

  end subroutine ls_reconstruct
  !
  ! ls_advance: advances the level set in time
  !
  subroutine ls_advance(div, time)
    use grid, only : levset, a_scr1, write_anal, nxp, nyp, nzp, a_phi, a_phitilde

    real, intent(in) :: div
    real, optional, intent(in) :: time

    if (levset>=1) then
        if (present(time)) then
            call ls_mom2nod(div, time)
        else
            call ls_mom2nod(div)
        end if
        call ls_advect !phi(n+1) has now been stored in grid.a_phitilde
        call ls_characterize(a_phi, a_phitilde)
        call ls_fractions
        call ls_crossing(ctime=a_scr1)
        call ls_abchange(ctime=a_scr1)
        call ls_ghost(conservative=.false.)
    end if
  end subroutine ls_advance
  !
  ! ----------------------------------------------------------------------
  ! ls_finalize: extracts level set information and reinitializes
  !
  subroutine ls_finalize
  
  use grid, only : a_phi, xm, ym, zm, dxi, dyi, dzi_t, nxp, nyp, nzp, nxyzp, levset, &
                   nstep, newphi
  use util, only : atob
    !
    ! reinitialize at the end of each time step dt
    !
    if (levset>=1) then
        !call ls_characterize
        !call ls_normals('bckw',ls_nx,ls_ny,ls_nz)
        if (nstep == 3) call reinitialize(rimthd, a_phi, xm, ym, zm, nxp, nyp, nzp, gkmin, gkmax, dxi, dyi, dzi_t)
    end if

  end subroutine ls_finalize
  !
  !----------------------------------------------------------------------
  ! subroutine ls_advect: This is the driver for advection of the level set
  ! function. A 2nd-order non-limited finite-differences upwind scheme is
  ! used that uses interpolated velocities at cell corners.
  !
  subroutine ls_advect
  
    use grid, only : a_up, a_vp, a_wp, a_sp, a_st, liquid, a_scr1, a_scr2,    &
         dn0 , nxp, nyp, nzp, nxyzp, dt, dzi_t, dzi_m, rz, zt, zm, dxi, dyi,      &
         level, nscl, newphi, nstep, a_phitilde

    use util, only      : atob,sclrset
    use advf, only      : advtnd

    call newphi(istep=nstep)
    call atob(nxyzp,a_sp,a_phitilde)
    call ls_fdupwind(a_sp,a_phitilde,dt)
    ! Need to set BCs here, a_phitilde is used in a following step (ls_abchange)
    call sclrset('cnst',nzp,nxp,nyp,a_phitilde)
    call advtnd(nzp,nxp,nyp,a_sp,a_phitilde,a_st,dt)

  end subroutine ls_advect
  !
  !---------------------------------------------------------------------- 
  ! subroutine ls_fdupwind: 2nd-order 3d finite-differences upwind scheme for
  ! uniform grids. 'scp0' points to the original field of phi. 'scp' is a copy
  ! of the original field of phi.
  !
  subroutine ls_fdupwind(scp0,scp,dt)
  
    use mpi_interface, only : myid, appl_abort
    use grid, only          : nzp,nxp,nyp,dxi,dyi,dzi_t,dzi_m,zm,rz,ym,ls_up, &
                              ls_vp,ls_wp, gkmin, gkmax
    use util, only : cyclics, cyclicc

    real, intent (in)    :: dt
    real, intent (in)    :: scp0(nzp,nxp,nyp)
    real, intent (inout) :: scp(nzp,nxp,nyp)

    integer :: i, j, k, kk, ii, jj, req(16)
    real    :: a,b,c,cx,cy,cz, dQ1x,dQ1y,dQ1z, dQ2x,dQ2y, dQ2z
    real    :: D1x(nzp,nxp,nyp), D2x(nzp,nxp,nyp), absD2x(nzp,nxp,nyp)
    real    :: D1y(nzp,nxp,nyp), D2y(nzp,nxp,nyp), absD2y(nzp,nxp,nyp)
    real    :: D1z(nzp,nxp,nyp), D2z(nzp,nxp,nyp), absD2z(nzp,nxp,nyp)
    real    :: absgradphi
   
    select case(fdmthd)

    case default ! first order upwind
    do j = 3, nyp-2
        do i = 3,nxp-2
            do k = gkmin+1, gkmax-1
                !
                ! x advection
                !
                scp(k,i,j) = scp(k,i,j) - dt*dxi*( &
                             max(ls_up(k,i,j),0.)*(scp0(k,i,j)   - scp0(k,i-1,j)) &
                           + min(ls_up(k,i,j),0.)*(scp0(k,i+1,j) - scp0(k,i,j)  ) )
                !
                ! y advection 
                !
                scp(k,i,j) = scp(k,i,j) - dt*dyi*( &
                             max(ls_vp(k,i,j),0.)*(scp0(k,i,j)   - scp0(k,i,j-1)) &
                           + min(ls_vp(k,i,j),0.)*(scp0(k,i,j+1) - scp0(k,i,j)  ) )
                !
                ! z advection
                !
                scp(k,i,j) = scp(k,i,j) - dt*( &
                             max(ls_wp(k,i,j),0.)*dzi_t(k)  *(scp0(k,i,j)   - scp0(k-1,i,j)) &
                           + min(ls_wp(k,i,j),0.)*dzi_t(k-1)*(scp0(k+1,i,j) - scp0(k,i,j)  ) )
            end do
        end do
    end do

    case(1) ! second-order polynomial upwind.
    do j = 3, nyp-2
        do i = 3,nxp-2
            do k = gkmin+2, gkmax-2
                !
                ! x advection (non-conservative)
                !
                scp(k,i,j) = scp(k,i,j) -  dt*0.5*dxi*( &
                             min(ls_up(k,i,j),0.)*(-3*scp0(k,i,j)+4*scp0(k,i+1,j)-1*scp0(k,i+2,j)) &
                           + max(ls_up(k,i,j),0.)*( 3*scp0(k,i,j)-4*scp0(k,i-1,j)+1*scp0(k,i-2,j)) )
                !
                ! y advection (non-conservative)
                !
                scp(k,i,j) = scp(k,i,j) -  dt*0.5*dyi*( &
                             min(ls_vp(k,i,j),0.)*(-3*scp0(k,i,j)+4*scp0(k,i,j+1)-1*scp0(k,i,j+2)) &
                           + max(ls_vp(k,i,j),0.)*( 3*scp0(k,i,j)-4*scp0(k,i,j-1)+1*scp0(k,i,j-2)) )
                !
                ! z advection (non-conservative)
                !
                if (ls_wp(k,i,j) >= 0.) then
                  a = ( 2*zm(k) - (zm(k-1) + zm(k-2)) ) / ( (zm(k)   - zm(k-1))*(zm(k)   - zm(k-2))  )
                  b = ( 2*zm(k) - (zm(k)   + zm(k-2)) ) / ( (zm(k-1) - zm(k)  )*(zm(k-1) - zm(k-2))  )
                  c = ( 2*zm(k) - (zm(k)   + zm(k-1)) ) / ( (zm(k-2) - zm(k)  )*(zm(k-2) - zm(k-1))  )
                else
                  a = ( 2*zm(k) - (zm(k+1) + zm(k+2)) ) / ( (zm(k)   - zm(k+1))*(zm(k)   - zm(k+2))  )
                  b = ( 2*zm(k) - (zm(k)   + zm(k+2)) ) / ( (zm(k+1) - zm(k)  )*(zm(k+1) - zm(k+2))  )
                  c = ( 2*zm(k) - (zm(k)   + zm(k+1)) ) / ( (zm(k+2) - zm(k)  )*(zm(k+2) - zm(k+1))  )
                end if
                scp(k,i,j) = scp(k,i,j) -  dt*( &
                             min(ls_wp(k,i,j),0.)*( a*scp0(k,i,j)+b*scp0(k+1,i,j)+c*scp0(k+2,i,j)) &
                           + max(ls_wp(k,i,j),0.)*( a*scp0(k,i,j)+b*scp0(k-1,i,j)+c*scp0(k-2,i,j)) )
            end do
        end do
    end do

    ! second-order ENO
    ! This method assumes equidistant grids in all three dimensions.
    ! On stretched grids, numerical order will be reduced.
    case(2)
    ! computing devided differences table
    do j = 1, nyp-1
        do i = 1,nxp-1
            do k = 1,nzp-1
                D1x(k,i,j) = scp0(k,i+1,j) - scp0(k,i,j) !D1_{k+1/2} phi
                D1y(k,i,j) = scp0(k,i,j+1) - scp0(k,i,j) !D1_{k+1/2} phi
                D1z(k,i,j) = scp0(k+1,i,j) - scp0(k,i,j) !D1_{k+1/2} phi
            end do
        end do
    end do

    do j = 2, nyp-1
        do i = 2,nxp-1
            do k = 2,nzp-1
                D2x(k,i,j)    = 0.5*(D1x(k,i,j) - D1x(k,i-1,j))
                absD2x(k,i,j) = abs(D2x(k,i,j))
                D2y(k,i,j)    = 0.5*(D1y(k,i,j) - D1y(k,i,j-1))
                absD2y(k,i,j) = abs(D2y(k,i,j))
                D2z(k,i,j)    = 0.5*(D1z(k,i,j) - D1z(k-1,i,j))
                absD2z(k,i,j) = abs(D2z(k,i,j))
            end do
        end do
    end do
    do j = 3, nyp-2
        do i = 3,nxp-2
            do k = 3,nzp-2
                if (ls_up(k,i,j) > 0.) then
                  ii = i-1
                else
                  ii = i
                end if
                if (ls_vp(k,i,j) > 0.) then
                  jj = j-1
                else
                  jj = j
                end if
                if (ls_wp(k,i,j) > 0.) then
                  kk = k-1
                else
                  kk = k
                end if
                if (absD2x(k,ii,j) < absD2x(k,ii+1,j)) then
                  cx = D2x(k,ii,j)
                else
                  cx = D2x(k,ii+1,j)
                end if
                if (absD2y(k,i,jj) < absD2y(k,i,jj+1)) then
                  cy = D2y(k,i,jj)
                else
                  cy = D2y(k,i,jj+1)
                end if
                if (absD2z(kk,i,j) < absD2z(kk+1,i,j)) then
                  cz = D2z(kk,i,j)
                else
                  cz = D2z(kk+1,i,j)
                end if
                dQ1x = D1x(k,ii,j)
                dQ2x = cx*(2*(i-ii)-1)
                dQ1y = D1y(k,i,jj)
                dQ2y = cy*(2*(j-jj)-1)
                dQ1z = D1z(kk,i,j)
                dQ2z = cz*(2*(k-kk)-1)
                absgradphi = sqrt( (0.5*dxi*(scp0(k,i+1,j) - scp0(k,i-1,j)))**2 +     &
                                   (0.5*dyi*(scp0(k,i,j+1) - scp0(k,i,j-1)))**2 +     &
                      (1./(zm(k+1)-zm(k-1))*(scp0(k+1,i,j) - scp0(k-1,i,j)))**2 )
                scp(k,i,j) = scp(k,i,j) - dt*( &
                     ls_up(k,i,j) * ( (dQ1x + dQ2x*dxi) * dxi )           + &
                     ls_vp(k,i,j) * ( (dQ1y + dQ2y*dyi) * dyi )           + &
                     ls_wp(k,i,j) * ( (dQ1z + dQ2z*dzi_t(k)) * dzi_t(k))  - &
                     1.*dxi*dxi*(scp0(k,i+1,j) - 2*scp0(k,i,j) + scp0(k,i-1,j) + &
                                 scp0(k,i,j+1) - 2*scp0(k,i,j) + scp0(k,i,j-1) + &
                                 scp0(k+1,i,j) - 2*scp0(k,i,j) + scp0(k-1,i,j) ) * absgradphi)
            end do
        end do
    end do
    !
    ! lin. extrapolation at top and bottom boundaries
    !
    do j=3,nyp-2
        do i=3,nxp-2
           scp(2,i,j) = 2*scp(3,i,j) - scp(4,i,j)
           scp(1,i,j) = 2*scp(2,i,j) - scp(3,i,j)

           scp(nzp-1,i,j) = 2*scp(nzp-2,i,j) - scp(nzp-3,i,j)
           scp(nzp,i,j)   = 3*scp(nzp-2,i,j) - 2*scp(nzp-3,i,j)
        end do
    end do

    end select
    !
    ! Boundary conditions
    !
    do j=3,nyp-2
        do i=3,nxp-2
           scp(gkmin+1,i,j) = 2*scp(gkmin+2,i,j) - scp(gkmin+3,i,j)
           scp(gkmin  ,i,j) = 2*scp(gkmin+1,i,j) - scp(gkmin+2,i,j)
           do k = 1, max(gkmin-1, 1)
               scp(k,i,j) = scp(gkmin,i,j)
           end do

           scp(gkmax-1,i,j) = 2*scp(gkmax-2,i,j) - scp(gkmax-3,i,j)
           scp(gkmax,i,j)   = 3*scp(gkmax-2,i,j) - 2*scp(gkmax-3,i,j)
           do k = min(gkmax+1, nzp), nzp
               scp(k,i,j) = scp(gkmax,i,j)
           end do
        end do
    end do
    call cyclics(nzp,nxp,nyp,scp,req)
    call cyclicc(nzp,nxp,nyp,scp,req)

  end subroutine ls_fdupwind

  !
  ! ls_normals: computes normal vectors of the level set. Normal vectors live
  !    on the grid nodes, as does phi itslef..
  !
  subroutine ls_normals(type,nx,ny,nz)
    use grid, only : nxp, nyp, nzp, dxi, dyi, zm, dzi_m, dzi_t, a_phi, ls_absdphi, &
                     ls_phidx,ls_phidy,ls_phidz
    use util, only : sclrset, cyclicc, cyclics
    use mpi_interface, only : myid, appl_abort

    !type = 'cntr' central; 'bckw' backward; 'upwd' upwind; 'thrm' for normals at thermal points
    character(len=4), intent(in) :: type 
    real,dimension(nzp,nxp,nyp), intent(inout) :: nx,ny,nz
    real    :: ai
    integer :: k,i,j,req(16)

    ls_phidx(:,:,:) = 1e10
    ls_phidy(:,:,:) = 1e10
    ls_phidz(:,:,:) = 1e10
    nx(:,:,:) = 0.0
    ny(:,:,:) = 0.0
    nz(:,:,:) = 0.0

    select case(type)
    case default
      if (myid==0) print *, '  ABORTING: normals discretization not specified' 
      call appl_abort(0)

    case('cntr')
    do j=3,nyp-2
      do i=3,nxp-2
        do k=gkmin+1, gkmax-1
          !not 2nd-order if dz not equidistant
          ls_phidx(k,i,j) = (a_phi(k  ,i+1,j  ) - a_phi(k  ,i-1,j  ))*0.5*dxi
          ls_phidy(k,i,j) = (a_phi(k  ,i  ,j+1) - a_phi(k  ,i  ,j-1))*0.5*dyi
          ls_phidz(k,i,j) = (a_phi(k+1,i  ,j  ) - a_phi(k-1,i  ,j  ))/(zm(k+1) - zm(k-1))
        end do
      end do
    end do

    case('bckw')
    do j=2,nyp
      do i=2,nxp
        do k=gkmin+1, gkmax-1
          ls_phidx(k,i,j) = (a_phi(k,i,j) - a_phi(k  ,i-1,j  ))*dxi
          ls_phidy(k,i,j) = (a_phi(k,i,j) - a_phi(k  ,i  ,j-1))*dyi
          ls_phidz(k,i,j) = (a_phi(k,i,j) - a_phi(k-1,i  ,j  ))*dzi_t(k) 
        end do
      end do
    end do

    case('upwd')
    do j=3,nyp-1
      do i=3,nxp-1
        do k=gkmin+1, gkmax-1
          if ( ((a_phi(k,i,j) - a_phi(k,i-1,j)).lt.0.).and. &
               ((a_phi(k,i,j) - a_phi(k,i+1,j)).lt.0.) ) then
            ls_phidx(k,i,j) = 0.
          else
            ls_phidx(k,i,j) = dxi*(sign(1.,a_phi(k,i+1,j)-a_phi(k,i-1,j))  &
              *max(a_phi(k,i,j)-a_phi(k,i+1,j), a_phi(k,i,j)-a_phi(k,i-1,j)) )
          end if
          if ( ((a_phi(k,i,j) - a_phi(k,i,j-1)).lt.0.).and. &
               ((a_phi(k,i,j) - a_phi(k,i,j-1)).lt.0.) ) then
            ls_phidy(k,i,j) = 0.
          else
            ls_phidy(k,i,j) = dyi*(sign(1.,a_phi(k,i,j+1)-a_phi(k,i,j-1))  &
              *max(a_phi(k,i,j)-a_phi(k,i,j+1), a_phi(k,i,j)-a_phi(k,i,j-1)) )
          end if
          if ( ((a_phi(k,i,j) - a_phi(k-1,i,j)).lt.0.).and. &
               ((a_phi(k,i,j) - a_phi(k+1,i,j)).lt.0.) ) then
            ls_phidz(k,i,j) = 0.
          else
            ls_phidz(k,i,j) = dzi_t(k)*(sign(1.,a_phi(k+1,i,j)-a_phi(k-1,i,j))  &
              *max(a_phi(k,i,j)-a_phi(k+1,i,j), a_phi(k,i,j)-a_phi(k-1,i,j)) )
          end if
        end do
      end do
    end do
    
    case('thrm')
    do j=2,nyp
      do i=2,nxp
        do k=gkmin+1, gkmax-1
          ls_phidx(k,i,j) = 0.25*dxi  * (a_phi(k  ,i,j  ) - a_phi(k  ,i-1,j  ) &
                                      +  a_phi(k  ,i,j-1) - a_phi(k  ,i-1,j-1) &
                                      +  a_phi(k-1,i,j  ) - a_phi(k-1,i-1,j  ) &
                                      +  a_phi(k-1,i,j-1) - a_phi(k-1,i-1,j-1) ) 
          ls_phidy(k,i,j) = 0.25*dyi  * (a_phi(k  ,i,j  ) - a_phi(k  ,i  ,j-1) &
                                      +  a_phi(k  ,i-1,j) - a_phi(k  ,i-1,j-1) &
                                      +  a_phi(k-1,i  ,j) - a_phi(k-1,i  ,j-1) &
                                      +  a_phi(k-1,i-1,j) - a_phi(k-1,i-1,j-1) ) 
          ls_phidz(k,i,j) = 0.25*dzi_m(k)*(a_phi(k  ,i  ,j)   - a_phi(k-1,i  ,j  ) &
                                         + a_phi(k  ,i  ,j-1) - a_phi(k-1,i  ,j-1) &
                                         + a_phi(k  ,i-1,j)   - a_phi(k-1,i-1,j  ) &
                                         + a_phi(k  ,i-1,j-1) - a_phi(k-1,i-1,j-1) )
        end do
      end do
    end do

    end select
    ls_phidx(gkmin,:,:) = ls_phidx(gkmin+1,:,:)
    ls_phidy(gkmin,:,:) = ls_phidy(gkmin+1,:,:)
    ls_phidz(gkmin,:,:) = ls_phidz(gkmin+1,:,:)
    ls_phidx(gkmax,:,:) = ls_phidx(gkmax-1,:,:)
    ls_phidy(gkmax,:,:) = ls_phidy(gkmax-1,:,:)
    ls_phidz(gkmax,:,:) = ls_phidz(gkmax-1,:,:)
    call cyclics(nzp,nxp,nyp,ls_phidx,req)
    call cyclicc(nzp,nxp,nyp,ls_phidx,req)
    call cyclics(nzp,nxp,nyp,ls_phidy,req)
    call cyclicc(nzp,nxp,nyp,ls_phidy,req)
    call cyclics(nzp,nxp,nyp,ls_phidz,req)
    call cyclicc(nzp,nxp,nyp,ls_phidz,req)

    do j=3,nyp-2
      do i=3,nxp-2
        do k=gkmin, gkmax
          ls_absdphi(k,i,j) = ( ls_phidx(k,i,j)**2.0 &
                              + ls_phidy(k,i,j)**2.0 &
                              + ls_phidz(k,i,j)**2.0 )**0.5

          ai = 1./(ls_absdphi(k,i,j)+tiny(ai))
          nx(k,i,j) = ls_phidx(k,i,j)*ai
          ny(k,i,j) = ls_phidy(k,i,j)*ai
          nz(k,i,j) = ls_phidz(k,i,j)*ai
        end do
      end do
    end do

    ! boundary conditions (sclrset includes call of cyclics)
    call sclrset('grad',nzp,nxp,nyp,ls_absdphi,dzi_m)
    call sclrset('grad',nzp,nxp,nyp,nx,dzi_m)
    call sclrset('grad',nzp,nxp,nyp,ny,dzi_m)
    call sclrset('grad',nzp,nxp,nyp,nz,dzi_m)

  end subroutine ls_normals
  !
  ! ls_mom2nod: interpolates bulk velocities (w,u,v) to grid nodes.
  !
  subroutine ls_mom2nod(div, time)
    use grid, only : ls_up,ls_vp,ls_wp, a_up, a_vp, a_wp, nxp, nyp, nzp, rzmt, &
      dzi_m, dzi_t, zm, xm, ls_UCN, ls_nUCN
    use defs, only : pi
    use util, only : velset
    use mpi_interface, only : myid, appl_abort,cyclics, cyclicc

    integer :: k,i,j,extravel,use_subsidence,req(16)
    real    :: z, ang, r, omega, xc, zc, U0, mydiv
    real, intent(in) :: div
    real,optional,intent(in) :: time
    !
    ! extravel = 1: After setting up level set velocities, extrapolate velocities
    !     in the normal direction.
    !
    extravel = 0
    use_subsidence = 1

    select case (mode)
    case default ! normal mode
      if (use_subsidence==1) then
         mydiv = div
      else
         mydiv = 0.0
      end if

      ls_up (:,:,:) = 0.0
      ls_vp (:,:,:) = 0.0
      ls_wp (:,:,:) = 0.0

      do j = 3, nyp-2
      do i = 3, nxp-2
      do k = gkmin+1, gkmax-1
          ls_up(k,i,j) = 0.5*(a_up(k,i,j) + a_up(k,i,j+1))             &
                       + 0.25*rzmt(k)*(a_up(k+1,i,j) + a_up(k+1,i,j+1) &
                                      -a_up(k,i,j)   - a_up(k,i,j+1)   )
          ls_vp(k,i,j) = 0.5*(a_vp(k,i,j) + a_vp(k,i+1,j))             &
                       + 0.25*rzmt(k)*(a_vp(k+1,i,j) + a_vp(k+1,i+1,j) &
                                      -a_vp(k,i,j)   - a_vp(k,i+1,j)   )
          ls_wp(k,i,j) = 0.25*(a_wp(k,i,j)     + a_wp(k,i+1,j)         &
                              +a_wp(k,i+1,j+1) + a_wp(k,i,j+1) )       &
                       - mydiv * zm(k)
      end do
      end do
      end do

    case(1) !test mode: solid body rotation
      xc = 0.0
      zc = 50.0
      omega = -2.*pi/628.
      ls_vp = 0.0
      do k=1,nzp
        z = zm(k)-zc
        do i=1,nxp
          r = ((xm(i)-xc)**2. + (z)**2.)**0.5
          if (xm(i)-xc > -tiny(1.)) then
            ang = atan(z/(xm(i) - xc + tiny(1.0) ))
          else
            ang = pi + atan(z/(xm(i) - xc + tiny(1.0) ))
          end if
          ls_up(k,i,:) = -sin(ang)*omega*r
          ls_wp(k,i,:) =  cos(ang)*omega*r
        end do
      end do

    case(2)
    if (present(time)) then
        U0 = 1.0
        ls_up = 0.
        ls_vp = 0.
        ls_wp = U0*sin(2.*pi*(1./40.)*time)
    else
        if (myid==0) print *,"ABORTING: time not passed for time-dependent level set test "
        call appl_abort(0)
    end if
    end select
    ! 
    ! Extrapolate level set transport velocities. Level set velocities are
    !   extrapolated only in the level set band, and only on nodes of cells not
    !   intersected by the zero level set. This is the set of nodes 'ls_UCN'.
    !
    if (extravel==1) then
        call ls_spread2(ls_up,ls_UCN,ls_nUCN,type='upwd')
        call ls_spread2(ls_vp,ls_UCN,ls_nUCN,type='updw')
        call ls_spread2(ls_wp,ls_UCN,ls_nUCN,type='updw')
    end if
    !
    ! Set periodic boundary conditions
    !
    call velset(nzp,nxp,nyp,ls_up,ls_vp,ls_wp)

  end subroutine ls_mom2nod
  !
  ! ls_characterize: This routine does three things:
  !         1) Find cells which are intersected by the current zero level set (cut
  !              cells)
  !         2) Based on the location of cut cells, define a (small) band of cells
  !              around cut cells, 2*bandw wide vertically.
  !         3) Integrate face and volume fractions, beta and alpha, resp., for the
  !              cut cells. alpha and beta indicate the the fractions immersed in
  !              dry fluid of the free atmospehre.
  !
  subroutine ls_characterize(phi, phitilde)
    use grid, only : nxp,nyp,nzp,ls_G, ls_nG,ls_B,ls_nB,  &
        ls_UCN,ls_nUCN,a_alpha,ls_Gtilde, ls_nGtilde, ls_cellkind,   &
        ls_nodekind,a_alpha,a_beta,ls_Gm,ls_Gp,ls_nGm,ls_nGp,ls_CN,ls_nCN, &
        a_alphatilde, a_betabar, ls_G_currently, ls_nG_currently, ls_Gpc, ls_nGpc
    use mpi_interface, only  : int_scalar_par_max
    use util, only : markerset
    
    real, dimension(nzp, nxp, nyp), intent(in) :: phi, phitilde
    real    :: m
    integer :: k, i, j, l, ll, kmin, kmax
    logical :: typ, typtilde
    
    ls_cellkind(:,:,:) = 42
    ls_nodekind(:,:,:) = 42
    ls_nG              = 0
    ls_nG_currently    = 0
    kmax               = 1
    kmin               = nzp
    !
    ! Note: ordering of the loop indices is important:
    !    1. It is efficient in Fortran to have the inner loop iterate trough the
    !       first index. (Arrays are indexed (k,i,j).)
    !    2. The entrainment routine relies on ls_G to carry indices with
    !       vertically neighbouring cells next to each other.
    !
    do j=3, nyp-2
    do i=3, nxp-2
    do k=2, nzp-1
        typ      = ls_celltype(nzp, nxp, nyp, k, i, j, phi)
        typtilde = ls_celltype(nzp, nxp, nyp, k, i, j, phitilde)
        if (typ.or.typtilde) then
            ! These cells are cut during the LES time step [t, t+dt].
            call packIndices(ls_G, ls_nG, k, i, j)
            kmin = max(2,     min(kmin, k-bandw)) 
            kmax = min(nzp-1, max(kmax, k+bandw)) 

            if (typ) then
                ! These cells are cut currently (at beginning of the time step)
                call packIndices(ls_G_currently, ls_nG_currently, k, i, j)
            end if
        end if
    end do

    end do
    end do

    if (ls_nG < 1) then
        stop "  ls_characterize: No cut cells found in this MPI task. I cannot handle this."
    end if 
    !
    ! define level set band
    !
    if (ls_nG == 0) then
        kmin = 0
        kmax = 0
    end if
    ! set kmin,kmax to match entire domain if no band is to be used
    if (bandw==0) then
        gkmin = 1
        gkmax = nzp
    else
        call int_scalar_par_max(kmax,gkmax)
        call int_scalar_par_max(-kmin,gkmin)
        gkmin = -gkmin
    end if
    !
    ! find band cells and store indices
    !
    ls_nB = 0 ! number of band cells
    do j=3,nyp-2
    do i=3,nxp-2
    do k=gkmin,gkmax
        call packIndices(ls_B, ls_nB, k, i, j)
    end do
    end do
    end do
    !
    ! Mark 'moist' and 'dry' domain by negative and positive cellkind, 
    ! respectively. Band cells have magnitude 1, non-band cells have magnitude
    ! 2. In summary: 
    !     -2, -1 -> moist (boundary layer)
    !     +2, +1 -> dry  (free atmosphere)
    !
    do k=2, nzp-1
      if ( (k<=gkmax).and.(k>=gkmin) ) then
        m = 1.
      else
        m = 2.
      end if
      do j=3,nyp-2
        do i=3,nxp-2
          ls_cellkind(k,i,j) = int( sign( m, phi(k,i,j) ) )
          ls_nodekind(k,i,j) = ls_cellkind(k,i,j)
        end do
      end do
    end do
    !
    ! Mark nodes of cells in the band
    !
    do l=1,ls_nB
        call unpackIndices(ls_B, l, k, i, j)
        ll = int(sign(1.,phi(k,i,j)))
        ls_nodekind(k,i,j)       = ll
        ls_nodekind(k,i-1,j)     = ll
        ls_nodekind(k,i,j-1)     = ll
        ls_nodekind(k,i-1,j-1)   = ll
        ls_nodekind(k-1,i,j)     = ll
        ls_nodekind(k-1,i-1,j)   = ll
        ls_nodekind(k-1,i,j-1)   = ll
        ls_nodekind(k-1,i-1,j-1) = ll
    end do
    !
    ! Now, correctly mark all cut cells '0'. Also, mark all 8 nodes of these
    ! cells '0'.
    !
    do l=1, ls_nG
        call unpackIndices(ls_G, l, k, i, j)
        ls_cellkind(k,i,j)       = 0
        ls_nodekind(k,i,j)       = 0
        ls_nodekind(k,i-1,j)     = 0
        ls_nodekind(k,i,j-1)     = 0
        ls_nodekind(k,i-1,j-1)   = 0
        ls_nodekind(k-1,i,j)     = 0
        ls_nodekind(k-1,i-1,j)   = 0
        ls_nodekind(k-1,i,j-1)   = 0
        ls_nodekind(k-1,i-1,j-1) = 0
    end do

    call markerset(nzp,nxp,nyp,ls_nodekind)
    call markerset(nzp,nxp,nyp,ls_cellkind)
    !
    ! find band nodes, excluding nodes of cut cells (ls_UCN). Also establish sets
    ! for ghost fluid extrapolation: lower band including cut cells (ls_Gm) and
    ! upper band excluding cut cells (ls_Gp).
    !
    ls_nCN  = 0
    ls_nUCN = 0
    ls_nGm  = 0
    ls_nGp  = 0
    ls_nGpc = 0
    do j=3,nyp-2
    do i=3,nxp-2
    do k=gkmin,gkmax
        if (abs(ls_nodekind(k,i,j))==1) then
            call packIndices(ls_UCN, ls_nUCN, k, i, j)
        end if
        if (ls_nodekind(k,i,j)==0) then !is a node of a cut cell
            call packIndices(ls_CN, ls_nCN, k, i, j)
        end if
        !
        ! Define band cells where the bottom ghost fluid is to be generated 
        ! (inside BL, currently cut cells and below, phi <= 0).
        !
        if ((ls_cellkind(k,i,j)==0).or.(ls_cellkind(k,i,j)==-1)) then
            call packIndices(ls_Gm, ls_nGm, k, i, j)
        end if
        !
        ! These are all band cells in the free atmosphere including cut cells
        ! (phi >= 0).
        !
        if ((ls_cellkind(k,i,j)==1).or.(ls_cellkind(k,i,j)==0)) then
            call packIndices(ls_Gp, ls_nGp, k, i, j)
            ! For conservative reconstruction: The following are all band cells
            ! in the free atmosphere excluding cut cells (phi > 0).
            if (ls_cellkind(k,i,j)==1) then
                call packIndices(ls_Gpc, ls_nGpc, k, i, j)
            end if
        end if
    end do
    end do
    end do
  end subroutine ls_characterize


  subroutine ls_fractions
  use grid, only : ls_G, ls_nG, a_phi, a_phitilde, a_alpha, a_alphatilde, &
    a_beta, nzp, nxp, nyp

    call ls_integrateAlpha(ls_G,ls_nG,a_phitilde,a_alphatilde,nzp,nxp,nyp)
    call ls_integrateAlpha(ls_G,ls_nG,a_phi,     a_alpha,     nzp,nxp,nyp)
    call ls_integrateBeta (ls_G,ls_nG,a_phi,     a_beta,      nzp,nxp,nyp)
  end subroutine ls_fractions
  !
  ! ls_crossing: find interface crossing times in all nodes of cut cells.
  !     Crossing times are measured relative to the current model time.
  !
  subroutine ls_crossing(ctime)
    use grid, only : dt, ls_nCN, ls_CN, a_phi, a_phitilde,nxp,nyp,nzp
    use util, only : sclrset
    
    real, intent(out), dimension(nzp,nxp,nyp) :: ctime
    integer :: k,i,j,l

    ctime(:,:,:) = -42.e16

    do l = 1, ls_nCN
      call unpackIndices(ls_CN, l, k, i, j)
      ctime(k,i,j) = a_phi(k,i,j) * dt &
                   / ( a_phi(k,i,j) - a_phitilde(k,i,j) + spacing(1.0) )
    end do

    call sclrset('cnst',nzp,nxp,nyp,ctime)
  end subroutine ls_crossing


  subroutine ls_integrateAlpha(cell_set,number_of_cells,phi,alpha,nzp,nxp,nyp)
    use util, only : sclrset

    integer, intent(in) :: number_of_cells, nzp, nxp, nyp
    integer, dimension(nxp*nyp*nzp,3), intent(in) :: cell_set
    real, dimension(nzp,nxp,nyp),   intent(in)    :: phi
    real, dimension(nzp,nxp,nyp),   intent(inout) :: alpha
    real, dimension(0:3) :: gb, gt
    integer              :: k, i, j, l
    !
    ! Initialize all alpha values according to the state of the rear-top-right
    ! cell node. This gives the correct values for uncut cells which are fully
    ! immersed in one phase. The cut (during the LES time step) cell alphas are
    ! overwritten below.
    !
    do j = 3, nyp - 2
    do i = 3, nxp - 2
    do k = 1, nzp
        ! alpha = 0 in moist phase, alpha = 1 in dry phase
        alpha(k,i,j) = 0.5*( 1. + sign(1., phi(k,i,j)) )
    end do
    end do
    end do

    do l = 1, number_of_cells
        call unpackIndices(cell_set, l, k, i, j)
        !
        ! level set function values on bottom face (z-)
        !
        gb(0) = phi(k-1,i-1,j-1)
        gb(1) = phi(k-1,i-1,j  )
        gb(2) = phi(k-1,i  ,j  )
        gb(3) = phi(k-1,i  ,j-1)
        !
        ! level set function values on top face (z+)
        !
        gt(0) = phi(k,i-1,j-1)
        gt(1) = phi(k,i-1,j  )
        gt(2) = phi(k,i  ,j  )
        gt(3) = phi(k,i  ,j-1)

        select case(vfmode)
        case(1)
            alpha(k,i,j) = NewVolumeFraction3D(gb, gt)
        case default
            alpha(k,i,j) = 1.-VolumeFraction3D(gb, gt)
        end select
    end do

    call sclrset('cnst',nzp,nxp,nyp,alpha)
  end subroutine ls_integrateAlpha

  subroutine ls_integrateBeta(cell_set,number_of_cells,phi,beta,nzp,nxp,nyp)
    use util, only : sclrset

    integer, intent(in) :: number_of_cells, nzp, nxp, nyp
    integer, dimension(nxp*nyp*nzp,3), intent(in) :: cell_set
    real, dimension(nzp,nxp,nyp),   intent(in)    :: phi
    real, dimension(nzp,nxp,nyp,3), intent(inout) :: beta
    real, dimension(0:3) :: gt
    integer              :: k, i, j, l
    !
    ! Initialize all beta values according to the state of the rear-top-right
    ! cell node. This gives the correct values for uncut cells which are fully
    ! immersed in one phase. The cut (during the LES time step) cell betas are
    ! overwritten below.
    !
    do j = 3, nyp - 2
    do i = 3, nxp - 2
    do k = 2, nzp
        ! beta = 0 in moist phase, beta = 1 in dry phase
        beta(k,  i,  j,  :) = 0.5*( 1. + sign(1., phi(k,i,j)) )
        beta(k-1,i,  j,  1) = beta(k,i,j,1) 
        beta(k,  i-1,j,  2) = beta(k,i,j,1)
        beta(k,  i,  j-1,3) = beta(k,i,j,1)
    end do
    end do
    end do

    do l = 1, number_of_cells
        call unpackIndices(cell_set, l, k, i, j)
        !
        ! +- z
        !
        gt(0) = phi(k-1,i-1,j-1)
        gt(1) = phi(k-1,i-1,j  )
        gt(2) = phi(k-1,i  ,j  )
        gt(3) = phi(k-1,i  ,j-1)
        beta(k-1,i,j,1) = NewVolumeFraction2D(gt)
        gt(0) = phi(k,i-1,j-1)
        gt(1) = phi(k,i-1,j  )
        gt(2) = phi(k,i  ,j  )
        gt(3) = phi(k,i  ,j-1)
        beta(k,i,j,1) = NewVolumeFraction2D(gt)
        !
        ! +- x
        !
        gt(0) = phi(k-1,i  ,j-1)
        gt(1) = phi(k  ,i  ,j-1)
        gt(2) = phi(k  ,i  ,j  )
        gt(3) = phi(k-1,i  ,j  )
        beta(k,i,j,2) = NewVolumeFraction2D(gt)
        gt(0) = phi(k-1,i-1,j-1)
        gt(1) = phi(k  ,i-1,j-1)
        gt(2) = phi(k  ,i-1,j  )
        gt(3) = phi(k-1,i-1,j  )
        beta(k,i-1,j,2) = NewVolumeFraction2D(gt)
        !
        ! +- y
        !
        gt(0) = phi(k-1,i-1,j  )
        gt(1) = phi(k  ,i-1,j  )
        gt(2) = phi(k  ,i  ,j  )
        gt(3) = phi(k-1,i  ,j  )
        beta(k,i,j,3) = NewVolumeFraction2D(gt)
        gt(0) = phi(k-1,i-1,j-1)
        gt(1) = phi(k  ,i-1,j-1)
        gt(2) = phi(k  ,i  ,j-1)
        gt(3) = phi(k-1,i  ,j-1)
        beta(k,i,j-1,3) = NewVolumeFraction2D(gt)
    end do
    call sclrset('cnst',nzp,nxp,nyp,beta(:,:,:,1))
    call sclrset('cnst',nzp,nxp,nyp,beta(:,:,:,2))
    call sclrset('cnst',nzp,nxp,nyp,beta(:,:,:,3))
  end subroutine ls_integrateBeta
  !
  ! ls_abchange: computes mean values of beta during the current time step based
  !     on crossing times at cut cell nodes.
  !
  subroutine ls_abchange(ctime)
    use grid, only : a_phi,a_phitilde, a_betabar,nxp,nyp,nzp,dt,ls_G,ls_nG,ls_nodekind,a_beta
    use util, only : sclrset
    real, dimension(nzp,nxp,nyp), intent(in) :: ctime
    real, dimension(0:5)    :: t,b
    real, dimension(1:4)    :: phi
    real                    :: dti,dti05,bsum
    integer, dimension(0:5) :: tseq
    integer, dimension(1:4) :: kk,ii,jj
    integer                 :: k,i,j,l,m,node,dir
    
    dti   = 1./dt
    dti05 = dti*0.5

    a_betabar(:,:,:,:) = a_beta(:,:,:,:)

    do dir=1,3 ! directions
      select case(dir)
      case(1) !z
        kk = (/0,0,0,0/)
        ii = (/1,1,0,0/)
        jj = (/1,0,0,1/)
      case(2) !x
        kk = (/1,0,0,1/)
        ii = (/0,0,0,0/)
        jj = (/1,1,0,0/)
      case(3) !y
        kk = (/1,0,0,1/)
        ii = (/1,1,0,0/)
        jj = (/0,0,0,0/)
      end select
      do l=1,ls_nG
        k = ls_G(l,1)
        i = ls_G(l,2)
        j = ls_G(l,3)
        ! crossing times at nodes of current face
        t(0)    = 0.0
        tseq(0) = 0
        do m=1,4
            t(m)    = min(dt, max(0.0, ctime(k-kk(m),i-ii(m),j-jj(m)) ) )
            tseq(m) = m
        end do
        t(5)    = dt
        tseq(5) = 5
        !
        ! sort crossing times. Note, 'sort' assumes indices range from 0..3,
        !     covert to 1..4.
        !
        call sort(tseq(1:4),t(1:4))
        tseq(1:4) = tseq(1:4)+1
        !
        ! Integrate the mean alphas and betas
        !
        m = 0
        do node=1,4
          phi(node) = a_phi(k-kk(node),i-ii(node),j-jj(node)) &
            + (a_phitilde(k-kk(node),i-ii(node),j-jj(node)) - a_phi(k-kk(node),i-ii(node),j-jj(node)) ) &
            * dti*t(tseq(m))
        end do
        b(m) = NewVolumeFraction2D(phi)
        !
        m = 5
        do node=1,4
          phi(node) = a_phi(k-kk(node),i-ii(node),j-jj(node)) &
            + (a_phitilde(k-kk(node),i-ii(node),j-jj(node)) - a_phi(k-kk(node),i-ii(node),j-jj(node)) ) &
            * dti*t(tseq(m))
        end do
        b(m) = NewVolumeFraction2D(phi)
        !
        ! checking for time collisions 
        !
        if(t(1)==t(2).and.t(2)==t(3).and.t(3)==t(4)) then
            ! In case crossing times are the same at each node, beta jumps, so
            !     we define the following. Note, in this case, betas are exactly
            !     0.0 or 1.0.
            b(1) = b(0)
            b(2) = b(0)
            b(3) = b(5)
            b(4) = b(5)
        else
            ! This is the general cas, beta does not jump.
            do m=1,4
                do node=1,4
                    phi(node) = a_phi(k-kk(node),i-ii(node),j-jj(node)) &
                      + (a_phitilde(k-kk(node),i-ii(node),j-jj(node)) - a_phi(k-kk(node),i-ii(node),j-jj(node)) ) &
                      * dti*t(tseq(m))
                end do
                if (t(tseq(m))>t(tseq(m-1))) then
                    !
                    ! All four nodes are imersed in one phase.
                    !
                    if      (phi(1).ge.0.0.and.phi(2).ge.0.0  .and. &
                             phi(3).ge.0.0.and.phi(4).ge.0.0) then
                        b(m) = 1.0
                    else if (phi(1).le.0.0.and.phi(2).le.0.0  .and. &
                             phi(3).le.0.0.and.phi(4).le.0.0) then
                        b(m) = 0.0
                    !
                    ! Normal case 
                    !
                    else 
                        b(m) = NewVolumeFraction2D(phi)
                    end if
                else 
                    b(m) = b(m-1)
                end if
            end do
        end if
        !
        ! Integratie beta(t) using the trapezoidal rule
        !
        bsum = 0.0
        do m=1,5
            bsum = bsum + (b(m)+b(m-1)) * (t(tseq(m))-t(tseq(m-1)))
        end do
        a_betabar(k,i,j,dir) = dti05*bsum
      end do
    end do
    call sclrset('cnst',nzp,nxp,nyp,a_betabar(:,:,:,1))
    call sclrset('cnst',nzp,nxp,nyp,a_betabar(:,:,:,2))
    call sclrset('cnst',nzp,nxp,nyp,a_betabar(:,:,:,3))
  end subroutine ls_abchange
  !
  ! ls_spread2: extrapolates scalar field S along the level set normals in
  !   either the \phi+ and the \phi- direction, based on the sign of 'sig'.
  !   
  subroutine ls_spread2(S,passed_cells, passed_ncells,type,sig)
    use mpi_interface, only : cyclics, cyclicc, appl_abort, myid
    use grid, only          : nxp,nyp,nzp,nxyzp,nstep,dxi,dyi,dzi_t,ls_nx,ls_ny,  &
      ls_nz,ls_phidx,ls_phidy,ls_phidz,a_phi,ls_nodekind,ls_cellkind

    integer :: k,i,j,l,ncells,kupw,iupw,jupw,iter,req(16)
    real    :: res, dtau, cfl, absn_i, nxfact, nyfact, nzfact
    real, dimension(nzp,nxp,nyp), intent(inout) :: S
    real, dimension(nzp,nxp,nyp) :: Sdx,Sdy,Sdz,dS,nx,ny,nz
    integer, dimension(nxyzp,3),optional,intent(in) :: passed_cells
    integer,                    optional,intent(in) :: passed_ncells,sig
    integer, dimension(nxyzp,3) :: cell_indices !internal index array
    character(4), intent(in) :: type
    real :: avg 
    logical :: not_two_arguments
    !
    ! check if both arguments present
    !
    not_two_arguments = xor(present(passed_cells), present(passed_ncells))
    if (not_two_arguments) then
        if (myid==0) print *,"ABORTING: wrong number of arguments in ls_spread2"
        call appl_abort(0)
    end if !both present or both not present
    !
    ! define the cells which are to be updated
    !
    if (present(passed_cells)) then
        cell_indices(:,:) = passed_cells(:,:) ! indices internal to this routine
        ncells            = passed_ncells     ! number of internal indices
    else
        !
        ! process the entire domain
        !
        l=0
        do j=3,nyp-2
        do i=3,nxp-2
        do k=2,nzp-1
            l = l+1
            cell_indices(l,1) = k
            cell_indices(l,2) = i
            cell_indices(l,3) = j
        end do
        end do
        end do
        ncells = l
    end if
    
    cfl       = 0.5
    call ls_normals(type,nx,ny,nz)
    dtau = cfl*abs(min(1./dxi/maxval(nx), 1./dyi/maxval(ny), &
                       1./maxval(dzi_t)/maxval(nz)) )
    !
    ! use average of source part of the level set band as as initial condition.
    !
    avg = 0.
    iter = 0
    do j=3,nyp-2
    do i=3,nxp-2
    do k=gkmin,gkmax
        if (ls_cellkind(k,i,j)==-sig) then
            avg = avg + S(k,i,j)
            iter = iter+1
        end if
    end do
    end do
    end do

    avg = avg/(iter+tiny(1.0))

    do l=1,ncells
        call unpackIndices(cell_indices, l, k, i, j)
        S(k,i,j) = avg
    end do
    call cyclics(nzp,nxp,nyp,S,req)
    call cyclicc(nzp,nxp,nyp,S,req)
    !print *,"phi = ",sig, "IC of q = ",avg
    !
    ! iterate extrapolation equation
    !
    res = 1.0
    iter = 0
    do while (iter<maxiter)
      iter = iter + 1
      res  = 0.0
      dS   = 0.0
      !
      ! Compute gradient of S
      !
      Sdx = 1e16
      Sdy = 1e16
      Sdz = 1e16
      do j=2,nyp-1
      do i=2,nxp-1
      do k=gkmin+1, gkmax
          Sdx(k,i,j) = (S(k,i,j) - S(k,i-1,j))*dxi
          Sdy(k,i,j) = (S(k,i,j) - S(k,i,j-1))*dyi
          Sdz(k,i,j) = (S(k,i,j) - S(k-1,i,j))*dzi_t(k)
      end do
      end do
      end do
      Sdx(gkmin,:,:) = Sdx(gkmin+1,:,:)
      Sdy(gkmin,:,:) = Sdy(gkmin+1,:,:)
      Sdz(gkmin,:,:) = Sdz(gkmin+1,:,:)
      Sdx(gkmax,:,:) = Sdx(gkmax-1,:,:)
      Sdy(gkmax,:,:) = Sdy(gkmax-1,:,:)
      Sdz(gkmax,:,:) = Sdz(gkmax-1,:,:)
      call cyclics(nzp,nxp,nyp,Sdx,req)
      call cyclicc(nzp,nxp,nyp,Sdx,req)
      call cyclics(nzp,nxp,nyp,Sdy,req)
      call cyclicc(nzp,nxp,nyp,Sdy,req)
      call cyclics(nzp,nxp,nyp,Sdz,req)
      call cyclicc(nzp,nxp,nyp,Sdz,req)
      !
      ! accumulate tendencies dS
      !
      do l=1,ncells
          call unpackIndices(cell_indices, l, k, i, j)
          select case(type)
          case('thrm')
            dS(k,i,j) = dS(k,i,j) + dtau*( min(sig*nx(k,i,j),0.)*Sdx(k,i+1,j) &
                                          +max(sig*nx(k,i,j),0.)*Sdx(k,i,  j))
            dS(k,i,j) = dS(k,i,j) + dtau*( min(sig*ny(k,i,j),0.)*Sdy(k,i,j+1) &
                                          +max(sig*ny(k,i,j),0.)*Sdy(k,i,j  ))
            dS(k,i,j) = dS(k,i,j) + dtau*( min(sig*nz(k,i,j),0.)*Sdz(k+1,i,j) &
                                          +max(sig*nz(k,i,j),0.)*Sdz(k,  i,j))
          case('upwd')
          !
          ! x direction
          !
          if (Sdx(k,i,j)<Sdx(k,i+1,j)) then
            dS(k,i,j) = dS(k,i,j) + dtau*min(sign(1.,a_phi(k,i,j))*nx(k,i,j)*Sdx(k,i,j),   &
                                             sign(1.,a_phi(k,i,j))*nx(k,i,j)*Sdx(k,i+1,j)  )
          else
            dS(k,i,j) = dS(k,i,j) + dtau*max(sign(1.,a_phi(k,i,j))*nx(k,i,j)*Sdx(k,i,j),   &
                                             sign(1.,a_phi(k,i,j))*nx(k,i,j)*Sdx(k,i+1,j)  )
          end if
          !
          ! y direction
          !
          if (Sdy(k,i,j)<Sdy(k,i,j+1)) then
            dS(k,i,j) = dS(k,i,j) + dtau*min(sign(1.,a_phi(k,i,j))*ny(k,i,j)*Sdy(k,i,j),   &
                                             sign(1.,a_phi(k,i,j))*ny(k,i,j)*Sdy(k,i,j+1)  )
          else
            dS(k,i,j) = dS(k,i,j) + dtau*max(sign(1.,a_phi(k,i,j))*ny(k,i,j)*Sdy(k,i,j),   &
                                             sign(1.,a_phi(k,i,j))*ny(k,i,j)*Sdy(k,i,j+1)  )
          end if
          !
          ! z direction
          !
          if (Sdz(k,i,j)<Sdz(k+1,i,j)) then
            dS(k,i,j) = dS(k,i,j) + dtau*min(sign(1.,a_phi(k,i,j))*nz(k,i,j)*Sdz(k,i,j),   &
                                             sign(1.,a_phi(k,i,j))*nz(k,i,j)*Sdz(k+1,i,j)  )
          else
            dS(k,i,j) = dS(k,i,j) + dtau*max(sign(1.,a_phi(k,i,j))*nz(k,i,j)*Sdz(k,i,j),   &
                                             sign(1.,a_phi(k,i,j))*nz(k,i,j)*Sdz(k+1,i,j)  )
          end if
          end select
      end do
      !
      ! update scalar field
      !
      do l=1,ncells
          call unpackIndices(cell_indices, l, k, i, j)
          res      = res + dS(k,i,j)**2.
          S(k,i,j) = S(k,i,j) - dS(k,i,j)
      end do
      res = (res/ncells)**.5
      !
      ! boundary conditions
      !
      call cyclics(nzp,nxp,nyp,S,req)
      call cyclicc(nzp,nxp,nyp,S,req)

    end do
  end subroutine ls_spread2

  !
  !
  ! ls_fixmixed: sets the correct average values for in cut cells. after UCLALES initialization
  !           values in cut cells are either the one of the lower or the upper layer, NOT
  !           the appropriate cell average in a finite-volume sense.
  !           
  !
  subroutine ls_fixmixed(nzp,nxp,nyp,nscl,a,alpha,dzi_t,dzi_m)
    
    use grid, only : ls_G, ls_nG, write_anal,a_phi, ls_q0, ls_q1, deltax, deltay, a_xp, dn0, th00
    use util, only : atob
    use mpi_interface, only : cyclics, cyclicc
    integer, intent(in)  :: nzp,nxp,nyp,nscl
    real, intent(in)     :: alpha(nzp,nxp,nyp),dzi_t(nzp),dzi_m(nzp)
    real, intent(inout)  :: a(nzp,nxp,nyp,nscl)
    integer              :: k,i,j,l,req(16),n

    do n = 4, nscl
        do l = 1, ls_nG
            call unpackIndices(ls_G, l, k, i, j)
            a(k,i,j,n) = (1.0-alpha(k,i,j)) * ls_q0(k,i,j,n) &
                       +      alpha(k,i,j)  * ls_q1(k,i,j,n)
        end do

        call cyclics(nzp,nxp,nyp,a(:,:,:,n),req)
        call cyclicc(nzp,nxp,nyp,a(:,:,:,n),req)
    end do
  end subroutine ls_fixmixed

  real function NewVolumeFraction2D(g)

    real, dimension(0:3), intent(in) :: g
    real, dimension(0:3)             :: b,t
    real                             :: VF,g01, g32, t01, t32
    integer                          :: i

    select case(vfmode)
    case default
        NewVolumeFraction2D = 1.-VolumeFraction2D(g)

    case(1)
    ! NewVolumeFraction3D may supply us with a face that is not cut, i.e.
    ! all values in g have the same sign. For these cases, we test for sign
    ! changes. If the sign does not change, we simply return 1.0 or 0.0 .
    if (g(0).ge.0.0.and.g(1).ge.0.0  .and. &
        g(2).ge.0.0.and.g(3).ge.0.0) then
        NewVolumeFraction2D = 1.0
    else if (g(0).le.0.0.and.g(1).le.0.0  .and. &
             g(2).le.0.0.and.g(3).le.0.0) then
        NewVolumeFraction2D = 0.0
    else
        t01 = intermediateSlicePosition(g(0), g(1))
        t32 = intermediateSlicePosition(g(3), g(2))

        t(0) = 0.
        t(1) = max(0.0, min(t01, t32) - eps);
        t(2) = min(1.0, max(t01, t32) + eps);
        t(3) = 1.

        b(0) = intermediateSliceFraction(g(0), g(3))
        b(3) = intermediateSliceFraction(g(1), g(2))
        do i = 1, 2
            g01  = g(0) + t(i) * (g(1) - g(0))
            g32  = g(3) + t(i) * (g(2) - g(3))
            b(i) = intermediateSliceFraction(g32, g01)
        end do

        VF = 0.0
        do i = 0, 2
            VF = VF + (t(i+1) - t(i)) * (b(i+1) + b(i))
        end do
        VF = 0.5 * VF
        
        if ((VF > 1.0).or.(VF < 0.0)) then
          print *, '  WARNING: Volume fraction greater than 1 or less than 0: ', VF
        end if

        NewVolumeFraction2D = max(0.0, min(1.0, VF))
    end if
    end select
  end function NewVolumeFraction2D


  real function NewVolumeFraction3D(gb, gt)

    use mpi_interface, only : appl_abort

    real, dimension(0:3),intent(in) :: gb, gt
    real, dimension(0:5)            :: b, t, g
    real, dimension(0:2)            :: bb, tt, c
    real, dimension(0:2,0:2)        :: TMtrx, TInv
    real                            :: VF, detT
    integer                         :: j, jj, k, kk
    integer, dimension(0:5)         :: i_seq

    b(0)     = NewVolumeFraction2D(gb)
    t(0)     = 0.0
    i_seq(0) = 0

    b(5)     = NewVolumeFraction2D(gt)
    t(5)     = 1.0
    i_seq(5) = 5

    ! potential passage of zero level set at corners
    do j = 0, 3
      t(j+1) = intermediateSlicePosition(gb(j), gt(j))
    end do

    ! establish time sequence
    call sort(i_seq(1:4),t(1:4))
    do j = 1, 4
      i_seq(j) = i_seq(j) + 1
    end do

    ! within intervals use quadratic interpolation of area fractions
    VF = 0.0
    do j = 0, 4
      if (t(i_seq(j+1)) - t(i_seq(j)) > 3.0e-6) then
        ! get times and area fractions within sub-interval
        tt(0) = t(i_seq(j)) + eps
        tt(1) = 0.5 * (t(i_seq(j)) + t(i_seq(j+1)))
        tt(2) = t(i_seq(j+1)) - eps

        ! time coefficient matrix and its inverse
        do k = 0, 2
          do kk = 0, 2
            TMtrx(k, kk) = tt(k)**(2-kk)
          end do
        end do

        detT = Invert_3x3_Matrix(Tinv, TMtrx)  ! Tinv is the inverse matrix of T up to a factor of 1/detT

        if (abs(detT) < 1e-12) then ! integrate linear interpolation
            do jj = 0, 3
                G(jj) = Gb(jj) + tt(0) * (Gt(jj) - Gb(jj))
            end do
            bb(0)  = NewVolumeFraction2D(G)

            do jj = 0, 3
                G(jj) = Gb(jj) + tt(2) * (Gt(jj) - Gb(jj))
            end do
            bb(2)  = NewVolumeFraction2D(G)

            VF = VF + (tt(2)-tt(0)) * 0.5 * (bb(0) + bb(2))
        else ! integrate quadratic interpolation
            do k = 0, 2
              do jj = 0, 3
                G(jj) = Gb(jj) + tt(k) * (Gt(jj) - Gb(jj))
              end do
              bb(k)  = NewVolumeFraction2D(G)
            end do

            ! interpolation coefficients
            do k = 0, 2
              c(k) = 0.0
              do kk = 0, 2
                c(k) = c(k) + Tinv(k, kk) * bb(kk)
              end do
              c(k) = c(k) / detT
            end do

            VF = VF +       c(0) * (tt(2)*TMtrx(2, 0) - tt(0)*TMtrx(0,0)) / 3 &
                    + 0.5 * c(1) * (TMtrx(2, 0)       - TMtrx(0, 0)     )     &
                    +       c(2) * (tt(2)             - tt(0)           )
        end if
      end if
    end do         
  
    NewVolumeFraction3D = clippedVolumeFraction(VF)
  end function NewVolumeFraction3D



  real function intermediateSlicePosition(phi0, phi1)
    real, intent(in) :: phi0, phi1

    if (phi0 * phi1 < 0.0) then
        intermediateSlicePosition = - phi0 / (phi1 - phi0)
    else 
        intermediateSlicePosition = 1.0
    end if
  end function intermediateSlicePosition



  real function intermediateSliceFraction(phi0, phi1)
    real, intent(in) ::  phi0, phi1

    if (phi1 * phi0 < 0.0) then
        if (phi1 < 0.0) then
            intermediateSliceFraction = - phi0 / ( phi1 - phi0)
        else
            intermediateSliceFraction = - phi1 / ( phi0 - phi1)
        end if
    else
        intermediateSliceFraction = heaviside( max(phi0, phi1) )
    end if
  end function intermediateSliceFraction
  !
  ! VolumeFraction2D: integrates up the area fraction of negative 'g' in a rectangular
  !                   2D grid cell having the corner values g(i = 0,1,2,3) at the
  !                   bottom left, top left, top right, and bottom right, resp.
  !                   
  real function VolumeFraction2D(g)

    real, dimension(0:3), intent(in) :: g
    real, dimension(0:3)             :: b,t
    real                             :: VF,g01, g32, t01, t32, diff
    integer                          :: i
 
    t(0) = 0.
    t(1) = 0.
    t(2) = 1.
    t(3) = 1.

    if (g(0).ge.0.0.and.g(1).ge.0.0  .and. &
        g(2).ge.0.0.and.g(3).ge.0.0) then
        VolumeFraction2D = 0.0
    else if (g(0).le.0.0.and.g(1).le.0.0  .and. &
             g(2).le.0.0.and.g(3).le.0.0) then
        VolumeFraction2D = 1.0
    else
        if (g(1)*g(0) < 0.0) then
          t01 = -g(0)/(g(1)-g(0))
        else 
          t01 = 1.0
        end if
        if (g(2)*g(3) < 0.0) then
          t32 = -g(3)/(g(2)-g(3))
        else 
          t32 = 1.0
        end if

        t(1) = MAX(0.0,MIN(t01,t32)-eps);
        t(2) = MIN(1.0,MAX(t01,t32)+eps);
        !
        ! b(0)
        !
        if (g(3)*g(0) < 0.0) then
          if (g(3) < 0.0) then
            b(0) = -g(3) / (g(0) - g(3))
          else
            b(0) = -g(0) / (g(3) - g(0))
          end if
        else
          b(0) = 1. - heaviside( max(g(0),g(3)) )
        end if
        !
        ! b(3)
        !
        if (g(2)*g(1) < 0.0) then
          if (g(2) < 0.0) then
            b(3) = -g(2) / (g(1) - g(2))
          else
            b(3) = -g(1) / (g(2) - g(1))
          end if
        else
          b(3) = 1. - heaviside( max(g(1),g(2)) )
        end if
        !
        ! intermediate slices
        !
        do i=1,2
          g01  = g(0) + t(i) * (g(1) - g(0))
          g32  = g(3) + t(i) * (g(2) - g(3))
          if (g32*g01<0.0) then
            if (g32<0.0) then 
              b(i) = -g32/(g01-g32)
            else
              b(i) = -g01/(g32-g01)
            end if
          else
            b(i) = 1.0 - heaviside(max(g32,g01))
          end if 
        end do

        VF = 0.0
        do i=0,2
          VF = VF + (t(i+1) - t(i)) * (b(i+1) + b(i))
        end do
        VF = 0.5*VF
        
        if ((VF > 1.0).or.(VF < 0.0)) then
          print *, '  WARNING: Volume fraction greater than 1 or less than 0: ', VF
        end if

        VolumeFraction2D = MAX(0.0,MIN(1.0,VF))

    end if

  end function VolumeFraction2D
  !
  ! VolumeFraction3D: Integrates the volume fraction of negative 'g' in a cubical
  !                   3D grid cell having the corner values... 
  !
  real function VolumeFraction3D(gb,gt)

    use mpi_interface, only : appl_abort

    real, dimension(0:3),intent(in) :: gb, gt
    real, dimension(0:5)     :: t, g
    real, dimension(0:2)     :: bb, tt, c
    real, dimension(0:2,0:2) :: TMtrx, TInv
    real                     :: VF, detT
    integer                  :: j, jj, k, kk
    integer, dimension(0:5)  :: i_seq

    t(0)     = 0.0
    i_seq(0) = 0
    t(5)     = 1.0
    i_seq(5) = 5

    ! potential passage of zero level set at corners
    do j = 0,3
      if (Gb(j) * Gt(j) < 0.0) then
        t(j+1) = -Gb(j)/(Gt(j) -Gb(j))
      else
        t(j+1) = 1.0
      end if
    end do

    ! establish time sequence
    call sort(i_seq(1:4),t(1:4))
    do j=1,4
      i_seq(j) = i_seq(j) + 1
    end do

    ! within intervals use quadratic interpolation of area fractions
    VF = 0.0
    do j=0,4
      if (t(i_seq(j+1))-t(i_seq(j)) > 3.0e-6) then
        ! get times and area fractions within sub-interval
        tt(0) = t(i_seq(j)) + eps!tiny(1.0)
        tt(1) = 0.5*(t(i_seq(j))+t(i_seq(j+1)))
        tt(2) = t(i_seq(j+1)) - eps!tiny(1.0)

        ! time coefficient matrix and its inverse
        do k = 0,2
          do kk = 0,2
            TMtrx(k,kk) = tt(k)**(2-kk)
          end do
        end do

        detT = Invert_3x3_Matrix(Tinv,TMtrx)  ! Tinv is the inverse matrix of T up to a factor of 1/detT

        if (abs(detT) < 1e-12) then
            do jj = 0, 3
                G(jj) = Gb(jj) + tt(0) * (Gt(jj) - Gb(jj))
            end do
            bb(0)  = VolumeFraction2D(G)

            do jj = 0, 3
                G(jj) = Gb(jj) + tt(2) * (Gt(jj) - Gb(jj))
            end do
            bb(2)  = VolumeFraction2D(G)

            VF = VF + (tt(2)-tt(0)) * 0.5 * (bb(0) + bb(2))
        else
            do k = 0,2
              do jj = 0,3
                G(jj) = Gb(jj) + tt(k) * (Gt(jj) - Gb(jj))
              end do
              bb(k)  = VolumeFraction2D(G)
            end do

            ! interpolation coefficients
            do k = 0,2
              c(k) = 0.0
              do kk = 0,2
                c(k) = c(k) + Tinv(k,kk) * bb(kk)
              end do
              c(k) = c(k)/detT
            end do

            VF = VF + c(0)*(tt(2)*TMtrx(2,0) - tt(0)*TMtrx(0,0) ) / 3   &
                    + c(1)*(      TMtrx(2,0) -       TMtrx(0,0) ) * 0.5 &
                    + c(2)*(      tt(2)      -       tt(0)      )
        end if
      end if
    end do         
  
    VolumeFraction3D = clippedVolumeFraction(VF)

  end function VolumeFraction3D

  real function clippedVolumeFraction(VF)
  real, intent(in) :: VF
  if (VF < -3.0 * eps) then
    print *, '  WARNING: unphysical VF, setting to 0.0 ',VF
  end if

  if(VF > 1.0 + 3.0 * eps) then
    print *, '  WARNING: unphysical VF, setting to 1.0 ',VF
  end if

  clippedVolumeFraction = max(0.0, min(1.0, VF))
  end function clippedVolumeFraction

  !
  ! Entrainment driver routine
  !
  subroutine entrainment
    use grid, only : ls_up, ls_vp, ls_wp, a_phitilde, a_alphatilde, ls_cellkind,&
        ls_q1, ls_q0, a_sp, a_st, a_xp, dzi_t, dn0, nxp, nyp, nzp, nxyzp, nstep,&
        nscl, dt, ls_G, ls_nG, newphi, zm, levset, a_scr1, a_alpha, entmodl
    use util, only : atob, sclrset
    use advf, only : advtnd
    
    if (levset>=1.and.entmodl>=1) then
        
        ! set entrainment velocities
        select case(entmodl)

        case(1)
            ls_up(:,:,:) = 0.0
            ls_vp(:,:,:) = 0.0
            ls_wp(:,:,:) = 0.04
        
        case default
            stop '[entrainment] entrainment model not supported'
  
        end select
  
        ! move interface
        if (entmodl >= 1) then
            call atob(nxyzp, a_phitilde, a_scr1)
            call ls_fdupwind(a_scr1, a_phitilde, dt)
            call sclrset('cnst', nzp, nxp, nyp, a_phitilde)
            call advtnd(nzp, nxp, nyp, a_scr1, a_phitilde, a_st, dt)
  
            ! characterize_entrainment_cells
            call atob(nxyzp, a_alphatilde, a_alpha)
            call ls_characterize(a_scr1, a_phitilde)
            call ls_integrateAlpha(ls_G, ls_nG, a_phitilde, a_alphatilde, nzp, nxp, nyp)
  
            call entrain_scalars(nscl, nxp, nyp, nzp, dzi_t, dn0, a_alphatilde, ls_cellkind, ls_q1, a_xp, ls_q0)
        end if
    end if
  
  end subroutine entrainment
  !
  ! Compute effect of entrainment mixing based on the input:
  !     - a_alphatilde: post-entrainment volume fraction
  !     - ls_cellkind:  marker that marks cell that are cut during entrainment process.
  !                     Note: this is different from the pre/post advection alphas!
  !
  subroutine entrain_scalars(nscl, nxp, nyp, nzp, dzi_t, dn0, a_alphatilde, &
                             ls_cellkind, ls_q1, a_xp, ls_q0)
    use grid, only : a_phi, a_alpha, a_phitilde, ls_G, ls_nG
  
    integer, intent(in)                                 :: nscl, nxp, nyp, nzp
    real, dimension(nzp), intent(in)                    :: dn0, dzi_t
    real, dimension(nzp, nxp, nyp), intent(in)          :: a_alphatilde
    real, dimension(nzp, nxp, nyp, nscl), intent(in)    :: ls_q1
    real, dimension(nzp, nxp, nyp, nscl), intent(inout) :: a_xp, ls_q0
    integer, dimension(nzp, nxp, nyp), intent(in)       :: ls_cellkind
  
    real    :: num, denom, q0, scalar_total_in_column
    integer :: n, k, i, j, k0, k1, k00, k11
    logical :: column_is_cut

    do n = 4, nscl
    
    ! memorizing total mass of species n
    scalar_total_in_column = arraysum(nzp, dn0(:)*a_xp(:,3,3,n))
  
    do j = 3, nyp - 2
    do i = 3, nxp - 2
        ! find vertical range of cut cells during the mixing event
        column_is_cut = cutCellsInColumn(i, j, ls_cellkind, nzp, nxp, nyp, gkmin, gkmax, k0, k1)
        if (column_is_cut) then
  
           num   = 0.0
           denom = 0.0
           do k = k0, k1
               denom = denom + (1.0 - a_alphatilde(k,i,j)) * dn0(k) / dzi_t(k)
               num   = num   + (a_xp(k,i,j,n) - a_alphatilde(k,i,j) * ls_q1(k,i,j,n)) &
                             * dn0(k) / dzi_t(k)
           end do

           q0 = num / (denom + tiny(1.))

           ! Update q0 and q = alpha q1 + (1-alpha) q0
           do k = k0, k1
               ls_q0(k,i,j,n) = q0
               a_xp (k,i,j,n) = a_alphatilde(k,i,j)       * ls_q1(k,i,j,n) &
                              + (1.- a_alphatilde(k,i,j)) * q0
           end do
        else
           print *, '    Entrainment: There are uncut columns. Stopping.'
           print *, '      Column (i,j) = ',i,j
           print *, '      nxp-4        = ',nxp-4
           print *, '      nyp-4        = ',nyp-4
           print *, '          k         a_phi(k,i,j)            a_phi(k,i-1,j)          a_phi(k,i-1,j-1)            a_phi(k,i,j-1)            a_alpha(k,i,j)       a_alphatilde(k,i,j) ls_cellkind(k,i,j)'
           do k = 33, 31, -1
              print *, k, a_phi(k,i,j), a_phi(k,i-1,j), a_phi(k,i-1,j-1), a_phi(k,i,j-1), a_alpha(k,i,j), a_alphatilde(k,i,j), ls_cellkind(k,i,j)
           end do
           print *, '          k         a_phi(k,i,j)            a_phi(k,i-1,j)          a_phi(k,i-1,j-1)            a_phi(k,i,j-1)            a_alpha(k,i,j)       a_alphatilde(k,i,j) ls_cellkind(k,i,j)'
           do k = 33, 31, -1
              print *, k, a_phitilde(k,i,j), a_phitilde(k,i-1,j), a_phitilde(k,i-1,j-1), a_phitilde(k,i,j-1), a_alpha(k,i,j), a_alphatilde(k,i,j), ls_cellkind(k,i,j)
           end do
           stop('entrain: no cut cells in this column')
        end if
  
    end do
    end do
  
    end do
  
  end subroutine entrain_scalars
  
  
  real function arraysum(n, a)
  
    integer, intent(in) :: n
    real, dimension(n), intent(in)    :: a
    integer :: k
  
    arraysum = 0.0
    do k = 1, n
        arraysum = arraysum + a(k)
    end do
  
  end function arraysum
  
  
  real function qbar(alpha, q0, q1)
    real, intent(in) :: alpha, q0, q1
  
    qbar = alpha * q1 + (1.0 - alpha) * q0
  end function qbar
  !
  ! Inverts a 3x3 matrix T. The inverse is Tinv, the function value
  ! is the determinant of T.
  !
  real function Invert_3x3_Matrix(Tinv, T)
  ! returns the determinant detT of T[][] and computes the inverse of it
  ! up to a factor of 1/detT
  
  real, dimension(0:2,0:2),intent(in)  :: T
  real, dimension(0:2,0:2),intent(out) :: Tinv
  integer                              :: i, ii
  integer, dimension(0:2)              :: next, prev
  parameter (next = (/1,2,0/))
  parameter (prev = (/2,0,1/))
  !detT =   T[0][0] * (T[1][1]*T[2][2] - T[1][2]*T[2][1]) 
  !       - T[1][0] * (T[0][1]*T[2][2] - T[0][2]*T[2][1]) 
  !       + T[2][0] * (T[0][1]*T[1][2] - T[0][2]*T[1][1]);
  !
  ! determinant of T
  !
  Invert_3x3_Matrix =   T(0,0) * (T(1,1)*T(2,2) - T(1,2)*T(2,1)) &
                      - T(1,0) * (T(0,1)*T(2,2) - T(0,2)*T(2,1)) &
                      + T(2,0) * (T(0,1)*T(1,2) - T(0,2)*T(1,1))

          
  do i = 0,2
    do ii = 0,2
       Tinv(ii,i) =   T(next(i),next(ii)) * T(prev(i),prev(ii))  &
                    - T(next(i),prev(ii)) * T(prev(i),next(ii))
    end do
  end do 
  
  end function Invert_3x3_Matrix
  !
  ! Heaviside function
  !
  real function heaviside(x)

    real,intent(in) :: x

    if (x.ge.0.0) then
      heaviside = 1.0
    else
      heaviside = 0.0
    end if

  end function heaviside
  !
  ! 
  !
  subroutine sort(i_seq, t)
  ! an auxiliary function doing some sorting of indices with respect to
  ! increasing values of t[i]
  real,    dimension(0:3), intent(in)    :: t
  integer, dimension(0:3), intent(inout) :: i_seq
  integer                                :: i, j
  integer, dimension(0:3)                :: ihelp

  do i =0,3
    ihelp(i) = 0.0
    do j=0,3
      if ( (t(i)>t(j)).or.((t(i) == t(j)).and.(i > j) ) ) then
        ihelp(i) = ihelp(i) + 1
      end if
    end do
  end do

  do j=0,3
      i_seq(ihelp(j)) = j
  end do

  end subroutine sort
  !
  ! function b(chi): Buoyancy mixing function
  !
  real function bfunc(chi)

    real, intent(in) :: chi
    real :: deltas,D,chis,b1, chi_limited

    D      = 0.031
    chis   = 0.09
    b1     = 0.247
    deltas = chis/16.
    
    chi_limited = min(1.0,max(0.0,chi))
    bfunc = b1 * ( - D/chis * chi_limited  &
                   + ((1.+D)/(1.-chis)+D/chis)*deltas*log( exp((chi_limited-chis)/deltas) + 1.)  )
  end function bfunc
  !
  ! logical funtion ls_celltype: Determines wheather a cell is cut or not. It
  !     returns True in case it is cut and False in case it is not.
  !     
  !     The cut state is determined by summing up products of sign(phi) along
  !     the 12 edges of the cell. If the sum is 12, all phi values have the
  !     same sign and the cell is not cut. If the sum is a value less than 12,
  !     phi changes sign along one or more edges and the cell is cut. 
  !
  logical function ls_celltype(nzp,nxp,nyp,k,i,j,phi)
    integer,intent(in)                     :: nzp, nxp, nyp,k,i,j
    real,dimension(nzp,nxp,nyp),intent(in) :: phi
    ls_celltype = (12 .ne. &
              sign(1.,phi(k,  i,  j  )*phi(k-1,i,  j  )) + &
              sign(1.,(phi(k,  i,  j  )*phi(k,  i-1,j  ))) + &
              sign(1.,(phi(k,  i,  j  )*phi(k,  i,  j-1))) + &
              !
              sign(1.,(phi(k,  i-1,j-1)*phi(k-1,i-1,j-1))) + &
              sign(1.,(phi(k,  i-1,j-1)*phi(k,  i,  j-1))) + &
              sign(1.,(phi(k,  i-1,j-1)*phi(k,  i-1,j  ))) + &
              !
              sign(1.,(phi(k-1,i-1,j  )*phi(k-1,i-1,j-1))) + &
              sign(1.,(phi(k-1,i-1,j  )*phi(k-1,i,  j  ))) + &
              sign(1.,(phi(k-1,i-1,j  )*phi(k,i-1,  j  ))) + &
              !
              sign(1.,(phi(k-1,i,  j-1)*phi(k-1,i-1,j-1))) + &
              sign(1.,(phi(k-1,i,  j-1)*phi(k-1,i,  j  ))) + &
              sign(1.,(phi(k-1,i,  j-1)*phi(k  ,i  ,j-1))) )
  end function ls_celltype


  subroutine ls_nonconservative_fix
    use grid, only : nxp, nyp, nzp, ls_G, ls_nG, ls_q0, ls_q1, levset, a_alpha, &
        a_rp, a_tp, qt_change, th_change, deltax, deltay, dzi_t, dn0, filprf
    
    integer :: k, i, j, l, typ
    real :: qt_temp, th_temp, dxy, dV
    real :: q_bottom = 0.0
    real :: q_top = 1.0
    
    typ = 3
    dxy = deltax * deltay
    if (trim(filprf)=='wobble') then
        q_bottom = 0.0
        q_top = 1.0
    else
        q_bottom = 1e-3
        q_top = 0.0
    end if
    
    if (levset>=1) then
    select case(typ)
    
    case default
      stop("This nonconservative reconstruction case is not implemented. Stopping.")
    
    case(1) ! set qt = alpha
    do l=1,ls_nG
        call unpackIndices(ls_G, l, k, i, j)
        a_rp(k,i,j) = a_alpha(k,i,j)
    end do
    
    case(2) ! alpha-based qt
    do l=1,ls_nG
        call unpackIndices(ls_G, l, k, i, j)
        qt_temp = a_alpha(k,i,j) * q_top + (1.0 - a_alpha(k,i,j)) &
            * q_bottom
        a_rp(k,i,j) = qt_temp
    end do
    
    case(3) ! alpha-based qt and theta_l
    do l=1,ls_nG
        call unpackIndices(ls_G, l, k, i, j)
         
        dV = dxy / dzi_t(k)
    
        qt_temp     = a_alpha(k,i,j) * q_top + (1.0 - a_alpha(k,i,j)) * q_bottom
        qt_change   = qt_change + dn0(k) * (qt_temp - a_rp(k,i,j)) * dV
        a_rp(k,i,j) = qt_temp
    
        th_temp = a_alpha(k,i,j) * ls_q1(k,i,j,4) + (1.0 - a_alpha(k,i,j)) &
                * ls_q0(k,i,j,4)
        th_change = th_change + dn0(k) * (th_temp - a_tp(k,i,j)) * dV
        a_tp(k,i,j) = th_temp
    end do
    
    end select
    
    print *, "total water mass change: ", qt_change, " kg"
    print *, "total energy/cp change:  ", th_change, " kg K m^-3"
    
    end if
  end subroutine ls_nonconservative_fix
  
  subroutine packIndices(cell_set, counter, k, i, j)
  use grid, only : nxyzp
  integer, intent(in)    :: k, i, j
  integer, intent(inout) :: counter, cell_set(nxyzp, 3)
    counter = counter + 1
    cell_set(counter, 1) = k
    cell_set(counter, 2) = i
    cell_set(counter, 3) = j
  end subroutine packIndices
  
  subroutine unpackIndices(cell_set, counter, k, i, j)
  use grid, only : nxyzp
  integer, intent(in)  :: counter, cell_set(nxyzp, 3)
  integer, intent(out) :: k, i, j
      k = cell_set(counter, 1) 
      i = cell_set(counter, 2) 
      j = cell_set(counter, 3)
  end subroutine unpackIndices
  !
  ! cutCellsInColumn: In a column (i,j), finds indices k of cut cells and returns
  !     kstart (the lowest cut cell) and kstop (the highest cut cell). This info
  !     is based on 'cellkind', which assumes all cut cells marked 0 and uncut
  !     cells marked something else.
  !     
  !     The return value of the function is either 0 (in case there are no cut
  !     cells in the column or or 1 (if there are).
  !
  !     Parameters:
  !        i, j             : real, double
  !            x, y indices of column to be examined
  !
  !        cellkind         : real, double array
  !            3D marker array signifying the kind of cell
  !
  !        nzp, nxp, nyp    : integer
  !            shape of cellkind
  !
  !        gkmin, gkmax     : integer
  !            minimal and maximal 'k' index of the level set band and the range
  !            the search operates on
  ! 
  !     Returns:
  !        kstart, kstop    : integer
  !            indicex of lowest and highest cut cell, respectively
  !
  !        cutCellsInColumn : logical
  !            boolean indicating whether the column is intersected by the level
  !            set. Often .true.
  !
  logical function cutCellsInColumn(i, j, cellkind, nzp, nxp, nyp, gkmin, gkmax, kstart, kstop)
    integer, dimension(nzp, nxp, nyp), intent(in) :: cellkind
    integer, intent(in) :: i, j, nzp, nxp, nyp, gkmin, gkmax
    integer, intent(inout) :: kstart, kstop
    integer :: k
  
    cutCellsInColumn = .false.
    kstart = 0
    kstop  = 0
    do k = gkmin+1, gkmax-1
        if (cellkind(k,i,j) == 0) then
            if (cutCellsInColumn .eqv. .false.) then
                kstart = k
            end if
            cutCellsInColumn = .true.
            kstop = k
        end if
    end do
  end function cutCellsInColumn

end module lset
