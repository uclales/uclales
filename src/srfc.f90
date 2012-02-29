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
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module srfc

  integer :: isfctyp = 0
 !irina 
 ! real    :: sst = 292.
  real    :: zrough =  0.1
  real    :: ubmin  =  0.20
  real    :: dthcon = 100.0
  real    :: drtcon = 0.0
! linda, b
  logical :: larms=.true.
! linda, e

contains 
  !
  ! --------------------------------------------------------------------------
  ! SURFACE: Calcualtes surface fluxes using an algorithm chosen by ISFCLYR
  ! and fills the appropriate 2D arrays
  !
  !     default: specified thermo-fluxes (drtcon, dthcon)
  !     isfclyr=1: specified surface layer gradients (drtcon, dthcon)
  !     isfclyr=2: fixed lower boundary of water at certain sst
  !     isfclyr=3: bulk aerodynamic law with coefficeints (drtcon, dthcon)
  !irina
! linda, b
!  subroutine surface(sst)
  subroutine surface(sst,time_in)
! linda, e

    use defs, only: vonk, p00, rcp, g, cp, alvl, ep2
    use grid, only: nzp, nxp, nyp, a_up, a_vp, a_theta, vapor, zt, psrf,   &
         th00, umean, vmean, dn0, level, a_ustar, a_tstar, a_rstar,        &
         uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc,                           &
! linda, b
          shls,lhls
    use step, only : strtim
! linda, e
    use thrm, only: rslf
    use stat, only: sfc_stat, sflg
    use mpi_interface, only : nypg, nxpg, double_array_par_sum

    implicit none
 
 !irina
! linda, b
!    real, optional, intent (inout) :: sst
    real, optional, intent (inout) :: sst,time_in
! linda, e
 !   
    real :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp)       &
         ,wspd(nxp,nyp), bfct(nxp,nyp)

    integer :: i, j, iterate
    real    :: zs, bflx0,bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum
    real (kind=8) :: bfl(2), bfg(2)
! linda,b
    real ::times,tfrac,tlag
    integer::tcnt,tcnt2
!linda,e


    select case(isfctyp)

       !
       ! use prescribed surface gradients dthcon, drton from NAMELIST
       ! use then similarity theory to compute the fluxes 
       !
    case(1)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j)=dthcon
             drdz(i,j)=drtcon
          end do
       end do
       zs = zrough
       call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar     &
            ,a_rstar)
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

       !
       ! use prescribed SST and assume qsurf=qsat (i.e. ocean) to compute
       ! gradients. Then use similarity theory to predict the fluxes. 
       !
    case(2)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       usum = 0.
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
             drdz(i,j) = vapor(2,i,j) - rslf(psrf,sst)
             bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
             usum = usum + a_ustar(i,j)
          end do
       end do
       usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
       zs = zrough
       if (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
       call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar     &
            ,a_rstar)
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
       !
       ! get fluxes from bulk formulae with coefficients given by dthcon and
       ! drtcon (wq=Ch*u*dth, Garrat p.55) and using prescribed sst 
       !  and qsurf=qsat (ocean); note that here zrough is not the roughness
       ! length but the drag coefficient
       !
   case(3)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
             drdz(i,j) = vapor(2,i,j) - rslf(psrf,sst)
   !cgils
            ! drdz(i,j) = vapor(2,i,j) - 0.98*rslf(psrf,sst)
             if (ubmin > 0.) then
                a_ustar(i,j) = sqrt(zrough)* wspd(i,j)
             else
                a_ustar(i,j) = abs(ubmin)
             end if
   !cgils  : drag coeff C = 0.0012*6.75 m/s = 0.0081 m/s
             !a_tstar(i,j) =  0.0081*dtdz(i,j)/a_ustar(i,j)
             !a_rstar(i,j) =  0.0081*drdz(i,j)/a_ustar(i,j)
             a_tstar(i,j) =  dthcon * wspd(i,j)*dtdz(i,j)/a_ustar(i,j)
             a_rstar(i,j) =  drtcon * wspd(i,j)*drdz(i,j)/a_ustar(i,j)
             bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
          end do
       end do
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
       !
       ! fix surface temperature to yield a constant surface buoyancy flux
       ! dthcon
       !
   case(4)

       Vzt   = 10.* (log(zt(2)/zrough)/log(10./zrough))       
       Vbulk = Vzt * (vonk/log(zt(2)/zrough))**2

       bfl(:) = 0.
       do j=3,nyp-2
          do i=3,nxp-2
             bfl(1) = bfl(1)+a_theta(2,i,j)
             bfl(2) = bfl(2)+vapor(2,i,j)
          end do
       end do

       call double_array_par_sum(bfl,bfg,2)

       bfg(2) = bfg(2)/real((nxpg-4)*(nypg-4))
       bfg(1) = bfg(1)/real((nxpg-4)*(nypg-4))

       do iterate=1,5
          bflx  = ((sst -bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst) -bfg(2))) &
               * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
          sst1 = sst + 0.1
          bflx1 = ((sst1-bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst1)-bfg(2))) &
               * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
          sst  = sst + 0.1* (dthcon - bflx) / (bflx1-bflx)
       end do

       do j=3,nyp-2
          do i=3,nxp-2
             wt_sfc(i,j) = Vbulk * (sst -a_theta(2,i,j))
             wq_sfc(i,j) = Vbulk * (rslf(psrf,sst) - vapor(2,i,j))
             wspd(i,j)    = max(0.1,                                    &
                  sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
             bflx         = wt_sfc(i,j)*g/bfg(1) + g*ep2*wq_sfc(i,j)
             a_ustar(i,j) = diag_ustar(zt(2),zrough,bflx,wspd(i,j))
             uw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                  *(a_up(2,i,j)+umean)/wspd(i,j)
             vw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                  *(a_vp(2,i,j)+vmean)/wspd(i,j)
             ww_sfc(i,j)  = 0.
             a_rstar(i,j) = wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = wt_sfc(i,j)/a_ustar(i,j)
          end do
       end do   
       !
       ! fix thermodynamic fluxes at surface given values in energetic 
       ! units and calculate  momentum fluxes from similarity theory
       !
    case default
       ffact = 1.

! linda, b
       if (larms) then ! linda
          tlag=3600.0
          times=(time_in-strtim)*24.0*3600.0
          tfrac=(times/tlag)-int(times/tlag)
          tcnt=int(times/tlag)+2
          tcnt2=tcnt-1
          if (tcnt.gt.24) then
             tcnt=tcnt-24
             if (tcnt.ne.1) tcnt2=tcnt-1
          end if
          wt_sfc(1,1)=(shls(tcnt2)+(shls(tcnt)-shls(tcnt2))*tfrac)/(0.5*(dn0(1)+dn0(2))*cp)
          wq_sfc(1,1)=(lhls(tcnt2)+(lhls(tcnt)-lhls(tcnt2))*tfrac)/(0.5*(dn0(1)+dn0(2))*alvl)
       else
          wt_sfc(1,1)  = ffact* dthcon/(0.5*(dn0(1)+dn0(2))*cp)
          wq_sfc(1,1)  = ffact* drtcon/(0.5*(dn0(1)+dn0(2))*alvl)
       endif
! linda, e

!       wt_sfc(1,1)  = ffact* dthcon/(0.5*(dn0(1)+dn0(2))*cp)
!       wq_sfc(1,1)  = ffact* drtcon/(0.5*(dn0(1)+dn0(2))*alvl)

       if (zrough <= 0.) then
          usum = 0.
          do j=3,nyp-2
             do i=3,nxp-2
                usum = usum + a_ustar(i,j)
             end do
          end do
          usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
          zs = max(0.0001,(0.016/g)*usum**2) !Charnock for flow over sea
       else
          zs = zrough
       end if

       do j=3,nyp-2
          do i=3,nxp-2
             wt_sfc(i,j)=wt_sfc(1,1)
             wq_sfc(i,j)=wq_sfc(1,1)

             wspd(i,j)    = max(0.1,                                    &
                  sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
             if (ubmin > 0.) then
                bflx = g*wt_sfc(1,1)/th00
                if (level >= 2) bflx = bflx + g*ep2*wq_sfc(i,j)
                a_ustar(i,j) = diag_ustar(zt(2),zs,bflx,wspd(i,j))
             else
                a_ustar(i,j) = abs(ubmin)
             end if

             ffact = a_ustar(i,j)*a_ustar(i,j)/wspd(i,j)
             uw_sfc(i,j)  = -ffact*(a_up(2,i,j)+umean)
             vw_sfc(i,j)  = -ffact*(a_vp(2,i,j)+vmean)
             ww_sfc(i,j)  = 0.
             a_rstar(i,j) = -wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = -wt_sfc(i,j)/a_ustar(i,j)
          end do
       end do

    end select

    if (sflg) call sfc_stat(nxp,nyp,wt_sfc,wq_sfc,a_ustar,sst)

    return
  end subroutine surface
  !
  ! -------------------------------------------------------------------
  ! GET_SWNDS: returns surface winds valid at cell centers
  !
  subroutine get_swnds(n1,n2,n3,usfc,vsfc,wspd,up,vp,umean,vmean)

    implicit none
    
    integer, intent (in) :: n1, n2, n3
    real, intent (in)    :: up(n1,n2,n3), vp(n1,n2,n3), umean, vmean
    real, intent (out)   :: usfc(n2,n3), vsfc(n2,n3), wspd(n2,n3)

    integer :: i, j, ii, jj

    do j=3,n3-2
       jj = j-1
       do i=3,n2-2
          ii = i-1
          usfc(i,j) = (up(2,i,j)+up(2,ii,j))*0.5+umean
          vsfc(i,j) = (vp(2,i,j)+vp(2,i,jj))*0.5+vmean
          wspd(i,j) = max(abs(ubmin),sqrt(usfc(i,j)**2+vsfc(i,j)**2))
       enddo
    enddo

  end subroutine get_swnds
  !
  ! ----------------------------------------------------------------------
  ! FUNCTION GET_USTAR:  returns value of ustar using the below 
  ! similarity functions and a specified buoyancy flux (bflx) given in
  ! kinematic units
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  real function diag_ustar(z,z0,bflx,wnd)

    use defs, only: vonk

    implicit none

    real, parameter      :: am   =  4.8   !   "          "         "
    real, parameter      :: bm   = 19.3   !   "          "         "
    real, parameter      :: eps  = 1.e-10 ! non-zero, small number

    real, intent (in)    :: z             ! height where u locates
    real, intent (in)    :: z0            ! momentum roughness height
    real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real, intent (in)    :: wnd           ! wind speed at z

    integer :: iterate
    real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

    lnz   = log(z/z0) 
    klnz  = vonk/lnz              
    c1    = 3.14159/2. - 3.*log(2.)

    ustar =  wnd*klnz
    if (bflx /= 0.0) then 
       do iterate=1,4
          lmo   = -(ustar**3)/(bflx*vonk + eps)
          zeta  = z/lmo
          if (zeta > 0.) then
             ustar =  vonk*wnd  /(lnz + am*zeta)
          else
             x     = sqrt( sqrt( 1.0 - bm*zeta ) )
             psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
             ustar = wnd*vonk/(lnz - psi1)
          end if
       end do
    end if

    diag_ustar = ustar

    return
  end function diag_ustar
  !
  ! ----------------------------------------------------------------------
  ! Subroutine srfcscls:  returns scale values based on Businger/Dye
  ! similarity functions.
  !
  ! phi_h (zeta > 0) =  Pr * (1 + ah * zeta)
  ! phi_h (zeta < 0) =  Pr * (1 - bh * zeta)^(-1/2)
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface 
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  subroutine srfcscls(n2,n3,z,z0,th00,u,dth,drt,ustar,tstar,rstar)

    use defs, only : vonk, g, ep2

    implicit none

    real, parameter     :: ah   =  7.8   ! stability function parameter
    real, parameter     :: bh   = 12.0   !   "          "         "
    real, parameter     :: am   =  4.8   !   "          "         "
    real, parameter     :: bm   = 19.3   !   "          "         "
    real, parameter     :: pr   = 0.74   ! prandlt number
    real, parameter     :: eps  = 1.e-10 ! non-zero, small number

    integer, intent(in) :: n2,n3         ! span of indicies covering plane
    real, intent(in)    :: z             ! height where u & T locate
    real, intent(in)    :: z0            ! momentum roughness height
    real, intent(in)    :: th00          ! reference temperature
    real, intent(in)    :: u(n2,n3)      ! velocities at z
    real, intent(in)    :: dth(n2,n3)    ! theta (th(z) - th(z0))
    real, intent(in)    :: drt(n2,n3)    ! qt(z) - qt(z0)
    real, intent(inout) :: ustar(n2,n3)  ! scale velocity
    real, intent(inout) :: tstar(n2,n3)  ! scale temperature
    real, intent(inout) :: rstar(n2,n3)  ! scale value of qt

    logical, save :: first_call=.True.
    integer :: i,j,iterate
    real    :: lnz, klnz, betg, cnst1, cnst2
    real    :: x, y, psi1, psi2, zeta, lmo, dtv

    lnz   = log(z/z0) 
    klnz  = vonk/lnz              
    betg  = th00/g
    cnst2 = -log(2.)
    cnst1 = 3.14159/2. + 3.*cnst2

    do j=3,n3-2
       do i=3,n2-2
          dtv = dth(i,j) + ep2*th00*drt(i,j)
          !
          ! stable case  ! Stable case is not tested!!
          !
          if (dtv > 0.) then
             x     = pr*(betg*u(i,j)**2)/dtv
             y     = (am*z - 0.5*x)/lnz
             x     = (x*ah*z - am**2*z*z)/(lnz**2)
             lmo   = -y + sqrt(x+y**2)
             zeta  = z/lmo
             ustar(i,j) =  vonk*u(i,j)  /(lnz + am*zeta)
             tstar(i,j) = (vonk*dtv/(lnz + ah*zeta))/pr
             !
             ! Neutral case
             ! 
          elseif (dtv == 0.) then
             ustar =  vonk*u(i,j)  /lnz
             tstar =  vonk*dtv/(pr*lnz)
             !
             ! ustable case, start iterations from values at previous tstep, 
             ! unless the sign has changed or if it is the first call, then 
             ! use neutral values.
             !
          else
             if (first_call .or. tstar(i,j)*dtv <= 0.) then
                ustar(i,j) = u(i,j)*klnz
                tstar(i,j) = (dtv*klnz/pr)
             end if

             do iterate = 1,3
                lmo   = betg*ustar(i,j)**2/(vonk*tstar(i,j))
                zeta  = z/lmo
                x     = sqrt( sqrt( 1.0 - bm*zeta ) )
                psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + cnst1
                y     = sqrt(1.0 - bh*zeta)
                psi2  = 2.*log(1.0 + y) +2.*cnst2
                ustar(i,j) = u(i,j)*vonk/(lnz - psi1)
                tstar(i,j) = (dtv*vonk/pr)/(lnz - psi2)
             end do
          end if

          rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
          tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)
       end do
    end do

    first_call = .False.

    return
  end subroutine srfcscls
  !
  ! ----------------------------------------------------------------------
  ! subroutine: sfcflxs:  this routine returns the surface fluxes based
  ! on manton-cotton algebraic surface layer equations. 
  !
  subroutine sfcflxs(n2,n3,vk,ubar,u,v,xx,us,ts,rs,uw,vw,tw,rw,ww)
    implicit none
    real, parameter      :: cc=4.7,eps=1.e-20

    integer, intent(in)  :: n2,n3
    real, intent(in)     :: ubar(n2,n3),u(n2,n3),v(n2,n3),xx(n2,n3),vk
    real, intent(in)     :: us(n2,n3),ts(n2,n3),rs(n2,n3)
    real, intent(out)    :: uw(n2,n3),vw(n2,n3),tw(n2,n3),rw(n2,n3),ww(n2,n3)

    real    :: x(n2,n3),y(n2,n3)
    integer i,j

    do j=3,n3-2
       do i=3,n2-2

          uw(i,j)=-(u(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          vw(i,j)=-(v(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          tw(i,j)=-ts(i,j)*us(i,j)
          rw(i,j)=-rs(i,j)*us(i,j)

          x(i,j) = xx(i,j)*vk*ts(i,j)*(ubar(i,j)/us(i,j))**2
          x(i,j) = x(i,j)*sqrt(sqrt(1.-15.*min(0.,x(i,j)))) &
               /(1.0+cc*max(0.,x(i,j)))
          y(i,j) =sqrt((1.-2.86*x(i,j))/(1.+x(i,j)* &
               (-5.39+x(i,j)*6.998 )))
          ww(i,j)=(0.27*max(6.25*(1.-x(i,j))*y(i,j),eps)-&
               1.18*x(i,j)*y(i,j))*us(i,j)**2
       enddo
    enddo
    return
  end subroutine sfcflxs

end module srfc



