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

use lsmdata

  integer :: isfctyp = 0
  real    :: zrough =  0.01
  real    :: ubmin  =  0.20
  real    :: dthcon = 100.0
  real    :: drtcon = 0.0

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
  subroutine surface(sst)

    use defs, only: vonk, p00, rcp, g, cp, alvl, ep2
    use grid, only: nzp, nxp, nyp, a_up, a_vp, a_theta, vapor, zt, psrf,   &
         th00, umean, vmean, dn0, level, a_ustar, a_tstar, a_rstar,        &
         uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc, nstep
    use thrm, only: rslf
    use stat, only: sfc_stat, sflg
    use mpi_interface, only : nypg, nxpg, double_array_par_sum

    implicit none
 
    real, optional, intent (inout) :: sst
    integer		:: i, j, iterate
    real		:: zs, bflx0,bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum
    real (kind=8)	:: bfl(2), bfg(2)

    real :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp) &
            ,wspd(nxp,nyp), bfct(nxp,nyp), ustar(nxp,nyp)

    select case(isfctyp)

  !
  ! ----------------------------------------------------------------------
  ! Set surface gradients
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
            ,a_rstar,obl)
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

  !
  ! ----------------------------------------------------------------------
  ! Get fluxes from profiles
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
            ,a_rstar,obl)
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

  !
  ! ----------------------------------------------------------------------
  ! Get fluxes from bulk formulae with coefficients given by
  ! dthcon and drtcon
  !
   case(3)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
             drdz(i,j) = vapor(2,i,j) - rslf(psrf,sst)
   	    !cgils:
            ! drdz(i,j) = vapor(2,i,j) - 0.98*rslf(psrf,sst)
             if (ubmin > 0.) then
                a_ustar(i,j) = sqrt(zrough)* wspd(i,j)
             else
                a_ustar(i,j) = abs(ubmin)
             end if
   	     !cgils: drag coeff C = 0.0012*6.75 m/s = 0.0081 m/s
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
  ! ----------------------------------------------------------------------
  ! Fix surface temperature to yield a constant surface buoyancy flux
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
  ! ----------------------------------------------------------------------
  ! Malte: Get surface fluxes using a land surface model
  !
   case(5)

       !Initialize Land Surface 
       if (init_lsm) then
          print*,"Start Case(5) - initialize LSM"
          call initlsm
          print*,"Land surface initialized..."
          init_lsm = .false.
       end if

       if (local) then
       u0bar(:,:,:)  = a_up(:,:,:)
       v0bar(:,:,:)  = a_vp(:,:,:)
       thetabar(:,:) = a_theta(2,:,:)
       vaporbar(:,:) = vapor(2,:,:)
       tskinbar(:,:) = tskin(:,:)
       qskinbar(:,:) = qskin(:,:)
       else
       u0av          = sum(a_up(2,3:nxp-2,3:nxp-2))/(nxp-4)/(nyp-4) + umean
       v0av          = sum(a_vp(2,3:nxp-2,3:nxp-2))/(nxp-4)/(nyp-4) + vmean
       thetabar(:,:) = sum(a_theta(2,3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)
       vaporbar(:,:) = sum(vapor(2,3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)
       tskinbar(:,:) = sum(tskin(3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)
       qskinbar(:,:) = sum(qskin(3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)
       end if

       !Calculate surface wind for flux calculation
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,u0bar,v0bar,umean,vmean)

       ! a) Calculate Monin Obuhkov Length iteratively
       !call getobl

       ! b) Calculate Monin Obuhkov Length from surface scalars
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j) = thetabar(i,j) - tskinbar(i,j)
             drdz(i,j) = vaporbar(i,j) - qskinbar(i,j)
          end do
       end do
       tskinavg = sum(tskin(3:(nxp-2),3:(nyp-2)))/(nxp-4)/(nyp-4)
       call srfcscls(nxp,nyp,zt(2),zrough,tskinavg,wspd,dtdz,drdz,a_ustar,a_tstar,a_rstar,obl)

       !Calculate the drag coefficients and aerodynamic resistance
       do j=3,nyp-2
          do i=3,nxp-2
             cm(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
			obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
             cs(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
			obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zt(2) / &
			z0h(i,j)) - psih(zt(2) / obl(i,j)) + psih(z0h(i,j) / &
			obl(i,j)))
             ra(i,j) =  1. / (cs(i,j)* wspd(i,j))
          end do
       end do

       !Get skin temperature and humidity from land surface model
       call lsm

       !Calculate the surface fluxes with bulk law (Fairall, 2003)
       !Fluxes in kinematic form with units [K m/s] and [kg/kg m/s]
       do j=3, nyp-2
          do i=3, nxp-2

             wt_sfc(i,j) = - (thetabar(i,j) - tskinbar(i,j)) / ra(i,j) 
             wq_sfc(i,j) = - (vaporbar(i,j) - qskinbar(i,j)) / ra(i,j)

             uw_sfc(i,j)  = - a_ustar(i,j)*a_ustar(i,j)                  &
                              *(a_up(2,i,j)+umean)/wspd(i,j)
             vw_sfc(i,j)  = - a_ustar(i,j)*a_ustar(i,j)                  &
                              *(a_vp(2,i,j)+vmean)/wspd(i,j)
             ww_sfc(i,j)  = 0.

             !uw_sfc(i,j)  = - 0.04
             !vw_sfc(i,j)  = - 0.005

             a_rstar(i,j) = - wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = - wt_sfc(i,j)/a_ustar(i,j)

          if  (wt_sfc(i,j) .lt. 0) then
             wt_sfc(i,j) = 0.
          end if

          if  (wq_sfc(i,j) .lt. 0) then
             wq_sfc(i,j) = 0.
          end if

          end do
       end do

  !
  ! ----------------------------------------------------------------------
  ! fix thermodynamic fluxes at surface given values in energetic 
  ! units and calculate momentum fluxes from winds
  !
   case default

       ffact = 1.
       wt_sfc(1,1)  = ffact* dthcon/(0.5*(dn0(1)+dn0(2))*cp)
       wq_sfc(1,1)  = ffact* drtcon/(0.5*(dn0(1)+dn0(2))*alvl)

       if (zrough <= 0.) then
          usum = 0.
          do j=3,nyp-2
             do i=3,nxp-2
                usum = usum + a_ustar(i,j)
             end do
          end do
          usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
          zs = max(0.0001,(0.016/g)*usum**2)
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
             a_rstar(i,j) = wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = wt_sfc(i,j)/a_ustar(i,j)
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

    real, parameter      :: am   =  4.7   !   "          "         "
    real, parameter      :: bm   = 16.0   !   "          "         "
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
  subroutine srfcscls(n2,n3,z,z0,th00,u,dth,drt,ustar,tstar,rstar,obl)

    use defs, only : vonk, g, ep2
    use grid ,only: nstep

    implicit none

    real, parameter     :: ah   =  4.7   ! stability function parameter
    real, parameter     :: bh   = 16.0   !   "          "         "
    real, parameter     :: am   =  4.7   !   "          "         "
    real, parameter     :: bm   = 16.0   !   "          "         "
    real, parameter     :: pr   =  0.74  ! prandlt number
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
    real, intent(inout) :: obl(n2,n3)    ! Obukhov Length

    logical, save :: first_call=.True.
    integer :: i,j,iterate
    real    :: lnz, klnz, betg
    real    :: x, y, zeta, lmo, dtv

    lnz   = log(z/z0) 
    klnz  = vonk/lnz              
    betg  = th00/g

    do j=3,n3-2
       do i=3,n2-2
          dtv = dth(i,j) + ep2*th00*drt(i,j)

          !
          ! STABLE CASE
          !
          if (dtv > 0.) then
             x     = (betg*u(i,j)**2)/dtv
             y     = (am - 0.5*x)/lnz
             x     = (x*ah - am**2)/(lnz**2)
             lmo   = -y + sqrt(x+y**2)
             zeta  = z/lmo
             ustar(i,j) =  vonk*u(i,j)/(lnz + am*zeta)
             tstar(i,j) = (vonk*dtv/(lnz + ah*zeta))/pr
          !print*,"STABLE ATMOSPHERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

          !
          ! NEUTRAL CASE
          ! 
          elseif (dtv == 0.) then
             ustar =  vonk*u(i,j)/lnz
             tstar =  vonk*dtv/(pr*lnz)
             lmo = -1.e10
          !print*,"NEUTRAL ATMOSPHERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

          !
          ! UNSTABLE CASE
          !
          ! start iterations from values at previous tstep, 
          ! unless the sign has changed or if it is the first call, then 
          ! use neutral values.
          !
          else
             if (first_call .or. tstar(i,j)*dtv <= 0.) then
             ustar(i,j) =  vonk*u(i,j)/lnz
             tstar(i,j) =  vonk*dtv/(pr*lnz)
             lmo = -1.e10
             end if

             do iterate = 1,3
                lmo   = betg*ustar(i,j)**2/(vonk*tstar(i,j))
                if (lmo .gt. -1) then
                lmo = -1.
                end if
                zeta  = z/lmo
                ustar(i,j) = u(i,j)*vonk/(lnz - psim(zeta))
                tstar(i,j) = (dtv*vonk/pr)/(lnz - psih(zeta))
             end do
          end if

          obl(i,j) = lmo
          rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
          tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)
          
       end do
    end do

       if (nstep == 1 .and. .false.) then
       print*,"********************************"
       print*,"*********SURFACE SCALARS********"
       print*,"********************************"
       print*,"obl (min)",minval(obl(3:n2-2,3:n3-2))
       print*,"obl (max)",maxval(obl(3:n2-2,3:n3-2))
       print*,"ustar (min)",minval(ustar(3:n2-2,3:n3-2))
       print*,"ustar (max)",maxval(ustar(3:n2-2,3:n3-2))
       print*,"tstar (min)",minval(tstar(3:n2-2,3:n3-2))
       print*,"tstar (max)",maxval(tstar(3:n2-2,3:n3-2))
       print*,"********************************"
       print*,"********************************"
       print*,"********************************"
       end if

    first_call = .False.

    return
  end subroutine srfcscls

  !
  ! ----------------------------------------------------------------------
  ! subroutine: sfcflxs:  this routine returns the surface fluxes based
  ! on manton-cotton algebraic surface layer equations. 
  !
  subroutine sfcflxs(n2,n3,vk,ubar,u,v,xx,us,ts,rstar,uw,vw,tw,rw,ww)
    implicit none
    real, parameter      :: cc=4.7,eps=1.e-20

    integer, intent(in)  :: n2,n3
    real, intent(in)     :: ubar(n2,n3),u(n2,n3),v(n2,n3),xx(n2,n3),vk
    real, intent(in)     :: us(n2,n3),ts(n2,n3),rstar(n2,n3)
    real, intent(out)    :: uw(n2,n3),vw(n2,n3),tw(n2,n3),rw(n2,n3),ww(n2,n3)

    real    :: x(n2,n3),y(n2,n3)
    integer i,j

    do j=3,n3-2
       do i=3,n2-2

          uw(i,j)=-(u(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          vw(i,j)=-(v(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          tw(i,j)=-ts(i,j)*us(i,j)
          rw(i,j)=-rstar(i,j)*us(i,j)

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

  !
  ! ----------------------------------------------------------------------
  ! subroutine: getobl - Calculates the Obuhkov length iteratively.
  !
  subroutine getobl
    use defs, only: g,ep2
    use grid, only: nzp,nxp,nyp,a_up,a_vp,a_theta,umean,vmean,vapor,zt,u0,v0
    implicit none

    integer        :: i,j,iter
    real           :: upcu, vpcv
    real           :: thv,thvs, thvsl
    real           :: L, Lstart, Lend, Lold
    real           :: Rib, fx, fxdif, wspd2
    real           :: thetavbar, tvskinbar
    !real          :: tskinav,thvskinbar,qskinav,thl0av,qt0av,u0av2,v0av2
     
    if (local) then

      do j=3,nyp-2
        do i=3,nxp-2
          thetavbar   = thetabar(i,j) * (1. + ep2 * vaporbar(i,j))
          tvskinbar   = tskinbar(i,j) * (1. + ep2 * qskinbar(i,j))

          upcu    = 0.5 * (u0bar(2,i,j) + u0bar(2,i+1,j)) + umean
          vpcv    = 0.5 * (v0bar(2,i,j) + v0bar(2,i,j+1)) + vmean
          wspd2   = max(abs(ubmin), upcu**2. + vpcv**2.)

          Rib     = g / thetavbar * zt(2) * (thetavbar - tvskinbar) / wspd2

          iter = 0
          L = obl(i,j)

          if (Rib * L < 0. .or. abs(L) == 1e5) then
             if(Rib > 0) L = 0.01
             if(Rib < 0) L = -0.01
          end if

          do while (.true.)
             iter    = iter + 1
             Lold    = L
             fx      = Rib - zt(2) / L * (log(zt(2) / z0h(i,j)) - psih(zt(2)  &
                       / L) + psih(z0h(i,j) / L)) / (log(zt(2) / z0m(i,j)) -  &
                       psim(zt(2) / L) + psim(z0m(i,j) / L)) ** 2.

             Lstart  = L - 0.001*L
             Lend    = L + 0.001*L

             fxdif   = ( (- zt(2) / Lstart * (log(zt(2) / z0h(i,j)) - psih(zt(2) / Lstart) + psih(z0h(i,j) / Lstart)) / (log(zt(2) / z0m(i,j)) - psim(zt(2) / Lstart) + psim(z0m(i,j) / Lstart)) ** 2.) - (-zt(2) / Lend * (log(zt(2) / z0h(i,j)) - psih(zt(2) / Lend) + psih(z0h(i,j) / Lend)) / (log(zt(2) / z0m(i,j)) - psim(zt(2) / Lend) + psim(z0m(i,j) / Lend)) ** 2.) ) / (Lstart - Lend)

             L       = L - fx / fxdif

             if(Rib * L < 0. .or. abs(L) == 1e5) then
               if(Rib > 0) L = 0.01
               if(Rib < 0) L = -0.01
             end if

             if(abs(L - Lold) < 0.0001) exit
             if(iter > 1000) stop 'Obukhov length calculation does not converge!'
           end do

           obl(i,j) = L

        end do
      end do
    end if

    ! CvH also do a global evaluation if local = .true.
    ! to get an appropriate local mean ??????????????????????????
    if (.not. local) then
       thetavbar   = thetabar(3,3) * (1. + ep2 * vaporbar(3,3))
       tvskinbar   = tskinbar(3,3) * (1. + ep2 * qskinbar(3,3))

       wspd2   = max(abs(ubmin), u0av**2. + v0av**2.)
       print*,wspd2,thetavbar,tvskinbar

       Rib     = g / thetavbar * zt(2) * (thetavbar - tvskinbar) / wspd2

       iter = 0
       L = oblav

       if(Rib * L < 0. .or. abs(L) == 1e5) then
         if(Rib > 0) L = 0.01
         if(Rib < 0) L = -0.01
       end if

       do while (.true.)
         iter    = iter + 1
         Lold    = L
         fx      = Rib - zt(2) / L * (log(zt(2) / z0hav) - psih(zt(2) / L) + psih(z0hav / L)) / (log(zt(2) / z0mav) - psim(zt(2) / L) + psim(z0mav / L)) ** 2.
         Lstart  = L - 0.001*L
         Lend    = L + 0.001*L
         fxdif   = ( (- zt(2) / Lstart * (log(zt(2) / z0hav) - psih(zt(2) / Lstart) + psih(z0hav / Lstart)) / (log(zt(2) / z0mav) - psim(zt(2) / Lstart) + psim(z0mav / Lstart)) ** 2.) - (-zt(2) / Lend * (log(zt(2) / z0hav) - psih(zt(2) / Lend) + psih(z0hav / Lend)) / (log(zt(2) / z0mav) - psim(zt(2) / Lend) + psim(z0mav / Lend)) ** 2.) ) / (Lstart - Lend)
         L       = L - fx / fxdif
         if(Rib * L < 0. .or. abs(L) == 1e5) then
            if(Rib > 0) L = 0.01
            if(Rib < 0) L = -0.01
         end if
         if(abs(L - Lold) < 0.0001) exit
         if(iter > 1000) stop 'Obukhov length calculation does not converge!'
       end do

       obl(:,:) = L
    end if

    oblav = L

    !print*,"**********************************************"
    !print*,"wspd2:",wspd2
    !print*,"Rib:",Rib
    !print*,"obl",minval(obl(3:nxp-2,3:nxp-2))
    !print*,"**********************************************"

    return
  end subroutine getobl

  !
  ! ----------------------------------------------------------
  ! Malte: Integrated stability function for momentum
  !
  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 15. * zeta) ** (0.25)
      !psim = 3.14159265/2. - 2. *atan(x) + log((1.+x)** 2. * (1. + x**2.)/ 8.)
      psim  = 3.14159265/2. - atan(x) + 2.*log((1+x)/2.) + log((1+x*x)/2.)
    else
      psim  = - 4.7 * zeta
    end if

    return
  end function psim

  !
  ! ----------------------------------------------------------------------
  ! Malte: Integrated stability function for heat
  !
  function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 15. * zeta) ** (0.25)
      psih  = 2. * log( (1. + x ** 2.) / 2. )
    else
      psih  = - 4.7 * zeta
    end if

    return
  end function psih

  !
  ! ----------------------------------------------------------------------
  ! Malte: Calculate the surface humidity assuming saturation. 
  ! DALES (vanHeerwarden)
  !
  subroutine qtsurf
    use defs, only : tmelt, R, Rm, rcp, p00
    use grid, only : nxp, nyp, vapor, psrf

    implicit none
    real       :: exner, tsurf, qsatsurf, surfwet, es, qtsl
    integer    :: i,j

    qtsl = 0.
    do j=3,nyp-2
      do i=3,nxp-2
        exner      = (psrf / p00)**rcp
        tsurf      = tskin(i,j) * exner
        es = 0.611e3 * exp(17.2694 * (tsurf - tmelt) / (tsurf - 35.86))

        qsatsurf   = R / Rm * es / psrf
        surfwet    = ra(i,j) / (ra(i,j) + rs(i,j))
        qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * vapor(2,i,j)
      end do
    end do

    return
  end subroutine qtsurf

  !
  ! ----------------------------------------------------------
  ! Malte: LSM to calculate temperature and moisture at the surface
  ! DALES (vanHeerwaarden)
  ! "PLEASE NOTE: LSM variables are defined and initialized in lsmdata.f90
  !
  subroutine lsm

    use defs, only: p00, stefan, rcp, cp, R, alvl, rowt, g
    use grid, only: nzp, nxp, nyp, a_up, a_vp, a_theta, vapor, liquid, zt, &
                    psrf, th00, umean, vmean, dn0, iradtyp, dt, &
		    a_lflxu, a_lflxd, a_sflxu, a_sflxd, nstep, press

    integer  :: i, j, k
    real     :: f1, f2, f3, f4, fsoil !Correction functions for Jarvis-Stewart
    real     :: lflxu_av, lflxd_av, sflxu_av, sflxd_av
    real     :: tsurfm, Tatm, qskinn, exnerair, exner, pair, thetaavg
    real     :: e, esat, qsat, desatdT, dqsatdT, Acoef, Bcoef
    real     :: fH, fLE, fLEveg, fLEsoil, fLEliq, LEveg, LEsoil, LEliq
    real     :: Wlmx, rk3coef=0


    thetaavg  = sum(a_theta(2,3:(nxp-2),3:(nyp-2)))/(nxp-4)/(nyp-4)
    pair      = psrf*exp(-1.*g*(zt(2)-zt(1))/2/R/thetaavg)
    exner     = (psrf/p00)**rcp
    exnerair  = (pair/p00)**rcp   

    !"1.0 - Compute water content per layer
    do j = 3,nyp-2
      do i = 3,nxp-2

        phitot(i,j) = 0.0
        do k = 1, ksoilmax
          phitot(i,j) = phitot(i,j) + phiw(i,j,k) * dzsoil(k)
        end do
        phitot(i,j) = phitot(i,j) / zsoil(ksoilmax)

        do k = 1, ksoilmax
          phifrac(i,j,k) = phiw(i,j,k)*dzsoil(k) / zsoil(ksoilmax) / phitot(i,j)
        end do

      end do
    end do

    do j = 3, nyp-2
      do i = 3, nxp-2

        !" 1.1 - Calculate net radiation (average in time (nradtime))
        if(iradtyp == 4) then
            if(nstep == 1) then
              sflxd_avn(2:nradtime,i,j) = sflxd_avn(1:nradtime-1,i,j)
              sflxu_avn(2:nradtime,i,j) = sflxu_avn(1:nradtime-1,i,j)
              lflxd_avn(2:nradtime,i,j) = lflxd_avn(1:nradtime-1,i,j)
              lflxu_avn(2:nradtime,i,j) = lflxu_avn(1:nradtime-1,i,j)

              sflxd_avn(1,i,j) = a_sflxd(1,i,j)
              sflxu_avn(1,i,j) = a_sflxu(1,i,j)
              lflxd_avn(1,i,j) = a_lflxd(1,i,j)
              lflxu_avn(1,i,j) = a_lflxu(1,i,j)
            end if

            !sflxd_av = sum(sflxd_avn(:,i,j))/nradtime
            !sflxu_av = sum(sflxu_avn(:,i,j))/nradtime
            !lflxd_av = sum(lflxd_avn(:,i,j))/nradtime
            !lflxu_av = sum(lflxu_avn(:,i,j))/nradtime

            sflxd_av = sum(sflxd_avn(:,:,:))/nradtime/(nxp-4)/(nyp-4)
            sflxu_av = sum(sflxu_avn(:,:,:))/nradtime/(nxp-4)/(nyp-4)
            lflxd_av = sum(lflxd_avn(:,:,:))/nradtime/(nxp-4)/(nyp-4)
            lflxu_av = sum(lflxu_avn(:,:,:))/nradtime/(nxp-4)/(nyp-4)

            Qnet(i,j) = (sflxd_av - sflxu_av + lflxd_av - lflxu_av)
            
        else
          !" Not using full radiation: use average radiation from Namelist
          Qnet(i,j) = Qnetav
        end if

        if ((nstep==1) .and. (i==10) .and. (j==10) .and. .false.) then
        print*,"Qnet***************sflxd*******************a_sflxd*******************lflxd*******************lflxu
          print*,Qnet(10,10),sflxd_av,a_sflxd(2,i,j),lflxd_av,lflxu_av
        end if

        !" 2.1 - Calculate the surface resistance with vegetation

        !" a) Stomatal opening as a function of incoming short wave radiation
        if (iradtyp == 4) then
          f1  = 1. / min(1., (0.004 * max(0.,sflxd_av) + 0.05) / &
                (0.81 * (0.004 * max(0.,sflxd_av) + 1.)))
        else
          f1  = 1.
        end if

        !" b) Soil moisture availability
        f2  = (phifc - phiwp) / (phitot(i,j) - phiwp)
        !" Put upper bound f2 in case of very dry soils and prevent less than 1
        f2  = max(f2, 1.)
        f2  = min(1.e8, f2)

        !" c) Response of stomata to vapor deficit of atmosphere
        Tatm    = a_theta(2,i,j)*exnerair
        esat = 0.611e3 * exp(17.2694 * (Tatm - 273.16) / (Tatm - 35.86))
        e    = vapor(2,i,j) * psrf / 0.622
        f3   = 1. / exp(-gD(i,j) * (esat - e) / 100.)

        !" d) Response to temperature
        f4      = 1./ (1. - 0.0016 * (298.0 - Tatm) ** 2.)
        rsveg(i,j)  = rsmin(i,j) / LAI(i,j) * f1 * f2 * f3 * f4

        !" 2.2 - Calculate soil resistance based on ECMWF method
        fsoil  = (phifc - phiwp) / (phiw(i,j,1) - phiwp)
        fsoil  = max(fsoil, 1.)
        fsoil  = min(1.e8, fsoil)
        rssoil(i,j) = rssoilmin(i,j) * fsoil

        !" 2.3 - Calculate the heat transport properties of the soil.
        !" Put in init function, as we don't have prognostic soil moisture atm.

        !"Save temperature and liquid water from previous timestep 
        if (nstep == 1) then
          tskinm(i,j) = tskin(i,j)
          Wlm(i,j)    = Wl(i,j)
        end if

        !"Solve the surface temperature implicitly including variations in LWout
        tsurfm  = tskinm(i,j) * exner
        esat    = 0.611e3 * exp(17.2694 * (tsurfm-273.16) / (tsurfm-35.86))
        qsat    = 0.622 * esat / psrf
        desatdT = esat * (17.2694 / (tsurfm-35.86) - 17.2694 *  &
                  (tsurfm-273.16) / (tsurfm-35.86)**2.)
        dqsatdT = 0.622 * desatdT / psrf

        !" First, remove LWup from Qnet calculation
        Qnet(i,j) = Qnet(i,j) + stefan * (tskinm(i,j)*exner)**4.
 
        !" Allow for dew fall and calculate dew water on leaves
        if(qsat - vapor(2,i,j) < 0.) then
          rsveg(i,j)  = 0.
          rssoil(i,j) = 0.
        end if

        Wlmx      = LAI(i,j) * Wmax
        Wl(i,j)   = min(Wl(i,j), Wlmx)
        cliq(i,j) = Wl(i,j) / Wlmx

        !" Calculate coefficients for surface fluxes
        fH      = 0.5*(dn0(1)+dn0(2)) * cp / ra(i,j) 
        fLEveg  = (1. - cliq(i,j)) *cveg(i,j) * (0.5* (dn0(1)+dn0(2))) * alvl &
                  /(ra(i,j) + rsveg(i,j))
        fLEsoil = (1. - cveg(i,j))            * (0.5* (dn0(1)+dn0(2))) * alvl &
                  /(ra(i,j) + rssoil(i,j))
        fLEliq  = (cliq(i,j) * cveg(i,j))     * (0.5* (dn0(1)+dn0(2))) * alvl &
                  /(ra(i,j))
        fLE     = fLEveg + fLEsoil + fLEliq

        !" Weighted timestep in runge-kutta scheme
        rk3coef = dt / (4- dble(nstep))

        !" Compute skin temperature from linarized surface energy balance
        Acoef   = Qnet(i,j) - stefan * (tskinm(i,j)*exner) **4. &
                  + 4. * stefan * (tskinm(i,j)*exner)**4. + fH*Tatm &
                  + fLE * (dqsatdT* (tskinm(i,j)*exner) - qsat + vapor(2,i,j)) &
                  + lambdaskin(i,j) * tsoil(i,j,1)
        Bcoef   = 4. * stefan * (tskinm(i,j)* exner) ** 3. + fH &
                  + fLE * dqsatdT + lambdaskin(i,j)

        if (Cskin(i,j) == 0.) then
           tskin(i,j) = Acoef * Bcoef ** (-1.) / exner
        else
           tskin(i,j) = (1. + rk3coef/Cskin(i,j) * Bcoef) ** (-1.) / exner &
                       * ((tskinm(i,j)*exner) + rk3coef/Cskin(i,j) * Acoef) 
        end if

        Qnet(i,j) = Qnet(i,j) - (stefan* (tskinm(i,j)*exner)**4.  &
                    + 4.*stefan * (tskinm(i,j)*exner)**3. *(tskin(i,j)*exner - &
                    tskinm(i,j)*exner))

        G0(i,j)   = lambdaskin(i,j) * ( tskin(i,j) * exner - tsoil(i,j,1) )

        qskinn    = (dqsatdT * (tskin(i,j)*exner - tskinm(i,j)*exner) + qsat)
        LE(i,j)   = - fLE     * ( vapor(2,i,j) - qskinn)

        LEveg     = - fLEveg  * ( vapor(2,i,j) - qskinn)
        LEsoil    = - fLEsoil * ( vapor(2,i,j) - qskinn)
        LEliq     = - fLEliq  * ( vapor(2,i,j) - qskinn)

        if(LE(i,j) == 0.) then
          rs(i,j) = 1.e8           ! "WHY???
        else
          rs(i,j) = -0.5*(dn0(1)+dn0(2)) * alvl * (vapor(2,i,j) - qskinn) &
                    / LE(i,j) - ra(i,j) 
        end if

        H(i,j)    = - fH  * ( a_theta(2,i,j)*exner - tskin(i,j)*exner )
        tendskin(i,j) = Cskin(i,j)*(tskin(i,j) - tskinm(i,j)) * exner / rk3coef

        !" In case of dew formation, allow all water to enter skin reservoir Wl
        if(qsat - vapor(2,i,j) < 0.) then
          Wl(i,j) =  Wlm(i,j) - rk3coef*((LEliq + LEsoil + LEveg)/(rowt * alvl))
        else
          Wl(i,j) =  Wlm(i,j) - rk3coef*(LEliq / (rowt * alvl))
        end if

        !" Save temperature and liquid water from previous timestep 
        if(nstep == 1) then
          tsoilm(i,j,:) = tsoil(i,j,:)
          phiwm(i,j,:)  = phiw(i,j,:)
        end if

        !" Calculate soil heat capacity and conductivity(based on water content)
        do k = 1, ksoilmax
          pCs(i,j,k)    = (1. - phi) * pCm + phiw(i,j,k) * pCw
          Ke            = log10(max(0.1,(phiw(i,j,k)/phi))) + 1.
          lambda(i,j,k) = Ke * (lambdasat - lambdadry) + lambdadry
        end do
     
        !" Calculate soil heat conductivity at half levels
        do k = 1, ksoilmax-1
          lambdah(i,j,k) = 0.5 * (lambda(i,j,k) * dzsoil(k) &
                            + lambda(i,j,k+1) * dzsoil(k+1)) / dzsoilh(k)
        end do
        lambdah(i,j,ksoilmax) = lambda(i,j,ksoilmax)

        !" Calculate soil moisture conductivity and difusivity 
        do k = 1, ksoilmax
          gammas(i,j,k)  = gammasat * (phiw(i,j,k) / phi) ** (2. * bc + 3.)
          lambdas(i,j,k) = bc * gammasat * (-1.) * psisat / phi &
                           * (phiw(i,j,k) / phi) ** (bc + 2.)
        end do

        !" Calculate soil moisture conductivity and difusivity at half levels
        do k = 1, ksoilmax-1
          lambdash(i,j,k) = 0.5 * (lambdas(i,j,k) * dzsoil(k) &
                            + lambdas(i,j,k+1) * dzsoil(k+1)) / dzsoilh(k)
          gammash(i,j,k)  = 0.5 * (gammas(i,j,k)  * dzsoil(k) &
                            + gammas(i,j,k+1)  * dzsoil(k+1)) / dzsoilh(k)
        end do
        lambdash(i,j,ksoilmax) = lambdas(i,j,ksoilmax)

        !" Solve the diffusion equation for the heat transport
        tsoil(i,j,1) = tsoilm(i,j,1) + rk3coef * ( lambdah(i,j,1) * (tsoil(i,j,2) - tsoil(i,j,1)) / dzsoilh(1) + G0(i,j) ) / dzsoil(1) / pCs(i,j,1)

        do k = 2, ksoilmax-1
          tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef / pCs(i,j,k) * ( lambdah(i,j,k) * (tsoil(i,j,k+1) - tsoil(i,j,k)) / dzsoilh(k) - lambdah(i,j,k-1) * (tsoil(i,j,k) - tsoil(i,j,k-1)) / dzsoilh(k-1) ) / dzsoil(k)
        end do

        tsoil(i,j,ksoilmax) = tsoilm(i,j,ksoilmax) + rk3coef / pCs(i,j,ksoilmax) * ( lambda(i,j,ksoilmax) * (tsoildeep(i,j) - tsoil(i,j,ksoilmax)) / dzsoil(ksoilmax) - lambdah(i,j,ksoilmax-1) * (tsoil(i,j,ksoilmax) - tsoil(i,j,ksoilmax-1)) / dzsoil(ksoilmax-1) ) / dzsoil(ksoilmax)

        !" Solve the diffusion equation for the moisture transport (closed bottom for now)
        phiw(i,j,1) = phiwm(i,j,1) + rk3coef * ( lambdash(i,j,1) * (phiw(i,j,2) - phiw(i,j,1)) / dzsoilh(1) - gammash(i,j,1) - (phifrac(i,j,1) * LEveg + LEsoil) / (rowt*alvl)) / dzsoil(1)

        do k = 2, ksoilmax-1
          phiw(i,j,k) = phiwm(i,j,k) + rk3coef * ( lambdash(i,j,k) * (phiw(i,j,k+1) - phiw(i,j,k)) / dzsoilh(k) - gammash(i,j,k) - lambdash(i,j,k-1) * (phiw(i,j,k) - phiw(i,j,k-1)) / dzsoilh(k-1) + gammash(i,j,k-1) - (phifrac(i,j,k) * LEveg) / (rowt*alvl)) / dzsoil(k)
        end do

        phiw(i,j,ksoilmax) = phiwm(i,j,ksoilmax) + rk3coef * (- lambdash(i,j,ksoilmax-1) * (phiw(i,j,ksoilmax) - phiw(i,j,ksoilmax-1)) / dzsoil(ksoilmax-1) + gammash(i,j,ksoilmax-1) - (phifrac(i,j,ksoilmax) * LEveg) / (rowt*alvl) ) / dzsoil(ksoilmax)

      end do
    end do
      
    call qtsurf

        if (nstep==1 .and. .false.) then
        print*,"********************************"
        print*,"timestep dt (0,...,dtlong)",dt
        print*,"********************************"
        print*,"Tskin (potential temp skin)",maxval(tskin(3:nxp-2,3:nyp-2))
        print*,"Tskin*exner (temp skin)",maxval(tskin(3:nxp-2,3:nyp-2))*exner
        print*,"Tskinm (potential temp skin, previous timestep)",maxval(tskinm(3:nxp-2,3:nyp-2))
        print*,"tskinm*exner (temp skin, previous timestep)",maxval(tskinm(3:nxp-2,3:nyp-2))*exner
        print*,"a_theta (potential temp first model layer)",maxval(a_theta(1,3:nxp-2,3:nyp-2))
        print*,"a_theta*exner (temp first model layer)",maxval(a_theta(1,3:nxp-2,3:nyp-2)*exner)
        print*,"delta T",minval((a_theta(1,3:nxp-2,3:nyp-2)*exner-tskin(3:nxp-2,3:nyp-2))*exner)
        print*,"********************************"
        print*,"lambdaskin (lambdaskinav=5.)",minval(lambdaskin(3:nxp-2,3:nyp-2))
        print*,"lambda (conductivity 1st layer, dry: 0.19)",minval(lambda(3:nxp-2,3:nyp-2,1))
        print*,"lambdah (conductivity at half levels)",minval(lambdah(3:nxp-2,3:nyp-2,1))
        print*,"lambda (conductivity 2nd layer, dry: 0.19)",minval(lambda(3:nxp-2,3:nyp-2,2))
        print*,"lambdas 1st layer",minval(lambdas(3:nxp-2,3:nyp-2,1))
        print*,"lambdash (at 1st half level)",minval(lambdash(3:nxp-2,3:nyp-2,1))
        print*,"lambdas 2nd layer",minval(lambdas(3:nxp-2,3:nyp-2,2))
        print*,"lambdasat",lambdasat
        print*,"pCs (2.2e6,....,4.2e6)",minval(pCs(3:nxp-2,3:nyp-2,1))
        print*,"Ke(0,...,1)",Ke
        print*,"********************************"
        print*,"phiw (water content 1st layer, 0.3)",minval(phiw(3:nxp-2,3:nyp-2,1))
        print*,"phiw (water content 2nd layer, 0.3)",minval(phiw(3:nxp-2,3:nyp-2,2))
        print*,"phiw (water content 2nd layer, 0.3)",minval(phiw(3:nxp-2,3:nyp-2,3))
        print*,"phiw (water content 2nd layer, 0.3)",minval(phiw(3:nxp-2,3:nyp-2,4))
        print*,"tsoil (soil temperature 1st layer, 290K)",minval(tsoil(3:nxp-2,3:nyp-2,1))
        print*,"tsoil (soil temperature 2nd layer, 287K)",minval(tsoil(3:nxp-2,3:nyp-2,2))
        print*,"tsoil (soil temperature 3rd layer, 285K)",minval(tsoil(3:nxp-2,3:nyp-2,3))
        print*,"tsoil (soil temperature 4th layer, 283K)",minval(tsoil(3:nxp-2,3:nyp-2,4))
        print*,"tsoil (soil temperature deeeeeeep, 283K)",minval(tsoildeep(3:nxp-2,3:nyp-2))
        print*,"********************************"
        print*,"Qnet(450 W/m2)",sum(Qnet(3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)
        print*,"G0 (min)",minval(G0(3:nxp-2,3:nyp-2))
        print*,"G0 (max)",maxval(G0(3:nxp-2,3:nyp-2))
        print*,"Tendskin (min)",minval(tendskin(3:nxp-2,3:nyp-2))
        print*,"Tendskin (max)",maxval(tendskin(3:nxp-2,3:nyp-2))
        print*,"********************************"
        print*,"Qnet-H-LE-G-tendSkin (min) !=! 0 -->",minval(Qnet(3:nxp-2,3:nyp-2)-(H(3:nxp-2,3:nyp-2)+LE(3:nxp-2,3:nyp-2)+G0(3:nxp-2,3:nyp-2))-tendskin(3:nxp-2,3:nyp-2))
        print*,"Qnet-H-LE-G-tendSkin (max) !=! 0 -->",maxval(Qnet(3:nxp-2,3:nyp-2)-(H(3:nxp-2,3:nyp-2)+LE(3:nxp-2,3:nyp-2)+G0(3:nxp-2,3:nyp-2))-tendskin(3:nxp-2,3:nyp-2))
        print*,"********************************"
        print*,"********************************"
        print*,"SH (min,lsm):",minval(H(3:nxp-2,3:nyp-2))
        print*,"SH (max,lsm):",maxval(H(3:nxp-2,3:nyp-2))
        print*,"LH (min,lsm):",minval(LE(3:nxp-2,3:nyp-2))
        print*,"LH (max,lsm):",maxval(LE(3:nxp-2,3:nyp-2))
        print*,"********************************"
        end if

  end subroutine lsm

end module srfc

