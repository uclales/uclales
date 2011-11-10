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
  real    :: zrough =  0.1
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
         uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc
    use thrm, only: rslf
    use stat, only: sfc_stat, sflg
    use mpi_interface, only : nypg, nxpg, double_array_par_sum

    implicit none
 
    real, optional, intent (inout) :: sst
    real :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp)       &
         ,wspd(nxp,nyp), bfct(nxp,nyp), cm(nxp,nyp), cs(nxp,nyp), ra(nxp,nyp) &
	 ,tskin(nxp,nyp), qskin(nxp,nyp), dudz(nxp,nyp), dvdz(nxp,nyp) &
	 ,dthldz(nxp,nyp), dqtdz(nxp,nyp), z0m(nxp,nyp), z0h(nxp,nyp) &
         ,ustar(nxp,nyp), phimzf(nxp,nyp), phihzf(nxp,nyp), obl(nxp,nyp)
    integer :: i, j, iterate
    real    :: zs, bflx0,bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum
    real (kind=8) :: bfl(2), bfg(2)


    select case(isfctyp)

  !
  ! ----------------------------------------------------------------------
  ! set surface gradients
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
  ! get fluxes from profiles
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
	print*,wt_sfc(55,3),wq_sfc(55,3)
	

  !
  ! ----------------------------------------------------------------------
  ! get fluxes from bulk formulae with coefficients given by
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
  ! fix surface temperature to yield a constant surface buoyancy flux
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
  ! xxxx
  !
   case(5)
       Vbulk = 0.01

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
       sst=289.
       iterate=0
       bflx0=0.0007
       bflx=g/bfg(1)*Vbulk*((sst-bfg(1))+ep2*bfg(1)*(rslf(psrf,sst)-bfg(2)))
       do while (abs(bflx-bflx0)>0.00001.and.iterate<800)
        sst=sst+0.01
        bflx=g/bfg(1)*Vbulk*((sst-bfg(1))+ep2*bfg(1)*(rslf(psrf,sst)-bfg(2)))
        iterate=iterate+1
       end do
       if (iterate.eq.799) print*,'WARNING'
       if (sst.lt.289) print*,'WRONG SST'
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
  ! Malte: Get surface fluxes using a land surface model (vanHeerwaarden)
  !
   case(6)

       !roughness length for momentum and heat (future: from NAMELIST)
       z0m(:,:) = 0.1
       z0h(:,:) = 0.1
       zs = zrough

       do j=3,nyp-2
          do i=3,nxp-2
             wspd(i,j) = max(0.1,                                            &
                         sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
             dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
             drdz(i,j) = vapor(2,i,j) - rslf(psrf,sst)
          end do
       end do
	!print*,wspd(55,3),dtdz(55,3),drdz(55,3)

       call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar    &
            ,a_rstar,obl)
	!print*,a_ustar(55,3),a_tstar(55,3),a_rstar(55,3),obl(55,3)

       !Calculate the drag coefficients and aerodynamic resistance
       do j=2,nyp-2
          do i=2,nxp-2

             cm(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
			obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
             cs(i,j) =  vonk**2. / (log(zt(2)/z0m(i,j)) - psim(zt(2) / &
			obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zt(2) / &
			z0h(i,j)) - psih(zt(2) / obl(i,j)) + psih(z0h(i,j) / &
			obl(i,j)))

             ra(i,j) = 1. / (cs(i,j)* wspd(i,j))

          end do
       end do

       !Get skin temperature and humidity from land surface model
       tskin(:,:) = sst*(p00/psrf)**rcp 	!call do_lsm
       qskin(:,:) = rslf(psrf,sst)		!call qtsurf
	!print*,tskin(55,3),qskin(55,3)

       !Calculate the surface fluxes with bulk law (Fairall, 2003)
       do j=3, nyp-2
          do i=3, nxp-2

	     wt_sfc(i,j) = -(a_theta(2,i,j) - tskin(i,j)) / ra(i,j) 
             wq_sfc(i,j) = -(vapor(2,i,j) - qskin(i,j)) / ra(i,j)
             uw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                  *(a_up(2,i,j)+umean)/wspd(i,j)
             vw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                  *(a_vp(2,i,j)+vmean)/wspd(i,j)
             ww_sfc(i,j)  = 0.
             a_rstar(i,j) = -wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = -wt_sfc(i,j)/a_ustar(i,j)

          end do
       end do

	!print*,a_rstar(55,3),a_tstar(55,3)
	print*,wt_sfc(55,3),wq_sfc(55,3)
	!print*,uw_sfc(55,3),vw_sfc(55,3)
	!print*,ra(55,3),cm(55,3),cs(55,3)

	
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
  subroutine srfcscls(n2,n3,z,z0,th00,u,dth,drt,ustar,tstar,rstar,obl)

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
    real, intent(inout) :: obl(n2,n3)    ! Obukhov Length

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
          ! Stable case
          !
          if (dtv > 0.) then
             x     = (betg*u(i,j)**2)/dtv
             y     = (am - 0.5*x)/lnz
             x     = (x*ah - am**2)/(lnz**2)
             lmo   = -y + sqrt(x+y**2)
             zeta  = z/lmo
             ustar(i,j) =  vonk*u(i,j)/(lnz + am*zeta)
             tstar(i,j) = (vonk*dtv/(lnz + ah*zeta))/pr
          !
          ! Neutral case
          ! 
          elseif (dtv == 0.) then
             ustar =  vonk*u(i,j)/lnz
             tstar =  vonk*dtv/(pr*lnz)
             lmo = -1.e10
          !
          ! Unstable case, start iterations from values at previous tstep, 
          ! unless the sign has changed or if it is the first call, then 
          ! use neutral values.
          !
          else
             if (first_call .or. tstar(i,j)*dtv <= 0.) then
                ustar(i,j) = u(i,j)*klnz
                tstar(i,j) = (dtv*klnz/pr)
		lmo = -1.e10

             end if

             do iterate = 1,3
                lmo   = betg*ustar(i,j)**2/(vonk*tstar(i,j))
                zeta  = z/lmo
                x     = sqrt( sqrt( 1.0 - bm*zeta ) )
                psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + cnst1
                y     = sqrt(1.0 - bh*zeta)
                psi2  = log(1.0 + y) + cnst2
                ustar(i,j) = u(i,j)*vonk/(lnz - psi1)
                tstar(i,j) = (dtv*vonk/pr)/(lnz - psi2)
             end do
          end if

          obl(i,j) = lmo
          rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
          tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)

		if (i .eq. 50) then		
		!print*,"dtv(50,3)",dtv
		!print*,"obl(50,3)",obl(50,3)
		!print*, zeta
		end if


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


  !
  ! ----------------------------------------------------------
  ! Malte: Integrated stability function psi_m
  !
  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    real    :: cnst1, cnst2
    cnst2 = -log(2.)
    cnst1 = 3.14159/2. + 3.*cnst2

    if(zeta <= 0) then
      x    = (1. - 19.3 * zeta) ** (0.25)
      psim = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + cnst1
      !psim = 3.14159265/2. - 2.*atan(x) + log( (1.+x)**2. * (1. + x ** 2.) / 8.)
    else
      psim  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    end if

    return
  end function psim


  !
  ! ----------------------------------------------------------------------
  ! Malte: Integrated stability function psi_h
  !
  function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    real    :: cnst1, cnst2
    cnst2 = -log(2.)
    cnst1 = 3.14159/2. + 3.*cnst2

    if(zeta <= 0) then
      x     = sqrt(1.0 - 12.*zeta)
      psih  = log(1.0 + x) + cnst2
      !x     = (1. - 16. * zeta) ** (0.25)
      !psih  = 2. * log( (1. + x ** 2.) / 2. )
    else
      psih  = -2./3.*(zeta - 5./0.35)*exp(-0.35 * zeta) - (1. + (2./3.)*zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    end if

    return
  end function psih

end module srfc



