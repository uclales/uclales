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
!!!!!!!!!!!!!!!!




module mcrp

  use wolken_driver,     only: cloud_type => wolke_typ

  use mpi_interface,  only : myid, double_scalar_par_max, double_scalar_par_sum
  use defs, only : tmelt,alvl, alvi,rowt,roice, pi, Rm, cp,t_hn
  use grid, only : dt,nstep,rkbeta,rkalpha, dxi, dyi ,dzi_t, nxp, nyp, nzp,nfpt, a_pexnr, pi0,pi1,a_rp, a_tp, th00, ccn,    &
       dn0, pi0,pi1, a_rt, a_tt,a_rpp, a_rpt, a_npp, a_npt, vapor, liquid, a_wp, zm,      &
       a_theta, a_scr1, a_scr2, a_scr7, rsi, &
       a_zpp ,   a_zpt    , & ! rain reflectivity
       a_ricep , a_ricet  , & ! ice mixing ratio
       a_nicep , a_nicet  , & ! ice number concentration
       a_rsnowp, a_rsnowt , & ! snow mass
       a_nsnowp, a_nsnowt , & ! snow number
       a_rgrp,   a_rgrt,    & ! graupel mass
       a_ngrp,   a_ngrt,    & ! graupel number
       a_rhailp, a_rhailt,  & ! hail mass
       a_nhailp, a_nhailt,  & ! hail number
       prc_c, prc_r, prc_i, prc_s, prc_g, prc_h, & 
       prc_acc, rev_acc, a_rct, cnd_acc, cev_acc, aut_acc, acc_acc, eps_max, prc_lev, lmptend, evap, &
       mom3

  USE parallele_umgebung, ONLY: isIO,double_global_maxval,global_maxval,global_minval,global_maxval_stdout,global_sumval_stdout
  USE modcross, ONLY: calcintpath

  use thrm, only : thermo, fll_tkrs, esl, esi
  use util, only : get_avg3, azero, sclrset
!  implicit none

  logical            :: lpartdrop      = .false.        !< Switch for rain drop like particles (namelist parameter)
  logical, parameter :: ldropfeedback  = .false.        !< Switch for feedback of LD's on bulk variables N, L and Z
                                                        !  (needs mom3=.true. in Namelist)

  logical, parameter :: droplet_sedim = .true., khairoutdinov = .false., &
                        kessler = .false., khairoutdinov_au = .false.

  real, dimension(:),allocatable :: convice,convliq,frho
  integer :: iturbkern ! 0) no turb. 1) Onishi-2015 2) Ayala-2015 
  real :: timenuc = 60
  real :: nin_set = 1.7e3
  real, dimension(:,:,:),allocatable :: a_npauto
  !
  ! drop sizes definition is based on vanZanten (2005)
  ! cloud droplets' diameter: 2-50 e-6 m
  ! drizzle drops' diameter: 50-1000 e-6 m
  !
  real, parameter    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
  real, parameter    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]

  real, parameter :: eps0 = 1e-20       ! small number
  real, parameter :: eps1 = 1e-9        ! small number
!   real, parameter :: rthres = 0.0          ! small number
  real, parameter :: rthres = 1e-20        ! small number
  real, parameter :: rho_0 = 1.21       ! air density at surface

  real, parameter :: prw = pi * rowt / 6.
  real, parameter :: nu_l = 1.460e-5     !..kinem. visc. of air
  real, parameter :: d_v  = 3.000e-5     !..diff. of water vapor
  real, parameter :: k_t  = 2.500e-2     !..heat conductivity

  real, parameter :: n_sc = 0.710        !..schmidt-number(pk, p.541)
  real, parameter :: n_f  = 0.333        !..exponent of n_sc in the vent-coeff. (pk, p.541)
  real, parameter :: m_f  = 0.500        !..exponent of n_re in the vent-coeff. (pk, p.541)
  !
  real, parameter :: a_e  = 2.18745584e1 !..const. in saturation pressure wrt ice
  real, parameter :: a_w  = 1.72693882e1 !..const. in saturation pressure wrt water
  real, parameter :: b_e  = 7.66000000e0 !..const. in saturation pressure wrt ice
  real, parameter :: b_w  = 3.58600000e1 !..const. in saturation pressure wrt water
  real, parameter :: e_3  = 6.10780000e2 !..saturation pressure at melting temp.

  !..Bjorn's mu-Dm-relation
  real, parameter :: cmur1 = 10.0    ! mu-Dm-relation for rain following
  real, parameter :: cmur2 = 1.20e+3 ! Milbrandt&Yau 2005, JAS, but with
  real, parameter :: cmur3 = 1.4e-3  ! revised constants

  !..Axel's mu-Dm-relation for raindrops based on 1d-bin model
  real, parameter :: rain_cmu0 = 6.0             
  real, parameter :: rain_cmu1 = 30.0             
  real, parameter :: rain_cmu2 = 1.00d+3         
  real, parameter :: rain_cmu3 = 1.10d-3   ! D_eq,breakup
  real, parameter :: rain_cmu4 = 1.0        
  logical, parameter :: mue_SB = .true.

  ! Milbrandt and Yau (2005, JAS) mu-Dm-relation for raindrops
  real, parameter :: cmy1 = 19.0    
  real, parameter :: cmy2 = 0.6e+3 
  real, parameter :: cmy3 = 1.8e-3  
  real, parameter :: cmy4 = 17.0  

  ! Milbrandt and McTaggart-Cowan (2010) mu-Dm-relation for raindrops
  real, parameter :: cmu_mc0 = 11.8
  real, parameter :: cmu_mc1 = 1.0d3
  real, parameter :: cmu_mc2 = 0.7
  real, parameter :: cmu_mc3 = 2.0

  ! look-up table for phillips et al. nucleation
  integer, parameter :: &
    ttmax  = 30,      &  ! sets limit for temperature in look-up table
    ssmax  = 60,      &  ! sets limit for ice supersaturation in look-up table
    ttstep = 2,       &  ! increment for temperature in look-up table
    ssstep = 1           ! increment for ice supersaturation in look-up table

  real, dimension(0:100,0:100), save :: &
    afrac_dust, &  ! look-up table of activated fraction of dust particles acting as ice nuclei
    afrac_soot, &  ! ... of soot particles
    afrac_orga     ! ... of organic material

  include 'phillips_nucleation_2010.incf'

  type particle
     character(10) :: name
     integer :: moments
     real :: nr
     real :: nu
     real :: mu
     real :: x_max !..maximum particle mass
     real :: x_min !..minimum particle mass
     real :: a_geo !..coeff. geometry
     real :: b_geo !..coeff. geometry = 1/3
     real :: a_vel !..coeff. fall velocity
     real :: b_vel !..coeff. fall velocity
     real :: s_vel !..Dispersion of fall velocity
     real :: a_ven !..coeff. ventilation param.
     real :: b_ven !..coeff. ventilation param.
     real :: cap   !..coeff. capacity

  end type particle
  type(particle) :: cldw,rain,ice,snow,graupel,hail

  integer :: istep = -999

  real :: rct_acc, rpt_acc

contains

  !
  ! ---------------------------------------------------------------------
  ! MICRO: sets up call to microphysics
  !
  subroutine micro(level,istp)

    integer, intent (in) :: level

    logical :: debug = .false.
    REAL(KIND=8) :: wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax, &
         &          precmax,rhomax,rhomin,ncmax,nimax,dzmin,tmax,tmin, &
         &          hlp1,hlp2,hlp3,hlp4,hlp5,hlp6,hlp7,hlp8,hsum
    integer, optional :: istp
    real, dimension(3:nxp-2,3:nyp-2) :: tmp
    integer :: n

    if (present(istp)) istep = istp


    if (mod(istep,100).eq.0) then
       debug = .true.
    else
       debug = .false.
    end if

    call fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,a_scr1,rs=a_scr2)
    select case (level)
    case(2)
       if (droplet_sedim) then
          call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,vapor,a_scr2,liquid,prc_c)
       end if
    case(3)
       if (.not.mom3) then
         call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,vapor,a_scr2,liquid,prc_c,a_rpp, &
              a_npp,a_rt,a_rpt,a_npt,a_scr7,prc_r)
       else
         call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,vapor,a_scr2,liquid,prc_c,a_rpp, &
              a_npp,a_rt,a_rpt,a_npt,a_scr7,prc_r,zp=a_zpp,zpt=a_zpt)
       end if
       do n =1,count(prc_lev>0)
         ! Christopher: implemented Axel's bug fix (also at the corresponding lines)
         prc_acc(:,:,n) = prc_acc(:,:,n) + (prc_c(prc_lev(n),:,:)+prc_r(prc_lev(n),:,:)) * dt / 3.
       end do
    case(5)
       call mcrph_sb(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,a_wp,vapor,liquid, &
            a_rpp,    a_npp    , & ! rain
            a_ricep , a_nicep  , & ! ice
            a_rsnowp, a_nsnowp , & ! snow
            a_rgrp,   a_ngrp,    & ! graupel
            a_rhailp, a_nhailp,  & ! hail
            a_rt,                & ! total water tendency
            a_rpt,    a_npt    , & ! rain
            a_ricet , a_nicet  , & ! ice
            a_rsnowt, a_nsnowt , & ! snow
            a_rgrt,   a_ngrt,    & ! graupel
            a_rhailt, a_nhailt,  & ! hail
            prc_c, prc_r, prc_i, prc_s, prc_g, prc_h)
      do n =1,count(prc_lev>0)
        prc_acc(:,:,n) = prc_acc(:,:,n) + (prc_c(prc_lev(n),:,:)+prc_r(prc_lev(n),:,:)+prc_i(prc_lev(n),:,:)          &
                           +prc_s(prc_lev(n),:,:)+prc_g(prc_lev(n),:,:)+prc_h(prc_lev(n),:,:)) * dt / 3.
      end do
   end select

  end subroutine micro
  !
  ! ---------------------------------------------------------------------
  ! MCRPH: calls microphysical parameterization
  !

  subroutine mcrph(level,n1,n2,n3,dn0,exner,pi0,pi1,thl,tlt,tk,vapor,rsat,rcld,prc_c, &
       rp, np, rtt,rpt,npt,dissip, prc_r,rsati, ricet,nicet,rsnowt,rgrpt,ricep,nicep,rsnowp,rgrpp, &
       prc_i, prc_s, prc_g, rct, zp, zpt)

    integer, intent (in) :: level,n1,n2,n3
    real, dimension(n1)      , intent (in)             :: dn0,pi0,pi1  !Density and mean pressures
    real, dimension(n1,n2,n3), intent (inout)          :: exner, & !exner function
         thl,   & ! Theta_l
         tk,    & ! temperature
         tlt,   & ! theta_l tendency
         vapor, & ! water vapor
         rsat,  & ! saturation mixing ration
         rcld,  & ! cloud water
         prc_c
    real, dimension(n1,n2,n3), intent (inout),optional :: np,&!rain drop number
         rp, &!rain water
         rtt,& !Total water tendency
         rpt,&!rain water tendency
         npt,&!rain droplet number tendency
         rct,&! cloud water tendency
         dissip, &
         prc_r, &
         zp,&! rain reflectivity
         zpt  ! rain reflectivity tendency   

    real, dimension(n1,n2,n3), intent (inout),optional :: rsati, ricet,nicet,rsnowt,rgrpt,&
         ricep,nicep,rsnowp,rgrpp

    real, dimension(n1,n2,n3), intent (inout),optional :: prc_i, prc_s, prc_g

    real, dimension(n1) :: tl,temp,rv,rc,nrain,rrain,zrain,nice, nin_active,rice,nsnow,rsnow,ngrp,rgrp,rs,rsi,s_i,r1
    integer :: i, j,n
    
    allocate(convice(n1),convliq(n1), frho(n1))

    if(lpartdrop .and. nstep==1) allocate(a_npauto(nzp,nxp,nyp))
    if(lpartdrop .and. nstep==1) a_npauto(:,:,:) = 0.

    evap = 0.
    do j=3,n3-2
       do i=3,n2-2
          if (.not.mom3) then
            call resetvar(rain,rp(1:n1,i,j),np(1:n1,i,j))
          else
            call resetvar(rain,rp(1:n1,i,j),np(1:n1,i,j),zp(1:n1,i,j))
            zrain = zp(1:n1,i,j)
          end if
          tl = thl(1:n1,i,j)
          temp = tk(1:n1,i,j)
          rv = vapor(1:n1,i,j)
          rc = rcld(1:n1,i,j)
          rs = rsat(1:n1,i,j)
          rrain = rp(1:n1,i,j)
          nrain = np(1:n1,i,j)
          convliq = alvl/(cp*(pi0+pi1+exner(1:n1,i,j))/cp)
          frho = (rho_0/dn0(1:n1))**0.35 * dn0(1:n1)

          ! autoconversion of cloud water
          call resetvar(cldw,rc)
          if(.not.mom3) then
             call resetvar(rain,rrain,nrain)             
             call auto_SB(n1,dn0,rc,rrain,nrain,tl,dissip(1:n1,i,j),i,j,au_acc=aut_acc(i,j),emax=eps_max(i,j),zm=zm)
          else
             call resetvar(rain,rrain,nrain,zrain)
             call auto_SB(n1,dn0,rc,rrain,nrain,tl,dissip(1:n1,i,j),i,j,zp=zrain,au_acc=aut_acc(i,j),emax=eps_max(i,j),zm=zm)
          end if

          if (.not.ldropfeedback) then ! if LDs feed back on dynamics
                                       ! there is no bulk raindrop phase anymore
                                       ! and accretion, evaporation and sedimentation of rain
                                       ! are done by the LDs

          ! accretion
          call resetvar(cldw,rc)
          if(.not.mom3) then
             call resetvar(rain,rrain,nrain)
             call accr_SB(n1,dn0,rc,rrain,nrain,tl,dissip(1:n1,i,j),ac_acc=acc_acc(i,j),zm=zm)
          else
             call resetvar(rain,rrain,nrain,zrain)
             call accr_SB(n1,dn0,rc,rrain,nrain,tl,dissip(1:n1,i,j),ac_acc=acc_acc(i,j),zm=zm,zp=zrain)
          end if

          ! sedimentation of rain
          if(.not.mom3) then
             call resetvar(rain,rrain,nrain)
             call sedim_rd(n1,dt,dn0,rrain,nrain,prc_r(1:n1,i,j))
          else
             call resetvar(rain,rrain,nrain,zrain)
             call sedim_rd(n1,dt,dn0,rrain,nrain,prc_r(1:n1,i,j),zp=zrain)
          end if

          ! evaporation of rain
          call resetvar(cldw,rc)
          if(.not.mom3) then
             call resetvar(rain,rrain,nrain)
             call wtr_dff_SB(n1,dn0,rrain,nrain,rc,rs,rv,tl,temp,ev=rev_acc(i,j),eva=evap(i,j,:),zm=zm)
          else
             call resetvar(rain,rrain,nrain,zrain)
             call wtr_dff_SB(n1,dn0,rrain,nrain,rc,rs,rv,tl,temp,ev=rev_acc(i,j),eva=evap(i,j,:),zm=zm,zp=zrain)
          end if

          end if ! end of IF(.not.ldropfeedback) 

          if (droplet_sedim) then
          ! sedimentation of cloud water
          call resetvar(cldw,rc)
          call sedim_cd(n1,dt,tl,rc,prc_c(1:n1,i,j))
          end if

          ! tendencies for Runge-Kutta
          tlt(2:n1,i,j) = tlt(2:n1,i,j)+(tl(2:n1) - thl(2:n1,i,j))/dt
          rtt(2:n1,i,j) = rtt(2:n1,i,j)+(rv(2:n1) - vapor(2:n1,i,j))/dt + (rc(2:n1) - rcld(2:n1,i,j))/dt
          rpt(2:n1,i,j) = max(rpt(2:n1,i,j) +(rrain(2:n1) - rp(2:n1,i,j))/dt,-rp(2:n1,i,j)/dt)
          npt(2:n1,i,j) = max(npt(2:n1,i,j) +(nrain(2:n1) - np(2:n1,i,j))/dt,-np(2:n1,i,j)/dt)
          if (mom3) then
            zpt(2:n1,i,j) = max(zpt(2:n1,i,j) +(zrain(2:n1) - zp(2:n1,i,j))/dt,-zp(2:n1,i,j)/dt)
          end if
          
          if (present(rct)) then
             rct(2:n1,i,j) = rct(2:n1,i,j) + (rc(2:n1) - rcld(2:n1,i,j))/dt
          end if
       end do
    end do

    deallocate(convice,convliq,frho)

  end subroutine mcrph

  subroutine wtr_dff_SB(n1,dn0,rp,np,rl,rs,rv,tl,tk,ev,eva,zm,zp)
    !
    ! ---------------------------------------------------------------------
    ! WTR_DFF_SB: calculates the evolution of the both number- and
    ! mass mixing ratio large drops due to evaporation in the absence of
    ! cloud water.
    !

    integer, intent (in) :: n1
    real, intent (in)    :: tk(n1),rs(n1),dn0(n1)
    real, intent (inout) :: rp(n1), np(n1),tl(n1),rv(n1),rl(n1)
    real, intent (inout), optional :: ev
    real, intent (inout), optional :: eva(n1)
    real, intent (in),    optional :: zm(n1)
    real, intent (inout), optional :: zp(n1)

    real, parameter     :: c_Nevap = 1.
    integer             :: k
    real                :: Xp, Dp, G, S, cerpt, cenpt, cezpt

    real, parameter :: a2 = 9.65       ! in SI [m/s]
    real, parameter :: c2 = 6e2        ! in SI [1/m]
    real, parameter :: Dv = 25.0e-6    ! in SI [m/s]

    logical, parameter :: oldevaporation = .false.

    real :: mue,lam,gfak,f_q,gamma_eva,b2

    b2 = a2*exp(c2*Dv)


    if (oldevaporation) then
      if (mom3) then
        print *, 'no evaporation with old sedimentation and 3 mom mcrp'
      end if
      do k=2,n1
         if (rp(k) > 0) then
            Xp = rp(k)/ (np(k)+eps0)
            Dp = ( Xp / prw )**(1./3.)
            G = 1. / (1. / (dn0(k)*rs(k)*Dv) + &
                 alvl*(alvl/(Rm*tk(k))-1.) / (Kt*tk(k)))
            S = rv(k)/rs(k) - 1.

            if (S < 0) then
               cerpt = 2. * pi * Dp * G * S * np(k) * dt
               cerpt = max (cerpt, -rp(k))
               cenpt = c_Nevap*cerpt * np(k) / rp(k)
               np(k)=np(k) + cenpt
               rp(k)=rp(k) + cerpt
               rv(k)=rv(k) - cerpt
               tl(k)=tl(k) + convliq(k)*cerpt
             
               if (present(ev)) ev=ev - cerpt * (zm(k)-zm(k-1))*dn0(k) / 3.
            end if
         end if
      end do
    else
       do k=2,n1

          S = rv(k)/rs(k) - 1.

          if (rp(k) > 0 .and. S < 0) then
             Xp  = rp(k)/ (np(k)+eps0)
             Dp  = ( Xp / prw )**(1./3.)
             G = 1.0 / ( alvl**2 / (Kt * Rm * tk(k)**2) + Rm * tk(k) / (Dv * esl(tk(k))) )

             if (mom3) then
                mue = rain_mue_z_inv(np(k), rp(k), zp(k))
             else
                mue = rain_mue_2mom(Dp)
             end if

             lam = (pi/6.* rowt &
                  &      * (mue+3.0)*(mue+2.0)*(mue+1.0) / Xp)**(1./3.)

             gfak =  0.1357940435E+01 &
                  &  + mue * ( +0.3033273220E+00  &
                  &  + mue * ( -0.1299313363E-01  &
                  &  + mue * ( +0.4002257774E-03  &
                  &  - mue * 0.4856703981E-05 ) ) )

             f_q  = rain%a_ven + rain%b_ven * 0.71**0.333 * sqrt(a2/nu_l)            &
                  &    * gfak / SQRT(lam)                                    &
                  &    * (1.0                                                &
                  &      - 1./2.  * (b2/a2)**1 * (lam/(1.*c2+lam))**(mue+5.0/2.0) &
                  &      - 1./8.  * (b2/a2)**2 * (lam/(2.*c2+lam))**(mue+5.0/2.0) &
                  &      - 1./16. * (b2/a2)**3 * (lam/(3.*c2+lam))**(mue+5.0/2.0) &
                  &      - 5./127.* (b2/a2)**4 * (lam/(4.*c2+lam))**(mue+5.0/2.0) &
                  &      )

             if (mom3) then
               ! three-moment scheme
               !gamma_eva = EXP(-0.1*mue) 
               gamma_eva = min((1.1e-3/Dp) * EXP(-0.2*mue),1.0)
               !gamma_eva = 1.0
             else
               if (mue_SB) then
                 ! two-moment scheme with S08 mue-lambda relation
                 !gamma_eva = EXP(-0.1*mue) 
                 gamma_eva = min((1.1e-3/Dp) * EXP(-0.2*mue),1.0)
                 !gamma_eva = 1.0
               else
                 ! two-moment scheme with MM10 mue-lambda relation
                 !gamma_eva = EXP(-0.1*mue) 
                 gamma_eva = min((1.1e-3/Dp) * EXP(-0.2*mue),1.0)
                 !gamma_eva = 1.0
               end if
             end if 
             !f_q = 1.0

             cerpt = 2. * pi * G * np(k) * (mue+1.0) / lam * f_q * S * dt
             cerpt = max (cerpt, -rp(k))
             cenpt = gamma_eva * cerpt / Xp
             if (mom3) then
               cezpt = (mue+4)/(mue+1) * cerpt * Xp
             end if

             np(k) = np(k) + cenpt
             rp(k) = rp(k) + cerpt
             rv(k) = rv(k) - cerpt
             tl(k) = tl(k) + convliq(k)*cerpt
             if (mom3) then
               zp(k) = zp(k) + cezpt
             end if
             
             if (present(ev)) ev = ev - cerpt * (zm(k)-zm(k-1))*dn0(k) / 3.

             if (present(eva)) eva(k) = cerpt
          end if
       end do
    end if

  end subroutine wtr_dff_SB

  subroutine auto_SB(n1,dn0,rc,rp,np,tl,diss,i,j,zp,au_acc,emax,zm)
    !
    ! ---------------------------------------------------------------------
    ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
    ! drizzle drops due to autoconversion. The autoconversion rate assumes
    ! f(x)=A*x**(nu_c)*exp(-Bx), an exponential in drop MASS x. It can
    ! be reformulated for f(x)=A*x**(nu_c)*exp(-Bx**(mu)), where formu=1/3
    ! one would get a gamma dist in drop diam -> faster rain formation.
    !
    
    integer, intent (in) :: n1,i,j
    real, intent (in)    :: dn0(n1), diss(n1)
    real, intent (inout) :: rc(n1), rp(n1), np(n1),tl(n1)
    real, intent (inout), optional :: zp(n1)
    real, intent (inout), optional :: au_acc, emax
    real, intent (in),    optional :: zm(n1)

    real            :: nu_c  = 0.           ! width parameter of cloud DSD
    real, parameter :: kc_0  = 9.44e+9      ! Long-Kernel
    real, parameter :: k_1  = 6.e+2        ! Parameter for phi function
    real, parameter :: k_2  = 0.68         ! Parameter for phi function

    real, parameter :: Cau = 4.1e-15 ! autoconv. coefficient in KK param.
    real, parameter :: Eau = 5.67    ! autoconv. exponent in KK param.
    real, parameter :: mmt = 1.e+6   ! transformation from m to \mu m

    ! Ayala kernel as in 2010 paper
    real, parameter :: kc0_a1 = 0.00489
    real, parameter :: kc0_a2 = -0.00042
    real, parameter :: kc0_a3 = -0.01400
    real, parameter :: kc0_b1 = 11.45
    real, parameter :: kc0_b2 = 9.68
    real, parameter :: kc0_b3 = 0.62
    real, parameter :: kc0_c1 = 4.82
    real, parameter :: kc0_c2 = 4.80
    real, parameter :: kc0_c3 = 0.76
    real, parameter :: kc0_bet = 0.00174

    ! Onishi kernel (of 29 July 2015)
    real, parameter :: kc1_a1 = 3.985e-03
    real, parameter :: kc1_a2 = 6.210e-03
    real, parameter :: kc1_a3 = 1.331e+00
    real, parameter :: kc1_b1 = 1.381e+01
    real, parameter :: kc1_b2 = 9.980e+00
    real, parameter :: kc1_b3 = 5.018e-01
    real, parameter :: kc1_c1 = 6.325e+00
    real, parameter :: kc1_c2 = -9.238e-01
    real, parameter :: kc1_c3 = -1.528e-01
    real, parameter :: kc1_bet = 2.026e-03

    ! Ayala kernel with 2009 eta and up to eps of 1000 (of 20150807)
    real, parameter :: kc2_a1 = 7.432e-04
    real, parameter :: kc2_a2 = -6.993e-05
    real, parameter :: kc2_a3 = -9.497e-02
    real, parameter :: kc2_b1 = 1.073e+01
    real, parameter :: kc2_b2 = 1.356e+01
    real, parameter :: kc2_b3 = 1.005e+00
    real, parameter :: kc2_c1 = 6.607e+00
    real, parameter :: kc2_c2 = 2.547e+00
    real, parameter :: kc2_c3 = 2.350e-01
    real, parameter :: kc2_bet = 3.480e-04

    real :: prey

    real, parameter :: csx = 0.23
    real, parameter :: ce = 0.93
    
    real, parameter :: alpha = 1.e-3 ! s-1 Kessler scheme
    real, parameter :: rc0 = 1.e-3   ! kg/kg Kessler scheme

    integer :: k
    real    :: k_au0,k_au, Xc, Dc, au, tau, phi, kc_alf, kc_rad, kc_sig, kc_bet, k_c, Re, epsilon, l
    real, dimension(n1) :: Resum, Recnt

    emax   = 0.0
!    if (present(au_acc)) au_acc = 0.0
    !
    ! Calculate the effect of turbulence on the autoconversion/
    ! accretion rate using equation (6) in Seifert et al. (2009)
    !

    ! 3mom paper uses nu=0 for dycoms, but nu=1 for RICO
    nu_c = 0.0 !! cldw%nu!1*dn0(k)*rc(k)+cldw%nu2

    k_au0  = kc_0 / (20.*cldw%x_max) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2

    if (iturbkern.ne.0) then
       if (iturbkern.eq.22) then
          ! Ayala as in 2010 QJ paper
          kc_alf = ( kc0_a1 + kc0_a2 * nu_c )/ ( 1. + kc0_a3 * nu_c )
          kc_rad = ( kc0_b1 + kc0_b2 * nu_c )/ ( 1. + kc0_b3 * nu_c )
          kc_sig = ( kc0_c1 + kc0_c2 * nu_c )/ ( 1. + kc0_c3 * nu_c )
          kc_bet = kc0_bet
          prey   = 0.25
       elseif (iturbkern.eq.1) then
          ! Onishi of August 2015
          kc_alf = ( kc1_a1 + kc1_a2 * nu_c )/ ( 1. + kc1_a3 * nu_c )
          kc_rad = ( kc1_b1 + kc1_b2 * nu_c )/ ( 1. + kc1_b3 * nu_c )
          kc_sig = ( kc1_c1 + kc1_c2 * nu_c )/ ( 1. + kc1_c3 * nu_c )
          kc_bet = kc1_bet
          prey   = -0.125
       elseif (iturbkern.eq.2) then
          ! Ayala kernel with 2009 eta and up to eps of 1000 (of 20150807)
          kc_alf = ( kc2_a1 + kc2_a2 * nu_c )/ ( 1. + kc2_a3 * nu_c )
          kc_rad = ( kc2_b1 + kc2_b2 * nu_c )/ ( 1. + kc2_b3 * nu_c )
          kc_sig = ( kc2_c1 + kc2_c2 * nu_c )/ ( 1. + kc2_c3 * nu_c )
          kc_bet = kc2_bet
          prey   = 0.25
       end if
    end if

    do k=2,n1-1
       if (rc(k) > 0.) then
          Xc = rc(k)/(cldw%nr+eps0)
          k_c = kc_0

          if (iturbkern.ne.0) then
             Dc = ( Xc / prw )**(1./3.)  ! mass mean diameter cloud droplets in m
             !
             ! Calculate the mixing length, dissipation rate and Taylor-Reynolds number
             ! (Axel: I think this is wrong, the equation can not be used for LES flow)
             !epsilon = diss(k)
             !l = csx*((1/dxi)*(1/dyi)*(1/dzi_t(k)))**(1./3.)
             !Re = (6./11.)*((l/ce)**(2./3))*((15./(1.5e-5))**0.5)*(epsilon**(1./6) )
             epsilon = min(diss(k)*1.e4,1e3)
             Re = max(10000. * (epsilon/100.)**(1./6.), 2000.) 
             !
             ! Dissipation needs to be converted to cm^2/s^3, which explains the factor
             ! 1.e4 with which epsilon is multiplied
             !
             k_au = k_au0 * (1. + epsilon * Re**prey &
                  * (kc_bet + kc_alf * exp( -1.* ((((Dc/2.)*1.e6-kc_rad)/kc_sig)**2) )))
             !
             if (present(emax)) emax = max(diss(k),emax)
          else
             k_au = k_au0
          end if

          ! which one is correct here?
          !au = k_au * rc(k)**2 * Xc**2
          ! neglecting density correction:
          au = k_au * dn0(k) * rc(k)**2 * Xc**2
          !
          ! small threshold that should not influence the result
          !
          if (rc(k) > 1.e-6) then
             tau = 1.0-rc(k)/(rc(k)+rp(k)+eps0)
             tau = MIN(MAX(tau,eps0),0.9)
             phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
             au  = au * (1.0 + phi/(1.0 - tau)**2)
          endif
          !
          ! Khairoutdinov and Kogan
          !
          if (khairoutdinov_au) then
             Dc = ( Xc / prw )**(1./3.)
             au = Cau * (Dc * mmt / 2.)**Eau
          end if
          if (kessler) then
             if (rc(k) > rc0) then
                au = alpha * (rc(k) -rc0)
             else
                au = 0.
             end if
          end if
          au    = au * dt
          au    = min(au,rc(k))
          rp(k) = rp(k) + au
          rc(k) = rc(k) - au
          tl(k) = tl(k) + convliq(k)*au
          if (.not.mom3) then
            !np(k) = np(k) + au/cldw%x_max    !old autoconversion, not consistent with mom3
            np(k) = np(k) + 2./3. / cldw%x_max*au
          else
            np(k) = np(k) + 2./3. / cldw%x_max*au
            zp(k) = zp(k) +14./9. * cldw%x_max*au
            !np(k) = np(k) + 4./5.  / cldw%x_max*au
            !zp(k) = zp(k) +19./15. * cldw%x_max*au
            !np(k) = np(k) + 8./9.  / cldw%x_max*au
            !zp(k) = zp(k) + 1.12963 * cldw%x_max*au
          end if
          
          if (present(au_acc)) au_acc = au_acc + au * (zm(k)-zm(k-1))*dn0(k) / 3.

          ! For particles: a_npauto in #/(kg*dt)
          if(lpartdrop .and. nstep==1) a_npauto(k,i,j) = a_npauto(k,i,j) + (rkalpha(1) + rkalpha(2))* au*dn0(k)
          if(lpartdrop .and. nstep==3) a_npauto(k,i,j) = a_npauto(k,i,j) + (rkalpha(3)             )* au*dn0(k)
          !if(lpartdrop .and. nstep==3) a_npauto(k,i,j) = au*dn0(k)
       end if
    end do

    !if (isIO().and.MAXVAL(rc).gt.0) THEN
    !if (isIO()) THEN
    !  WRITE(0,'(a,15e15.4)') '  k_au = ',k_au0
    !end if

  end subroutine auto_SB

  subroutine accr_SB(n1,dn0,rc,rp,np,tl,diss,zp,ac_acc,zm)
    !
    ! ---------------------------------------------------------------------
    ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
    ! and self collection following Seifert & Beheng (2001).  Included is
    ! an alternative formulation for accretion only, following
    ! Khairoutdinov and Kogan
    !

    integer, intent (in) :: n1
    real, intent (inout)    :: rc(n1), rp(n1), np(n1),tl(n1)
    real, intent (in)    :: dn0(n1),diss(n1)
    real, intent (inout), optional :: zp(n1)
    real, intent (inout), optional :: ac_acc
    real, intent (in),    optional :: zm(n1)

    real, parameter :: k_r0 = 5.78
    real, parameter :: k_rr = 4.33
    real, parameter :: k_1 = 5.e-4
    real, parameter :: Cac = 67.     ! accretion coefficient in KK param.
    real, parameter :: Eac = 1.15    ! accretion exponent in KK param.
    real, parameter :: gam = 6e2     ! [1/m]
    real, parameter :: bet = 9.65     ! [m/s]
    logical, parameter :: longkernel = .false.
    logical, parameter :: breakup = .true.

    integer :: k
    real    :: tau, phi, ac, sc, k_r, epsilon, k_turb, xp
    real    :: zac, zsc, mue, lam, x_r, d_r, cscn, cscz

!    if (present(ac_acc)) ac_acc = 0.0

    do k=2,n1-1
       if (rp(k) > 0.) then

          k_r = k_r0

          ! Simulate the effect of turbulence on the collision kernel
          ! (dissipation needs to be converted to cm^2/s^3)
          if (iturbkern.ne.0) then
             epsilon = min(diss(k)*1.e4,1e3)   ! 1e4 for CGS
             if (iturbkern.eq.2.or.iturbkern.eq.22) then
                k_turb = 1. + (0.05*epsilon**0.25)  ! Ayala-Wang kernel
             elseif (iturbkern.eq.1) then
                xp = rp(k) / np(k)
                xp = MIN(MAX(xp,rain%x_min),rain%x_max)  ! Onishi kernel
                k_turb = 1. + 0.8e-3 * epsilon * (rain%x_min/xp)**(2./3.)
             end if
          !!!   k_turb = 1.0
          else
             k_turb = 1.0
          end if
          if (rc(k) > 0.) then
             ! accretion
             tau = 1.0-rc(k)/(rc(k)+rp(k)+eps0)
             tau = MIN(MAX(tau,eps0),1.)
             phi = (tau/(tau+k_1))**4

             ! including density correction:
             ac  = k_r * rc(k) * rp(k) * phi * frho(k) * k_turb
             if (mom3) then
               zac = 2.0 * k_r * rc(k) * zp(k) * phi * frho(k) * k_turb
             end if
             !
             ! Khairoutdinov and Kogan
             !
             if (khairoutdinov) then
                ac = Cac * (rc(k) * rp(k))**Eac
             end if
             !
             ac    = ac * dt
             ac    = min(ac, rc(k))
             rp(k) = rp(k) + ac
             rc(k) = rc(k) - ac
             tl(k) = tl(k) + convliq(k)*ac
             if (mom3) then
               zac   = zac * dt
               zp(k) = zp(k) + zac
             end if
             if (present(ac_acc)) ac_acc = ac_acc + ac * (zm(k)-zm(k-1))*dn0(k) / 3.
          end if

          !selfcollection
          if (.true.) then 
            if (longkernel) then        ! using Long's kernel
              sc = k_rr * np(k) * rp(k) * frho(k) *dt
              zsc  = 2.* k_rr * zp(k) *np(k)* frho(k) * dt
            else                        ! using variance approximation

              x_r = rp(k)/(np(k)+eps0) 
              d_r = ( x_r / prw )**(1./3.)
              if (mom3) then
                 mue = rain_mue_z_inv(np(k), rp(k), zp(k))
              else
                 mue = rain_mue_2mom(d_r)
              end if
              
              cscn = 0.5 * (mue+1.)*((mue+1.)+(mue+2.)) 
              lam  = (prw*(mue+3.)*(mue+2.)*(mue+1.)/x_r)**(1./3.)
                    
              sc = pi/sqrt(8.0)*cscn*bet*np(k)*np(k)/lam**2 * frho(k) & 
                         *( 1./(1. + 2.*gam/lam)**(1.0*(mue+3.0)) & 
                          - 1./(1. + 1.*gam/lam)**(2.0*(mue+3.0)) )**0.5 * dt
              if (mom3) then
                cscz = 0.5 * (mue+4.)*(mue+3.)*(mue+2.)*(mue+1.)/(mue+6.)/(mue+5.) * (1.0+(mue+5.)/(mue+4.))
                zsc = pi/sqrt(2.0)*cscz*bet*np(k)*zp(k)/lam**2 * frho(k) & 
                         *( 1./(1. + 2.*gam/lam)**(1.0*(mue+6.0)) & 
                          - 1./(1. + 1.*gam/lam)**(2.0*(mue+6.0)) )**0.5 & 
                         * min(max(0.53*(1.-0.69e3*d_r),0.0),0.5) * dt
              end if
            end if
          end if

          if (breakup) then
            x_r = rp(k)/(np(k)+eps0) 
            d_r = ( x_r / prw )**(1./3.) * 1e3  ! D in mm
            if (d_r > 0.3) then 
               sc = sc * ( 1.0 - (1.0*(d_r-1.1) + 1.0))
            endif
          endif
          
          sc = min(sc, np(k))
          np(k) = np(k) - sc
          if (mom3) then
            zp(k) = zp(k) + zsc
          end if

       end if
    end do

  end subroutine accr_SB

  subroutine sedim_rd(n1,dt,dn0,rp,np,rrate,zp)
    !
    ! ---------------------------------------------------------------------
    ! SEDIM_RD: calculates the sedimentation of the rain drops and its
    ! effect on the evolution of theta_l and r_t.  This is expressed in
    ! terms of Dp the mean diameter, not the mass weighted mean diameter
    ! as is used elsewhere.  This is just 1/lambda in the exponential
    ! distribution
    !

    integer, intent (in)                      :: n1
    real, intent (in)                         :: dt
    real, intent (in),    dimension(n1)       :: dn0
    real, intent (inout), dimension(n1) :: rp,np
    real, intent (out)  , dimension(n1) :: rrate
    real, intent (inout), optional :: zp(n1)
    real, parameter :: a2 = 9.65       ! in SI [m/s]
    real, parameter :: c2 = 6e2        ! in SI [1/m]
    real, parameter :: Dv = 25.0e-6    ! in SI [m/s]
    real, parameter :: aq = 6.0e3
    real, parameter :: bq = -0.2
    real, parameter :: an = 3.5e3
    real, parameter :: bn = -0.1


    integer ::  k, kp1, kk, km1
    real    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz, cmax
    real, dimension(n1) :: nslope,rslope,dn,dr, rfl, nfl, vn, vr, cn, cr
    real, dimension(n1) :: zslope, dzz, zfl, vz, cz

    logical, parameter :: oldsedimentation = .false.

    b2 = a2*exp(c2*Dv)

    nfl(n1) = 0.
    rfl(n1) = 0.
    vn = 0.
    vr = 0.
    do k=n1-1,2,-1
      if (rp(k) > rthres) then
        Xp = rp(k) / np(k)
        xp = MIN(MAX(xp,rain%x_min),rain%x_max)
     
        !
        ! Adjust Dm and mu-Dm and Dp=1/lambda following Milbrandt & Yau
        !
        Dm = ( 6. / (rowt*pi) * Xp )**(1./3.)
        if (mom3) then
          mu = rain_mue_z_inv(np(k),rp(k),zp(k))
        else
          if (oldsedimentation) then
            mu = cmur1*(1.+tanh(cmur2*(Dm-cmur3)))    !MY05 revised constants
          else
            mu = rain_mue_2mom(Dm) 
          end if
        end if
        Dp = (Dm**3/((mu+3.)*(mu+2.)*(mu+1.)))**(1./3.)
        
        ! including density corrcetion
        vn(k) = frho(k)/dn0(k) *(a2 - b2*(1.+c2*Dp)**(-(1.+mu)))
        vr(k) = frho(k)/dn0(k) *(a2 - b2*(1.+c2*Dp)**(-(4.+mu)))
        if (mom3) then
          vz(k) = frho(k)/dn0(k) *(a2 - b2*(1.+c2*Dp)**(-(7.+mu)))
        end if
        !
        ! Set fall speeds following Khairoutdinov and Kogan

        !              if (khairoutdinov) then
        !                 vn(k) = max(0.,an * Dp + bn)
        !                 vr(k) = max(0.,aq * Dp + bq)
        !              end if
        !irina-olivier
        if (khairoutdinov) then
            vn(k) = max(0.,an * Dp + bn)
            vr(k) = max(0.,aq * Dp + bq)
        end if
        vn(k) = min(max(vn(k),0.1),20.)
        vr(k) = min(max(vr(k),0.1),20.)
        if (mom3) then
          vz(k) = min(max(vz(k),0.1),20.)
        end if
      end if
      
    end do

    do k=2,n1-1
       kp1 = min(k+1,n1-1)
       km1 = max(k,2)
       cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzi_t(k)*dt
       cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzi_t(k)*dt
       if (mom3) then
         cz(k) = 0.25*(vz(kp1)+2.*vz(k)+vz(km1))*dzi_t(k)*dt
       end if
    end do

    !cmax = MAXVAL(cr)
    !if (cmax.gt.1.0) then
    !  WRITE(0,'(a,i5,a,f8.2)') 'sedim: myid = ',myid,' meteor = rain,  cmax = ',cmax
    !end if

    !...piecewise linear method: get slopes
    do k=n1-1,2,-1
       dn(k) = np(k+1)-np(k)
       dr(k) = rp(k+1)-rp(k)
       if (mom3) then
         dzz(k) = zp(k+1)-zp(k)
       end if
    enddo
    dn(1)  = dn(2)
    dn(n1) = dn(n1-1)
    dr(1)  = dr(2)
    dr(n1) = dr(n1-1)
    if (mom3) then
      dzz(1) = dzz(2)
      dzz(n1) = dzz(n1-1)
    end if
    do k=n1-1,2,-1
       !...slope with monotone limiter for np
       sk = 0.5 * (dn(k-1) + dn(k))
       mini = min(np(k-1),np(k),np(k+1))
       maxi = max(np(k-1),np(k),np(k+1))
       nslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(np(k)-mini), &
            &                                     2.*(maxi-np(k)))
       !...slope with monotone limiter for rp
       sk = 0.5 * (dr(k-1) + dr(k))
       mini = min(rp(k-1),rp(k),rp(k+1))
       maxi = max(rp(k-1),rp(k),rp(k+1))
       rslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(rp(k)-mini), &
            &                                     2.*(maxi-rp(k)))
       if (mom3) then
         !...slope with monotone limiter for zp
         sk = 0.5 * (dzz(k-1) + dzz(k))
         mini = min(zp(k-1),zp(k),zp(k+1))
         maxi = max(zp(k-1),zp(k),zp(k+1))
         zslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(zp(k)-mini), &
            &                                       2.*(maxi-zp(k)))
       end if
    enddo

    rfl = 0.
    nfl = 0.
    zfl = 0.
    do k=n1-2,2,-1

       kk = k
       tot = 0.0
       zz  = 0.0
       cc  = min(1.,cn(k))
       do while (cc > 0 .and. kk <= n1-1)
          tot = tot + dn0(kk)*(np(kk)+nslope(kk)*(1.-cc))*cc/dzi_t(kk)
          zz  = zz + 1./dzi_t(kk)
          kk  = kk + 1
          cc  = min(1.,cn(kk) - zz*dzi_t(kk))
       enddo
       tot = min(tot,dn0(k)/dzi_t(k) * np(k) - nfl(k+1) * dt - rthres)  
       nfl(k) = -tot /dt

       kk = k
       tot = 0.0
       zz  = 0.0
       cc  = min(1.,cr(k))
       do while (cc > 0 .and. kk <= n1-1)
          tot = tot + dn0(kk)*(rp(kk)+rslope(kk)*(1.-cc))*cc/dzi_t(kk)
          zz  = zz + 1./dzi_t(kk)
          kk  = kk + 1
          cc  = min(1.,cr(kk) - zz*dzi_t(kk))
       enddo
       tot = min(tot,dn0(k)/dzi_t(k) * rp(k) - rfl(k+1) * dt - rthres)
       rfl(k) = -tot /dt

       if (mom3) then
         kk = k
         tot = 0.0
         zz  = 0.0
         cc  = min(1.,cz(k))
         do while (cc > 0 .and. kk <= n1-1)
            tot = tot + dn0(kk)*(zp(kk)+zslope(kk)*(1.-cc))*cc/dzi_t(kk)
            zz  = zz + 1./dzi_t(kk)
            kk  = kk + 1
            cc  = min(1.,cz(kk) - zz*dzi_t(kk))
         enddo
         tot = min(tot,dn0(k)/dzi_t(k) * zp(k) - zfl(k+1) * dt - rthres)
         zfl(k) = -tot /dt
       end if

       kp1=k+1
       flxdiv = (rfl(kp1)-rfl(k))*dzi_t(k)/dn0(k)
       rp(k) = rp(k)-flxdiv*dt
       np(k) = np(k)-(nfl(kp1)-nfl(k))*dzi_t(k)/dn0(k)*dt
       if (mom3) then
         zp(k) = zp(k)-(zfl(kp1)-zfl(k))*dzi_t(k)/dn0(k)*dt
       end if

       rrate(k)    = -rfl(k)

    end do

  end subroutine sedim_rd

  subroutine sedim_cd(n1,dt,tl,rc,rrate)
    !
    ! ---------------------------------------------------------------------
    ! SEDIM_CD: calculates the cloue-droplet sedimentation flux and its effect
    ! on the evolution of r_t and theta_l assuming a log-normal distribution
    !

    integer, intent (in):: n1
    real, intent (in)                        :: dt
    !irina
    !     real, intent (in),   dimension(n1) :: th,tk
    !irina
    !
    real, intent (inout),dimension(n1) :: rc,tl
    real, intent (out)  , dimension(n1) :: rrate

    real, parameter :: c = 1.19e8 ! Stokes fall velocity coef [m^-1 s^-1]
    real, parameter :: sgg = 1.2  ! geometric standard dev of cloud droplets

    integer ::  k, kp1
    real    :: Dc, Xc, vc, flxdiv
    real    :: rfl(n1)


    !
    ! calculate the precipitation flux and its effect on r_t and theta_l
    !
    rfl = 0.
    do k=n1-1,2,-1
      if (rc(k) > 0.) then
        Xc = rc(k) / cldw%nr
        Dc = ( Xc / prw )**(1./3.)
        vc = min(c*(Dc*0.5)**2 * exp(4.5*(log(sgg))**2),1./(dzi_t(k)*dt))
        !               vc = min(c*(Dc*0.5)**2 * exp(5*(log(1.3))**2),1./(dzi_t(k)*dt))
        rfl(k) = - rc(k) * vc
        !
      end if
      kp1=k+1
      flxdiv = (rfl(kp1)-rfl(k))*dzi_t(k)
      rc(k) = rc(k)-flxdiv*dt
      tl(k) = tl(k)+flxdiv*convliq(k)*dt
      rrate(k)    = -rfl(k)
    end do

  end subroutine sedim_cd

  real function rain_mue_2mom(Dp) result(mue)
    real, intent(in)   :: Dp

    if (mue_SB) then
       IF (Dp.LE.rain_cmu3) THEN ! see Seifert (2008)            
          mue = rain_cmu0*TANH((4.*rain_cmu2*(Dp-rain_cmu3))**2) &
               & + rain_cmu4
       ELSE
          mue = rain_cmu1*TANH((1.*rain_cmu2*(Dp-rain_cmu3))**2) &
                         & + rain_cmu4
       ENDIF
    else
       mue=cmu_mc0*(cmu_mc1*Dp-cmu_mc2)**2+cmu_mc3  !Milbrandt and McTaggart-Cowan
       !mue = cmur1*(1.+tanh(cmur2*(Dp-cmur3)))     !Bjorn's revised constants
       !mue = cmy1*tanh(cmy2*(Dp-cmy3))+cmy4        !MY05 original constants
       !mue = 3.
    end if
    
  end function rain_mue_2mom

  real function rain_mue_z_inv(N,L,Z) 

    real, intent(in)   :: N,L,Z
    real               :: z1
    real               :: gg,gg1,aa0,aa1,aa2,qq,rr,dd,ss,tt
    real, dimension(4), parameter :: cc = (/0.5569344,0.03052414,-1.078084,-0.000936101/)
    logical, parameter :: Cardano = .true.

    real, parameter :: gg20 = 1.4681 ! gg(mue=20)
    real, parameter :: gg15 = 1.6299 ! gg(mue=15)
    real, parameter :: gg01 = 8.75   ! gg(mue=1)
    real, parameter :: gg00 = 20.0   ! gg(mue=0)

    ! gg(mue) = (mue+6.)*(mue+5.)*(mue+4.)/((mue+3.)*(mue+2.)*(mue+1.))

    if (L.gt.1e-10) then
       gg  = N*Z/(L*L)        ! dimensionless parameter, depends only on mue         
       if (Cardano) then
          ! inversion of gg(mue) using Cardano's formula (cubic root)
          ! (see e.g. http://mathworld.wolfram.com/CubicFormula.html)          
          !if (gg.gt.gg00) then
          !   z1 = 0.0
          if (gg.gt.gg01) then
             z1 = 1.0
          elseif (gg.gt.gg20) then 
             gg1 = 1.0/(gg-1.0)
             aa0 = (6.0*gg-120.)*gg1
             aa1 = (11.*gg-74.0)*gg1
             aa2 = (6.0*gg-15.0)*gg1
             qq = (3.0*aa1-aa2**2)/9.
             rr = (9.0*aa1*aa2-27.*aa0-2.*aa2**3)/54.0
             dd = qq*qq*qq+rr*rr
             if (dd.gt.0) then
                ss = (rr+sqrt(dd))**(1./3.)
                tt = (rr-sqrt(dd))**(1./3.)         
                z1 = -aa2/3. + (ss+tt)
             else
                print *, 'dd < 0 in rain_mue_z_inv, gg = ',gg
                z1 = 1.0
             end if
          else
             z1 = 20.0
          end if
       else
          ! rational function fit (not much faster than Cardano's formula)
          if (gg.gt.gg20) then
             z1 = (gg-20.)*(cc(1)+cc(2)*gg+cc(4)*gg*gg)/(1+cc(3)*gg)
          else
             z1 = 20.
          end if
       end if
    else
       z1 = 5.0
    end if      
    rain_mue_z_inv = max(z1,0.0)
         
  end function rain_mue_z_inv

  real elemental function gfct(x)
    !*******************************************************************************
    !                                                                              *
    !       gammafunction from numerical recipes (f77)                              *
    !                                                                              *
    !*******************************************************************************
    real, intent(in) :: x
    real, parameter :: cof(6) = (/76.18009173e0,-86.50532033e0,24.01409822e0, -1.231739516e0,.120858003e-2,-.536382e-5/)
    real, parameter :: stp    = 2.50662827465e0
    real :: xx,tmp,ser,gamma
    integer j

    xx  = x  - 1.0
    tmp = xx + 5.5
    tmp = (xx + 0.5) * log(tmp) - tmp
    ser = 1.0
    do j = 1,6
       xx  = xx  + 1.0
       ser = ser + cof(j) / xx
    enddo
    gamma = tmp + log(stp*ser)
    gamma = exp(gamma)

    gfct = gamma
  end function gfct

  real elemental function incgfct_upper(a,x)
    !*******************************************************************************
    !
    ! upper incomplete gamma function
    !
    !              int(x)(oo) exp(-t) t^(a-1) dt
    !
    !*******************************************************************************

    real, intent(in) :: a, x
    real :: gam, gln

    call gammq(a,x,gam,gln)
    incgfct_upper = exp(gln) * gam

  end function incgfct_upper


  real elemental function incgfct_lower(a,x)
    !*******************************************************************************
    !
    ! lower incomplete gamma function
    !

    !              int(0)(x) exp(-t) t^(a-1) dt
    !
    !   !*******************************************************************************

    real, intent(in) :: a, x
    real :: gam, gln

    call gammp(a,x,gam,gln)
    incgfct_lower = exp(gln) * gam

  end function incgfct_lower


  real elemental function incgfct(a,x1,x2)
    !*******************************************************************************
    !
    ! actual incomplete gamma-function
    !
    !*******************************************************************************

    real, intent(in) :: a, x1, x2

    incgfct = incgfct_lower(a,x2) - incgfct_lower(a,x1)

  end function incgfct

  elemental subroutine gammp(a,x,gamma,gln)


    real, intent(in) :: a, x
    real, intent(out) :: gamma, gln

    real :: gammcf, gamser

    if (x.lt.0.0e0 .or. a .le. 0.0e0) then
       gamma = 0.0e0
       return
    end if

    if (x .lt. a+1.)then
       call gser(gamser,a,x,gln)
       gamma = gamser
    else
       call gcf(gammcf,a,x,gln)
       gamma = 1.0e0 - gammcf
    endif
  end subroutine gammp

  elemental subroutine gammq(a,x,gamma,gln)

    real, intent(in) :: a, x
    real, intent(out) :: gamma,gln
    real :: gammcf, gamser

    if (x.lt.0.0e0 .or. a .le. 0.0e0) then
!        write(*,*) 'error in gammq: bad arguments (module gamma_functions)'
       gamma = 0.0e0
       return
    end if

    if (x.lt.a+1.) then
       call gser(gamser,a,x,gln)
       gamma = 1.0e0 - gamser
    else
       call gcf(gammcf,a,x,gln)
       gamma = gammcf
    endif
  end subroutine gammq

  real elemental function gammln(x)
    !*******************************************************************************
    !                                                                              *
    !       log(gamma function) taken from press et al.,  numerical recipes (f77)
    !
    !*******************************************************************************

    real, intent(in) :: x
    real, parameter :: cof(6) = (/76.18009172947146e0,-86.50532032941677e0, &
         24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2, -.5395239384953e-5/)
    real, parameter :: stp = 2.5066282746310005e0
    real :: xx,tmp,ser
    integer :: j

    xx  = x
    tmp = xx + 5.5e0
    tmp = (xx + 0.5e0) * log(tmp) - tmp
    ser = 1.000000000190015e0
    do j = 1,6
       xx  = xx  + 1.0e0
       ser = ser + cof(j) / xx
    enddo
    gammln = tmp + log(stp*ser/x)
  end function gammln

  real elemental function gfct2(x)
    !*******************************************************************************
    !                                                                              *
    !       gamma function taken from press et al.,  numerical recipes (f77)
    !
    !       (other formulation, same results)
    !*******************************************************************************
    real, intent(in) :: x
    real, parameter  :: cof(6) =(/76.18009173e0,-86.50532033e0,24.01409822e0,  &
         &     -1.231739516e0,.120858003e-2,-.536382e-5/)
    real, parameter  :: stp = 2.50662827465e0
    real :: xx,tmp,ser,gamma
    integer j

    xx  = x  - 1.
    tmp = xx + 5.5
    tmp = (xx + 0.5) * log(tmp) - tmp
    ser = 1.
    do j = 1,6
       xx  = xx  + 1.
       ser = ser + cof(j) / xx
    enddo
    gamma = tmp + log(stp*ser)
    gamma = exp(gamma)

    gfct2 = gamma
  end function gfct2

  elemental subroutine gcf(gammcf,a,x,gln)
    integer, parameter :: itmax = 100
    real, parameter :: eps = 3.e-7, fpmin = 1.e-30
    real, intent(in) :: a, x
    real, intent(out) :: gammcf, gln

    integer :: i
    real :: an,b,c,d,del,h

    gln=gammln(a)
    b=x+1.-a
    c=1./fpmin
    d=1./b
    h=d
    do i=1,itmax
       an=-i*(i-a)
       b=b+2.0e0
       d=an*d+b
       if (abs(d).lt.fpmin) d=fpmin
       c=b+an/c
       if (abs(c).lt.fpmin) c=fpmin
       d=1./d
       del=d*c
       h=h*del
       if (abs(del-1.).lt.eps) exit
    end do

    if (abs(del-1.).ge.eps) then
!        write (*,*) 'error in gcf: a too large, itmax too small (module gamma_functions)'
       gammcf = 0.0e0
       return
    end if

    gammcf=exp(-x+a*log(x)-gln)*h
  end subroutine gcf

  elemental subroutine gser(gamser,a,x,gln)
    integer, parameter :: itmax = 100
    real, parameter :: eps=3.e-7
    real, intent(in) :: a, x
    real, intent(out) :: gamser, gln

    integer :: n
    real :: ap,del,sum

    gln=gammln(a)
    if (x.le.0.) then
       if (x.lt.0.) then
!           write (*,*) 'error in gser: x < 0 (module gamma_functions)'
       end if
       gamser=0.0e0
       return
    endif

    ap=a
    sum=1./a
    del=sum
    do n=1,itmax
       ap=ap+1.
       del=del*x/ap
       sum=sum+del
       if (abs(del).lt.abs(sum)*eps) exit
    end do

    if (abs(del).ge.abs(sum)*eps) then
!        write (*,*) 'error in gser: a too large, itmax too small'
!        write (*,*) '  (module gamma_functions)'
       gamser = 0.0e0
       return
    end if

    gamser = sum*exp(-x+a*log(x)-gln)

  end subroutine gser

  subroutine resetvar(meteor,mass,num,refl)
    type(particle),intent(in)        :: meteor
    real, dimension(:), intent(inout) :: mass
    real, dimension(:), intent(inout), optional :: num
    real, dimension(:), intent(inout), optional :: refl

 !   if (any(mass < 0.)) then
 !     print *, trim(meteor%name), 'below zero', mass
 !   end if
    where (mass < rthres)
       mass = 0.
    end where
    mass(1) = 0.
    if (present(num)) then
       where (meteor%x_max*num < mass)
          num = mass/meteor%x_max
       end where

       where (meteor%x_min*num > mass)
          num = mass/meteor%x_min
       end where
    end if
    if (present(refl)) then !restrict mu to [0,20]
      where (refl< 1.4681*mass**2/num)
        refl = 1.4681*mass**2/num
      end where
      where (refl> 20.*mass**2/num)
        refl = 20.*mass**2/num
      end where
    end if
  end subroutine resetvar

  subroutine initmcrp(level,firsttime)
    integer , intent(in) :: level
    logical, optional, intent(inout) :: firsttime

    if (present(firsttime)) then
       firsttime = .false.
    end if

    if (level<3)  mom3 = .false.

    if (level==5) then
       CALL init_sb
    else
       !Cloudwater properties
       cldw = PARTICLE('cloudw', &! HN: made variable for seeding experiments
            1, & !number of moments
            CCN, & !Number of droplets
            1.000000, & !.nu.....Width parameter of the distribution
            1.000000, & !.mu.....exponential parameter of the distribution
            2.60e-10, & !.x_max..maximum particle mass D=80e-6m
            4.20e-15, & !.x_min..minimale particler mass D=2.e-6m
            1.24e-01, & !.a_geo..coefficient of meteor geometry
            0.333333, & !.b_geo..coefficient of meteor geometry = 1/3
            3.75e+05, & !.a_vel..coefficient of fall velocity
            0.666667, & !.b_vel..coefficient of fall velocity
            0.25     , &!.s_vel...dispersion of the fall velocity
            0.780000, & !.a_ven..ventilation coefficient (PK, p.541)
            0.308000, & !.b_ven..ventilation coefficient (PK, p.541)
            2.0)        !.cap....capacity coefficient

       !Rainwater properties
       rain = PARTICLE('rain', &! HN: made variable for seeding experiments
            2, & !number of moments
            0, & !Number of droplets
            1.000000, & !.nu.....Width parameter of the distribution
            0.333333, & !.mu.....exponential parameter of the distribution
            3.00e-06, & !.x_max..maximum particle mass
            2.60e-10, & !.x_min..minimale particler mass D=80.e-6m
            1.24e-01, & !.a_geo..coefficient of meteor geometry
            0.333333, & !.b_geo..coefficient of meteor geometry = 1/3
            1.59e+02, & !.a_vel..coefficient of fall velocity
            0.266667, & !.b_vel..coefficient of fall velocity
            0.25    , & !.s_vel...dispersion of the fall velocity
            0.780000, & !.a_ven..ventilation coefficient (PK, p.541)
            0.308000, & !.b_ven..ventilation coefficient (PK, p.541)
            2.0)        !.cap....capacity coefficient

    end if

  end subroutine initmcrp

  subroutine sedimentation(n1,dn0,meteor,rp,nr,rrate,dtopt)
    integer, intent(in) :: n1
    real, dimension(n1), intent(inout) :: rp
    type(particle),  intent(in) :: meteor
    real, dimension(n1), intent(inout),optional :: nr
    real, dimension(n1), intent(in),optional :: dn0
    real, intent(in), optional :: dtopt
    real, intent (out)  , dimension(n1) :: rrate

    integer  :: k,kk,km1,kp1
    real,dimension(4),save      :: alfn,alfq,c_lam
    real :: lam,sk,tot,zz,xp
    real, dimension(n1)   :: nfl,rfl,vn,vr,dn,dr,rslope,cn,cr,nslope,np
    real :: cc, flxdiv,maxi,mini,vlimit,vmin,cmax,dtsedi
    logical, save :: firsttime(4) = .true.
    integer :: metnr = 0
    logical :: plm = .false.

    if(present(dtopt)) then
       dtsedi = dtopt
    else
       dtsedi = dt
    end if

    select case (meteor%name)
    case('ice')
       metnr = 1
       vmin   = 0.1
       vlimit = 2.0
    case('snow')
       metnr = 2
       vmin   = 0.1
       vlimit = 2.0
    case('graupel')
       metnr = 3
       vmin   = 0.1
       vlimit = 30.0
    case('hail')
       metnr = 4
       vmin   = 0.1
       vlimit = 30.0
    case default
       WRITE (0,*) 'icrmcrp: stopped in sedimentation ',meteor%name
       stop
    end select
    if (present(nr)) then
       np = nr
    else
       np = meteor%nr
    end if

    if (firsttime(metnr)) then
       alfn(metnr) = meteor%a_vel * gfct((meteor%nu+meteor%b_vel+1.0)/meteor%mu)&
            &                / gfct((meteor%nu+1.0)/meteor%mu)
       alfq(metnr) = meteor%a_vel * gfct((meteor%nu+meteor%b_vel+2.0)/meteor%mu)&
            &                / gfct((meteor%nu+2.0)/meteor%mu)
       c_lam(metnr) = gfct((meteor%nu+1.0)/meteor%mu)/gfct((meteor%nu+2.0)/meteor%mu)
       WRITE(0,'(a,i5,3a,e11.3)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  alfn = ',alfn(metnr)
       WRITE(0,'(a,i5,3a,e11.3)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  alfq = ',alfq(metnr)
       WRITE(0,'(a,i5,3a,e11.3)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  clam = ',c_lam(metnr)
       firsttime(metnr) =.false.
    end if

    do k=n1-1,2,-1
      if ((rp(k) > rthres).and.(np(k).ne.0.)) then
       Xp = rp(k) / np(k)
       xp = MIN(MAX(xp,meteor%x_min),meteor%x_max)
       lam = ( c_lam(metnr) * xp )**(meteor%b_vel)
       vr(k) = alfq(metnr) * lam
       vr(k) = max(vr(k),vmin)
       vr(k) = min(vr(k),vlimit)
       !      vr(k) = -vr(k)

       vn(k) = alfn(metnr) * lam * (dn0(1)/dn0(k))**0.5
       vn(k) = max(vn(k),vmin)
       vn(k) = min(vn(k),vlimit)
       !      vn(k) = -vn(k)
      else
        vr(k) = 0.
        vn(k) = 0.
      end if
    end do
    do k=2,n1-1
       kp1 = min(k+1,n1-1)
       km1 = max(k,2)
       cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzi_t(k)*dtsedi
       cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzi_t(k)*dtsedi
    end do

    !cmax = MAXVAL(cr)
    !if (cmax.gt.0.8) then
    !  WRITE(0,'(a,i5,3a,f8.2)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  cmax = ',cmax
    !end if

    if (plm) then
    !...piecewise linear method: get slopes
    do k=n1-1,2,-1
       dn(k) = np(k+1)-np(k)
       dr(k) = rp(k+1)-rp(k)
    enddo
    dn(1)  = dn(2)
    dn(n1) = dn(n1-1)
    dr(1)  = dr(2)
    dr(n1) = dr(n1-1)
    do k=n1-1,2,-1
       !...slope with monotone limiter for np
       sk = 0.5 * (dn(k-1) + dn(k))
       mini = min(np(k-1),np(k),np(k+1))
       maxi = max(np(k-1),np(k),np(k+1))
       nslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(np(k)-mini), &
            &                                     2.*(maxi-np(k)))
       !...slope with monotone limiter for rp
       sk = 0.5 * (dr(k-1) + dr(k))
       mini = min(rp(k-1),rp(k),rp(k+1))
       maxi = max(rp(k-1),rp(k),rp(k+1))
       rslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(rp(k)-mini), &
            &                                     2.*(maxi-rp(k)))
    enddo
    else
       nslope = 0.0
       rslope = 0.0
    end if

    rfl(n1-1) = 0.
    nfl(n1-1) = 0.
    do k=n1-2,2,-1

       kk = k
       tot = 0.0
       zz  = 0.0
       cc  = min(1.,cn(k))
       do while (cc > 0 .and. kk <= n1-1)
          tot = tot + dn0(kk)*(np(kk)+nslope(kk)*(1.-cc))*cc/dzi_t(kk)
          zz  = zz + 1./dzi_t(kk)
          kk  = kk + 1
          cc  = min(1.,cn(kk) - zz*dzi_t(kk))
       enddo
       tot = min(tot,dn0(k)/dzi_t(k) * np(k) - nfl(k+1) * dtsedi - rthres)  
       nfl(k) = -tot /dtsedi

       kk = k
       tot = 0.0
       zz  = 0.0
       cc  = min(1.,cr(k))

       do while (cc > 0 .and. kk <= n1-1)
          tot = tot + dn0(kk)*(rp(kk)+rslope(kk)*(1.-cc))*cc/dzi_t(kk)
          zz  = zz + 1./dzi_t(kk)
          kk  = kk + 1
          cc  = min(1.,cr(kk) - zz*dzi_t(kk))
       enddo
       tot = min(tot,dn0(k)/dzi_t(k) * rp(k) - rfl(k+1) * dtsedi - rthres)
       rfl(k) = -tot /dtsedi

       kp1=k+1
       flxdiv = (rfl(kp1)-rfl(k))*dzi_t(k)/dn0(k)
       rp(k) = rp(k)-flxdiv*dtsedi
       if(meteor%moments==2) then
          np(k) = np(k)-(nfl(kp1)-nfl(k))*dzi_t(k)/dn0(k)*dtsedi
       end if
       rrate(k)    = -rfl(k)
    end do


  end subroutine sedimentation



  !==============================================================================
  !
  ! Two-moment mixed-phase bulk microphysics
  !
  ! original version by Axel Seifert, May 2003
  ! with modifications by Ulrich Blahak, August 2007
  !
  ! Description:
  !
  ! The subroutine is the interface to the original KAMM2 modules,
  ! which are provided by src_seifert.f90. A major difference between UCLA-LES and
  ! KAMM2 is that KAMM2 uses mass densities instead of mixing ratios, thus
  ! q_cloud=rho*qc.
  ! Temporary allocation of memory to the KAMM2 variables is done by the
  ! subroutines ALLOC_DRIVER and ALLOC_WOLKEN.
  ! All microphysical source terms e.g. nucleation, condensation, coagulation,
  ! freezing and melting are then calculated and time integrated within the
  ! subroutine CLOUDS.
  !
  !==============================================================================

  SUBROUTINE mcrph_sb(level,ke,je,ie,dn0,exner,pi0,pi1,thl,tlt,tk,w,qvin,qcin, &
       qrin, qnrin, qiin, qniin, qsin, qnsin, qgin, qngin, qhin, qnhin,  &
       qttend, qrtend, qnrtend, qitend, qnitend, qstend, qnstend, qgtend, qngtend, qhtend, qnhtend,  &
       prec_c, prec_r, prec_i, prec_s, prec_g, prec_h)

    ! Note that in F90 local variables overwrite global ones, nevertheless
    ! some SB variables are renamed here

    ! UCLA-LES modules

    use defs, only : cpr,p00,alvl,alvi
    use grid, only : dt,dzi_t,zm, &
!LINDA
         nstep, mp_qt, mp_qr, mp_qi, mp_qs, mp_qg, mp_qh, &
         mp_nqr, mp_nqi, mp_nqs, mp_nqg, mp_nqh, mp_tlt
    use thrm, only : rsif

    ! KAMM2 modules

    USE wolken_driver,     ONLY: loc_ix, loc_iy, loc_iz,           &
         &                       dt_kamm2 => dt,                   &
         &                       dz_kamm2 => dz3d,                 &
         &                       w_kamm2 => w,                     &
         &                       T_kamm2 => T,                     &
         &                       p_0, T_0, rho_k=>rho_0, q,        &
         &                       q_cloud, n_cloud,                 &
         &                       q_ice, q_rain, q_snow, q_graupel, &
         &                       n_ice, n_rain, n_snow, n_graupel, &
         &                       prec_ice, prec_rain,              &
         &                       prec_snow, prec_graupel, prec,    &
         &                       n_hail, q_hail, prec_hail,        &
         &                       alloc_driver, dealloc_driver, cp, &
         &                       S_w, S_i,dSwdz, dSidz, dT0dz,     &
         &                       zml_k,                            &
         &                       w_cb,                             &
         &                       speichere_umwandlungsraten,       &
         &                       speichere_precipstat,             &
         &                       dmin_wg_g, pvec_wg_g, Tvec_wg_g,  &
         &                       qwvec_wg_g, qivec_wg_g, &
         &                       anzp_wg, anzT_wg, anzi_wg, anzw_wg, &
         &                       ltabdminwgg

!    USE parallele_umgebung, ONLY:                                  &
!         &                       isIO,double_global_maxval

    USE wolken,            ONLY: alloc_wolken, dealloc_wolken,     &
         &                       clouds,                           &
         &                       rrho_04, rrho_c

    USE wolken_eis, ONLY : init_dmin_wetgrowth, &
         init_dmin_wg_gr_ltab_equi

    USE wolken_konstanten, ONLY: init_seifert,                     &
         L_ew,L_ed,R_d,R_l,                &
         myparticle=>particle,             &
         mycloud=>cloud,myrain=>rain,myice=>ice,&
         mysnow=>snow,mygraupel=>graupel,myhail=>hail, &
         e_ws,e_es,                        &
         e_ws_vec,e_es_vec,                &
         T_nuc, T_f, satad_nach_mikrophysik, & 
         graupel_shedding, hail_shedding,    &
         qnc_const

    IMPLICIT NONE

    ! Declare variables in argument list

    integer, intent (in) :: level,ie,je,ke

    real, dimension(ke)      , intent (in)             :: &
         dn0,   & ! Density
         pi0,   & ! base state pressure 1
         pi1      ! base state pressure 2

    real, dimension(ke,je,ie), intent (inout)          :: &
         exner, & ! exner function
         thl,   & ! liquid water potential temperature
         tk,    & ! temperature
         w        ! vertical velocity

    real, dimension(ke,je,ie), intent (in) ::  & 
         qvin,qcin,qrin,qnrin,qiin,qniin,qsin,qnsin,qgin,qngin,qhin,qnhin

    real, dimension(ke,je,ie), intent (inout) ::  & 
         tlt,   & ! tendency of liquid water potential temperature
         qttend,qrtend,qnrtend,qitend,qnitend,qstend,qnstend,qgtend,qngtend,qhtend,qnhtend

    real, dimension(ke,je,ie), intent (inout) :: &
         prec_c, prec_r, prec_i, prec_s, prec_g, prec_h

    ! ... Local Variables

    real, dimension(ke,je,ie)                 ::  &
         qv,qc,qr,qnr,qi,qni,qs,qns,qg,qng,qh,qnh

    TYPE(particle) :: pice,psnow,pgraupel,phail

    INTEGER        :: its,ite,jts,jte,kts,kte
    INTEGER        :: ims,ime,jms,jme,kms,kme
    INTEGER        :: i,j,k,ii,jj,kk,n,nt_kamm2,ntsedi,igridpoints,izstat,kkk
    INTEGER, SAVE  :: firstcall, firstcall_init
    REAL(KIND=8)   :: wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,&
         &            precmax,rhomax,rhomin,ncmax,nimax,dzmin,tmax,tmin
    REAL(KIND=8) :: q_vap_new,q_vap_old
    REAL(KIND=8) :: q_liq_new,q_liq_old
    REAL(KIND=8) :: q_ice_new,q_ice_old
    REAL(KIND=8) :: q_v, rho_v, hlp, hlp2, hlp3, dt_sedi,e_v,T_a,qerr, mrho, zxi, zxc
    REAL(KIND=8) :: convliq,convice
    REAL(KIND=8), DIMENSION(10) :: timing

    INTEGER, DIMENSION(:), AllOCATABLE     :: ilm,jlm,klm
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: ssi 

    ! using eps=0 increases runtime, but may remove artifacts
    REAL,    PARAMETER :: eps = 0.0
    !REAL,    PARAMETER :: eps = 1e-10
    LOGICAL, PARAMETER :: debug = .false.
    LOGICAL            :: debug_maxval = .true.
    LOGICAL            :: debug_timing = .true.
    LOGICAL, PARAMETER :: CGP_SEARCH = .true.

    INTEGER :: izerror,unitnr,nt

    LOGICAL, PARAMETER :: lgshed_ascm = .false.
    LOGICAL, PARAMETER :: lhshed_ascm = .false.

    character (len=4)             :: bcstr

    IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: start"

    if (mod(istep,100).eq.0) then
       debug_maxval = .true.
    else
       debug_maxval = .false.
    end if
    if (mod(istep,1234).eq.0.and.nstep.eq.3) then
       debug_timing = .true.
    else
       debug_timing = .false.
    end if

    qnc_const = ccn

    qv = qvin
    qc = qcin
    qr = qrin
    qi = qiin
    qs = qsin
    qg = qgin
    qh = qhin
    qnr = qnrin
    qni = qniin
    qns = qnsin
    qng = qngin
    qnh = qnhin

    ! memory domain (with halo)
    ims = 1
    ime = ie
    jms = 1
    jme = je
    kms = 1
    kme = ke

    ! computational domain (without halo)
    its = 3
    ite = ie-2
    jts = 3
    jte = je-2
    kts = 2
    kte = ke-1
    loc_ix = (jte-jts+1)*(kte-kts+1)*(ite-its+1)

    IF (isIO().and.debug_maxval) &
         & WRITE (*,'(1X,A,I4,A,E12.4)') "mcrph_sb: cloud_type = ",cloud_type,", CCN = ",qnc_const

    IF (firstcall_init /= 1) THEN

       !CALL init_seifert( cloud_type )

       ! initialisiere die lookup-table fuer die graupel-wet-growth parametrisierung:
       IF (isIO()) WRITE (*,*) "mcrph_sb: init_dmin_wetgrowth"

       unitnr = 11
       CALL init_dmin_wetgrowth('datafiles/dmin_wetgrowth_lookup.dat', unitnr)
       CALL init_dmin_wg_gr_ltab_equi(&
            'datafiles/dmin_wetgrowth_lookup.dat', &
            unitnr, 61, ltabdminwgg)
       IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: finished init_dmin_wetgrowth"

       graupel_shedding = lgshed_ascm
       hail_shedding    = lhshed_ascm
       speichere_precipstat = .FALSE.
       speichere_umwandlungsraten = .FALSE.

       firstcall_init = 1
    END IF

    IF (debug_maxval) THEN
       wmax  = double_global_maxval(w(kts:kte,jts:jte,its:ite))
       qvmax = double_global_maxval(qv(kts:kte,jts:jte,its:ite))
       qcmax = double_global_maxval(qc(kts:kte,jts:jte,its:ite))
       qrmax = double_global_maxval(qr(kts:kte,jts:jte,its:ite))
       qimax = double_global_maxval(qi(kts:kte,jts:jte,its:ite))
       qsmax = double_global_maxval(qs(kts:kte,jts:jte,its:ite))
       qgmax = double_global_maxval(qg(kts:kte,jts:jte,its:ite))
       qhmax = double_global_maxval(qh(kts:kte,jts:jte,its:ite))
       nimax = double_global_maxval(qni(kts:kte,jts:jte,its:ite))
       tmax  = double_global_maxval(tk(kts:kte,jts:jte,its:ite))
       IF (isIO()) THEN
          WRITE (*,*) "mcrph_sb: output before microphysics"
          WRITE(*,'(10x,a)') 'Maximum Values:'
          WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','ni','tmax'
          WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,nimax,tmax
       ENDIF
    ENDIF
    IF (debug_timing) THEN       
!    IF (isIO().and.debug_timing) THEN       
       nt=1
       timing = 0.0
       call cpu_time(timing(nt))
    END IF

    allocate( ilm(0:loc_ix), jlm(0:loc_ix), klm(0:loc_ix))
    IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: ilm allocated, i = ", loc_ix
    allocate( ssi(kts:kte,jts:jte,its:ite))
    IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: ssi allocated"

    !..calculate supersaturations after dynamics/advection
    IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: calculate supersaturations"
    DO jj = jts, jte
       DO ii = its, ite
          DO kk = kts, kte
             ! supersaturation w.r.t. ice
             hlp           = p00 * ((pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp)**cpr
             ssi(kk,jj,ii) = qv(kk,jj,ii)/rsif(hlp,tk(kk,jj,ii)) - 1.0

!             ssi(kk,jj,ii) =  R_d & 
!                  & * dn0(kk) * qv(kk,jj,ii) &
!                  & * tk(kk,jj,ii) &
!                  & / e_es(dble(tk(kk,jj,ii))) - 1.0
          ENDDO
       ENDDO
    ENDDO

    !..search for cloudy grid points and store locations
    i = -1
    IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: search for cloudy grid points"
    IF (CGP_SEARCH) THEN
       DO jj = jts, jte
          DO ii = its, ite
             DO kk = kts, kte
                IF (  ssi(kk,jj,ii)  > eps .or. &
                     & qc(kk,jj,ii)  > eps .or. &
                     & qr(kk,jj,ii)  > eps .or. &
                     & qi(kk,jj,ii)  > eps .or. &
                     & qs(kk,jj,ii)  > eps .or. &
                     & qg(kk,jj,ii)  > eps .or. &
                     & qh(kk,jj,ii)  > eps ) THEN
                   i = i+1
                   ilm(i) = ii     ! they start with i=0
                   jlm(i) = jj
                   klm(i) = kk
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO jj = jts, jte
          DO ii = its, ite
             DO kk = kts, kte
                i = i+1
                ilm(i) = ii
                jlm(i) = jj
                klm(i) = kk
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    igridpoints = i

    !..return now, if no clouds are found
    IF (igridpoints == -1) THEN

       IF (debug) WRITE (*,*) "mcrph_sb: no clouds found!"
       DEALLOCATE(ssi)
       DEALLOCATE(ilm,jlm,klm)

       IF (debug) WRITE (*,*) "mcrph_sb: deallocated"

    ELSE ! cloudy points have been found

       IF (debug.and.isIO()) WRITE (0,'(a,i5,a,i10,a,f6.2)') &
            & " mcrph_sb: myid = ",myid,", number of cloudy grid points = ",i+1, &
            & ", percent of domain = ", 100.0*(i+1.0)/loc_ix

       loc_ix = igridpoints
       loc_iy = 1
       loc_iz = 1
       IF (debug.and.isIO()) THEN
          WRITE (*,*) "mcrph_sb: grid"
          WRITE (*,*) "      its = ",its
          WRITE (*,*) "      ite = ",ite
          WRITE (*,*) "      jts = ",jts
          WRITE (*,*) "      jte = ",jte
          WRITE (*,*) "      kts = ",kts
          WRITE (*,*) "      kte = ",kte
          WRITE (*,*) "      loc_ix = ",loc_ix
          WRITE (*,*) "      loc_iy = ",loc_iy
          WRITE (*,*) "      loc_iz = ",loc_iz
       ENDIF
       IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: calling allocation"
       ! ... Allocate memory to temporary KAMM2 variables
       CALL alloc_driver()
       CALL alloc_wolken()
       IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: loc_ix = ", loc_ix
       ! ... transpose to one-dimensional array and variables used in cloud module
       j = 1
       k = 1
       IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: copy to one-dimensional array"
       do i=0,loc_ix
          ! ... grid points
          ii = ilm(i)
          jj = jlm(i)
          kk = klm(i)
           
          hlp = R_l * (1.0 +  (R_d/R_l-1.0) * qv(kk,jj,ii)) 

          ! ... dynamics
          T_0(i,j,k)      = tk(kk,jj,ii)
          p_0(i,j,k)      = p00 * ((pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp)**cpr
          rho_k(i,j,k)    = p_0(i,j,k)/(hlp*tk(kk,jj,ii))
!         rho_k(i,j,k)    = dn0(kk)

          ! .. the ice supersaturation
          S_i(i,j,k)   = ssi(kk,jj,ii)

          ! ... concentrations --> number densities
          n_rain(i,j,k)    = rho_k(i,j,k) * dble( qnr(kk,jj,ii) )
          n_ice(i,j,k)     = rho_k(i,j,k) * dble( qni(kk,jj,ii) )
          n_snow(i,j,k)    = rho_k(i,j,k) * dble( qns(kk,jj,ii) )
          n_graupel(i,j,k) = rho_k(i,j,k) * dble( qng(kk,jj,ii) )
          n_hail(i,j,k)    = rho_k(i,j,k) * dble( qnh(kk,jj,ii) )

          ! ... mixing ratios -> mass densities
          q(i,j,k)         = rho_k(i,j,k) * dble( qv(kk,jj,ii) )
          q_cloud(i,j,k)   = rho_k(i,j,k) * dble( qc(kk,jj,ii) )
          q_rain(i,j,k)    = rho_k(i,j,k) * dble( qr(kk,jj,ii) )
          q_ice(i,j,k)     = rho_k(i,j,k) * dble( qi(kk,jj,ii) )
          q_snow(i,j,k)    = rho_k(i,j,k) * dble( qs(kk,jj,ii) )
          q_graupel(i,j,k) = rho_k(i,j,k) * dble( qg(kk,jj,ii) )
          q_hail(i,j,k)    = rho_k(i,j,k) * dble( qh(kk,jj,ii) )

          w_kamm2(i,j,k)   = w(kk,jj,ii) !+ 0.5 * tke(kk,jj,ii)

       enddo
       IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: finished copy to one-dimensional array"

       ! ... timestep
       dt_kamm2 = dt

       IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: calling clouds"


       IF (debug.and.debug_maxval) THEN
          wmax  = global_maxval(w_kamm2)
          qvmax = global_maxval(q)
          qcmax = global_maxval(q_cloud)
          qrmax = global_maxval(q_rain)
          qimax = global_maxval(q_ice)
          qsmax = global_maxval(q_snow)
          qgmax = global_maxval(q_graupel)
          qhmax = global_maxval(q_hail)
          nimax = global_maxval(n_ice)
          ncmax = global_maxval(n_cloud)
          tmax  = global_maxval(T_0)
          IF (isIO()) THEN
             WRITE (*,*) "mcrph_sb: output before clouds (cloudy grid points only)"
             WRITE(*,'(10x,a)') 'Maximum Values:'
             WRITE(*,'(A10,11A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','ni','nc','tmax'
             WRITE(*,'(A10,11E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,nimax,ncmax,tmax
          ENDIF
       ENDIF

       IF (debug_timing) THEN
!       IF (isIO().and.debug_timing) THEN
          nt=nt+1
          call cpu_time(timing(nt))
       END IF

       ! .. this subroutine calculates all the microphysical sources and sinks
       CALL clouds ()

       IF (debug.and.debug_maxval) THEN
          wmax  = global_maxval(w_kamm2)
          qvmax = global_maxval(q)
          qcmax = global_maxval(q_cloud)
          qrmax = global_maxval(q_rain)
          qimax = global_maxval(q_ice)
          qsmax = global_maxval(q_snow)
          qgmax = global_maxval(q_graupel)
          qhmax = global_maxval(q_hail)
          nimax = global_maxval(n_ice)
          tmax  = global_maxval(T_0)
          IF (isIO()) THEN
             WRITE (*,*) "mcrph_sb: output after clouds (cloudy grid points only)"
             WRITE(*,'(10x,a)') 'Maximum Values:'
             WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','ni','tmax'
             WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,nimax,tmax
          ENDIF
       ENDIF

       IF (debug_timing) THEN
!       IF (isIO().and.debug_timing) THEN
          nt=nt+1
          call cpu_time(timing(nt))
       END IF


       IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: back from clouds"

       ! ... output coeffs and density correction for fall velocity of precip particles
       IF (isIO().AND.firstcall.NE.1.AND.debug) THEN
          WRITE (*,*) "mcrph_sb: particles and PSD shapes"
          WRITE (*,*) "          type = ", mycloud%name
          WRITE (*,*) "          nu   = ", mycloud%nu
          WRITE (*,*) "          mu   = ", mycloud%mu
          WRITE (*,*) "          type = ", myrain%name
          WRITE (*,*) "          nu   = ", myrain%nu
          WRITE (*,*) "          mu   = ", myrain%mu
          WRITE (*,*) "          type = ", myice%name
          WRITE (*,*) "          nu   = ", myice%nu
          WRITE (*,*) "          mu   = ", myice%mu
          WRITE (*,*) "          type = ", mysnow%name
          WRITE (*,*) "          nu   = ", mysnow%nu
          WRITE (*,*) "          mu   = ", mysnow%mu
          WRITE (*,*) "          type = ", mygraupel%name
          WRITE (*,*) "          nu   = ", mygraupel%nu
          WRITE (*,*) "          mu   = ", mygraupel%mu
          WRITE (*,*) "          type = ", myhail%name
          WRITE (*,*) "          nu   = ", myhail%nu
          WRITE (*,*) "          mu   = ", myhail%mu
          WRITE (*,*)
          firstcall = 1
       ENDIF

       IF (debug) THEN
          IF (MINVAL(q) < 0.0) THEN
             WRITE (*,*) ' mcrph_sb: q < 0 ENCOUNTERED, FILLED WITH 0.0'
             WHERE (q < 0.0) q = 0.0
          ENDIF
          IF (MINVAL(q_cloud) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_cloud < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_rain) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_rain < 0, STOPPED AFTER CLOUDS 1', MINVAL(q_rain)
             stop
          ENDIF
          IF (MINVAL(q_ice) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_ice < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_snow) < 0.0) THEN
             write (*,*) ' mcrph_sb: q_snow < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_graupel) < 0.0) THEN
             WRITE (*,*) ' mcrph_sb: q_graupel < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(q_hail) < 0.0) THEN
             WRITE (*,*) ' mcrph_sb: q_hail < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_cloud) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_cloud < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_rain) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_rain < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_ice) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_ice < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_snow) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_snow < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_graupel) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_graupel < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
          IF (MINVAL(n_hail) < 0.0) THEN
             write (*,*) ' mcrph_sb: n_hail < 0, STOPPED AFTER CLOUDS 1'
             stop
          ENDIF
       ENDIF

       ! ... Transformation of variables back to driving model and latent heat equation
       IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: trafo back"

       j = 1
       k = 1
       DO i=0,loc_ix
          ii = ilm(i)
          jj = jlm(i)
          kk = klm(i)

          hlp = 1.0 / rho_k(i,j,k)
          q_v = hlp * q(i,j,k)

          ! ... latent heat for temperature equation
          !     in UCLA-LES the qc is part of theta_l, but qr is not

          ! ... ... new variables
          q_vap_new = q_v + hlp*q_cloud(i,j,k)
          q_liq_new = hlp *  q_rain(i,j,k)
          q_ice_new = hlp * ( q_ice(i,j,k)+q_snow(i,j,k)+q_graupel(i,j,k)+q_hail(i,j,k) )
          ! ... ... old variables
          q_vap_old = qv(kk,jj,ii)+qc(kk,jj,ii)
          q_liq_old = qr(kk,jj,ii)
          q_ice_old = qi(kk,jj,ii) + qs(kk,jj,ii) + qg(kk,jj,ii) + qh(kk,jj,ii)

          ! ... temperature tendency (for liquid water potential temperature)
          convice = alvi/(cp*(pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp)
          convliq = (alvl-alvi)/(cp*(pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp)
          tlt(kk,jj,ii) = tlt(kk,jj,ii) &
               &      - convice * (q_vap_new - q_vap_old) / dt  &
               &      + convliq * (q_liq_new - q_liq_old) / dt

          ! LINDA, b, write out tendencies
          if (lmptend) then
            if (nstep==1) mp_tlt(kk,jj,ii)=0.0
            mp_tlt(kk,jj,ii) = mp_tlt(kk,jj,ii)                          & 
                              - convice * (q_vap_new - q_vap_old) / dt  &
                              + convliq * (q_liq_new - q_liq_old) / dt
          end if
          ! LINDA, e
          ! ... mass densities to mixing ratios with actual density:
          qv(kk,jj,ii) = q_v
          qc(kk,jj,ii) = hlp * q_cloud(i,j,k)
          qr(kk,jj,ii) = hlp * q_rain(i,j,k)
          qi(kk,jj,ii) = hlp * q_ice(i,j,k)
          qs(kk,jj,ii) = hlp * q_snow(i,j,k)
          qg(kk,jj,ii) = hlp * q_graupel(i,j,k)
          qh(kk,jj,ii) = hlp * q_hail(i,j,k)

          ! ... number concentrations
          qnr(kk,jj,ii) = hlp * n_rain(i,j,k)
          qni(kk,jj,ii) = hlp * n_ice(i,j,k)
          qns(kk,jj,ii) = hlp * n_snow(i,j,k)
          qng(kk,jj,ii) = hlp * n_graupel(i,j,k)
          qnh(kk,jj,ii) = hlp * n_hail(i,j,k)

       ENDDO

       deallocate(ilm,jlm,klm)

       if (debug) then
          WHERE (qv < 0.0) qv = 0.0
          WHERE (qc < 0.0) qc = 0.0
          WHERE (qr < 0.0) qr = 0.0
          WHERE (qi < 0.0) qi = 0.0
          WHERE (qs < 0.0) qs = 0.0
          WHERE (qg < 0.0) qg = 0.0
          WHERE (qh < 0.0) qh = 0.0
          WHERE (qnr < 0.0) qnr = 0.0
          WHERE (qni < 0.0) qni = 0.0
          WHERE (qns < 0.0) qns = 0.0
          WHERE (qng < 0.0) qng = 0.0
          WHERE (qnh < 0.0) qnh = 0.0
       end if

       ! ... deallocation
       IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: calling deallocation_driver"
       CALL dealloc_driver()
       IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: calling deallocation_wolken"
       CALL dealloc_wolken()
       IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: deallocation"
       DEALLOCATE(ssi)

       IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: end"

    END IF ! ... This ends the loooong if-block 'cloudy points present'

    pice     = copy_particle(myice,0.00)
    psnow    = copy_particle(mysnow,0.25)
    pgraupel = copy_particle(mygraupel,0.00)
    phail    = copy_particle(myhail,0.00)

    ! need exactly those names to use UCLA-LES sedimentation code
    pice%name     = 'ice'
    psnow%name    = 'snow'
    pgraupel%name = 'graupel'
    phail%name    = 'hail'

    IF (debug_timing) THEN
!    IF (isIO().and.debug_timing) THEN
       nt=nt+1
       call cpu_time(timing(nt))
    END IF

    IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: BCs for sedimentation"

    !bcstr = 'cnst'   ! 'mixd' or 'cnst'
    !call sclrset(bcstr,ke,je,ie,qr ,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qnr,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qi ,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qni,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qs ,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qns,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qg ,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qng,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qh ,dzi_t)
    !call sclrset(bcstr,ke,je,ie,qnh,dzi_t)

    IF (debug.and.isIO()) WRITE (0,*) "mcrph_sb: sedimentation"

    prec_c = 0.0
    prec_r = 0.0
    DO j=3,je-2
      DO i=3,ie-2
         !    call sedim_cd(ke,dt,thl(1:ke,j,i),qc(1:ke,j,i),prec_c(1:ke,j,i))
         call sedim_rd(ke,dt,dn0,qr(1:ke,j,i),qnr(1:ke,j,i),prec_r(1:ke,j,i))
      end do
    end do

    if (cloud_type.ge.1000) then
    prec_i = 0.0
    prec_s = 0.0
    prec_g = 0.0
    DO j=3,je-2
      DO i=3,ie-2
         IF (ANY(qi(1:ke,j,i).gt.0.0)) call sedimentation (ke,dn0,pice, qi(1:ke,j,i),qni(1:ke,j,i),rrate = prec_i(1:ke,j,i))
         IF (ANY(qs(1:ke,j,i).gt.0.0)) call sedimentation (ke,dn0,psnow,qs(1:ke,j,i),qns(1:ke,j,i),rrate = prec_s(1:ke,j,i))
      end do
    end do
    ntsedi = 3
    DO j=3,je-2
      DO i=3,ie-2
        IF (ANY(qg(1:ke,j,i).gt.0.0)) THEN
          DO ii=1,ntsedi
            call sedimentation (ke,dn0,pgraupel,qg(1:ke,j,i),qng(1:ke,j,i),rrate = prec_g(1:ke,j,i),dtopt=dt/ntsedi)
          end do
        ENDIF
      end do
    end do
    end if
    if (cloud_type.ge.2000) then
    prec_h = 0.0
    ntsedi = 3
    DO j=3,je-2
      DO i=3,ie-2
        IF (ANY(qh(1:ke,j,i).gt.0.0)) THEN
          DO ii=1,ntsedi
            call sedimentation (ke,dn0,phail,qh(1:ke,j,i),qnh(1:ke,j,i),rrate = prec_h(1:ke,j,i),dtopt=dt/ntsedi)
          end do
        ENDIF
      end do
    end do
    end if

    IF (debug_timing) THEN
!    IF (isIO().and.debug_timing) THEN
       nt=nt+1
       call cpu_time(timing(nt))
    END IF

    IF (debug) THEN
       IF (MINVAL(qc) < 0.0) THEN
          write (*,*) ' mcrph_sb: qc < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qi) < 0.0) THEN
          write (*,*) ' mcrph_sb: qi < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qs) < 0.0) THEN
          write (*,*) ' mcrph_sb: qs < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qg) < 0.0) THEN
          write (*,*) ' mcrph_sb: qg < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qh) < 0.0) THEN
          write (*,*) ' mcrph_sb: qh < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qr) < 0.0) THEN
          write (*,*) ' mcrph_sb: qr < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qni) < 0.0) THEN
          write (*,*) ' mcrph_sb: qni < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qns) < 0.0) THEN
          write (*,*) ' mcrph_sb: qns < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qng) < 0.0) THEN
          write (*,*) ' mcrph_sb: qng < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qnh) < 0.0) THEN
          write (*,*) ' mcrph_sb: qnh < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
       IF (MINVAL(qnr) < 0.0) THEN
          write (*,*) ' mcrph_sb: qnr < 0, STOPPED AFTER SEDIMENTATION'
          stop
       ENDIF
    ENDIF

    WHERE (qv < 0.0) qv = 0.0
    WHERE (qc < 0.0) qc = 0.0
    WHERE (qr < 0.0) qr = 0.0
    WHERE (qi < 0.0) qi = 0.0
    WHERE (qs < 0.0) qs = 0.0
    WHERE (qg < 0.0) qg = 0.0
    WHERE (qh < 0.0) qh = 0.0
    WHERE (qnr < 0.0) qnr = 0.0
    WHERE (qni < 0.0) qni = 0.0
    WHERE (qns < 0.0) qns = 0.0
    WHERE (qng < 0.0) qng = 0.0
    WHERE (qnh < 0.0) qnh = 0.0

    DO i=1,ie
       DO j=1,je
          DO k=2,ke
             qttend(k,j,i) = qttend(k,j,i) + (qv(k,j,i) - qvin(k,j,i))/dt &
                                           + (qc(k,j,i) - qcin(k,j,i))/dt

             qrtend(k,j,i)  = max(qrtend(k,j,i) + (qr(k,j,i) - qrin(k,j,i))/dt,-qrin(k,j,i)/dt)
             qitend(k,j,i)  = max(qitend(k,j,i) + (qi(k,j,i) - qiin(k,j,i))/dt,-qiin(k,j,i)/dt)
             qstend(k,j,i)  = max(qstend(k,j,i) + (qs(k,j,i) - qsin(k,j,i))/dt,-qsin(k,j,i)/dt)
             qgtend(k,j,i)  = max(qgtend(k,j,i) + (qg(k,j,i) - qgin(k,j,i))/dt,-qgin(k,j,i)/dt)
             qhtend(k,j,i)  = max(qhtend(k,j,i) + (qh(k,j,i) - qhin(k,j,i))/dt,-qhin(k,j,i)/dt)

             qnrtend(k,j,i) = qnrtend(k,j,i) + (qnr(k,j,i) - qnrin(k,j,i))/dt
             qnitend(k,j,i) = qnitend(k,j,i) + (qni(k,j,i) - qniin(k,j,i))/dt
             qnstend(k,j,i) = qnstend(k,j,i) + (qns(k,j,i) - qnsin(k,j,i))/dt
             qngtend(k,j,i) = qngtend(k,j,i) + (qng(k,j,i) - qngin(k,j,i))/dt
             qnhtend(k,j,i) = qnhtend(k,j,i) + (qnh(k,j,i) - qnhin(k,j,i))/dt
! LINDA, b, write out tendenties
             if (lmptend) then
                if (nstep==1) then
                    mp_qt(k,j,i) = 0.0
                    mp_qr(k,j,i) = 0.0
                    mp_qi(k,j,i) = 0.0
                    mp_qs(k,j,i) = 0.0
                    mp_qg(k,j,i) = 0.0
                    mp_qh(k,j,i) = 0.0

                    mp_nqr(k,j,i) = 0.0
                    mp_nqi(k,j,i) = 0.0
                    mp_nqs(k,j,i) = 0.0
                    mp_nqg(k,j,i) = 0.0
                    mp_nqh(k,j,i) = 0.0
                endif

                mp_qt  (k,j,i) = mp_qt(k,j,i)                 &
                                + (qv(k,j,i) - qvin(k,j,i))/dt &
                                + (qc(k,j,i) - qcin(k,j,i))/dt
                mp_qr  (k,j,i) = mp_qr(k,j,i) + (qr(k,j,i) - qrin(k,j,i))/dt
                mp_qi  (k,j,i) = mp_qi(k,j,i) + (qi(k,j,i) - qiin(k,j,i))/dt
                mp_qs  (k,j,i) = mp_qs(k,j,i) + (qs(k,j,i) - qsin(k,j,i))/dt
                mp_qg  (k,j,i) = mp_qg(k,j,i) + (qg(k,j,i) - qgin(k,j,i))/dt
                mp_qh  (k,j,i) = mp_qh(k,j,i) + (qh(k,j,i) - qhin(k,j,i))/dt

                mp_nqr (k,j,i) = mp_nqr(k,j,i) + (qnr(k,j,i) - qnrin(k,j,i))/dt
                mp_nqi (k,j,i) = mp_nqi(k,j,i) + (qni(k,j,i) - qniin(k,j,i))/dt
                mp_nqs (k,j,i) = mp_nqs(k,j,i) + (qns(k,j,i) - qnsin(k,j,i))/dt
                mp_nqg (k,j,i) = mp_nqg(k,j,i) + (qng(k,j,i) - qngin(k,j,i))/dt
                mp_nqh (k,j,i) = mp_nqh(k,j,i) + (qnh(k,j,i) - qnhin(k,j,i))/dt
              end if
              ! LINDA, e
          END DO
       END DO
    END DO

    IF (debug_maxval) THEN
       wmax  = double_global_maxval(w)
       qvmax = double_global_maxval(qv)
       qcmax = double_global_maxval(qc)
       qrmax = double_global_maxval(qr)
       qimax = double_global_maxval(qi)
       qsmax = double_global_maxval(qs)
       qgmax = double_global_maxval(qg)
       qhmax = double_global_maxval(qh)
       nimax = double_global_maxval(qni)
       tmax  = double_global_maxval(tk)
       IF (isIO()) THEN
          WRITE (*,*) "mcrph_sb: output after microphysics and sedimentation"
          WRITE(*,'(10x,a)') 'Maximum Values:'
          WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','ni','tmax'
          WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,nimax,tmax
       ENDIF
    ENDIF


    IF (debug_timing) THEN
!    IF (isIO().and.debug_timing) THEN
       nt=nt+1
       call cpu_time(timing(nt))
       WRITE(*,'(1x,a,i3,a,f8.3,a)') 'mcrph_sb: CPU ',myid,' time  t2-t1(sec) = ',timing(2)-timing(1),' (trafo1)'
       WRITE(*,'(1x,a,i3,a,f8.3,a)') '          CPU ',myid,' time  t3-t2(sec) = ',timing(3)-timing(2),' (clouds)'
       WRITE(*,'(1x,a,i3,a,f8.3,a)') '          CPU ',myid,' time  t4-t3(sec) = ',timing(4)-timing(3),' (trafo2)'
       WRITE(*,'(1x,a,i3,a,f8.3,a)') '          CPU ',myid,' time  t5-t4(sec) = ',timing(5)-timing(4),' (sedi)'
       WRITE(*,'(1x,a,i3,a,f8.3,a)') '          CPU ',myid,' time  t6-t5(sec) = ',timing(6)-timing(5),' (final)'
       WRITE(*,'(1x,a,i3,a,f8.3)')   '          CPU ',myid,' time  total(sec) = ',timing(nt)-timing(1)
    END IF

    IF (debug.and.isIO()) WRITE (*,*) "mcrph_sb: end"

    RETURN

  END SUBROUTINE mcrph_sb

  SUBROUTINE init_sb

    USE wolken_konstanten, ONLY: init_seifert,                     &
         myparticle=>particle,             &
         mycloud=>cloud,myrain=>rain,myice=>ice,&
         mysnow=>snow,mygraupel=>graupel,myhail=>hail

    call init_seifert(cloud_type)

    cldw    = copy_particle(mycloud,0.00)
    rain    = copy_particle(myrain,0.00)
    ice     = copy_particle(myice,0.00)
    snow    = copy_particle(mysnow,0.25)
    graupel = copy_particle(mygraupel,0.00)
    hail    = copy_particle(myhail,0.00)

  END SUBROUTINE init_sb

  FUNCTION copy_particle(p,s_vel) result(q)

    USE wolken_konstanten, ONLY: myparticle=>particle

    IMPLICIT NONE
    TYPE(myparticle) :: p  ! axel's particle
    TYPE(particle)   :: q  ! thijs' particle
    REAL             :: s_vel

    q%name  = p%name
    q%moments = 2
    q%nr    = 0.0
    q%nu    = p%nu
    q%mu    = p%mu
    q%x_max = p%x_max
    q%x_min = p%x_min
    q%a_geo = p%a_geo
    q%b_geo = p%b_geo
    q%a_vel = p%a_vel
    q%b_vel = p%b_vel
    q%s_vel = s_vel
    q%a_ven = p%a_ven
    q%b_ven = p%b_ven
    q%cap   = p%cap

  END FUNCTION COPY_PARTICLE

end module mcrp
