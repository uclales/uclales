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
! TODOLIST ICEMICRO
! *Include irina - alterations !DONE
! *Remove irreversible water frm q_t !DONE
! *Work do-loops outward   !DONE
! *Include Oliviers modifications  !DONE
! *Do serialization of warmmicrophysics !DONE
! *Include partitioning function between ice and water !DONE
! *Include Statistics
! *Resetvars
! *Include riming
! *Include selfcollection
! *Include collection



module mcrp

  use defs, only : tmelt,alvl, alvi,rowt, pi, Rm, cp,t_hn
  use grid, only : dt, dxi, dyi ,dzt, nxp, nyp, nzp, a_pexnr, a_rp, a_tp, th00, ccn,    &
       dn0, pi0,pi1, a_rt, a_tt,a_rpp, a_rpt, a_npp, a_npt, vapor, liquid,       &
       a_theta, a_scr1, a_scr2, a_scr7, precip, &
       a_ninuct,a_ricet,a_nicet,a_rsnowt,a_nsnowt,a_rgrt,a_ngrt,&
       a_ninucp,a_ricep,a_nicep,a_rsnowp,a_nsnowp,a_rgrp,a_ngrp,rsup
  use thrm, only : thermo, fll_tkrs, thetal_noprecip
  use stat, only : sflg, updtst
  use util, only : get_avg3, azero
  implicit none

  logical, parameter :: droplet_sedim = .False., khairoutdinov = .False., turbulence = .False.
  logical            :: firsttime = .true.
  integer            :: nprocess

  integer,parameter  :: iwtrdff = 1,iauto = 2,iaccr = 3,isedimrd = 4,isedimcd = 5,iicenucnr = 6,iicenuc = 7,ifreez=8,idep=9,imelt=10,ised_ice=11,iauto_ice=12,iaggr_ice=13
  real, dimension(:),allocatable :: convice,convliq
  integer, dimension(:), allocatable :: seq
  ! 
  ! drop sizes definition is based on vanZanten (2005)
  ! cloud droplets' diameter: 2-50 e-6 m
  ! drizzle drops' diameter: 50-1000 e-6 m
  !
  real, parameter    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
  real, parameter    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]

  real, parameter :: eps0 = 1e-20       ! small number
  real, parameter :: eps1 = 1e-9        ! small number
  real, parameter :: rho_0 = 1.21       ! air density at surface

  real, parameter :: prw = pi * rowt / 6.
  
  real, parameter :: t_nuc     = 273.2d+0 ! maximum temperature for ice nucleation
  real, parameter :: t_freeze  = 273.2d+0 ! maximum temperature for freezing
  real, parameter :: nu_l = 1.460d-5     !..kinem. visc. von luft
  real, parameter :: d_v  = 3.000d-5     !..diffusivitaet von wasserdampf
  real, parameter :: k_t  = 2.500d-2     !..waermeleitfaehigkeit

!   real, parameter :: k_w  = 0.930d+0     !..dielektrizitaetsfaktor wasser
!   real, parameter :: k_i  = 0.176d+0     !..dielektrizitaetsfaktor eis
!   real, parameter :: rho0 = 1.225d+0     !..norm-luftdichte
!   real, parameter :: t_f  = 2.330d+2     !..bei t < t_f kein fl.wasser
!   real, parameter :: r_l  = 2.870d+2     !..gaskonstante trockener luft
!   real, parameter :: c_w  = 4.218d+3     !..spezifische waerme von wasser
!   real, parameter :: sigma_wa  = 7.1d-2  !..oberflaechenspannung wasser-luft
!   real, parameter :: rowt   = 1000.0    !..materialdichte von fluessigwasser
!   real, parameter :: rho_ice = 900.0     !..materialdichte von eis
  ! .. wolkenphysikalische konstanten und parameter ..
! 
  real, parameter :: n_sc = 0.710        !..schmidt-zahl (pk, s.541)
  real, parameter :: n_f  = 0.333        !..exponent von n_sc im vent-koeff. (pk, s.541)
  real, parameter :: m_f  = 0.500        !..exponent von n_re im vent-koeff. (pk, s.541)
! 
  real, parameter :: a_e  = 2.18745584d1 !..konst. saettigungsdamppfdruck - eis
  real, parameter :: a_w  = 1.72693882d1 !..konst. saettigungsdamppfdruck - wasser
  real, parameter :: b_e  = 7.66000000d0 !..konst. saettigungsdamppfdruck - eis
  real, parameter :: b_w  = 3.58600000d1 !..konst. saettigungsdamppfdruck - wasser
  real, parameter :: e_3  = 6.10780000d2 !..saettigungsdamppfdruck bei t = tmelt
! 
!   real, parameter :: c_mult     = 3.5d8    !..koeff. fuer splintering
!   real, parameter :: t_mult_min = 265.0    !..minimale temp. splintering
!   real, parameter :: t_mult_max = 270.0    !..maximale temp. splintering
!   real, parameter :: t_mult_opt = 268.0    !..optimale temp. splintering
! 
!   ! ... spezielle parameter des kamm2-wolkenmoduls
! 
!   real, parameter :: r_c     = 12.0e-6         !..mittlerer radius (bei 1-moment)
!   real, parameter :: rho_vel    = 0.5d0        !..exponent in dichtekorrektur
!   real, parameter :: rho_vel_c  = 1.0d0        !..fuer wolkentropfen
! 
!   ! ... spezielle parameter des kamm2-wolkenmoduls (eisphase)
! 
!   real, parameter :: e_ii  = 0.00              !..min. eff.
!   real, parameter :: e_ic  = 0.80              !..max. eff. fuer ice_cloud_riming
!   real, parameter :: e_sc  = 0.80              !..max. eff. fuer snow_cloud_riming
!   real, parameter :: e_gc  = 1.00              !..max. eff. fuer graupel_cloud_riming
!   real, parameter :: e_hc  = 1.00              !..max. eff. fuer hail_cloud_riming
!   real, parameter :: e_min = 0.01              !..min. eff. fuer gc,ic,sc
!   real, parameter :: alpha_spacefilling = 0.1  !..raumerfuellungskoeff (max. 0.68)
!   real, parameter :: ice_s_vel  = 0.25         !..dispersion der fallgeschw.
!   real, parameter :: snow_s_vel = 0.25         !..dispersion der fallgeschw.
!   real, parameter :: r_shedding = 500.0e-6     !..mittlerer radius shedding
!   real, parameter :: t_shed = 263.2
! 
!   ! rain_freeze: der teil des regenspektrums kleiner als d_rainfrz_ig
!   ! wird nach gefrieren dem eis zugeschlagen, der teil von dort bis zu d_rainfrz_gh dem graupel
!   ! und der rest dem hagel.
!   real, parameter :: d_rainfrz_ig = 0.50d-3 !  rain --> ice oder graupel
!   real, parameter :: d_rainfrz_gh = 1.25d-3 ! rain --> graupel oder hail
! 
!   real, parameter :: q_krit_ii = 1.000d-6 ! q-schwellenwert fuer ice_selfcollection
!   real, parameter :: d_krit_ii = 50.00d-6 ! d-schwellenwert fuer ice_selfcollection
!   real, parameter :: d_conv_ii = 75.00d-6 ! d-schwellenwert fuer ice_selfcollection
!   real, parameter :: q_krit_ic = 1.000d-5 ! q-schwellenwert fuer ice_cloud_riming
!   real, parameter :: d_krit_ic = 150.0d-6 ! d-schwellenwert fuer ice_cloud_riming
!   real, parameter :: q_krit_ir = 1.000d-5 ! q-schwellenwert fuer ice_rain_riming
!   real, parameter :: d_krit_ir = 100.0d-6 ! d-schwellenwert fuer ice_rain_riming
!   real, parameter :: q_krit_sc = 1.000d-5 ! q-schwellenwert fuer snow_cloud_riming
!   real, parameter :: d_krit_sc = 150.0d-6 ! d-schwellenwert fuer snow_cloud_riming
!   real, parameter :: q_krit_sr = 1.000d-5 ! q-schwellenwert fuer snow_rain_riming
!   real, parameter :: d_krit_sr = 100.0d-6 ! d-schwellenwert fuer snow_rain_riming
!   real, parameter :: q_krit_gc = 1.000d-6 ! q-schwellenwert fuer graupel_cloud_riming
!   real, parameter :: d_krit_gc = 100.0d-6 ! d-schwellenwert fuer graupel_cloud_riming
!   real, parameter :: q_krit_hc = 1.000d-6 ! q-schwellenwert fuer hail_cloud_riming
!   real, parameter :: d_krit_hc = 100.0d-6 ! d-schwellenwert fuer hail_cloud_riming
!   real, parameter :: q_krit_fr = 1.000d-6 ! q-schwellenwert fuer rain_freeze
!   real, parameter :: q_krit_c  = 1.000d-6 ! q-schwellenwert sonst
!   real, parameter :: q_krit    = 1.000d-9 ! q-schwellenwert sonst
!   real, parameter :: d_conv_sg = 200.0d-6 ! d-schwellenwert
!   real, parameter :: d_conv_ig = 200.0d-6 ! d-schwellenwert
!   real, parameter :: x_conv    = 0.100d-9 ! minimale graupel-/hagelmasse riming
!   real, parameter :: d_shed_g  = 3.000d-3 ! d-schwellenwert fuer graupel_shedding
!   real, parameter :: d_shed_h  = 5.000d-3 ! d-schwellenwert fuer hagel_shedding
!   real, parameter :: d_krit_c  = 10.00d-6 ! d-schwellenwert fuer cloud_collection
!   real, parameter :: d_coll_c  = 40.00d-6 ! oberer wert fuer cloud_coll_eff
!   real, parameter :: t_nuc     = 273.2d+0 ! temperatur ab der eisnukleation einsetzt
!   real, parameter :: t_freeze  = 273.2d+0 ! temperatur ab der gefrieren einsetzt
! 
!   real, parameter :: q_krit_aus = 1.00d-5 ! q-schwellenwert fuer ausgabe von d und z


  type particle
    integer :: moments
    real :: nr
    real :: nu1    !..width =  nu1*mass+nu2
    real :: nu2
    real :: nu
    real :: mu
    real :: x_max !..maximum particle mass
    real :: x_min !..minimum particle mass
    real :: a_geo !..koeff. geometrie
    real :: b_geo !..koeff. geometrie = 1/3
    real :: a_vel !..koeff. fallgesetz
    real :: b_vel !..koeff. fallgesetz
    real :: a_ven !..koeff. ventilationsparam.
    real :: b_ven !..koeff. ventilationsparam.
    real :: cap   !..koeff. kapazitaet
  end type particle
  type(particle) :: cldw,rain,ice,snow,graupel


contains

  !
  ! ---------------------------------------------------------------------
  ! MICRO: sets up call to microphysics
  !

  subroutine micro(level)

    integer, intent (in) :: level

    select case (level) 
    case(2)
!        if (droplet_sedim) call sedim_cd(level,nzp,nxp,nyp,dn0,dt,a_theta,a_scr1,      &
!             liquid,precip,a_rt,a_tt)     
    case(3)
       call mcrph(level,nzp,nxp,nyp,dn0,a_tp,a_theta,a_scr1,vapor,a_scr2,liquid,a_rpp, &
              a_npp,precip,a_rt,a_tt,a_rpt,a_npt,a_scr7)
    case(4)
       call mcrph(level,nzp,nxp,nyp,dn0,a_tp,a_theta,a_scr1,vapor,a_scr2,liquid,a_rpp, &
              a_npp,precip,a_rt,a_tt,a_rpt,a_npt,a_scr7,&
              a_ninuct,a_ricet,a_nicet,a_rsnowt,a_nsnowt,a_rgrt,a_ngrt,&
              a_ninucp,a_ricep,a_nicep,a_rsnowp,a_nsnowp,a_rgrp,a_ngrp,rsup)
    end select

  end subroutine micro
  !
  ! ---------------------------------------------------------------------
  ! MCRPH: calls microphysical parameterization 
  !

      subroutine mcrph(level,n1,n2,n3,dn0,thl,th,tk,rv,rs,rc,rp,np,         &
       rrate, rtt,tlt,rpt,npt,dissip,    &
       ninuct,ricet,nicet,rsnowt,nsnowt,rgrt,ngrt,ninucp,ricep,nicep,rsnowp,nsnowp,rgrp,ngrp,rsup)

    integer, intent (in) :: level,n1,n2,n3
    real, dimension(n1)      , intent (in)             :: dn0  !Density
    real, dimension(n1,n2,n3), intent (inout)          :: thl,& !Theta_l
                                                          th, & !theta
                                                          tk, & !temperature
                                                          rv, & !water vapor
                                                          rs, & !saturation mixing ration
                                                          rc,& !Condensate/cloud water
                                                          np,&!rain drop number
                                                          rp, &!rain water
                                                          rtt,& !Total water tendency
                                                          tlt,&! theta_l tendency
                                                          rpt,&!rain water tendency
                                                          npt,&!rain droplet number tendency
                                                          dissip
    real, intent (out)                                 :: rrate(n1,n2,n3)
    real, dimension(n1,n2,n3), intent (inout),optional :: ninuct,ricet,nicet,rsnowt,nsnowt,rgrt,ngrt,&
                                                          ninucp,ricep,nicep,rsnowp,nsnowp,rgrp,ngrp,&
                                                          rsup

    real, dimension(n1) :: q1,q2,q3,q4,q5 !dummy-arrays for conversion between mixing rate 
    integer :: i, j,n
    rrate = 0.    
    if(firsttime) call initmcrp(level)
    
    tlt    = tlt    - thl/dt
    rtt    = rtt    - (rv+rc)/dt
    rpt    = rpt    - rp/dt
    npt    = npt    - np/dt
    if (level==4) then
      ninuct = ninuct - ninucp/dt
      ricet  = ricet  - ricep/dt
      nicet  = nicet  - nicep/dt
      rsnowt = rsnowt - rsnowp/dt
      nsnowt = nsnowt - nsnowp/dt
      rgrt   = rgrt   - rgrp/dt
      ngrt   = ngrt   - ngrp/dt
    end if

    seq =  gen_sequence(nprocess)


    do j=3,n3-2
      do i=3,n2-2
        convliq = alvl*th(:,i,j)/(cp*tk(:,i,j))
        convice = alvi*th(:,i,j)/(cp*tk(:,i,j))
        do n=1,nprocess
          select case(seq(n))
          case(iwtrdff)
            call resetvar(cldw,rc(:,i,j))
            call wtr_dff_SB(n1,dn0,rp(:,i,j),np(:,i,j),rc(:,i,j),rs(:,i,j),rv(:,i,j),thl(:,i,j),tk(:,i,j))
          case(iauto)
            call resetvar(cldw,rc(:,i,j))
            call auto_SB(n1,dn0,rc(:,i,j),rp(:,i,j),np(:,i,j),thl(:,i,j),dissip(:,i,j))
          case(iaccr)
            call resetvar(cldw,rc(:,i,j))
            call accr_SB(n1,dn0,rc(:,i,j),rp(:,i,j),np(:,i,j),thl(:,i,j),dissip(:,i,j))
          case(isedimrd)
            call sedim_rd(n1,dt,dn0,rp(:,i,j),np(:,i,j))
          case(isedimcd)
            call sedim_cd(n1,dt,thl(:,i,j),th(:,i,j),tk(:,i,j),rc(:,i,j))
          case(iicenucnr)
            q1 = dn0*rsup(:,i,j)
            call n_icenuc(n1,ninucp(:,i,j),tk(:,i,j),q1)
          case(iicenuc)
            q1 = dn0*ricep(:,i,j)
            q2 = dn0*rsup(:,i,j)
            call ice_nucleation(n1,ninucp,q1,nicep(:,i,j),q2,thl(:,i,j),tk(:,i,j))
            ricep(:,i,j) = q1/dn0
          case(ifreez)
            call resetvar(cldw,rc(:,i,j))
            q1 = dn0*rc(:,i,j)
            q2 = dn0*rp(:,i,j)
            q3 = dn0*ricep(:,i,j)
            q4 = dn0*rgrp(:,i,j)
            call cloud_freeze(n1,q1,q3,nicep(:,i,j),thl(:,i,j),tk(:,i,j))
            call rain_freeze(n1,q2,q3,np(:,i,j),q4,nicep(:,i,j),thl(:,i,j),tk(:,i,j))
            rc(:,i,j)    = q1/dn0
            rp(:,i,j)    = q2/dn0
            ricep(:,i,j) = q3/dn0
            rgrp(:,i,j)  = q4/dn0
          case(idep)
          case(imelt)
            call resetvar(cldw,rc(:,i,j))
            q1 = dn0*rc(:,i,j)
            q2 = dn0*rp(:,i,j)
            q3 = dn0*ricep(:,i,j)
            q4 = dn0*rgrp(:,i,j)
            q5 = dn0*rsnowp(:,i,j)
            call melting(n1,ice,q3,thl(:,i,j),nicep(:,i,j),q1,q2,np(:,i,j),tk(:,i,j))
            call melting(n1,graupel,q4,thl(:,i,j),nicep(:,i,j),q1,q2,np(:,i,j),tk(:,i,j))
            call melting(n1,snow,q5,thl(:,i,j),nicep(:,i,j),q1,q2,np(:,i,j),tk(:,i,j))
            rc(:,i,j)    = q1/dn0
            rp(:,i,j)    = q2/dn0
            ricep(:,i,j) = q3/dn0
            rgrp(:,i,j)  = q4/dn0
            rsnowp(:,i,j)= q5/dn0
          case(ised_ice)
            call sedimentation (n1,ice, ricep(:,i,j),nicep(:,i,j))
            call sedimentation (n1,graupel, rgrp(:,i,j))
            call sedimentation (n1,snow, rsnowp(:,i,j))
          case(iauto_ice)
          case(iaggr_ice)
          case default
          end select
        end do
      end do
    end do

    tlt    = tlt    + thl/dt
    rtt    = rtt    + (rv+rc)/dt
    rpt    = rpt    + rp/dt
    npt    = npt    + np/dt
    if (level==4) then
      ninuct = ninuct + ninucp/dt
      ricet  = ricet  + ricep/dt
      nicet  = nicet  + nicep/dt
      rsnowt = rsnowt + rsnowp/dt
      nsnowt = nsnowt + nsnowp/dt
      rgrt   = rgrt   + rgrp/dt
      ngrt   = ngrt   + ngrp/dt
    end if
  end subroutine mcrph
  subroutine wtr_dff_SB(n1,dn0,rp,np,rl,rs,rv,tl,tk)
 !
 ! ---------------------------------------------------------------------
 ! WTR_DFF_SB: calculates the evolution of the both number- and
 ! mass mixing ratio large drops due to evaporation in the absence of
 ! cloud water.
 !

    integer, intent (in) :: n1
    real, intent (in)    :: tk(n1),rs(n1), dn0(n1)
    real, intent (inout) :: rp(n1), np(n1),tl(n1),rv(n1),rl(n1)

    real, parameter     :: c_Nevap = 0.7
    integer             :: k
    real                :: Xp, Dp, G, S, cerpt, cenpt
    do k=2,n1
      if (rp(k) > 0 .and. rl(k)<=0.) then
          Xp = rp(k)/ (np(k)+eps0)
          Dp = ( Xp / prw )**(1./3.)
          G = 1. / (1. / (dn0(k)*rs(k)*Dv) + &
                alvl*(alvl/(Rm*tk(k))-1.) / (Kt*tk(k)))
          S = rv(k)/rs(k) - 1.

          if (S < 0) then
            cerpt = 2. * pi * Dp * G * S * np(k)
            cenpt = c_Nevap*cerpt * np(k) / rp(k)
            np(k)=np(k) + cenpt*dt
            rl(k)=rl(k) + cerpt*dt
            rv(k)=rv(k) - cerpt*dt
            tl(k)=tl(k) + convliq(k)*cerpt*dt
          end if
        end if
    end do
  
  end subroutine wtr_dff_SB
  subroutine auto_SB(n1,dn0,rc,rp,np,tl,diss)
  !
  ! ---------------------------------------------------------------------
  ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
  ! drizzle drops due to autoconversion. The autoconversion rate assumes
  ! f(x)=A*x**(nu_c)*exp(-Bx), an exponential in drop MASS x. It can
  ! be reformulated for f(x)=A*x**(nu_c)*exp(-Bx**(mu)), where formu=1/3
  ! one would get a gamma dist in drop diam -> faster rain formation.
  !

    integer, intent (in) :: n1
    real, intent (in)    :: dn0(n1), diss(n1)
    real, intent (inout) :: rc(n1), rp(n1), np(n1),tl(n1)

    real            :: nu_c  = 0.           ! width parameter of cloud DSD
    real, parameter :: kc_0  = 9.44e+9      ! Long-Kernel
    real, parameter :: k_1  = 6.e+2        ! Parameter for phi function
    real, parameter :: k_2  = 0.68         ! Parameter for phi function

    real, parameter :: Cau = 4.1e-15 ! autoconv. coefficient in KK param.
    real, parameter :: Eau = 5.67    ! autoconv. exponent in KK param.
    real, parameter :: mmt = 1.e+6   ! transformation from m to \mu m

    real, parameter :: kc_a1 = 0.00489
    real, parameter :: kc_a2 = -0.00042
    real, parameter :: kc_a3 = -0.01400
    real, parameter :: kc_b1 = 11.45
    real, parameter :: kc_b2 = 9.68
    real, parameter :: kc_b3 = 0.62
    real, parameter :: kc_c1 = 4.82
    real, parameter :: kc_c2 = 4.80
    real, parameter :: kc_c3 = 0.76
    real, parameter :: kc_bet = 0.00174
    real, parameter :: csx = 0.23
    real, parameter :: ce = 0.93

    integer :: k
    real    :: k_au, Xc, Dc, au, tau, phi, kc_alf, kc_rad, kc_sig, k_c, Re, epsilon, l
    real, dimension(n1) :: Resum, Recnt

    !
    ! Calculate the effect of turbulence on the autoconversion/
    ! accretion rate using equation (6) in Seifert et al. (2009)
    !


    do k=2,n1-1
        Xc = rc(k)/(cldw%nr+eps0)
        if (Xc > 0.) then
          k_c = kc_0
          nu_c = cldw%nu1*dn0(k)*rc(k)+cldw%nu2

          if (turbulence) then
            kc_alf = ( kc_a1 + kc_a2 * nu_c )/ ( 1. + kc_a3 * nu_c )
            kc_rad = ( kc_b1 + kc_b2 * nu_c )/ ( 1. + kc_b3 * nu_c )
            kc_sig = ( kc_c1 + kc_c2 * nu_c )/ ( 1. + kc_c3 * nu_c )
            call azero(n1,Recnt,Resum)
            Dc = ( Xc / prw )**(1./3.)  ! mass mean diameter cloud droplets in m
              !
              ! Calculate the mixing length, dissipation rate and Taylor-Reynolds number
              !
              l = csx*((1/dxi)*(1/dyi)*(1/dzt(k)))**(1./3.)
              epsilon = min(diss(k),0.06)
              Re = (6./11.)*((l/ce)**(2./3))*((15./(1.5e-5))**0.5)*(epsilon**(1./6) )
              !
              ! Dissipation needs to be converted to cm^2/s^3, which explains the factor
              ! 1.e4 with which diss(k) is multiplied
              !
              k_c = k_c * (1. + epsilon *1.e4* (Re*1.e-3)**(0.25) &
                    * (kc_bet + kc_alf * exp( -1.* ((((Dc/2.)*1.e+6-kc_rad)/kc_sig)**2) )))
              !print *,'enhancement factor = ', k_c/(9.44e+9)
              !
              ! Calculate conditional average of Re i.e., conditioned on cloud/rain water
              !
              Resum(k) = Resum(k)+ Re
              Recnt(k) = Recnt(k)+ 1
          end if

          k_au  = k_c / (20.*cldw%x_max) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2
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
          if (khairoutdinov) then
              Dc = ( Xc / prw )**(1./3.)
              au = Cau * (Dc * mmt / 2.)**Eau
          end if

          rp(k) = rp(k) + au*dt
          rc(k) = rc(k) - au*dt
          tl(k) = tl(k) + convliq(k)*au*dt
          np(k) = np(k) + au/cldw%x_max*dt
          !
        end if
    end do

  end subroutine auto_SB
  subroutine accr_SB(n1,dn0,rc,rp,np,tl,diss)
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

    real, parameter :: k_r0 = 4.33
    real, parameter :: k_1 = 5.e-5
    real, parameter :: Cac = 67.     ! accretion coefficient in KK param.
    real, parameter :: Eac = 1.15    ! accretion exponent in KK param.

    integer :: k
    real    :: tau, phi, ac, sc, k_r, epsilon

      do k=2,n1-1
          if (rc(k) > 0. .and. rp(k) > 0.) then
            tau = 1.0-rc(k)/(rc(k)+rp(k)+eps0)
            tau = MIN(MAX(tau,eps0),1.)
            phi = (tau/(tau+k_1))**4

            k_r = k_r0
            !
            ! Simulate the effect of turbulence on the collision kernel
            ! (dissipation needs to be converted to cm^2/s^3)
            !
            if (turbulence) then
                epsilon = min(600.,diss(k)*1.e4)   ! put an upper limit to the dissipation rate used
                k_r = k_r*(1+(0.05*epsilon**0.25))
            end if

            ac  = k_r * rc(k) * rp(k) * phi * sqrt(rho_0*dn0(k))
            !
            ! Khairoutdinov and Kogan
            !
            if (khairoutdinov) then
            ac = Cac * (rc(k) * rp(k))**Eac
            end if
            !
            rp(k) = rp(k) + ac*dt
            rc(k) = rc(k) - ac*dt
            tl(k) = tl(k) + convliq(k)*ac*dt

          end if
          sc = k_r * np(k) * rp(k) * sqrt(rho_0*dn0(k))
          np(k) = np(k) - sc*dt
      end do

  end subroutine accr_SB
  subroutine sedim_rd(n1,dt,dn0,rp,np)
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
     real, parameter :: a2 = 9.65       ! in SI [m/s]
     real, parameter :: c2 = 6e2        ! in SI [1/m]
     real, parameter :: Dv = 25.0e-6    ! in SI [m/s]
     real, parameter :: cmur1 = 10.0    ! mu-Dm-relation for rain following
     real, parameter :: cmur2 = 1.20e+3 ! Milbrandt&Yau 2005, JAS, but with
     real, parameter :: cmur3 = 1.5e-3  ! revised constants
     real, parameter :: aq = 6.0e3
     real, parameter :: bq = -0.2
     real, parameter :: an = 3.5e3
     real, parameter :: bn = -0.1


     integer ::  k, kp1, kk, km1
     real    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz
     real, dimension(n1) :: nslope,rslope,dn,dr, rfl, nfl, vn, vr, cn, cr

     b2 = a2*exp(c2*Dv)

        nfl(n1) = 0.
        rfl(n1) = 0.
        do k=n1-1,2,-1
          Xp = rp(k) / (np(k)+eps0)
          !
          ! Adjust Dm and mu-Dm and Dp=1/lambda following Milbrandt & Yau
          !
          Dm = ( 6. / (rowt*pi) * Xp )**(1./3.)
          mu = cmur1*(1.+tanh(cmur2*(Dm-cmur3)))
          Dp = (Dm**3/((mu+3.)*(mu+2.)*(mu+1.)))**(1./3.)

          vn(k) = sqrt(dn0(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(1.+mu)))
          vr(k) = sqrt(dn0(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(4.+mu)))
          !
          ! Set fall speeds following Khairoutdinov and Kogan

!              if (khairoutdinov) then
!                 vn(k) = max(0.,an * Dp + bn)
!                 vr(k) = max(0.,aq * Dp + bq)
!              end if
!irina-olivier
          if (khairoutdinov) then
              vn(k) = max(0.,an * Dm + bn)
              vr(k) = max(0.,aq * Dm + bq)
          end if

        end do

        do k=2,n1-1
          kp1 = min(k+1,n1-1)
          km1 = max(k,2)
          cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzt(k)*dt
          cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzt(k)*dt
        end do

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

        rfl(n1-1) = 0.
        nfl(n1-1) = 0.
        do k=n1-2,2,-1

          kk = k
          tot = 0.0
          zz  = 0.0
          cc  = min(1.,cn(k))
          do while (cc > 0 .and. kk <= n1-1)
              tot = tot + dn0(kk)*(np(kk)+nslope(kk)*(1.-cc))*cc/dzt(kk)
              zz  = zz + 1./dzt(kk)
              kk  = kk + 1
              cc  = min(1.,cn(kk) - zz*dzt(kk))
          enddo
          nfl(k) = -tot /dt

          kk = k
          tot = 0.0
          zz  = 0.0
          cc  = min(1.,cr(k))
          do while (cc > 0 .and. kk <= n1-1)
              tot = tot + dn0(kk)*(rp(kk)+rslope(kk)*(1.-cc))*cc/dzt(kk)
              zz  = zz + 1./dzt(kk)
              kk  = kk + 1
              cc  = min(1.,cr(kk) - zz*dzt(kk))
          enddo
          rfl(k) = -tot /dt

          kp1=k+1
          flxdiv = (rfl(kp1)-rfl(k))*dzt(k)/dn0(k)
          rp(k) = rp(k)-flxdiv*dt
          np(k) = np(k)-(nfl(kp1)-nfl(k))*dzt(k)/dn0(k)*dt

        end do

   end subroutine sedim_rd
  subroutine sedim_cd(n1,dt,tl,th,tk,rc)
  !
  ! ---------------------------------------------------------------------
  ! SEDIM_CD: calculates the cloud-droplet sedimentation flux and its effect
  ! on the evolution of r_t and theta_l assuming a log-normal distribution
  !

    integer, intent (in):: n1
    real, intent (in)                        :: dt
    !irina
    real, intent (in),   dimension(n1) :: th,tk
    !irina
    !
    real, intent (inout),dimension(n1) :: rc,tl

    real, parameter :: c = 1.19e8 ! Stokes fall velocity coef [m^-1 s^-1]

    integer ::  k, kp1
    real    :: Dc, Xc, vc, flxdiv
    real    :: rfl(n1)


    !
    ! calculate the precipitation flux and its effect on r_t and theta_l
    !
          rfl(n1) = 0.
          do k=n1-1,2,-1
             Xc = rc(k) / (cldw%nr+eps0)
             Dc = ( Xc / prw )**(1./3.)
             vc = min(c*(Dc*0.5)**2 * exp(5*(log(1.3))**2),1./(dzt(k)*dt))
             rfl(k) = - rc(k) * vc
             !
             kp1=k+1
             flxdiv = (rfl(kp1)-rfl(k))*dzt(k)
             rc(k) = rc(k)-flxdiv*dt
             tl(k) = tl(k)+flxdiv*(alvl/cp)*th(k)/tk(k)*dt
          end do

  end subroutine sedim_cd
  subroutine n_icenuc(n1,nin,tk,qsup)
    integer, intent(in) :: n1
    real, intent(inout), dimension (n1) :: nin
    real, intent(in) , dimension(n1) :: tk,qsup
    integer :: k
    real :: fact
    real, parameter :: t_nuc = 60.
    fact = (1-exp(-dt/t_nuc))
    do k=2,n1
      nin(k)  = nin(k) + (n_ice_meyers_contact(tk(k),min(qsup(k),0.25))-nin(k))*fact
    end do
  end subroutine n_icenuc

  subroutine ice_nucleation(n1,nin,qice,nice,qsup,tl,tk)
    integer,intent(in) :: n1
    real, intent(inout),dimension(n1) :: nin,nice, qice, tl
    real,intent(in), dimension(n1) :: tk,qsup
    real :: nuc_n, nuc_q
    integer :: k
    
    do k=2,n1
      if (tk(k) < t_nuc .and. qsup(k) > 0.0) then
        nuc_n = 0d0
        nuc_n = max(nin(k) - (nice(k)),0.)

        nuc_q = min(nuc_n * ice%x_min, qsup(k))
        !nuc_n = nuc_q / ice%x_min                !axel 20040416
        nin(k)  = nin(k)  - nuc_n
        nice(k) = nice(k) + nuc_n
        qice(k) = qice(k) + nuc_q
        tl(k)   = tl(k)   + convice(k)*nuc_q
      endif
    end do
  end subroutine
  subroutine cloud_freeze (n1,qcloud,qice,nice,tl,tk)
!
    integer, intent(in) :: n1
    real, intent(inout), dimension(n1) :: qice,nice,qcloud,tl
    real, intent(in), dimension(n1) :: tk
!     real, intent(in) :: nc
!         ! .. local variables ..
        integer          ::  k
        real :: frq, frn,qc,xc,jhet,jhom,jtot,tc
        real, parameter :: ahet = 6.5d-1 ! 1/k,      messung nach barklie and gokhale
        real, parameter :: bhet = 2.0d+2 ! 1/(m3 s), messung nach barklie and gokhale
        real            :: facg


        facg = moment_gamma(cldw,2)    ! <hn

        !..test auf schmelzen oder gefrieren von wolkenteilchen
        do k = 1, n1
          if (tk(k) < tmelt) then
            qc = qcloud(k)
            tc = tk(k) - tmelt
            if (tc < -50.0) then
              frq = qc                                               !..komplettes hom. gefrieren
              frn = cldw%nr                                               !..unterhalb -50 c
            else
              xc = min(max(qc/(cldw%nr+eps0),cldw%x_min),cldw%x_max)    !..mittlere masse

              !..hom. gefrieren nach jeffrey und austin (1997), siehe auch cotton und field (2001)
              if (tc > -30.0) then
                jhom = 1.0d6 * exp(-7.63-2.996*(tc+30.0))           !..j in 1/(m3 s)
              else
                jhom = 1.0d6 * exp(-243.4-14.75*tc-0.307*tc**2-0.00287*tc**3-0.0000102*tc**4)
              endif

              !..het. gefrieren: stochastisches modell nach bigg, daten nach barklie and gokhale
              jhet = bhet * ( exp( - ahet * tc) - 1.0 )            !..j in 1/(m3 s)
              !jhet = 0.0 ! neglected for cloud droplets

              !..umwandlungsraten fuer anzahl- und massendichte
              jtot = (jhom + jhet) / rowt * dt                     !..j*dt in 1/kg
              frn  = jtot * qc
              frq  = jtot * qc * xc * facg
            end if
            frq  = min(frq,qc)

            !..berechnung der h2o-komponenten
            qcloud(k) = qcloud(k) - frq
            tl(k) = tl(k)+ (convice(k)-convliq(k))*frq
!                 ncloud(k) = ncloud(k) - frn

            frn  = max(frn,frq/cldw%x_max)
            qice(k)   = qice(k)   + frq
            nice(k)   = nice(k)   + frn

            ! ub<<

          end if
        end do

      end subroutine cloud_freeze

      subroutine rain_freeze (n1,qrain,nrain,qice,nice,qgrp,tl,tk)
        integer, intent(in) :: n1
        real, dimension(n1), intent(inout) :: qrain,nrain,qice,nice,qgrp,tl
        real, dimension(n1), intent(in) :: tk
        ! .. local variables ..
        integer                     :: k
        real            :: fr_q,fr_n,q_r,x_r,n_r,j_het,&
            &  fr_q_i,fr_n_i,fr_q_g,fr_n_g,n_0,lam,xmax_ice,fr_q_tmp,fr_n_tmp
        real, save      :: coeff_z
        real, parameter :: a_het = 6.5d-1 ! messung nach barklie and gokhale (pk s.350)
        real, parameter :: b_het = 2.0d+2 ! messung nach barklie and gokhale (pk s.350)
        real, parameter :: q_crit_fr = 1.000d-6 ! q-schwellenwert fuer rain_freeze
        real, parameter :: d_rainfrz_ig = 0.50d-3 !  rain --> ice oder graupel
!         real, parameter :: d_rainfrz_gh = 1.25d-3 ! rain --> graupel oder hail

! 
! 
!           !..koeff. fuer reflektivitaet z (2. moment)
          coeff_z = moment_gamma(rain,2)
! 
        xmax_ice = (d_rainfrz_ig/rain%a_geo)**(1.0d0/rain%b_geo)
!         xmax_gr  = (d_rainfrz_gh/rain%a_geo)**(1.0d0/rain%b_geo)
        fr_q_g = 0.0
        !..test auf schmelzen oder gefrieren von regentropfen
        do k = 1, n1
              q_r = qrain(k)
              n_r = nrain(k)

              if (tk(k) < tmelt) then
                if (q_r <= q_crit_fr) then
                  if (tk(k) < t_hn) then
                    fr_q = q_r                  !  ausfrieren unterhalb t_hn \approx -40 c
                    fr_n = n_r
                    fr_n_i= n_r
                    fr_q_i= q_r
                    fr_n_g= 0.0
                    fr_q_g= 0.0
    ! ub_20080220>
                    fr_n_tmp = 1.0
                    fr_q_tmp = 1.0
    ! <ub_20080220
                  else
                    fr_q = 0.0
                    fr_n = 0.0
                    fr_n_i= 0.0
                    fr_q_i= 0.0
                    fr_n_g= 0.0
                    fr_q_g= 0.0
    ! ub_20080220>
                    fr_n_tmp = 0.0
                    fr_q_tmp = 0.0
    ! <ub_20080220
                  end if
                else
                  x_r = min(max(q_r/(n_r+eps0),rain%x_min),rain%x_max)
                  n_r = q_r / x_r
                  if (tk(k) < t_hn) then            !..nur eis
                    fr_q = q_r                  !  ausfrieren unterhalb t_hn \approx -40 c
                    fr_n = n_r

                    ! ub>> je nach groesse werden die gefrorenen regentropfen dem wolkeneis zugeschlagen
                    !      oder dem graupel oder hagel. hierzu erfolgt eine partielle integration des spektrums von 0
                    !      bis zu einer ersten trennmasse xmax_ice (--> eis), von dort bis zu xmax_gr (--> graupel)
                    !      und von xmax_gr bis unendlich (--> hagel).

                    lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                    n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                    fr_n_i = n_0/(rain%mu*lam**((rain%nu+1.0)/rain%mu))* &
                        incgfct_lower((rain%nu+1.0)/rain%mu, lam*xmax_ice**rain%mu)
                    fr_q_i = n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                        incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)

                    fr_n_g = fr_n - fr_n_i
                    fr_q_g = fr_q - fr_q_i
                    fr_n_tmp = n_r/max(fr_n,n_r)
                    fr_q_tmp = q_r/max(fr_q,q_r)

                  else                           !..heterogenes gefrieren
                    j_het = max(b_het * ( exp( a_het * (tmelt - tk(k))) - 1.0 ),0.d0) / rowt * dt

                    ! ub>> je nach groesse werden die gefrorenen regentropfen dem wolkeneis zugeschlagen
                    !      oder dem graupel oder hagel. hierzu erfolgt eine partielle integration des spektrums von 0
                    !      bis zu einer ersten trennmasse xmax_ice (--> eis), von dort bis zu xmax_gr (--> graupel)
                    !      und von xmax_gr bis unendlich (--> hagel).

                    if (j_het >= 1d-20) then
                      fr_n  = j_het * q_r
                      fr_q  = j_het * q_r * x_r * coeff_z

                      lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                      n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                      fr_n_i = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                          incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                      fr_q_i = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                          incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_ice**rain%mu)
                      fr_n_tmp = n_r/max(fr_n,n_r)
                      fr_q_tmp = q_r/max(fr_q,q_r)
                    else
                      fr_n= 0.0
                      fr_q= 0.0
                      fr_n_i= 0.0
                      fr_q_i= 0.0
                      fr_q_g= 0.0
    ! <ub_20080212
                      fr_n_tmp = 0.0
                      fr_q_tmp = 0.0
                    end if

                  end if

                  fr_n = fr_n * fr_n_tmp
                  fr_q = fr_q * fr_q_tmp

                  fr_n_i = fr_n_i * fr_n_tmp
                  fr_q_i = fr_q_i * fr_q_tmp
                  fr_q_g = fr_q_g * fr_q_tmp
                end if

                !..berechnung der h2o-komponenten

                qrain(k) = qrain(k) - fr_q
                tl(k) = tl(k)+ (convice(k)-convliq(k))*fr_q
                nrain(k) = n_r - fr_n

                if (qrain(k) < 0.0) then
                  !write (*,*) 'seifert rain_freeze: qrain < 0.0, ', k, tk(k), q_r, j_het, fr_q, fr_q_tmp
                  qrain(k) = 0.0d0
                end if
                if (nrain(k) < 0.0) then
                  !write (*,*) 'seifert rain_freeze: nrain < 0.0, ', k, tk(k), n_r, j_het, fr_n, fr_n_tmp
                  nrain(k) = 0.0d0
                end if

              qice(k) = qice(k)  + fr_q_i
              nice(k) = nice(k)  + fr_n_i
              qgrp(k) = qgrp(k)  + fr_q_g

              end if
        end do

      end subroutine rain_freeze

      subroutine deposition(n1,meteor,qice,nice,qv,tl,tk,qsup)
      integer, intent(in) :: n1
      type(particle), intent(in) :: meteor
      real, dimension(n1), intent(inout) :: qice,nice,qv,tl
      
      real, dimension(n1), intent(in)    :: tk,qsup

    
      ! locale variablen
      integer                     :: k
      real        :: q_g,n_g,x_g,d_g,v_g,f_v,f_n,n_re,f_v_fakt,vent_fakt
      real        :: c_g                 !..koeff. fuer mittlere kapazitaet
      real        :: a_f,b_f,a_n,b_n     !..koeff. fuer mittleren ventilationkoeff.
      real :: dep,gi
        c_g = 1.0 / meteor%cap
        a_n = vent_coeff_a(meteor,0)
        b_n = vent_coeff_b(meteor,0)
        a_f = vent_coeff_a(meteor,1)
        b_f = vent_coeff_b(meteor,1)

      f_v_fakt = n_sc**n_f
      vent_fakt = b_n / b_f
      do k = 1, n1

        ! hn: in case q_garupel=0, dep_graupel has to be zero too

        if (qice(k)> 0.0) then
          n_g = nice(k)                                     !..anzahldichte
          q_g = qice(k)                                     !..massendichte

          x_g = min(max(q_g/(n_g+eps0),max(meteor%x_min,eps0)),meteor%x_max)  !..mittlere masse
          ! x_g**b_geo durch exp(b_geo*log(x_g)) ersetzen, da x_g == 0 nicht mehr vorkommen kann; 20 % schneller:
          d_g = meteor%a_geo * exp(meteor%b_geo*log(x_g))          !..mittlerer durchmesser
          v_g = meteor%a_vel * exp(meteor%b_vel*log(x_g)) * dn0(k)  !..mittlere sedimentationsgeschw.

          n_re = v_g * d_g / nu_l                                    !..mittlere reynoldszahl
          ! n_re**m_f durch exp(m_f*log(n_re)) ersetzen, da n_re == 0 nicht mehr vorkommen kann:
          f_v  = a_f + b_f * f_v_fakt * exp(m_f*log(n_re))           !..mittlerer vent.koeff.
          ! ub>> schnellere berechnung wg. verzicht auf "**":
          !              f_n  = a_n + b_n * f_v_fakt * n_re**m_f                    !..mittlerer vent.koeff.
          f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer vent.koeff.
          ! ub<<
          f_v  = max(f_v,1.d0) !unnoetig??
          f_n  = max(f_n,1.d0) !unnoetig??
!           e_si = e_es(tk(k))
          gi = 4.0*pi / ( alvi**2 / (K_T * Rm * tk(k)**2) + Rm * tk(k) / (D_v * e_es(tk(k))) )
          dep = gi * n_g * c_g * d_g * f_v * qsup(k) * dt

          qice(k) = qice(k) + dep
          qv(k) = qv(k) - dep
          tl(k) = tl(k) + convice(k)*dep
          if (meteor%moments==2) then
            nice(k) = nice(k) + dep * f_n/f_v / x_g
          end if
        endif

      enddo

    end subroutine deposition


  subroutine melting(n1,meteor,q,thl,nr,qcld,qrain,nrain,tk)
    integer, intent(in) :: n1
    real,dimension(n1),intent(in) :: tk
    type(particle), intent(in) :: meteor
    real,dimension(n1),intent(inout)  ::q,nr,qrain,nrain,qcld,thl

    integer                     :: k

    real            :: x_m,d_m,v_m,t_a,n_re,d_t,e_a
    real            :: melt,melt_v,melt_h,melt_n,melt_q
    real            :: fh_q,fv_q
    real, save      :: a_melt_n,b_melt_n
    real, save      :: a_melt_q,b_melt_q

    real, parameter :: eps   = 1.d-20


      a_melt_n = vent_coeff_a(meteor,0)
      b_melt_n = vent_coeff_b(meteor,0)
      a_melt_q = vent_coeff_a(meteor,1)
      b_melt_q = vent_coeff_b(meteor,1)

    do k = 2,n1
          t_a = tk(k) !wrf!+ t(k) + t_g(k)
          e_a = e_ws(T_a)                                     !..Saettigungsdampfdruck
                   
          if (t_a > tmelt .and. q(k) > 0.0) then

            x_m = min(max(q(k)/(nr(k)+eps),meteor%x_min),meteor%x_max)  !..mittlere qe in si

            D_m = meteor%a_geo * x_m**meteor%b_geo                   !..mittlerer Durchmesser
            v_m = meteor%a_vel * x_m**meteor%b_vel * rho_0  !..mittlere Sedimentationsgeschw.

            N_re = v_m * D_m / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
  !          fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

            D_T  = Kt / (cp * rho_0)!WRF!+rho_g(i,j,k)))
            fh_q = D_T / D_v * fv_q
 !           fh_n = D_T / D_v * fv_n

            melt   = 2.0*pi / (alvi-alvl) * D_m * nr(k) * dt

            melt_h = melt * Kt * (T_a - tmelt)
            melt_v = melt * D_v*alvl/Rm * (e_a/T_a - e_3/tmelt)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_m

! ub>> setzte melt_n so, dass x_h beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q(k)) / x_m + nr(k), 0.0d0), nr(k))

            melt_q = MIN(q(k),melt_q)
            melt_n = MIN(nr(k),melt_n)

            melt_q = MAX(0.d0,melt_q)
            melt_n = MAX(0.d0,melt_n)


            q(k) = q(k) - melt_q
            nr(k) = nr(k) - melt_n
            if(x_m<cldw%x_max) then
              qcld(k)    = qcld(k)    + melt_q
              thl(k)      = thl(k)      - convice(k)*melt_q

            else
              qrain(k)    = qrain(k)    + melt_q
              nrain(k)    = nrain(k)    + melt_n
              thl(k)      = thl(k)      - (convice(k)-convliq(k))*melt_q
            end if
          endif
        enddo



  end subroutine melting

  subroutine auto_ice
  end subroutine auto_ice

! 
  subroutine aggr_ice()
! !     !*******************************************************************************
! !     !                                                                              *
! !     !       berechnung der koagulation von graupel und schnee                      *
! !     !                                                                              *
! !     !*******************************************************************************
! ! 
! !     use globale_variablen,  only: loc_ix, loc_iy, loc_iz, t, p, q, p_g, t_g, rho_g, &
! !          &                        q_graupel, n_graupel, q_ice, n_ice, dt, dqdt, speichere_umwandlungsraten
! !     use konstanten,         only: pi,cv,cp,wolke_typ
! !     use parallele_umgebung, only: isio
! !     use initialisierung,    only: t_0,p_0,rho_0
! ! 
! !     implicit none
! ! 
! !     ! locale variablen
! !     integer                     :: i,j,k
! !     integer, save               :: firstcall
! ! 
! !     double precision            :: t_a
! !     double precision            :: q_g,n_g,x_g,d_g,v_g
! !     double precision            :: q_i,n_i,x_i,d_i,v_i
! !     double precision            :: coll_n,coll_q,e_coll
! !     double precision, save      :: delta_n_gg,delta_n_gi,delta_n_ii
! !     double precision, save      :: delta_q_gg,delta_q_gi,delta_q_ii
! !     double precision, save      :: theta_n_gg,theta_n_gi,theta_n_ii
! !     double precision, save      :: theta_q_gg,theta_q_gi,theta_q_ii
! ! 
! !     double precision, parameter :: eps  = 1.d-20
! !     double precision, parameter :: e_gi = 0.05
! 
!       delta_n_gg = coll_delta_11(graupel,ice,0)
!       delta_n_gi = coll_delta_12(graupel,ice,0)
!       delta_n_ii = coll_delta_22(graupel,ice,0)
!       delta_q_gg = coll_delta_11(graupel,ice,0)
!       delta_q_gi = coll_delta_12(graupel,ice,1)
!       delta_q_ii = coll_delta_22(graupel,ice,1)
! 
!       theta_n_gg = coll_theta_11(graupel,ice,0)
!       theta_n_gi = coll_theta_12(graupel,ice,0)
!       theta_n_ii = coll_theta_22(graupel,ice,0)
!       theta_q_gg = coll_theta_11(graupel,ice,0)
!       theta_q_gi = coll_theta_12(graupel,ice,1)
!       theta_q_ii = coll_theta_22(graupel,ice,1)
! 
!     do k = 1, loc_iz
!       q_i = q_ice(i,j,k)                                     !..fluessigwassergehalt in si
!       q_g = q_graupel(i,j,k)                                 !..fluessigwassergehalt in si
!       if (q_i > q_krit .and. q_g > q_krit) then
!         t_a = t_0(i,j,k) !wrf!+ t_0(i,j,k) + t_g(i,j,k)
! 
!         !.. temperaturabhaengige sticking efficiency nach lin (1983)
!         if (t_a > t_3) then
!           e_coll = 1.0
!         else
!           e_coll = min(exp(0.09*(t_a-t_3)),1.0d0)
!         end if
! 
!         n_i = n_ice(i,j,k)                                        !..anzahldichte
!         n_g = n_graupel(i,j,k)                                    !..anzahldichte
! 
!         x_i = min(max(q_i/(n_i+eps),ice%x_min),ice%x_max)         !..mittlere masse
!         x_g = min(max(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere masse
! 
!         d_i = ice%a_geo * x_i**ice%b_geo                          !..mittlerer durchmesser
!         v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)         !..mittlere sedimentationsgeschw.
! 
!         d_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer durchmesser
!         v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere sedimentationsgeschw.
! 
!         coll_n = pi/4.0 * n_g * n_i * e_coll * dt &
!               &   * (delta_n_gg * d_g**2 + delta_n_gi * d_g*d_i + delta_n_ii * d_i**2) &
!               &   * (theta_n_gg * v_g**2 - theta_n_gi * v_g*v_i + theta_n_ii * v_i**2)**0.5
! 
!         coll_q = pi/4.0 * n_g * q_i * e_coll * dt &
!               &   * (delta_q_gg * d_g**2 + delta_q_gi * d_g*d_i + delta_q_ii * d_i**2) &
!               &   * (theta_q_gg * v_g**2 - theta_q_gi * v_g*v_i + theta_q_ii * v_i**2)**0.5
! 
!         coll_n = min(n_i,coll_n)
!         coll_q = min(q_i,coll_q)
! 
!         q_graupel(i,j,k) = q_graupel(i,j,k) + coll_q
!         q_ice(i,j,k)     = q_ice(i,j,k)     - coll_q
!         n_ice(i,j,k)     = n_ice(i,j,k)     - coll_n
! 
!       endif
!     enddo
! 
  end subroutine aggr_ice

  subroutine sedimentation(n1,meteor,rp,np)
    integer, intent(in) :: n1
    real, dimension(n1), intent(inout) :: rp
    type(particle),  intent(in) :: meteor
    real, dimension(n1), intent(inout),optional :: np
! 
!     ! .. local variables ..
    integer  :: k,kk,km1,kp1
    real      :: alfn,alfq,c_lam,lam,sk,tot,zz,xp
    real, dimension(n1)   :: nfl,rfl,vn,vr,dn,dr,rslope,cn,cr,nslope
    real :: cc, flxdiv,maxi,mini

    alfn = meteor%a_vel * gfct((meteor%nu+meteor%b_vel+1.0)/meteor%mu)&
          &                / gfct((meteor%nu+1.0)/meteor%mu)
    alfq = meteor%a_vel * gfct((meteor%nu+meteor%b_vel+2.0)/meteor%mu)&
          &                / gfct((meteor%nu+2.0)/meteor%mu)
    c_lam = gfct((meteor%nu+1.0)/meteor%mu)/gfct((meteor%nu+2.0)/meteor%mu)

    where (rp < 0.0d0) rp = 0.0d0
    where (np < 0.0d0) np = 0.0d0

    do k=n1-1,2,-1

      Xp = rp(k) / (np(k)+eps0)
      lam = ( c_lam * xp )**(meteor%b_vel)
      vr(k) = alfq * lam
      vr(k) = max(vr(k),0.1d+0)
      vr(k) = min(vr(k),30.d0)
      vr(k) = -vr(k)

      vn(k) = alfn * lam
      vn(k) = max(vn(k),0.1d+0)
      vn(k) = min(vn(k),30.d0)
      vn(k) = -vn(k)
    end do
    do k=2,n1-1
      kp1 = min(k+1,n1-1)
      km1 = max(k,2)
      cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzt(k)*dt
      cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzt(k)*dt
    end do

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

    rfl(n1-1) = 0.
    nfl(n1-1) = 0.
    do k=n1-2,2,-1

      kk = k
      tot = 0.0
      zz  = 0.0
      cc  = min(1.,cn(k))
      do while (cc > 0 .and. kk <= n1-1)
          tot = tot + dn0(kk)*(np(kk)+nslope(kk)*(1.-cc))*cc/dzt(kk)
          zz  = zz + 1./dzt(kk)
          kk  = kk + 1
          cc  = min(1.,cn(kk) - zz*dzt(kk))
      enddo
      nfl(k) = -tot /dt

      kk = k
      tot = 0.0
      zz  = 0.0
      cc  = min(1.,cr(k))
      do while (cc > 0 .and. kk <= n1-1)
          tot = tot + dn0(kk)*(rp(kk)+rslope(kk)*(1.-cc))*cc/dzt(kk)
          zz  = zz + 1./dzt(kk)
          kk  = kk + 1
          cc  = min(1.,cr(kk) - zz*dzt(kk))
      enddo
      rfl(k) = -tot /dt

      kp1=k+1
      flxdiv = (rfl(kp1)-rfl(k))*dzt(k)/dn0(k)
      rp(k) = rp(k)-flxdiv*dt
      if(meteor%moments==2) then
        np(k) = np(k)-(nfl(kp1)-nfl(k))*dzt(k)/dn0(k)*dt
      end if

    end do


  end subroutine sedimentation
  
  real function n_ice_meyers_contact(t_a,s)
    ! diagnostische beziehung fuer anzahldichte der eisteilchen nach meyers (1992)
!
    real, intent(in):: s,t_a
    real, parameter :: n_0 = 1.0d+3
    real, parameter :: n_m = 1.0d+3
    real, parameter :: a_d = -0.639
    real, parameter :: b_d = 12.960
    real, parameter :: c_d = -2.8
    real, parameter :: d_d = 0.262
    real, parameter :: tmelt  = 2.732d+2     !..tripelpunkt wasser

    n_ice_meyers_contact = n_0 * exp( a_d + b_d * s )         &
          &               + n_m * exp( c_d + d_d * (t_a - tmelt) )

!     return
  end function n_ice_meyers_contact
  real function vent_coeff_a(parti,n)
    implicit none

    integer        :: n
    type(particle) :: parti

    vent_coeff_a = parti%a_ven * gfct((parti%nu+n+parti%b_geo)/parti%mu)              &
         &                  / gfct((parti%nu+1.0)/parti%mu)                        &
         &                * ( gfct((parti%nu+1.0)/parti%mu)                        &
         &                  / gfct((parti%nu+2.0)/parti%mu) )**(parti%b_geo+n-1.0)

  end function vent_coeff_a

  real function vent_coeff_b(parti,n)
    implicit none

    integer        :: n
    type(particle) :: parti


    vent_coeff_b = parti%b_ven                                                  &
         & * gfct((parti%nu+n+(m_f+1.0)*parti%b_geo+m_f*parti%b_vel)/parti%mu)  &
         &             / gfct((parti%nu+1.0)/parti%mu)                          &
         &           * ( gfct((parti%nu+1.0)/parti%mu)                          &
         &             / gfct((parti%nu+2.0)/parti%mu)                          &
         &             )**((m_f+1.0)*parti%b_geo+m_f*parti%b_vel+n-1.0)

  end function vent_coeff_b

  real function moment_gamma(p,n)
    integer          :: n
    type(particle)   :: p

    moment_gamma  = gfct((n+p%nu+1.0)/p%mu) / gfct((p%nu+1.0)/p%mu)        &
         &     * ( gfct((  p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) )**n
  end function moment_gamma
  real function gfct(x)
    !*******************************************************************************
    !                                                                              *
    !       gammafunktion aus numerical recipes (f77)                              *
    !                                                                              *
    !*******************************************************************************
    implicit none

    real cof(6)
    real stp,half,one,x,xx,fpf,tmp,ser,gamma
    integer j

    data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,  &
         &     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
    data half,one,fpf/0.5d0,1.0d0,5.5d0/

    xx  = x  - one
    tmp = xx + fpf
    tmp = (xx + half) * log(tmp) - tmp
    ser = one
    do j = 1,6
      xx  = xx  + one
      ser = ser + cof(j) / xx
    enddo
    gamma = tmp + log(stp*ser)
    gamma = exp(gamma)

    gfct = gamma
    return
  end function gfct
    !*******************************************************************************
  !
  ! upper incomplete gamma function
  !
  ! eigentliche obere unvollstaendige gamma-funktion, hier direkt
  ! das integral
  !              int(x)(oo) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  real function incgfct_upper(a,x)

    implicit none

    real, intent(in) :: a, x
    real :: gam, gln

    gam = gammq(a,x,gln)
    incgfct_upper = exp(gln) * gam

  end function incgfct_upper

  !*******************************************************************************
  !
  ! lower incomplete gamma function
  !
  ! eigentliche untere unvollstaendige gamma-funktion, hier direkt
  ! das integral
  !              int(0)(x) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  real function incgfct_lower(a,x)

    implicit none

    real, intent(in) :: a, x
    real :: gam, gln

    gam = gammp(a,x,gln)
    incgfct_lower = exp(gln) * gam

  end function incgfct_lower

  !*******************************************************************************
  !
  ! eigentliche unvollstaendige gamma-funktion, hier direkt
  ! das integral
  !              int(x1)(x2) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  real function incgfct(a,x1,x2)

    implicit none

    real, intent(in) :: a, x1, x2

    incgfct = incgfct_lower(a,x2) - incgfct_lower(a,x1)

  end function incgfct

  real function gammp(a,x,gln)

    implicit none

    real, intent(in) :: a, x
    real, intent(out) :: gln

    real :: gammcf, gamser

    if (x.lt.0.0d0 .or. a .le. 0.0d0) then
      write(*,*) 'error in gammp: bad arguments'
      write(*,*) '  (module gamma_functions, src_seifert.f90)'
      gammp = 0.0d0
      return
    end if

    if (x .lt. a+1.)then
      call gser(gamser,a,x,gln)
      gammp = gamser
    else
      call gcf(gammcf,a,x,gln)
      gammp = 1.0d0 - gammcf
    endif
    return
  end function gammp

  real function gammq(a,x,gln)

    implicit none

    real, intent(in) :: a, x
    real, intent(out) :: gln
    real :: gammcf, gamser

    if (x.lt.0.0d0 .or. a .le. 0.0d0) then
      write(*,*) 'error in gammq: bad arguments (module gamma_functions, src_seifert.f90)'
      gammq = 0.0d0
      return
    end if

    if (x.lt.a+1.) then
      call gser(gamser,a,x,gln)
      gammq = 1.0d0 - gamser
    else
      call gcf(gammcf,a,x,gln)
      gammq = gammcf
    endif
    return
  end function gammq

  real function gammln(x)
  !*******************************************************************************
  !                                                                              *
  !       log(gamma function) taken from press et al.,  numerical recipes (f77)
  !
  !       log(gammafunktion) aus numerical recipes (f77)                         *
  !       (original)                                                             *
  !*******************************************************************************
    implicit none

    real, intent(in) :: x
    real, save :: cof(6), stp
    real :: xx,tmp,ser
    integer :: j
    data cof /76.18009172947146d0,-86.50532032941677d0, &
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5/
    data stp /2.5066282746310005d0/

    xx  = x
    tmp = xx + 5.5d0
    tmp = (xx + 0.5d0) * log(tmp) - tmp
    ser = 1.000000000190015d0
    do j = 1,6
       xx  = xx  + 1.0d0
       ser = ser + cof(j) / xx
    enddo
    gammln = tmp + log(stp*ser/x)
    return
  end function gammln

  real function gfct2(x)
  !*******************************************************************************
  !                                                                              *
  !       gamma function taken from press et al.,  numerical recipes (f77)
  !
  !       gammafunktion aus numerical recipes (f77)                              *
  !       (etwas umformuliert, aber dieselben ergebnisse wie obige originalfunktion)
  !*******************************************************************************
    implicit none

    real, intent(in) :: x
    real, save :: cof(6), stp, half, one, fpf
    real :: xx,tmp,ser,gamma
    integer j

    data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,  &
          &     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
    data half,one,fpf/0.5d0,1.0d0,5.5d0/

    xx  = x  - one
    tmp = xx + fpf
    tmp = (xx + half) * log(tmp) - tmp
    ser = one
    do j = 1,6
       xx  = xx  + one
       ser = ser + cof(j) / xx
    enddo
    gamma = tmp + log(stp*ser)
    gamma = exp(gamma)

    gfct2 = gamma
    return
  end function gfct2
  
  subroutine gcf(gammcf,a,x,gln)

    implicit none

    integer, parameter :: itmax = 100
    real, parameter :: eps = 3.d-7, fpmin = 1.d-30
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
      b=b+2.0d0
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
      write (*,*) 'error in gcf: a too large, itmax too small (module gamma_functions, src_seifert.f90)'
      gammcf = 0.0d0
      return
    end if

    gammcf=exp(-x+a*log(x)-gln)*h

    return
  end subroutine gcf

  subroutine gser(gamser,a,x,gln)

    implicit none

    integer, parameter :: itmax = 100
    real, parameter :: eps=3.d-7
    real, intent(in) :: a, x
    real, intent(out) :: gamser, gln

    integer :: n
    real :: ap,del,sum

    gln=gammln(a)
    if (x.le.0.) then
      if (x.lt.0.) then
        write (*,*) 'error in gser: x < 0 (module gamma_functions, src_seifert.f90)'
      end if
      gamser=0.0d0
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
      write (*,*) 'error in gser: a too large, itmax too small'
      write (*,*) '  (module gamma_functions, src_seifert.f90)'
      gamser = 0.0d0
      return
    end if

    gamser = sum*exp(-x+a*log(x)-gln)

    return
  end subroutine gser

  real function e_es (t_)
    !*******************************************************************************
    !                        saettigungsdampfdruck ueber eis                       *
    !*******************************************************************************

    real, intent (in) :: t_

    e_es  = e_3 * exp (a_e * (t_ - tmelt) / (t_ - b_e))

  end function e_es

  real function e_ws (t_)
    !*******************************************************************************
    !                      saettigungsdampfdruck ueber wasser                      *
    !*******************************************************************************

    real, intent (in) :: t_

    e_ws  = e_3 * exp (a_w * (t_ - tmelt) / (t_ - b_w))

  end function e_ws

  function gen_sequence(nprocess)
  integer, intent(in) :: nprocess
  integer, dimension(nprocess) :: gen_sequence
  integer :: i
    gen_sequence = (/(i,i=1,nprocess)/)
    if (.not. droplet_sedim) then
      where (gen_sequence==isedimcd)
        gen_sequence = 0
      end where
    end if

  end function gen_sequence

  subroutine resetvar(meteor,mass,number)
    type(particle),intent(in)        :: meteor
    real, dimension(:), intent(inout) :: mass
    real, dimension(:), intent(inout), optional :: number
    mass = max(0.,mass)
    if (present(number)) then
      number = max(min(number,mass/meteor%x_min),mass/meteor%x_max)
    end if
  end subroutine resetvar
  subroutine initmcrp(level)
  integer , intent(in) :: level

  firsttime = .false.
  nprocess = 0
  if (droplet_sedim) nprocess = nprocess + 1
  if (level==3)      nprocess = 4
  if (level==4)      nprocess = nprocess + 0

  allocate(seq(nprocess))

!Cloudwater properties
  cldw%nr        = CCN
  cldw%x_min     = 4.2e-15
  cldw%x_max     = 2.6e-10
!Rainwater properties
  rain%x_min     = 5.e-10
  rain%x_max     = 2.6e-10
!Cloud ice properties
!Snow properties
!Graupel properties
  
  end subroutine initmcrp
end module mcrp
