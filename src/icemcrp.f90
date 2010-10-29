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

  use defs, only : tmelt,alvl, alvi,rowt,roice, pi, Rm, cp,t_hn 
  use grid, only : dt,nstep,rkbeta,rkalpha, dxi, dyi ,dzt, nxp, nyp, nzp, a_pexnr, pi0,pi1,a_rp, a_tp, th00, ccn,    &
       dn0, pi0,pi1, a_rt, a_tt,a_rpp, a_rpt, a_npp, a_npt, vapor, liquid,       &
       a_theta, a_scr1, a_scr2, a_scr7, precip, &
       a_ninuct,a_ricet,a_nicet,a_rsnowt,a_rgrt,&
       a_ninucp,a_ricep,a_nicep,a_rsnowp,a_rgrp,rsup
  use thrm, only : thermo, fll_tkrs
  use stat, only : sflg, updtst
  use util, only : get_avg3, azero
  implicit none
 
  logical, parameter :: droplet_sedim = .False., khairoutdinov = .False., turbulence = .False.,ice_multiplication = .TRUE.
  integer            :: nprocess,nprocwarm=5,nprocice=18
  
  integer,parameter  :: iwtrdff = 1,iauto = 2,iaccr = 3,isedimrd = 4,isedimcd = 5, &
                        iicenucnr = 6, iicenuc = 7, ifreez=8, idep=9, imelt_ice=10,imelt_snow=11,imelt_grp=12, ised_ice=13, &
                        iself_ice=14, icoll_ice_snow=15, icoll_ice_grp=16, icoll_snow_grp=17, &
                        iriming_ice_cloud=18, iriming_snow_cloud=19, iriming_grp_cloud=20, &
                        iriming_ice_rain=21, iriming_snow_rain=22,iriming_grp_rain=23
  real, dimension(:),allocatable :: convice,convliq
  integer, dimension(23) :: microseq = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,22/)
  logical :: lrandommicro = .false.
  real :: timenuc = 60
  real :: nin_set = 1.7e3
  ! 
  ! drop sizes definition is based on vanZanten (2005)
  ! cloud droplets' diameter: 2-50 e-6 m
  ! drizzle drops' diameter: 50-1000 e-6 m
  !
  real, parameter    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
  real, parameter    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]

  real, parameter :: eps0 = 1e-20       ! small number
  real, parameter :: eps1 = 1e-9        ! small number
  real, parameter :: qthres = 1e-14        ! small number
  real, parameter :: rho_0 = 1.21       ! air density at surface

  real, parameter :: prw = pi * rowt / 6.
  real, parameter :: nu_l = 1.460e-5     !..kinem. visc. von luft
  real, parameter :: d_v  = 3.000e-5     !..diffusivitaet von wasserdampf
  real, parameter :: k_t  = 2.500e-2     !..waermeleitfaehigkeit

  real, parameter :: n_sc = 0.710        !..schmidt-zahl (pk, s.541)
  real, parameter :: n_f  = 0.333        !..exponent von n_sc im vent-koeff. (pk, s.541)
  real, parameter :: m_f  = 0.500        !..exponent von n_re im vent-koeff. (pk, s.541)
! 
  real, parameter :: a_e  = 2.18745584e1 !..konst. saettigungsdamppfdruck - eis
  real, parameter :: a_w  = 1.72693882e1 !..konst. saettigungsdamppfdruck - wasser
  real, parameter :: b_e  = 7.66000000e0 !..konst. saettigungsdamppfdruck - eis
  real, parameter :: b_w  = 3.58600000e1 !..konst. saettigungsdamppfdruck - wasser
  real, parameter :: e_3  = 6.10780000e2 !..saettigungsdamppfdruck bei t = tmelt
! 
  real, parameter :: c_mult     = 3.5d8    !..koeff. fuer splintering
  real, parameter :: t_mult_min = 265.0    !..minimale temp. splintering
  real, parameter :: t_mult_max = 270.0    !..maximale temp. splintering
  real, parameter :: t_mult_opt = 268.0    !..optimale temp. splintering

  real, parameter :: e_ic  = 0.80              !..max. eff. fuer ice_cloud_riming
  real, parameter :: e_sc  = 0.80              !..max. eff. fuer snow_cloud_riming
  real, parameter :: e_gc  = 1.00              !..max. eff. fuer graupel_cloud_riming
  real, parameter :: alpha_spacefilling = 0.1  !..raumerfuellungskoeff (max. 0.68)

!   real, parameter :: q_crit_ic = 1.000e-5 ! q-schwellenwert fuer ice_cloud_riming
!   real, parameter :: d_crit_ic = 150.0e-6 ! e-schwellenwert fuer ice_cloud_riming
!   real, parameter :: q_crit_sc = 1.000e-5 ! q-schwellenwert fuer snow_cloud_riming
!   real, parameter :: d_crit_sc = 150.0e-6 ! e-schwellenwert fuer snow_cloud_riming
!   real, parameter :: q_crit_ir = 1.000e-5 ! q-schwellenwert fuer ice_rain_riming
!   real, parameter :: d_crit_ir = 100.0e-6 ! e-schwellenwert fuer ice_rain_riming
!   real, parameter :: q_crit_sr = 1.000e-5 ! q-schwellenwert fuer snow_rain_riming
!   real, parameter :: d_crit_sr = 100.0e-6 ! e-schwellenwert fuer snow_rain_riming
!   real, parameter :: q_crit_gc = 1.000e-6 ! q-schwellenwert fuer graupel_cloud_riming
!   real, parameter :: d_crit_gc = 100.0e-6 ! e-schwellenwert fuer graupel_cloud_riming
!   real, parameter :: q_crit_c  = 1.000e-6 ! q-schwellenwert sonst
!   real, parameter :: q_crit    = 1.000e-9 ! q-schwellenwert sonst
!   real, parameter :: d_crit_c  = 10.00e-6 ! e-schwellenwert fuer cloud_collection
!   real, parameter :: d_crit_r  = 10.00e-5 ! e-schwellenwert fuer cloud_collection
!   real, parameter :: d_coll_c  = 40.00e-6 ! oberer wert fuer cloud_coll_eff
!   real, parameter :: q_crit_is = 1.000e-4 ! q-schwellenwert fuer ice_selfcollection
!   real, parameter :: d_crit_is = 50.00e-6 ! e-schwellenwert fuer ice_selfcollection
! 
!   real, parameter :: d_conv_is = 75.00e-6 ! e-schwellenwert fuer ice_selfcollection
!   real, parameter :: d_conv_sg = 200.0e-5 ! e-schwellenwert
!   real, parameter :: d_conv_ig = 200.0e-6 ! e-schwellenwert

  real, parameter :: q_crit_ic = 1.000e-5 ! q-schwellenwert fuer ice_cloud_riming
  real, parameter :: d_crit_ic = 150.0e-6 ! e-schwellenwert fuer ice_cloud_riming
  real, parameter :: q_crit_ir = 1.000e-5 ! q-schwellenwert fuer ice_rain_riming
  real, parameter :: d_crit_ir = 100.0e-6 ! e-schwellenwert fuer ice_rain_riming
  real, parameter :: q_crit_sc = 1.000e-5 ! q-schwellenwert fuer snow_cloud_riming
  real, parameter :: d_crit_sc = 150.0e-6 ! e-schwellenwert fuer snow_cloud_riming
  real, parameter :: q_crit_sr = 1.000e-5 ! q-schwellenwert fuer snow_rain_riming
  real, parameter :: d_crit_sr = 100.0e-6 ! e-schwellenwert fuer snow_rain_riming
  real, parameter :: q_crit_gc = 1.000e-6 ! q-schwellenwert fuer graupel_cloud_riming
  real, parameter :: d_crit_gc = 100.0e-6 ! e-schwellenwert fuer graupel_cloud_riming
  real, parameter :: q_crit_c  = 1.000e-6 ! q-schwellenwert sonst
  real, parameter :: q_crit    = 1.000e-9 ! q-schwellenwert sonst
  real, parameter :: d_conv_sg = 200.0e-5 ! e-schwellenwert
  real, parameter :: d_conv_ig = 200.0e-6 ! e-schwellenwert
  real, parameter :: d_crit_c  = 10.00e-6 ! e-schwellenwert fuer cloud_collection
  real, parameter :: d_crit_r  = 10.00e-6 ! e-schwellenwert fuer cloud_collection
  real, parameter :: d_coll_c  = 40.00e-6 ! oberer wert fuer cloud_coll_eff
  real, parameter :: q_crit_is = 1.000e-4 ! q-schwellenwert fuer ice_selfcollection
  real, parameter :: d_crit_is = 50.00e-6 ! e-schwellenwert fuer ice_selfcollection
  real, parameter :: d_conv_is = 75.00e-6 ! e-schwellenwert fuer ice_selfcollection

  type particle
    character(10) :: name
    integer :: moments
    real :: nr
    real :: nu
    real :: mu
    real :: x_max !..maximum particle mass
    real :: x_min !..minimum particle mass
    real :: a_geo !..koeff. geometrie
    real :: b_geo !..koeff. geometrie = 1/3
    real :: a_vel !..koeff. fallgesetz
    real :: b_vel !..koeff. fallgesetz
    real :: s_vel !..Dispersion der Fallgeschw.
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
    call fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,dn0,th00,a_scr1,rs=a_scr2)
    select case (level) 
    case(2)
!        if (droplet_sedim) call sedim_cd(level,nzp,nxp,nyp,dn0,dt,a_theta,a_scr1,      &
!             liquid,precip,a_rt,a_tt)     
    case(3)
       call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_scr1,vapor,a_scr2,liquid,a_rpp, &
              a_npp,precip,a_rt,a_tt,a_rpt,a_npt,a_scr7)
    case(4)
       call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_scr1,vapor,a_scr2,liquid,a_rpp, &
              a_npp,precip,a_rt,a_tt,a_rpt,a_npt,a_scr7,&
              a_ninuct,a_ricet,a_nicet,a_rsnowt,a_rgrt,&
              a_ninucp,a_ricep,a_nicep,a_rsnowp,a_rgrp,rsup)
    end select

  end subroutine micro
  !
  ! ---------------------------------------------------------------------
  ! MCRPH: calls microphysical parameterization 
  !

  subroutine mcrph(level,n1,n2,n3,dn0,exner,pi0,pi1,thl,tk,rv,rs,rc,rp,np,         &
       rrate, rtt,tlt,rpt,npt,dissip,    &
       ninuct,ricet,nicet,rsnowt,rgrt,ninucp,ricep,nicep,rsnowp,rgrp,rsup)

    integer, intent (in) :: level,n1,n2,n3
    real, dimension(n1)      , intent (in)             :: dn0,pi0,pi1  !Density and mean pressures
    real, dimension(n1,n2,n3), intent (inout)          :: exner, & !exner function
                                                          thl,& !Theta_l
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
    real, dimension(n1,n2,n3), intent (inout),optional :: ninuct,ricet,nicet,rsnowt,rgrt,&
                                                          ninucp,ricep,nicep,rsnowp,rgrp,&
                                                          rsup

    real, dimension(n1) :: tl,temp,qv,qc,nrain,qrain,ninuc,nice, nin_active,si,qice,nsnow,qsnow,ngrp,qgrp,qsat,qsup,q1 !dummy-arrays for conversion between mixing rate
    integer :: i, j,n,k
    logical,save  :: firsttime = .true.
    real :: dtrk
    rrate = 0.

    if(firsttime) call initmcrp(level,firsttime)
    allocate(convice(n1),convliq(n1))

    do j=3,n3-2
      do i=3,n2-2
        call resetvar(rain,rp(1:n1,i,j),np(1:n1,i,j))
        tl = thl(1:n1,i,j)
        temp = tk(1:n1,i,j)
        qv = rv(1:n1,i,j)/dn0
        qc = rc(1:n1,i,j)/dn0
        qsat = rs(1:n1,i,j)/dn0
        qrain = rp(1:n1,i,j)/dn0
        nrain = np(1:n1,i,j)
        convliq = alvl/cp*(pi0+pi1+exner(1:n1,i,j))/cp
        if (level == 4) then
          call resetvar(ice,ricep(1:n1,i,j),nicep(1:n1,i,j))
          call resetvar(snow,rsnowp(1:n1,i,j))
          call resetvar(graupel,rgrp(1:n1,i,j))

          qsup = rsup(1:n1,i,j)/dn0

          ninuc = ninucp(1:n1,i,j)

          qice = ricep(1:n1,i,j)/dn0
          nice = nicep(1:n1,i,j)
          qsnow = rsnowp(1:n1,i,j)/dn0
          where (qsnow>0)
            nsnow = snow%nr
          elsewhere
            nsnow = 0
          end where 
          qgrp = rgrp(1:n1,i,j)/dn0
          where (qgrp>0)
            ngrp = graupel%nr
          elsewhere
            ngrp = 0
          end where 
          nin_active = nice + nsnow + ngrp

          convice = alvi/cp*(pi0+pi1+exner(1:n1,i,j))/cp
        end if
        if (lrandommicro) call shuffle(microseq)
         do n=1,nprocess
          select case(microseq(n))
          case(iwtrdff)
            call resetvar(cldw,qc)
            call resetvar(rain,qrain,nrain)
            qv = qv*dn0
            qc = qc*dn0
            qrain = qrain*dn0
            call wtr_dff_SB(n1,dn0,qrain,nrain,qc,qsat,qv,tl,temp)
            qv = qv/dn0
            qc = qc/dn0
            qrain = qrain/dn0
          case(iauto)
            call resetvar(cldw,qc)
            call resetvar(rain,qrain,nrain)
            qc = qc*dn0
            qrain = qrain*dn0
            call auto_SB(n1,dn0,qc,qrain,nrain,tl,dissip(1:n1,i,j))
            qc = qc/dn0
            qrain = qrain/dn0
          case(iaccr)
            call resetvar(cldw,qc)
            call resetvar(rain,qrain,nrain)
            qc = qc*dn0
            qrain = qrain*dn0
            call accr_SB(n1,dn0,qc,qrain,nrain,tl,dissip(1:n1,i,j))
            qc = qc/dn0
            qrain = qrain/dn0
          case(isedimrd)
            call resetvar(rain,qrain,nrain)
            qrain = qrain*dn0
            call sedim_rd(n1,dt,dn0,qrain,nrain)
            qrain = qrain/dn0
          case(isedimcd)
            call resetvar(cldw,qc)
            qc = qc*dn0
            call sedim_cd(n1,dt,tl,qc)
            qc = qc/dn0
          case(iicenucnr)
            where (qsup<0.) qsup = 0.
            call n_icenuc(n1,ninuc,nin_active,temp,qv,qsup)
          case(iicenuc)
            where (ninuc<0.) ninuc = 0.
            where (qsup<0.) qsup = 0.
            call ice_nucleation(n1,ninuc,qice,nice,qsup,tl,temp)

          case(ifreez)
            call resetvar(cldw,qc)
            call resetvar(rain,qrain,nrain)

            call cloud_freeze(n1,qc,ninuc,qice,nice,tl,temp)
            call rain_freeze(n1,qrain,nrain,ninuc,qice,nice,qgrp,temp)
          case(idep)
            where (qsup<0.) qsup = 0.
            call deposition(n1,ice,ninuc,qice,nice,qv,tl,temp,qsup)
            where (qsup<0.) qsup = 0.
            call deposition(n1,snow,ninuc,qsnow,nsnow,qv,tl,temp,qsup)
            where (qsup<0.) qsup = 0.
            call deposition(n1,graupel,ninuc,qgrp,ngrp,qv,tl,temp,qsup)
          case(imelt_ice)
            call resetvar(ice,qice,nice)
            call melting(n1,ice,qice,tl,nice,qc,qrain,nrain,temp)
          case(imelt_snow)
            call resetvar(snow,qsnow)
            call melting(n1,snow,qsnow,tl,nsnow,qc,qrain,nrain,temp)
          case(imelt_grp)
            call resetvar(graupel,qgrp)
            call melting(n1,graupel,qgrp,tl,ngrp,qc,qrain,nrain,temp)
          case(iself_ice)
            call resetvar(ice,qice,nice)
            call ice_selfcollection(n1,ice,snow,qice,qsnow,nice,temp,d_conv_is,q_crit_is,d_crit_is)
          case(icoll_ice_snow)
            call resetvar(ice,qice,nice)
            call resetvar(snow,qsnow)
            call ice_collection(n1,ice,snow,qice,nice,qsnow,temp,q_crit_is)
          case(icoll_ice_grp)
            call resetvar(ice,qice,nice)
            call resetvar(graupel,qgrp)
            call ice_collection(n1,ice,graupel,qice,nice,qgrp,temp,q_crit_is)
      
          case(icoll_snow_grp)
            call resetvar(ice,qice,nice)
            call resetvar(snow,qsnow)
            call resetvar(graupel,qgrp)
            call ice_collection(n1,snow,graupel,qsnow,nsnow,qgrp,temp,q_crit_is)
          case(iriming_ice_cloud)
            call resetvar(cldw,qc)
            call resetvar(ice,qice,nice)
            call ice_cloud_riming(n1,cldw,ice,qc,qice,nice,qgrp,tl,temp,d_coll_c,q_crit_c,d_crit_c,q_crit_ic,d_crit_ic,d_conv_ig,e_ic)
          case(iriming_snow_cloud)
            call resetvar(snow,qsnow)
            call resetvar(cldw,qc)
            call ice_cloud_riming(n1,cldw,snow,qc,qsnow,nsnow,qgrp,tl,temp,d_coll_c,q_crit_c,d_crit_c,q_crit_sc,d_crit_sc,d_conv_sg,e_sc)
          case(iriming_grp_cloud)
            call resetvar(graupel,qgrp)
            call resetvar(cldw,qc)
            q1 = 0.
            call ice_cloud_riming(n1,cldw,graupel,qc,qgrp,ngrp,q1,tl,temp,d_coll_c,q_crit_c,d_crit_c,q_crit_gc,d_crit_gc,d_conv_sg,e_gc)
            qgrp = qgrp + q1
          case(iriming_ice_rain)
            call resetvar(ice,qice,nice)
            call resetvar(rain,qrain,nrain)
            call ice_rain_riming(n1,rain,ice    ,qrain,nrain,qice ,nice ,qgrp,temp,q_crit,d_crit_r, q_crit_ir,d_crit_ir)
!             call ice_rain_riming(n1,rain,ice    ,qrain,nrain,qice ,nice ,qgrp,temp,d_crit_r,1e-2,1e-2,1e-2)
          case(iriming_snow_rain)
            call resetvar(snow,qsnow)
            call resetvar(rain,qrain,nrain)
            call ice_rain_riming(n1,rain,snow   ,qrain,nrain,qsnow,nsnow,qgrp,temp,q_crit,d_crit_r, q_crit_sr,d_crit_sr)
          case(iriming_grp_rain)
            call resetvar(graupel,qgrp)
            call resetvar(rain,qrain,nrain)
            q1 = 0.
            call ice_rain_riming(n1,rain,graupel,qrain,nrain,qgrp ,ngrp,q1,temp,q_crit,d_crit_r, q_crit,0.)
            qgrp = qgrp + q1
          case(ised_ice)
            call resetvar(ice,qice,nice)
            call resetvar(snow,qsnow)
            call resetvar(graupel,qgrp)
            call sedimentation (n1,ice, qice,nice)
            call sedimentation (n1,snow, qsnow)
            call sedimentation (n1,graupel, qgrp)
          case default
          end select
        end do
        call resetvar(cldw,qc)
        call resetvar(rain,qrain,nrain)
        if (level==4) then
          call resetvar(ice,qice,nice)
          call resetvar(snow,qsnow)
          call resetvar(graupel,qgrp)
        end if

        tlt(2:n1,i,j) =  tlt(2:n1,i,j)+      (tl(2:n1) - thl(2:n1,i,j))/dt
        rtt(2:n1,i,j) = rtt(2:n1,i,j) +(dn0(2:n1)*qv(2:n1) - rv(2:n1,i,j))/dt + (dn0(2:n1)*qc(2:n1) - rc(2:n1,i,j))/dt
        rpt(2:n1,i,j) = max(rpt(2:n1,i,j) +(dn0(2:n1)*qrain - rp(2:n1,i,j))/dt,-rp(2:n1,i,j)/dt)
        npt(2:n1,i,j) = max(npt(2:n1,i,j) +(nrain(2:n1) - np(2:n1,i,j))/dt,-np(2:n1,i,j)/dt)
        
        if (level == 4) then
          ninuct(2:n1,i,j) = ninuct(2:n1,i,j) + (ninuc(2:n1)-ninucp(2:n1,i,j))/dt !max(ninuct(2:n1,i,j) + (ninuc(2:n1) - ninucp(2:n1,i,j))/dtrk,-ninucp(2:n1,i,j)/dt)
          ricet(2:n1,i,j)  = max(ricet(2:n1,i,j)  + (dn0(2:n1)*qice(2:n1) - ricep(2:n1,i,j))/dt,-ricep(2:n1,i,j)/dt)
          nicet(2:n1,i,j)  = nicet(2:n1,i,j)  + (nice(2:n1) - nicep(2:n1,i,j))/dt
          rsnowt(2:n1,i,j) = max(rsnowt(2:n1,i,j) + (dn0(2:n1)*qsnow(2:n1) - rsnowp(2:n1,i,j))/dt,-rsnowp(2:n1,i,j)/dt)
          rgrt(2:n1,i,j)   = max(rgrt(2:n1,i,j)   + (dn0(2:n1)*qgrp(2:n1) - rgrp(2:n1,i,j))/dt,-rp(2:n1,i,j)/dt)
        end if

     end do
    end do
    deallocate(convice,convliq)
!print *,maxval(rsup),maxloc(rsup)
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
        if (rc(k) > 0.) then
          Xc = rc(k)/(cldw%nr+eps0)
          k_c = kc_0
          nu_c = cldw%nu!1*dn0(k)*rc(k)+cldw%nu2

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
        else
          au = 0.
        end if
          rp(k) = rp(k) + au*dt
          rc(k) = rc(k) - au*dt
          tl(k) = tl(k) + convliq(k)*au*dt
          np(k) = np(k) + au/cldw%x_max*dt
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
    real, parameter :: k_1 = 5.e-4
    real, parameter :: Cac = 67.     ! accretion coefficient in KK param.
    real, parameter :: Eac = 1.15    ! accretion exponent in KK param.

    integer :: k
    real    :: tau, phi, ac, sc, k_r, epsilon

      do k=2,n1-1
          if (rp(k) > 0.) then

           k_r = k_r0

            !
            ! Simulate the effect of turbulence on the collision kernel
            ! (dissipation needs to be converted to cm^2/s^3)
            !
            if (turbulence) then
                epsilon = min(600.,diss(k)*1.e4)   ! put an upper limit to the dissipation rate used
                k_r = k_r*(1+(0.05*epsilon**0.25))
            end if
          if (rc(k) > 0.) then
            tau = 1.0-rc(k)/(rc(k)+rp(k)+eps0)
            tau = MIN(MAX(tau,eps0),1.)
            phi = (tau/(tau+k_1))**4


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
          end if
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
              vn(k) = max(0.,an * Dp + bn)
              vr(k) = max(0.,aq * Dp + bq)
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
  subroutine sedim_cd(n1,dt,tl,rc)
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

    real, parameter :: c = 1.19e8 ! Stokes fall velocity coef [m^-1 s^-1]
    real, parameter :: sgg = 1.2  ! geometric standard dev of cloud droplets

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
             vc = min(c*(Dc*0.5)**2 * exp(4.5*(log(sgg))**2),1./(dzt(k)*dt))
!               vc = min(c*(Dc*0.5)**2 * exp(5*(log(1.3))**2),1./(dzt(k)*dt))
             rfl(k) = - rc(k) * vc
             !
             kp1=k+1
             flxdiv = (rfl(kp1)-rfl(k))*dzt(k)
             rc(k) = rc(k)-flxdiv*dt
             tl(k) = tl(k)+flxdiv*convliq(k)*dt
          end do

  end subroutine sedim_cd
  subroutine n_icenuc(n1,nin,nin_active,tk,qv,qsup)
    integer, intent(in) :: n1
    real, intent(inout), dimension (n1) :: nin
    real, intent(in) , dimension(n1) :: tk,qsup,qv,nin_active
    integer :: k
!     real :: fact
    real :: nineq
!     fact = (1-exp(-dt/timenuc))

    do k=1,n1
!       nineq = n_ice_meyers_contact(tk(k),min(qsup(k),0.25))
!       nineq = 
!       nin(k)  = nin(k) + (n_ice_meyers_contact(tk(k),min(qsup(k),0.25))-nin(k))*fact
!       nin(k) = min(nin(k) + nineq*dt/timenuc,nineq)
!       nin(k)  = n_ice_meyers_contact(tk(k),min(qsup(k),0.25))
      if (qsup(k)/qv(k)>0.05) then
        nin(k) = max(0.,nin_set -max(nin_active(k),0.))
      else
        nin(k) = 0.
      end if
    end do

  end subroutine n_icenuc

  subroutine ice_nucleation(n1,nin,qice,nice,qsup,tl,tk)
    integer,intent(in) :: n1
    real, intent(inout),dimension(n1) :: nin,nice, qice, tl,qsup
    real,intent(in), dimension(n1) :: tk
    real :: nuc_n, nuc_q
    integer :: k
    
    do k=2,n1
      if (tk(k) < tmelt .and. qsup(k) > 0.0) then
        nuc_n = nin(k) 

        if (nuc_n>eps0) then
!          nuc_q = min(nuc_n * ice%x_min, qsup(k))
          nuc_n = nuc_q / ice%x_min                !axel 20040416
          nin(k)  = nin(k)  - nuc_n
          qsup(k) = qsup(k) - nuc_q
          nice(k) = nice(k) + nuc_n
          qice(k) = qice(k) + nuc_q
          tl(k)   = tl(k)   + convice(k)*nuc_q
        end if
      endif
    end do
  end subroutine
  subroutine cloud_freeze (n1,qcloud,ninuc,qice,nice,tl,tk)
!

    integer, intent(in) :: n1
    real, intent(inout), dimension(n1) :: qice,nice,ninuc,qcloud,tl
    real, intent(in), dimension(n1) :: tk
!     real, intent(in) :: nc
!         ! .. local variables ..
        integer          ::  k
        real :: frq, frn,qc,xc,jhet,jhom,jtot,tc
        real, parameter :: ahet = 6.5e-1 ! 1/k,      messung nach barklie and gokhale
        real, parameter :: bhet = 2.0e+2 ! 1/(m3 s), messung nach barklie and gokhale
        real,save       :: facg
          facg = moment_gamma(cldw,2)    ! <hn
        do k = 1, n1
          if (tk(k) < tmelt .and. qcloud(k) > 0.) then


            qc = qcloud(k)
            tc = tk(k) - tmelt
            if (tc < -50.0) then
              frq = qc                                               !..komplettes hom. gefrieren
              frn = cldw%nr                                               !..unterhalb -50 c
            else
              xc = min(max(qc/(cldw%nr+eps0),cldw%x_min),cldw%x_max)    !..mean masse

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
              frn  = min(jtot * qc,ninuc(k))
              frq  = frn * xc * facg
            end if
            frq  = min(frq,qc)
!             if (frq>eps0)  frn  = cldw%nr*frq/qc
            !..berechnung der h2o-komponenten
            qcloud(k) = qcloud(k) - frq
            tl(k) = tl(k)+ (convice(k)-convliq(k))*frq

! 
            frn  = max(frn,frq/cldw%x_max)
            ninuc(k) = ninuc(k)   - frn
            qice(k)   = qice(k)   + frq
            nice(k)   = nice(k)   + frn

            ! ub<<

          end if
        end do

      end subroutine cloud_freeze

      subroutine rain_freeze (n1,qrain,nrain,ninuc,qice,nice,qgrp,tk)
        integer, intent(in) :: n1
        real, dimension(n1), intent(inout) :: qrain,nrain,qice,nice,qgrp,ninuc
        real, dimension(n1), intent(in) :: tk
        ! .. local variables ..
        integer                     :: k
        real            :: fr_q,fr_n,q_r,x_r,n_r,j_het,&
            &  fr_q_i,fr_n_i,fr_q_g,fr_n_g,n_0,lam!,fr_q_tmp,fr_n_tmp
        real, save      :: coeff_z,xmax_ice
        real, parameter :: a_het = 6.5e-1 ! messung nach barklie and gokhale (pk s.350)
        real, parameter :: b_het = 2.0e+2 ! messung nach barklie and gokhale (pk s.350)
        real, parameter :: q_crit_fr = 1.000e-6 ! q-schwellenwert fuer rain_freeze
        real, parameter :: d_rainfrz_ig = 0.50e-3 !  rain --> ice oder graupel
          coeff_z = moment_gamma(rain,2)
  !
          xmax_ice = (d_rainfrz_ig/rain%a_geo)**(1.0e0/rain%b_geo)
        fr_q_g = 0.0
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
                  else
                    fr_q = 0.0
                    fr_n = 0.0
                    fr_n_i= 0.0
                    fr_q_i= 0.0
                    fr_n_g= 0.0
                    fr_q_g= 0.0
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

                  else                           !..heterogenes gefrieren
                    j_het = max(b_het * ( exp( a_het * (tmelt - tk(k))) - 1.0 ),0.) / rowt * dt
                    ! ub>> je nach groesse werden die gefrorenen regentropfen dem wolkeneis zugeschlagen
                    !      oder dem graupel oder hagel. hierzu erfolgt eine partielle integration des spektrums von 0
                    !      bis zu einer ersten trennmasse xmax_ice (--> eis), von dort bis zu xmax_gr (--> graupel)
                    !      und von xmax_gr bis unendlich (--> hagel).

                    if (j_het >= 1-20) then
                      fr_n  = min(j_het * q_r,ninuc(k))
                      fr_q  = fr_n * x_r * coeff_z

                      lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                      n_0 = min(ninuc(k),rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu))
                      fr_n_i = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                          incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                      fr_q_i = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                          incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_ice**rain%mu)
                      fr_n = min(fr_n,n_r)
                      fr_q = min(fr_q,q_r)
                    else
                      fr_n= 0.0
                      fr_q= 0.0
                      fr_n_i= 0.0
                      fr_q_i= 0.0
                      fr_q_g= 0.0
                    end if

                  end if

                  fr_n = fr_n
                  fr_q = fr_q

                  fr_n_i = fr_n_i
                  fr_q_i = fr_q_i
                  fr_q_g = fr_q_g 
                end if
                !..berechnung der h2o-komponenten
                qrain(k) = qrain(k) - fr_q
                nrain(k) = n_r - fr_n

                if (qrain(k) < 0.0) then
                  write (*,*) 'seifert rain_freeze: qrain < 0.0, ', k, tk(k), q_r, j_het, fr_q
                  qrain(k) = 0.0e0
                  stop
                end if
                if (nrain(k) < 0.0) then
                  write (*,*) 'seifert rain_freeze: nrain < 0.0, ', k, tk(k), n_r, j_het, fr_n
                  nrain(k) = 0.0e0
                  stop
                end if
              qice(k) = qice(k)  + fr_q_i
              nice(k) = nice(k)  + fr_n_i
              qgrp(k) = qgrp(k)  + fr_q_g

              end if
        end do

      end subroutine rain_freeze

      subroutine deposition(n1,meteor,ninuc,qice,nice,qv,tl,tk,qsup)
      integer, intent(in) :: n1
      type(particle), intent(in) :: meteor
      real, dimension(n1), intent(inout) :: ninuc,qice,nice,qv,tl,qsup
      
      real, dimension(n1), intent(in)    :: tk

    
      ! locale variablen
      integer                     :: k
      real        :: q_g,n_g,x_g,d_g,v_g,f_v,f_n,n_re,f_v_fakt,vent_fakt
      real        :: c_g                 !..koeff. fuer mean kapazitaet
      real        :: a_f,b_f,a_n,b_n     !..koeff. fuer meann ventilationkoeff.
      real :: dep,ndep,gi
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

          x_g = min(max(q_g/(n_g+eps0),max(meteor%x_min,eps0)),meteor%x_max)  !..mean masse
          ! x_g**b_geo durch exp(b_geo*log(x_g)) ersetzen, da x_g == 0 nicht mehr vorkommen kann; 20 % schneller:
          d_g = meteor%a_geo * exp(meteor%b_geo*log(x_g))          !..mean durchmesser
          v_g = meteor%a_vel * exp(meteor%b_vel*log(x_g)) * dn0(k)  !..mean sedimentationsgeschw.

          n_re = v_g * d_g / nu_l                                    !..mean reynoldszahl
          ! n_re**m_f durch exp(m_f*log(n_re)) ersetzen, da n_re == 0 nicht mehr vorkommen kann:
          f_v  = a_f + b_f * f_v_fakt * exp(m_f*log(n_re))           !..mean vent.koeff.
          ! ub>> schnellere berechnung wg. verzicht auf "**":
          !              f_n  = a_n + b_n * f_v_fakt * n_re**m_f                    !..mean vent.koeff.
          f_n = a_n + vent_fakt * (f_v - a_f)                        !..mean vent.koeff.
          ! ub<<
          f_v  = max(f_v,1.e0) !unnoetig??
          f_n  = max(f_n,1.e0) !unnoetig??
!           e_si = e_es(tk(k))
        
          gi = 4.0*pi / ( alvi**2 / (K_T * Rm * tk(k)**2) + Rm * tk(k) / (D_v * e_es(tk(k))) )
          !si = 1-(qv(k)+qsup(k))/rsif(p,tk(k))
	  !ndep  = min(gi * n_g * c_g * d_g * qsup(k)/qv(k) * dt * f_n / x_g,ninuc(k))
          dep = gi * n_g * c_g * d_g * f_v * qsup(k)/qv(k) * dt 
          qice(k) = qice(k) + dep
          qsup(k) = qsup(k) - dep
          qv(k) = qv(k) - dep
          tl(k) = tl(k) + convice(k)*dep
         ! if (meteor%moments==2) then
         !   nice(k) = nice(k) + ndep 
	 !   ninuc(k) = ninuc(k) - ndep
         ! end if
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
    real            :: a_melt_n,b_melt_n
    real            :: a_melt_q,b_melt_q



      a_melt_n = vent_coeff_a(meteor,0)
      b_melt_n = vent_coeff_b(meteor,0)
      a_melt_q = vent_coeff_a(meteor,1)
      b_melt_q = vent_coeff_b(meteor,1)

    do k = 2,n1
          t_a = tk(k) !wrf!+ t(k) + t_g(k)
          e_a = e_ws(T_a)                                     !..Saettigungsdampfdruck
                   
          if (t_a > tmelt .and. q(k) > 0.0) then

            x_m = min(max(q(k)/(nr(k)+eps0),meteor%x_min),meteor%x_max)  !..mean qe in si

            D_m = meteor%a_geo * x_m**meteor%b_geo                   !..mean Durchmesser
            v_m = meteor%a_vel * x_m**meteor%b_vel * rho_0  !..mean Sedimentationsgeschw.

            N_re = v_m * D_m / nu_l                             !..mean Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mean Vent.Koeff. Wasserdampf
  !          fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mean Vent.Koeff. Wasserdampf

            D_T  = Kt / (cp * rho_0)!WRF!+rho_g(i,j,k)))
            fh_q = D_T / D_v * fv_q
 !           fh_n = D_T / D_v * fv_n

            melt   = 2.0*pi / (alvi-alvl) * D_m * nr(k) * dt

            melt_h = melt * Kt * (T_a - tmelt)
            melt_v = melt * D_v*alvl/Rm * (e_a/T_a - e_3/tmelt)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_m

! ub>> setzte melt_n so, dass x_h beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q(k)) / x_m + nr(k), 0.0e0), nr(k))

            melt_q = MIN(q(k),melt_q)
            melt_n = MIN(nr(k),melt_n)

            melt_q = MAX(0.e0,melt_q)
            melt_n = MAX(0.e0,melt_n)
            q(k) = q(k) - melt_q
            if (meteor%nr==2) then
              nr(k) = nr(k) - melt_n
            end if
            if(x_m<cldw%x_max) then
              qcld(k)    = qcld(k)    + melt_q
              thl(k)      = thl(k)      - convice(k)*melt_q
! 
            else
              qrain(k)    = qrain(k)    + melt_q
              nrain(k)    = nrain(k)    + melt_n
            end if
          endif
        enddo



  end subroutine melting

  subroutine ice_selfcollection(n1,metin,metout,q_in,q_out,n_in,tk,d_conv,q_crit,d_crit)
    integer, intent(in) :: n1
    real, dimension(n1), intent(in) :: tk
    type(particle), intent(in) :: metin, metout
    real, dimension(n1), intent(inout) :: q_in,q_out,n_in
    real, intent(in) :: d_conv, q_crit, d_crit

  
!     !*******************************************************************************
!     !                                                                              *
!     !       berechnung des selbsteinfangs der eispartikel                          *
!     !                                                                              *
!     !*******************************************************************************
! 
!     use globale_variablen,  only: loc_ix, loc_iy, loc_iz, t, p, q, p_g, t_g, rho_g, &
!          &                        q_ice, n_ice, q_snow, n_snow, dt, dqdt, speichere_umwandlungsraten
!     use konstanten,         only: pi,cv,cp,wolke_typ
!     use parallele_umgebung, only: isio
!     use initialisierung,    only: t_0,p_0,rho_0
! 
!     implicit none
! 
!     ! locale variablen
    integer                     :: k
    real       :: q_i,n_i,x_i,d_i,v_i,e_coll
    real, save :: x_conv
    real       :: self_n,self_q
      real, dimension(3),save :: delta_n,delta_q, theta_n,theta_q
      logical, save :: firsttime(3) = .true.
      integer :: metnr = 0
      select case (metin%name)
      case('ice')
        metnr = 1
      case('snow')
        metnr = 2
      case('graupel')
        metnr = 3
       case default
        stop
     end select
      if (firsttime(metnr)) then
        firsttime(metnr) = .false.
        delta_n(metnr) = 2.0*coll_delta_11(metin,0) + coll_delta_12(metin,metin,0)
        delta_q(metnr) = 2.0*coll_delta_11(metin,1) + coll_delta_12(metin,metin,1)
        theta_n(metnr) = 2.0*coll_theta_11(metin,0) + coll_theta_12(metin,metin,0)
        theta_q(metnr) = 2.0*coll_theta_11(metin,1) + coll_theta_12(metin,metin,1)

        x_conv = (d_conv/metout%a_geo)**(1./metout%b_geo)

      end if
! 
!       delta_n(metnr) = 2.0*coll_delta_11(metin,0) + coll_delta_12(metin,metin,0)
!       delta_q(metnr) = 2.0*coll_delta_11(metin,1) + coll_delta_12(metin,metin,1)
!       theta_n(metnr) = 2.0*coll_theta_11(metin,0) + coll_theta_12(metin,metin,0)
!       theta_q(metnr) = 2.0*coll_theta_11(metin,1) + coll_theta_12(metin,metin,1)
! 

    do k = 1, n1
          q_i = q_in(k)                                   !..fluessigwassergehalt in si
          n_i = n_in(k)                                   !..anzahldichte in si
          x_i = min(max(q_i/(n_i+eps0),metin%x_min),metin%x_max)    !..mean masse in si
          d_i = metin%a_geo * x_i**metin%b_geo                     !..mean durchmesser
          if ( n_i > 0.0 .and. q_i > q_crit .and. d_i > d_crit ) then
            if (tk(k)>tmelt) then
              e_coll = 1.0
            else
              !.. temperaturabhaengige efficiency nach cotton et al. (1986)
              !   (siehe auch straka, 1989; s. 53)
              e_coll = min(10**(0.035*(tk(k)-tmelt)-0.7),0.2e0)
              !.. temperaturabhaengige efficiency nach lin et al. (1983)
!               e_coll = 0.001*min(exp(0.09*(tk(k)-tmelt)),1.0e0)
              !e_coll = max(e_ii,min(exp(0.09*(t_a-t_3)),1.0e0))
            end if


            v_i = metin%a_vel * x_i**metin%b_vel * dn0(k)      !..mean sedimentationsgeschw.


            self_q = pi * 0.25e0 * e_coll * delta_q(metnr) * n_i * q_i * d_i * d_i &
                  * ( theta_q(metnr) * v_i * v_i + 2.0*metin%s_vel**2 )**0.5 * dt

            self_q = min(self_q,q_i)

            q_in(k)  = q_in(k)  - self_q
            q_out(k) = q_out(k) + self_q

              self_n = pi * 0.25e0 * e_coll * delta_n(metnr) * n_i * n_i * d_i * d_i &
                      * ( theta_n(metnr) * v_i * v_i + 2.0*metin%s_vel**2 )**0.5 * dt
              self_n = min(min(self_n,self_q/x_conv),n_i)
              n_in(k)  = n_in(k)  - self_n


          endif
        enddo

  end subroutine ice_selfcollection
  subroutine ice_cloud_riming(n1,cloud,ice,q_c,q_i,n_i,q_g,thl,tk,d_coll,q_crit_c,d_crit_c,q_crit_i,d_crit_i,d_conv,e_ic)
      integer, intent(in) :: n1
      type(particle), intent(in) :: cloud,ice
      real, dimension(n1), intent(inout) :: q_c,q_i,n_i,q_g,thl
      real, dimension(n1), intent(in) :: tk
      real, intent(in) :: d_coll,q_crit_c,d_crit_c,q_crit_i,d_crit_i,e_ic
      real, intent(in), optional :: d_conv
      
      integer                     :: k
  real, parameter :: e_min = 0.01              !..min. eff. fuer gc,ic,sc
!   real, parameter :: x_conv    = 0.100e-9 ! minimale graupel-/hagelmasse riming

      real     :: x_i,d_i,v_i
      real     :: x_c,d_c,v_c,e_coll,x_coll_c
      real     :: rime_q
      real     :: conv_n,conv_q
      real     :: mult_n,mult_q,mult_1,mult_2
!       real     :: delta_n_ii,delta_n_ic,delta_n_cc
      real     :: delta_q_ii,delta_q_ic,delta_q_cc
!       real     :: theta_n_ii,theta_n_ic,theta_n_cc
      real     :: theta_q_ii,theta_q_ic,theta_q_cc
      real     :: const1,const2,const3,const4,const5
      real, dimension(2,2,3),save :: delta, theta
      logical, save :: firsttime(3) = .true.
      integer :: metnr = 0
      select case (ice%name)
      case('ice')
        metnr = 1
      case('snow')
        metnr = 2
      case('graupel')
        metnr = 3
       case default
        stop
     end select
      if (firsttime(metnr)) then
        firsttime(metnr) = .false.
!         delta_n_ii = coll_delta_11(ice,0)
!         delta_n_ic = coll_delta_12(ice,cloud,0)
!         delta_n_cc = coll_delta_22(cloud,0)
        delta(1,1,metnr) = coll_delta_11(ice,1)
        delta(1,2,metnr) = coll_delta_12(ice,cloud,1)
        delta(2,2,metnr) = coll_delta_22(cloud,1)

!         theta_n_ii = coll_theta_11(ice,0)
!         theta_n_ic = coll_theta_12(ice,cloud,0)
!         theta_n_cc = coll_theta_22(cloud,0)
        theta(1,1,metnr) = coll_theta_11(ice,1)
        theta(1,2,metnr) = coll_theta_12(ice,cloud,1)
        theta(2,2,metnr) = coll_theta_22(cloud,1)
      end if
      delta_q_ii =  delta(1,1,metnr)
      delta_q_ic =  delta(1,2,metnr)
      delta_q_cc =  delta(2,2,metnr)
      theta_q_ii =  theta(1,1,metnr)
      theta_q_ic =  theta(1,2,metnr)
      theta_q_cc =  theta(2,2,metnr)

      x_coll_c = (d_coll/cldw%a_geo)**3          !..mindestmasse fuer collection, begrenzt rime_n

      const1 = e_ic/(d_coll - d_crit_c)
      const2 = 1/x_coll_c
      const3 = 1/(t_mult_opt - t_mult_min)
      const4 = 1/(t_mult_opt - t_mult_max)
      const5 = alpha_spacefilling * rowt/roice
      do k = 1, n1

            x_c = min(max(q_c(k)/(cloud%nr+eps0),cloud%x_min),cloud%x_max)  !..mean masse
            d_c = cloud%a_geo * x_c**cloud%b_geo                   !..mean durchmesser

            x_i = min(max(q_i(k)/(n_i(k)+eps0),ice%x_min),ice%x_max)      !..mean masse
            d_i = ice%a_geo * x_i**ice%b_geo                       !..mean durchmesser


            if (q_c(k) > q_crit_c .and. q_i(k) > q_crit_i .and. d_i > d_crit_i .and. d_c > d_crit_c) then

              v_c = cloud%a_vel * x_c**cloud%b_vel * dn0(k)  !..mean sedimentationsgeschw.
              v_i = ice%a_vel   * x_i**ice%b_vel   * dn0(k) !..mean sedimentationsgeschw.

              e_coll = min(e_ic, max(const1*(d_c - d_crit_c), e_min))

              rime_q = pi/4.0 * e_coll * n_i(k) * q_c(k) * dt &
                  &   *     (delta_q_ii * d_i*d_i + delta_q_ic * d_i*d_c + delta_q_cc * d_c*d_c) &
                  &   * sqrt(theta_q_ii * v_i*v_i - theta_q_ic * v_i*v_c + theta_q_cc * v_c*v_c  &
                  &          +ice%s_vel**2)
                rime_q = min(q_c(k),rime_q)

                q_i(k)   = q_i(k)  + rime_q
                q_c(k) = q_c(k) - rime_q
                thl(k) = thl(k) + convice(k)*rime_q
!                 cloud%nrloud(k) = cloud%nrloud(k) - rime_n
                ! ub>>
                ! ub<<

                ! eismultiplikation nach hallet und mossop

!                 mult_q = 0.0
                if (tk(k) < tmelt .and. ice_multiplication .and. ice%moments == 2) then
                  mult_1 = (tk(k) - t_mult_min)*const3
                  mult_2 = (tk(k) - t_mult_max)*const4
                  mult_1 = max(0.e0,min(mult_1,1.e0))
                  mult_2 = max(0.e0,min(mult_2,1.e0))
                  mult_n = c_mult * mult_1 * mult_2 * rime_q

                  n_i(k)  = n_i(k)  + mult_n
                endif

                ! umwandlung ice -> graupel
                  conv_q = 0.0
                  conv_n = 0.0
                  
                if(present(d_conv)) then
                  if (d_i > d_conv) then
                    conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*roice*d_i**3/x_i - 1.0) )
                    ! d_i darf durch conv nicht kleiner werden als d_conv_ig
                    conv_q = min(q_i(k)-n_i(k)*(d_conv_ig/ice%a_geo)**(1.0/ice%b_geo),conv_q)
                    conv_q = min(q_i(k),conv_q)
                    ! ub >>
                    x_i = min(max((q_i(k))/(n_i(k)+eps0),ice%x_min),ice%x_max)    !..mean masse incl. riming
                    ! ub <<
                    conv_n = conv_q / x_i
                    conv_n = min(n_i(k),conv_n)
                  end if
                endif
              
                q_i(k)     = q_i(k)     - conv_q
                q_g(k) = q_g(k) + conv_q

                if (ice%moments == 2) n_i(k)     = n_i(k)     - conv_n
!                 n_graupel(k) = n_graupel(k) + conv_n
                ! ub>>

!               else
! 
!                 rimeqcrate_ice(k) = rimeqcrate_ice(k) + rime_q
!                 rimencrate_ice(k) = rimencrate_ice(k) + rime_n
! 
!               end if
              ! ub<<

              ! temperatur und druck

              !wrf! if (tk(k) < tmelt) then
              !wrf!   t(k)     = t(i, j, k)   + l_ew / cp * rime_q / (rho_0(k)+rho_g(k))
              !wrf!   p(k)     = p(i, j, k)   + l_ew / cv * rime_q * r_l
              !wrf! endif
            endif
          enddo

    end subroutine ice_cloud_riming
  subroutine ice_rain_riming(n1,rain,ice,q_r,n_r,q_i,n_i,q_g,tk,q_crit_r, d_crit_r,q_crit_i,d_crit_i)
      integer, intent(in) :: n1
      type(particle), intent(in) :: rain,ice
      real, dimension(n1), intent(inout) :: q_r,n_r,q_i,n_i,q_g
      real, dimension(n1), intent(in) :: tk
      real, intent(in) :: q_crit_r,d_crit_r,q_crit_i,d_crit_i

      integer                     :: k

    real            :: x_i,d_i,v_i,d_id
    real            :: x_r,d_r,v_r,d_rd
    real            :: rime_n,rime_qi,rime_qr
    real            :: mult_n,mult_q,mult_1,mult_2
    real      :: delta_n_ii,delta_n_ir,           delta_n_rr
    real      :: delta_q_ii,delta_q_ir,delta_q_ri,delta_q_rr
    real      :: theta_n_ii,theta_n_ir,           theta_n_rr
    real      :: theta_q_ii,theta_q_ir,theta_q_ri,theta_q_rr
    real, save :: d_av_fakt_r
    real,dimension(3),save :: d_av_fakt_i

      real, dimension(2,2,3),save :: delta_n,delta_q, theta_n,theta_q
      logical, save :: firsttime(3) = .true.
      integer :: metnr = 0
      select case (ice%name)
      case('ice')
        metnr = 1
      case('snow')
        metnr = 2
      case('graupel')
        metnr = 3
       case default
        stop
     end select
      if (firsttime(metnr)) then
        firsttime(metnr) = .false.
!         delta_n_ii = coll_delta_11(ice,0)
!         delta_n_ic = coll_delta_12(ice,cloud,0)
!         delta_n_cc = coll_delta_22(cloud,0)
        delta_n(1,1,metnr) = coll_delta_11(ice,0)
        delta_n(1,2,metnr) = coll_delta_12(ice,rain,0)
        delta_n(2,2,metnr) = coll_delta_22(rain,0)
        delta_q(1,1,metnr) = coll_delta_11(ice,1)
        delta_q(1,2,metnr) = coll_delta_12(ice,rain,1)
        delta_q(2,2,metnr) = coll_delta_22(rain,1)

!         theta_n_ii = coll_theta_11(ice,0)
!         theta_n_ic = coll_theta_12(ice,rain,0)
!         theta_n_cc = coll_theta_22(rain,0)
        theta_n(1,1,metnr) = coll_theta_11(ice,0)
        theta_n(1,2,metnr) = coll_theta_12(ice,rain,0)
        theta_n(2,2,metnr) = coll_theta_22(rain,0)
        theta_q(1,1,metnr) = coll_theta_11(ice,1)
        theta_q(1,2,metnr) = coll_theta_12(ice,rain,1)
        theta_q(2,2,metnr) = coll_theta_22(rain,1)

        d_av_fakt_r = d_average_factor(rain)
        d_av_fakt_i(metnr) = d_average_factor(ice)
      end if
      delta_n_ii =  delta_n(1,1,metnr)
      delta_n_ir =  delta_n(1,2,metnr)
      delta_n_rr =  delta_n(2,2,metnr)
      theta_n_ii =  theta_n(1,1,metnr)
      theta_n_ir =  theta_n(1,2,metnr)
      theta_n_rr =  theta_n(2,2,metnr)
      delta_q_ii =  delta_q(1,1,metnr)
      delta_q_ir =  delta_q(1,2,metnr)
      delta_q_rr =  delta_q(2,2,metnr)
      theta_q_ii =  theta_q(1,1,metnr)
      theta_q_ir =  theta_q(1,2,metnr)
      theta_q_rr =  theta_q(2,2,metnr)


    do k = 1, n1

            x_r = min(max(q_r(k)/(n_r(k)+eps0),rain%x_min),rain%x_max)  !..mean masse
            d_r = rain%a_geo * x_r**rain%b_geo                   !..mean durchmesser

          x_i = min(max(q_i(k)/(n_i(k)+eps0),ice%x_min),ice%x_max)    !..mean masse in si
          d_i = ice%a_geo * x_i**ice%b_geo                     !..mean durchmesser
          d_id = d_av_fakt_i(metnr) * d_i
          
          if (q_r(k) > q_crit_r  .and. d_r > d_crit_r .and. q_i(k) > q_crit_i .and. d_i > d_crit_i) then

            x_r = min(max(q_r(k)/(n_r(k)+eps0),rain%x_min),rain%x_max)  !..mean masse in si

            v_i = ice%a_vel * x_i**ice%b_vel * dn0(k)    !..mean sedimentationsgeschw.

            d_r = rain%a_geo * x_r**rain%b_geo                   !..mean durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * dn0(k)  !..mean sedimentationsgeschw.
            d_rd = d_av_fakt_r * d_r

            rime_n  = pi/4.0 * n_i(k) * n_r(k) * dt &
                 &   * (delta_n_ii * d_i*d_i + delta_n_ir * d_i*d_r + delta_n_rr * d_r*d_r) &
                 &   * sqrt(theta_n_ii * v_i*v_i - theta_n_ir * v_i*v_r + theta_n_rr * v_r*v_r  &
                 &     +ice%s_vel**2)

            rime_qr = pi/4.0 * n_i(k) * q_r(k) * dt &
                 &   * (delta_n_ii * d_i*d_i + delta_q_ir * d_i*d_r + delta_q_rr * d_r*d_r) &
                 &   * sqrt(theta_n_ii * v_i*v_i - theta_q_ir * v_i*v_r + theta_q_rr * v_r*v_r  &
                 &     +ice%s_vel**2)

            rime_qi = pi/4.0 * n_r(k) * q_i(k) * dt &
                 &   * (delta_q_ii * d_i*d_i + delta_q_ri * d_i*d_r + delta_n_rr * d_r*d_r) &
                 &   * sqrt(theta_q_ii * v_i*v_i - theta_q_ri * v_i*v_r + theta_n_rr * v_r*v_r  &
                 &     +ice%s_vel**2)


              rime_n  = min(min(n_r(k),n_i(k)),rime_n)
              rime_qr = min(q_r(k),rime_qr)
              rime_qi = min(q_i(k),rime_qi)

              if (ice%moments==2) n_i(k)  = n_i(k)  - rime_n
              n_r(k) = n_r(k) - rime_n
              q_i(k)  = q_i(k)  - rime_qi
              q_r(k) = q_r(k) - rime_qr


              ! eismultiplikation nach hallet und mossop

              mult_q = 0.0
              if (tk(k) < tmelt .and. ice_multiplication .and. ice%moments==2) then
                mult_1 = (tk(k) - t_mult_min) / (t_mult_opt - t_mult_min)
                mult_2 = (tk(k) - t_mult_max) / (t_mult_opt - t_mult_max)
                mult_1 = max(0.e0,min(mult_1,1.e0))
                mult_2 = max(0.e0,min(mult_2,1.e0))
                mult_n = c_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = min(rime_qr,mult_q)

                n_i(k) = n_i(k)  + mult_n
                q_i(k) = q_i(k)  + mult_q
              endif
              if (tk(k) >= tmelt) then
                if (d_id > d_rd) then
                   ! regen wird wieder abgeshedded
                  if (ice%moments==2) n_i(k) = n_i(k) + rime_n
                  n_r(k) = n_r(k) + rime_qr / x_r
                  q_i(k) = q_i(k) + rime_qi - mult_q
                  q_r(k) = q_r(k) + rime_qr
                else
                   ! eisteilchen schmelzen instantan zu regen.
                  n_r(k) = n_r(k) + rime_n
                  q_r(k) = q_r(k) + rime_qr + rime_qi - mult_q
                end if
              else
                  ! eis + angefrorenes regenwasser ergibt graupel:
!                   n_graupel(k) = n_graupel(k) + rime_n
                  q_g(k) = q_g(k) + rime_qi + rime_qr - mult_q
              end if
              ! ub<<

            !wrf!if (tk(k) < tmelt) then
            !wrf!   t(k)         = t(i, j, k)   + l_ew / cp * rime_qr / (rho_0(k)+rho_g(k))
            !wrf!   p(k)         = p(i, j, k)   + l_ew / cv * rime_qr * r_l
            !wrf!endif
          endif
    enddo

  end subroutine ice_rain_riming
  subroutine ice_collection(n1,met1,met2,q_i,n_i,q_g,tk,q_crit)
      integer, intent(in) :: n1
      type(particle), intent(in) :: met1,met2
      real, dimension(n1), intent(inout) :: q_i,q_g,n_i
      real, dimension(n1), intent(in) :: tk
      real, intent(in) :: q_crit

    integer                     :: k

    real            :: x_g,d_g,v_g
    real            :: x_i,d_i,v_i
    real      :: delta_n_ss,delta_n_si,delta_n_ii
    real      :: delta_q_ss,delta_q_si,delta_q_ii
    real      :: theta_n_ss,theta_n_si,theta_n_ii
    real      :: theta_q_ss,theta_q_si,theta_q_ii

    real            :: e_coll,coll_n,coll_q
      real, dimension(2,2,3,3),save :: delta_n,delta_q, theta_n,theta_q
      logical, save :: firsttime(3,3) = .true.
      integer :: metnr1 = 0,metnr2=0
      select case (met1%name)
      case('ice')
        metnr1 = 1
      case('snow')
        metnr1 = 2
      case('graupel')
        metnr1 = 3
       case default
        stop
     end select
      select case (met2%name)
      case('ice')
        metnr2 = 1
      case('snow')
        metnr2 = 2
      case('graupel')
        metnr2 = 3
      case default
        stop
      end select
      if (firsttime(metnr1,metnr2)) then
        firsttime(metnr1,metnr2) = .false.
!         delta_n_ii = coll_delta_11(met2,0)
!         delta_n_ic = coll_delta_12(met2,cloud,0)
!         delta_n_cc = coll_delta_22(cloud,0)
        delta_n(1,1,metnr1,metnr2) = coll_delta_11(met2,0)
        delta_n(1,2,metnr1,metnr2) = coll_delta_12(met2,met1,0)
        delta_n(2,2,metnr1,metnr2) = coll_delta_22(met1,0)
        delta_q(1,1,metnr1,metnr2) = coll_delta_11(met2,1)
        delta_q(1,2,metnr1,metnr2) = coll_delta_12(met2,met1,1)
        delta_q(2,2,metnr1,metnr2) = coll_delta_22(met1,1)

!         theta_n_ii = coll_theta_11(met2,0)
!         theta_n_ic = coll_theta_12(met2,met1,0)
!         theta_n_cc = coll_theta_22(met1,0)
        theta_n(1,1,metnr1,metnr2) = coll_theta_11(met2,0)
        theta_n(1,2,metnr1,metnr2) = coll_theta_12(met2,met1,0)
        theta_n(2,2,metnr1,metnr2) = coll_theta_22(met1,0)
        theta_q(1,1,metnr1,metnr2) = coll_theta_11(met2,1)
        theta_q(1,2,metnr1,metnr2) = coll_theta_12(met2,met1,1)
        theta_q(2,2,metnr1,metnr2) = coll_theta_22(met1,1)

      end if
      delta_n_ss =  delta_n(1,1,metnr1,metnr2)
      delta_n_si =  delta_n(1,2,metnr1,metnr2)
      delta_n_ii =  delta_n(2,2,metnr1,metnr2)
      theta_n_ss =  theta_n(1,1,metnr1,metnr2)
      theta_n_si =  theta_n(1,2,metnr1,metnr2)
      theta_n_ii =  theta_n(2,2,metnr1,metnr2)
      delta_q_ss =  delta_q(1,1,metnr1,metnr2)
      delta_q_si =  delta_q(1,2,metnr1,metnr2)
      delta_q_ii =  delta_q(2,2,metnr1,metnr2)
      theta_q_ss =  theta_q(1,1,metnr1,metnr2)
      theta_q_si =  theta_q(1,2,metnr1,metnr2)
      theta_q_ii =  theta_q(2,2,metnr1,metnr2)

!       delta_n_ss = coll_delta_11(met2,0)
!       delta_n_si = coll_delta_12(met2,met1,0)
!       delta_n_ii = coll_delta_22(met1,0)
!       delta_q_ss = coll_delta_11(met2,0)
!       delta_q_si = coll_delta_12(met2,met1,1)
!       delta_q_ii = coll_delta_22(met1,1)
! 
!       theta_n_ss = coll_theta_11(met2,0)
!       theta_n_si = coll_theta_12(met2,met1,0)
!       theta_n_ii = coll_theta_22(met1,0)
!       theta_q_ss = coll_theta_11(met2,0)
!       theta_q_si = coll_theta_12(met2,met1,1)
!       theta_q_ii = coll_theta_22(met1,1)


    do k = 1, n1

          if (q_i(k) > q_crit .and. q_g(k) > q_crit) then

            if (tk(k) > tmelt) then
              e_coll = 1.0
            else
              !.. temperaturabhaengige sticking efficiency nach lin (1983)
              e_coll = max(0.1e0,min(exp(0.09*(tk(k)-tmelt)),1.0e0))
              !.. temperaturabhaengige efficiency nach cotton et al. (1986)
              !   (siehe auch straka, 1989; s. 53)
              !e_coll = min(10**(0.035*(tk(k)-tmelt)-0.7),0.2e0)
            endif

            x_i = min(max(q_i(k)/(n_i(k)+eps0),met1%x_min),met1%x_max)     !..mean masse in si
            x_g = min(max(q_g(k)/(met2%nr+eps0),met2%x_min),met2%x_max)   !..mean masse in si

            d_i = met1%a_geo * x_i**met1%b_geo                      !..mean durchmesser
            v_i = met1%a_vel * x_i**met1%b_vel * dn0(k)     !..mean sedimentationsgeschw.

            d_g = met2%a_geo * x_g**met2%b_geo                    !..mean durchmesser
            v_g = met2%a_vel * x_g**met2%b_vel * dn0(k)   !..mean sedimentationsgeschw.


            coll_q = pi/4.0 * met2%nr * q_i(k) * e_coll * dt &
                 &   * (delta_q_ss * d_g**2 + delta_q_si * d_g*d_i + delta_q_ii * d_i**2) &
                 &   * (theta_q_ss * v_g**2 - theta_q_si * v_g*v_i + theta_q_ii * v_i**2  &
                 &     +met2%s_vel**2 + met1%s_vel**2)**0.5

            coll_q = min(q_i(k),coll_q)

            q_g(k) = q_g(k) + coll_q
            q_i(k)  = q_i(k)  - coll_q
            if (met1%nr ==2) then
              coll_n = pi/4.0 * met2%nr * n_i(k) * e_coll * dt &
                  &   * (delta_n_ss * d_g**2 + delta_n_si * d_g*d_i + delta_n_ii * d_i**2) &
                  &   * (theta_n_ss * v_g**2 - theta_n_si * v_g*v_i + theta_n_ii * v_i**2  &
                  &     +met2%s_vel**2 + met1%s_vel**2)**0.5
              coll_n = min(n_i(k),coll_n)
              n_i(k)  = n_i(k)  - coll_n
            end if


          endif
        enddo

  end subroutine ice_collection
  
  subroutine sedimentation(n1,meteor,rp,nr)
    integer, intent(in) :: n1
    real, dimension(n1), intent(inout) :: rp
    type(particle),  intent(in) :: meteor
    real, dimension(n1), intent(inout),optional :: nr
! 
!     ! .. local variables ..
    integer  :: k,kk,km1,kp1
    real,dimension(3),save      :: alfn,alfq,c_lam
    real :: lam,sk,tot,zz,xp
    real, dimension(n1)   :: nfl,rfl,vn,vr,dn,dr,rslope,cn,cr,nslope,np
    real :: cc, flxdiv,maxi,mini
      logical, save :: firsttime(3) = .true.
      integer :: metnr = 0
      select case (ice%name)
      case('ice')
        metnr = 1
      case('snow')
        metnr = 2
      case('graupel')
        metnr = 3
      case default
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
      firsttime(metnr) =.false.
    end if

    do k=n1-1,2,-1

      Xp = rp(k) / (np(k)+eps0)
      lam = ( c_lam(metnr) * xp )**(meteor%b_vel)      
      vr(k) = alfq(metnr) * lam
      vr(k) = max(vr(k),0.1e+0)
      vr(k) = min(vr(k),30.e0)
      vr(k) = -vr(k)

      vn(k) = alfn(metnr) * lam
      vn(k) = max(vn(k),0.1e+0)
      vn(k) = min(vn(k),30.e0)
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
    real, parameter :: n_0 = 1.0e+3
    real, parameter :: n_m = 1.0e+3
    real, parameter :: a_d = -0.639
    real, parameter :: b_d = 12.960
    real, parameter :: c_d = -2.8
    real, parameter :: d_d = 0.262
    real, parameter :: tmelt  = 2.732e+2     !..tripelpunkt wasser
    real :: temp
    temp = t_a - tmelt
    n_ice_meyers_contact = n_0 * exp( a_d + b_d * s )
    n_ice_meyers_contact = n_ice_meyers_contact  &
          &               + n_m * exp( c_d + d_d * temp )

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

    data cof,stp/76.18009173e0,-86.50532033e0,24.01409822e0,  &
         &     -1.231739516e0,.120858003e-2,-.536382e-5,2.50662827465e0/
    data half,one,fpf/0.5e0,1.0e0,5.5e0/

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

  real function incgfct_upper(a,x)
    !*******************************************************************************
  !
  ! upper incomplete gamma function
  !
  ! actual obere unvollstaendige gamma-funktion, hier direkt
  ! das integral
  !              int(x)(oo) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

    implicit none

    real, intent(in) :: a, x
    real :: gam, gln

    gam = gammq(a,x,gln)
    incgfct_upper = exp(gln) * gam

  end function incgfct_upper


  real function incgfct_lower(a,x)
  !*******************************************************************************
  !
  ! lower incomplete gamma function
  !
  ! actual untere unvollstaendige gamma-funktion, hier direkt
  ! das integral
  !              int(0)(x) exp(-t) t^(a-1) dt
  !
!   !*******************************************************************************

    implicit none

    real, intent(in) :: a, x
    real :: gam, gln

    gam = gammp(a,x,gln)
    incgfct_lower = exp(gln) * gam

  end function incgfct_lower


  real function incgfct(a,x1,x2)
  !*******************************************************************************
  !
  ! actual unvollstaendige gamma-funktion, hier direkt
  ! das integral
  !              int(x1)(x2) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

    implicit none

    real, intent(in) :: a, x1, x2

    incgfct = incgfct_lower(a,x2) - incgfct_lower(a,x1)

  end function incgfct

  real function gammp(a,x,gln)

    implicit none

    real, intent(in) :: a, x
    real, intent(out) :: gln

    real :: gammcf, gamser

    if (x.lt.0.0e0 .or. a .le. 0.0e0) then
      write(*,*) 'error in gammp: bad arguments'
      write(*,*) '  (module gamma_functions, src_seifert.f90)'
      gammp = 0.0e0
      return
    end if

    if (x .lt. a+1.)then
      call gser(gamser,a,x,gln)
      gammp = gamser
    else
      call gcf(gammcf,a,x,gln)
      gammp = 1.0e0 - gammcf
    endif
    return
  end function gammp

  real function gammq(a,x,gln)

    implicit none

    real, intent(in) :: a, x
    real, intent(out) :: gln
    real :: gammcf, gamser

    if (x.lt.0.0e0 .or. a .le. 0.0e0) then
      write(*,*) 'error in gammq: bad arguments (module gamma_functions, src_seifert.f90)'
      gammq = 0.0e0
      return
    end if

    if (x.lt.a+1.) then
      call gser(gamser,a,x,gln)
      gammq = 1.0e0 - gamser
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
    data cof /76.18009172947146e0,-86.50532032941677e0, &
         24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2, &
         -.5395239384953e-5/
    data stp /2.5066282746310005e0/

    xx  = x
    tmp = xx + 5.5e0
    tmp = (xx + 0.5e0) * log(tmp) - tmp
    ser = 1.000000000190015e0
    do j = 1,6
       xx  = xx  + 1.0e0
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

    data cof,stp/76.18009173e0,-86.50532033e0,24.01409822e0,  &
          &     -1.231739516e0,.120858003e-2,-.536382e-5,2.50662827465e0/
    data half,one,fpf/0.5e0,1.0e0,5.5e0/

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
      write (*,*) 'error in gcf: a too large, itmax too small (module gamma_functions, src_seifert.f90)'
      gammcf = 0.0e0
      return
    end if

    gammcf=exp(-x+a*log(x)-gln)*h

    return
  end subroutine gcf

  subroutine gser(gamser,a,x,gln)

    implicit none

    integer, parameter :: itmax = 100
    real, parameter :: eps=3.e-7
    real, intent(in) :: a, x
    real, intent(out) :: gamser, gln

    integer :: n
    real :: ap,del,sum

    gln=gammln(a)
    if (x.le.0.) then
      if (x.lt.0.) then
        write (*,*) 'error in gser: x < 0 (module gamma_functions, src_seifert.f90)'
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
      write (*,*) 'error in gser: a too large, itmax too small'
      write (*,*) '  (module gamma_functions, src_seifert.f90)'
      gamser = 0.0e0
      return
    end if

    gamser = sum*exp(-x+a*log(x)-gln)

    return
  end subroutine gser

  real function coll_delta(p1,n)
    implicit none
    type(particle) :: p1
    integer        :: n

    coll_delta = gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)         &
         &                  / gfct((p1%nu+1.0  )/p1%mu)         &
         &        * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_geo+n)   &
         &        / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_geo+n)

    return
  end function coll_delta

  real function coll_delta_11(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_delta_11 = coll_delta(p1,n)

    return
  end function coll_delta_11

  real function coll_delta_22(p2,n)

    type(particle) :: p2
    integer        :: n

    coll_delta_22 = coll_delta(p2,n)

    return
  end function coll_delta_22

  real function coll_delta_12(p1,p2,n)

    type(particle) :: p1,p2
    integer        :: n

    coll_delta_12 = 2.0 * gfct((p1%b_geo+p1%nu+1.0)/p1%mu)               &
         &                       / gfct((p1%nu+1.0)/p1%mu)               &
         &                * gfct((p1%nu+1.0)/p1%mu)**(p1%b_geo)          &
         &                / gfct((p1%nu+2.0)/p1%mu)**(p1%b_geo)          &
         &              * gfct((p2%b_geo+p2%nu+1.0+n)/p2%mu)             &
         &                        /gfct((p2%nu+1.0  )/p2%mu)             &
         &                * gfct((p2%nu+1.0)/p2%mu)**(p2%b_geo+n)        &
         &                / gfct((p2%nu+2.0)/p2%mu)**(p2%b_geo+n)

    return
  end function coll_delta_12

  real function coll_theta(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_theta = gfct((2.0*p1%b_vel+2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &               / gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &          * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_vel)        &
         &          / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_vel)

    return
  end function coll_theta

  real function coll_theta_11(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_theta_11 = coll_theta(p1,n)

    return
  end function coll_theta_11

  real function coll_theta_22(p2,n)

    type(particle) :: p2
    integer        :: n

    coll_theta_22 = coll_theta(p2,n)

    return
  end function coll_theta_22

  real function coll_theta_12(p1,p2,n)

    type(particle) :: p1,p2
    integer        :: n

    coll_theta_12 = 2.0 * gfct((p1%b_vel+2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &
         &                    / gfct((2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &
         &              * gfct((p1%nu+1.0)/p1%mu)**(p1%b_vel)              &
         &              / gfct((p1%nu+2.0)/p1%mu)**(p1%b_vel)              &
         &           * gfct((p2%b_vel+2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &
         &                    / gfct((2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &
         &              * gfct((p2%nu+1.0)/p2%mu)**(p2%b_vel)              &
         &              / gfct((p2%nu+2.0)/p2%mu)**(p2%b_vel)

    return
  end function coll_theta_12

  function d_average_factor (parti)

  ! faktor zur berechnung des meann durchmessers von verallg. gammaverteilten hydrometeoren:
    ! gueltig fuer d = a_geo * x^b_geo
    ! berechnung des meann durchmessers: d_average = parti%b_geo * d_average_factor * (q/qn)**parti%b_geo

    implicit none

    real :: d_average_factor
    type(particle), intent(in) :: parti

    d_average_factor = &
         ( gfct( (parti%b_geo+parti%nu+1.0)/parti%mu ) / &
         gfct( (parti%nu+1.0)/parti%mu ) ) * &
         (gfct( (parti%nu+1.0)/parti%mu ) / gfct( (parti%nu+2.0)/parti%mu ) ) ** parti%b_geo

  end function d_average_factor

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
  i=0
    gen_sequence = (/(i,i=1,nprocess)/)
    if (.not. droplet_sedim) then
      where (gen_sequence==isedimcd)
        gen_sequence = 0
      end where
    end if
    call shuffle(gen_sequence)
    print *,gen_sequence
  end function gen_sequence
  
  subroutine shuffle(a)
    integer, intent(in out) :: a(:)
    integer :: i, rand, temp
    real :: x

    do i = size(a), 1, -1
       call random_number(x)
       rand = int(x * i) + 1
       temp = a(rand)
       a(rand) = a(i)
       a(i) = temp
    end do
  end subroutine

  subroutine resetvar(meteor,mass,number)
    type(particle),intent(in)        :: meteor
    real, dimension(:), intent(inout) :: mass
    real, dimension(:), intent(inout), optional :: number
    where (mass < qthres)
      mass = 0.
    end where
    if (present(number)) then
      number = max(0.,max(min(number,mass/meteor%x_min),mass/meteor%x_max))
    end if
  end subroutine resetvar
  subroutine initmcrp(level,firsttime)
  integer , intent(in) :: level
  logical, intent(inout) :: firsttime

  firsttime = .false.
 if (level>=3)      nprocess = nprocwarm
!  if (droplet_sedim) nprocess = nprocess + 1
 if (level==4)      nprocess = nprocess + nprocice


!Cloudwater properties
 cldw = PARTICLE('cloudw', &! HN: made variable for seeding experiments
         1, & !number of moments
         CCN, & !Number of droplets
         0.333333, & !.nu.....Breiteparameter der Verteil.
         0.666666, & !.mu.....Exp.-parameter der Verteil.
         2.60e-10, & !.x_max..maximale Teilchenmasse D=80e-6m
         4.20e-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
         1.24e-01, & !.a_geo..Koeff. Geometrie
         0.333333, & !.b_geo..Koeff. Geometrie = 1/3
         3.75e+05, & !.a_vel..Koeff. Fallgesetz
         0.666667, & !.b_vel..Koeff. Fallgesetz
         0.25     , &!.s_vel...Dispersion der Fallgeschw.
         0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
         0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
         2.0)        !.cap....Koeff. Kapazitaet

!Rainwater properties
 rain = PARTICLE('rain', &! HN: made variable for seeding experiments
         2, & !number of moments
         0, & !Number of droplets
         0.333333, & !.nu.....Breiteparameter der Verteil.
         0.666666, & !.mu.....Exp.-parameter der Verteil.
         2.60e-10, & !.x_max..maximale Teilchenmasse D=80e-6m
         4.20e-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
         1.24e-01, & !.a_geo..Koeff. Geometrie
         0.333333, & !.b_geo..Koeff. Geometrie = 1/3
         3.75e+05, & !.a_vel..Koeff. Fallgesetz
         0.666667, & !.b_vel..Koeff. Fallgesetz
         0.25     , &!.s_vel...Dispersion der Fallgeschw.
         0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
         0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
         2.0)        !.cap....Koeff. Kapazitaet
!Cloud ice properties
 ice = PARTICLE('ice', &
        2, & !number of moments
        0, & !
        0.000000, & !.nu.....Breiteparameter der Verteil.
        0.333333, & !.mu.....Exp.-parameter der Verteil.
        1.00e-07, & !.x_max..maximale Teilchenmasse D=???e-2m
        1.00e-12, & !.x_min..minimale Teilchenmasse D=200e-6m
        3.303633, & !.a_geo..Koeff. Geometrie
        0.476191, & !.b_geo..Koeff. Geometrie = 1/2.1
        2.77e+01, & !.a_vel..Koeff. Fallgesetz
        0.215790, & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
         0.25     , &!.s_vel...Dispersion der Fallgeschw.
        0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
        0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
        2.0)        !.cap....Koeff. Kapazitaet

!Snow properties

 snow =  PARTICLE('snow', & ! nach Andy Heymsfield (CRYSTAL-FACE)
         1,  & !number of moments
         1e2, &
        0.000000, & !.nu.....Breiteparameter der Verteil.
        0.333333, & !.mu.....Exp.-parameter der Verteil.
        1.00e-07, & !.x_max..maximale Teilchenmasse D=???e-2m
        1.00e-12, & !.x_min..minimale Teilchenmasse D=200e-6m
        3.303633, & !.a_geo..Koeff. Geometrie
        0.476191, & !.b_geo..Koeff. Geometrie = 1/2.1
        2.47e+02, & !.a_vel..Koeff. Fallgesetz
        0.333333, & !.b_vel..Koeff. Fallgesetz
         0.25     , &!.s_vel...Dispersion der Fallgeschw.
        0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
        0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
        2.0)        !.cap....Koeff. Kapazitaet

!Graupel properties
graupel = PARTICLE('graupel', & ! 'graupel im standard'
         1,  & !number of moments
         1e2, &
         1.000000, & !.nu.....Breiteparameter der Verteil.
         0.166666, & !.mu.....Exp.-parameter der Verteil.
         1.00e-04, & !.x_max..maximale Teilchenmasse
         2.60e-10, & !.x_min..minimale Teilchenmasse
         1.10e-01, & !.a_geo..Koeff. Geometrie
         0.300000, & !.b_geo..Koeff. Geometrie = 1/3.10
         7.64e+01,& !.a_vel..Koeff. Fallgesetz
         0.255200, & !.b_vel..Koeff. Fallgesetz
         0.25     , &!.s_vel...Dispersion der Fallgeschw.
         0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
         0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
         2.0)        !.cap....Koeff. Kapazitaet
         
  end subroutine initmcrp
end module mcrp
