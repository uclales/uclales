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

  use mpi_interface,  only : myid, double_scalar_par_max
  use defs, only : tmelt,alvl, alvi,rowt,roice, pi, Rm, cp,t_hn 
  use grid, only : dt,nstep,rkbeta,rkalpha, dxi, dyi ,dzi_t, nxp, nyp, nzp,nfpt, a_pexnr, pi0,pi1,a_rp, a_tp, th00, ccn,    &
       dn0, pi0,pi1, a_rt, a_tt,a_rpp, a_rpt, a_npp, a_npt, vapor, liquid, a_wp,      &
       a_theta, a_scr1, a_scr2, a_scr7, &
       a_ninucp, a_ninuct , & ! ice nuclei concentration 
       a_ricep , a_ricet  , & ! ice mixing ratio
       a_nicep , a_nicet  , & ! ice number concentration
       a_rsnowp, a_rsnowt , & ! snow mass
       a_nsnowp, a_nsnowt , & ! snow number
       a_rgrp,   a_rgrt,    & ! graupel mass
       a_ngrp,   a_ngrt,    & ! graupel number
       a_rhailp, a_rhailt,  & ! hail mass
       a_nhailp, a_nhailt,  & ! hail number
       prc_c, prc_r, prc_i, prc_s, prc_g, prc_h

  USE parallele_umgebung, ONLY: isIO,double_global_maxval,global_maxval, global_minval

  use thrm, only : thermo, fll_tkrs
!axel!  use stat, only : sflg, updtst
  use util, only : get_avg3, azero, sclrset
  implicit none

  logical, parameter :: droplet_sedim = .False., khairoutdinov = .False., turbulence = .False.,ice_multiplication = .TRUE.
  integer            :: nprocess,nprocwarm=5,nprocice=18

  integer,parameter  :: iwtrdff = 3,iauto = 1,iaccr = 2,isedimrd = 4,isedimcd = 5, &
       iicenucnr = 6, iicenuc = 7, ifreez=8, idep=9, imelt_ice=10,imelt_snow=11,imelt_grp=12, ised_ice=13, &
       iself_ice=14, icoll_ice_snow=15, icoll_ice_grp=16, icoll_snow_grp=17, &
       iriming_ice_cloud=18, iriming_snow_cloud=19, iriming_grp_cloud=20, &
       iriming_ice_rain=21, iriming_snow_rain=22,iriming_grp_rain=23
  real, dimension(:),allocatable :: convice,convliq
  integer, dimension(23) :: microseq = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23/)
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
  real, parameter :: rthres = 0.0          ! small number
  !real, parameter :: rthres = 1e-14        ! small number
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
  ! 
  real, parameter :: c_mult     = 3.5d8    !..splintering coefficient
  real, parameter :: t_mult_min = 265.0    !..min temp. splintering
  real, parameter :: t_mult_max = 270.0    !..max temp. splintering
  real, parameter :: t_mult_opt = 268.0    !..opt temp. splintering

  real, parameter :: e_ic  = 0.80              !..max. eff. for ice_cloud_riming
  real, parameter :: e_sc  = 0.80              !..max. eff. for snow_cloud_riming
  real, parameter :: e_gc  = 1.00              !..max. eff. for graupel_cloud_riming
  real, parameter :: alpha_spacefilling = 0.1  !..space filling coef (max. 0.68)

  real, parameter :: r_crit_ic = 1.000e-5 ! r-critical value for ice_cloud_riming
  real, parameter :: d_crit_ic = 150.0e-6 ! e-critical value for ice_cloud_riming
  real, parameter :: r_crit_ir = 1.000e-5 ! r-critical value for ice_rain_riming
  real, parameter :: d_crit_ir = 100.0e-6 ! e-critical value for ice_rain_riming
  real, parameter :: r_crit_sc = 1.000e-5 ! r-critical value for snow_cloud_riming
  real, parameter :: d_crit_sc = 150.0e-6 ! e-critical value for snow_cloud_riming
  real, parameter :: r_crit_sr = 1.000e-5 ! r-critical value for snow_rain_riming
  real, parameter :: d_crit_sr = 100.0e-6 ! e-critical value for snow_rain_riming
  real, parameter :: r_crit_gc = 1.000e-6 ! r-critical value for graupel_cloud_riming
  real, parameter :: d_crit_gc = 100.0e-6 ! e-critical value for graupel_cloud_riming
  real, parameter :: r_crit_c  = 1.000e-6 ! r-critical value else
  real, parameter :: r_crit    = 1.000e-9 ! r-critical value else
  real, parameter :: d_conv_sg = 200.0e-5 ! e-critical value
  real, parameter :: d_conv_ig = 200.0e-6 ! e-critical value
  real, parameter :: d_crit_c  = 10.00e-6 ! e-critical value for cloud_collection
  real, parameter :: d_crit_r  = 10.00e-6 ! e-critical value for cloud_collection
  real, parameter :: d_coll_c  = 40.00e-6 ! max value for cloud_coll_eff
  real, parameter :: r_crit_is = 1.000e-4 ! r-critical value for ice_selfcollection
  real, parameter :: d_crit_is = 50.00e-6 ! e-critical value for ice_selfcollection
  real, parameter :: d_conv_is = 75.00e-6 ! e-critical value for ice_selfcollection
  real :: &
    na_dust    = 162.e3,   & ! initial number density of dust [1/m³], phillips08 value 162e3
    na_soot    =  15.e6,   & ! initial number density of soot [1/m³], phillips08 value 15e6
    na_orga    = 177.e6,   & ! initial number density of organics [1/m3], phillips08 value 177e6
    ni_het_max = 100.0d3,   & ! max number of in between 1-10 per liter, i.e. 1d3-10d3
    ni_hom_max = 5000.0d3     ! number of liquid aerosols between 100-5000 per liter

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

contains

  !
  ! ---------------------------------------------------------------------
  ! MICRO: sets up call to microphysics
  !
  subroutine micro(level,istp)

    integer, intent (in) :: level

    logical :: debug = .false.
    REAL(KIND=8) :: wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,&
         &                  precmax,rhomax,rhomin,ncmax,nimax,dzmin,tmax,tmin
    integer, optional :: istp

    if (present(istp)) istep = istp


    if (mod(istep,100).eq.0) then
       debug = .true.
    else
       debug = .false.
    end if


    call fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,a_scr1,rs=a_scr2)
    select case (level) 
    case(2)
       if (droplet_sedim) call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,vapor,a_scr2,liquid,prc_c)
    case(3)
    if (mod(istep,100).eq.0.and.nstep.eq.3) then
    !IF (debug) THEN
      wmax  = double_global_maxval(a_wp)
      qvmax = double_global_maxval(vapor)
      qcmax = double_global_maxval(liquid)
      qrmax = double_global_maxval(a_rpp)
      nimax = double_global_maxval(a_npp)
      tmax  = double_global_maxval(a_scr1)
      IF (isIO()) THEN
      WRITE (*,*) "micro: output before mcrph"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','nr','tmax'
      WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,nimax,tmax
      ENDIF
    ENDIF
       call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,vapor,a_scr2,liquid,prc_c,a_rpp, &
            a_npp,a_rt,a_rpt,a_npt,a_scr7, prc_r)
    if (mod(istep,100).eq.0.and.nstep.eq.3) then
    !IF (debug) THEN
      wmax  = double_global_maxval(a_wp)
      qvmax = double_global_maxval(vapor)
      qcmax = double_global_maxval(liquid)
      qrmax = double_global_maxval(a_rpp)
      nimax = double_global_maxval(a_npp)
      tmax  = double_global_maxval(a_scr1)
      IF (isIO()) THEN
      WRITE (*,*) "micro: output after mcrph"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','nr','tmax'
      WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,nimax,tmax
      ENDIF
    ENDIF
    case(4)
    IF (debug) THEN
      wmax  = double_global_maxval(a_wp)
      qvmax = double_global_maxval(vapor)
      qcmax = double_global_maxval(liquid)
      qrmax = double_global_maxval(a_rpp)
      qimax = double_global_maxval(a_ricep)
      qsmax = double_global_maxval(a_rsnowp)
      qgmax = double_global_maxval(a_rgrp)
      nimax = double_global_maxval(a_nicep)
      tmax  = double_global_maxval(a_scr1)
      IF (isIO()) THEN
      WRITE (*,*) "micro: output before mcrph"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','ni','tmax'
      WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,nimax,tmax
      ENDIF
    ENDIF
       call mcrph(level,nzp,nxp,nyp,dn0,a_pexnr,pi0,pi1,a_tp,a_tt,a_scr1,vapor,a_scr2,liquid,prc_c,a_rpp, &
            a_npp,a_rt,a_rpt,a_npt,a_scr7, prc_r, &
            a_ninuct,a_ricet,a_nicet,a_rsnowt,a_rgrt,&
            a_ninucp,a_ricep,a_nicep,a_rsnowp,a_rgrp, &
            prc_i, prc_s, prc_g)
    IF (debug) THEN
      wmax  = double_global_maxval(a_wp)
      qvmax = double_global_maxval(vapor)
      qcmax = double_global_maxval(liquid)
      qrmax = double_global_maxval(a_rpp)
      qimax = double_global_maxval(a_ricep)
      qsmax = double_global_maxval(a_rsnowp)
      qgmax = double_global_maxval(a_rgrp)
      nimax = double_global_maxval(a_nicep)
      tmax  = double_global_maxval(a_scr1)
      IF (isIO()) THEN
      WRITE (*,*) "micro: output after mcrph"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','ni','tmax'
      WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,nimax,tmax
      ENDIF
    ENDIF
    case(5)
    IF (debug) THEN
      wmax  = double_global_maxval(a_wp)
      qvmax = double_global_maxval(vapor)
      qcmax = double_global_maxval(liquid)
      qrmax = double_global_maxval(a_rpp)
      qimax = double_global_maxval(a_ricep)
      qsmax = double_global_maxval(a_rsnowp)
      qgmax = double_global_maxval(a_rgrp)
      qhmax = double_global_maxval(a_rhailp)
      nimax = double_global_maxval(a_nicep)
      tmax  = double_global_maxval(a_scr1)
      IF (isIO()) THEN
      WRITE (*,*) "micro: output before mcrph_sb"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','ni','tmax'
      WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,nimax,tmax
      ENDIF
      tmax  = double_global_maxval(a_tt)
      qcmax = double_global_maxval(a_rt)
      qrmax = double_global_maxval(a_rpt)
      qimax = double_global_maxval(a_ricet)
      qsmax = double_global_maxval(a_rsnowt)
      qgmax = double_global_maxval(a_rgrt)
      qhmax = double_global_maxval(a_rhailt)
      IF (isIO()) THEN
      WRITE(*,'(A10,10A11)')   '   ', 'ttend','qttend','qrtend','qitend','qstend','qgtend', 'qhtend'
      WRITE(*,'(A10,10E11.3)') '   ', tmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax
      ENDIF
    ENDIF
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
    IF (debug) THEN
      wmax  = double_global_maxval(a_wp)
      qvmax = double_global_maxval(vapor)
      qcmax = double_global_maxval(liquid)
      qrmax = double_global_maxval(a_rpp)
      qimax = double_global_maxval(a_ricep)
      qsmax = double_global_maxval(a_rsnowp)
      qgmax = double_global_maxval(a_rgrp)
      qhmax = double_global_maxval(a_rhailp)
      nimax = double_global_maxval(a_nicep)
      tmax  = double_global_maxval(a_scr1)
      IF (isIO()) THEN
      WRITE (*,*) "micro: output after mcrph_sb"
      WRITE(*,'(10x,a)') 'Maximum Values:'
      WRITE(*,'(A10,10A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','ni','tmax'
      WRITE(*,'(A10,10E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,nimax,tmax
      ENDIF
      tmax  = double_global_maxval(a_tt)
      qcmax = double_global_maxval(a_rt)
      qrmax = double_global_maxval(a_rpt)
      qimax = double_global_maxval(a_ricet)
      qsmax = double_global_maxval(a_rsnowt)
      qgmax = double_global_maxval(a_rgrt)
      qhmax = double_global_maxval(a_rhailt)
      IF (isIO()) THEN
      WRITE(*,'(A10,10A11)')   '   ', 'ttend','qttend','qrtend','qitend','qstend','qgtend', 'qhtend'
      WRITE(*,'(A10,10E11.3)') '   ', tmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax
      ENDIF
    ENDIF
    end select

  end subroutine micro
  !
  ! ---------------------------------------------------------------------
  ! MCRPH: calls microphysical parameterization 
  !

  subroutine mcrph(level,n1,n2,n3,dn0,exner,pi0,pi1,thl,tlt,tk,vapor,rsat,rcld,prc_c, &
       rp, np, rtt,rpt,npt,dissip, prc_r, &
       ninuct,ricet,nicet,rsnowt,rgrpt,ninucp,ricep,nicep,rsnowp,rgrpp,si, &
       prc_i, prc_s, prc_g)

    integer, intent (in) :: level,n1,n2,n3
    real, dimension(n1)      , intent (in)             :: dn0,pi0,pi1  !Density and mean pressures
    real, dimension(n1,n2,n3), intent (inout)          :: exner, & !exner function
         thl,& !Theta_l
         tk, & !temperature
         tlt,&! theta_l tendency
         vapor, & !water vapor
         rsat, & !saturation mixing ration
         rcld, &  !Cloud water
         prc_c
    real, dimension(n1,n2,n3), intent (inout),optional :: np,&!rain drop number
         rp, &!rain water
         rtt,& !Total water tendency
         rpt,&!rain water tendency
         npt,&!rain droplet number tendency
         dissip, &
         prc_r

    real, dimension(n1,n2,n3), intent (inout),optional :: ninuct,ricet,nicet,rsnowt,rgrpt,&
         ninucp,ricep,nicep,rsnowp,rgrpp,&
         si
    real, dimension(n1,n2,n3), intent (inout),optional :: prc_i, prc_s, prc_g

    real, dimension(n1) :: tl,temp,rv,rc,nrain,rrain,ninuc,nice, nin_active,rice,nsnow,rsnow,ngrp,rgrp,rs,s_i,r1
    integer :: i, j,n
!as    logical,save  :: firsttime = .true.
    !     real :: dtrk
    ! 
!as    if(firsttime) call initmcrp(level,firsttime)
    allocate(convice(n1),convliq(n1))
    do j=3,n3-2
       do i=3,n2-2
          call resetvar(rain,rp(1:n1,i,j),np(1:n1,i,j))
          tl = thl(1:n1,i,j)
          temp = tk(1:n1,i,j)
          rv = vapor(1:n1,i,j)
          rc = rcld(1:n1,i,j)
          rs = rsat(1:n1,i,j)
          rrain = rp(1:n1,i,j)
          nrain = np(1:n1,i,j)
          convliq = alvl/cp*(pi0+pi1+exner(1:n1,i,j))/cp
          if (level == 4) then
             call resetvar(ice,ricep(1:n1,i,j),nicep(1:n1,i,j))
             call resetvar(snow,rsnowp(1:n1,i,j))
             call resetvar(graupel,rgrpp(1:n1,i,j))
             s_i(k) =  R_d &
                  & * dn0(k) * rv(k) &
                  & * tk(k) &
                  & / e_es(dble(tk(k))) - 1.0
 

             ninuc = ninucp(1:n1,i,j)

             rice = ricep(1:n1,i,j)
             nice = nicep(1:n1,i,j)
             rsnow = rsnowp(1:n1,i,j)
             where (rsnow>0)
                nsnow = snow%nr
             elsewhere
                nsnow = 0
             end where
             rgrp = rgrpp(1:n1,i,j)
             where (rgrp>0)
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
                call resetvar(cldw,rc)
                call resetvar(rain,rrain,nrain)
                call wtr_dff_SB(n1,dn0,rrain,nrain,rc,rs,rv,tl,temp)
             case(iauto)
                call resetvar(cldw,rc)
                call resetvar(rain,rrain,nrain)
                call auto_SB(n1,dn0,rc,rrain,nrain,tl,dissip(1:n1,i,j))
             case(iaccr)
                call resetvar(cldw,rc)
                call resetvar(rain,rrain,nrain)
                call accr_SB(n1,dn0,rc,rrain,nrain,tl,dissip(1:n1,i,j))
             case(isedimrd)
                call resetvar(rain,rrain,nrain)
                call sedim_rd(n1,dt,dn0,rrain,nrain,prc_r(1:n1,i,j))
             case(isedimcd)
                call resetvar(cldw,rc)
                call sedim_cd(n1,dt,tl,rc,prc_c(1:n1,i,j))
             case(iicenucnr)
!                 call n_icenuc(n1,ninuc,nin_active,temp,rv,s_i)
             case(iicenuc)
                where (ninuc<0.) ninuc = 0.
!                 call ice_nucleation(n1,ninuc,rice,nice,s_i,tl,temp)
                 call ice_nucleation_homhet(n1,ice,rv,rc,rice,nice,rsnow, nsnow,s_i,tl,temp)
                 

             case(ifreez)
                call resetvar(cldw,rc)
                call resetvar(rain,rrain,nrain)

                call cloud_freeze(n1,rc,ninuc,rice,nice,tl,temp)
                call rain_freeze(n1,rrain,nrain,ninuc,rice,nice,rgrp,temp)
             case(idep)
! if (i==15 .and. j==15) then
!   print *, s_i(10:17)
! end if
                call deposition(n1,ice,ninuc,rice,nice,rv,tl,temp,s_i)
                call deposition(n1,snow,ninuc,rsnow,nsnow,rv,tl,temp,s_i)
                call deposition(n1,graupel,ninuc,rgrp,ngrp,rv,tl,temp,s_i)
             case(imelt_ice)
                call resetvar(ice,rice,nice)
                call melting(n1,ice,rice,tl,nice,rc,rrain,nrain,temp)
             case(imelt_snow)
                call resetvar(snow,rsnow)
                call melting(n1,snow,rsnow,tl,nsnow,rc,rrain,nrain,temp)
             case(imelt_grp)
                call resetvar(graupel,rgrp)
                call melting(n1,graupel,rgrp,tl,ngrp,rc,rrain,nrain,temp)
             case(iself_ice)
                call resetvar(ice,rice,nice)
                call ice_selfcollection(n1,ice,snow,rice,rsnow,nice,dn0,temp,d_conv_is,r_crit_is,d_crit_is)
             case(icoll_ice_snow)
                call resetvar(ice,rice,nice)
                call resetvar(snow,rsnow)
                call ice_collection(n1,ice,snow,rice,nice,rsnow,temp,r_crit_is)
             case(icoll_ice_grp)
                call resetvar(ice,rice,nice)
                call resetvar(graupel,rgrp)
                call ice_collection(n1,ice,graupel,rice,nice,rgrp,temp,r_crit_is)

             case(icoll_snow_grp)
                call resetvar(ice,rice,nice)
                call resetvar(snow,rsnow)
                call resetvar(graupel,rgrp)
                call ice_collection(n1,snow,graupel,rsnow,nsnow,rgrp,temp,r_crit_is)
             case(iriming_ice_cloud)
                call resetvar(cldw,rc)
                call resetvar(ice,rice,nice)
                call ice_cloud_riming(n1,cldw,ice,rc,rice,nice,rgrp,tl,temp,d_coll_c,r_crit_c,d_crit_c,r_crit_ic,d_crit_ic,d_conv_ig,e_ic)
             case(iriming_snow_cloud)
                call resetvar(snow,rsnow)
                call resetvar(cldw,rc)
                call ice_cloud_riming(n1,cldw,snow,rc,rsnow,nsnow,rgrp,tl,temp,d_coll_c,r_crit_c,d_crit_c,r_crit_sc,d_crit_sc,d_conv_sg,e_sc)
             case(iriming_grp_cloud)
                call resetvar(graupel,rgrp)
                call resetvar(cldw,rc)
                r1 = 0.
                call ice_cloud_riming(n1,cldw,graupel,rc,rgrp,ngrp,r1,tl,temp,d_coll_c,r_crit_c,d_crit_c,r_crit_gc,d_crit_gc,d_conv_sg,e_gc)
                rgrp = rgrp + r1
             case(iriming_ice_rain)
                call resetvar(ice,rice,nice)
                call resetvar(rain,rrain,nrain)
                call ice_rain_riming(n1,rain,ice    ,rrain,nrain,rice ,nice ,rgrp,dn0,temp,r_crit,d_crit_r, r_crit_ir,d_crit_ir)
                !             call ice_rain_riming(n1,rain,ice    ,rrain,nrain,rice ,nice ,rgrp,temp,d_crit_r,1e-2,1e-2,1e-2)
             case(iriming_snow_rain)
                call resetvar(snow,rsnow)
                call resetvar(rain,rrain,nrain)
                call ice_rain_riming(n1,rain,snow   ,rrain,nrain,rsnow,nsnow,rgrp,dn0,temp,r_crit,d_crit_r, r_crit_sr,d_crit_sr)
             case(iriming_grp_rain)
                call resetvar(graupel,rgrp)
                call resetvar(rain,rrain,nrain)
                r1 = 0.
                call ice_rain_riming(n1,rain,graupel,rrain,nrain,rgrp ,ngrp,r1,dn0,temp,r_crit,d_crit_r, r_crit,0.)
                rgrp = rgrp + r1
             case(ised_ice)
                call resetvar(ice,rice,nice)
                call resetvar(snow,rsnow)
                call resetvar(graupel,rgrp)
                call sedimentation (n1,ice, rice,nice,rrate = prc_i(1:n1,i,j))
                call sedimentation (n1,snow, rsnow, rrate = prc_s(1:n1,i,j))
                call sedimentation (n1,graupel, rgrp, rrate = prc_g(1:n1,i,j))
                prc_i(1:n1,i,j) = prc_i(1:n1,i,j)
                prc_s(1:n1,i,j) = prc_s(1:n1,i,j)
                prc_g(1:n1,i,j) = prc_g(1:n1,i,j)            
             case default
             end select
          end do
          call resetvar(cldw,rc)
          call resetvar(rain,rrain,nrain)
          if (level==4) then
             call resetvar(ice,rice,nice)
             call resetvar(snow,rsnow)
             call resetvar(graupel,rgrp)
          end if

          tlt(2:n1,i,j) =  tlt(2:n1,i,j)+      (tl(2:n1) - thl(2:n1,i,j))/dt
          rtt(2:n1,i,j) = rtt(2:n1,i,j) +(rv(2:n1) - vapor(2:n1,i,j))/dt + (rc(2:n1) - rcld(2:n1,i,j))/dt
          rpt(2:n1,i,j) = max(rpt(2:n1,i,j) +(rrain(2:n1) - rp(2:n1,i,j))/dt,-rp(2:n1,i,j)/dt)
          npt(2:n1,i,j) = max(npt(2:n1,i,j) +(nrain(2:n1) - np(2:n1,i,j))/dt,-np(2:n1,i,j)/dt)

          if (level == 4) then
             ninuct(2:n1,i,j) = ninuct(2:n1,i,j) + (ninuc(2:n1)-ninucp(2:n1,i,j))/dt !max(ninuct(2:n1,i,j) + (ninuc(2:n1) - ninucp(2:n1,i,j))/dtrk,-ninucp(2:n1,i,j)/dt)
             ricet(2:n1,i,j)  = max(ricet(2:n1,i,j)  + (rice(2:n1) - ricep(2:n1,i,j))/dt,-ricep(2:n1,i,j)/dt)
             nicet(2:n1,i,j)  = nicet(2:n1,i,j)  + (nice(2:n1) - nicep(2:n1,i,j))/dt
             rsnowt(2:n1,i,j) = max(rsnowt(2:n1,i,j) + (rsnow(2:n1) - rsnowp(2:n1,i,j))/dt,-rsnowp(2:n1,i,j)/dt)
             rgrpt(2:n1,i,j)   = max(rgrpt(2:n1,i,j)   + (rgrp(2:n1) - rgrpp(2:n1,i,j))/dt,-rgrpp(2:n1,i,j)/dt)
          end if

       end do
    end do
    deallocate(convice,convliq)
    !print *,maxval(s_i),maxloc(s_i)
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

    real, parameter     :: c_Nevap = 1.
    integer             :: k
    real                :: Xp, Dp, G, S, cerpt, cenpt
    do k=2,n1
       if (rp(k) > 0) then
          Xp = rp(k)/ (np(k)+eps0)
          Dp = ( Xp / prw )**(1./3.)
          G = 1. / (1. / (dn0(k)*rs(k)*Dv) + &
               alvl*(alvl/(Rm*tk(k))-1.) / (Kt*tk(k)))
          S = rv(k)/rs(k) - 1.

          if (S < 0) then
             cerpt = 2. * pi * Dp * G * S * np(k)
             cerpt = min (cerpt, rv(k)/dt)
             cenpt = c_Nevap*cerpt * np(k) / rp(k)
             np(k)=np(k) + cenpt*dt
             rp(k)=rp(k) + cerpt*dt
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
    real    :: k_au0,k_au, Xc, Dc, au, tau, phi, kc_alf, kc_rad, kc_sig, k_c, Re, epsilon, l
    real, dimension(n1) :: Resum, Recnt

    !
    ! Calculate the effect of turbulence on the autoconversion/
    ! accretion rate using equation (6) in Seifert et al. (2009)
    !

    nu_c = 1.0 !! cldw%nu!1*dn0(k)*rc(k)+cldw%nu2

    k_au0  = kc_0 / (20.*cldw%x_max) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2

    do k=2,n1-1
       if (rc(k) > 0.) then
          Xc = rc(k)/(cldw%nr+eps0)
          k_c = kc_0

          if (turbulence) then
             kc_alf = ( kc_a1 + kc_a2 * nu_c )/ ( 1. + kc_a3 * nu_c )
             kc_rad = ( kc_b1 + kc_b2 * nu_c )/ ( 1. + kc_b3 * nu_c )
             kc_sig = ( kc_c1 + kc_c2 * nu_c )/ ( 1. + kc_c3 * nu_c )
             call azero(n1,Recnt,Resum)
             Dc = ( Xc / prw )**(1./3.)  ! mass mean diameter cloud droplets in m
             !
             ! Calculate the mixing length, dissipation rate and Taylor-Reynolds number
             !
             l = csx*((1/dxi)*(1/dyi)*(1/dzi_t(k)))**(1./3.)
             epsilon = min(diss(k),0.06)
             Re = (6./11.)*((l/ce)**(2./3))*((15./(1.5e-5))**0.5)*(epsilon**(1./6) )
             !
             ! Dissipation needs to be converted to cm^2/s^3, which explains the factor
             ! 1.e4 with which diss(k) is multiplied
             !
             k_au = k_au0 * (1. + epsilon *1.e4* (Re*1.e-3)**(0.25) &
                  * (kc_bet + kc_alf * exp( -1.* ((((Dc/2.)*1.e+6-kc_rad)/kc_sig)**2) )))
             !print *,'enhancement factor = ', k_c/(9.44e+9)
             !
             ! Calculate conditional average of Re i.e., conditioned on cloud/rain water
             !
             Resum(k) = Resum(k)+ Re
             Recnt(k) = Recnt(k)+ 1
          else
             k_au = k_au0
          end if

          ! which one is correct here?
          !au = k_au * rc(k)**2 * Xc**2
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
          au    = min(au,rc(k)/dt)
          rp(k) = rp(k) + au*dt
          rc(k) = rc(k) - au*dt
          tl(k) = tl(k) + convliq(k)*au*dt
          np(k) = np(k) + au/cldw%x_max*dt
       end if
    end do
    
    !if (isIO().and.MAXVAL(rc).gt.0) THEN
    !if (isIO()) THEN
    !  WRITE(0,'(a,15e15.4)') '  k_au = ',k_au0
    !end if

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

    real, parameter :: k_r0 = 5.78
    real, parameter :: k_rr = 4.33
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

             ! correct??
             ac  = k_r * rc(k) * rp(k) * phi * sqrt(rho_0*dn0(k))

             !
             ! Khairoutdinov and Kogan
             !
             if (khairoutdinov) then
                ac = Cac * (rc(k) * rp(k))**Eac
             end if
             !
             ac    = min(ac, rc(k)/dt)
             rp(k) = rp(k) + ac*dt
             rc(k) = rc(k) - ac*dt
             tl(k) = tl(k) + convliq(k)*ac*dt

          end if
          sc = k_rr * np(k) * rp(k) * sqrt(rho_0*dn0(k))
          sc = min(sc, np(k)/dt)
          np(k) = np(k) - sc*dt
       end if
    end do

  end subroutine accr_SB
  subroutine sedim_rd(n1,dt,dn0,rp,np,rrate)
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
    real    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz, cmax
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
       cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzi_t(k)*dt
       cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzi_t(k)*dt
    end do

    !cmax = MAXVAL(cr)
    !if (cmax.gt.1.0) then
    !  WRITE(0,'(a,i5,a,f8.2)') 'sedim: myid = ',myid,' meteor = rain,  cmax = ',cmax
    !end if

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
          tot = tot + dn0(kk)*(np(kk)+nslope(kk)*(1.-cc))*cc/dzi_t(kk)
          zz  = zz + 1./dzi_t(kk)
          kk  = kk + 1
          cc  = min(1.,cn(kk) - zz*dzi_t(kk))
       enddo
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
       rfl(k) = -tot /dt

       kp1=k+1
       flxdiv = (rfl(kp1)-rfl(k))*dzi_t(k)/dn0(k)
       rp(k) = rp(k)-flxdiv*dt
       np(k) = np(k)-(nfl(kp1)-nfl(k))*dzi_t(k)/dn0(k)*dt
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
    rfl(n1) = 0.
    do k=n1-1,2,-1
       Xc = rc(k) / (cldw%nr+eps0)
       Dc = ( Xc / prw )**(1./3.)
       vc = min(c*(Dc*0.5)**2 * exp(4.5*(log(sgg))**2),1./(dzi_t(k)*dt))
       !               vc = min(c*(Dc*0.5)**2 * exp(5*(log(1.3))**2),1./(dzi_t(k)*dt))
       rfl(k) = - rc(k) * vc
       !
       kp1=k+1
       flxdiv = (rfl(kp1)-rfl(k))*dzi_t(k)
       rc(k) = rc(k)-flxdiv*dt
       tl(k) = tl(k)+flxdiv*convliq(k)*dt
       rrate(k)    = -rfl(k)
    end do

  end subroutine sedim_cd
  subroutine n_icenuc(n1,nin,nin_active,tk,rv,s_i)
    integer, intent(in) :: n1
    real, intent(inout), dimension (n1) :: nin
    real, intent(in) , dimension(n1) :: tk,s_i,rv,nin_active
    integer :: k
    real :: fact, supsat
    real :: niner
    fact = (1-exp(-dt/timenuc))

    do k=1,n1
      if (s_i(k) > 0) then
!         niner = n_ice_meyers_contact(tk(k),min(supsat,0.25))
!         nin(k)  = nin(k) + (n_ice_meyers_contact(tk(k),min(supsat,0.25))-nin(k))*fact
!         nin(k) = min(nin(k) + niner*dt/timenuc,niner)
        nin(k)  = n_ice_meyers_contact(tk(k),min(s_i(k)*100.,0.25))
! print *, 'nicenuc',k,nin(k)
!         if (s_i(k)/rv(k)>0.05) then
!           nin(k) = max(0.,nin_set -max(nin_active(k),0.))
!         else
!           nin(k) = 0.
!         end if
      end if
    end do

  end subroutine n_icenuc

  subroutine ice_nucleation(n1,nin,rice,nice,s_i,tl,tk)
    integer,intent(in) :: n1
    real, intent(inout),dimension(n1) :: nin,nice, rice, tl,s_i
    real,intent(in), dimension(n1) :: tk
    real :: nuc_n, nuc_r
    integer :: k

    do k=2,n1
       if (tk(k) < tmelt .and. s_i(k) > 0.0) then
          nuc_n = nin(k) 

          if (nuc_n>eps0) then
         nuc_r = nuc_n * ice%x_min
             nuc_n   = nuc_r / ice%x_min        
             nin(k)  = nin(k)  - nuc_n
             nice(k) = nice(k) + nuc_n
             rice(k) = rice(k) + nuc_r
             tl(k)   = tl(k)   + convice(k)*nuc_r
! print *, 'icen',k, tk(k), nin(k), s_i(k), nuc_r, nuc_n
          end if
       endif
    end do
  end subroutine ice_nucleation

  subroutine ice_nucleation_homhet(n1,ice,rv,rcloud,rice,nice,rsnow, nsnow,s_i,tl,tk)
  !*******************************************************************************
  !                                                                              *
  ! berechnung der nukleation der wolkeneispartikel                              *
  !                                                                              *
  ! nucleation scheme is based on the papers:                                    *
  !                                                                              *
  ! "a parametrization of cirrus cloud formation: homogenous                     *
  ! freezing of supercooled aerosols" by b. kaercher and                         *
  ! u. lohmann 2002 (kl02 hereafter)                                             *
  !                                                                              *
  ! "physically based parameterization of cirrus cloud formation                 *
  ! for use in global atmospheric models" by b. kaercher, j. hendricks           *
  ! and u. lohmann 2006 (khl06 hereafter)                                        *
  !                                                                              *
  !*******************************************************************************
! 
!   use globale_variablen,  only: loc_ix, loc_iy, loc_iz, t, p, q,                &
!     &                           q_cloud,q_ice,q_rain,q_snow,q_graupel,          &
!     &                           n_ice, n_snow, p_g, t_g, rho_g,w,rho,dz3d,dt,   &
!     &                           s_i,s_w,dsidz,dt0dz,                            &
!     &                           dqdt,speichere_dqdt,                &
!     &                           cond_neu_sb, evap_neu_sb, speichere_precipstat
!   use konstanten,         only: pi,r,cv,cp,wolke_typ
!   use parallele_umgebung, only: isio,abortparallel
!   use initialisierung,    only: t_0,p_0,rho_0,dichte
!   use geometrie,          only: x3_x3

  implicit none
    integer,intent(in) :: n1
    type(particle), intent(in) :: ice
    real, intent(inout),dimension(n1) :: rv,rcloud, nice, rice, nsnow, rsnow, tl,s_i
    real,intent(in), dimension(n1) :: tk
    real :: nuc_n, nuc_r
    integer :: k

  ! locale variablen
  real            :: t_a,p_a,rho_a,q_d,e_d
  real            :: q_i,n_i,x_i,r_i
  real            :: ndiag, n_m, n_f

  real, parameter :: eps  = 1.d-20


  ! coeffs of meyer formula
  real, parameter :: a_md = -0.639
  real, parameter :: b_md = 12.960
  real, parameter :: n_m0 = 1.0d+3

  ! some constants needed for kaercher and lohmann approach
  real, parameter :: &
    r_0     = 0.25e-6          , &    ! aerosol particle radius prior to freezing
    alpha_d = 0.5              , &    ! deposition coefficient (kl02; spichtinger & gierens 2009)
    k_b     = 1.38065e-23      , &    ! boltzmann constant [j/k]
    m_w     = 18.01528e-3      , &    ! molecular mass of water [kg/mol]
    m_a     = 28.96e-3         , &    ! molecular mass of air [kg/mol]
    n_avo   = 6.022e23         , &    ! avogadro number [1/mol]
    ma_w    = m_w / n_avo      , &    ! mass of water molecule [kg]
    grav    = 9.81             , &    ! acceleration of gravity [m/s2]
    p0      = 1013.25e2        , &
    svol    = ma_w / roice          ! specific volume of a water molecule in ice

  real(kind=8)  :: si, e_si
  real(kind=8)  :: ni_hom,ri_hom,mi_hom
  real(kind=8)  :: v_th,n_sat,flux,phi,d_v,cool,tau,delta,w_p,scr
  real(kind=8)  :: ctau, tau_g,acoeff(3),bcoeff(2), ri_dot
  real(kind=8)  :: kappa,sqrtkap,ren,r_imfc,r_im,r_ik,ri_0
  real(kind=8)  :: tgrow,ri_max,beta,xj,dxj,xmid,fmid,nimax
  real(kind=8)  :: xt,xs
  integer       :: jj,ss,tt

  logical :: use_hetnuc, use_homnuc, use_wp

  real(kind=8), dimension(3) :: infrac

!   real, allocatable :: nuc_n(:), nuc_q(:)
!   allocate (nuc_n(1:n1))
!   allocate (nuc_q(1:n1))

!   nuc_n = 0.0d0
!   nuc_q = 0.0d0

  ! with an upper limit of n_het_max

!         if (nuc_typ.eq.4) then
          na_dust    = 162.e4 ! number density of dust [1/m³], phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], phillips08 value 15e6
          na_orga    = 177.e5 ! number density of organics [1/m3], phillips08 value 177e6
!         if (nuc_typ.eq.6) then
!           ! increase dust by factor 100, reduce organics by factor of 10
!           na_dust    = 162.e5 ! number density of dust [1/m³], phillips08 value 162e3
!           na_soot    =  15.e6 ! number density of soot [1/m³], phillips08 value 15e6
!           na_orga    = 177.e5 ! number density of organics [1/m3], phillips08 value 177e6
      do k = 2, n1
        if (tk(k) < tmelt .and. s_i(k) > 0.0  &
          & .and. ( nice(k) < ni_het_max ) )then
          if (rcloud(k) > 0.0) then
            ! immersion freezing at water saturation
            xt = (274. - tk(k))  / ttstep
            xt = min(xt,ttmax-1.)
            tt = int(xt)
            infrac(1) = (tt+1-xt) * afrac_dust(tt,99) + (xt-tt) * afrac_dust(tt+1,99)
            infrac(2) = (tt+1-xt) * afrac_soot(tt,99) + (xt-tt) * afrac_soot(tt+1,99)
            infrac(3) = (tt+1-xt) * afrac_orga(tt,99) + (xt-tt) * afrac_orga(tt+1,99)

          else

            ! calculate indices used for look-up tables
            xt = (274. - tk(k))  / ttstep
            xs = (100*s_i(k)) / ssstep
            xt = min(xt,ttmax-1.)
            xs = min(xs,ssmax-1.)
            tt = int(xt)
            ss = int(xs)

            ! bi-linear interpolation in look-up tables
            infrac(1) = (tt+1-xt) * (ss+1-xs) * afrac_dust(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_dust(tt+1,ss  ) &
                      + (tt+1-xt) * (xs-ss)   * afrac_dust(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_dust(tt+1,ss+1)
            infrac(2) = (tt+1-xt) * (ss+1-xs) * afrac_soot(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_soot(tt+1,ss  ) &
                      + (tt+1-xt) * (xs-ss)   * afrac_soot(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_soot(tt+1,ss+1)
            infrac(3) = (tt+1-xt) * (ss+1-xs) * afrac_orga(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_orga(tt+1,ss  ) &
                      + (tt+1-xt) * (xs-ss)   * afrac_orga(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_orga(tt+1,ss+1)
          endif !end if ql>0

          ndiag = na_dust * infrac(1) + na_soot * infrac(2) + na_orga * infrac(3)
          ndiag = min(ndiag,ni_het_max)

          nuc_n = max( ndiag - (nice(k)+nsnow(k)),0.)

          nuc_r = min(nuc_n * ice%x_min, rv(k))

          nice(k) = nice(k) + nuc_n
          rice(k) = rice(k) + nuc_r
          rv(k)    = rv(k)     - nuc_r
          tl(k) = tl(k)+ convice(k) * nuc_r

        endif !if temp & supersaturation OK

      end do

  ! homogeneous nucleation using khl06 approach
!   if (use_homnuc) then
!     do k = 1, loc_iz
!       do j = 1, loc_iy
!         do i = 0, loc_ix
! 
!           p_a  = p_0(i,j,k)
!           t_a  = t_0(i,j,k)
!           e_si = e_es(t_a)
!           si   = q(i,j,k) * r_d * t_a / e_si
! 
!           ! critical supersaturation for homogenous nucleation, cf. eq. (1) of kb08
!           scr  = 2.349 - t_a / 259.00
! 
!           if (si > scr  &
!             .and. n_ice(i,j,k) < ni_hom_max ) then
! 
!             n_i = n_ice(i,j,k) ! + n_snow(i,j,k)
!             q_i = q_ice(i,j,k) ! + q_snow(i,j,k)
!             x_i = min(max(q_i/(n_i+eps),ice%x_min),ice%x_max)
!             !r_i = 0.5 * ice%a_geo * x_i**ice%b_geo
!             r_i = (x_i/(4./3.*pi*rho_ice))**(1./3.)
! 
!             d_v     = 0.211e-4 * (t_a/t_3)**1.94 * (p0/p_a)
!             v_th    = sqrt( 8.*k_b*t_a/(pi*ma_w) )
!             flux    = alpha_d * v_th/4.
!             n_sat   = e_si / (k_b*t_a)
! 
!             ! coeffs of supersaturation equation
!             acoeff(1) = (l_ed * grav) / (cp * r_d * t_a**2) - grav/(r_l * t_a)
!             acoeff(2) = 1.0/n_sat
!             acoeff(3) = (l_ed**2 * m_w * ma_w)/(cp * p_a * t_a * m_a)
! 
!             ! coeffs of depositional growth equation
!             bcoeff(1) = flux * svol * n_sat * (si - 1.)
!             bcoeff(2) = flux / d_v
! 
!             if (use_wp) then
!               ! pre-existing ice crystals included as reduced updraft speed
!               ri_dot = bcoeff(1) / (1. + bcoeff(2) * r_i)
!               r_ik   = (4 * pi) / svol * n_i * r_i**2 * ri_dot
!               w_p    = (acoeff(2) + acoeff(3) * si)/(acoeff(1) * si) * r_ik  ! khl06 eq. 19
!               w_p    = max(dble(w_p),0.d0)
!             else
!               w_p    = 0.d0
!             end if
! 
!             if (w(i,j,k) > w_p) then   ! homogenous nucleation event
! 
!               ! timescales of freezing event (see kl02, rm05, khl06)
!               cool    = grav / cp * w(i,j,k)
!               ctau    = t_a * ( 0.004*t_a - 2. ) + 304.4
!               tau     = 1.0 / (ctau * cool)                       ! freezing timescale, eq. (5)
!               delta   = (bcoeff(2) * r_0)                         ! dimless aerosol radius, eq.(4)
!               tau_g   = (bcoeff(1) / r_0) / (1 + delta)           ! timescale for initial growth, eq.(4)
!               phi     = acoeff(1)*si / ( acoeff(2) + acoeff(3)*si) * (w(i,j,k) - w_p)
! 
!               ! monodisperse approximation following khl06
!               kappa   = 2. * bcoeff(1) * bcoeff(2) * tau / (1.+ delta)**2  ! kappa, eq. 8 khl06
!               sqrtkap = sqrt(kappa)                                        ! root of kappa
!               ren     = 3. * sqrtkap / ( 2. + sqrt(1.+9.*kappa/pi) )       ! analy. approx. of erfc by rm05
!               r_imfc  = 4. * pi * bcoeff(1)/bcoeff(2)**2 / svol
!               r_im    = r_imfc / (1.+ delta) * ( delta**2 - 1. &
!                 & + (1.+0.5*kappa*(1.+ delta)**2) * ren/sqrtkap)           ! rim eq. 6 khl06
! 
!               ! number concentration and radius of ice particles
!               ni_hom  = phi / r_im                                         ! ni eq.9 khl06
!               ri_0    = 1. + 0.5 * sqrtkap * ren                           ! for eq. 3 khl06
!               ri_hom  = (ri_0 * (1. + delta) - 1. ) / bcoeff(2)            ! eq. 3 khl06 * ren = eq.23 khl06
!               mi_hom  = (4./3. * pi * rho_ice) * ni_hom * ri_hom**3
!               mi_hom  = max(dble(mi_hom),ice%x_min)
! 
!               nuc_n(i,j,k) = max(min(dble(ni_hom), ni_hom_max), 0.d0)
!               nuc_q(i,j,k) = min(nuc_n(i,j,k) * mi_hom, q(i,j,k))
! 
!               n_ice(i,j,k) = n_ice(i,j,k) + nuc_n(i,j,k)
!               q_ice(i,j,k) = q_ice(i,j,k) + nuc_q(i,j,k)
!               q(i,j,k)     = q(i,j,k)     - nuc_q(i,j,k)
! 
!             end if
!           end if
! 
!         enddo
!       enddo
!     enddo
!   end if
! 
!   deallocate(nuc_n, nuc_q)

  end subroutine ice_nucleation_homhet

  subroutine cloud_freeze (n1,rcloud,ninuc,rice,nice,tl,tk)
    integer, intent(in) :: n1
    real, intent(inout), dimension(n1) :: rice,nice,ninuc,rcloud,tl
    real, intent(in), dimension(n1) :: tk
    !     real, intent(in) :: nc
    !         ! .. local variables ..
    integer          ::  k
    real :: frr, frn,rc,xc,jhet,jhom,jtot,tc
    real, parameter :: ahet = 6.5e-1 ! 1/k,      after barklie and gokhale
    real, parameter :: bhet = 2.0e+2 ! 1/(m3 s), after barklie and gokhale
    real,save       :: facg

    facg = moment_gamma(cldw,2)    ! <hn
    do k = 1, n1
       if (tk(k) < tmelt .and. rcloud(k) > 0.) then


          rc = rcloud(k)
          tc = tk(k) - tmelt
          if (tc < -50.0) then
             frr = rc                                                    !..completely homogeneous freezing
             frn = cldw%nr                                              
          else
             xc = min(max(rc/(cldw%nr+eps0),cldw%x_min),cldw%x_max)    !..mean mass

             !..hom. freezing after jeffrey and austin (1997), see also cotton and field (2001)
             if (tc > -30.0) then
                jhom = 1.0d6 * exp(-7.63-2.996*(tc+30.0))           !..j in 1/(m3 s)
             else
                jhom = 1.0d6 * exp(-243.4-14.75*tc-0.307*tc**2-0.00287*tc**3-0.0000102*tc**4)
             endif

             !..het. freezing: stochastic model after bigg, numbers after barklie and gokhale
             jhet = bhet * ( exp( - ahet * tc) - 1.0 )            !..j in 1/(m3 s)
             !jhet = 0.0 ! neglected for cloud droplets


             jtot = (jhom + jhet) / rowt * dt                     !..j*dt in 1/kg
             frn  = jtot * rc
             frr  = frn * xc * facg
          end if
          frr  = min(frr,rc)
          !         if (frq>eps0)  frn  = cldw%nr*frq/rc
          rcloud(k) = rcloud(k) - frr
          tl(k) = tl(k)+ (convice(k)-convliq(k))*frr

          frn  = min(frn,frr/cldw%x_max)
!         ninuc(k) = ninuc(k)   - frn
          rice(k)   = rice(k)   + frr
          nice(k)   = nice(k)   + frn
! print *, 'cfr',k, tk(k), rcloud(k), rice(k),frr

       end if
    end do

  end subroutine cloud_freeze

  subroutine rain_freeze (n1,rrain,nrain,ninuc,rice,nice,rgrp,tk)
    integer, intent(in) :: n1
    real, dimension(n1), intent(inout) :: rrain,nrain,rice,nice,rgrp,ninuc
    real, dimension(n1), intent(in) :: tk
    ! .. local variables ..
    integer                     :: k
    real            :: fr_r,fr_n,r_r,x_r,n_r,j_het,&
         &  fr_r_i,fr_n_i,fr_r_g,fr_n_g,n_0,lam!,fr_r_tmp,fr_n_tmp
    real, save      :: coeff_z,xmax_ice
    real, parameter :: a_het = 6.5e-1 !  after barklie and gokhale (pk p.350)
    real, parameter :: b_het = 2.0e+2 !  after barklie and gokhale (pk p.350)
    real, parameter :: r_crit_fr = 1.000e-6 ! r-critical value for rain_freeze
    real, parameter :: d_rainfrz_ig = 0.50e-3 !  rain --> ice or graupel distinction
    coeff_z = moment_gamma(rain,2)

    xmax_ice = (d_rainfrz_ig/rain%a_geo)**(1.0e0/rain%b_geo)
    fr_r_g = 0.0
    do k = 1, n1
       r_r = rrain(k)
       n_r = nrain(k)

       if (tk(k) < tmelt) then
          if (r_r <= r_crit_fr) then
             if (tk(k) < t_hn) then
                fr_r = r_r                  !   below t_hn, everything freezes
                fr_n = n_r
                fr_n_i= n_r
                fr_r_i= r_r
                fr_n_g= 0.0
                fr_r_g= 0.0
             else
                fr_r = 0.0
                fr_n = 0.0
                fr_n_i= 0.0
                fr_r_i= 0.0
                fr_n_g= 0.0
                fr_r_g= 0.0
             end if
          else
             x_r = min(max(r_r/(n_r+eps0),rain%x_min),rain%x_max)
             n_r = r_r / x_r
             if (tk(k) < t_hn) then            !..only ice
                fr_r = r_r                  
                fr_n = n_r

                ! depending on their size, raindrops freeze to become either cloud ice or graupel

                lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                fr_n_i = n_0/(rain%mu*lam**((rain%nu+1.0)/rain%mu))* &
                     incgfct_lower((rain%nu+1.0)/rain%mu, lam*xmax_ice**rain%mu)
                fr_r_i = n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                     incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)

                fr_n_g = fr_n - fr_n_i
                fr_r_g = fr_r - fr_r_i

             else                           !..heterogenic freezing
                j_het = max(b_het * ( exp( a_het * (tmelt - tk(k))) - 1.0 ),0.) / rowt * dt
                ! depending on their size, raindrops freeze to become either cloud ice or graupel


                if (j_het >= 1-20) then
!               fr_n  = min(j_het * r_r,ninuc(k))
                   fr_r  = fr_n * x_r * coeff_z

                   lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
              n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                   fr_n_i = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                        incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                   fr_r_i = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                        incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_ice**rain%mu)
! print *,'a'   ,fr_r, r_r,ninuc(k)
                   fr_n = min(fr_n,n_r)
                   fr_r = min(fr_r,r_r)
                else
                   fr_n= 0.0
                   fr_r= 0.0
                   fr_n_i= 0.0
                   fr_r_i= 0.0
                   fr_r_g= 0.0
                end if

             end if

             fr_n = fr_n
             fr_r = fr_r

             fr_n_i = fr_n_i
             fr_r_i = fr_r_i
             fr_r_g = fr_r_g 
          end if
          rrain(k) = rrain(k) - fr_r
          nrain(k) = n_r - fr_n

          if (rrain(k) < 0.0) then
             write (*,*) 'rain_freeze: rrain < 0.0, ', k, tk(k), r_r, j_het, fr_r
             rrain(k) = 0.0e0
             stop
          end if
          if (nrain(k) < 0.0) then
             write (*,*) 'rain_freeze: nrain < 0.0, ', k, tk(k), n_r, j_het, fr_n
             nrain(k) = 0.0e0
             stop
          end if
          rice(k) = rice(k)  + fr_r_i
          nice(k) = nice(k)  + fr_n_i
          rgrp(k) = rgrp(k)  + fr_r_g
! if (r_r > r_crit_fr) print *, 'rfr',k, tk(k), rrain(k), rice(k),rgrp(k), j_het,fr_r

       end if
    end do
  end subroutine rain_freeze

  subroutine deposition(n1,meteor,ninuc,rice,nice,rv,tl,tk,s_i)
    use thrm, only : esi
    integer, intent(in) :: n1
    type(particle), intent(in) :: meteor
    real, dimension(n1), intent(inout) :: ninuc,rice,nice,rv,tl,s_i
    real, dimension(n1), intent(in)    :: tk

    integer                     :: k
    real        :: r_g,n_g,x_g,d_g,v_g,f_v,f_n,n_re,f_v_fakt,vent_fakt
    real        :: c_g                 !..coeff. for av. capacity
    real        :: a_f,b_f,a_n,b_n     !..coeff. for av. ventilation coef.
    real :: dep,ndep,gi
    c_g = 1.0 / meteor%cap
    a_n = vent_coeff_a(meteor,0)
    b_n = vent_coeff_b(meteor,0)
    a_f = vent_coeff_a(meteor,1)
    b_f = vent_coeff_b(meteor,1)

    f_v_fakt = n_sc**n_f
    vent_fakt = b_n / b_f
    do k = 1, n1

       ! hn: in case r_garupel=0, dep_graupel has to be zero too

       if (rice(k)> 0.0) then
          n_g = nice(k)                                     !..number density
          r_g = rice(k)                                     !..mass density

          x_g = min(max(r_g/(n_g+eps0),max(meteor%x_min,eps0)),meteor%x_max)  !..mean mass
          d_g = meteor%a_geo * exp(meteor%b_geo*log(x_g))          !..mean diameter
          v_g = meteor%a_vel * exp(meteor%b_vel*log(x_g)) * dn0(k)  !..mean sedimentation velocity

          n_re = v_g * d_g / nu_l                                    !..mean reynolds number
          f_v  = a_f + b_f * f_v_fakt * exp(m_f*log(n_re))           !..mean vent.coeff.
          f_n = a_n + vent_fakt * (f_v - a_f)                        !..mean vent.coeff.
          ! ub<<
          !         f_v  = max(f_v,1.e0) 
          !         f_n  = max(f_n,1.e0) 
  !       e_si = esi(tk(k))

          gi = 4.0*pi / ( alvi**2 / (K_T * Rm * tk(k)**2) + Rm * tk(k) / (D_v * esi(tk(k))) )
          ndep  = gi * n_g * c_g * d_g * s_i(k) * dt * f_n / x_g
        if (ndep < 0.) then
          ndep = max(max(ndep, -nice(k)),-rice(k) / x_g * f_n / f_v)
        end if
          dep   = ndep *x_g*f_v/f_n
          rice(k) = rice(k) + dep
          rv(k) = rv(k) - dep
          tl(k) = tl(k) + convice(k)*dep
        if (meteor%moments==2 .and. ndep < 0.) then
          nice(k) = nice(k) + ndep
!   	      ninuc(k) = ninuc(k) - ndep
          end if
       endif
    enddo

  end subroutine deposition


  subroutine melting(n1,meteor,r,thl,nr,rcld,rrain,nrain,tk)
    use thrm, only : esl
    integer, intent(in) :: n1
    real,dimension(n1),intent(in) :: tk
    type(particle), intent(in) :: meteor
    real,dimension(n1),intent(inout)  ::r,nr,rrain,nrain,rcld,thl

    integer                     :: k

    real            :: x_m,d_m,v_m,t_a,n_re,d_t,e_a
    real            :: melt,melt_v,melt_h,melt_n,melt_r
    real            :: fh_r,fv_r
    real            :: a_melt_n,b_melt_n
    real            :: a_melt_r,b_melt_r



    a_melt_n = vent_coeff_a(meteor,0)
    b_melt_n = vent_coeff_b(meteor,0)
    a_melt_r = vent_coeff_a(meteor,1)
    b_melt_r = vent_coeff_b(meteor,1)

    do k = 2,n1
       t_a = tk(k) !wrf!+ t(k) + t_g(k)
      e_a = esl(T_a)                                     !..Saturation pressure

       if (t_a > tmelt .and. r(k) > 0.0) then

          x_m = min(max(r(k)/(nr(k)+eps0),meteor%x_min),meteor%x_max)  !..mean re in si

          D_m = meteor%a_geo * x_m**meteor%b_geo                   !..mean diameter
          v_m = meteor%a_vel * x_m**meteor%b_vel * rho_0  !..mean fall velocity

          N_re = v_m * D_m / nu_l                             !..mean Reynolds nr
          fv_r = a_melt_r + b_melt_r * N_sc**n_f * N_re**m_f  !..mean Vent.coeff. water vapor
          !      fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mean Vent.coeff. water vapor

          D_T  = Kt / (cp * rho_0)!WRF!+rho_g(i,j,k)))
          fh_r = D_T / D_v * fv_r
          !       fh_n = D_T / D_v * fv_n

          melt   = 2.0*pi / (alvi-alvl) * D_m * nr(k) * dt

          melt_h = melt * Kt * (T_a - tmelt)
          melt_v = melt * D_v*alvl/Rm * (e_a/T_a - e_3/tmelt)

          melt_r = (melt_h * fh_r + melt_v * fv_r)
          !        melt_n = (melt_h * fh_n + melt_v * fv_n) / x_m

          melt_n = MIN(MAX( (melt_r - r(k)) / x_m + nr(k), 0.0e0), nr(k))

          melt_r = MIN(r(k),melt_r)
          melt_n = MIN(nr(k),melt_n)

          melt_r = MAX(0.e0,melt_r)
          melt_n = MAX(0.e0,melt_n)
          r(k) = r(k) - melt_r
          if (meteor%nr==2) then
             nr(k) = nr(k) - melt_n
          end if
          if(x_m<cldw%x_max) then
             rcld(k)    = rcld(k)    + melt_r
             thl(k)      = thl(k)      - convice(k)*melt_r
          else
             rrain(k)    = rrain(k)    + melt_r
             nrain(k)    = nrain(k)    + melt_n
          end if
       endif
    enddo
  end subroutine melting

  subroutine ice_selfcollection(n1,metin,metout,r_in,r_out,n_in,dn0,tk,d_conv,r_crit,d_crit)
    integer, intent(in) :: n1
    real, dimension(n1), intent(in) :: dn0,tk
    type(particle), intent(in) :: metin, metout
    real, dimension(n1), intent(inout) :: r_in,r_out,n_in
    real, intent(in) :: d_conv, r_crit, d_crit

    integer                     :: k
    real       :: r_i,n_i,x_i,d_i,v_i,e_coll
    real, save :: x_conv
    real       :: self_n,self_r
    real, dimension(3),save :: delta_n,delta_r, theta_n,theta_r
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
       delta_r(metnr) = 2.0*coll_delta_11(metin,1) + coll_delta_12(metin,metin,1)
       theta_n(metnr) = 2.0*coll_theta_11(metin,0) + coll_theta_12(metin,metin,0)
       theta_r(metnr) = 2.0*coll_theta_11(metin,1) + coll_theta_12(metin,metin,1)

       x_conv = (d_conv/metout%a_geo)**(1./metout%b_geo)

    end if

    do k = 1, n1
       r_i = r_in(k)                                   
       n_i = n_in(k)                                   
       x_i = min(max(r_i/(n_i+eps0),metin%x_min),metin%x_max)    
       d_i = metin%a_geo * x_i**metin%b_geo                     
       if ( n_i > 0.0 .and. r_i > r_crit .and. d_i > d_crit ) then
          if (tk(k)>tmelt) then
             e_coll = 1.0
          else
             !.. efficiency as a function of temperature, after cotton et al. (1986)
             !   (see also straka, 1989; s. 53)
             e_coll = min(10**(0.035*(tk(k)-tmelt)-0.7),0.2e0)
          end if
          v_i = metin%a_vel * x_i**metin%b_vel * dn0(k)      !..mean fall velocity
          self_r = pi * 0.25e0 * e_coll * delta_r(metnr) / dn0(k) * n_i * r_i * d_i * d_i &
               * ( theta_r(metnr) * v_i * v_i + 2.0*metin%s_vel**2 )**0.5 * dt
          self_r = min(self_r,r_i)

          r_in(k)  = r_in(k)  - self_r
          r_out(k) = r_out(k) + self_r

          self_n = pi * 0.25e0 * e_coll * delta_n(metnr) / dn0(k) * n_i * n_i * d_i * d_i &
               * ( theta_n(metnr) * v_i * v_i + 2.0*metin%s_vel**2 )**0.5 * dt
          self_n = min(min(self_n,self_r/x_conv),n_i)
          n_in(k)  = n_in(k)  - self_n
       endif
    enddo

  end subroutine ice_selfcollection

  subroutine ice_cloud_riming(n1,cloud,ice,r_c,r_i,n_i,r_g,thl,tk,d_coll,r_crit_c,d_crit_c,r_crit_i,d_crit_i,d_conv,e_ic)
    integer, intent(in) :: n1
    type(particle), intent(in) :: cloud,ice
    real, dimension(n1), intent(inout) :: r_c,r_i,n_i,r_g,thl
    real, dimension(n1), intent(in) :: tk
    real, intent(in) :: d_coll,r_crit_c,d_crit_c,r_crit_i,d_crit_i,e_ic
    real, intent(in), optional :: d_conv

    integer                     :: k
    real, parameter :: e_min = 0.01              !..min. eff. fuer gc,ic,sc
    !   real, parameter :: x_conv    = 0.100e-9 ! min. graupel-/hagelmass riming
    real     :: x_i,d_i,v_i
    real     :: x_c,d_c,v_c,e_coll,x_coll_c
    real     :: rime_r
    real     :: conv_n,conv_r
    real     :: mult_n,mult_r,mult_1,mult_2
    !     real     :: delta_n_ii,delta_n_ic,delta_n_cc
    real     :: delta_r_ii,delta_r_ic,delta_r_cc
    !     real     :: theta_n_ii,theta_n_ic,theta_n_cc
    real     :: theta_r_ii,theta_r_ic,theta_r_cc
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
       delta(1,1,metnr) = coll_delta_11(ice,1)
       delta(1,2,metnr) = coll_delta_12(ice,cloud,1)
       delta(2,2,metnr) = coll_delta_22(cloud,1)

       theta(1,1,metnr) = coll_theta_11(ice,1)
       theta(1,2,metnr) = coll_theta_12(ice,cloud,1)
       theta(2,2,metnr) = coll_theta_22(cloud,1)
    end if
    delta_r_ii =  delta(1,1,metnr)
    delta_r_ic =  delta(1,2,metnr)
    delta_r_cc =  delta(2,2,metnr)
    theta_r_ii =  theta(1,1,metnr)
    theta_r_ic =  theta(1,2,metnr)
    theta_r_cc =  theta(2,2,metnr)

    x_coll_c = (d_coll/cldw%a_geo)**3          !..minimal mass for collection, limits rime_n

    const1 = e_ic/(d_coll - d_crit_c)
    const2 = 1/x_coll_c
    const3 = 1/(t_mult_opt - t_mult_min)
    const4 = 1/(t_mult_opt - t_mult_max)
    const5 = alpha_spacefilling * rowt/roice
    do k = 1, n1
       x_c = min(max(r_c(k)/(cloud%nr+eps0),cloud%x_min),cloud%x_max) 
       d_c = cloud%a_geo * x_c**cloud%b_geo                   

       x_i = min(max(r_i(k)/(n_i(k)+eps0),ice%x_min),ice%x_max)      
       d_i = ice%a_geo * x_i**ice%b_geo                       


       if (r_c(k) > r_crit_c .and. r_i(k) > r_crit_i .and. d_i > d_crit_i .and. d_c > d_crit_c) then

          v_c = cloud%a_vel * x_c**cloud%b_vel * dn0(k)  !..mean fall velocity
          v_i = ice%a_vel   * x_i**ice%b_vel   * dn0(k) !..mean fall velocity

          e_coll = min(e_ic, max(const1*(d_c - d_crit_c), e_min))

          rime_r = pi/4.0 * e_coll / dn0(k) * n_i(k) * r_c(k) * dt &
               &   *     (delta_r_ii * d_i*d_i + delta_r_ic * d_i*d_c + delta_r_cc * d_c*d_c) &
               &   * sqrt(theta_r_ii * v_i*v_i - theta_r_ic * v_i*v_c + theta_r_cc * v_c*v_c  &
               &          +ice%s_vel**2)
          rime_r = min(r_c(k),rime_r)

          r_i(k)   = r_i(k)  + rime_r
          r_c(k) = r_c(k) - rime_r
          thl(k) = thl(k) + convice(k)*rime_r
          ! ice multiplication after hallet and mossop

          !         mult_r = 0.0
          if (tk(k) < tmelt .and. ice_multiplication .and. ice%moments == 2) then
             mult_1 = (tk(k) - t_mult_min)*const3
             mult_2 = (tk(k) - t_mult_max)*const4
             mult_1 = max(0.e0,min(mult_1,1.e0))
             mult_2 = max(0.e0,min(mult_2,1.e0))
             mult_n = c_mult * mult_1 * mult_2 * rime_r

             n_i(k)  = n_i(k)  + mult_n
          endif

          ! conversion ice -> graupel
          conv_r = 0.0
          conv_n = 0.0

          if(present(d_conv)) then
             if (d_i > d_conv) then
                conv_r = (rime_r - mult_r) / ( const5 * (pi/6.0*roice*d_i**3/x_i - 1.0) )
                ! d_i can't become smaller than d_conv_ig
                conv_r = min(r_i(k)-n_i(k)*(d_conv_ig/ice%a_geo)**(1.0/ice%b_geo),conv_r)
                conv_r = min(r_i(k),conv_r)
                ! ub >>
                x_i = min(max((r_i(k))/(n_i(k)+eps0),ice%x_min),ice%x_max)    !..mean mass incl. riming
                ! ub <<
                conv_n = conv_r / x_i
                conv_n = min(n_i(k),conv_n)
             end if
          endif

          r_i(k) = r_i(k) - conv_r
          r_g(k) = r_g(k) + conv_r

          if (ice%moments == 2) n_i(k)     = n_i(k)     - conv_n
       endif
    enddo

  end subroutine ice_cloud_riming

  subroutine ice_rain_riming(n1,rain,ice,r_r,n_r,r_i,n_i,r_g,dn0,tk,r_crit_r, d_crit_r,r_crit_i,d_crit_i)
    integer, intent(in) :: n1
    type(particle), intent(in) :: rain,ice
    real, dimension(n1), intent(inout) :: r_r,n_r,r_i,n_i,r_g
    real, dimension(n1), intent(in) :: dn0,tk
    real, intent(in) :: r_crit_r,d_crit_r,r_crit_i,d_crit_i

    integer                     :: k

    real            :: x_i,d_i,v_i,d_id
    real            :: x_r,d_r,v_r,d_rd
    real            :: rime_n,rime_ri,rime_rr
    real            :: mult_n,mult_r,mult_1,mult_2
    real      :: delta_n_ii,delta_n_ir,           delta_n_rr
    real      :: delta_r_ii,delta_r_ir,delta_r_ri,delta_r_rr
    real      :: theta_n_ii,theta_n_ir,           theta_n_rr
    real      :: theta_r_ii,theta_r_ir,theta_r_ri,theta_r_rr
    real, save :: d_av_fakt_r
    real,dimension(3),save :: d_av_fakt_i

    real, dimension(2,2,3),save :: delta_n,delta_r, theta_n,theta_r
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
       delta_n(1,1,metnr) = coll_delta_11(ice,0)
       delta_n(1,2,metnr) = coll_delta_12(ice,rain,0)
       delta_n(2,2,metnr) = coll_delta_22(rain,0)
       delta_r(1,1,metnr) = coll_delta_11(ice,1)
       delta_r(1,2,metnr) = coll_delta_12(ice,rain,1)
       delta_r(2,2,metnr) = coll_delta_22(rain,1)

       theta_n(1,1,metnr) = coll_theta_11(ice,0)
       theta_n(1,2,metnr) = coll_theta_12(ice,rain,0)
       theta_n(2,2,metnr) = coll_theta_22(rain,0)
       theta_r(1,1,metnr) = coll_theta_11(ice,1)
       theta_r(1,2,metnr) = coll_theta_12(ice,rain,1)
       theta_r(2,2,metnr) = coll_theta_22(rain,1)

       d_av_fakt_r = d_average_factor(rain)
       d_av_fakt_i(metnr) = d_average_factor(ice)
    end if
    delta_n_ii =  delta_n(1,1,metnr)
    delta_n_ir =  delta_n(1,2,metnr)
    delta_n_rr =  delta_n(2,2,metnr)
    theta_n_ii =  theta_n(1,1,metnr)
    theta_n_ir =  theta_n(1,2,metnr)
    theta_n_rr =  theta_n(2,2,metnr)
    delta_r_ii =  delta_r(1,1,metnr)
    delta_r_ir =  delta_r(1,2,metnr)
    delta_r_rr =  delta_r(2,2,metnr)
    theta_r_ii =  theta_r(1,1,metnr)
    theta_r_ir =  theta_r(1,2,metnr)
    theta_r_rr =  theta_r(2,2,metnr)

    do k = 1, n1

       x_r = min(max(r_r(k)/(n_r(k)+eps0),rain%x_min),rain%x_max)  
       d_r = rain%a_geo * x_r**rain%b_geo                   

       x_i = min(max(r_i(k)/(n_i(k)+eps0),ice%x_min),ice%x_max)    
       d_i = ice%a_geo * x_i**ice%b_geo                     
       d_id = d_av_fakt_i(metnr) * d_i

       if (r_r(k) > r_crit_r  .and. d_r > d_crit_r .and. r_i(k) > r_crit_i .and. d_i > d_crit_i) then
          x_r = min(max(r_r(k)/(n_r(k)+eps0),rain%x_min),rain%x_max)  
          v_i = ice%a_vel * x_i**ice%b_vel * dn0(k)    
          d_r = rain%a_geo * x_r**rain%b_geo           
          v_r = rain%a_vel * x_r**rain%b_vel * dn0(k)  
          d_rd = d_av_fakt_r * d_r

          rime_n  = pi/4.0 / dn0(k) * n_i(k) * n_r(k) * dt &
               &   * (delta_n_ii * d_i*d_i + delta_n_ir * d_i*d_r + delta_n_rr * d_r*d_r) &
               &   * sqrt(theta_n_ii * v_i*v_i - theta_n_ir * v_i*v_r + theta_n_rr * v_r*v_r  &
               &     +ice%s_vel**2)

          rime_rr = pi/4.0 / dn0(k) * n_i(k) * r_r(k) * dt &
               &   * (delta_n_ii * d_i*d_i + delta_r_ir * d_i*d_r + delta_r_rr * d_r*d_r) &
               &   * sqrt(theta_n_ii * v_i*v_i - theta_r_ir * v_i*v_r + theta_r_rr * v_r*v_r  &
               &     +ice%s_vel**2)

          rime_ri = pi/4.0 / dn0(k) * n_r(k) * r_i(k) * dt &
               &   * (delta_r_ii * d_i*d_i + delta_r_ri * d_i*d_r + delta_n_rr * d_r*d_r) &
               &   * sqrt(theta_r_ii * v_i*v_i - theta_r_ri * v_i*v_r + theta_n_rr * v_r*v_r  &
               &     +ice%s_vel**2)


          rime_n  = min(min(n_r(k),n_i(k)),rime_n)
          rime_rr = min(r_r(k),rime_rr)
          rime_ri = min(r_i(k),rime_ri)

          if (ice%moments==2) n_i(k)  = n_i(k)  - rime_n
          n_r(k) = n_r(k) - rime_n
          r_i(k)  = r_i(k)  - rime_ri
          r_r(k) = r_r(k) - rime_rr

          ! ice multiplication after hallet and mossop
          mult_r = 0.0
          if (tk(k) < tmelt .and. ice_multiplication .and. ice%moments==2) then
             mult_1 = (tk(k) - t_mult_min) / (t_mult_opt - t_mult_min)
             mult_2 = (tk(k) - t_mult_max) / (t_mult_opt - t_mult_max)
             mult_1 = max(0.e0,min(mult_1,1.e0))
             mult_2 = max(0.e0,min(mult_2,1.e0))
             mult_n = c_mult * mult_1 * mult_2 * rime_rr
             mult_r = mult_n * ice%x_min
             mult_r = min(rime_rr,mult_r)

             n_i(k) = n_i(k)  + mult_n
             r_i(k) = r_i(k)  + mult_r
          endif
          if (tk(k) >= tmelt) then
             if (d_id > d_rd) then
                ! rain is shedded from the ice
                if (ice%moments==2) n_i(k) = n_i(k) + rime_n
                n_r(k) = n_r(k) + rime_rr / x_r
                r_i(k) = r_i(k) + rime_ri - mult_r
                r_r(k) = r_r(k) + rime_rr
             else
                ! small ice particles melt immediately to rain
                n_r(k) = n_r(k) + rime_n
                r_r(k) = r_r(k) + rime_rr + rime_ri - mult_r
             end if
          else
             ! ice + the frozen water becomes graupel
             r_g(k) = r_g(k) + rime_ri + rime_rr - mult_r
          end if
       endif
    enddo

  end subroutine ice_rain_riming

  subroutine ice_collection(n1,met1,met2,r_i,n_i,r_g,tk,r_crit)
    integer, intent(in) :: n1
    type(particle), intent(in) :: met1,met2
    real, dimension(n1), intent(inout) :: r_i,r_g,n_i
    real, dimension(n1), intent(in) :: tk
    real, intent(in) :: r_crit

    integer                     :: k

    real            :: x_g,d_g,v_g
    real            :: x_i,d_i,v_i
    real      :: delta_n_ss,delta_n_si,delta_n_ii
    real      :: delta_r_ss,delta_r_si,delta_r_ii
    real      :: theta_n_ss,theta_n_si,theta_n_ii
    real      :: theta_r_ss,theta_r_si,theta_r_ii

    real            :: e_coll,coll_n,coll_r
    real, dimension(2,2,3,3),save :: delta_n,delta_r, theta_n,theta_r
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
       delta_n(1,1,metnr1,metnr2) = coll_delta_11(met2,0)
       delta_n(1,2,metnr1,metnr2) = coll_delta_12(met2,met1,0)
       delta_n(2,2,metnr1,metnr2) = coll_delta_22(met1,0)
       delta_r(1,1,metnr1,metnr2) = coll_delta_11(met2,1)
       delta_r(1,2,metnr1,metnr2) = coll_delta_12(met2,met1,1)
       delta_r(2,2,metnr1,metnr2) = coll_delta_22(met1,1)

       theta_n(1,1,metnr1,metnr2) = coll_theta_11(met2,0)
       theta_n(1,2,metnr1,metnr2) = coll_theta_12(met2,met1,0)
       theta_n(2,2,metnr1,metnr2) = coll_theta_22(met1,0)
       theta_r(1,1,metnr1,metnr2) = coll_theta_11(met2,1)
       theta_r(1,2,metnr1,metnr2) = coll_theta_12(met2,met1,1)
       theta_r(2,2,metnr1,metnr2) = coll_theta_22(met1,1)

    end if
    delta_n_ss =  delta_n(1,1,metnr1,metnr2)
    delta_n_si =  delta_n(1,2,metnr1,metnr2)
    delta_n_ii =  delta_n(2,2,metnr1,metnr2)
    theta_n_ss =  theta_n(1,1,metnr1,metnr2)
    theta_n_si =  theta_n(1,2,metnr1,metnr2)
    theta_n_ii =  theta_n(2,2,metnr1,metnr2)
    delta_r_ss =  delta_r(1,1,metnr1,metnr2)
    delta_r_si =  delta_r(1,2,metnr1,metnr2)
    delta_r_ii =  delta_r(2,2,metnr1,metnr2)
    theta_r_ss =  theta_r(1,1,metnr1,metnr2)
    theta_r_si =  theta_r(1,2,metnr1,metnr2)
    theta_r_ii =  theta_r(2,2,metnr1,metnr2)

    do k = 1, n1

       if (r_i(k) > r_crit .and. r_g(k) > r_crit) then

          if (tk(k) > tmelt) then
             e_coll = 1.0
          else
             !.. Temp. dependent sticking efficiency after lin (1983)
             e_coll = max(0.1e0,min(exp(0.09*(tk(k)-tmelt)),1.0e0))
          endif

          x_i = min(max(r_i(k)/(n_i(k)+eps0),met1%x_min),met1%x_max)     !..mean mass in si
          x_g = min(max(r_g(k)/(met2%nr+eps0),met2%x_min),met2%x_max)   !..mean mass in si

          d_i = met1%a_geo * x_i**met1%b_geo                      !..mean diameter
          v_i = met1%a_vel * x_i**met1%b_vel * dn0(k)     !..mean velocity

          d_g = met2%a_geo * x_g**met2%b_geo                    !..mean diameter
          v_g = met2%a_vel * x_g**met2%b_vel * dn0(k)   !..mean fall velocity


          coll_r = pi/4.0 * met2%nr * r_i(k) * e_coll * dt &
               &   * (delta_r_ss * d_g**2 + delta_r_si * d_g*d_i + delta_r_ii * d_i**2) &
               &   * (theta_r_ss * v_g**2 - theta_r_si * v_g*v_i + theta_r_ii * v_i**2  &
               &     +met2%s_vel**2 + met1%s_vel**2)**0.5

          coll_r = min(r_i(k),coll_r)

          r_g(k) = r_g(k) + coll_r
          r_i(k)  = r_i(k)  - coll_r
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

  subroutine sedimentation(n1,meteor,rp,nr,rrate,dtopt)
    integer, intent(in) :: n1
    real, dimension(n1), intent(inout) :: rp
    type(particle),  intent(in) :: meteor
    real, dimension(n1), intent(inout),optional :: nr
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
       vmin   = 1.0
       vlimit = 10.0
    case('hail')
       metnr = 4
       vmin   = 1.0
       vlimit = 10.0
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
       WRITE(0,'(a,i5,3a,e9.3)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  alfn = ',alfn(metnr)
       WRITE(0,'(a,i5,3a,e9.3)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  alfq = ',alfq(metnr)
       WRITE(0,'(a,i5,3a,e9.3)') 'sedim: myid = ',myid,' meteor = ',meteor%name,',  clam = ',c_lam(metnr)
       firsttime(metnr) =.false.
    end if

    do k=n1-1,2,-1

       Xp = rp(k) / (np(k)+eps0)
       lam = ( c_lam(metnr) * xp )**(meteor%b_vel)      
       vr(k) = alfq(metnr) * lam
       vr(k) = max(vr(k),vmin)
       vr(k) = min(vr(k),vlimit)
       !      vr(k) = -vr(k)

       vn(k) = alfn(metnr) * lam
       vn(k) = max(vn(k),vmin)
       vn(k) = min(vn(k),vlimit)
       !      vn(k) = -vn(k)
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

  real function n_ice_meyers_contact(t_a,s)
    ! diagnostic calculation of number of ice particles after meyers (1992)
    !
    real, intent(in):: s,t_a
    real, parameter :: n_0 = 1.0e+3
    real, parameter :: n_m = 1.0e+3
    real, parameter :: a_d = -0.639
    real, parameter :: b_d = 12.960
    real, parameter :: c_d = -2.8
    real, parameter :: d_d = 0.262
    real, parameter :: tmelt  = 2.732e+2     !..triple point water
    real :: temp
    temp = t_a - tmelt
    n_ice_meyers_contact = n_0 * exp( a_d + b_d * s )
    n_ice_meyers_contact = n_ice_meyers_contact  &
         &               + n_m * exp( c_d + d_d * temp )

  end function n_ice_meyers_contact

  real function vent_coeff_a(parti,n)

    integer        :: n
    type(particle) :: parti

    vent_coeff_a = parti%a_ven * gfct((parti%nu+n+parti%b_geo)/parti%mu)              &
         &                  / gfct((parti%nu+1.0)/parti%mu)                        &
         &                * ( gfct((parti%nu+1.0)/parti%mu)                        &
         &                  / gfct((parti%nu+2.0)/parti%mu) )**(parti%b_geo+n-1.0)

  end function vent_coeff_a

  real function vent_coeff_b(parti,n)

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
    !       gammafunction from numerical recipes (f77)                              *
    !                                                                              *
    !*******************************************************************************

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
  end function gfct

  real function incgfct_upper(a,x)
    !*******************************************************************************
    !
    ! upper incomplete gamma function
    !
    !              int(x)(oo) exp(-t) t^(a-1) dt
    !
    !*******************************************************************************

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

    !              int(0)(x) exp(-t) t^(a-1) dt
    !
    !   !*******************************************************************************

    real, intent(in) :: a, x
    real :: gam, gln

    gam = gammp(a,x,gln)
    incgfct_lower = exp(gln) * gam

  end function incgfct_lower


  real function incgfct(a,x1,x2)
    !*******************************************************************************
    !
    ! actual incomplete gamma-function
    !
    !*******************************************************************************

    real, intent(in) :: a, x1, x2

    incgfct = incgfct_lower(a,x2) - incgfct_lower(a,x1)

  end function incgfct

  real function gammp(a,x,gln)


    real, intent(in) :: a, x
    real, intent(out) :: gln

    real :: gammcf, gamser

    if (x.lt.0.0e0 .or. a .le. 0.0e0) then
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
  end function gammp

  real function gammq(a,x,gln)

    real, intent(in) :: a, x
    real, intent(out) :: gln
    real :: gammcf, gamser

    if (x.lt.0.0e0 .or. a .le. 0.0e0) then
       write(*,*) 'error in gammq: bad arguments (module gamma_functions)'
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
  end function gammq

  real function gammln(x)
    !*******************************************************************************
    !                                                                              *
    !       log(gamma function) taken from press et al.,  numerical recipes (f77)
    !
    !*******************************************************************************

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
  end function gammln

  real function gfct2(x)
    !*******************************************************************************
    !                                                                              *
    !       gamma function taken from press et al.,  numerical recipes (f77)
    !
    !       (other formulation, same results)
    !*******************************************************************************
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
  end function gfct2

  subroutine gcf(gammcf,a,x,gln)
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
       write (*,*) 'error in gcf: a too large, itmax too small (module gamma_functions)'
       gammcf = 0.0e0
       return
    end if

    gammcf=exp(-x+a*log(x)-gln)*h
  end subroutine gcf

  subroutine gser(gamser,a,x,gln)
    integer, parameter :: itmax = 100
    real, parameter :: eps=3.e-7
    real, intent(in) :: a, x
    real, intent(out) :: gamser, gln

    integer :: n
    real :: ap,del,sum

    gln=gammln(a)
    if (x.le.0.) then
       if (x.lt.0.) then
          write (*,*) 'error in gser: x < 0 (module gamma_functions)'
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
       write (*,*) '  (module gamma_functions)'
       gamser = 0.0e0
       return
    end if

    gamser = sum*exp(-x+a*log(x)-gln)

  end subroutine gser

  real function coll_delta(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_delta = gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)         &
         &                  / gfct((p1%nu+1.0  )/p1%mu)         &
         &        * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_geo+n)   &
         &        / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_geo+n)

  end function coll_delta

  real function coll_delta_11(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_delta_11 = coll_delta(p1,n)

  end function coll_delta_11

  real function coll_delta_22(p2,n)

    type(particle) :: p2
    integer        :: n

    coll_delta_22 = coll_delta(p2,n)

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

  end function coll_delta_12

  real function coll_theta(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_theta = gfct((2.0*p1%b_vel+2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &               / gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &          * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_vel)        &
         &          / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_vel)

  end function coll_theta

  real function coll_theta_11(p1,n)

    type(particle) :: p1
    integer        :: n

    coll_theta_11 = coll_theta(p1,n)

  end function coll_theta_11

  real function coll_theta_22(p2,n)

    type(particle) :: p2
    integer        :: n

    coll_theta_22 = coll_theta(p2,n)

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

  end function coll_theta_12

  function d_average_factor (parti)

    ! factor for the calculation of av. diameter for gamma-distributed hydrometeors
    ! valid for d = a_geo * x^b_geo

    real :: d_average_factor
    type(particle), intent(in) :: parti

    d_average_factor = &
         ( gfct( (parti%b_geo+parti%nu+1.0)/parti%mu ) / &
         gfct( (parti%nu+1.0)/parti%mu ) ) * &
         (gfct( (parti%nu+1.0)/parti%mu ) / gfct( (parti%nu+2.0)/parti%mu ) ) ** parti%b_geo

  end function d_average_factor

  real function e_es (t_)
    !*******************************************************************************
    !                        saturation pressure over ice                       *
    !*******************************************************************************

    real, intent (in) :: t_

    e_es  = e_3 * exp (a_e * (t_ - tmelt) / (t_ - b_e))

  end function e_es

  real function e_ws (t_)
    !*******************************************************************************
    !                      saturation pressure over water                      *
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
  end subroutine shuffle

  subroutine resetvar(meteor,mass,num)
    type(particle),intent(in)        :: meteor
    real, dimension(:), intent(inout) :: mass
    real, dimension(:), intent(inout), optional :: num

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
  end subroutine resetvar

  subroutine initmcrp(level,firsttime)
    integer , intent(in) :: level
    logical, optional, intent(inout) :: firsttime

    if (present(firsttime)) then
       firsttime = .false.
    end if

    if (level==2) then
       nprocess = 1
       microseq = isedimcd
    end if
    if (level>=3)      nprocess = nprocwarm
    if (level==4)      nprocess = nprocess + nprocice

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
       !Cloud ice properties
       ice = PARTICLE('ice', &
            2, & !number of moments
            0, & !Number of droplets
            0.000000, & !.nu.....Width parameter of the distribution
            0.333333, & !.mu.....exponential parameter of the distribution
            1.00e-07, & !.x_max..maximum particle mass D=???e-2m
            1.00e-12, & !.x_min..minimale particler mass D=200e-6m
            3.303633, & !.a_geo..coefficient of meteor geometry
            0.476191, & !.b_geo..coefficient of meteor geometry = 1/2.1
            2.77e+01, & !.a_vel..coefficient of fall velocity
            0.215790, & !.b_vel..coefficient of fall velocity = 0.41/1.9
            0.25    , & !.s_vel...dispersion of the fall velocity
            0.780000, & !.a_ven..ventilation coefficient (PK, p.541)
            0.308000, & !.b_ven..ventilation coefficient (PK, p.541)
            2.0)        !.cap....capacity coefficient
       
       !Snow properties
       
       snow =  PARTICLE('snow', & ! after Andy Heymsfield (CRYSTAL-FACE)
            1,  & !number of moments
            1e2, & !Number of droplets
            0.000000, & !.nu.....Width parameter of the distribution
            0.333333, & !.mu.....exponential parameter of the distribution
            1.00e-9, & !.x_max..maximum particle mass D=???e-2m
            1.00e-9, & !.x_min..minimale particler mass D=200e-6m
            3.303633, & !.a_geo..coefficient of meteor geometry
            0.476191, & !.b_geo..coefficient of meteor geometry = 1/2.1
            2.47e+02, & !.a_vel..coefficient of fall velocity
            0.333333, & !.b_vel..coefficient of fall velocity
            0.25    , & !.s_vel...dispersion of the fall velocity
            0.780000, & !.a_ven..ventilation coefficient (PK, p.541)
            0.308000, & !.b_ven..ventilation coefficient (PK, p.541)
            2.0)        !.cap....capacity coefficient
       
       !Graupel properties
       graupel = PARTICLE('graupel', & ! 'graupel'
            1,  & !number of moments
            1e2, &      !Number of droplets
            1.000000, & !.nu.....Width parameter of the distribution
            0.166666, & !.mu.....exponential parameter of the distribution
            1.00e-06, & !.x_max..maximum particle mass
            1.00e-06, & !.x_min..minimale particler mass
            1.10e-01, & !.a_geo..coefficient of meteor geometry
            0.300000, & !.b_geo..coefficient of meteor geometry = 1/3.10
            7.64e+01, & !.a_vel..coefficient of fall velocity
            0.255200, & !.b_vel..coefficient of fall velocity
            0.25    , & !.s_vel...dispersion of the fall velocity
            0.780000, & !.a_ven..ventilation coefficient (PK, p.541)
            0.308000, & !.b_ven..ventilation coefficient (PK, p.541)
            2.0)        !.cap....capacity coefficient

       hail = graupel
    end if

  end subroutine initmcrp

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
    use grid, only : dt,dzi_t,zm

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
             ssi(kk,jj,ii) =  R_d & 
                  & * dn0(kk) * qv(kk,jj,ii) &
                  & * tk(kk,jj,ii) &
                  & / e_es(dble(tk(kk,jj,ii))) - 1.0
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
           
          ! ... dynamics
          T_0(i,j,k)      = tk(kk,jj,ii)
          p_0(i,j,k)      = p00 * ((pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp)**cpr
          rho_k(i,j,k)    = dn0(kk)

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
             write (*,*) ' mcrph_sb: q_rain < 0, STOPPED AFTER CLOUDS 1'
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
          convice = alvi/cp*(pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp
          convliq = (alvl-alvi)/cp*(pi0(kk)+pi1(kk)+exner(kk,jj,ii))/cp
          tlt(kk,jj,ii) = tlt(kk,jj,ii) &
               &      - convice * (q_vap_new - q_vap_old) / dt  &
               &      + convliq * (q_liq_new - q_liq_old) / dt

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
         IF (ANY(qi(1:ke,j,i).gt.0.0)) call sedimentation (ke,pice, qi(1:ke,j,i),qni(1:ke,j,i),rrate = prec_i(1:ke,j,i))
         IF (ANY(qs(1:ke,j,i).gt.0.0)) call sedimentation (ke,psnow,qs(1:ke,j,i),qns(1:ke,j,i),rrate = prec_s(1:ke,j,i))
      end do
    end do
    ntsedi = 2
    DO j=3,je-2
      DO i=3,ie-2 
        IF (ANY(qg(1:ke,j,i).gt.0.0)) THEN
          DO ii=1,ntsedi
            call sedimentation (ke,pgraupel,qg(1:ke,j,i),qng(1:ke,j,i),rrate = prec_g(1:ke,j,i),dtopt=dt/ntsedi)
          end do
        ENDIF
      end do
    end do
    end if
    if (cloud_type.ge.2000) then
    prec_h = 0.0
    ntsedi = 2
    DO j=3,je-2
      DO i=3,ie-2 
        IF (ANY(qh(1:ke,j,i).gt.0.0)) THEN
          DO ii=1,ntsedi
            call sedimentation (ke,phail,qh(1:ke,j,i),qnh(1:ke,j,i),rrate = prec_g(1:ke,j,i),dtopt=dt/ntsedi)
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
