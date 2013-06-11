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
module grid

!axel
!  use ncio, only : open_nc, define_nc
!irina
!  use step, only : case_name
  implicit none
  !
  integer           :: nxp = 132           ! number of x points
  integer           :: nyp = 132           ! number of y points
  integer           :: nzp = 105           ! number of z points

  logical           :: nxpart = .true.     ! number of processors in x

  real              :: deltax = 35.        ! dx for basic grid
  real              :: deltay = 35.        ! dy for basic grid
  real              :: deltaz = 17.5       ! dz for basic grid
  real              :: dzrat  = 1.02       ! grid stretching ratio
  real              :: dzmax  = 1200.      ! height to start grid-stretching
  real              :: dtlong = 10.0       ! long timestep
  real              :: th00   = 288.       ! basic state temperature

  real              :: CCN = 150.e6        ! Number of CCN per kg
  real              :: umean = 0.          ! Galilean transformation
  real              :: vmean = 0.          ! Galilean transformation
  real              :: sfc_albedo = 0.05   ! surface albedo

  integer           :: isfctyp = 0         ! surface flux parameterization type
  integer           :: igrdtyp = 1         ! vertical grid type
  integer           :: iradtyp = 0         ! radiation model type
  integer           :: level   = 0         ! thermodynamic level
  integer           :: naddsc  = 0         ! number of additional scalars

  logical           :: lmptend = .false.   ! Write out microphysical 3D fields.
  logical           :: lrad_ca = .false.   ! Perform clear air radiation calculations
  logical           :: lcouvreux = .false.  ! switch for 'radioactive' scalar
  logical           :: lwaterbudget = .false.  ! switch for liquid water budget diagnostics
  integer           :: ncvrx               ! Number of Couvreux scalar
  integer           :: ncld               ! Number of Couvreux scalar

  integer           :: nfpt = 10           ! number of rayleigh friction points
  real              :: distim = 300.0      ! dissipation timescale

  character (len=7), allocatable, save :: sanal(:)
  character (len=80):: expnme = 'Default' ! Experiment name
  character (len=80):: filprf = 'x'       ! File Prefix
  character (len=7) :: runtype = 'INITIAL'! Run Type Selection

  real, parameter   ::  rkalpha(3) = (/ 8./15., -17./60.,  3./4. /), &
                        rkbeta(3)  = (/    0.0,   5./12., -5./12./)

  integer           :: nz, nxyzp, nxyp, nstep
  real              :: dxi, dyi, dt, psrf
  real, dimension(:), allocatable :: xt, xm, yt, ym, zt, zm, dzi_t, dzi_m,        &
       u0, v0, pi0, pi1, th0, dn0, rt0, spngt, spngm, wsavex, wsavey,         &
       wfls, dthldtls,dqtdtls                                       !cgils
  !
  ! 2D Arrays (surface fluxes)
  !
  real, dimension (:,:), allocatable :: albedo, a_ustar, a_tstar, a_rstar,    &
       uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc, trac_sfc, sflxu_toa,sflxd_toa,lflxu_toa,lflxd_toa, sflxu_toa_ca,sflxd_toa_ca,lflxu_toa_ca,lflxd_toa_ca, &
       cnd_acc, &  ! accumulated condensation  [kg/m2] (diagnostic for 2D output)
       cev_acc, &  ! accumulated evaporation of cloud water [kg/m2] (diagnostic for 2D output)
       rev_acc     ! accumulated evaporation of rainwater   [kg/m2] (diagnostic for 2D output)
  integer, dimension(10) :: prc_lev = -1
  integer :: nv1, nv2, nsmp = 0

  !Malte: variables to restart the land surface
  real, dimension (:,:,:),  allocatable :: a_tsoil, a_phiw,                   &
                            a_sflxd_avn, a_sflxu_avn, a_lflxd_avn, a_lflxu_avn
  real, dimension (:,:),    allocatable :: a_tskin, a_qskin, a_Wl, a_Qnet, a_G0

  !Malte: variables to read homogeneous fluxes from nc file
  real, dimension (:),      allocatable :: shls, lhls, usls, timels

  !
  ! 3D Arrays
  !irina
  real, dimension (:,:,:), allocatable ::                                     &
       a_theta, a_pexnr, press, vapor, a_rflx, a_sflx, liquid, rsi,           &
       a_scr1, a_scr2, a_scr3, a_scr4, a_scr5, a_scr6, a_scr7,                &
       a_lflxu, a_lflxd, a_sflxu, a_sflxd,a_km, &
       prc_c, prc_r, prc_i, prc_s, prc_g, prc_h , prc_acc,               &
       a_lflxu_ca, a_lflxd_ca, a_sflxu_ca, a_sflxd_ca

  real, dimension (:,:), allocatable :: svctr
  real, dimension (:)  , allocatable :: ssclr
  !
  ! Named pointers (to 3D arrays)
  !
  real, dimension (:,:,:), pointer :: a_up, a_ut, a_vp, a_vt, a_wp, a_wt,     &
       a_sp, a_st, a_tp, a_tt, a_rp, a_rt, a_rpp, a_rpt, a_npp, a_npt,        &
       a_ricep , a_ricet  , & ! ice mixing ratio
       a_nicep , a_nicet  , & ! ice number concentration
       a_rsnowp, a_rsnowt , & ! snow
       a_nsnowp, a_nsnowt , &
       a_rgrp,   a_rgrt,    & ! graupel
       a_ngrp,   a_ngrt,    &
       a_rhailp, a_rhailt,  & ! hail
       a_nhailp, a_nhailt, a_rct, a_cld, a_cvrxp, a_cvrxt
 ! linda,b, output of tendencies
  real, dimension (:,:,:), allocatable :: &
        mp_qt, mp_qr, mp_qi, mp_qs, mp_qg, mp_qh, &
        mp_nqt, mp_nqr, mp_nqi, mp_nqs, mp_nqg, mp_nqh, &
        mp_tlt
! linda, e 

  character(40)      :: zname      = 'zt'
  character(40)      :: zhname     = 'zm'
  character(40)      :: zlongname  = 'Vertical position of cell centers'
  character(40)      :: zhlongname = 'Vertical position of cell faces'
  character(40)      :: zunit      = 'm'
  character(40)      :: xname      = 'xt'
  character(40)      :: xhname     = 'xm'
  character(40)      :: xlongname  = 'Longitudinal position of cell centers'
  character(40)      :: xhlongname = 'Longitudinal position of cell faces'
  character(40)      :: xunit      = 'm'
  character(40)      :: yname      = 'yt'
  character(40)      :: yhname     = 'ym'
  character(40)      :: ylongname  = 'Lateral position of cell centers'
  character(40)      :: yhlongname = 'Lateral position of cell faces'
  character(40)      :: yunit      = 'm'
  character(40)      :: tname      = 'time'
  character(40)      :: tlongname  = 'Time'
  character(40)      :: tunit      = 's'
  integer, parameter :: ictr = 0
  integer, parameter :: ihlf = 1
  !
  ! Memory for prognostic variables
  !
  real, dimension (:,:,:,:), allocatable, target :: a_xp, a_xt1, a_xt2
  !
  !
  integer :: nscl = 4
  !
contains
  !
  !----------------------------------------------------------------------
  !
  subroutine define_vars

    use mpi_interface, only :myid
    integer :: memsize

    allocate (u0(nzp),v0(nzp),pi0(nzp),pi1(nzp),th0(nzp),dn0(nzp),rt0(nzp))
!cgils
    allocate (wfls(nzp),dthldtls(nzp),dqtdtls(nzp))
    wfls(:) = 0.
    dthldtls(:) = 0.
    dqtdtls(:) = 0.

    memsize = 2*nxyzp ! complex array in pressure solver

    allocate (a_theta(nzp,nxp,nyp),a_pexnr(nzp,nxp,nyp),press(nzp,nxp,nyp),a_km(nzp,nxp,nyp))
    a_theta(:,:,:) = 0.
    a_pexnr(:,:,:) = 0.
    press(:,:,:) = 0.
    a_km(:,:,:) = 0.
    memsize = memsize + nxyzp*14 !

    if (level > 0) then
       allocate (vapor(nzp,nxp,nyp))
       vapor(:,:,:) = 0.
       memsize = memsize + nxyzp
       if (level > 1) then
          allocate (liquid(nzp,nxp,nyp))
          liquid(:,:,:) = 0.
          memsize = memsize + nxyzp
       end if
    end if

    if (iradtyp > 0 ) then
       allocate (a_rflx(nzp,nxp,nyp))
       a_rflx(:,:,:) = 0.
       memsize = memsize + nxyzp
       !irina
       allocate (a_lflxu(nzp,nxp,nyp))
       a_lflxu(:,:,:) = 0.
       memsize = memsize + nxyzp
       allocate (a_lflxd(nzp,nxp,nyp))
       a_lflxd(:,:,:) = 0.
       memsize = memsize + nxyzp
       if (lrad_ca) then
        allocate (a_lflxu_ca(nzp,nxp,nyp))
        a_lflxu_ca(:,:,:) = 0.
        memsize = memsize + nxyzp
        allocate (a_lflxd_ca(nzp,nxp,nyp))
        a_lflxd_ca(:,:,:) = 0.
        memsize = memsize + nxyzp
       end if
    end if
    !irina
    if (iradtyp > 2 .or. (iradtyp == 2 .and. level > 1)) then
       allocate (a_sflx(nzp,nxp,nyp),albedo(nxp,nyp),sflxu_toa(nxp,nyp),sflxd_toa(nxp,nyp),lflxu_toa(nxp,nyp),lflxd_toa(nxp,nyp))
       a_sflx(:,:,:) = 0.
       albedo(:,:) = 0.
       sflxu_toa(:,:) = 0.
       sflxd_toa(:,:) = 0.
       lflxu_toa(:,:) = 0.
       lflxd_toa(:,:) = 0.
       memsize = memsize + nxyzp + nxyp
       !irina
       allocate (a_sflxu(nzp,nxp,nyp))
       a_sflxu(:,:,:) = 0.
       memsize = memsize + nxyzp
       allocate (a_sflxd(nzp,nxp,nyp))
       a_sflxd(:,:,:) = 0.
       memsize = memsize + nxyzp
       if (lrad_ca) then 
        allocate (sflxu_toa_ca(nxp,nyp),sflxd_toa_ca(nxp,nyp),lflxu_toa_ca(nxp,nyp),lflxd_toa_ca(nxp,nyp))
        sflxu_toa_ca(:,:) = 0.
        sflxd_toa_ca(:,:) = 0.
        lflxu_toa_ca(:,:) = 0.
        lflxd_toa_ca(:,:) = 0.
        allocate (a_sflxu_ca(nzp,nxp,nyp))
        a_sflxu_ca(:,:,:) = 0.
        memsize = memsize + nxyzp
        allocate (a_sflxd_ca(nzp,nxp,nyp))
        a_sflxd_ca(:,:,:) = 0.
        memsize = memsize + nxyzp
      end if
    end if
 !
    allocate (a_scr1(nzp,nxp,nyp),a_scr2(nzp,nxp,nyp),a_scr3(nzp,nxp,nyp))
    allocate (a_scr4(nzp,nxp,nyp),a_scr5(nzp,nxp,nyp),a_scr6(nzp,nxp,nyp),a_scr7(nzp,nxp,nyp))
    a_scr1(:,:,:) = 0.
    a_scr2(:,:,:) = 0.
    a_scr3(:,:,:) = 0.
    a_scr4(:,:,:) = 0.
    a_scr5(:,:,:) = 0.
    a_scr6(:,:,:) = 0.
    a_scr7(:,:,:) = 0.
    memsize = memsize + 7*nxyzp

    nscl = nscl+naddsc
    if (level   > 0) nscl = nscl+1  ! qt only
    if (level   > 2) nscl = nscl+2  ! nr,qr
    if (level   > 3) nscl = nscl+4  ! ni,qi,qs,qg
    if (level   > 4) nscl = nscl+4  ! ns,ng,qh,nh (for Axel's two-moment scheme)

    if (lwaterbudget) then
      nscl = nscl+1 ! additional cloud water a_cld in the tracer array
      ncld = nscl
    end if
    if (lcouvreux) then
      nscl = nscl+1 ! Additional radioactive scalar
      ncvrx = nscl
    end if

    allocate (a_xp(nzp,nxp,nyp,nscl), a_xt1(nzp,nxp,nyp,nscl),        &
         a_xt2(nzp,nxp,nyp,nscl))

    a_xp(:,:,:,:)  = 0.
    a_xt1(:,:,:,:) = 0.
    a_xt2(:,:,:,:) = 0.

    a_up =>a_xp (:,:,:,1)
    a_vp =>a_xp (:,:,:,2)
    a_wp =>a_xp (:,:,:,3)
    a_tp =>a_xp (:,:,:,4)

    if (level > 0) a_rp =>a_xp (:,:,:,5)

    ! warm rain with number and mass of rain
    if (level >= 3) then
      a_rpp =>a_xp(:,:,:,6)
      a_npp =>a_xp(:,:,:,7)
      allocate (prc_acc(nxp,nyp,count( prc_lev>=0)))
      prc_acc(:,:,:) = 0.   ! accumulated precipitation for 2D output  [kg/m2]
      allocate (rev_acc(nxp,nyp))
      rev_acc(:,:) = 0.   ! accumulated evaporation of rain water    [kg/m2]
    else
      a_rpp => NULL()
      a_npp => NULL()
    end if
    if (lwaterbudget) then
      ! for liquid water budget and precipitation efficiency diagnostic
      a_cld=>a_xp(:,:,:,ncld)
      a_cld(:,:,:) = 0.
      allocate (cnd_acc(nxp,nyp),cev_acc(nxp,nyp))
      cnd_acc(:,:) = 0.   ! accumulated condensation                 [kg/m2]
      cev_acc(:,:) = 0.   ! accumulated evaporation of cloud water   [kg/m2]
    else
      a_cld => NULL()
    end if
    ! ice microphysics
    if (level >= 4) then
      a_ricep  =>a_xp(:,:,:, 8)
      a_nicep  =>a_xp(:,:,:, 9)
      a_rsnowp =>a_xp(:,:,:,10)
      a_rgrp   =>a_xp(:,:,:,11)
    else
      a_ricep  => NULL()
      a_nicep  => NULL()
      a_rsnowp => NULL()
      a_rgrp   => NULL()
    end if
    ! SB2006 two-moment scheme with hail
    if (level == 5) then
      a_nsnowp =>a_xp(:,:,:,12)
      a_ngrp   =>a_xp(:,:,:,13)
      a_rhailp =>a_xp(:,:,:,14)
      a_nhailp =>a_xp(:,:,:,15)
    else
      a_nsnowp => NULL()
      a_ngrp   => NULL()
      a_rhailp => NULL()
      a_nhailp => NULL()
    end if
    if (lcouvreux) then
      a_cvrxp=>a_xp(:,:,:,ncvrx)
      a_cvrxp(:,:,:) = 0.
      allocate (trac_sfc(nxp,nyp))
      trac_sfc = 1.
    else
      a_cvrxp => NULL()
    end if

    allocate (a_ustar(nxp,nyp),a_tstar(nxp,nyp))
    allocate (uw_sfc(nxp,nyp),vw_sfc(nxp,nyp),ww_sfc(nxp,nyp))
    allocate (wt_sfc(nxp,nyp))

    if (level > 0) allocate(wq_sfc(nxp,nyp),a_rstar(nxp,nyp))

    !Malte: allocate Land surface variables for restart
    if (isfctyp == 5) then
       allocate (a_tsoil(4,nxp,nyp))
       allocate (a_phiw (4,nxp,nyp))
       allocate (a_tskin (nxp,nyp))
       allocate (a_qskin (nxp,nyp))
       allocate (a_Wl    (nxp,nyp))
       allocate (a_Qnet  (nxp,nyp))
       allocate (a_G0    (nxp,nyp))
       allocate (a_sflxd_avn(100,nxp,nyp))
       allocate (a_sflxu_avn(100,nxp,nyp))
       allocate (a_lflxd_avn(100,nxp,nyp))
       allocate (a_lflxu_avn(100,nxp,nyp))
       memsize = memsize + 2*nxp*nyp*4 + 5*nxp*nyp + 4*nxp*nyp*100
    end if

    !Malte: allocate variables for homogeneous fluxes (no lsm used)
    if (isfctyp == 0) then
       allocate(shls(1740))
       allocate(lhls(1740))
       allocate(usls(1740))
       allocate(timels(1740))
       memsize = memsize + 4*1740
    end if
    !End Malte

    if (level >= 2) then
       allocate(prc_c(nzp,nxp,nyp))
       prc_c = 0.
       memsize = memsize + nxyzp
    end if
    if (level >= 3) then
       allocate(prc_r(nzp,nxp,nyp))
       prc_r = 0.
       memsize = memsize + nxyzp
    end if
    if (level >= 4) then
       allocate(prc_i(nzp,nxp,nyp))
       allocate(prc_s(nzp,nxp,nyp))
       allocate(prc_g(nzp,nxp,nyp))
       allocate(rsi(nzp,nxp,nyp))
       rsi = 0.
       prc_i = 0.
       prc_s = 0.
       prc_g = 0.
       memsize = memsize + 4*nxyzp
    end if
    if (level >= 5) then
       allocate(prc_h(nzp,nxp,nyp))
       prc_h = 0.
       memsize = memsize + nxyzp
    end if

    a_ustar(:,:) = 0.
    a_tstar(:,:) = 0.
    uw_sfc(:,:)  = 0.
    vw_sfc(:,:)  = 0.
    ww_sfc(:,:)  = 0.
    wt_sfc(:,:) = 0.
    
    if (level > 0) then
      wq_sfc(:,:) = 0.
      a_rstar(:,:) = 0.
    end if

    memsize = memsize +  nxyzp*nscl*2 + 3*nxyp + nxyp*10

!linda,b

    if (lmptend) then
       allocate(mp_qt(nzp,nxp,nyp))
       allocate(mp_qr(nzp,nxp,nyp))
       allocate(mp_qi(nzp,nxp,nyp))
       allocate(mp_qs(nzp,nxp,nyp))
       allocate(mp_qg(nzp,nxp,nyp))
       allocate(mp_qh(nzp,nxp,nyp))
       allocate(mp_nqt(nzp,nxp,nyp))
       allocate(mp_nqr(nzp,nxp,nyp))
       allocate(mp_nqi(nzp,nxp,nyp))
       allocate(mp_nqs(nzp,nxp,nyp))
       allocate(mp_nqg(nzp,nxp,nyp))
       allocate(mp_nqh(nzp,nxp,nyp))
       allocate(mp_tlt(nzp,nxp,nyp))

       memsize = memsize + 13*nxyzp

       mp_qt(:,:,:) = 0.
       mp_qr(:,:,:) = 0.
       mp_qi(:,:,:) = 0.
       mp_qs(:,:,:) = 0.
       mp_qh(:,:,:) = 0.
       mp_qg(:,:,:) = 0.
       mp_nqt(:,:,:) = 0.
       mp_nqr(:,:,:) = 0.
       mp_nqi(:,:,:) = 0.
       mp_nqs(:,:,:) = 0.
       mp_nqg(:,:,:) = 0.
       mp_nqh(:,:,:) = 0.
       mp_tlt(:,:,:) = 0.

    end if
!linda,e
    if(myid == 0) then
       print "(//' ',49('-')/,' ',/3x,i3.3,' prognostic scalars')", nscl
       print "('   memory to be allocated  -  ',f8.3,' mbytes')", &
            memsize*1.e-6*kind(0.0)
    end if

  end subroutine define_vars
  !
  !----------------------------------------------------------------------
  !
  subroutine define_grid

    use mpi_interface, only: xoffset, yoffset, wrxid, wryid, nxpg, nypg,   &
         appl_abort, myid

    integer :: i,j,k,kmax,nchby
    real    :: dzrfm,dz,zb,dzmin
    real    :: zmnvc(-1:nzp+1)
    character (len=51) :: &
         fm1 = '(//" ",49("-")/,"   grid dimensions:"/)            ',      &
         fm2 = '("   nxp-4 = ",i4,", dx, dx = ",f8.1,",",f9.1," m")',      &
         fm3 = '("   nyp-4 = ",i4,", dy, dy = ",f8.1,",",f9.1," m")',      &
         fm4 = '("   nzp   = ",i4,", dz, dz = ",f8.1,",",f9.1," m")',      &
         fm5 = '("   thermo level: ",i3)                        '

    nxyzp  = nxp*nyp*nzp
    nxyp   = nxp*nyp

    nz= nzp-1
    dzmin = 0.
    dxi=1./deltax
    dyi=1./deltay
    allocate(wsavex(4*nxpg+100),wsavey(4*nypg+100))
    wsavex=0.0
    wsavey=0.0

    !
    ! define xm array for grid 1 from deltax
    !
    allocate (xm(nxp))
    xm(1)=-float(max(nxpg-2,1))*.5*deltax+xoffset(wrxid)*deltax
    do i=2,nxp-1
       xm(i)=xm(i-1)+deltax
    end do
    xm(nxp)=2*xm(nxp-1)-xm(nxp-2)
    !
    ! define ym array for grid 1 from deltay
    !
    allocate (ym(nyp))
    ym(1)=-float(max(nypg-2,1))*.5*deltay+yoffset(wryid)*deltay
    do j=2,nyp-1
       ym(j)=ym(j-1)+deltay
    end do
    ym(nyp)=2*ym(nyp-1)-ym(nyp-2)

    !
    !      define where the momentum points will lie in vertical
    !
  allocate (zm(nzp))
  select case (abs(igrdtyp))
     !
     ! Read in grid spacings from a file
     !
  case(3)
     open (1,file='zm_grid_in',status='old',form='formatted')
     do k=1,nzp
        read (1,*) zm(k)
     end do
     close (1)
     if (zm(1) /= 0.) then
       if (myid == 0) print *, 'ABORTING:  Error in input grid'
       call appl_abort(0)
    end if
     !
     ! Tschebyschev Grid with vertical size given by dzmax
     !
  case(2)
     zm(1) = 0.
     nchby = nzp-3
     do k=1,nzp-2
        zm(k+1) = cos( ((2.*nchby - 1. - 2.*(k-1))*2.*asin(1.))/(2.*nchby))
        zm(k+1) = (zm(k+1)+1.)*dzmax/2.
     end do
     zm(nzp-1) = dzmax
     zm(nzp) = dzmax + (zm(nzp-1)-zm(nzp-2))
     !
     ! define zm array for grid 1 from deltaz and dzrat, if dzrat is
     ! negative compress grid so that dzmin is the grid spacing in a 100m
     ! interval below dzmax.  In both cases stretcvh grid uniformly by the
     ! ration |dzrat| above dzmax
     !
  case(1)
     zm(1)=0.
     zm(2)=deltaz
     zb=dzmax+100.
     if (dzrat.lt.0.) then
        dzmin = -float(int(dzrat))
        dzrat =  dzrat+dzmin-1
        kmax = int(log(deltaz/dzmin)/log(abs(dzrat)))
        zb=dzmax-100.
        do k=1,kmax
           zb=zb-dzmin*abs(dzrat)**k
        end do
     end if

     dz=deltaz
     do k=3,nzp
        if(zm(k-1) > zb .and. zm(k-1) < dzmax)then
           dz=max(dzmin,dz/abs(dzrat))
        else if (zm(k-1) >= dzmax) then
           dz=dz*abs(dzrat)
        end if
        zm(k)=zm(k-1)+dz
     end do
  case default
     zm(1)=0.
     do k=1,nzp
        zm(k)=zm(k-1)+deltaz
     end do
  end select
  !
  ! Grid Points for Thermal Points (T-Grid):
  !
  allocate (xt(nxp))
  do i=2,nxp
     xt(i)=.5*(xm(i)+xm(i-1))
  end do
  xt(1)=1.5*xm(1)-.5*xm(2)
  !
  allocate (yt(nyp))
  do j=2,nyp
     yt(j)=.5*(ym(j)+ym(j-1))
  end do
  yt(1)=1.5*ym(1)-.5*ym(2)
  !
  allocate (zt(nzp))
  if (igrdtyp .lt. 0) then
     !
     ! Read in grid spacings from a file
     !
     open (2,file='zt_grid_in',status='old',form='formatted')
     do k=1,nzp
        read (2,*) zt(k)
     end do
     close (2)
   else
     !
     ! calculate where the thermo points will lie based on geometric
     ! interpolation from the momentum points
     !
     do k=1,nzp
        zmnvc(k)=zm(k)
     end do
     zmnvc(0)=-(zmnvc(2)-zmnvc(1))**2 /(zmnvc(3)-zmnvc(2))
     zmnvc(-1)=zmnvc(0)-(zmnvc(1)-zmnvc(0))**2 /(zmnvc(2)-zmnvc(1))
     zmnvc(nzp+1)=zmnvc(nzp)+(zmnvc(nzp)-zmnvc(nzp-1))**2              &
                  /(zmnvc(nzp-1)-zmnvc(nzp-2))

     do k=1,nzp
       dzrfm=sqrt(sqrt((zmnvc(k+1)-zmnvc(k)) /(zmnvc(k-1)-zmnvc(k-2))))
       zt(k)=zmnvc(k-1)+(zmnvc(k)-zmnvc(k-1))/(1.+dzrfm)
     end do
  end if
  !
  ! compute other arrays based on the vertical grid.
  !   dzi_m: inverse of distance between thermal points k+1 and k
  !   dzi_t: inverse of distance between momentum points k and k-1
  !
  allocate (dzi_m(nzp))
  do k=1,nzp-1
     dzi_m(k)=1./(zt(k+1)-zt(k))
  end do
  dzi_m(nzp)=dzi_m(nzp-1)*dzi_m(nzp-1)/dzi_m(nzp-2)

  allocate (dzi_t(nzp))
  do k=2,nzp
     dzi_t(k)=1./(zm(k)-zm(k-1))
  end do
  dzi_t(1)=dzi_t(2)*dzi_t(2)/dzi_t(3)
  !
  ! set timesteps
  !
  if(myid == 0) then
     write(6,fm1)
     write(6,fm2) nxpg-4, deltax, 2.*xt(nxp-2)
     write(6,fm3) nypg-4, deltay, 2.*yt(nyp-2)

     write(6,fm4) nzp,zm(2)-zm(1),zm(nzp)
     write(6,fm5) level
  endif

  end subroutine define_grid
  !
  ! ----------------------------------------------------------------------
  ! Subroutine write_hist:  This subroutine writes a binary history file
  !
  subroutine write_hist(htype, time)

    use mpi_interface, only : appl_abort, myid, wrxid, wryid

    integer :: errcode=-17

    integer, intent (in) :: htype
    real, intent (in)    :: time

    character (len=80) :: hname

    integer :: n, iblank, nseed
    integer, allocatable, dimension(:) :: seed
    !
    ! create and open a new output file.
    !
    write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
    hname = trim(hname)//'.'//trim(filprf)

    select case(htype)
    case default
       hname = trim(hname)//'.iflg'
    case(0)
       hname = trim(hname)//'.R'
    case(1)
       hname = trim(hname)//'.rst'
    case(2)
       iblank=index(hname,' ')
       write (hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(time), 's'
    end select

    call random_seed(size=nseed)
    allocate (seed(nseed))
    call random_seed(get=seed)

    !
    ! Write fields
    !
    if (myid == 0) print "(//' ',49('-')/,' ',/,'   History write to: ',A30)" &
         ,hname
    open(10,file=trim(hname), form='unformatted')

    write(10) time,th00,umean,vmean,dt,level,iradtyp,nzp,nxp,nyp,nscl
    write(10) nseed
    write(10) seed
    write(10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf
    write(10) a_ustar, a_tstar, a_rstar
    write(10) a_pexnr

    !Malte: Restart land surface
    if (isfctyp == 5) then
       write(10) a_tsoil
       write(10) a_phiw
       write(10) a_tskin
       write(10) a_qskin
       write(10) a_Wl
       write(10) a_sflxd
       write(10) a_sflxu
       write(10) a_lflxd
       write(10) a_lflxu
       write(10) a_sflxd_avn
       write(10) a_sflxu_avn
       write(10) a_lflxd_avn
       write(10) a_lflxu_avn
    end if
    !End Malte

    do n=1,nscl
       call newvar(n)
       write(10) a_sp
    end do
    if(level>=3) then
      write(10) prc_acc, rev_acc
    end if
    if(lwaterbudget) then
      write(10) cnd_acc, cev_acc
    end if
    write(10) nv2, nsmp
    write(10) svctr
    close(10)

    if (myid == 0 .and. htype < 0) then
       print *, 'CFL Violation'
       call appl_abort(errcode)
    end if

    return
  end subroutine write_hist
  !
  ! ----------------------------------------------------------------------
  ! Subroutine read_hist:  This subroutine reads a binary history file
  !
  subroutine read_hist(time, hfilin)

    use mpi_interface, only : appl_abort, myid, wrxid, wryid

    character(len=80), intent(in) :: hfilin
    real, intent(out)             :: time

    character (len=80) :: hname
    integer :: n, nseed, nxpx, nypx, nzpx, nsclx, iradx, lvlx
    integer, dimension(:), allocatable :: seed
    logical :: exans
    real :: umx, vmx, thx
    !
    ! open input file.
    !

    write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
    hname = trim(hname)//'.'//trim(hfilin)

    inquire(file=trim(hname),exist=exans)
    if (.not.exans) then
       print *,'ABORTING: History file', trim(hname),' not found'
       call appl_abort(0)
    else
       open (10,file=trim(hname),status='old',form='unformatted')
       read (10) time,thx,umx,vmx,dt,lvlx,iradx,nzpx,nxpx,nypx,nsclx
       read (10) nseed
       allocate(seed(nseed))
       read(10) seed
       call random_seed(put=seed)
       if (nxpx /= nxp .or. nypx /= nyp .or. nzpx /= nzp)  then
          if (myid == 0) print *, nxp, nyp, nzp, nxpx, nypx, nzpx
          call appl_abort(-1)
       end if

       read (10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf
       read (10) a_ustar, a_tstar, a_rstar
       read (10) a_pexnr

       !Malte: Restart land surface
       if (isfctyp == 5) then
          read(10) a_tsoil
          read(10) a_phiw
          read(10) a_tskin
          read(10) a_qskin
          read(10) a_Wl
          read(10) a_sflxd
          read(10) a_sflxu
          read(10) a_lflxd
          read(10) a_lflxu
          read(10) a_sflxd_avn
          read(10) a_sflxu_avn
          read(10) a_lflxd_avn
          read(10) a_lflxu_avn
       end if
       !End Malte

       do n=1,nscl
          call newvar(n)
          if (n <= nsclx) read (10) a_sp
       end do
       do n=nscl+1,nsclx
          read (10)
       end do
      if(level>=3) then
        read(10) prc_acc, rev_acc
      end if
      if(lwaterbudget) then
        read(10) cnd_acc, cev_acc
      end if
      read(10) nv2, nsmp
      allocate (svctr(nzp,nv2))
      read(10) svctr

       close(10)
       !
       ! adjust namelist and basic state appropriately
       !
       if (thx /= th00) then
          if (myid == 0) print "('  th00 changed  -  ',2f8.2)",th00,thx
          a_tp(:,:,:) = a_tp(:,:,:) + thx - th00
       end if
       if (umx /= umean) then
          if (myid == 0) print "('  umean changed  -  ',2f8.2)",umean,umx
          a_up = a_up + umx - umean
          u0 = u0 + umx - umean
       end if
       if (vmx /= vmean) then
          if (myid == 0) print "('  vmean changed  -  ',2f8.2)",vmean,vmx
          a_vp = a_vp + vmx - vmean
          v0 = v0 +vmx - vmean
       end if

    end if

  end subroutine read_hist
  !
  ! ----------------------------------------------------------------------
  ! Subroutine newvar:  This routine updates the variabe pointer to the
  ! value corresponding to the next variable in the table
  !
  subroutine newvar(inum,istep)

    integer, intent(in) :: inum
    integer, optional, intent(in)   :: istep
    a_sp =>a_xp (:,:,:,inum)

    if (present(istep)) then
       select case (istep)
       case(1)
          a_st=>a_xt1(:,:,:,inum)
       case(2)
          a_st=>a_xt2(:,:,:,inum)
       case(3)
          a_st=>a_xt1(:,:,:,inum)
       end select
    end if

  end subroutine newvar

end module grid

