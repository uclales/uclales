! Module to define and initialize all land surface variables
! Adopted from DALES (vanHeerwarden)
!----------------------------------------------------------------------------
!

module lsmdata

  SAVE

  integer           :: nradtime  = 200

  ! Flags
  ! ----------------------------------------------------------
  logical          :: init_lsm   = .true.   !<   Flag for initializing LSM
  logical          :: local      = .true.   !<  Switch to MOST locally to get local Obukhov length
  logical          :: smoothflux = .false.  !<  Create uniform sensible & latent heat over domain
  logical          :: neutral    = .false.  !<  Disable stability corrections
  logical          :: filter     = .false.  !<  Filter variables at 3 dx to  prevent peak in
                                            !<  dimensionless wind profile if MOST-local enabled

  ! Soil properties
  ! ----------------------------------------------------------

  ! Domain-uniform properties
  integer,parameter :: ksoilmax = 4       !<  Number of soil layers [-]

  real              :: lambdasat          !<  heat conductivity saturated soil [W/m/K]
  real              :: Ke                 !<  Kersten number [-]

  real, allocatable :: zsoil  (:)         !<  Height of bottom soil layer from surface [m]
  real, allocatable :: zsoilc (:)         !<  Height of center soil layer from surface [m]
  real, allocatable :: dzsoil (:)         !<  Depth of soil layer [m]
  real, allocatable :: dzsoilh(:)         !<  Depth of soil layer between center of layers [m]

  ! Spatially varying properties
  real, allocatable :: lambda  (:,:,:)    !<  Heat conductivity soil layer [W/m/K]
  real, allocatable :: lambdah (:,:,:)    !<  Heat conductivity soil layer half levels [W/m/K]

  real, allocatable :: lambdas (:,:,:)    !<  Soil moisture diffusivity soil layer 
  real, allocatable :: lambdash(:,:,:)    !<  Soil moisture diffusivity soil half levels 

  real, allocatable :: gammas  (:,:,:)    !<  Soil moisture conductivity soil layer 
  real, allocatable :: gammash (:,:,:)    !<  Soil moisture conductivity soil half levels 

  real, allocatable :: rootf   (:,:,:)    !<  Root fraction per soil layer [-]
  real              :: rootfav (ksoilmax) !<  Average root fraction per soil layer [-]

  real, allocatable :: phiw    (:,:,:)    !<  Water content soil matrix [-]
  real              :: phiwav  (ksoilmax) !<  Average water content soil matrix [-]

  real, allocatable :: phiwm   (:,:,:)    !<  Water content soil matrix previous time step [-]
  real, allocatable :: phifrac (:,:,:)    !<  Relative water content per layer [-]
  real, allocatable :: phitot  (:,:)      !<  Total soil water content [-]
  real, allocatable :: pCs     (:,:,:)    !<  Volumetric heat capacity [J/m3/K]
  real, allocatable :: Dh      (:,:,:)    !<  Heat diffusivity

  real, allocatable :: sflxd_avn (:,:,:)  !<  Radiation of nradtime timesteps
  real, allocatable :: sflxu_avn (:,:,:)  !<  Radiation of nradtime timesteps
  real, allocatable :: lflxd_avn (:,:,:)  !<  Radiation of nradtime timesteps
  real, allocatable :: lflxu_avn (:,:,:)  !<  Radiation of nradtime timesteps

  real, allocatable :: tsoil   (:,:,:)    !<  Soil temperature [K]
  real, allocatable :: tsoilm  (:,:,:)    !<  Soil temperature previous time step [K]
  real              :: tsoilav (ksoilmax) !<  Average Soil temperature [K]

  real, allocatable :: tsoildeep (:,:)    !<  Deep soil temperature [K]
  real              :: tsoildeepav = 283. !< Average deep soil temperature [K]

  ! Soil related constants [adapted from ECMWF]
  real, parameter   :: phi       = 0.472  !<  volumetric soil porosity [-]
  real, parameter   :: phifc     = 0.323  !<  volumetric moisture at field capacity [-]
  real, parameter   :: phiwp     = 0.171  !<  volumetric moisture at wilting point [-]

  real, parameter   :: pCm       = 2.19e6 !<  Volumetric soil heat capacity [J/m3/K]
  real, parameter   :: pCw       = 4.2e6  !<  Volumetric water heat capacity [J/m3/K]

  real, parameter   :: lambdadry = 0.190  !<  Heat conductivity dry soil [W/m/K]
  real, parameter   :: lambdasm  = 3.11   !<  Heat conductivity soil matrix [W/m/K]
  real, parameter   :: lambdaw   = 0.57   !<  Heat conductivity water [W/m/K]

  real, parameter   :: bc        = 6.04     !< Clapp and Hornberger non-dimensional exponent [-]
  real, parameter   :: gammasat  = 0.57e-6  !< Hydraulic conductivity at saturation [m s-1]
  real, parameter   :: psisat    = -0.388   !< Matrix potential at saturation [m]

  ! Land surface properties
  ! ----------------------------------------------------------

  ! Surface properties
  real, allocatable :: z0m        (:,:) !<  Roughness length for momentum [m]
  real              :: z0mav    = 0.1

  real, allocatable :: z0h        (:,:) !<  Roughness length for heat [m]
  real              :: z0hav    = 0.025

  real, allocatable :: albedo     (:,:) !<  Surface albedo [-]
  real              :: albedoav = 0.25

  real, allocatable :: LAI        (:,:) !<  Leaf area index vegetation [-]
  real              :: LAIav    = 2.

  real, allocatable :: Cskin      (:,:) !<  Heat capacity skin layer [J]
  real              :: Cskinav  = 20000.

  real, allocatable :: cveg       (:,:) !<  Vegetation cover [-]
  real              :: cvegav   = 0.9

  real, allocatable :: lambdaskin (:,:) !<  Heat conductivity skin layer [W/m/K]
  real              :: lambdaskinav = 5.

  real, allocatable :: Wl         (:,:) !<  Liquid water reservoir [m]
  real              :: Wlav    = 0.     
                             
  real, parameter   :: Wmax    = 0.0002 !<  Maximum layer of liquid water on surface [m]
  real, allocatable :: Wlm        (:,:) !<  Liquid water reservoir previous timestep [m]

  real, allocatable :: cliq       (:,:) !<  Fraction of vegetated surface covered with liquid water 
  real, allocatable :: tskin      (:,:) !<  Skin temperature [K]
  real, allocatable :: tskinm     (:,:) !<  Skin temperature previous timestep [K]
  real, allocatable :: qskin      (:,:) !<  Skin specific humidity [kg/kg]
  real              :: tskinavg = 0.    !<  Slab average of tskin used by srfcsclrs 

  ! Surface energy balance
  real, allocatable :: Qnet     (:,:)   !<  Net radiation [W/m2]
  real              :: Qnetav   = 300.
  real, allocatable :: Qnetm    (:,:)   !<  Net radiation previous timestep [W/m2]
  real, allocatable :: Qnetn    (:,:)   !<  Net radiation dummy [W/m2]

  real, allocatable :: rsmin    (:,:)   !<  Minimum vegetation resistance [s/m]
  real              :: rsminav = 110.

  real, allocatable :: rssoilmin(:,:)   !<  Minimum soil evaporation resistance [s/m]
  real              :: rssoilminav = 50.

  real, allocatable :: gD       (:,:)   !<  Response factor vegetation to vapor pressure deficit [-]
  real              :: gDav	= 0.

  real, allocatable :: LE       (:,:)   !<  Latent heat flux [W/m2]
  real, allocatable :: H        (:,:)   !<  Sensible heat flux [W/m2]
  real, allocatable :: G0       (:,:)   !<  Ground heat flux [W/m2]
  real, allocatable :: ra       (:,:)   !<  Aerodynamic resistance [s/m]
  real, allocatable :: rs       (:,:)   !<  Composite resistance [s/m]
  real, allocatable :: rsveg    (:,:)   !<  Vegetation resistance [s/m]
  real, allocatable :: rsvegm   (:,:)   !<  Vegetation resistance previous timestep [s/m]
  real, allocatable :: rssoil   (:,:)   !<  Soil evaporation resistance [s/m]
  real, allocatable :: rssoilm  (:,:)   !<  Soil evaporation resistance previous timestep [s/m]
  real, allocatable :: tndskin  (:,:)   !<  Tendency of skin [W/m2]
  !real              :: rsisurf2 = 0.   !<  Vegetation resistance [s/m] if isurf2 is used

  ! Turbulent exchange variables
  real, allocatable :: obl     (:,:)    !<  local obuhkov length [m]
  real              :: oblav            !<  Spatially averaged obukhov length [m]
  real, allocatable :: cm      (:,:)    !<  Drag coefficient for momentum [-]
  real, allocatable :: cs      (:,:)    !<  Drag coefficient for scalars [-]

  real              :: u0av             !<  Mean u-wind component
  real              :: v0av             !<  Mean v-wind component
  real, allocatable :: u0bar    (:,:,:) !<  Filtered u-wind component
  real, allocatable :: v0bar    (:,:,:) !<  Filtered v-wind component
  real, allocatable :: thetaav (:,:)   !<  Filtered liquid water pot temp at first level
  real, allocatable :: vaporav (:,:)   !<  Filtered specific humidity at first full level
  real, allocatable :: tskinav (:,:)   !<  Filtered surface temperature
  real, allocatable :: qskinav (:,:)   !<  Filtered surface specific humidity

  contains

  !
  ! ----------------------------------------------------------
  ! Malte: initialize LSM and surface layer
  ! Adopted from DALES (vanHeerwaarden)
  !
  subroutine initlsm
  use grid, only : nzp, nxp, nyp, th00, vapor, iradtyp
  !use forc, only : sfc_albedo

    integer :: k,ierr

    ! --------------------------------------------------------
    ! Read LSM-specific NAMELIST (SURFNAMELIST)
    !
    namelist/SURFNAMELIST/ & 
    !< Switches
    local, filter, smoothflux, neutral, &
    !< Soil related variables
    tsoilav, tsoildeepav, phiwav, rootfav, &
    !< Land surface related variables
    z0mav, z0hav, rsisurf2, &
    Cskinav, lambdaskinav, albedoav, Qnetav, cvegav, Wlav, &
    !< Jarvis-Steward related variables
    rsminav, rssoilminav, LAIav, gDav

    open(17,file='SURFNAMELIST',status='old',iostat=ierr)
    read (17,SURFNAMELIST,iostat=ierr)
    if (ierr > 0) then
      print *, 'Problem in namoptions SURFNAMELIST'
      print *, 'iostat error: ', ierr
      stop 'ERROR: Problem in namoptions SURFNAMELIST'
    endif
    write(6 ,SURFNAMELIST)
    close(17)

    ! --------------------------------------------------------
    ! Allocate surface arrays
    ! 
    allocate(obl(nxp,nyp))
    allocate(ra(nxp,nyp))
    allocate(rs(nxp,nyp))
    allocate(z0m(nxp,nyp))
    allocate(z0h(nxp,nyp))
    allocate(tskin(nxp,nyp))
    allocate(qskin(nxp,nyp))
    allocate(cm(nxp,nyp))
    allocate(cs(nxp,nyp))
    allocate(albedo(nxp,nyp))

    if(iradtyp == 4) then
      allocate(sflxd_avn(nradtime,nxp,nyp))
      allocate(sflxu_avn(nradtime,nxp,nyp))
      allocate(lflxd_avn(nradtime,nxp,nyp))
      allocate(lflxu_avn(nradtime,nxp,nyp))
      
      sflxd_avn (:,:,:) = 0.
      sflxu_avn (:,:,:) = 0.
      lflxd_avn (:,:,:) = 0.
      lflxu_avn (:,:,:) = 0.
    end if

    ! Allocate LSM arrays
    allocate(zsoil(ksoilmax))
    allocate(zsoilc(ksoilmax))
    allocate(dzsoil(ksoilmax))
    allocate(dzsoilh(ksoilmax))
    allocate(lambda(nxp,nyp,ksoilmax))
    allocate(lambdah(nxp,nyp,ksoilmax))
    allocate(lambdas(nxp,nyp,ksoilmax))
    allocate(lambdash(nxp,nyp,ksoilmax))
    allocate(gammas(nxp,nyp,ksoilmax))
    allocate(gammash(nxp,nyp,ksoilmax))
    allocate(Dh(nxp,nyp,ksoilmax))
    allocate(phiw(nxp,nyp,ksoilmax))
    allocate(phiwm(nxp,nyp,ksoilmax))
    allocate(phifrac(nxp,nyp,ksoilmax))
    allocate(pCs(nxp,nyp,ksoilmax))
    allocate(rootf(nxp,nyp,ksoilmax))
    allocate(tsoil(nxp,nyp,ksoilmax))
    allocate(tsoilm(nxp,nyp,ksoilmax))
    allocate(tsoildeep(nxp,nyp))
    allocate(phitot(nxp,nyp))

    allocate(Qnet(nxp,nyp))
    allocate(Qnetm(nxp,nyp))
    allocate(Qnetn(nxp,nyp))
    allocate(LE(nxp,nyp))
    allocate(H(nxp,nyp))
    allocate(G0(nxp,nyp))

    allocate(rsveg(nxp,nyp))
    allocate(rsvegm(nxp,nyp))
    allocate(rsmin(nxp,nyp))
    allocate(rssoil(nxp,nyp))
    allocate(rssoilm(nxp,nyp))
    allocate(rssoilmin(nxp,nyp))
    allocate(cveg(nxp,nyp))
    allocate(cliq(nxp,nyp))
    allocate(tndskin(nxp,nyp))
    allocate(tskinm(nxp,nyp))
    allocate(Cskin(nxp,nyp))
    allocate(lambdaskin(nxp,nyp))
    allocate(LAI(nxp,nyp))
    allocate(gD(nxp,nyp))
    allocate(Wl(nxp,nyp))
    allocate(Wlm(nxp,nyp))

    ! Allocate filtered variables
    allocate(u0bar(nzp,nxp,nyp))
    allocate(v0bar(nzp,nxp,nyp))
    allocate(thetaav(nxp,nyp))
    allocate(vaporav(nxp,nyp))
    allocate(tskinav(nxp,nyp))
    allocate(qskinav(nxp,nyp))

    ! --------------------------------------------------------
    ! Initialize arrays
    ! 
    tskinm	= th00
    tskin	= th00
    qskin       = sum(vapor(2,3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)
    albedo	= albedoav
    z0m		= z0mav
    z0h		= z0hav
    ra		= 1.
    obl 	= -10e10
    oblav	= -10e10
    qskin	= sum(vapor(2,3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)

    dzsoil(1) = 0.07    ! z = 0.07 m  COSMO config from Linda
    dzsoil(2) = 0.27    ! z = 0.34 m
    dzsoil(3) = 1.13    ! z = 1.47 m
    dzsoil(4) = 1.39    ! z = 2.86 m 

    !dzsoil(1) = 0.07   ! ECMWF config from Chiel
    !dzsoil(2) = 0.21
    !dzsoil(3) = 0.72
    !dzsoil(4) = 1.89
 
    ! Calculate vertical layer properties
    zsoil(1)  = dzsoil(1)
    do k = 2, ksoilmax
      zsoil(k) = zsoil(k-1) + dzsoil(k)
    end do
    zsoilc = -(zsoil-0.5*dzsoil)
    do k = 1, ksoilmax-1
      dzsoilh(k) = 0.5 * (dzsoil(k+1) + dzsoil(k))
    end do
    dzsoilh(ksoilmax) = 0.5 * dzsoil(ksoilmax)

    ! Set evaporation related properties
    ! Set water content of soil - constant in this scheme
    phiw(:,:,1) = phiwav(1)
    phiw(:,:,2) = phiwav(2)
    phiw(:,:,3) = phiwav(3)
    phiw(:,:,4) = phiwav(4)
    
    phiwm(:,:,1) = phiwav(1)
    phiwm(:,:,2) = phiwav(2)
    phiwm(:,:,3) = phiwav(3)
    phiwm(:,:,4) = phiwav(4)    

    phitot = 0.0
    do k = 1, ksoilmax
      phitot(:,:) = phitot(:,:) + phiw(:,:,k) * dzsoil(k)
    end do
    phitot(:,:) = phitot(:,:) / zsoil(ksoilmax)
    do k = 1, ksoilmax
      phifrac(:,:,k) = phiw(:,:,k) * dzsoil(k) / zsoil(ksoilmax) / phitot(:,:)
    end do

    ! Set root fraction per layer for short grass
    rootf(:,:,1) = rootfav(1)
    rootf(:,:,2) = rootfav(2)
    rootf(:,:,3) = rootfav(3)
    rootf(:,:,4) = rootfav(4)

    tsoil(:,:,1)   = tsoilav(1)
    tsoil(:,:,2)   = tsoilav(2)
    tsoil(:,:,3)   = tsoilav(3)
    tsoil(:,:,4)   = tsoilav(4)
    tsoilm(:,:,1)   = tsoilav(1)
    tsoilm(:,:,2)   = tsoilav(2)
    tsoilm(:,:,3)   = tsoilav(3)
    tsoilm(:,:,4)   = tsoilav(4)
    tsoildeep(:,:) = tsoildeepav

    ! Calculate conductivity saturated soil
    lambdasat = lambdasm ** (1. - phi) * lambdaw ** (phi)

    Qnet       = Qnetav
    Qnetm      = Qnetav
    Qnetn      = Qnetav
    Cskin      = Cskinav
    lambdaskin = lambdaskinav
    rsmin      = rsminav
    rssoilmin  = rssoilminav
    LAI        = LAIav
    gD         = gDav
    cveg       = cvegav
    cliq       = 0.
    Wl         = Wlav
    !sfc_albedo = albedoav

  end subroutine initlsm

end module lsmdata
