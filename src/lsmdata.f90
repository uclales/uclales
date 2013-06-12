! Module to define and initialize all land surface variables
! Adopted from DALES (vanHeerwarden)
!
! Output is written in both ts and analysis files.
! For correct 3D output level=5 is required (see ncio.f90).
!
! Malte Rieck, June 2012
!
!----------------------------------------------------------------------------
!

module lsmdata

  SAVE

  integer           :: nradtime  = 100

  ! Flags
  ! ----------------------------------------------------------
  logical          :: init_lsm   = .true.   !<  Flag for initializing LSM
  logical          :: local      = .true.   !<  Switch to MOST locally to get local Obukhov length
  logical          :: smoothflux = .false.  !<  Create uniform sensible & latent heat over domain
  logical          :: neutral    = .false.  !<  Disable stability corrections
  logical          :: hetero     = .false.  !<  Switch to heteogeneous surface conditions
  logical          :: filter     = .false.  !<  Filter variables at 3 dx to  prevent peak in
                                            !<  dimensionless wind profile if MOST-local enabled

  ! Soil properties
  ! ----------------------------------------------------------

  ! Domain-uniform properties
  integer,parameter :: ksoilmax = 4       !<  Number of soil layers [-]

  ! Surface heterogeneity variables
  integer           :: hetper = 1         !<  Set number of heterogeneity periods in the domain
  integer           :: hetlen             !<  Length of a patch [number of grid cells]
  integer           :: xhet               !<  Loop variable for patch in x-dir
  integer           :: yhet               !<  Loop variable for patch in y-dir

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

 !real, allocatable :: phiw    (:,:,:)    !<  Water content soil matrix [-]
  real              :: phiwav  (ksoilmax) !<  Average water content soil matrix [-]

  real, allocatable :: phiwm   (:,:,:)    !<  Water content soil matrix previous time step [-]
  real, allocatable :: phifrac (:,:,:)    !<  Relative water content per layer [-]
  real, allocatable :: phitot  (:,:)      !<  Total soil water content [-]
  real, allocatable :: pCs     (:,:,:)    !<  Volumetric heat capacity [J/m3/K]
  real, allocatable :: Dh      (:,:,:)    !<  Heat diffusivity

 !real, allocatable :: sflxd_avn (:,:,:)  !<  Radiation of nradtime timesteps
 !real, allocatable :: sflxu_avn (:,:,:)  !<  Radiation of nradtime timesteps
 !real, allocatable :: lflxd_avn (:,:,:)  !<  Radiation of nradtime timesteps
 !real, allocatable :: lflxu_avn (:,:,:)  !<  Radiation of nradtime timesteps

 !real, allocatable :: tsoil   (:,:,:)    !<  Soil temperature [K]
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

  !real, allocatable :: albedo    (:,:) !<  Surface albedo [-]

  real, allocatable :: LAI        (:,:) !<  Leaf area index vegetation [-]
  real, allocatable :: LAIG       (:,:) !<  Global leaf area index vegetation [-]
  real              :: LAIav    = 4.    !<  Average leaf area index [-]
  real              :: LAImin   = 2.    !<  Minimum leaf area index [-]  
  real              :: LAImax   = 6.    !<  Maximum leaf area index [-]

  real, allocatable :: Cskin      (:,:) !<  Heat capacity skin layer [J]
  real              :: Cskinav  = 20000.

  real, allocatable :: cveg       (:,:) !<  Vegetation cover [-]
  real              :: cvegav   = 0.9

  real, allocatable :: lambdaskin (:,:) !<  Heat conductivity skin layer [W/m/K]
  real              :: lambdaskinav = 5.

 !real, allocatable :: Wl         (:,:) !<  Liquid water reservoir [m]
  real              :: Wlav     = 0.     
                             
  real, parameter   :: Wmax     = 0.0002 !<  Maximum layer of liquid water on surface [m]
  real, allocatable :: Wlm        (:,:) !<  Liquid water reservoir previous timestep [m]

  real, allocatable :: cliq       (:,:) !<  Fraction of vegetated surface covered with liquid water 
 !real, allocatable :: tskin      (:,:) !<  Skin temperature [K]
  real, allocatable :: tskinm     (:,:) !<  Skin temperature previous timestep [K]
 !real, allocatable :: qskin      (:,:) !<  Skin specific humidity [kg/kg]
  real              :: tskinavg = 0.    !<  Slab average of tskin used by srfcsclrs 

  ! Surface energy balance
  !real, allocatable :: Qnet     (:,:)   !<  Net radiation [W/m2]
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
  !real, allocatable :: G0       (:,:)   !<  Ground heat flux [W/m2]
  real, allocatable :: ra       (:,:)   !<  Aerodynamic resistance [s/m]
  real, allocatable :: rsurf    (:,:)   !<  Composite resistance [s/m]
  real, allocatable :: rsveg    (:,:)   !<  Vegetation resistance [s/m]
  real, allocatable :: rsvegm   (:,:)   !<  Vegetation resistance previous timestep [s/m]
  real, allocatable :: rssoil   (:,:)   !<  Soil evaporation resistance [s/m]
  real, allocatable :: rssoilm  (:,:)   !<  Soil evaporation resistance previous timestep [s/m]
  real, allocatable :: tndskin  (:,:)   !<  Tendency of skin [W/m2]

  ! Turbulent exchange variables
  !real, allocatable :: obl     (:,:)    !<  local obuhkov length [m] BvS: moved to grid, used beyond lsm
  real              :: oblav            !<  Spatially averaged obukhov length [m]
  real, allocatable :: cm      (:,:)    !<  Drag coefficient for momentum [-]
  real, allocatable :: cs      (:,:)    !<  Drag coefficient for scalars [-]

  real              :: u0av             !<  Mean u-wind component
  real              :: v0av             !<  Mean v-wind component
  real, allocatable :: u0bar    (:,:,:) !<  Filtered u-wind component
  real, allocatable :: v0bar    (:,:,:) !<  Filtered v-wind component
  real, allocatable :: thetaav  (:,:)   !<  Filtered liquid water pot temp at first level
  real, allocatable :: vaporav  (:,:)   !<  Filtered specific humidity at first full level
  real, allocatable :: tskinav  (:,:)   !<  Filtered surface temperature
  real, allocatable :: qskinav  (:,:)   !<  Filtered surface specific humidity

  contains

  !
  ! ----------------------------------------------------------
  ! Malte: initialize LSM and surface layer
  ! Adopted from DALES (vanHeerwaarden)
  !
  subroutine initlsm(sst,time_in)
    
    use grid, only : nzp, nxp, nyp, th00, vapor, iradtyp, a_G0, a_Qnet, &
                     a_tskin, a_qskin, a_phiw, a_tsoil, a_Wl, dt

    use mpi_interface, only: myid, xoffset, yoffset, wrxid, wryid, nxpg, nypg

    integer :: k,ierr

    real, intent(in) :: sst,time_in

    ! --------------------------------------------------------
    ! Read LSM-specific NAMELIST (SURFNAMELIST)
    !
    namelist/SURFNAMELIST/ & 
    !< Switches
    local, filter, smoothflux, neutral, hetero, &
    !< Soil related variables
    tsoilav, tsoildeepav, phiwav, rootfav, &
    !< Land surface related variables
    z0mav, z0hav, Cskinav, albedoav, &
    lambdaskinav, Qnetav, cvegav, Wlav, &
    !< Jarvis-Steward related variables
    rsminav, rssoilminav, LAIav, gDav, &
    !< Heterogeneity related variables
    hetper, LAImin, LAImax

    open(17,file='SURFNAMELIST',status='old',iostat=ierr)
    read (17,SURFNAMELIST,iostat=ierr)
    if (ierr > 0) then
      print *, 'Problem in namoptions SURFNAMELIST'
      print *, 'iostat error: ', ierr
      stop 'ERROR: Problem in namoptions SURFNAMELIST'
    endif
    !write(6 ,SURFNAMELIST)
    close(17)

    ! --------------------------------------------------------
    ! Allocate land surface arrays
    ! 

    !Allocate surface scheme arrays
    !allocate(obl(nxp,nyp))
    allocate(ra(nxp,nyp))
    allocate(rsurf(nxp,nyp))
    allocate(z0m(nxp,nyp))
    allocate(z0h(nxp,nyp))
    allocate(cm(nxp,nyp))
    allocate(cs(nxp,nyp))
    !allocate(albedo(nxp,nyp))

    ! Allocate LSM arrays
    allocate(zsoil(ksoilmax))
    allocate(zsoilc(ksoilmax))
    allocate(dzsoil(ksoilmax))
    allocate(dzsoilh(ksoilmax))
    allocate(lambda(ksoilmax,nxp,nyp))
    allocate(lambdah(ksoilmax,nxp,nyp))
    allocate(lambdas(ksoilmax,nxp,nyp))
    allocate(lambdash(ksoilmax,nxp,nyp))
    allocate(gammas(ksoilmax,nxp,nyp))
    allocate(gammash(ksoilmax,nxp,nyp))
    allocate(Dh(ksoilmax,nxp,nyp))
    allocate(phitot(nxp,nyp))
    allocate(phiwm(ksoilmax,nxp,nyp))
    allocate(phifrac(ksoilmax,nxp,nyp))
    allocate(pCs(ksoilmax,nxp,nyp))
    allocate(rootf(ksoilmax,nxp,nyp))
    allocate(tsoilm(ksoilmax,nxp,nyp))
    allocate(tsoildeep(nxp,nyp))

    !allocate(Qnet(nxp,nyp))
    allocate(Qnetm(nxp,nyp))
    allocate(Qnetn(nxp,nyp))
    allocate(LE(nxp,nyp))
    allocate(H(nxp,nyp))
    !allocate(G0(nxp,nyp))

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
    allocate(Wlm(nxp,nyp))

    ! Allocate global variables (nxpg,nypg)
    allocate(LAIG(nxpg,nypg))

    ! Allocate filtered variables
    allocate(u0bar(nzp,nxp,nyp))
    allocate(v0bar(nzp,nxp,nyp))
    allocate(thetaav(nxp,nyp))
    allocate(vaporav(nxp,nyp))
    allocate(tskinav(nxp,nyp))
    allocate(qskinav(nxp,nyp))

    !Variables are now defined in grid
    !allocate(Wl(nxp,nyp))
    !allocate(tskin(nxp,nyp))
    !allocate(qskin(nxp,nyp))
    !allocate(phiw(ksoilmax,nxp,nyp))
    !allocate(tsoil(ksoilmax,nxp,nyp))

    ! --------------------------------------------------------
    ! Initialize arrays
    ! 

    if ((time_in*86400.) .le. dt) then

       a_tskin          = sst
       a_qskin          = sum(vapor(2,3:nxp-2,3:nyp-2))/(nxp-4)/(nyp-4)

       a_phiw(1,:,:)    = phiwav(1)
       a_phiw(2,:,:)    = phiwav(2)
       a_phiw(3,:,:)    = phiwav(3)
       a_phiw(4,:,:)    = phiwav(4)

       a_tsoil(1,:,:)   = tsoilav(1)
       a_tsoil(2,:,:)   = tsoilav(2)
       a_tsoil(3,:,:)   = tsoilav(3)
       a_tsoil(4,:,:)   = tsoilav(4)

       a_Wl             = Wlav
    end if 

    z0m       = z0mav
    z0h       = z0hav
    
    dzsoil(1) = 0.07    !z = 0.07 m COSMO config from Linda
    dzsoil(2) = 0.27    !z = 0.34 m
    dzsoil(3) = 1.13    !z = 1.47 m
    dzsoil(4) = 1.39    !z = 2.86 m 

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
    phitot = 0.0
    do k = 1, ksoilmax
       phitot(:,:) = phitot(:,:) + a_phiw(k,:,:) * dzsoil(k)
    end do
    phitot(:,:)    = phitot(:,:) / zsoil(ksoilmax)

    do k = 1, ksoilmax
      phifrac(k,:,:) = a_phiw(k,:,:) * dzsoil(k) / zsoil(ksoilmax) / phitot(:,:)
    end do

    ! Set root fraction per layer for short grass (not used for now!!!)
    rootf(1,:,:) = rootfav(1)
    rootf(2,:,:) = rootfav(2)
    rootf(3,:,:) = rootfav(3)
    rootf(4,:,:) = rootfav(4)

    tsoildeep(:,:) = tsoildeepav

    ! Calculate conductivity saturated soil
    lambdasat  = lambdasm ** (1. - phi) * lambdaw ** (phi)

    Cskin      = Cskinav
    lambdaskin = lambdaskinav
    rsmin      = rsminav
    rssoilmin  = rssoilminav
    LAI        = LAIav
    gD         = gDav
    cveg       = cvegav

    ! --------------------------------------------------------http://journals.ametsoc.org/doi/abs/10.1175/2011MWR3599.1
    ! Static Heterogeneity of Vegetation
    ! 
    ! std:    LAI=4. albedo=0.20 z0mav=0.035 z0hav=0.035
    ! grass:  LAI=2. albedo=0.25 z0mav=0.035 z0hav=0.035
    ! forest: LAI=6. albedo=0.15 z0mav=0.500 z0hav=0.500
    ! hetper sets the number of heterogeneity periods in the domain
    ! hetper*2 is the number of patches in x or y direction
    
    if (hetero) then
      hetlen = nxpg/(2*hetper)

      do yhet=0,(hetper*2)
         do xhet=0,(hetper*2) 

            if (mod(yhet,2) .eq. 0) then 
               if (mod(xhet,2) .eq. 0) then
                   LAIG(xhet*hetlen+3:xhet*hetlen+hetlen+3,yhet*hetlen+3:yhet*hetlen+hetlen+3) = LAImin
               else if (mod(xhet,2) .ne. 0) then
                   LAIG(xhet*hetlen+3:xhet*hetlen+hetlen+3,yhet*hetlen+3:yhet*hetlen+hetlen+3) = LAImax
               end if
            end if              

            if (mod(yhet,2) .ne. 0) then
               if (mod(xhet,2) .eq. 0) then
                   LAIG(xhet*hetlen+3:xhet*hetlen+hetlen+3,yhet*hetlen+3:yhet*hetlen+hetlen+3) = LAImax
               else if (mod(xhet,2) .ne. 0) then
                   LAIG(xhet*hetlen+3:xhet*hetlen+hetlen+3,yhet*hetlen+3:yhet*hetlen+hetlen+3) = LAImin
               end if 
            end if   

         end do
      end do

      LAI(3:(nxp-2),3:(nyp-2)) = LAIG(3+xoffset(wrxid):nxp+xoffset(wrxid)-2, &
                                  3+yoffset(wryid):nyp+yoffset(wryid)-2 )

    end if

    deallocate(LAIG)

  end subroutine initlsm

end module lsmdata
