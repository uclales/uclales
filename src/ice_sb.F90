!===============================================================================!
! Ulrich Blahak, 18.3.2009:
!
! Newly vectorized routines (the originals are still left in the code for comparisons and tests):
!
! - ice_nucleation_vec         (untested)
!   - n_ice_huffmann_vali_vec  (untested)
!
! - complete_ice_snow_rim_vec  (untested)
!
! - rain_sedi_lm_vec           (untested)
! - ice_sedi_lm_vec            (untested)
! - snow_sedi_lm_vec           (untested)
! - graupel_sedi_lm_vec        (untested)
! - hail_sedi_lm_vec           (untested)
!
! - rain_freeze_gamlook   (spectral partitioning of freezing rain; inc_gfct replaced by equidistant LUT)
! - vapor_dep_simple      (the complicated limiting of evaporation is replaced by something very simple)
! - clnuc_sk_4D           (Segal/Khain cloud activation; transformed into an equidistant LUT)
! - graupel_hail_conv_wet_gamlook (calls to inc_gfct replaced by equidistant LUT)
! - dmin_wg_gr_ltab_equi  (LUT for wet growth diameter of graupel made equidistant)
!
!===============================================================================!
!
!===============================================================================!
! Ulrich Blahak, 18.3.2009:
!
! Further changes and Bugfixes (partly in german, sorry ...):
!
!!! - in Segal/Khain cloud activation: linear extrapolation for NCN from 50 down to 0 m**-3!
!!!
!!! - in init_dmin_graupel_wetgrowth_lookup(): Fehler in der Einlesereihenfolge von
!!!   anzp_wg_g, anzT_wg_g,anzw_wg_g,anzi_wg_g korrigiert.
!!!
!!! - von speichere_umwandlungsraten auf speichere_dqdt umgestellt. Nur verwenden
!!!   in Verbindung mit src_gscp_wetgrowth_20081202.F90.isobar.incloudnuc
!!!
!!! - bei snow_melting/graupel_melting/hail_melting falsches fh_q: Nach Auswertung der
!!!   entsprechenden Formeln fuer fh_q und fv_q aus Rasmussen and Heymsfield (1987, Part I)
!!!   fuer einen weiten T- und p-Bereich ergibt sich fh_q ungefaehr gleich 1.05 * fv_q, was
!!!   im Code jetzt beruecksichtigt wird. Die alte Formel ergab dagegen
!!!   eine weite Streuung des Verh. fh_q / fv_q von unter 0.5 bis ueber 1.5, je nach
!!!   T und p, was an einer falschen Definition (fh_q ist nicht D_T/D_v*fv_q!) und 
!!!   an der anschliessenden Vernachlaessigung der Temperaturabhaengigkeit von D_T und D_v
!!!   liegt. Mit Temperaturabh. waere nach der falschen Definition fh_q ca. 0.83*fv_q.
!!!
!!! - ice_rain_riming / snow_rain_riming:  bei rime_qi, rime_qr z.T. falsche deltas und thetas verbessert
!!!   Ausserdem dort die Speicherung von Umwandlungsraten (Prozessindices) korrigiert.
!!!   Das ganze musste auch in complete_ice_snow_riming und complete_ice_snow_rim_vec angepasst werden.
!!!
!!! - graupel_rain_riming:  bei rime_q, rime_qr ebenfalls 
!!!   1) deltas und thetas in der ersten IF-Verzweigung konsistent gemacht
!!!      (diese waren nicht falsch)
!!!   und 2) in der zweiten IF-Verzweigung Fehler korrigiert.
!!!
!!! - BUG in iceCRY2test behoben! Dort waren versehentlich noch die ageo- und bgeo-Werte 
!!!   von iceCRY drin, welche viel zu hohe Bulk-Dichten bei kleinen Partikeln produziert haben.
!!!   Nach ersten Erfahrungen fallen dadurch im Endeffekt die Reflektivitaeten in 
!!!   den Aufwindschlaeuchen in der Hoehe etwas kleiner aus, d.h. kleinere und/oder mehr 
!!!   feste Teilchen.
!
!
!
!===============================================================================!
!
!
!===============================================================================!
!
! Two-moment bulk microphysics by Axel Seifert
!
! (with contributions from Uli Blahak and Heike Noppel)
!
! Description:
! Provides various modules and subroutines for two-moment bulk microphysics
!
! Current Code Owner: Axel Seifert, DWD
!                     axel.seifert@dwd.de
!
! Language: Fortran 90.
!
! Some code standards or recommendations, at least:
!
! - All new variables/subroutines should be named in English
! - In the future also comments should be written in English, 
!   but temporary use of some German is OK, too.
! - Length of names of subroutines should be <= 20
! - Length of names of variables should be <= 15
! - Length of lines has to be <= 100 including comments,
!   recommended is <= 80 for code without comments. 
! - Temporary modifications for experiments should be marked by
!
!     AS_YYYYMMDD>
!         ! Change for / Bugfix ....
!     <AS_YYYYMMDD
!
!   until they are validated as improvments and made official
!   (with AS, or whatever, being the initials of the guy doing the stuff
!   and YYYYMMDD=year/month/day).
!
!===============================================================================!
!
! This module can be used for LM or WRF. If it is used for LM, set
! the following preprocessor switch to 1, otherwise to 0:
! (either uncomment the following line or do it with the appropriate compiler options)
!

#define FOR_LM    0

!!!===============================================================================================
!!!
!!! ATTENTION: For the LM, we need 3D rainrate variables (rain_r, ...) in the interface routine!!!
!!!            (For WRF, these are still 2D)
!!!
!!!===============================================================================================


!!! Vorsicht: hier satad_nach_mikrophysik = .true. !!!
!!!
!!! Temperaturlimit bei cloud_nucleation_SK herausgenommen (nach Erkenntnissen von Heike)!!!
!!!
!!! Bei snow_rain_riming und ice_rain_riming kann bei T_a > T_3 nun kein Graupel mehr entstehen.
!!! Stattdessen werden alle Kollisionsteilchen in Regen umgewandelt oder bleiben unveraendert.


!*******************************************************************************!
!
! If used with LM:
!   use with interface module
!     src_gscp_wetgrowth_20080212.F90.isobar.incloudnuc
!
! If used with WRF:
!   use with interface 
!     module module_mp_seifert_wrfiface_20080212.F
!
!*******************************************************************************!
!
! ub: Now exp. decrease of N_ccn with height for all kinds of cloud nucleation types!
!

MODULE wolken_konstanten

  IMPLICIT NONE 

  !INTERFACE OPERATOR(**)
  !   module procedure power_explog_dd     
  !   module procedure power_explog_dr     
  !END INTERFACE

  TYPE PARTICLE
    CHARACTER(20)    :: name  !..Bezeichnung der Partikelklasse
    DOUBLE PRECISION :: nu    !..Breiteparameter der Verteil.
    DOUBLE PRECISION :: mu    !..Exp.-parameter der Verteil.
    DOUBLE PRECISION :: x_max !..maximale Teilchenmasse
    DOUBLE PRECISION :: x_min !..minimale Teilchenmasse
    DOUBLE PRECISION :: a_geo !..Koeff. Geometrie
    DOUBLE PRECISION :: b_geo !..Koeff. Geometrie = 1/3
    DOUBLE PRECISION :: a_vel !..Koeff. Fallgesetz
    DOUBLE PRECISION :: b_vel !..Koeff. Fallgesetz
    DOUBLE PRECISION :: a_ven !..Koeff. Ventilationsparam.
    DOUBLE PRECISION :: b_ven !..Koeff. Ventilationsparam.
    DOUBLE PRECISION :: cap   !..Koeff. Kapazitaet
  END TYPE PARTICLE

! UB_20090227>>
  ! Following are two new type declarations to hold the values
  ! of equidistand lookup tables. Up to now, such lookup tables are
  ! used in the Segal-Khain parameterization of CCN-activation and
  ! for determining the wet growth diameter of graupel.

  ! Type declaration for a general 2D equidistant lookup table:
  TYPE lookupt_2D
    INTEGER :: n1  ! number of grid points in x1-direction
    INTEGER :: n2  ! number of grid points in x2-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x1  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x2  ! grid vector in x1-direction
    DOUBLE PRECISION                     :: dx1          ! dx1   (grid distance w.r.t. x1)
    DOUBLE PRECISION                     :: dx2          ! dx2   (grid distance w.r.t. x2)
    DOUBLE PRECISION                     :: odx1         ! one over dx 1
    DOUBLE PRECISION                     :: odx2         ! one over dx 2
    DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: ltable
  END TYPE lookupt_2D
  
  ! Type declaration for a general 4D equidistant lookup table:
  TYPE lookupt_4D
    INTEGER :: n1  ! number of grid points in x1-direction
    INTEGER :: n2  ! number of grid points in x2-direction
    INTEGER :: n3  ! number of grid points in x3-direction
    INTEGER :: n4  ! number of grid points in x4-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x1  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x2  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x3  ! grid vector in x1-direction
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x4  ! grid vector in x1-direction
    DOUBLE PRECISION                     :: dx1          ! dx1   (grid distance w.r.t. x1)
    DOUBLE PRECISION                     :: dx2          ! dx2   (grid distance w.r.t. x2)
    DOUBLE PRECISION                     :: dx3          ! dx3   (grid distance w.r.t. x3)
    DOUBLE PRECISION                     :: dx4          ! dx4   (grid distance w.r.t. x4)
    DOUBLE PRECISION                     :: odx1         ! one over dx 1
    DOUBLE PRECISION                     :: odx2         ! one over dx 2
    DOUBLE PRECISION                     :: odx3         ! one over dx 3
    DOUBLE PRECISION                     :: odx4         ! one over dx 4
    DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: ltable
  END TYPE lookupt_4D
! UB_20090227<<

  ! .. allgemeine physikalische Konstanten und Parameter ..

  DOUBLE PRECISION, PARAMETER :: nu_l = 1.460d-5     !..Kinem. Visc. von Luft
  DOUBLE PRECISION, PARAMETER :: D_v  = 3.000d-5     !..Diffusivitaet von Wasserdampf
  DOUBLE PRECISION, PARAMETER :: K_T  = 2.500d-2     !..Waermeleitfaehigkeit     
  DOUBLE PRECISION, PARAMETER :: L_ew = 3.330d+5     !..Schmelzwaerme     
  DOUBLE PRECISION, PARAMETER :: L_ed = 2.830d+6     !..Sublimationswaerme     
  DOUBLE PRECISION, PARAMETER :: L_wd = 2.500d+6     !..Verdampfungswaerme     
  DOUBLE PRECISION, PARAMETER :: K_w  = 0.930d+0     !..Dielektrizitaetsfaktor Wasser
  DOUBLE PRECISION, PARAMETER :: K_i  = 0.176d+0     !..Dielektrizitaetsfaktor Eis
  DOUBLE PRECISION, PARAMETER :: rho0 = 1.225d+0     !..Norm-Luftdichte
  DOUBLE PRECISION, PARAMETER :: T_3  = 2.732d+2     !..Tripelpunkt Wasser
  DOUBLE PRECISION, PARAMETER :: T_f  = 2.330d+2     !..Bei T < T_f kein Fl.wasser
  DOUBLE PRECISION, PARAMETER :: R_l  = 2.870d+2     !..Gaskonstante trockener Luft 
  DOUBLE PRECISION, PARAMETER :: R_d  = 4.615d+2     !..Gaskonstante von Wasserdampf
  DOUBLE PRECISION, PARAMETER :: c_w  = 4.218d+3     !..spezifische Waerme von Wasser
  DOUBLE PRECISION, PARAMETER :: sigma_wa  = 7.1d-2  !..Oberflaechenspannung Wasser-Luft
  DOUBLE PRECISION, PARAMETER :: rho_w   = 1000.0    !..Materialdichte von Fluessigwasser
  DOUBLE PRECISION, PARAMETER :: rho_ice = 900.0     !..Materialdichte von Eis

  ! .. wolkenphysikalische Konstanten und Parameter ..

  DOUBLE PRECISION, PARAMETER :: N_sc = 0.710        !..Schmidt-Zahl (PK, S.541)
  DOUBLE PRECISION, PARAMETER :: n_f  = 0.333        !..Exponent von N_sc im Vent-koeff. (PK, S.541)
  DOUBLE PRECISION, PARAMETER :: m_f  = 0.500        !..Exponent von N_re im Vent-koeff. (PK, S.541)

  DOUBLE PRECISION, PARAMETER :: A_e  = 2.18745584d1 !..Konst. Saettigungsdamppfdruck - Eis
  DOUBLE PRECISION, PARAMETER :: A_w  = 1.72693882d1 !..Konst. Saettigungsdamppfdruck - Wasser
  DOUBLE PRECISION, PARAMETER :: B_e  = 7.66000000d0 !..Konst. Saettigungsdamppfdruck - Eis
  DOUBLE PRECISION, PARAMETER :: B_w  = 3.58600000d1 !..Konst. Saettigungsdamppfdruck - Wasser
  DOUBLE PRECISION, PARAMETER :: e_3  = 6.10780000d2 !..Saettigungsdamppfdruck bei T = T_3

  DOUBLE PRECISION, PARAMETER :: C_mult     = 3.5d8    !..Koeff. fuer Splintering
  DOUBLE PRECISION, PARAMETER :: T_mult_min = 265.0    !..Minimale Temp. Splintering
  DOUBLE PRECISION, PARAMETER :: T_mult_max = 270.0    !..Maximale Temp. Splintering
  DOUBLE PRECISION, PARAMETER :: T_mult_opt = 268.0    !..Optimale Temp. Splintering

  ! ... spezielle Parameter des KAMM2-Wolkenmoduls

  DOUBLE PRECISION, PARAMETER :: r_c     = 12.0e-6         !..mittlerer Radius (bei 1-Moment)
  DOUBLE PRECISION, PARAMETER :: rho_vel    = 0.5d0        !..Exponent in Dichtekorrektur
  DOUBLE PRECISION, PARAMETER :: rho_vel_c  = 1.0d0        !..fuer Wolkentropfen

  ! ... spezielle Parameter des KAMM2-Wolkenmoduls (Eisphase)

  DOUBLE PRECISION, PARAMETER :: e_ii  = 0.00              !..min. Eff.
  DOUBLE PRECISION, PARAMETER :: e_ic  = 0.80              !..max. Eff. fuer ice_cloud_riming
  DOUBLE PRECISION, PARAMETER :: e_sc  = 0.80              !..max. Eff. fuer snow_cloud_riming
  DOUBLE PRECISION, PARAMETER :: e_gc  = 1.00              !..max. Eff. fuer graupel_cloud_riming
  DOUBLE PRECISION, PARAMETER :: e_hc  = 1.00              !..max. Eff. fuer hail_cloud_riming
  DOUBLE PRECISION, PARAMETER :: e_min = 0.01              !..min. Eff. fuer gc,ic,sc
  DOUBLE PRECISION, PARAMETER :: alpha_spacefilling = 0.1  !..Raumerfuellungskoeff (max. 0.68)
  DOUBLE PRECISION, PARAMETER :: ice_s_vel  = 0.00         !..Dispersion der Fallgeschw. 
  DOUBLE PRECISION, PARAMETER :: snow_s_vel = 0.25         !..Dispersion der Fallgeschw.
  DOUBLE PRECISION, PARAMETER :: r_shedding = 500.0e-6     !..mittlerer Radius Shedding
  DOUBLE PRECISION, PARAMETER :: T_shed = 263.2

  DOUBLE PRECISION :: &
    na_dust    = 162.e3,   & ! initial number density of dust [1/m³], Phillips08 value 162e3
    na_soot    =  15.e6,   & ! initial number density of soot [1/m³], Phillips08 value 15e6
    na_orga    = 177.e6,   & ! initial number density of organics [1/m3], Phillips08 value 177e6
    ni_het_max = 100.0d3,   & ! max number of IN between 1-10 per liter, i.e. 1d3-10d3 
    ni_hom_max = 5000.0d3     ! number of liquid aerosols between 100-5000 per liter

  ! Look-up table for Phillips et al. nucleation
  INTEGER, PARAMETER :: &
    ttmax  = 30,      &  ! sets limit for temperature in look-up table
    ssmax  = 60,      &  ! sets limit for ice supersaturation in look-up table
    ttstep = 2,       &  ! increment for temperature in look-up table
    ssstep = 1           ! increment for ice supersaturation in look-up table

  REAL(KIND=8), DIMENSION(0:100,0:100), SAVE :: &
    afrac_dust, &  ! look-up table of activated fraction of dust particles acting as ice nuclei
    afrac_soot, &  ! ... of soot particles
    afrac_orga     ! ... of organic material

  INCLUDE 'phillips_nucleation_2010.incf'

  ! rain_freeze: der Teil des Regenspektrums kleiner als D_rainfrz_ig
  ! wird nach Gefrieren dem Eis zugeschlagen, der Teil von dort bis zu D_rainfrz_gh dem Graupel
  ! und der Rest dem Hagel.
  DOUBLE PRECISION, PARAMETER :: D_rainfrz_ig = 0.50d-3 !  rain --> ice oder graupel
  DOUBLE PRECISION, PARAMETER :: D_rainfrz_gh = 1.25d-3 ! rain --> graupel oder hail

  DOUBLE PRECISION, PARAMETER :: q_krit_ii = 1.000d-6 ! q-Schwellenwert fuer ice_selfcollection 
  DOUBLE PRECISION, PARAMETER :: D_krit_ii = 100.0d-6 ! D-Schwellenwert fuer ice_selfcollection
  DOUBLE PRECISION, PARAMETER :: D_conv_ii = 75.00d-6 ! D-Schwellenwert fuer ice_selfcollection
  DOUBLE PRECISION, PARAMETER :: q_krit_ic = 1.000d-5 ! q-Schwellenwert fuer ice_cloud_riming
  DOUBLE PRECISION, PARAMETER :: D_krit_ic = 150.0d-6 ! D-Schwellenwert fuer ice_cloud_riming
  DOUBLE PRECISION, PARAMETER :: q_krit_ir = 1.000d-5 ! q-Schwellenwert fuer ice_rain_riming
  DOUBLE PRECISION, PARAMETER :: D_krit_ir = 100.0d-6 ! D-Schwellenwert fuer ice_rain_riming
  DOUBLE PRECISION, PARAMETER :: q_krit_sc = 1.000d-5 ! q-Schwellenwert fuer snow_cloud_riming
  DOUBLE PRECISION, PARAMETER :: D_krit_sc = 150.0d-6 ! D-Schwellenwert fuer snow_cloud_riming
  DOUBLE PRECISION, PARAMETER :: q_krit_sr = 1.000d-5 ! q-Schwellenwert fuer snow_rain_riming
  DOUBLE PRECISION, PARAMETER :: D_krit_sr = 100.0d-6 ! D-Schwellenwert fuer snow_rain_riming
  DOUBLE PRECISION, PARAMETER :: q_krit_gc = 1.000d-6 ! q-Schwellenwert fuer graupel_cloud_riming
  DOUBLE PRECISION, PARAMETER :: D_krit_gc = 100.0d-6 ! D-Schwellenwert fuer graupel_cloud_riming
  DOUBLE PRECISION, PARAMETER :: q_krit_hc = 1.000d-6 ! q-Schwellenwert fuer hail_cloud_riming
  DOUBLE PRECISION, PARAMETER :: D_krit_hc = 100.0d-6 ! D-Schwellenwert fuer hail_cloud_riming
  DOUBLE PRECISION, PARAMETER :: q_krit_fr = 1.000d-6 ! q-Schwellenwert fuer rain_freeze
  DOUBLE PRECISION, PARAMETER :: q_krit_c  = 1.000d-6 ! q-Schwellenwert sonst
  DOUBLE PRECISION, PARAMETER :: q_krit    = 1.000d-9 ! q-Schwellenwert sonst
  DOUBLE PRECISION, PARAMETER :: D_conv_sg = 200.0d-6 ! D-Schwellenwert
  DOUBLE PRECISION, PARAMETER :: D_conv_ig = 200.0d-6 ! D-Schwellenwert
  DOUBLE PRECISION, PARAMETER :: x_conv    = 0.100d-9 ! minimale Graupel-/Hagelmasse riming
  DOUBLE PRECISION, PARAMETER :: D_shed_g  = 3.000d-3 ! D-Schwellenwert fuer graupel_shedding
  DOUBLE PRECISION, PARAMETER :: D_shed_h  = 5.000d-3 ! D-Schwellenwert fuer hagel_shedding
  DOUBLE PRECISION, PARAMETER :: D_krit_c  = 10.00d-6 ! D-Schwellenwert fuer cloud_collection
  DOUBLE PRECISION, PARAMETER :: D_coll_c  = 40.00d-6 ! oberer Wert fuer cloud_coll_eff
  DOUBLE PRECISION, PARAMETER :: T_nuc     = 273.2d+0 ! Temperatur ab der Eisnukleation einsetzt
  DOUBLE PRECISION, PARAMETER :: T_freeze  = 273.2d+0 ! Temperatur ab der Gefrieren einsetzt

  DOUBLE PRECISION, PARAMETER :: q_krit_aus = 1.00d-5 ! q-Schwellenwert fuer Ausgabe von D und Z 

  !..Jason's mu-Dm-relation for snow (see also Milbrandt & Yau 2005)
  DOUBLE PRECISION, PARAMETER :: snow_cmu1 = 4.5      ! a 
  DOUBLE PRECISION, PARAMETER :: snow_cmu2 = 0.5e+3   ! b 
  DOUBLE PRECISION, PARAMETER :: snow_cmu3 = 5.0e-3   ! Dnue
  DOUBLE PRECISION, PARAMETER :: snow_cmu4 = 5.5      ! a
  INTEGER,          PARAMETER :: snow_cmu5 = 1        ! exponent

  LOGICAL, PARAMETER          :: use_mu_Dm_snow_sedi = .FALSE.
  LOGICAL, PARAMETER          :: use_mu_Dm_rain_sedi = .TRUE.
  LOGICAL, PARAMETER          :: use_mu_Dm_rain_evap = .TRUE.
  LOGICAL, PARAMETER          :: use_mu_orig_rain_sedi = .FALSE.
  LOGICAL, PARAMETER          :: use_mu_orig_rain_evap = .FALSE.


  ! ub: neue Parametrisierung des Gefrierens von Regen 
  ! (Aufteilung in Eis und Hagel anstelle nur Hagel):
  LOGICAL, PARAMETER          :: use_rain_freeze_uli = .TRUE.

  ! ub: Schalter, ob Konversion von ice und snow nur dann zu graupel 
  ! (gesteuert durch alpha_spacefilling), 
  ! wenn die riming-rate groesser als die Depositionsrate ist 
  ! (ansonsten findet dies immer statt und wird nur durch alpha_spacefilling reguliert):
  LOGICAL, PARAMETER          :: use_ice_graupel_conv_uli = .TRUE.

  ! ... Parameter der Partikelklassen (Geometrie, Fallgeschwindigkeit, Verteilungfunktion, ...)

  ! new Graupel and Hail categories for comparison with/without hail
  ! nu and mu were chosen by comparison of size distributions from the bin-model by HUJI
  ! with generalised gamma-distribution
  ! avel, bvel derived from fit to kurves by Heymsfield and Kajikawa
  ! for hail-classes:
  !     ageo, bgeo in order to get constant density of 910, 750 and 600 kg/m**3
  !  

  TYPE(PARTICLE), PARAMETER :: graupelmedium = PARTICLE( & ! 'graupel im standard' 
       &                                 'graupelmedium' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.166666, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.10d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.300000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 7.64d+01,& !.a_vel..Koeff. Fallgesetz
       &                                 0.255200, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet  

  TYPE(PARTICLE), PARAMETER :: graupelleicht = PARTICLE( & ! 'graupel wenn hail
       &                                 'graupelleicht' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.166666, & !.mu.....Exp.-parameter der Verteil.
 !      &                                 1.000000, & !.mu.....Exp.-parameter der Verteil.
       &                                 2.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.30d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.300000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.66d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.257800, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hail900 = PARTICLE( & ! 'hail mit 920 kg/m**3'
       &                                 'hail900' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.166666, & !.mu.....Exp.-parameter der Verteil.
 !      &                                 1.000000, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-02, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.28d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 0.618d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.2141800, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hail750 = PARTICLE( & ! 'hail mit 750 kg/m**2'
       &                                 'hail750' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
!       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 0.166666, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-02, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.36d-01,& !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 0.580d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.2141500, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hail600 = PARTICLE( & ! 'hail mit 600 kg/m**3'
       &                                 'hail600' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.166666, & !.mu.....Exp.-parameter der Verteil.
 !      &                                 1.0000000, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-02, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.47d-01,& !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 0.535d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.2142000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

!----------------------------------------------------------------------------------------
! <<hn

  TYPE(PARTICLE), PARAMETER :: graupelhailA = PARTICLE( & ! 'graupel/hail Axel' 
       &                                 'graupel/hail_A' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.10d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.300000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.70d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.250000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet  

  TYPE(PARTICLE), PARAMETER :: graupelhail = PARTICLE( & ! 'graupel/hail Uli'
       &                                 'graupel/hail' ,& !.name...Bezeichnung
       &                                 0.666666, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 2.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.10d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.300000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.70d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.250000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelhail2 = PARTICLE( & ! 'graupel/hail Uli2'
       &                                 'graupel/hail2' ,& !.name...Bezeichnung
       &                                 0.333333, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 2.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.30d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.300000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.65d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.250000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelsb2006 = PARTICLE( & ! 'graupel/hail Uli2'
       &                                 'graupelsb2006' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.90d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.323000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.40d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.230000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelhail2test = PARTICLE( & ! 'graupel/hail Uli2'
       &                                 'graupel/hail2_test' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-03, & !.x_max..maximale Teilchenmasse
       &                                 1.00d-09, & !.x_min..minimale Teilchenmasse 
       &                                 1.50d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.323000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.33d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.165000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelhail2test2 = PARTICLE( & ! 'graupel/hail Uli2_2'
       &                                 'graupel/hail2_test2' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
!       &                                 1.00d-03, & !.x_max..maximale Teilchenmasse
       &                                 5.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 1.00d-09, & !.x_min..minimale Teilchenmasse 
       &                                 1.50d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.323000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.32d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.180000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelhail2test3 = PARTICLE( & ! 'graupel/hail Uli2_3'
       &                                 'graupel/hail2_test3' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
!       &                                 1.00d-03, & !.x_max..maximale Teilchenmasse
       &                                 5.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 1.00d-09, & !.x_min..minimale Teilchenmasse 
       &                                 1.42d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.314000, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.33d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.1870000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hailstandard = PARTICLE( & ! 'hail'
       &                                 'hailstandard' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 2.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.29d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 0.56d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.212700, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hailNOPPEL = PARTICLE( & ! 'hail'
       &                                 'hailNOPPEL' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-02, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.29d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 0.56d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.212700, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hailNOPPELlarge = PARTICLE( & ! 'hail'
       &                                 'hailNOPPELlarge' ,& !.name...Bezeichnung
       &                                 2.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-02, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.29d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 0.56d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.212700, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hailULI = PARTICLE( & ! 'hail'
       &                                 'hailULI' ,& !.name...Bezeichnung
       &                                 0.666667, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 4.00d-03, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-9, & !.x_min..minimale Teilchenmasse 
       &                                 0.127567, & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 42.31069, & !.a_vel..Koeff. Fallgesetz
       &                                 0.166667, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: hailULItest = PARTICLE( & ! 'hail'
       &                                 'hailULItest' ,& !.name...Bezeichnung
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
!       &                                 4.00d-03, & !.x_max..maximale Teilchenmasse
       &                                 5.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-9, & !.x_min..minimale Teilchenmasse 
       &                                 0.1366 , & !.a_geo..Koeff. Geometrie
       &                                 0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                                 39.3    , & !.a_vel..Koeff. Fallgesetz
       &                                 0.166667, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: cloudfastautoconv = PARTICLE( &
       &                               'cloudfastautoconv',  & !.name...Bezeichnung der Partikelklasse
       &                               0.0,      & !.nu.....Breiteparameter der Verteil.
       &                               0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: cloudstandard = PARTICLE( &! HN: made variable for seeding experiments
       &                               'cloudsstandard',  & !.name...Bezeichnung der Partikelklasse
       &                               0.333333, & !.nu.....Breiteparameter der Verteil.
       &                               0.666666, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: cloud_nue1mue1 = PARTICLE( &
       &                               'cloud_nue1mue1',  & !.name...Bezeichnung der Partikelklasse
       &                               1.000000, & !.nu.....Breiteparameter der Verteil.
       &                               1.000000, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: cloudveryslowautoconv = PARTICLE( &
       &                               'cloudveryslowautoconv',  & !.name...Bezeichnung der Partikelklasse
       &                               6.000000,      & !.nu.....Breiteparameter der Verteil.
       &                               1.000000, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: cloudslowautoconv = PARTICLE( &
       &                               'cloudslowautoconv',  & !.name...Bezeichnung der Partikelklasse
       &                               0.333333,      & !.nu.....Breiteparameter der Verteil.
       &                               1.000000, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: cloudheike = PARTICLE( & ! nu und mu aus Vergl. mit bin-Modell
       &                               'cloudheike',  & !.name...Bezeichnung der Partikelklasse
       &                               6.000000, & !.nu.....Breiteparameter der Verteil.
       &                               0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                               2.60d-10, & !.x_max..maximale Teilchenmasse D=80e-6m
       &                               4.20d-15, & !.x_min..minimale Teilchenmasse D=2.e-6m
       &                               1.24d-01, & !.a_geo..Koeff. Geometrie
       &                               0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                               3.75d+05, & !.a_vel..Koeff. Fallgesetz
       &                               0.666667, & !.b_vel..Koeff. Fallgesetz
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: iceCRY =  PARTICLE( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'iceCRY', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              5.00d-07, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-12, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              3.303633, & !.a_geo..Koeff. Geometrie
       &                              0.476191, & !.b_geo..Koeff. Geometrie = 1/2.1
       &                              2.77d+01, & !.a_vel..Koeff. Fallgesetz 
       &                              0.215790, & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: iceCRY2test =  PARTICLE( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'iceCRY2test', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu...e..Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              1.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-12, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              0.835000, & !.a_geo..Koeff. Geometrie
       &                              0.390000, & !.b_geo..Koeff. Geometrie = 1/2.1
       &                              2.77d+01, & !.a_vel..Koeff. Fallgesetz 
       &                              0.215790, & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: iceCirrus =  PARTICLE( & ! single bullets nach Heymsfield and Iaquinta
       &                              'iceCirrus', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu...e..Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              1.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-12, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              32.68941, & !.a_geo..Koeff. Geometrie, m(D)=7.08e-4*D**2.08
       &                              0.480769, & !.b_geo..Koeff. Geometrie = 1/2.08
       &                              3.447d+3, & !.a_vel..Koeff. Fallgesetz v(D)=3577*D**1.31
       &                              0.629808, & !.b_vel..Koeff. Fallgesetz = 1.31/2.08
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowCRYSTAL =  PARTICLE( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'snowCRYSTAL', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              1.00d-07, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-12, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              3.303633, & !.a_geo..Koeff. Geometrie
       &                              0.476191, & !.b_geo..Koeff. Geometrie = 1/2.1
       &                              2.47d+02, & !.a_vel..Koeff. Fallgesetz 
       &                              0.333333, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowCRYSTALuli =  PARTICLE( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'snowCRYSTALuli', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.500000, & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-05, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-10, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              2.400000, & !.a_geo..Koeff. Geometrie
       &                              0.455000, & !.b_geo..Koeff. Geometrie = 1/2.1
       &                              4.200000, & !.a_vel..Koeff. Fallgesetz 
       &                              0.092000, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowCRYSTALuli2 =  PARTICLE( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'snowCRYSTALuli2', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.500000, & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-05, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-10, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              2.400000, & !.a_geo..Koeff. Geometrie
       &                              0.455000, & !.b_geo..Koeff. Geometrie = 1/2.1
       &                              8.800000, & !.a_vel..Koeff. Fallgesetz 
       &                              0.150000, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowULI =  PARTICLE( & ! nach Andy Heymsfield (CRYSTAL-FACE)
       &                              'snowULI', & !.name...Bezeichnung der Partikelklasse
       &                              0.333333, & !.nu.....Breiteparameter der Verteil.
       &                              0.5     , & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.00d-12, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              3.303633, & !.a_geo..Koeff. Geometrie
       &                              0.476191, & !.b_geo..Koeff. Geometrie = 1/2.1
       &                              3.9259,   & !.a_vel..Koeff. Fallgesetz 
       &                              0.08,     & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: iceHK87 = PARTICLE( & ! HK87 'hex plates' + 700 D**1.2
       &                             'iceHK87',    & !.name...Bezeichnung der Partikelklasse
       &                            -0.333333, & !.nu.....Breiteparameter der Verteil.
       &                             0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                             1.00d-07, & !.x_max..maximale Teilchenmasse D=13.6e-3m
       &                             1.00d-12, & !.x_min..minimale Teilchenmasse D=60.0e-6m
       &                             2.17d-01, & !.a_geo..Koeff. Geometrie
       &                             0.302115, & !.b_geo..Koeff. Geometrie = 1/3.31
       &                             3.17d+02, & !.a_vel..Koeff. Fallgesetz
       &                             0.362500, & !.b_vel..Koeff. Fallgesetz = 1.2/3.31
       &                             0.860000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                             0.280000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                             3.141593)   !.cap....Koeff. Kapazitaet = pi

  TYPE(PARTICLE), PARAMETER :: iceHK = PARTICLE( & ! HK87 'hex plates'
       &                             'iceHK',    & !.name...Bezeichnung der Partikelklasse
       &                            -0.333333,      & !.nu.....Breiteparameter der Verteil.
       &                             0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                             1.00d-07, & !.x_max..maximale Teilchenmasse D=1.e-3m
       &                             1.00d-12, & !.x_min..minimale Teilchenmasse D=60e-6m
       &                             2.17d-01, & !.a_geo..Koeff. Geometrie
       &                             0.302115, & !.b_geo..Koeff. Geometrie = 1/3.31
       &                             4.19d+01, & !.a_vel..Koeff. Fallgesetz
       &                             0.260000, & !.b_vel..Koeff. Fallgesetz
       &                             0.860000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                             0.280000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                             3.141593)   !.cap....Koeff. Kapazitaet = pi

  TYPE(PARTICLE), PARAMETER :: rainREISNER = PARTICLE( & ! Reisner98
       &                              'rainREISNER',   & !.name...Bezeichnung der Partikelklasse
       &                             -0.666666, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              3.00d-06, & !.x_max..maximale Teilchenmasse
       &                              2.60d-10, & !.x_min..minimale Teilchenmasse D=80.e-6m
       &                              1.24d-01, & !.a_geo..Koeff. Geometrie
       &                              0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                              1.59d+02, & !.a_vel..Koeff. Fallgesetz
       &                              0.266667, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: rainAS = PARTICLE( & ! 
       &                              'rainAS',   & !.name...Bezeichnung der Partikelklasse
       &                              1.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              3.00d-06, & !.x_max..maximale Teilchenmasse
       &                              2.60d-10, & !.x_min..minimale Teilchenmasse D=80.e-6m
       &                              1.24d-01, & !.a_geo..Koeff. Geometrie
       &                              0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                              1.59d+02, & !.a_vel..Koeff. Fallgesetz
       &                              0.266667, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: rainUCLA = PARTICLE('rain', &
         1.000000, & !.nu.....Width parameter of the distribution
         0.333333, & !.mu.....exponential parameter of the distribution
         3.00e-06, & !.x_max..maximum particle mass
         2.60e-10, & !.x_min..minimale particler mass D=80.e-6m
         1.24e-01, & !.a_geo..coefficient of meteor geometry
         0.333333, & !.b_geo..coefficient of meteor geometry = 1/3
         1.59e+02, & !.a_vel..coefficient of fall velocity
         0.266667, & !.b_vel..coefficient of fall velocity
         0.780000, & !.a_ven..ventilation coefficient (PK, p.541)
         0.308000, & !.b_ven..ventilation coefficient (PK, p.541)
         2.0)        !.cap....capacity coefficient

  TYPE(PARTICLE), PARAMETER :: rainsb2006 = PARTICLE( & ! 
       &                              'rainsb2006',   & !.name...Bezeichnung der Partikelklasse
       &                             -0.666666, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              3.00d-06, & !.x_max..maximale Teilchenmasse
       &                              2.60d-10, & !.x_min..minimale Teilchenmasse D=80.e-6m
       &                              1.24d-01, & !.a_geo..Koeff. Geometrie
       &                              0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                              1.59d+02, & !.a_vel..Koeff. Fallgesetz
       &                              0.266667, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: rainULI = PARTICLE( & ! Blahak, v=v(x) gefittet 6.9.2005
       &                              'rainULI',   & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              3.00d-06, & !.x_max..maximale Teilchenmasse
       &                              2.60d-10, & !.x_min..minimale Teilchenmasse D=80.e-6m
       &                              1.24d-01, & !.a_geo..Koeff. Geometrie
       &                              0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                              114.0137, & !.a_vel..Koeff. Fallgesetz
       &                              0.234370, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: rainULI2 = PARTICLE( & ! Blahak, v=v(x) gefittet 6.9.2005
       &                              'rainULI2',   & !.name...Bezeichnung der Partikelklasse
       &                              1.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                              3.00d-06, & !.x_max..maximale Teilchenmasse
       &                              2.60d-10, & !.x_min..minimale Teilchenmasse D=80.e-6m
       &                              1.24d-01, & !.a_geo..Koeff. Geometrie
       &                              0.333333, & !.b_geo..Koeff. Geometrie = 1/3
       &                              114.0137, & !.a_vel..Koeff. Fallgesetz
       &                              0.234370, & !.b_vel..Koeff. Fallgesetz
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelHKVT = PARTICLE( & ! HK87 'lump graupel' mit mod. VT 
       &                                 'graupelHKVT',& !.name...Bezeichnung der Partikelklasse
       &                                 1.0     , & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.50d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.90d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.322580, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 0.40d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.230000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet


  TYPE(PARTICLE), PARAMETER :: graupelLH = PARTICLE( & ! LH74 'lump graupel 3' 
       &                                 'graupelLH',& !.name...Bezeichnung der Partikelklasse
       &                                 1.0     ,      & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 3.46d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.370370, & !.b_geo..Koeff. Geometrie = 1/2.7
       &                                 9.45d+00, & !.a_vel..Koeff. Fallgesetz
       &                                 0.120000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelHK = PARTICLE( & ! HK87 'cold lump graupel' 
       &                                 'graupelHK',& !.name...Bezeichnung der Partikelklasse
       &                                 1.0,      & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.77d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.311527, & !.b_geo..Koeff. Geometrie = 1/3.21
       &                                 4.64d+01, & !.a_vel..Koeff. Fallgesetz
       &                                 0.260000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelHKULI = PARTICLE( & ! HK87 'cold lump graupel' 
       &                                 'graupelHKULI',& !.name...Bezeichnung der Partikelklasse
       &                                 1.0,      & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.77d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.311527, & !.b_geo..Koeff. Geometrie = 1/3.21
       &                                 79.41,    & !.a_vel..Koeff. Fallgesetz
       &                                 0.27726,  & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: graupelHK87 = PARTICLE( & ! HK87 'lump graupel' 
       &                                 'graupelHK87',& !.name...Bezeichnung der Partikelklasse
       &                                 1.000000, & !.nu.....Breiteparameter der Verteil.
       &                                 0.333333, & !.mu.....Exp.-parameter der Verteil.
       &                                 1.00d-04, & !.x_max..maximale Teilchenmasse
       &                                 2.60d-10, & !.x_min..minimale Teilchenmasse 
       &                                 1.90d-01, & !.a_geo..Koeff. Geometrie
       &                                 0.322580, & !.b_geo..Koeff. Geometrie = 1/3.10
       &                                 1.01d+02, & !.a_vel..Koeff. Fallgesetz
       &                                 0.290000, & !.b_vel..Koeff. Fallgesetz
       &                                 0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                                 0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                                 2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowLHm = PARTICLE( & ! nach LHmod
       &                              'snowLHm', & !.name...Bezeichnung der Partikelklasse
       &                              0.000000, & !.nu.....Breiteparameter der Verteil.
       &                              0.5     , & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.73d-09, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              8.16d+00, & !.a_geo..Koeff. Geometrie
       &                              0.526316, & !.b_geo..Koeff. Geometrie = 1/1.9
       &                              1.50d+01, & !.a_vel..Koeff. Fallgesetz 
       &                              0.215790, & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowLH = PARTICLE( & ! nach Locatelli und Hobbs
       &                              'snowLH', & !.name...Bezeichnung der Partikelklasse
       &                              0.500000, & !.nu.....Breiteparameter der Verteil.
       &                              0.5     , & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.73d-09, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              8.16d+00, & !.a_geo..Koeff. Geometrie
       &                              0.526316, & !.b_geo..Koeff. Geometrie = 1/1.9
       &                              2.77d+01, & !.a_vel..Koeff. Fallgesetz 
       &                              0.215790, & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowLHULI = PARTICLE( & ! nach Locatelli und Hobbs
       &                              'snowLHULI', & !.name...Bezeichnung der Partikelklasse
       &                              0.5     , & !.nu.....Breiteparameter der Verteil.
       &                              0.5     , & !.mu.....Exp.-parameter der Verteil.
       &                              2.00d-06, & !.x_max..maximale Teilchenmasse D=???e-2m
       &                              1.73d-09, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                              8.16d+00, & !.a_geo..Koeff. Geometrie
       &                              0.526316, & !.b_geo..Koeff. Geometrie = 1/1.9
       &                              10.33,    & !.a_vel..Koeff. Fallgesetz 
       &                              0.15,     & !.b_vel..Koeff. Fallgesetz = 0.41/1.9
       &                              0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                              0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                              2.0)        !.cap....Koeff. Kapazitaet

  TYPE(PARTICLE), PARAMETER :: snowC = PARTICLE( & ! nach Cotton et. al (1986)
       &                               'snowC',& !.name...Bezeichnung der Partikelklasse
       &                               0.333333, & !.nu.....Breiteparameter der Verteil.
       &                               0.5     , & !.mu.....Exp.-parameter der Verteil.
       &                               5.00d-06, & !.x_max..maximale Teilchenmasse D=1.0e-2m
       &                               8.30d-11, & !.x_min..minimale Teilchenmasse D=200e-6m
       &                               1.33398 , & !.a_geo..Koeff. Geometrie
       &                               0.416667, & !.b_geo..Koeff. Geometrie = 1/2.4
       &                               2.876292, & !.a_vel..Koeff. Fallgesetz 
       &                               0.083333, & !.b_vel..Koeff. Fallgesetz = 0.2/2.4
       &                               0.780000, & !.a_ven..Koeff. Ventilation (PK, S.541)
       &                               0.308000, & !.b_ven..Koeff. Ventilation (PK, S.541)
       &                               2.0)        !.cap....Koeff. Kapazitaet


  TYPE(PARTICLE)                       :: cloud, rain, ice, snow, graupel, hail

  LOGICAL, PARAMETER                   :: isdebug   = .FALSE.
  LOGICAL, PARAMETER                   :: ice_multiplication = .TRUE.

  ! Graupel-Shedding: Das angefrorene Wasser wird bei T > T_shed, 
  ! ggf. nach enhanced-melting, wieder abgeworfen und zu Regen.
  ! Default is .false. for both parameter (see src_gscp)
  LOGICAL                              :: graupel_shedding
  LOGICAL                              :: hail_shedding

  LOGICAL, PARAMETER                   :: enhanced_melting   = .TRUE.

  ! Saturation adjustment within subroutine clouds() (.false.) or external after clouds() (.true.)
  ! (The latter requires extra care by the calling program!!!)
  LOGICAL, PARAMETER                   :: satad_nach_mikrophysik = .TRUE.

  INTEGER                              :: ice_typ,nuc_i_typ,nuc_c_typ,cloud_typ,mu_Dm_rain_typ

  DOUBLE PRECISION, POINTER, DIMENSION (:,:,:) :: rrho_04
  DOUBLE PRECISION, POINTER, DIMENSION (:,:,:) :: rrho_c

  DOUBLE PRECISION :: rain_cmu0,rain_cmu1,rain_cmu2,rain_cmu3,rain_cmu4,rain_gfak
  INTEGER          :: rain_cmu5

  DOUBLE PRECISION :: qnc_const, c_reffc, c_reffi

CONTAINS

  ! This subroutine has to be called once at the start of the model run by
  ! the main program. It properly sets the parameters for the different hydrometeor
  ! classes according to predefined parameter sets (see above).

  SUBROUTINE init_seifert( cloud_typ )
    IMPLICIT NONE
    INTEGER, INTENT(in) :: cloud_typ
    INTEGER :: wolke_typ

    wolke_typ = MOD(cloud_typ/10,100)

    IF (cloud_typ < 2000) THEN
      ! Without hail class:
      cloud = cloud_nue1mue1
      rain  = rainUCLA
      !rain  = rainULI
      ice   = iceCRY2test
      snow  = snowCRYSTALuli2
      graupel = graupelhail2test2
      ! Dummy value for hail:
      hail    = hailNOPPEL

      ! old settings
!      cloud = cloudstandard
!      rain  = rainAS
!      ice   = iceHK
!      snow  = snowCRYSTAL
!      graupel = graupelhail
    ELSE
      ! Including hail class:
      cloud = cloud_nue1mue1
      rain  = rainULI
      ice   = iceCirrus
      snow  = snowCRYSTALuli2
      graupel = graupelhail2test3
      hail    = hailULItest
    END IF

    !mu_Dm_rain_typ = mod(cloud_typ,10)
    mu_Dm_rain_typ = 2

    IF (mu_Dm_rain_typ.EQ.0) THEN
      !..constant mue value
      rain_cmu0 = 0.0
      rain_cmu1 = 0.0
      rain_cmu2 = 1.0
      rain_cmu3 = 1.0
      rain_cmu4 = (rain%nu+1.0d0)/rain%b_geo - 1.0d0 ! <-- this is the (constant) mue value 
      rain_cmu5 = 1
      rain_gfak = -1.0  ! In this case gamma = 1 in rain_evaporation    
    ELSEIF (mu_Dm_rain_typ.EQ.1) THEN
      !..Axel's mu-Dm-relation for raindrops based on 1d-bin model
      rain_cmu0 = 6.0             ! Axel's 2007 relation 
      rain_cmu1 = 30.0            ! 
      rain_cmu2 = 1.00d+3         ! 
      rain_cmu3 = 1.10d-3         ! D_eq,break
      rain_cmu4 = 1.0             ! 
      rain_cmu5 = 2               ! exponent
      rain_gfak = 1.0      
    ELSEIF (mu_Dm_rain_typ.EQ.2) THEN
      !..Modifikation of mu-Dm-relation for experiments with increased evaporation
      rain_cmu0 = 11.0            ! instead of 6.0      
      rain_cmu1 = 30.0            ! 
      rain_cmu2 = 1.00d+3         ! 
      rain_cmu3 = 1.10d-3         ! 
      rain_cmu4 = 4.0             ! instead of 1.0  
      rain_cmu5 = 2               ! 
      rain_gfak = 0.5             ! instead of 1.0      
    ELSEIF (mu_Dm_rain_typ.EQ.3) THEN
      !..Jason Milbrandts mu-Dm-relation'
      rain_cmu0 = 19.0           !
      rain_cmu1 = 19.0           ! Jason Milbrandt's mu-Dm-relation for rain
      rain_cmu2 = 0.60d+3        ! (Milbrandt&Yau 2005, JAS, Table 1)
      rain_cmu3 = 1.80d-3        !
      rain_cmu4 = 17.0           !
      rain_cmu5 = 1              !
      rain_gfak = -1.0  ! In this case gamma = 1 in rain_evaporation    
    ENDIF
      
!<AS20080304

  END SUBROUTINE init_seifert

  DOUBLE PRECISION FUNCTION power_explog_dd(a,b)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in)  :: a,b
    power_explog_dd = EXP(b*LOG(a))
  END FUNCTION power_explog_dd

  DOUBLE PRECISION FUNCTION power_explog_dr(a,b)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in)  :: a
    REAL, INTENT(in)  :: b
    power_explog_dr = EXP(b*LOG(a))
  END FUNCTION power_explog_dr

  DOUBLE PRECISION FUNCTION vent_coeff_a(parti,n)
    IMPLICIT NONE

    INTEGER        :: n
    TYPE(PARTICLE) :: parti

    vent_coeff_a = parti%a_ven * gfct((parti%nu+n+parti%b_geo)/parti%mu)              &
         &                  / gfct((parti%nu+1.0)/parti%mu)                        & 
         &                * ( gfct((parti%nu+1.0)/parti%mu)                        & 
         &                  / gfct((parti%nu+2.0)/parti%mu) )**(parti%b_geo+n-1.0) 

  END FUNCTION vent_coeff_a

  DOUBLE PRECISION FUNCTION vent_coeff_b(parti,n)
    IMPLICIT NONE

    INTEGER        :: n
    TYPE(PARTICLE) :: parti

    DOUBLE PRECISION, PARAMETER :: m_f   = 0.500   ! Koeff. Ventilationskoeff. (PK, S.541)

    vent_coeff_b = parti%b_ven                                                  & 
         & * gfct((parti%nu+n+(m_f+1.0)*parti%b_geo+m_f*parti%b_vel)/parti%mu)  &
         &             / gfct((parti%nu+1.0)/parti%mu)                          & 
         &           * ( gfct((parti%nu+1.0)/parti%mu)                          &
         &             / gfct((parti%nu+2.0)/parti%mu)                          &
         &             )**((m_f+1.0)*parti%b_geo+m_f*parti%b_vel+n-1.0)

  END FUNCTION vent_coeff_b

  DOUBLE PRECISION FUNCTION moment_gamma(p,n) 
    IMPLICIT NONE

    INTEGER          :: n
    TYPE(PARTICLE)   :: p

    moment_gamma  = gfct((n+p%nu+1.0)/p%mu) / gfct((p%nu+1.0)/p%mu)        &
         &     * ( gfct((  p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) )**n
  END FUNCTION moment_gamma

  DOUBLE PRECISION FUNCTION fracmoment_gamma(p,f) 
    IMPLICIT NONE

    DOUBLE PRECISION :: f
    TYPE(PARTICLE)   :: p

    fracmoment_gamma  = gfct((f+p%nu+1.0)/p%mu) / gfct((p%nu+1.0)/p%mu)        &
         &     * ( gfct((  p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) )**f
  END FUNCTION fracmoment_gamma

  DOUBLE PRECISION FUNCTION lambda_gamma(p,x) 
    IMPLICIT NONE

    DOUBLE PRECISION :: x
    TYPE(PARTICLE)   :: p

    lambda_gamma  = ( gfct((p%nu+1.0)/p%mu) / gfct((p%nu+2.0)/p%mu) * x)**(-p%mu)
  END FUNCTION lambda_gamma

  DOUBLE PRECISION FUNCTION e_es (t_)
    !*******************************************************************************
    !                        Saettigungsdampfdruck ueber Eis                       *
    !*******************************************************************************

    DOUBLE PRECISION, INTENT (IN) :: t_

    e_es  = e_3 * EXP (A_e * (t_ - T_3) / (t_ - B_e))

  END FUNCTION e_es

  DOUBLE PRECISION FUNCTION e_ws (t_)
    !*******************************************************************************
    !                      Saettigungsdampfdruck ueber Wasser                      *
    !*******************************************************************************

    DOUBLE PRECISION, INTENT (IN) :: t_

    e_ws  = e_3 * EXP (A_w * (t_ - T_3) / (t_ - B_w))

  END FUNCTION e_ws

  FUNCTION e_ws_vec (ta,idim,jdim)

    INTEGER :: idim, jdim
    DOUBLE PRECISION               :: e_ws_vec(idim,jdim)
    DOUBLE PRECISION, INTENT (IN)  :: ta(idim,jdim)

    e_ws_vec  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))

  END FUNCTION e_ws_vec

  FUNCTION e_es_vec (ta,idim,jdim)

    INTEGER :: idim, jdim
    DOUBLE PRECISION               :: e_es_vec(idim,jdim)
    DOUBLE PRECISION, INTENT (IN)  :: ta(idim,jdim)

    e_es_vec  = e_3 * EXP (A_e * (ta - T_3) / (ta - B_e))

  END FUNCTION e_es_vec

  DOUBLE PRECISION FUNCTION gfct(x)
    !*******************************************************************************
    !                                                                              *
    !       Gammafunktion aus Numerical Recipes (F77)                              *
    !                                                                              *
    !*******************************************************************************
    IMPLICIT NONE

    DOUBLE PRECISION cof(6)
    DOUBLE PRECISION stp,half,one,x,xx,fpf,tmp,ser,gamma
    INTEGER j

    DATA cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,  &
         &     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
    DATA half,one,fpf/0.5d0,1.0d0,5.5d0/

    xx  = x  - one
    tmp = xx + fpf
    tmp = (xx + half) * LOG(tmp) - tmp
    ser = one
    DO j = 1,6
      xx  = xx  + one
      ser = ser + cof(j) / xx
    ENDDO
    gamma = tmp + LOG(stp*ser)
    gamma = EXP(gamma)

    gfct = gamma
    RETURN
  END FUNCTION gfct

END MODULE wolken_konstanten

MODULE gamma_functions_mp_seifert

  IMPLICIT NONE

! UB_20081118>> structure for holding the data of a lookup table for the
!               incomplete gamma function:

  INTEGER, PARAMETER                     :: nlookup   = 2000    ! Internal number of bins (low res part)

  INTEGER, PARAMETER                     :: nlookuphr = 10000   ! Internal number of bins (high res part)

  ! dummy of internal number of bins (high res part) in case the high resolution part is not really needed:
  INTEGER, PARAMETER                     :: nlookuphr_dummy = 10 

  ! Type to hold the lookup table for the incomplete gamma functions.
  ! The table is divided into a low resolution part, which spans the
  ! whole range of x-values up to the 99.5 % x-value, and a high resolution part for the
  ! smallest 1 % of these x-values, where the incomplete gamma function may increase
  ! very rapidly and nonlinearily, depending on paramter a.
  ! For some applications (e.g., Newtons Method in future subroutine 
  ! graupel_hail_conv_wetgrowth_Dg_gamlook() ), this rapid change requires a much higher
  ! accuracy of the table lookup as compared to be achievable with the low resolution table.
  TYPE gamlookuptable
    ! Number of bins in the tables:
    INTEGER                              :: n           ! Internal number of bins (low res part)
    INTEGER                              :: nhr         ! Internal number of bins (high res part)
    DOUBLE PRECISION                     :: a           ! a-parameter
    DOUBLE PRECISION, DIMENSION(:), POINTER :: x        ! vector of x-parameters (limit of integration) - 
                                                        ! always starts at 0 and has equidistant dx (low resolution part)
    DOUBLE PRECISION, DIMENSION(:), POINTER :: xhr      ! vector of x-parameters (limit of integration) - 
                                                        ! always starts at 0 and has equidistant dxhr (high resolution part) 
    DOUBLE PRECISION                     :: dx          ! dx   (low resolution part)
    DOUBLE PRECISION                     :: dxhr        ! dxhr (high resolution part) 
    DOUBLE PRECISION                     :: odx         ! one over dx 
    DOUBLE PRECISION                     :: odxhr       ! one over dxhr 
    DOUBLE PRECISION, DIMENSION(:), POINTER :: igf      ! value of the inc. gamma function at (a,x) (low res)
    DOUBLE PRECISION, DIMENSION(:), POINTER :: igfhr    ! value of the inc. gamma function at (a,x) (high res)
  END TYPE gamlookuptable
! UB_20081118<<

  PRIVATE :: gcf, gser

CONTAINS

  DOUBLE PRECISION FUNCTION gammln(x)
  !*******************************************************************************
  !                                                                              *
  !       LOG(Gamma function) taken from Press et al.,  Numerical Recipes (F77)
  !
  !       Log(Gammafunktion) aus Numerical Recipes (F77)                         *
  !       (original)                                                             *
  !*******************************************************************************
    IMPLICIT NONE
      
    DOUBLE PRECISION, INTENT(in) :: x
    DOUBLE PRECISION, SAVE :: cof(6), stp
    DOUBLE PRECISION :: xx,tmp,ser
    INTEGER :: j
    DATA cof /76.18009172947146d0,-86.50532032941677d0, &
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5/
    DATA stp /2.5066282746310005d0/

    xx  = x
    tmp = xx + 5.5d0
    tmp = (xx + 0.5d0) * LOG(tmp) - tmp
    ser = 1.000000000190015d0
    DO j = 1,6
       xx  = xx  + 1.0d0
       ser = ser + cof(j) / xx
    ENDDO
    gammln = tmp + LOG(stp*ser/x)
    RETURN
  END FUNCTION gammln

  DOUBLE PRECISION FUNCTION gfct2(x)
  !*******************************************************************************
  !                                                                              *
  !       Gamma function taken from Press et al.,  Numerical Recipes (F77)
  !
  !       Gammafunktion aus Numerical Recipes (F77)                              *
  !       (etwas umformuliert, aber dieselben Ergebnisse wie obige Originalfunktion)
  !*******************************************************************************
    IMPLICIT NONE
      
    DOUBLE PRECISION, INTENT(in) :: x
    DOUBLE PRECISION, SAVE :: cof(6), stp, half, one, fpf
    DOUBLE PRECISION :: xx,tmp,ser,gamma
    INTEGER j
      
    DATA cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,  &
          &     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
    DATA half,one,fpf/0.5d0,1.0d0,5.5d0/
      
    xx  = x  - one
    tmp = xx + fpf
    tmp = (xx + half) * LOG(tmp) - tmp
    ser = one
    DO j = 1,6
       xx  = xx  + one
       ser = ser + cof(j) / xx
    ENDDO
    gamma = tmp + LOG(stp*ser)
    gamma = EXP(gamma)

    gfct2 = gamma
    RETURN
  END FUNCTION gfct2
  
  !*******************************************************************************
  !
  !       Incomplete Gamma function taken from Press et al.,  Numerical Recipes (F77)
  !
  !       Unvollstaendige Gammafunktion aus Numerical Recipes (F77)              *
  !       (etwas umformuliert, aber dieselben Ergebnisse wie obige Originalfunktion)
  !*******************************************************************************

  !*******************************************************************************
  ! 1) diverse Hilfsfunktionen:

  SUBROUTINE gcf(gammcf,a,x,gln)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ITMAX = 100
    DOUBLE PRECISION, PARAMETER :: EPS = 3.d-7, FPMIN = 1.d-30
    DOUBLE PRECISION, INTENT(in) :: a, x
    DOUBLE PRECISION, INTENT(out) :: gammcf, gln

    INTEGER :: i
    DOUBLE PRECISION :: an,b,c,d,del,h

    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    DO i=1,ITMAX
      an=-i*(i-a)
      b=b+2.0d0
      d=an*d+b
      IF (ABS(d).LT.FPMIN) d=FPMIN
      c=b+an/c
      IF (ABS(c).LT.FPMIN) c=FPMIN
      d=1./d
      del=d*c
      h=h*del
      IF (ABS(del-1.).LT.EPS) EXIT
    END DO

    IF (ABS(del-1.).GE.EPS) THEN
      WRITE (*,*) 'ERROR in GCF: a too large, ITMAX too small (MODULE gamma_functions, src_seifert.f90)'
      gammcf = 0.0d0
      RETURN
    END IF

    gammcf=EXP(-x+a*LOG(x)-gln)*h

    RETURN
  END SUBROUTINE gcf

  SUBROUTINE gser(gamser,a,x,gln)

    IMPLICIT NONE

    INTEGER, PARAMETER :: ITMAX = 100
    DOUBLE PRECISION, PARAMETER :: EPS=3.d-7
    DOUBLE PRECISION, INTENT(in) :: a, x
    DOUBLE PRECISION, INTENT(out) :: gamser, gln
      
    INTEGER :: n
    DOUBLE PRECISION :: ap,del,sum
    
    gln=gammln(a)
    IF (x.LE.0.) THEN
      IF (x.LT.0.) THEN
        WRITE (*,*) 'ERROR in GSER: x < 0 (MODULE gamma_functions, src_seifert.f90)'
      END IF
      gamser=0.0d0
      RETURN
    ENDIF

    ap=a
    sum=1./a
    del=sum
    DO n=1,ITMAX
      ap=ap+1.
      del=del*x/ap
      sum=sum+del
      IF (ABS(del).LT.ABS(sum)*EPS) EXIT
    END DO

    IF (ABS(del).GE.ABS(sum)*EPS) THEN
      WRITE (*,*) 'ERROR in GSER: a too large, ITMAX too small'
      WRITE (*,*) '  (MODULE gamma_functions, src_seifert.f90)'
      gamser = 0.0d0
      RETURN
    END IF

    gamser = sum*EXP(-x+a*LOG(x)-gln)

    RETURN
  END SUBROUTINE gser

  DOUBLE PRECISION FUNCTION gammp(a,x,gln)
  
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: a, x
    DOUBLE PRECISION, INTENT(out) :: gln
   
    DOUBLE PRECISION :: gammcf, gamser
    
    IF (x.LT.0.0d0 .OR. a .LE. 0.0d0) THEN
      WRITE(*,*) 'ERROR in GAMMP: bad arguments'
      WRITE(*,*) '  (MODULE gamma_functions, src_seifert.f90)'
      gammp = 0.0d0
      RETURN
    END IF

    IF (x .LT. a+1.)THEN
      CALL gser(gamser,a,x,gln)
      gammp = gamser
    ELSE
      CALL gcf(gammcf,a,x,gln)
      gammp = 1.0d0 - gammcf
    ENDIF
    RETURN
  END FUNCTION gammp

  DOUBLE PRECISION FUNCTION gammq(a,x,gln)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: a, x
    DOUBLE PRECISION, INTENT(out) :: gln
    DOUBLE PRECISION :: gammcf, gamser

    IF (x.LT.0.0d0 .OR. a .LE. 0.0d0) THEN
      WRITE(*,*) 'ERROR in GAMMQ: bad arguments (MODULE gamma_functions, src_seifert.f90)'
      gammq = 0.0d0
      RETURN
    END IF

    IF (x.LT.a+1.) THEN
      CALL gser(gamser,a,x,gln)
      gammq = 1.0d0 - gamser
    ELSE
      CALL gcf(gammcf,a,x,gln)
      gammq = gammcf
    ENDIF
    RETURN
  END FUNCTION gammq

  ! Ende diverse Hilfsfunktionen
  !*******************************************************************************

  !*******************************************************************************
  !
  ! Upper incomplete gamma function
  !
  ! Eigentliche obere unvollstaendige Gamma-Funktion, hier direkt
  ! das Integral
  !              int(x)(oo) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  DOUBLE PRECISION FUNCTION incgfct_upper(a,x)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: a, x
    DOUBLE PRECISION :: gam, gln
    
    gam = gammq(a,x,gln)
    incgfct_upper = EXP(gln) * gam

  END FUNCTION incgfct_upper

  !*******************************************************************************
  !
  ! Lower incomplete gamma function
  !
  ! Eigentliche untere unvollstaendige Gamma-Funktion, hier direkt
  ! das Integral
  !              int(0)(x) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  DOUBLE PRECISION FUNCTION incgfct_lower(a,x)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: a, x
    DOUBLE PRECISION :: gam, gln
    
    gam = gammp(a,x,gln)
    incgfct_lower = EXP(gln) * gam

  END FUNCTION incgfct_lower

  !*******************************************************************************
  !
  ! Eigentliche unvollstaendige Gamma-Funktion, hier direkt
  ! das Integral
  !              int(x1)(x2) exp(-t) t^(a-1) dt
  !
  !*******************************************************************************

  DOUBLE PRECISION FUNCTION incgfct(a,x1,x2)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: a, x1, x2
    
    incgfct = incgfct_lower(a,x2) - incgfct_lower(a,x1)

  END FUNCTION incgfct
  
  !*******************************************************************************
  !
  ! Create Lookup-table vectors for the lower incomplete gamma function,
  !              int(0)(x) exp(-t) t^(a-1) dt
  ! as function of x at constant a. 
  ! The table runs from x=0 to the 99.5 % - value of the normalized 
  ! incomplete gamma function. This 99.5 % - value has been fitted
  ! with high accuracy as function of a in the range a in [0;20], but can
  ! safely be applied also to higher values of a. (Fit created with the
  ! matlab-program "gamma_unvoll_lower_lookup.m" by Ulrich Blahak, 2008/11/13).
  !
  ! The last value in the table corresponds to x = infinity, so that
  ! during the reconstruction of incgfct-values from the table,
  ! the x-value can safely be truncated at the maximum table x-value.
  !
  !*******************************************************************************

  SUBROUTINE incgfct_lower_lookupcreate(a,ltable,nl,nlhr)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: a  ! value of a
    TYPE(gamlookuptable), INTENT(inout) :: ltable
    INTEGER, INTENT(in) :: nl, nlhr

    INTEGER :: i, err

    DOUBLE PRECISION, PARAMETER ::   &
         c1 =  36.629433904824623d0, &
         c2 = -0.119475603955226d0,  &
         c3 =  0.339332937820052d0,  &
         c4 =  1.156369000458310d0

    ! Store parameters in the structure ltable:
    ltable%a = a
    ltable%n = nl
    ltable%nhr = nlhr

    ! Allocate Memory for the table vectors:
    NULLIFY(ltable%x)
    NULLIFY(ltable%xhr)
    NULLIFY(ltable%igf)
    NULLIFY(ltable%igfhr)

    ALLOCATE(ltable%x(nl), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error x'
      STOP
    END IF
    ALLOCATE(ltable%xhr(nlhr), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error xhr'
      STOP
    END IF
    ALLOCATE(ltable%igf(nl), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error igf'
      STOP
    END IF
    ALLOCATE(ltable%igfhr(nlhr), STAT=err)
    IF (err /= 0) THEN
      WRITE (*,*) 'INCGFCT_LOWER_LOOKUPCREATE: Allocation error igfhr'
      STOP
    END IF

    !==================================================================
    ! low resolution part of the table:
    !==================================================================

    ! maximum x-value of the lookup table (99.5-%-value):
    ltable%x(ltable%n-1) = c1 * ( 1.0d0 - EXP(c2*a**c3) ) + c4*a

    ! create lookup table vectors:
    ltable%dx = ltable%x(ltable%n-1) / (ltable%n-2.0d0)
    ltable%odx = 1.0d0 / ltable%dx
!!! Diese Schleife vektorisiert nicht wg. incgfct_lower():
    DO i = 1, ltable%n - 1
      ltable%x(i) = (i-1) * ltable%dx
      ltable%igf(i) = incgfct_lower(a,ltable%x(i))
    END DO

    ! The last value is for x = infinity:
    ltable%x(ltable%n) = (ltable%n-1) * ltable%dx
    ltable%igf(ltable%n) = gfct2(a)

    !==================================================================
    ! high resolution part of the table (lowest 2 % of the X-values):
    !==================================================================

    ! create lookup table vectors:
    ltable%dxhr = ltable%x(NINT(0.01*(ltable%n-1))) / (ltable%nhr-1.0d0)
    ltable%odxhr = 1.0d0 / ltable%dxhr
!!! Diese Schleife vektorisiert nicht wg. incgfct_lower():
    DO i = 1, ltable%nhr
      ltable%xhr(i) = (i-1) * ltable%dxhr
      ltable%igfhr(i) = incgfct_lower(a,ltable%xhr(i))
    END DO

    RETURN
  END SUBROUTINE incgfct_lower_lookupcreate

  !*******************************************************************************
  !
  ! Retrieve values from a lookup table of the lower incomplete gamma function,
  ! as function of x at a constant a, for which the lookup table has been 
  ! created.
  !
  ! The last value in the table has to correspond to x = infinity, so that
  ! during the reconstruction of incgfct-values from the table,
  ! the x-value can safely be truncated at the maximum table x-value:
  !
  ! ltable%igf( ltable%x(ltable%n),...) = gfct(a)
  !
  ! Profiling with ifort on a Linux-PC shows, that table lookup for the
  ! incompl. gamma-Funktion is faster by a factor of about 15 compared
  ! to the original function without optimization (-O0). Using optimization
  ! could change this ratio (we encoutered up to 300 depending on function inlining).
  !
  ! Concerning the accuracy, comparisons show that the results of table lookup
  ! are accurate to within better than 0.1 % or even much less, except for
  ! very small values of X, for which the absolute values are however very
  ! close to 0. For X -> infinity (X > 99.5 % - value), accuracy may be 
  ! somewhat reduced up to about 0.5 % ,
  ! because the table is truncated at the 99.5 % value (second-last value)
  ! and the last value is set to the ordinary gamma function.
  !
  ! This function only uses the low resolution part of the table!
  !
  !*******************************************************************************

  DOUBLE PRECISION FUNCTION incgfct_lower_lookup(x, ltable)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: x  ! value of x for table lookup
    TYPE(gamlookuptable), INTENT(in) :: ltable

    INTEGER :: iu, io
    DOUBLE PRECISION :: xt

    ! Trunkcate x to the range of the table:
    xt = MAX(MIN(x, ltable%x(ltable%n)), 0.0d0)

    ! calculate indices of the neighbouring regular x-values
    ! in the table:
    iu = MIN(FLOOR(xt * ltable%odx) + 1, ltable%n-1)
    io = iu + 1

    ! interpolate linearily and subtract from the ordinary gamma function to get the upper
    ! incomplete gamma function:
    incgfct_lower_lookup = ltable%igf(iu) + &
         (ltable%igf(io) - ltable%igf(iu)) * ltable%odx * (xt-ltable%x(iu))

  END FUNCTION incgfct_lower_lookup

  ! Statt linearer Interpolation wird eine quadratische Interpolation gemacht
  ! und hierfuer am "kleinen" Ende der lookup table der 50-fach hoeher aufgeloeste
  ! high-resolution Ast benutzt.
  ! (an benachbarte 3 Punkte eine Parabel interpolieren und interpolierten Wert von der Parabel nehmen -- 
  !  weil es jeweils 2 moegliche 3-Punkte-Nachbarschaften gibt, wird aus Stetigkeitsgruenden der Mittelwert
  ! von beiden genommen):
  DOUBLE PRECISION FUNCTION incgfct_lower_lookup_parabolic(x, ltable)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x  ! value of x for table lookup
    TYPE(gamlookuptable), INTENT(in) :: ltable

    INTEGER :: iu, im, io
    DOUBLE PRECISION :: xt, f12, f23, f123, yn1, yn2

    ! If x is within the high-resolution part of the table:
    IF (x <= ltable%xhr(ltable%nhr)) THEN

      ! Truncate x to the range of the table:
      xt = MAX(x, 0.0d0)

      ! Neighbourhood 1:
      ! calculate indices of the neighbouring regular x-values
      ! in the table:
      im = MIN(MAX(FLOOR(xt * ltable%odxhr) + 1, 2),  ltable%nhr-1)
      iu = im - 1
      io = im + 1

      ! interpolate parabolicly:
      ! coefficients of the interpolating parabola wrsp. to Newton's polynome form with Newton's tableau ("divided differences"):
      f12  = (ltable%igfhr(im) - ltable%igfhr(iu)) * ltable%odxhr
      f23  = (ltable%igfhr(io) - ltable%igfhr(im)) * ltable%odxhr
      f123 = 0.5d0 * (f23 - f12) * ltable%odxhr
      ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
      yn1 = f123 * (xt - ltable%xhr(im)) + f12
      yn1 = yn1  * (xt - ltable%xhr(iu)) + ltable%igfhr(iu)

      IF (im < ltable%nhr - 1) THEN
        ! Neighbourhood 2:
        ! calculate indices of the neighbouring regular x-values
        ! in the table:

        im = MIN(MAX(CEILING(xt * ltable%odxhr) + 1, 2),  ltable%nhr-1)
        iu = im - 1
        io = im + 1

        ! interpolate parabolicly:
        ! coefficients of the interpolating parabola wrsp. to Newton's polynome form with Newton's tableau ("divided differences"):
        f12  = (ltable%igfhr(im) - ltable%igfhr(iu)) * ltable%odxhr
        f23  = (ltable%igfhr(io) - ltable%igfhr(im)) * ltable%odxhr
        f123 = 0.5d0 * (f23 - f12) * ltable%odxhr
        ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
        yn2 = f123 * (xt - ltable%xhr(im)) + f12
        yn2 = yn2  * (xt - ltable%xhr(iu)) + ltable%igfhr(iu)

        incgfct_lower_lookup_parabolic = 0.5d0 * ( yn1 + yn2 )

      ELSE

        incgfct_lower_lookup_parabolic = yn1

      END IF

    ELSE

      ! Truncate x to the range of the table:
      xt = MAX(MIN(x, ltable%x(ltable%n)), 0.0d0)

      ! Neighbourhood 1:
      ! calculate indices of the neighbouring regular x-values
      ! in the table:
      im = MIN(MAX(FLOOR(xt * ltable%odx) + 1, 2),  ltable%n-1)
      iu = im - 1
      io = im + 1

      ! interpolate parabolicly:
      ! coefficients of the interpolating parabola wrsp. to Newton's polynome form with Newton's tableau ("divided differences"):
      f12  = (ltable%igf(im) - ltable%igf(iu)) * ltable%odx
      f23  = (ltable%igf(io) - ltable%igf(im)) * ltable%odx
      f123 = 0.5d0 * (f23 - f12) * ltable%odx
      ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
      yn1 = f123 * (xt - ltable%x(im)) + f12
      yn1 = yn1  * (xt - ltable%x(iu)) + ltable%igf(iu)

      IF (im < ltable%n - 2) THEN
        ! Neighbourhood 2:
        ! calculate indices of the neighbouring regular x-values
        ! in the table:

        im = MIN(MAX(CEILING(xt * ltable%odx) + 1, 2),  ltable%n-2)
        iu = im - 1
        io = im + 1

        ! interpolate parabolicly:
        ! coefficients of the interpolating parabola wrsp. to Newton's polynome form with Newton's tableau ("divided differences"):
        f12  = (ltable%igf(im) - ltable%igf(iu)) * ltable%odx
        f23  = (ltable%igf(io) - ltable%igf(im)) * ltable%odx
        f123 = 0.5d0 * (f23 - f12) * ltable%odx
        ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
        yn2 = f123 * (xt - ltable%x(im)) + f12
        yn2 = yn2  * (xt - ltable%x(iu)) + ltable%igf(iu)

        incgfct_lower_lookup_parabolic = 0.5d0 * ( yn1 + yn2 )

      ELSE

        incgfct_lower_lookup_parabolic = yn1

      END IF

    END IF

    incgfct_lower_lookup_parabolic = MAX(incgfct_lower_lookup_parabolic, 0.0d0)

  END FUNCTION incgfct_lower_lookup_parabolic

  !*******************************************************************************
  !
  ! Retrieve values of the upper incomplete gamma function 
  ! from a lookup table of the lower incomplete gamma function,
  ! as function of x at a constant a, for which the lookup table has been 
  ! created.
  !
  ! The last value in the table has to correspond to x = infinity 
  ! (the ordinary gamma function of a), so that
  ! during the reconstruction of incgfct-values from the table,
  ! the x-value can safely be truncated at the maximum table x-value:
  !
  ! ltable%igf( ltable%x(ltable%n),...) = gfct(a)
  !
  ! This function only uses the low resolution part of the table!
  !
  !*******************************************************************************

  DOUBLE PRECISION FUNCTION incgfct_upper_lookup(x, ltable)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: x  ! value of x for table lookup
    TYPE(gamlookuptable), INTENT(in) :: ltable

    INTEGER :: iu, io
    DOUBLE PRECISION :: xt

    ! Trunkcate x to the range of the table:
    xt = MAX(MIN(x, ltable%x(ltable%n)), 0.0d0)

    ! calculate indices of the neighbouring regular x-values
    ! in the table:
    iu = MIN(FLOOR(xt * ltable%odx) + 1, ltable%n-1)
    io = iu + 1

    ! interpolate lower inc. gamma function linearily and subtract from 
    ! the ordinary gamma function to get the upper
    ! incomplete gamma function:
    incgfct_upper_lookup = ltable%igf(ltable%n) - ltable%igf(iu) -  &
         (ltable%igf(io) - ltable%igf(iu)) * ltable%odx * (xt-ltable%x(iu))

    ! Aufgrund von Rundungsfehlern (Differenz von 2 fast gleichen Zahlen) kann es beim table lookup passieren,
    ! dass incgfct_upper_lookup(x, ltable) kleiner 0 wird, wenn eigentlich nahezu 0.0 herauskommen muesste.
    ! Dies kommt vor allem dann vor, wenn x sehr gross ist.
    ! Deswegen Begrenzung:

    incgfct_upper_lookup = MAX(incgfct_upper_lookup, 0.0d0)

    RETURN
  END FUNCTION incgfct_upper_lookup

  ! Statt linearer Interpolation wird eine quadratische Interpolation gemacht 
  ! und hierfuer am "kleinen" Ende der lookup table der 50-fach hoeher aufgeloeste
  ! high-resolution Ast benutzt.
  ! (an benachbarte 3 Punkte eine Parabel interpolieren und interpolierten Wert von der Parabel nehmen -- 
  !  weil es jeweils 2 moegliche 3-Punkte-Nachbarschaften gibt, wird aus Stetigkeitsgruenden der Mittelwert
  ! von beiden genommen):
  DOUBLE PRECISION FUNCTION incgfct_upper_lookup_parabolic(x, ltable)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x  ! value of x for table lookup
    TYPE(gamlookuptable), INTENT(in) :: ltable

    INTEGER :: iu, im, io
    DOUBLE PRECISION :: xt, f12, f23, f123, yn1, yn2

    ! If x is within the high-resolution part of the table:
    IF (x <= ltable%xhr(ltable%nhr)) THEN

      ! Truncate x to the range of the table:
      xt = MAX(x, 0.0d0)

      ! Neighbourhood 1:
      ! calculate indices of the neighbouring regular x-values
      ! in the table:
      im = MIN(MAX(FLOOR(xt * ltable%odxhr) + 1, 2),  ltable%nhr-1)
      iu = im - 1
      io = im + 1

      ! interpolate parabolicly (lower incomplete gamma function -- 
      ! will be converted to upper function later):
      ! coefficients of the interpolating parabola wrsp. to Newton's polynome form
      ! with Newton's tableau ("divided differences"):
      f12  = (ltable%igfhr(im) - ltable%igfhr(iu)) * ltable%odxhr
      f23  = (ltable%igfhr(io) - ltable%igfhr(im)) * ltable%odxhr
      f123 = 0.5d0 * (f23 - f12) * ltable%odxhr
      ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
      yn1 = f123 * (xt - ltable%xhr(im)) + f12
      yn1 = yn1  * (xt - ltable%xhr(iu)) + ltable%igfhr(iu)

      IF (im < ltable%nhr - 1) THEN
        ! Neighbourhood 2:
        ! calculate indices of the neighbouring regular x-values
        ! in the table:

        im = MIN(MAX(CEILING(xt * ltable%odxhr) + 1, 2),  ltable%nhr-1)
        iu = im - 1
        io = im + 1

        ! interpolate parabolicly:
        ! coefficients of the interpolating parabola wrsp. to Newton's polynome form
        ! with Newton's tableau ("divided differences"):
        f12  = (ltable%igfhr(im) - ltable%igfhr(iu)) * ltable%odxhr
        f23  = (ltable%igfhr(io) - ltable%igfhr(im)) * ltable%odxhr
        f123 = 0.5d0 * (f23 - f12) * ltable%odxhr
        ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
        yn2 = f123 * (xt - ltable%xhr(im)) + f12
        yn2 = yn2  * (xt - ltable%xhr(iu)) + ltable%igfhr(iu)

        incgfct_upper_lookup_parabolic = 0.5d0 * ( yn1 + yn2 )

      ELSE

        incgfct_upper_lookup_parabolic = yn1

      END IF

    ELSE

      ! Truncate x to the range of the table:
      xt = MAX(MIN(x, ltable%x(ltable%n)), 0.0d0)

      ! Neighbourhood 1:
      ! calculate indices of the neighbouring regular x-values
      ! in the table:
      im = MIN(MAX(FLOOR(xt * ltable%odx) + 1, 2),  ltable%n-1)
      iu = im - 1
      io = im + 1

      ! interpolate parabolicly (lower incomplete gamma function --
      ! will be converted to upper FUNCTION later):
      ! coefficients of the interpolating parabola wrsp. to Newton's polynome form
      ! with Newton's tableau ("divided differences"):
      f12  = (ltable%igf(im) - ltable%igf(iu)) * ltable%odx
      f23  = (ltable%igf(io) - ltable%igf(im)) * ltable%odx
      f123 = 0.5d0 * (f23 - f12) * ltable%odx
      ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
      yn1 = f123 * (xt - ltable%x(im)) + f12
      yn1 = yn1  * (xt - ltable%x(iu)) + ltable%igf(iu)

      IF (im < ltable%n - 2) THEN
        ! Neighbourhood 2:
        ! calculate indices of the neighbouring regular x-values
        ! in the table:

        im = MIN(MAX(CEILING(xt * ltable%odx) + 1, 2),  ltable%n-2)
        iu = im - 1
        io = im + 1

        ! interpolate parabolicly:
        ! coefficients of the interpolating parabola wrsp. to Newton's polynome form
        ! with Newton's tableau ("divided differences"):
        f12  = (ltable%igf(im) - ltable%igf(iu)) * ltable%odx
        f23  = (ltable%igf(io) - ltable%igf(im)) * ltable%odx
        f123 = 0.5d0 * (f23 - f12) * ltable%odx
        ! Use Horner's scheme to calculate the interpolated value from the interpolating parabola:
        yn2 = f123 * (xt - ltable%x(im)) + f12
        yn2 = yn2  * (xt - ltable%x(iu)) + ltable%igf(iu)

        incgfct_upper_lookup_parabolic = 0.5d0 * ( yn1 + yn2 )

      ELSE

        incgfct_upper_lookup_parabolic = yn1

      END IF

    END IF

    ! Convert to upper incomplete gamma function:
    incgfct_upper_lookup_parabolic = ltable%igf(ltable%n) - incgfct_upper_lookup_parabolic

    incgfct_upper_lookup_parabolic = MAX(incgfct_upper_lookup_parabolic, 0.0d0)

    RETURN
  END FUNCTION incgfct_upper_lookup_parabolic

END MODULE gamma_functions_mp_seifert


!=======================================================================

MODULE wolken_driver

  USE wolken_konstanten

  IMPLICIT NONE

  ! ... Parameter fuer Wolken ...
  LOGICAL          :: wolke
  INTEGER          :: wolke_typ=2603,ccn_typ

  ! ub>>
  ! Schalter fuer die Ausgabe von horizontal gemittelten Umwandlungsraten 
  ! fuer jeden einzelnen mikrophysikalischen Prozess: Initialisierung mit .FALSE.
  LOGICAL          :: speichere_umwandlungsraten = .FALSE.
  LOGICAL          :: speichere_dqdt = .FALSE.
  LOGICAL          :: speichere_precipstat = .FALSE.
  ! ub<<

  ! ... Gitter ...
  INTEGER          :: loc_ix, loc_iy, loc_iz
  DOUBLE PRECISION :: dz

  ! ... Zeitschritt ....
  DOUBLE PRECISION :: dt

  ! ... Pi ...
  DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793238462643383d0

  DOUBLE PRECISION, PARAMETER :: r  = 287.04D0
  DOUBLE PRECISION, PARAMETER :: r1 = 461.50D0
  DOUBLE PRECISION, PARAMETER :: cp = 1005.7D0
  DOUBLE PRECISION, PARAMETER :: cv = 718.66D0


  ! ... Felder fuer die Ausgabe der Umwandlungsraten ...
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: dqdt
  
  ! Anzahl der Umwandlungsraten + 18 
  ! (1 fuer die Heizrate, 1 fuer Temperaturprofil, 4 fuer Uebersaett. vor/nach Mikrophysik , 
  !  6 fuer die Q's, 6 fuer die QN's); 
  ! Die Umwandlungsraten beginnen bei 1 und gehen bis nrates-18;
  ! nrates-17 und nrates-16 ist fuer die Uebersaettigung geg. Wasser vor und nach der Mikrophysik;
  ! nrates-15 und nrates-14 ist fuer die Uebersaettigung geg. Eis vor und nach der Mikrophysik;
  ! nrates-13 und nrates-12 ist fuer die Heizrate und fuers T-Profil;
  ! nrates-11 bis nrates-6 ist fuer QC,QR,QI,QS,QG und QH;
  ! nrates-5 bis nrates ist fuer QNC,QNR,QNI,QNS,QNG und QNH
  INTEGER, PARAMETER :: nmicrorates = 63
  INTEGER, PARAMETER :: nsedifluxdiv = 12
  INTEGER, PARAMETER :: nrates = nmicrorates + nsedifluxdiv + 18

  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: cond_neu_sb, evap_neu_sb

  ! ... Felder ...
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: &
       & w, p, t, rho, q,                      &
       & q_cloud, q_ice, q_rain, q_graupel, q_snow, q_hail, &
       & n_cloud, n_ice, n_rain, n_graupel, n_snow, n_hail, &
       & S_w, S_i, dSwdz, dSidz, dT0dz

  ! ... Grundzustandsvariablen ...
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::   &
       & p_0, t_0, rho_0

  ! ... Referenzzustand ...
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::   &
       & p_g, t_g, rho_g

  ! ... Metrikkoeff ...
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::   &
       & w_g, x3_x3

  ! ... vertical grid ...
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::   &
       & dz3d

  ! w at cloud base and height of model level for cloud_nucleation_SK
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::   &
       & zml_k, w_cb

  ! ub>>
  ! Speicherfelder fuer die Depositionsraten von Eis und Schnee fuer die
  ! Entscheidung, ob Eis durch Riming Eis bleibt oder zu Graupel werden kann:
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::   &
       deprate_ice, deprate_snow, &
       rimeqcrate_ice, rimencrate_ice, rimeqrrate_ice, rimeqirate_ice, rimenrrate_ice, &
       rimeqcrate_snow, rimencrate_snow, rimeqrrate_snow, rimeqirate_snow, rimenrrate_snow, &
       d_id_sp, d_sd_sp, d_rd_sp_ice, d_rd_sp_snow
  ! ub<<

  ! ... Regenraten ...
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: &
       & prec, prec_s, prec_cloud, prec_ice, &
       & prec_rain, prec_snow, prec_graupel, prec_hail

  ! ... Variables for wet growth diameter lookup tables:
  DOUBLE PRECISION, DIMENSION(:,:,:,:),   ALLOCATABLE :: &
       & dmin_wg_g
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: &
       & pvec_wg_g, Tvec_wg_g, qwvec_wg_g, qivec_wg_g
  INTEGER :: anzp_wg, anzT_wg, anzi_wg, anzw_wg
! UB_20090227>>
! Struct to hold the new equidistant lookup table for graupel wetgrowth diameter:
  TYPE(lookupt_4d) :: ltabdminwgg
! UB_20090227<<



CONTAINS

  ! ... some subroutines for bulk cloud scheme

  ! Equation of state

  FUNCTION dichte(p,t,q,q_cloud,q_ice,q_rain,q_snow,q_graupel,q_hail)         
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(0:loc_ix,1:loc_iy,1:loc_iz), INTENT(in) :: &
         & p, t, q
    DOUBLE PRECISION, DIMENSION(0:loc_ix,1:loc_iy,1:loc_iz), OPTIONAL, INTENT(in) :: &
         & q_cloud, q_ice, q_rain, q_snow, q_graupel, q_hail
    DOUBLE PRECISION, DIMENSION(0:loc_ix,1:loc_iy,1:loc_iz) :: &
         & dichte

    dichte = (p / p_0 - T / T_0) * rho_0 - (R_d / R_l - 1.0d0) * q &
         &                 + q_cloud + q_ice + q_rain + q_snow + q_graupel + q_hail 


  END FUNCTION dichte

  SUBROUTINE alloc_driver ()

    IMPLICIT NONE

    INTEGER       :: stat_var = 0


    ALLOCATE(     w(0:loc_ix,1:loc_iy,1:loc_iz), &
         &        p(0:loc_ix,1:loc_iy,1:loc_iz), &
         &        T(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      rho(0:loc_ix,1:loc_iy,1:loc_iz), &
         &        q(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      p_0(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      T_0(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    rho_0(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      p_g(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      T_g(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    rho_g(0:loc_ix,1:loc_iy,1:loc_iz), &
         &     dz3d(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      S_w(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      S_i(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    dSwdz(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    dSidz(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    dT0dz(0:loc_ix,1:loc_iy,1:loc_iz), &
         &      w_g(0:loc_ix,1:loc_iy,1:loc_iz), &
         &    x3_x3(0:loc_ix,1:loc_iy,1:loc_iz), STAT = stat_var )



    ! >>hn height zml_k and updraft at cloud base w_cb needed for
    !      cloud drople nucleation
    ALLOCATE( zml_k(0:loc_ix,1:loc_iy,1:loc_iz), STAT = stat_var)        
    IF (stat_var /= 0) THEN
      WRITE (*,*) 'Fehler bei Allokierung zml_k ', stat_var
      STOP
    END IF
    ALLOCATE( w_cb(0:loc_ix,1:loc_iy,1:loc_iz), STAT = stat_var)        
    IF (stat_var /= 0) THEN
      WRITE (*,*) 'Fehler bei Allokierung w_cb ', stat_var
      STOP
    END IF
    ! <<hn

    ! ub>>
    IF (use_ice_graupel_conv_uli) THEN
      ALLOCATE( deprate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           deprate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           STAT = stat_var)        
      IF (stat_var /= 0) THEN
        WRITE (*,*) 'Fehler bei Allokierung deprate_ice, deprate_snow ', stat_var
        STOP
      END IF
      ALLOCATE( rimeqcrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimencrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqirate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqrrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimenrrate_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqcrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimencrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqirate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimeqrrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           rimenrrate_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_id_sp(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_sd_sp(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_rd_sp_ice(0:loc_ix,1:loc_iy,1:loc_iz), &
           d_rd_sp_snow(0:loc_ix,1:loc_iy,1:loc_iz), &
           STAT = stat_var)        
      IF (stat_var /= 0) THEN
        WRITE (*,*) 'Fehler bei Allokierung rimerate_x ', stat_var
        STOP
      END IF
    END IF
    ! ub<<

    wolke = .TRUE.

    w_g   = 1.0
    x3_x3 = 1.0

    w     = 0.0
    p     = 0.0
    T     = 0.0
    rho   = 0.0
    q     = 0.0
    p_0   = 1d5
    T_0   = 3d2
    rho_0 = 1.2
    p_g   = 0.0
    T_g   = 0.0
    rho_g = 0.0
    dz3d  = 1.0

    S_w   = 0.0
    S_i   = 0.0
    dSwdz = 0.0
    dSidz = 0.0
    dT0dz = 0.0

    zml_k = 0.0
    w_cb = 0.0

    IF (use_ice_graupel_conv_uli) THEN
      deprate_ice = 0.0d0
      deprate_snow = 0.0d0
      rimeqcrate_ice = 0.0d0
      rimencrate_ice = 0.0d0
      rimeqirate_ice = 0.0d0
      rimeqrrate_ice = 0.0d0
      rimenrrate_ice = 0.0d0
      rimeqcrate_snow = 0.0d0
      rimencrate_snow = 0.0d0
      rimeqirate_snow = 0.0d0
      rimeqrrate_snow = 0.0d0
      rimenrrate_snow = 0.0d0
      d_id_sp = 0.0d0
      d_sd_sp = 0.0d0
      d_rd_sp_ice = 0.0d0
      d_rd_sp_snow = 0.0d0
    END IF

  END SUBROUTINE alloc_driver

  SUBROUTINE dealloc_driver ()
    IMPLICIT NONE

    DEALLOCATE(w,p,T,rho,q,p_0,T_0,rho_0,p_g,T_g,rho_g,w_g,x3_x3,dz3d,&
         & S_w,S_i,dSwdz,dSidz,dT0dz,rrho_04,rrho_c)
    ! >> hn
    DEALLOCATE(zml_k, w_cb)
    ! << hn
    ! ub>>
    IF (use_ice_graupel_conv_uli) THEN
      DEALLOCATE(deprate_ice,deprate_snow,&
           rimeqcrate_ice,&
           rimencrate_ice,&
           rimeqirate_ice,&
           rimeqrrate_ice,&
           rimenrrate_ice,&
           rimeqcrate_snow,&
           rimencrate_snow,&
           rimeqirate_snow,&
           rimeqrrate_snow,&
           rimenrrate_snow,&
           d_id_sp,&
           d_sd_sp,&
           d_rd_sp_ice,&
           d_rd_sp_snow &
           )
    END IF
    ! ub>>
  END SUBROUTINE dealloc_driver

END MODULE wolken_driver


!****************************************************************************
!
! Some dummy modules to mimic a KAMM2 environment
!
!****************************************************************************

MODULE globale_variablen
  USE wolken_driver, ONLY: dt, loc_ix, loc_iy, loc_iz, &
       & w, T, p, q, rho, T_g, p_g, rho_g,       &
       & q_cloud, n_cloud, prec_cloud,   &
       & q_rain, n_rain, prec_rain,      &
       & q_ice, n_ice, prec_ice,         &
       & q_snow, n_snow, prec_snow,      &
       & q_graupel, n_graupel, prec_graupel,  &
       & q_hail, n_hail, prec_hail,  &
       & prec,prec_s, dz3d, S_w, dSwdz, S_i, dSidz, dT0dz, zml_k, w_cb, &
       & dqdt, speichere_umwandlungsraten, speichere_dqdt, nrates, nmicrorates, &
       & nsedifluxdiv, speichere_precipstat, cond_neu_sb, evap_neu_sb, &
       & dmin_wg_g, pvec_wg_g, Tvec_wg_g, qwvec_wg_g, qivec_wg_g, &
       & anzp_wg, anzT_wg, anzi_wg, anzw_wg, &
       & ltabdminwgg, &
       & deprate_ice,deprate_snow, rimeqcrate_ice, rimencrate_ice, rimeqirate_ice, &
       & rimeqrrate_ice, rimenrrate_ice, rimeqcrate_snow, rimencrate_snow, rimeqirate_snow, &
       & rimeqrrate_snow, rimenrrate_snow, d_id_sp, d_sd_sp, d_rd_sp_ice, d_rd_sp_snow
! ub<<      
END MODULE globale_variablen

MODULE konstanten
  USE wolken_driver, ONLY: dz, wolke_typ, pi, r, r1, cv, cp
END MODULE konstanten

MODULE geometrie
  USE wolken_driver, ONLY: w_g,x3_x3
END MODULE geometrie

MODULE parallele_umgebung

#if FOR_LM == 1
  USE data_parallel, ONLY : my_cart_id, num_compute
#else
  use mpi_interface, only : my_cart_id=>myid  ! proc id
#endif

  IMPLICIT NONE

  INTEGER :: mpi_max,mpi_min,mpi_sum

CONTAINS

  LOGICAL FUNCTION isIO()

    IF (my_cart_id == 0) THEN
      isIO = .TRUE.
    ELSE
      isIO = .FALSE.
    ENDIF

  END FUNCTION isIO

  SUBROUTINE abortparallel(text, ierror)
    INTEGER, INTENT(in) :: ierror
    CHARACTER*(*), INTENT(in) :: text

    WRITE (*,'(A,I6)') text,ierror 
    STOP
  END SUBROUTINE abortparallel

  SUBROUTINE reduce_1_int(loc_var,erg,Op)
    INTEGER :: loc_var,erg
    INTEGER:: Op
    erg=loc_var
  END SUBROUTINE reduce_1_int

  SUBROUTINE reduce_1(loc_var,erg,Op)
    DOUBLE PRECISION:: loc_var,erg
    INTEGER:: Op
    erg=loc_var
  END SUBROUTINE reduce_1

#if FOR_LM == 1

  DOUBLE PRECISION FUNCTION global_maxval(feld)
    USE parallel_utilities, ONLY :   &
         & global_values               ! collects values from all nodes
    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE data_parallel , ONLY:        &
         & icomm_cart,               & ! communicator for the virtual cartesian topology
         & imp_reals                   ! determines the correct REAL type used in the model
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:,:),INTENT(in):: feld

    INTEGER (KIND=iintegers) :: ierror  ! error status variable
    CHARACTER (LEN=80) :: yerrmsg ! for error message
    REAL (KIND=ireals) :: loc_max

    loc_max = MAXVAL(feld)

    IF (num_compute > 1) THEN
      CALL global_values (loc_max,1,'MAX',imp_reals,icomm_cart,-1,yerrmsg,ierror)
    END IF

    global_maxval = loc_max

  END FUNCTION global_maxval

  DOUBLE PRECISION FUNCTION global_minval(feld)
    USE parallel_utilities, ONLY :   &
         & global_values               ! collects values from all nodes
    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE data_parallel , ONLY:        &
         & icomm_cart,               & ! communicator for the virtual cartesian topology
         & imp_reals                   ! determines the correct REAL type used in the model
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:,:),INTENT(in):: feld

    INTEGER (KIND=iintegers) :: ierror  ! error status variable
    CHARACTER (LEN=80) :: yerrmsg ! for error message
    REAL (KIND=ireals) :: loc_min

    loc_min = MINVAL(feld)

    IF (num_compute > 1) THEN
      CALL global_values (loc_min,1,'MIN',imp_reals,icomm_cart,-1,yerrmsg,ierror)
    END IF

    global_minval = loc_min

  END FUNCTION global_minval

  DOUBLE PRECISION FUNCTION global_maxval_2d(feld)
    USE parallel_utilities, ONLY :   &
         & global_values               ! collects values from all nodes
    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE data_parallel , ONLY:        &
         & icomm_cart,               & ! communicator for the virtual cartesian topology
         & imp_reals                   ! determines the correct REAL type used in the model
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:),INTENT(in):: feld

    INTEGER (KIND=iintegers) :: ierror  ! error status variable
    CHARACTER (LEN=80) :: yerrmsg ! for error message
    REAL (KIND=ireals) :: loc_max

    loc_max = MAXVAL(feld)

    IF (num_compute > 1) THEN
      CALL global_values (loc_max,1,'MAX',imp_reals,icomm_cart,-1,yerrmsg,ierror)
    END IF

    global_maxval_2d = loc_max

  END FUNCTION global_maxval_2d

  DOUBLE PRECISION FUNCTION global_minval_2d(feld)
    USE parallel_utilities, ONLY :   &
         & global_values               ! collects values from all nodes
    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE data_parallel , ONLY:        &
         & icomm_cart,               & ! communicator for the virtual cartesian topology
         & imp_reals                   ! determines the correct REAL type used in the model
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(:,:),INTENT(in):: feld

    INTEGER (KIND=iintegers) :: ierror  ! error status variable
    CHARACTER (LEN=80) :: yerrmsg ! for error message
    REAL (KIND=ireals) :: loc_min

    loc_min = MINVAL(feld)

    IF (num_compute > 1) THEN
      CALL global_values (loc_min,1,'MIN',imp_reals,icomm_cart,-1,yerrmsg,ierror)
    END IF

    global_minval_2d = loc_min

  END FUNCTION global_minval_2d

#else 

  FUNCTION global_maxval(feld)
    use mpi_interface,  only : myid, double_scalar_par_max
    IMPLICIT NONE

    double precision, DIMENSION(:,:,:),INTENT(in):: feld
    real(kind=8)  :: loc_max, global_maxval

    loc_max = dble(MAXVAL(feld))

    call double_scalar_par_max(loc_max,global_maxval)

  END FUNCTION global_maxval

  FUNCTION global_minval(feld)
    use mpi_interface,  only : myid, double_scalar_par_min
    IMPLICIT NONE

    double precision, DIMENSION(:,:,:),INTENT(in):: feld
    real(kind=8)  :: loc_min, global_minval

    loc_min = MINVAL(feld)

    call double_scalar_par_min(loc_min,global_minval)

  END FUNCTION global_minval

  FUNCTION double_global_maxval(feld)
    use mpi_interface,  only : myid, double_scalar_par_max
    IMPLICIT NONE

    real(kind=8), DIMENSION(:,:,:),INTENT(in):: feld
    real(kind=8) :: loc_max, double_global_maxval

    loc_max = MAXVAL(feld)

    call double_scalar_par_max(loc_max,double_global_maxval)

  END FUNCTION double_global_maxval

  SUBROUTINE global_maxval_stdout(text1,text2,feld)
    IMPLICIT NONE
    character(len=*), intent (in) :: text1,text2
    real(kind=8), DIMENSION(:,:,:),INTENT(in):: feld
    real(kind=8) :: maxi

    maxi = double_global_maxval(feld)

    if (isIO()) WRITE(*,'(3x,A,2x,A,E11.3)') text1,text2,maxi    

  END SUBROUTINE global_maxval_stdout 
#endif

END MODULE parallele_umgebung

MODULE initialisierung
  USE wolken_driver, ONLY: p_0, T_0, rho_0, dichte
END MODULE initialisierung

!=======================================================================

MODULE wolken_sedi
  !*******************************************************************************
  !                                                                              *
  !   Sedimentation der Wolken- und Niederschlagspartikel                        *
  !                                                                              *
  !*******************************************************************************

  USE wolken_konstanten 

  IMPLICIT NONE

CONTAINS

#if FOR_LM == 1

  SUBROUTINE rain_sedimentation_lm (qr,qnr,qc,rainrate,rho,rhocorr,adz,dt, &
      &             ims,ime,jms,jme,kms,kme,    & ! memory dims
      &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
        & iintegers,                & ! KIND-type parameters for integer variables
        & ireals                      ! KIND-type parameters for real variables

    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
        INTENT(INOUT) :: qr,qnr,qc
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme),    &
        INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme),    &
        INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    REAL(KIND=ireals)     :: n_r,q_r,x_r,n0,lam,v_n,v_q,n,G1,G4,c_n,c_q,s_n,s_q,D_m,D_r,mue

    REAL(KIND=ireals), PARAMETER :: alf = 9.65e+00     ! in SI [m/s]
    REAL(KIND=ireals), PARAMETER :: bet = 1.03e+01     ! in SI [m/s]
    REAL(KIND=ireals), PARAMETER :: gam = 6.00e+02     ! in SI [1/m]
    REAL(KIND=ireals), PARAMETER :: n_0 = 1.00e+07     ! in SI [1/m^4]
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-15  

    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme) :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_rain,v_q_rain,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax

    WHERE (qnr < 0.0) qnr = 0.0
    WHERE (qr < 0.0) qr = 0.0

    v_n_rain = 0.0_ireals
    v_q_rain = 0.0_ireals
    q_fluss  = 0.0_ireals
    n_fluss  = 0.0_ireals

    n  = 0.0
    !G1 = gfct(n+1) 
    !G4 = gfct(n+4) * pi/6.*rho_w    
    G1 = 1.0 
    G4 = pi *rho_w

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite

          IF (qr(i,j,k,nnew) > eps) THEN

          qr(i,j,k,nnew) = rho(i,j,k)*qr(i,j,k,nnew)
          qnr(i,j,k,nnew) = rho(i,j,k)*qnr(i,j,k,nnew)
          n_r = qnr(i,j,k,nnew) + eps        !..Anzahldichte in SI
          q_r = qr(i,j,k,nnew)  + eps        !..Fluessigwassergehalt in SI
          IF (cloud_typ <= 1) THEN
            lam = (pi*rho_w*n_0/q_r)**(1./4.)         !..Nur q_r  (n_0 = const) 
            D_r = 1.0/lam
            mue = 0.0
          ELSEIF (use_mu_Dm_rain_sedi) THEN
            x_r = q_r / n_r                           !..mittlere Masse in SI
            x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
            D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
            IF (qc(i,j,k,nnew) >= eps) THEN ! Seifert (2007)            
              ! UB_20080212>
              !              mue = 2.0
              mue = (rain%nu+1.0d0)/rain%b_geo - 1.0d0
              ! <UB_20080212
            ELSEIF (D_m.LE.rain_cmu3) THEN    
              mue = rain_cmu0*TANH((4.*rain_cmu2*(D_m-rain_cmu3))**2) &
                  & + rain_cmu4
            ELSE
              mue = rain_cmu1*TANH((1.*rain_cmu2*(D_m-rain_cmu3))**2) &
                  & + rain_cmu4
            ENDIF
            D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
          ELSEIF (use_mu_orig_rain_sedi) THEN
            x_r = q_r / n_r                           !..mittlere Masse in SI
            x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
            D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
            mue = (rain%nu+1.0d0)/rain%b_geo -1.0d0;
            D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
          ELSE
! UB_20090316>>
!!$            x_r = q_r / n_r                           !..mittlere Masse in SI
!!$            x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
!!$            !n0  = n_r * (pi*rho_w/x_r)**(1./3.)      !..n0   
!!$            !n0  = n_r/G1 * (G4/G1/x_r)**(1./3.)      !..n0
!!$            n0  = n_r * (G4/x_r)**(1./3.)             !..n0
!!$            n0  = min(max(250.0d+03,n0),20000.0d+03)
!!$            !lam = (pi*rho_w*n0/q_r)**(1./4.)         !..lambda von n0 und q
!!$            !lam = (G4*n0/q_r)**(1./(n+4))            !..lambda
!!$            lam = (G4*n0/q_r)**(0.25)                 !..lambda
!!$            lam = min(max(1d+03,lam),1d+04)
!!$            D_r = 1.0/lam
!!$            mue = 0.0
            ! use exponential DSD with limited N0 outside cloudwater present regions
            ! and original mue otherwise:
            IF (qc(i,j,k,nnew) >= q_krit) THEN            
              x_r = q_r / n_r                           !..mittlere Masse in SI
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
              D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
              mue = (rain%nu+1.0d0)/rain%b_geo - 1.0d0
              D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
            ELSE
              x_r = q_r / n_r                           !..mittlere Masse in SI
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
              n0  = n_r * (G4/x_r)**(1./3.)             !..n0
              n0  = MIN(MAX(250.0d+03,n0),20000.0d+03)
              lam = (G4*n0/q_r)**(0.25)                 !..lambda
! UB 20080223              lam = MIN(MAX(1d+03,lam),1d+04)
              D_r = 1.0/lam
              mue = 0.0
            ENDIF
! UB_20090316<<
          ENDIF
          v_n = alf - bet / (1.0 + gam*D_r)**(mue+1.)
          v_q = alf - bet / (1.0 + gam*D_r)**(mue+4.)
          v_n = v_n * rhocorr(i,j,k)
          v_q = v_q * rhocorr(i,j,k)
          v_n = MAX(v_n,1.d-1)
          v_q = MAX(v_q,1.d-1)
          v_n = MIN(v_n,2.d+1)
          v_q = MIN(v_q,2.d+1)
          v_n_rain(i,j,k) = - v_n
          v_q_rain(i,j,k) = - v_q
        ELSE
          v_n_rain(i,j,k) = 0.0_ireals
          v_q_rain(i,j,k) = 0.0_ireals
        ENDIF
        END DO
      END DO
    END DO

    v_n_rain(:,:,kte+1) = v_n_rain(:,:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_rain(:,:,kte+1) = v_q_rain(:,:,kte)   !
    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite

          v_n = 0.5 * (v_n_rain(i,j,k+1)+v_n_rain(i,j,k))
          v_q = 0.5 * (v_q_rain(i,j,k+1)+v_q_rain(i,j,k))
          ! Formulierung unter der Annahme, dass v_n, v_q stets negativ
          c_n = -v_n * adz(i,j,k) * dt 
          c_q = -v_q * adz(i,j,k) * dt 
          IF (c_n > 1) THEN
            kk = k
            s_n = 0.0
            DO WHILE (c_n > 1 .AND. kk.GT.2)
              s_n = s_n + qnr(i,j,kk,nnew)/adz(i,j,kk)
              c_n = (c_n - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_n = s_n + qnr(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_n,1.0d0)
            s_n = -s_n / dt
          ELSE
            s_n = v_n * qnr(i,j,k,nnew)
          ENDIF
          IF (c_q > 1) THEN
            kk = k
            s_q = 0.0
            DO WHILE (c_q > 1 .AND. kk.GT.2)
              s_q = s_q + qr(i,j,kk,nnew)/adz(i,j,kk)
              c_q = (c_q - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_q = s_q + qr(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_q,1.0d0)
            s_q = -s_q / dt
          ELSE
            s_q = v_q * qr(i,j,k,nnew)
          ENDIF

          ! Flux-limiter to avoid negative values
          n_fluss(i,j,k) = MAX(s_n,n_fluss(i,j,k-1)-qnr(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_q,q_fluss(i,j,k-1)-qr(i,j,k,nnew) /(adz(i,j,k)*dt))

        END DO
      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          qnr(i,j,k,nnew) = qnr(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qr(i,j,k,nnew)  = qr(i,j,k,nnew)  + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qr(i,j,k,nnew)  = qr(i,j,k,nnew)  / rho(i,j,k)
          qnr(i,j,k,nnew) = qnr(i,j,k,nnew) / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE rain_sedimentation_lm

  SUBROUTINE rain_sedi_lm_vec (qr,qnr,qc,rainrate,rho,rhocorr,adz,dt, &
      &             ims,ime,jms,jme,kms,kme,    & ! memory dims
      &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
        & iintegers,                & ! KIND-type parameters for integer variables
        & ireals                      ! KIND-type parameters for real variables

    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
        INTENT(INOUT) :: qr,qnr,qc
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme),    &
        INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme),    &
        INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    REAL(KIND=ireals)     :: n_r,q_r,x_r,n0,lam,n,G1,G4,D_m,D_r,mue,v_n,v_q

    REAL(KIND=ireals), PARAMETER :: alf = 9.65e+00     ! in SI [m/s]
    REAL(KIND=ireals), PARAMETER :: bet = 1.03e+01     ! in SI [m/s]
    REAL(KIND=ireals), PARAMETER :: gam = 6.00e+02     ! in SI [1/m]
    REAL(KIND=ireals), PARAMETER :: n_0 = 1.00e+07     ! in SI [1/m^4]
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-15  

!!! Maybe convert the following arrays to ALLOCATABLE to save stack memory!
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme) :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_rain,v_q_rain,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax
    REAL(KIND=ireals), DIMENSION(ims:ime) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL, DIMENSION(ims:ime) :: cflag

    WHERE (qnr < 0.0) qnr = 0.0
    WHERE (qr < 0.0) qr = 0.0

    v_n_rain = 0.0_ireals
    v_q_rain = 0.0_ireals
    q_fluss  = 0.0_ireals
    n_fluss  = 0.0_ireals

    n  = 0.0
    G1 = 1.0 
    G4 = pi *rho_w

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite

          IF (qr(i,j,k,nnew) > eps) THEN

          qr(i,j,k,nnew) = rho(i,j,k)*qr(i,j,k,nnew)
          qnr(i,j,k,nnew) = rho(i,j,k)*qnr(i,j,k,nnew)
          n_r = qnr(i,j,k,nnew) + eps        !..Anzahldichte in SI
          q_r = qr(i,j,k,nnew)  + eps        !..Fluessigwassergehalt in SI
          IF (cloud_typ <= 1) THEN
            lam = (pi*rho_w*n_0/q_r)**(1./4.)         !..Nur q_r  (n_0 = const) 
            D_r = 1.0/lam
            mue = 0.0
          ELSEIF (use_mu_Dm_rain_sedi) THEN
            x_r = q_r / n_r                           !..mittlere Masse in SI
            x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
            D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
            IF (qc(i,j,k,nnew) >= eps) THEN ! Seifert (2007)            
              ! UB_20080212>
              !              mue = 2.0
              mue = (rain%nu+1.0d0)/rain%b_geo - 1.0d0
              ! <UB_20080212
            ELSEIF (D_m.LE.rain_cmu3) THEN    
              mue = rain_cmu0*TANH((4.*rain_cmu2*(D_m-rain_cmu3))**2) &
                  & + rain_cmu4
            ELSE
              mue = rain_cmu1*TANH((1.*rain_cmu2*(D_m-rain_cmu3))**2) &
                  & + rain_cmu4
            ENDIF
            D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
          ELSEIF (use_mu_orig_rain_sedi) THEN
            x_r = q_r / n_r                           !..mittlere Masse in SI
            x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
            D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
            mue = (rain%nu+1.0d0)/rain%b_geo -1.0d0;
            D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
          ELSE
! UB_20090316>>
            ! use exponential DSD with limited N0 outside cloudwater present regions
            ! and original mue otherwise:
            IF (qc(i,j,k,nnew) >= q_krit) THEN            
              x_r = q_r / n_r                           !..mittlere Masse in SI
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
              D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
              mue = (rain%nu+1.0d0)/rain%b_geo - 1.0d0
              D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
            ELSE
              x_r = q_r / n_r                           !..mittlere Masse in SI
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
              n0  = n_r * (G4/x_r)**(1./3.)             !..n0
              n0  = MIN(MAX(250.0d+03,n0),20000.0d+03)
              lam = (G4*n0/q_r)**(0.25)                 !..lambda
! UB 20080223              lam = MIN(MAX(1d+03,lam),1d+04)
              D_r = 1.0/lam
              mue = 0.0
            ENDIF
! UB_20090316<<
          ENDIF
          v_n = alf - bet / (1.0 + gam*D_r)**(mue+1.)
          v_q = alf - bet / (1.0 + gam*D_r)**(mue+4.)
          v_n = v_n * rhocorr(i,j,k)
          v_q = v_q * rhocorr(i,j,k)
          v_n = MAX(v_n,1.d-1)
          v_q = MAX(v_q,1.d-1)
          v_n = MIN(v_n,2.d+1)
          v_q = MIN(v_q,2.d+1)
          v_n_rain(i,j,k) = - v_n
          v_q_rain(i,j,k) = - v_q
        ELSE
          v_n_rain(i,j,k) = 0.0_ireals
          v_q_rain(i,j,k) = 0.0_ireals
        ENDIF
        END DO
      END DO
    END DO

    v_n_rain(:,:,kte+1) = v_n_rain(:,:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_rain(:,:,kte+1) = v_q_rain(:,:,kte)   ! lower BC for terminal veloc.

    DO k = kts,kte
      DO j = jts,jte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_rain(i,j,k+1)+v_n_rain(i,j,k))
          v_qv(i) = 0.5 * (v_q_rain(i,j,k+1)+v_q_rain(i,j,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,j,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,j,k) * dt
        END DO

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * qnr(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + qnr(i,j,kk,nnew)/adz(i,j,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + qnr(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_nv(i),1.0d0)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qr(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qr(i,j,kk,nnew)/adz(i,j,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qr(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_qv(i),1.0d0)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,j,k) = MAX(s_nv(i),n_fluss(i,j,k-1)-qnr(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_qv(i),q_fluss(i,j,k-1)-qr(i,j,k,nnew) /(adz(i,j,k)*dt))
        END DO

      END DO
    END DO

    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0 ! upper BC condition

    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          qnr(i,j,k,nnew) = qnr(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qr(i,j,k,nnew)  = qr(i,j,k,nnew)  + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qr(i,j,k,nnew)  = qr(i,j,k,nnew)  / rho(i,j,k)
          qnr(i,j,k,nnew) = qnr(i,j,k,nnew) / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE rain_sedi_lm_vec

  SUBROUTINE graupel_sedimentation_lm (qg,ng,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Graupel                               *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qg,ng
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_g,q_g,x_g,n0,lam,v_n,v_q,n,c_n,c_q,s_n,s_q
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax

    !IF (isIO()) WRITE (6,*) 'graupel_sedimentation_lm: start'

    IF (firstcall.NE.1) THEN
      alf_n = graupel%a_vel * gfct((graupel%nu+graupel%b_vel+1.0)/graupel%mu)&
           &                / gfct((graupel%nu+1.0)/graupel%mu)
      alf_q = graupel%a_vel * gfct((graupel%nu+graupel%b_vel+2.0)/graupel%mu)&
           &                / gfct((graupel%nu+2.0)/graupel%mu)
      c_lam = gfct((graupel%nu+1.0)/graupel%mu)/gfct((graupel%nu+2.0)/graupel%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD graupel_sedi:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qg < 0.0d0) qg = 0.0d0
    WHERE (ng < 0.0d0) ng = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qg(i,j,k,nnew) = rho(i,j,k)*qg(i,j,k,nnew)
          ng(i,j,k,nnew) = rho(i,j,k)*ng(i,j,k,nnew)
          x_g = qg(i,j,k,nnew) / (ng(i,j,k,nnew)+eps)
          x_g = MIN(MAX(x_g,graupel%x_min),graupel%x_max)
          lam = ( c_lam * x_g )**(graupel%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.1d+0)
          v_q = MAX(v_q,0.1d+0)
          v_n = MIN(v_n,30.d0)
          v_q = MIN(v_q,30.d0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO
    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! unter Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   !
! UB_20080313> changed order of k and j loops
    DO k = kts,kte
      DO j = jts,jte
! UB_20080313<
        DO i = its,ite
          v_n = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_q = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_n, v_q stets negativ!!!
          c_n = -v_n * adz(i,j,k) * dt 
          c_q = -v_q * adz(i,j,k) * dt 
          IF (c_n > 1) THEN
            kk = k
            s_n = 0.0
            DO WHILE (c_n > 1 .AND. kk.GT.2)
              s_n = s_n + ng(i,j,kk,nnew)/adz(i,j,kk)
              c_n = (c_n - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_n = s_n + ng(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_n,1.0d0)
            s_n = -s_n / dt
          ELSE
            s_n = v_n * ng(i,j,k,nnew)
          ENDIF
          IF (c_q > 1) THEN
            kk = k
            s_q = 0.0
            DO WHILE (c_q > 1 .AND. kk.GT.2)
              s_q = s_q + qg(i,j,kk,nnew)/adz(i,j,kk)
              c_q = (c_q - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_q = s_q + qg(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_q,1.0d0)
            s_q = -s_q / dt
          ELSE
            s_q = v_q * qg(i,j,k,nnew)
          ENDIF
          ! Flux-limiter to avoid negative values
          n_fluss(i,j,k) = MAX(s_n,n_fluss(i,j,k-1)-ng(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_q,q_fluss(i,j,k-1)-qg(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO
      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          ng(i,j,k,nnew) = ng(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qg(i,j,k,nnew) = qg(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qg(i,j,k,nnew) = qg(i,j,k,nnew)  / rho(i,j,k)
          ng(i,j,k,nnew) = ng(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE graupel_sedimentation_lm

  SUBROUTINE graupel_sedi_lm_vec (qg,ng,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Graupel                               *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qg,ng
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_g,q_g,x_g,n0,lam,n,v_n,v_q
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

!!! Maybe convert the following arrays to ALLOCATABLE to save stack memory!
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax
    REAL(KIND=ireals), DIMENSION(ims:ime) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL, DIMENSION(ims:ime) :: cflag

    IF (firstcall.NE.1) THEN
      alf_n = graupel%a_vel * gfct((graupel%nu+graupel%b_vel+1.0)/graupel%mu)&
           &                / gfct((graupel%nu+1.0)/graupel%mu)
      alf_q = graupel%a_vel * gfct((graupel%nu+graupel%b_vel+2.0)/graupel%mu)&
           &                / gfct((graupel%nu+2.0)/graupel%mu)
      c_lam = gfct((graupel%nu+1.0)/graupel%mu)/gfct((graupel%nu+2.0)/graupel%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD graupel_sedi_lm_vec:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qg < 0.0d0) qg = 0.0d0
    WHERE (ng < 0.0d0) ng = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qg(i,j,k,nnew) = rho(i,j,k)*qg(i,j,k,nnew)
          ng(i,j,k,nnew) = rho(i,j,k)*ng(i,j,k,nnew)
          x_g = qg(i,j,k,nnew) / (ng(i,j,k,nnew)+eps)
          x_g = MIN(MAX(x_g,graupel%x_min),graupel%x_max)
          lam = ( c_lam * x_g )**(graupel%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.1d+0)
          v_q = MAX(v_q,0.1d+0)
          v_n = MIN(v_n,30.d0)
          v_q = MIN(v_q,30.d0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO

    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   ! lower BC for the terminal veloc.

    DO k = kts,kte
      DO j = jts,jte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_qv(i) = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,j,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,j,k) * dt
        END DO

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * ng(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + ng(i,j,kk,nnew)/adz(i,j,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + ng(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_nv(i),1.0d0)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qg(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qg(i,j,kk,nnew)/adz(i,j,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qg(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_qv(i),1.0d0)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,j,k) = MAX(s_nv(i),n_fluss(i,j,k-1)-ng(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_qv(i),q_fluss(i,j,k-1)-qg(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO

      END DO
    END DO

    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0 ! upper BC

    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          ng(i,j,k,nnew) = ng(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qg(i,j,k,nnew) = qg(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qg(i,j,k,nnew) = qg(i,j,k,nnew)  / rho(i,j,k)
          ng(i,j,k,nnew) = ng(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE graupel_sedi_lm_vec

  SUBROUTINE hail_sedimentation_lm (qh,nh,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Hail                               *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qh,nh
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_g,q_g,x_g,n0,lam,v_n,v_q,n,c_n,c_q,s_n,s_q
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax

    !IF (isIO()) WRITE (6,*) 'hail_sedimentation_lm: start'

    IF (firstcall.NE.1) THEN
      alf_n = hail%a_vel * gfct((hail%nu+hail%b_vel+1.0)/hail%mu)&
           &                / gfct((hail%nu+1.0)/hail%mu)
      alf_q = hail%a_vel * gfct((hail%nu+hail%b_vel+2.0)/hail%mu)&
           &                / gfct((hail%nu+2.0)/hail%mu)
      c_lam = gfct((hail%nu+1.0)/hail%mu)/gfct((hail%nu+2.0)/hail%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD hail_sedi:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qh < 0.0d0) qh = 0.0d0
    WHERE (nh < 0.0d0) nh = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qh(i,j,k,nnew) = rho(i,j,k)*qh(i,j,k,nnew)
          nh(i,j,k,nnew) = rho(i,j,k)*nh(i,j,k,nnew)
          x_g = qh(i,j,k,nnew) / (nh(i,j,k,nnew)+eps)
          x_g = MIN(MAX(x_g,hail%x_min),hail%x_max)
          lam = ( c_lam * x_g )**(hail%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.1d+0)
          v_q = MAX(v_q,0.1d+0)
          v_n = MIN(v_n,30.d0)
          v_q = MIN(v_q,30.d0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO
    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! unter Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   !
! UB_20080313> changed order of k and j loops
    DO k = kts,kte
      DO j = jts,jte
! UB_20080313<
        DO i = its,ite
          v_n = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_q = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_n, v_q stets negativ
          c_n = -v_n * adz(i,j,k) * dt 
          c_q = -v_q * adz(i,j,k) * dt 
          IF (c_n > 1) THEN
            kk = k
            s_n = 0.0
            DO WHILE (c_n > 1 .AND. kk.GT.2)
              s_n = s_n + nh(i,j,kk,nnew)/adz(i,j,kk)
              c_n = (c_n - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_n = s_n + nh(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_n,1.0d0)
            s_n = -s_n / dt
          ELSE
            s_n = v_n * nh(i,j,k,nnew)
          ENDIF
          IF (c_q > 1) THEN
            kk = k
            s_q = 0.0
            DO WHILE (c_q > 1 .AND. kk.GT.2)
              s_q = s_q + qh(i,j,kk,nnew)/adz(i,j,kk)
              c_q = (c_q - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_q = s_q + qh(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_q,1.0d0)
            s_q = -s_q / dt
          ELSE
            s_q = v_q * qh(i,j,k,nnew)
          ENDIF
          ! Flux-limiter to avoid negative values
          n_fluss(i,j,k) = MAX(s_n,n_fluss(i,j,k-1)-nh(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_q,q_fluss(i,j,k-1)-qh(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO
      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          nh(i,j,k,nnew) = nh(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qh(i,j,k,nnew) = qh(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qh(i,j,k,nnew) = qh(i,j,k,nnew)  / rho(i,j,k)
          nh(i,j,k,nnew) = nh(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE hail_sedimentation_lm

  SUBROUTINE hail_sedi_lm_vec (qh,nh,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Hail                               *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qh,nh
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_g,q_g,x_g,n0,lam,v_n,v_q,n
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

!!! Maybe convert the following arrays to ALLOCATABLE to save stack memory!
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax
    REAL(KIND=ireals), DIMENSION(ims:ime) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL, DIMENSION(ims:ime) :: cflag

    IF (firstcall.NE.1) THEN
      alf_n = hail%a_vel * gfct((hail%nu+hail%b_vel+1.0)/hail%mu)&
           &                / gfct((hail%nu+1.0)/hail%mu)
      alf_q = hail%a_vel * gfct((hail%nu+hail%b_vel+2.0)/hail%mu)&
           &                / gfct((hail%nu+2.0)/hail%mu)
      c_lam = gfct((hail%nu+1.0)/hail%mu)/gfct((hail%nu+2.0)/hail%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD hail_sedi_lm_vec:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qh < 0.0d0) qh = 0.0d0
    WHERE (nh < 0.0d0) nh = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qh(i,j,k,nnew) = rho(i,j,k)*qh(i,j,k,nnew)
          nh(i,j,k,nnew) = rho(i,j,k)*nh(i,j,k,nnew)
          x_g = qh(i,j,k,nnew) / (nh(i,j,k,nnew)+eps)
          x_g = MIN(MAX(x_g,hail%x_min),hail%x_max)
          lam = ( c_lam * x_g )**(hail%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.1d+0)
          v_q = MAX(v_q,0.1d+0)
          v_n = MIN(v_n,30.d0)
          v_q = MIN(v_q,30.d0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO
    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! unter Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   ! upper BC for terminal veloc.

    DO k = kts,kte
      DO j = jts,jte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_qv(i) = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,j,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,j,k) * dt
        END DO

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * nh(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + nh(i,j,kk,nnew)/adz(i,j,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + nh(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_nv(i),1.0d0)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qh(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qh(i,j,kk,nnew)/adz(i,j,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qh(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_qv(i),1.0d0)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,j,k) = MAX(s_nv(i),n_fluss(i,j,k-1)-nh(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_qv(i),q_fluss(i,j,k-1)-qh(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO

      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0 ! upper BC
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          nh(i,j,k,nnew) = nh(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qh(i,j,k,nnew) = qh(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qh(i,j,k,nnew) = qh(i,j,k,nnew)  / rho(i,j,k)
          nh(i,j,k,nnew) = nh(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE hail_sedi_lm_vec

  SUBROUTINE snow_sedimentation_lm (qs,ns,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Snow                               *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qs,ns
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_s,q_s,x_s,n0,lam,v_n,v_q,n,c_n,c_q,s_n,s_q,mue,nue,D_m
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax

    !IF (isIO()) WRITE (6,*) 'snow_sedimentation_lm: start'

    IF (firstcall.NE.1) THEN
      alf_n = snow%a_vel * gfct((snow%nu+snow%b_vel+1.0)/snow%mu) / gfct((snow%nu+1.0)/snow%mu)
      alf_q = snow%a_vel * gfct((snow%nu+snow%b_vel+2.0)/snow%mu) / gfct((snow%nu+2.0)/snow%mu)
      c_lam = gfct((snow%nu+1.0)/snow%mu)/gfct((snow%nu+2.0)/snow%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD snow_sedi:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qs < 0.0d0) qs = 0.0d0
    WHERE (ns < 0.0d0) ns = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qs(i,j,k,nnew) = rho(i,j,k)*qs(i,j,k,nnew)
          ns(i,j,k,nnew) = rho(i,j,k)*ns(i,j,k,nnew)
          x_s = qs(i,j,k,nnew) / (ns(i,j,k,nnew)+eps)
          x_s = MIN(MAX(x_s,snow%x_min),snow%x_max)
          IF (use_mu_Dm_snow_sedi) THEN
            D_m = snow%a_geo * x_s**snow%b_geo 
            mue = snow_cmu1*TANH((snow_cmu2*(D_m-snow_cmu3))**snow_cmu5) + snow_cmu4
            nue = (mue-2.0)/3.0
            alf_n = snow%a_vel * gfct((nue+snow%b_vel+1.0)/snow%mu)/gfct((nue+1.0)/snow%mu)
            alf_q = snow%a_vel * gfct((nue+snow%b_vel+2.0)/snow%mu)/gfct((nue+2.0)/snow%mu)
            c_lam = gfct((nue+1.0)/snow%mu)/gfct((nue+2.0)/snow%mu)
          ENDIF
          lam = ( c_lam * x_s )**(snow%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.1d+0)
          v_q = MAX(v_q,0.1d+0)
          v_n = MIN(v_n,3.0d+0)
          v_q = MIN(v_q,3.0d+0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO
    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! unter Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   !
! UB_20080313> changed order of k and j loops
    DO k = kts,kte
      DO j = jts,jte
! UB_20080313<
        DO i = its,ite
          v_n = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_q = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_n, v_q stets negativ
          c_n = -v_n * adz(i,j,k) * dt 
          c_q = -v_q * adz(i,j,k) * dt 
          IF (c_n > 1) THEN
            kk = k
            s_n = 0.0
            DO WHILE (c_n > 1 .AND. kk.GT.2)
              s_n = s_n + ns(i,j,kk,nnew)/adz(i,j,kk)
              c_n = (c_n - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_n = s_n + ns(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_n,1.0d0)
            s_n = -s_n / dt
          ELSE
            s_n = v_n * ns(i,j,k,nnew)
          ENDIF
          IF (c_q > 1) THEN
            kk = k
            s_q = 0.0
            DO WHILE (c_q > 1 .AND. kk.GT.2)
              s_q = s_q + qs(i,j,kk,nnew)/adz(i,j,kk)
              c_q = (c_q - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_q = s_q + qs(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_q,1.0d0)
            s_q = -s_q / dt
          ELSE
            s_q = v_q * qs(i,j,k,nnew)
          ENDIF
          ! Flux-limiter to avoid negative values
          n_fluss(i,j,k) = MAX(s_n,n_fluss(i,j,k-1)-ns(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_q,q_fluss(i,j,k-1)-qs(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO
      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          ns(i,j,k,nnew) = ns(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qs(i,j,k,nnew) = qs(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qs(i,j,k,nnew) = qs(i,j,k,nnew)  / rho(i,j,k)
          ns(i,j,k,nnew) = ns(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE snow_sedimentation_lm

  SUBROUTINE snow_sedi_lm_vec (qs,ns,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Snow                               *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qs,ns
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_s,q_s,x_s,n0,lam,v_n,v_q,n,mue,nue,D_m
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

!!! Maybe convert the following arrays to ALLOCATABLE to save stack memory!
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax
    REAL(KIND=ireals), DIMENSION(ims:ime) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL, DIMENSION(ims:ime) :: cflag

    IF (firstcall.NE.1) THEN
      alf_n = snow%a_vel * gfct((snow%nu+snow%b_vel+1.0)/snow%mu) / gfct((snow%nu+1.0)/snow%mu)
      alf_q = snow%a_vel * gfct((snow%nu+snow%b_vel+2.0)/snow%mu) / gfct((snow%nu+2.0)/snow%mu)
      c_lam = gfct((snow%nu+1.0)/snow%mu)/gfct((snow%nu+2.0)/snow%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD snow_sedi_lm_vec:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qs < 0.0d0) qs = 0.0d0
    WHERE (ns < 0.0d0) ns = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qs(i,j,k,nnew) = rho(i,j,k)*qs(i,j,k,nnew)
          ns(i,j,k,nnew) = rho(i,j,k)*ns(i,j,k,nnew)
          x_s = qs(i,j,k,nnew) / (ns(i,j,k,nnew)+eps)
          x_s = MIN(MAX(x_s,snow%x_min),snow%x_max)
          IF (use_mu_Dm_snow_sedi) THEN
            D_m = snow%a_geo * x_s**snow%b_geo 
            mue = snow_cmu1*TANH((snow_cmu2*(D_m-snow_cmu3))**snow_cmu5) + snow_cmu4
            nue = (mue-2.0)/3.0
            alf_n = snow%a_vel * gfct((nue+snow%b_vel+1.0)/snow%mu)/gfct((nue+1.0)/snow%mu)
            alf_q = snow%a_vel * gfct((nue+snow%b_vel+2.0)/snow%mu)/gfct((nue+2.0)/snow%mu)
            c_lam = gfct((nue+1.0)/snow%mu)/gfct((nue+2.0)/snow%mu)
          ENDIF
          lam = ( c_lam * x_s )**(snow%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.1d+0)
          v_q = MAX(v_q,0.1d+0)
          v_n = MIN(v_n,3.0d+0)
          v_q = MIN(v_q,3.0d+0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO

    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   ! lower BC for sedimentation veloc.

    DO k = kts,kte
      DO j = jts,jte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_qv(i) = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,j,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,j,k) * dt
        END DO

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * ns(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + ns(i,j,kk,nnew)/adz(i,j,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + ns(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_nv(i),1.0d0)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qs(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qs(i,j,kk,nnew)/adz(i,j,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qs(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_qv(i),1.0d0)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,j,k) = MAX(s_nv(i),n_fluss(i,j,k-1)-ns(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_qv(i),q_fluss(i,j,k-1)-qs(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO

      END DO
    END DO

    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0 ! upper BC

    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          ns(i,j,k,nnew) = ns(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qs(i,j,k,nnew) = qs(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qs(i,j,k,nnew) = qs(i,j,k,nnew)  / rho(i,j,k)
          ns(i,j,k,nnew) = ns(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE snow_sedi_lm_vec

  SUBROUTINE ice_sedimentation_lm (qi,ni,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Ice                                   *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qi,ni
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_i,q_i,x_i,n0,lam,v_n,v_q,n,c_n,c_q,s_n,s_q
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax

    !IF (isIO()) WRITE (6,*) 'ice_sedimentation_lm: start'

    IF (firstcall.NE.1) THEN
      alf_n = ice%a_vel * gfct((ice%nu+ice%b_vel+1.0)/ice%mu) / gfct((ice%nu+1.0)/ice%mu)
      alf_q = ice%a_vel * gfct((ice%nu+ice%b_vel+2.0)/ice%mu) / gfct((ice%nu+2.0)/ice%mu)
      c_lam = gfct((ice%nu+1.0)/ice%mu)/gfct((ice%nu+2.0)/ice%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD ice_sedi:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qi < 0.0d0) qi = 0.0d0
    WHERE (ni < 0.0d0) ni = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qi(i,j,k,nnew) = rho(i,j,k)*qi(i,j,k,nnew)
          ni(i,j,k,nnew) = rho(i,j,k)*ni(i,j,k,nnew)
          x_i = qi(i,j,k,nnew) / (ni(i,j,k,nnew)+eps)
          x_i = MIN(MAX(x_i,ice%x_min),ice%x_max)
          lam = ( c_lam * x_i )**(ice%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.0d+0)
          v_q = MAX(v_q,0.0d+0)
          v_n = MIN(v_n,3.0d+0)
          v_q = MIN(v_q,3.0d+0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO
    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   !
! UB_20080313> changed order of k and j loops
    DO k = kts,kte
      DO j = jts,jte
! UB_20080313<
        DO i = its,ite
          v_n = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_q = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_n, v_q stets negativ
          c_n = -v_n * adz(i,j,k) * dt 
          c_q = -v_q * adz(i,j,k) * dt 
          IF (c_n > 1) THEN
            kk = k
            s_n = 0.0
            DO WHILE (c_n > 1 .AND. kk.GT.2)
              s_n = s_n + ni(i,j,kk,nnew)/adz(i,j,kk)
              c_n = (c_n - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_n = s_n + ni(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_n,1.0d0)
            s_n = -s_n / dt
          ELSE
            s_n = v_n * ni(i,j,k,nnew)
          ENDIF
          IF (c_q > 1) THEN
            kk = k
            s_q = 0.0
            DO WHILE (c_q > 1 .AND. kk.GT.2)
              s_q = s_q + qi(i,j,kk,nnew)/adz(i,j,kk)
              c_q = (c_q - 1) * adz(i,j,kk-1)/adz(i,j,kk) 
              kk  = kk - 1
            ENDDO
            s_q = s_q + qi(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_q,1.0d0)
            s_q = -s_q / dt
          ELSE
            s_q = v_q * qi(i,j,k,nnew)
          ENDIF
          ! Flux-limiter to avoid negative values
          n_fluss(i,j,k) = MAX(s_n,n_fluss(i,j,k-1)-ni(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_q,q_fluss(i,j,k-1)-qi(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO
      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          ni(i,j,k,nnew) = ni(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qi(i,j,k,nnew) = qi(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qi(i,j,k,nnew) = qi(i,j,k,nnew)  / rho(i,j,k)
          ni(i,j,k,nnew) = ni(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE ice_sedimentation_lm

  SUBROUTINE ice_sedi_lm_vec (qi,ni,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime,jms,jme,kms,kme,    & ! memory dims
       &             its,ite,jts,jte,kts,kte,nnew) ! dims
    !*******************************************************************************
    !       Berechnung der Sedimentation von Ice                                   *
    !*******************************************************************************

    USE data_parameters , ONLY :     &
         & iintegers,                & ! KIND-type parameters for integer variables
         & ireals                      ! KIND-type parameters for real variables
    USE konstanten
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
    INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
    INTEGER(KIND=iintegers), INTENT(IN) :: nnew 
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme, 2 ),    &
         INTENT(INOUT) :: qi,ni
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(IN)    :: rho,adz,rhocorr
    REAL(KIND=ireals), DIMENSION( ims:ime , jms:jme , kms:kme ),    &
         INTENT(INOUT) :: rainrate
    REAL(KIND=ireals) :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL(KIND=ireals)     :: n_i,q_i,x_i,n0,lam,v_n,v_q,n
    REAL(KIND=ireals), SAVE      :: alf_n,alf_q,c_lam
    REAL(KIND=ireals), PARAMETER :: eps = 1.00e-20   

!!! Maybe convert the following arrays to ALLOCATABLE to save stack memory!
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme)   :: q_fluss,n_fluss
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme,kms-1:kme+1) :: v_n_sedi,v_q_sedi,zr
    REAL(KIND=ireals), DIMENSION(ims:ime,jms:jme) :: zmax
    REAL(KIND=ireals), DIMENSION(ims:ime) :: v_nv, v_qv, s_nv, s_qv, c_nv, c_qv
    LOGICAL, DIMENSION(ims:ime) :: cflag

    IF (firstcall.NE.1) THEN
      alf_n = ice%a_vel * gfct((ice%nu+ice%b_vel+1.0)/ice%mu) / gfct((ice%nu+1.0)/ice%mu)
      alf_q = ice%a_vel * gfct((ice%nu+ice%b_vel+2.0)/ice%mu) / gfct((ice%nu+2.0)/ice%mu)
      c_lam = gfct((ice%nu+1.0)/ice%mu)/gfct((ice%nu+2.0)/ice%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD ice_sedi_lm_vec:"
        WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qi < 0.0d0) qi = 0.0d0
    WHERE (ni < 0.0d0) ni = 0.0d0

    q_fluss = 0.0
    n_fluss = 0.0

    DO k = kts,kte
      DO j = jts,jte
        DO i = its,ite
          qi(i,j,k,nnew) = rho(i,j,k)*qi(i,j,k,nnew)
          ni(i,j,k,nnew) = rho(i,j,k)*ni(i,j,k,nnew)
          x_i = qi(i,j,k,nnew) / (ni(i,j,k,nnew)+eps)
          x_i = MIN(MAX(x_i,ice%x_min),ice%x_max)
          lam = ( c_lam * x_i )**(ice%b_vel) * rhocorr(i,j,k)
          v_n = alf_n * lam
          v_q = alf_q * lam
          v_n = MAX(v_n,0.0d+0)
          v_q = MAX(v_q,0.0d+0)
          v_n = MIN(v_n,3.0d+0)
          v_q = MIN(v_q,3.0d+0)
          v_n_sedi(i,j,k) = -v_n
          v_q_sedi(i,j,k) = -v_q
        END DO
      END DO
    END DO

    v_n_sedi(:,:,kte+1) = v_n_sedi(:,:,kte)   ! untere Randbedingung fuer Fallgeschw.
    v_q_sedi(:,:,kte+1) = v_q_sedi(:,:,kte)   ! lower BC for sedimentation veloc.

    DO k = kts,kte
      DO j = jts,jte

        DO i = its,ite
          v_nv(i) = 0.5 * (v_n_sedi(i,j,k+1)+v_n_sedi(i,j,k))
          v_qv(i) = 0.5 * (v_q_sedi(i,j,k+1)+v_q_sedi(i,j,k))
          ! Formulierung unter der Annahme, dass v_nv, v_qv stets negativ
          c_nv(i) = -v_nv(i) * adz(i,j,k) * dt 
          c_qv(i) = -v_qv(i) * adz(i,j,k) * dt
        END DO

        kk = k
        s_nv = 0.0
        DO i = its,ite
          IF (c_nv(i) <= 1) THEN
            s_nv(i) = v_nv(i) * ni(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_nv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_nv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_nv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_nv(i) = s_nv(i) + ni(i,j,kk,nnew)/adz(i,j,kk)
                c_nv(i) = (c_nv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_nv(i) = s_nv(i) + ni(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_nv(i),1.0d0)
              s_nv(i) = -s_nv(i) / dt
            END IF
          END DO
        ENDIF

        kk = k
        s_qv = 0.0
        DO i = its,ite
          IF (c_qv(i) <= 1) THEN
            s_qv(i) = v_qv(i) * qi(i,j,k,nnew)
          END IF
        END DO
        IF (ANY(c_qv(its:ite) > 1)) THEN
          cflag = .FALSE.
          DO WHILE (ANY(c_qv(its:ite) > 1) .AND. kk > 2)
            DO i = its,ite
              IF (c_qv(i) > 1) THEN
                cflag(i) = .TRUE.
                s_qv(i) = s_qv(i) + qi(i,j,kk,nnew)/adz(i,j,kk)
                c_qv(i) = (c_qv(i) - 1) * adz(i,j,kk-1)/adz(i,j,kk)
              END IF
            END DO
            kk  = kk - 1
          ENDDO
          DO i = its,ite
            IF (cflag(i)) THEN
              s_qv(i) = s_qv(i) + qi(i,j,kk,nnew)/adz(i,j,kk)*MIN(c_qv(i),1.0d0)
              s_qv(i) = -s_qv(i) / dt
            END IF
          END DO
        ENDIF

        ! Flux-limiter to avoid negative values
        DO i = its,ite
          n_fluss(i,j,k) = MAX(s_nv(i),n_fluss(i,j,k-1)-ni(i,j,k,nnew)/(adz(i,j,k)*dt))
          q_fluss(i,j,k) = MAX(s_qv(i),q_fluss(i,j,k-1)-qi(i,j,k,nnew)/(adz(i,j,k)*dt))
        END DO

      END DO
    END DO
    n_fluss(:,:,kts-1) = 0.0 ! obere Randbedingung
    q_fluss(:,:,kts-1) = 0.0
    DO k = kts,kte
      DO j = jts, jte
        DO i = its,ite
          ni(i,j,k,nnew) = ni(i,j,k,nnew) + ( n_fluss(i,j,k) - n_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qi(i,j,k,nnew) = qi(i,j,k,nnew) + ( q_fluss(i,j,k) - q_fluss(i,j,k-1) )*adz(i,j,k)*dt
          qi(i,j,k,nnew) = qi(i,j,k,nnew)  / rho(i,j,k)
          ni(i,j,k,nnew) = ni(i,j,k,nnew)  / rho(i,j,k)
          rainrate(i,j,k) = - q_fluss(i,j,k) ! Regenrate
        ENDDO
      ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE ice_sedi_lm_vec

#else
   ! FOR_LM == 0

  SUBROUTINE rain_sedimentation_wrf (qr,qnr,qc,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime, jms,jme, kms,kme,               & ! memory dims
       &             its,ite, jts,jte, kts,kte                ) ! tile   dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE konstanten

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme , &
                           its,ite, jts,jte, kts,kte 
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(INOUT) :: qr,qnr,qc,adz
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(IN)    :: rho,rhocorr
    REAL, DIMENSION( ims:ime , jms:jme ),    &
         INTENT(INOUT) :: rainrate
    REAL :: dt

    ! .. Local Variables ..
    INTEGER  :: i, j, k, kk
    REAL     :: n_r,q_r,x_r,n0,lam,v_n,v_q,n,G1,G4,s_n,s_q,D_r,D_m,mue
    REAL :: c_n,c_q

    REAL, PARAMETER :: alf = 9.65e+00     ! in SI [m/s]
    REAL            :: bet = 1.03e+01     ! in SI [m/s]
    REAL, PARAMETER :: gam = 6.00e+02     ! in SI [1/m]
    REAL, PARAMETER :: n_0 = 1.00e+07     ! in SI [1/m^4]
    REAL, PARAMETER :: eps = 1.00e-20   
    REAL, PARAMETER :: D_v = 25.e-6   

    REAL, DIMENSION(its:ite,kts:kte+1,jts:jte) :: q_fluss,n_fluss
    REAL, DIMENSION(its:ite,kts-1:kte,jts:jte) :: v_n_rain,v_q_rain

    WHERE (qnr < 0.0) qnr = 0.0
    WHERE (qr < 0.0) qr = 0.0

!    bet = alf*exp(gam*D_v)  ! Bjorn's idea to adjust bet with D_v

    v_n_rain = 0.0
    v_q_rain = 0.0
    q_fluss = 0.0
    n_fluss = 0.0

    G1 = 1.0 
    G4 = pi *rho_w

    DO j = jts,jte
       DO k = kts,kte
          DO i = its,ite

!             IF (qr(i,k,j) < eps) CYCLE

             qr(i,k,j) = rho(i,k,j)*qr(i,k,j)
             qnr(i,k,j) = rho(i,k,j)*qnr(i,k,j)

             n_r = qnr(i,k,j) + eps                    !..Anzahldichte in SI
             q_r = qr(i,k,j)  + eps          !..Fluessigwassergehalt in SI
             IF (cloud_typ <= 1) THEN
                lam = (pi*rho_w*n_0/q_r)**(1./4.)         !..Nur q_r  (n_0 = const) 
                D_r = 1.0/lam
                mue = 0.0
             ELSEIF (use_mu_Dm_rain_sedi) THEN
               x_r = q_r / n_r                           !..mittlere Masse in SI
               x_r = MIN(MAX(x_r,REAL(rain%x_min)),REAL(rain%x_max))
               D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.)                               
               IF (qc(i,k,j) >= eps) THEN ! Seifert (2007)            
! UB_20080212>
!                 mue = 2.0
                 mue = (rain%nu+1.0d0)/rain%b_geo - 1.0d0
! <UB_20080212
               ELSEIF (D_m.LE.rain_cmu3) THEN    
                 mue = rain_cmu0*TANH((4.*rain_cmu2*(D_m-rain_cmu3))**2) &
                     & + rain_cmu4
               ELSE
                 mue = rain_cmu1*TANH((1.*rain_cmu2*(D_m-rain_cmu3))**2) &
                     & + rain_cmu4
               ENDIF
               D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)

             ELSEIF (use_mu_orig_rain_sedi) THEN
               x_r = q_r / n_r                           !..mittlere Masse in SI
               x_r = MIN(MAX(x_r,REAL(rain%x_min)),REAL(rain%x_max))
               D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
               mue = (rain%nu+1.0d0)/rain%b_geo -1.0d0;
               D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)

             ELSE
!!$               ! use exponential DSD with limited N0 and lambda:
!!$               x_r = q_r / n_r                           !..mittlere Masse in SI
!!$               x_r = MIN(MAX(x_r,REAL(rain%x_min)),REAL(rain%x_max))
!!$               n0  = n_r * (G4/x_r)**(1./3.)             !..n0
!!$               n0  = MIN(MAX(250.0e+03,n0),20000.0e+03)
!!$               lam = (G4*n0/q_r)**(0.25)                 !..lambda
!!$               lam = MIN(MAX(1e+03,lam),1e+04)
!!$               D_r = 1.0 / lam
!!$               mue = 0.0
!!$               ! ub<<
               ! use exponential DSD with limited N0 outside cloudwater present regions:
               IF (qc(i,k,j) >= q_krit) THEN ! Seifert (2007)            
                 x_r = q_r / n_r                           !..mittlere Masse in SI
                 x_r = MIN(MAX(x_r,REAL(rain%x_min)),REAL(rain%x_max))
                 D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.) 
                 mue = (rain%nu+1.0d0)/rain%b_geo - 1.0d0
                 D_r = (D_m**3/((mue+3.)*(mue+2.)*(mue+1.)))**(1./3.)
               ELSE
                 x_r = q_r / n_r                           !..mittlere Masse in SI
                 x_r = MIN(MAX(x_r,REAL(rain%x_min)),REAL(rain%x_max))
                 n0  = n_r * (G4/x_r)**(1./3.)             !..n0
                 n0  = MIN(MAX(250.0e+03,n0),20000.0e+03)
                 lam = (G4*n0/q_r)**(0.25)                 !..lambda
! UB 20080223              lam = MIN(MAX(1d+03,lam),1d+04)
                 D_r = 1.0/lam
                 mue = 0.0
               ENDIF
             ENDIF
             v_n = alf - bet / (1.0 + gam*D_r)**(mue+1.)    !..Fallgeschw. fuer n, MP-Vert,
             v_q = alf - bet / (1.0 + gam*D_r)**(mue+4.)    !..Fallgeschw. fuer q, MP-Vert.
             v_n = v_n * rhocorr(i,k,j)
             v_q = v_q * rhocorr(i,k,j)
             v_n = MAX(v_n,1.e-1)
             v_q = MAX(v_q,1.e-1)
             v_n = MIN(v_n,2.e+1)
             v_q = MIN(v_q,2.e+1)
             v_n_rain(i,k,j) = v_n
             v_q_rain(i,k,j) = v_q
          END DO
       END DO
    END DO
    ! Unteren Wert verdoppeln wg. spaeterer Differenzberechnung (quasi untere Randbedingung)
    v_n_rain(:,kts-1,:) = v_n_rain(:,kts,:)
    v_q_rain(:,kts-1,:) = v_q_rain(:,kts,:)
    DO j = jts,jte
! UB_20080313> Bugfix: k-loop has to be from model top to bottom for the flux-limiter!
       DO k = kte,kts,-1
!       DO k = kts,kte
! UB_20080313<
          DO i = its,ite

            v_n = 0.5 * (v_n_rain(i,k-1,j)+v_n_rain(i,k,j))
            v_q = 0.5 * (v_q_rain(i,k-1,j)+v_q_rain(i,k,j))

            ! Formulierung unter der Annahme, dass v_n, v_q stets positiv
            c_n = v_n * adz(i,k,j) * dt 
            c_q = v_q * adz(i,k,j) * dt 
            IF (c_n > 1) THEN
              kk = k
              s_n = 0.0
              DO WHILE (c_n > 1 .AND. kk.LE.kte-1)
                s_n = s_n + qnr(i,kk,j)/adz(i,kk,j)
                c_n = (c_n - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_n = s_n + qnr(i,kk,j)/adz(i,kk,j)*MIN(c_n,1.0e0)
              s_n = s_n / dt
            ELSE
              s_n = v_n * qnr(i,k,j)
            ENDIF
            IF (c_q > 1) THEN
              kk = k
              s_q = 0.0
              DO WHILE (c_q > 1 .AND. kk.LE.kte-1)
                s_q = s_q + qr(i,kk,j)/adz(i,kk,j)
                c_q = (c_q - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_q = s_q + qr(i,kk,j)/adz(i,kk,j)*MIN(c_q,1.0e0)
              s_q = s_q / dt
            ELSE
              s_q = v_q * qr(i,k,j)
            ENDIF

            n_fluss(i,k,j) = - s_n
            q_fluss(i,k,j) = - s_q

            ! Flux-limiter to avoid negative values
            n_fluss(i,k,j) = MAX(n_fluss(i,k,j),n_fluss(i,k+1,j)-qnr(i,k,j)/(adz(i,k,j)*dt))
            q_fluss(i,k,j) = MAX(q_fluss(i,k,j),q_fluss(i,k+1,j)- qr(i,k,j)/(adz(i,k,j)*dt))

          END DO
       END DO
    END DO
    DO j = jts,jte
       DO i = its,ite
          n_fluss(i,kte+1,j) = 0.0
          q_fluss(i,kte+1,j) = 0.0
          rainrate(i,j) = -q_fluss(i,kts,j)
       END DO
    END DO
    DO j = jts, jte
       DO k = kts,kte
          DO i = its,ite
             qnr(i,k,j) = qnr(i,k,j) + ( n_fluss(i,k,j) - n_fluss(i,k+1,j) )*adz(i,k,j)*dt
             qnr(i,k,j)  = qnr(i,k,j)  / rho(i,k,j)
             qr(i,k,j)  = qr(i,k,j)  + ( q_fluss(i,k,j) - q_fluss(i,k+1,j) )*adz(i,k,j)*dt
             qr(i,k,j)  = qr(i,k,j)  / rho(i,k,j)
          ENDDO
       ENDDO
    ENDDO
    
    RETURN

  END SUBROUTINE rain_sedimentation_wrf


  SUBROUTINE graupel_sedimentation_wrf (qg,qng,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime, jms,jme, kms,kme,               & ! memory dims
       &             its,ite, jts,jte, kts,kte                ) ! tile   dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE konstanten

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme , &
                           its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(INOUT) :: qg,qng,adz
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(IN)    :: rho,rhocorr
    REAL, DIMENSION( ims:ime , jms:jme ),    &
         INTENT(INOUT) :: rainrate
    REAL :: dt

    ! .. Local Variables ..
    INTEGER          :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL :: x,n0,lam,v_n,v_q,n,s_n,s_q
    REAL, SAVE      :: alf_n,alf_q,c_lam
    REAL, PARAMETER :: eps = 1.00e-20

    REAL :: c_n,c_q    
    REAL, DIMENSION(its:ite,kts:kte+1,jts:jte) :: q_fluss,n_fluss
    REAL, DIMENSION(its:ite,kts-1:kte,jts:jte) :: v_n_sedi,v_q_sedi

    IF (firstcall.NE.1) THEN
       alf_n = graupel%a_vel * gfct((graupel%nu+graupel%b_vel+1.0)/graupel%mu)&
            &                / gfct((graupel%nu+1.0)/graupel%mu)
       alf_q = graupel%a_vel * gfct((graupel%nu+graupel%b_vel+2.0)/graupel%mu)&
            &                / gfct((graupel%nu+2.0)/graupel%mu)
       c_lam = gfct((graupel%nu+1.0)/graupel%mu)/gfct((graupel%nu+2.0)/graupel%mu)
       IF(isdebug) THEN
          WRITE (6,*)           "CLOUD graupel_sedi:"
          WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
          WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
          WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
       ENDIF
       firstcall = 1
    ENDIF

    WHERE (qg < 0.0) qg = 0.0
    WHERE (qng < 0.0) qng = 0.0

    q_fluss = 0.0
    n_fluss = 0.0

    DO j = jts,jte
       DO k = kts,kte
          DO i = its,ite
             qg(i,k,j) = rho(i,k,j)*qg(i,k,j)
             qng(i,k,j) = rho(i,k,j)*qng(i,k,j)
             x = qg(i,k,j) / (qng(i,k,j)+eps)
             x = MIN(MAX(x,REAL(graupel%x_min)),REAL(graupel%x_max))
             lam = ( c_lam * x )**(graupel%b_vel) * rhocorr(i,k,j)
             v_n = alf_n * lam
             v_q = alf_q * lam
             v_n = MAX(v_n,0.1e+0)
             v_q = MAX(v_q,0.1e+0)
             v_n = MIN(v_n,30.e0)
             v_q = MIN(v_q,30.e0)
             v_n_sedi(i,k,j) = v_n
             v_q_sedi(i,k,j) = v_q
          END DO
       END DO
    END DO
    v_n_sedi(:,kts-1,:) = v_n_sedi(:,kts,:)
    v_q_sedi(:,kts-1,:) = v_q_sedi(:,kts,:)

    DO j = jts,jte
! UB_20080313> Bugfix: k-loop has to be from model top to bottom for the flux-limiter!
       DO k = kte,kts,-1
!       DO k = kts,kte
! UB_20080313<
          DO i = its,ite

            v_n = 0.5 * (v_n_sedi(i,k-1,j)+v_n_sedi(i,k,j))
            v_q = 0.5 * (v_q_sedi(i,k-1,j)+v_q_sedi(i,k,j))

            ! Formulierung unter der Annahme, dass v_n, v_q stets positiv
            c_n = v_n * adz(i,k,j) * dt 
            c_q = v_q * adz(i,k,j) * dt 
            IF (c_n > 1) THEN
              kk = k
              s_n = 0.0
              DO WHILE (c_n > 1 .AND. kk.LE.kte-1)
                s_n = s_n + qng(i,kk,j)/adz(i,kk,j)
                c_n = (c_n - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_n = s_n + qng(i,kk,j)/adz(i,kk,j)*MIN(c_n,1.0e0)
              s_n = s_n / dt
            ELSE
              s_n = v_n * qng(i,k,j)
            ENDIF
            IF (c_q > 1) THEN
              kk = k
              s_q = 0.0
              DO WHILE (c_q > 1 .AND. kk.LE.kte-1)
                s_q = s_q + qg(i,kk,j)/adz(i,kk,j)
                c_q = (c_q - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_q = s_q + qg(i,kk,j)/adz(i,kk,j)*MIN(c_q,1.0e0)
              s_q = s_q / dt
            ELSE
              s_q = v_q * qg(i,k,j)
            ENDIF

            n_fluss(i,k,j) = - s_n
            q_fluss(i,k,j) = - s_q

            ! Flux-limiter to avoid negative values
            n_fluss(i,k,j) = MAX(n_fluss(i,k,j),n_fluss(i,k+1,j)-qng(i,k,j)/(adz(i,k,j)*dt))
            q_fluss(i,k,j) = MAX(q_fluss(i,k,j),q_fluss(i,k+1,j)- qg(i,k,j)/(adz(i,k,j)*dt))

          END DO
       END DO
    END DO
    DO j = jts,jte
      DO i = its,ite
        n_fluss(i,kte+1,j) = 0.0
        q_fluss(i,kte+1,j) = 0.0
        rainrate(i,j) = -q_fluss(i,kts,j)
      END DO
    END DO
    DO j = jts, jte
       DO k = kts,kte
          DO i = its,ite
             qng(i,k,j) = qng(i,k,j) + ( n_fluss(i,k,j) - n_fluss(i,k+1,j) )*adz(i,k,j)*dt
             qng(i,k,j) = qng(i,k,j) / rho(i,k,j)
             qg(i,k,j) = qg(i,k,j) + ( q_fluss(i,k,j) - q_fluss(i,k+1,j) )*adz(i,k,j)*dt
             qg(i,k,j) = qg(i,k,j) / rho(i,k,j)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE graupel_sedimentation_wrf

  SUBROUTINE hail_sedimentation_wrf (qh,qnh,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime, jms,jme, kms,kme,               & ! memory dims
       &             its,ite, jts,jte, kts,kte                ) ! tile   dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE konstanten

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme , &
                           its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(INOUT) :: qh,qnh,adz
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(IN)    :: rho,rhocorr
    REAL, DIMENSION( ims:ime , jms:jme ),    &
         INTENT(INOUT) :: rainrate
    REAL :: dt

    ! .. Local Variables ..
    INTEGER          :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL :: x,n0,lam,v_n,v_q,n,s_n,s_q,cc,zz
    REAL, SAVE      :: alf_n,alf_q,c_lam
    REAL, PARAMETER :: eps = 1.00e-20

    REAL :: c_n,c_q    
    REAL, DIMENSION(its:ite,kts:kte+1,jts:jte) :: q_fluss,n_fluss
    REAL, DIMENSION(its:ite,kts-1:kte,jts:jte) :: v_n_sedi,v_q_sedi

    IF (firstcall.NE.1) THEN
       alf_n = hail%a_vel * gfct((hail%nu+hail%b_vel+1.0)/hail%mu)&
            &                / gfct((hail%nu+1.0)/hail%mu)
       alf_q = hail%a_vel * gfct((hail%nu+hail%b_vel+2.0)/hail%mu)&
            &                / gfct((hail%nu+2.0)/hail%mu)
       c_lam = gfct((hail%nu+1.0)/hail%mu)/gfct((hail%nu+2.0)/hail%mu)
       IF(isdebug) THEN
          WRITE (6,*)           "CLOUD hail_sedi:"
          WRITE (6,'(A,D10.3)') "  Coeff. lambda:          c_lam = ",c_lam
          WRITE (6,'(A,D10.3)') "  Coeff. sedimentation n: alf_n = ",alf_n
          WRITE (6,'(A,D10.3)') "  Coeff. sedimentation q: alf_q = ",alf_q
       ENDIF
       firstcall = 1
    ENDIF

    WHERE (qh < 0.0) qh = 0.0
    WHERE (qnh < 0.0) qnh = 0.0

    q_fluss = 0.0
    n_fluss = 0.0

    DO j = jts,jte
       DO k = kts,kte
          DO i = its,ite
             qh(i,k,j) = rho(i,k,j)*qh(i,k,j)
             qnh(i,k,j) = rho(i,k,j)*qnh(i,k,j)
             x = qh(i,k,j) / (qnh(i,k,j)+eps)
             x = MIN(MAX(x,REAL(hail%x_min)),REAL(hail%x_max))
             lam = ( c_lam * x )**(hail%b_vel) * rhocorr(i,k,j)
             v_n = alf_n * lam
             v_q = alf_q * lam
             v_n = MAX(v_n,0.1e+0)
             v_q = MAX(v_q,0.1e+0)
             v_n = MIN(v_n,30.e0)
             v_q = MIN(v_q,30.e0)
             v_n_sedi(i,k,j) = v_n
             v_q_sedi(i,k,j) = v_q
          END DO
       END DO
    END DO
    v_n_sedi(:,kts-1,:) = v_n_sedi(:,kts,:)
    v_q_sedi(:,kts-1,:) = v_q_sedi(:,kts,:)
    DO j = jts,jte
! UB_20080313> Bugfix: k-loop has to be from model top to bottom for the flux-limiter!
       DO k = kte,kts,-1
!       DO k = kts,kte
! UB_20080313<
          DO i = its,ite

            v_n = 0.5 * (v_n_sedi(i,k-1,j)+v_n_sedi(i,k,j))
            v_q = 0.5 * (v_q_sedi(i,k-1,j)+v_q_sedi(i,k,j))
            ! Formulierung unter der Annahme, dass v_n, v_q stets positiv!!!
            c_n = v_n * adz(i,k,j) * dt 
            c_q = v_q * adz(i,k,j) * dt 
            IF (c_n > 1) THEN
              kk = k
              s_n = 0.0
              DO WHILE (c_n > 1 .AND. kk.LE.kte-1)
                s_n = s_n + qnh(i,kk,j)/adz(i,kk,j)
                c_n = (c_n - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_n = s_n + qnh(i,kk,j)/adz(i,kk,j)*MIN(c_n,1.0e0)
              s_n = s_n / dt
            ELSE
              s_n = v_n * qnh(i,k,j)
            ENDIF
            IF (c_q > 1) THEN
              kk = k
              s_q = 0.0
              DO WHILE (c_q > 1 .AND. kk.LE.kte-1)
                s_q = s_q + qh(i,kk,j)/adz(i,kk,j)
                c_q = (c_q - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_q = s_q + qh(i,kk,j)/adz(i,kk,j)*MIN(c_q,1.0e0)
              s_q = s_q / dt
            ELSE
              s_q = v_q * qh(i,k,j)
            ENDIF

            n_fluss(i,k,j) = - s_n
            q_fluss(i,k,j) = - s_q

            ! Flux-limiter to avoid negative values
            n_fluss(i,k,j) = MAX(n_fluss(i,k,j),n_fluss(i,k+1,j)-qnh(i,k,j)/(adz(i,k,j)*dt))
            q_fluss(i,k,j) = MAX(q_fluss(i,k,j),q_fluss(i,k+1,j)- qh(i,k,j)/(adz(i,k,j)*dt))

          END DO
       END DO
    END DO
    DO j = jts,jte
       DO i = its,ite
          n_fluss(i,kte+1,j) = 0.0
          q_fluss(i,kte+1,j) = 0.0
          rainrate(i,j) = -q_fluss(i,kts,j)
       END DO
    END DO
    DO j = jts, jte
       DO k = kts,kte
          DO i = its,ite
             qnh(i,k,j) = qnh(i,k,j) + ( n_fluss(i,k,j) - n_fluss(i,k+1,j) )*adz(i,k,j)*dt
             qnh(i,k,j) = qnh(i,k,j) / rho(i,k,j)
             qh(i,k,j) = qh(i,k,j) + ( q_fluss(i,k,j) - q_fluss(i,k+1,j) )*adz(i,k,j)*dt
             qh(i,k,j) = qh(i,k,j) / rho(i,k,j)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE hail_sedimentation_wrf

  SUBROUTINE snow_sedimentation_wrf (qs,qns,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime, jms,jme, kms,kme,               & ! memory dims
       &             its,ite, jts,jte, kts,kte                ) ! tile   dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE konstanten

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme , &
                           its,ite, jts,jte, kts,kte 
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(INOUT) :: qs,qns,adz
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(IN)    :: rho,rhocorr
    REAL, DIMENSION( ims:ime , jms:jme ),    &
         INTENT(INOUT) :: rainrate
    REAL :: dt

    ! .. Local Variables ..
    INTEGER          :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL :: x,n0,lam,v_n,v_q,n,mue,nue,D_m,s_n,s_q
    REAL, SAVE      :: alf_n,alf_q,c_lam   
    REAL, PARAMETER :: eps = 1.00e-20   

    REAL :: c_n, c_q
    REAL, DIMENSION(its:ite,kts:kte+1,jts:jte) :: q_fluss,n_fluss
    REAL, DIMENSION(its:ite,kts-1:kte,jts:jte) :: v_n_sedi,v_q_sedi

    IF (firstcall.NE.1) THEN
       alf_n = snow%a_vel * gfct((snow%nu+snow%b_vel+1.0)/snow%mu)/gfct((snow%nu+1.0)/snow%mu)
       alf_q = snow%a_vel * gfct((snow%nu+snow%b_vel+2.0)/snow%mu)/gfct((snow%nu+2.0)/snow%mu)
       c_lam = gfct((snow%nu+1.0)/snow%mu)/gfct((snow%nu+2.0)/snow%mu)
       IF(isdebug) THEN
          WRITE (6,*)           "CLOUD snow_sedi:"
          WRITE (6,'(A,D10.3)') "  Koeff. fuer lambda:          c_lam = ",c_lam
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation n: alf_n = ",alf_n
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation q: alf_q = ",alf_q
       ENDIF
       firstcall = 1
    ENDIF

    WHERE (qs < 0.0) qs = 0.0
    WHERE (qns < 0.0) qns = 0.0

    q_fluss = 0.0
    n_fluss = 0.0

    DO j = jts,jte
       DO k = kts,kte
          DO i = its,ite
            qs(i,k,j) = rho(i,k,j)*qs(i,k,j)
            qns(i,k,j) = rho(i,k,j)*qns(i,k,j)
            x = qs(i,k,j) / (qns(i,k,j)+eps)
            x = MIN(MAX(x,REAL(snow%x_min)),REAL(snow%x_max))
            IF (use_mu_Dm_snow_sedi) THEN
              D_m = snow%a_geo * x**snow%b_geo
              mue = snow_cmu1*TANH((snow_cmu2*(D_m-snow_cmu3))**snow_cmu5) + snow_cmu4
              nue = (mue-2.0)/3.0
              alf_n = snow%a_vel * gfct((nue+snow%b_vel+1.0)/snow%mu)/gfct((nue+1.0)/snow%mu)
              alf_q = snow%a_vel * gfct((nue+snow%b_vel+2.0)/snow%mu)/gfct((nue+2.0)/snow%mu)
              c_lam = gfct((nue+1.0)/snow%mu)/gfct((nue+2.0)/snow%mu)
            ENDIF
            lam = ( c_lam * x )**(snow%b_vel) * rhocorr(i,k,j)
            v_n = alf_n * lam                                 !..Fallgeschw. fuer n
            v_q = alf_q * lam                                 !..Fallgeschw. fuer q
            v_n = MAX(v_n,0.1e+0)
            v_q = MAX(v_q,0.1e+0)
            v_n = MIN(v_n,3.0e+0)
            v_q = MIN(v_q,3.0e+0)
            v_n_sedi(i,k,j) = v_n
            v_q_sedi(i,k,j) = v_q
          END DO
       END DO
    END DO
    v_n_sedi(:,kts-1,:) = v_n_sedi(:,kts,:)
    v_q_sedi(:,kts-1,:) = v_q_sedi(:,kts,:)
    DO j = jts,jte
! UB_20080313> Bugfix: k-loop has to be from model top to bottom for the flux-limiter!
       DO k = kte,kts,-1
!       DO k = kts,kte
! UB_20080313<
          DO i = its,ite

            v_n = 0.5 * (v_n_sedi(i,k-1,j)+v_n_sedi(i,k,j))
            v_q = 0.5 * (v_q_sedi(i,k-1,j)+v_q_sedi(i,k,j))
            ! Formulierung unter der Annahme, dass v_n, v_q stets positiv!!!
            c_n = v_n * adz(i,k,j) * dt 
            c_q = v_q * adz(i,k,j) * dt 
            IF (c_n > 1) THEN
              kk = k
              s_n = 0.0
              DO WHILE (c_n > 1 .AND. kk.LE.kte-1)
                s_n = s_n + qns(i,kk,j)/adz(i,kk,j)
                c_n = (c_n - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_n = s_n + qns(i,kk,j)/adz(i,kk,j)*MIN(c_n,1.0e0)
              s_n = s_n / dt
            ELSE
              s_n = v_n * qns(i,k,j)
            ENDIF
            IF (c_q > 1) THEN
              kk = k
              s_q = 0.0
              DO WHILE (c_q > 1 .AND. kk.LE.kte-1)
                s_q = s_q + qs(i,kk,j)/adz(i,kk,j)
                c_q = (c_q - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_q = s_q + qs(i,kk,j)/adz(i,kk,j)*MIN(c_q,1.0e0)
              s_q = s_q / dt
            ELSE
              s_q = v_q * qs(i,k,j)
            ENDIF

            n_fluss(i,k,j) = - s_n
            q_fluss(i,k,j) = - s_q

            ! Flux-limiter to avoid negative values
            n_fluss(i,k,j) = MAX(n_fluss(i,k,j),n_fluss(i,k+1,j)-qns(i,k,j)/(adz(i,k,j)*dt))
            q_fluss(i,k,j) = MAX(q_fluss(i,k,j),q_fluss(i,k+1,j)- qs(i,k,j)/(adz(i,k,j)*dt))

          END DO
       END DO
    END DO

    DO j = jts,jte
       DO i = its,ite
          n_fluss(i,kte+1,j) = 0.0 ! rho(i,kte,j) * v_n_sedi(i,kte,j) * qns(i,kte,j) 
          q_fluss(i,kte+1,j) = 0.0 ! rho(i,kte,j) * v_q_sedi(i,kte,j) * qs(i,kte,j)
          rainrate(i,j) = -q_fluss(i,kts,j)
       END DO
    END DO
    DO j = jts, jte
       DO k = kts,kte
          DO i = its,ite
             qns(i,k,j) = qns(i,k,j) + ( n_fluss(i,k,j) - n_fluss(i,k+1,j) ) * adz(i,k,j) * dt
             qns(i,k,j) = qns(i,k,j) / rho(i,k,j)
             qs(i,k,j)  = qs(i,k,j)  + ( q_fluss(i,k,j) - q_fluss(i,k+1,j) ) * adz(i,k,j) * dt
             qs(i,k,j)  = qs(i,k,j)  / rho(i,k,j)
          ENDDO
       ENDDO
    ENDDO             

  END SUBROUTINE snow_sedimentation_wrf

  SUBROUTINE ice_sedimentation_wrf (qi,qni,rainrate,rho,rhocorr,adz,dt, &
       &             ims,ime, jms,jme, kms,kme,               & ! memory dims
       &             its,ite, jts,jte, kts,kte                ) ! tile   dims
    !*******************************************************************************
    !       Berechnung der Sedimentation der Regentropfen                          *
    !*******************************************************************************

    USE konstanten

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme , &
         its,ite, jts,jte, kts,kte 
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(INOUT) :: qi,qni,adz
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),    &
         INTENT(IN)    :: rho,rhocorr
    REAL, DIMENSION( ims:ime , jms:jme ),    &
         INTENT(INOUT) :: rainrate
    REAL :: dt

    ! .. Local Variables ..
    INTEGER          :: i, j, k, kk
    INTEGER, SAVE    :: firstcall
    REAL :: x,n0,lam,v_n,v_q,n,s_n,s_q
    REAL, SAVE      :: alf_n,alf_q,c_lam   
    REAL, PARAMETER :: eps = 1.00e-20   

    REAL :: c_n, c_q
    REAL, DIMENSION(its:ite,kts:kte+1,jts:jte) :: q_fluss,n_fluss
    REAL, DIMENSION(its:ite,kts-1:kte,jts:jte) :: v_n_sedi,v_q_sedi

    IF (firstcall.NE.1) THEN
      alf_n = ice%a_vel * gfct((ice%nu+ice%b_vel+1.0)/ice%mu)/gfct((ice%nu+1.0)/ice%mu)
      alf_q = ice%a_vel * gfct((ice%nu+ice%b_vel+2.0)/ice%mu)/gfct((ice%nu+2.0)/ice%mu)
      c_lam = gfct((ice%nu+1.0)/ice%mu)/gfct((ice%nu+2.0)/ice%mu)
      IF(isdebug) THEN
        WRITE (6,*)           "CLOUD snow_sedi:"
        WRITE (6,'(A,D10.3)') "  Koeff. fuer lambda:          c_lam = ",c_lam
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation n: alf_n = ",alf_n
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation q: alf_q = ",alf_q
      ENDIF
      firstcall = 1
    ENDIF

    WHERE (qi < 0.0) qi = 0.0
    WHERE (qni < 0.0) qni = 0.0

    q_fluss = 0.0
    n_fluss = 0.0

    DO j = jts,jte
      DO k = kts,kte
        DO i = its,ite
          qi(i,k,j) = rho(i,k,j)*qi(i,k,j)
          qni(i,k,j) = rho(i,k,j)*qni(i,k,j)
          x = qi(i,k,j) / (qni(i,k,j)+eps)
          x = MIN(MAX(x,REAL(ice%x_min)),REAL(ice%x_max))
          lam = ( c_lam * x )**(ice%b_vel) * rhocorr(i,k,j)
          v_n = alf_n * lam                                 !..Fallgeschw. fuer n
          v_q = alf_q * lam                                 !..Fallgeschw. fuer q
          v_n = MAX(v_n,0.0e+0)
          v_q = MAX(v_q,0.0e+0)
          v_n = MIN(v_n,3.0e+0)
          v_q = MIN(v_q,3.0e+0)
          v_n_sedi(i,k,j) = v_n
          v_q_sedi(i,k,j) = v_q
        END DO
      END DO
    END DO
    v_n_sedi(:,kts-1,:) = v_n_sedi(:,kts,:)
    v_q_sedi(:,kts-1,:) = v_q_sedi(:,kts,:)
    DO j = jts,jte
! UB_20080313> Bugfix: k-loop has to be from model top to bottom for the flux-limiter!
       DO k = kte,kts,-1
!       DO k = kts,kte
! UB_20080313<
          DO i = its,ite

            v_n = 0.5 * (v_n_sedi(i,k-1,j)+v_n_sedi(i,k,j))
            v_q = 0.5 * (v_q_sedi(i,k-1,j)+v_q_sedi(i,k,j))
            ! Formulierung unter der Annahme, dass v_n, v_q stets positiv!!!
            c_n = v_n * adz(i,k,j) * dt 
            c_q = v_q * adz(i,k,j) * dt 
            IF (c_n > 1) THEN
              kk = k
              s_n = 0.0
              DO WHILE (c_n > 1 .AND. kk.LE.kte-1)
                s_n = s_n + qni(i,kk,j)/adz(i,kk,j)
                c_n = (c_n - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_n = s_n + qni(i,kk,j)/adz(i,kk,j)*MIN(c_n,1.0e0)
              s_n = s_n / dt
            ELSE
              s_n = v_n * qni(i,k,j)
            ENDIF
            IF (c_q > 1) THEN
              kk = k
              s_q = 0.0
              DO WHILE (c_q > 1 .AND. kk.LE.kte-1)
                s_q = s_q + qi(i,kk,j)/adz(i,kk,j)
                c_q = (c_q - 1) * adz(i,kk+1,j)/adz(i,kk,j) 
                kk  = kk + 1
              ENDDO
              s_q = s_q + qi(i,kk,j)/adz(i,kk,j)*MIN(c_q,1.0e0)
              s_q = s_q / dt
            ELSE
              s_q = v_q * qi(i,k,j)
            ENDIF

            n_fluss(i,k,j) = - s_n
            q_fluss(i,k,j) = - s_q

            ! Flux-limiter to avoid negative values
            n_fluss(i,k,j) = MAX(n_fluss(i,k,j),n_fluss(i,k+1,j)-qni(i,k,j)/(adz(i,k,j)*dt))
            q_fluss(i,k,j) = MAX(q_fluss(i,k,j),q_fluss(i,k+1,j)- qi(i,k,j)/(adz(i,k,j)*dt))

          END DO
       END DO
    END DO

    DO j = jts,jte
      DO i = its,ite
        n_fluss(i,kte+1,j) = 0.0 ! rho(i,kte,j) * v_n_sedi(i,kte,j) * qni(i,kte,j)
        q_fluss(i,kte+1,j) = 0.0 ! rho(i,kte,j) * v_q_sedi(i,kte,j) * qi(i,kte,j)
        rainrate(i,j) = -q_fluss(i,kts,j)
      END DO
    END DO
    DO j = jts, jte
      DO k = kts,kte
        DO i = its,ite
          qni(i,k,j) = qni(i,k,j) + ( n_fluss(i,k,j) - n_fluss(i,k+1,j) ) * adz(i,k,j) * dt
          qni(i,k,j) = qni(i,k,j) / rho(i,k,j)
          qi(i,k,j)  = qi(i,k,j)  + ( q_fluss(i,k,j) - q_fluss(i,k+1,j) ) * adz(i,k,j) * dt
          qi(i,k,j)  = qi(i,k,j)  / rho(i,k,j)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_sedimentation_wrf

#endif

END MODULE wolken_sedi

!==============================================================================

MODULE wolken_eis
  !*******************************************************************************
  !                                                                              *
  !     Diverse Unterprogramme zur Berechnung der Wolkenphysik mit Eisphase!     *
  !                                                                              *
  !*******************************************************************************

  USE wolken_konstanten 
! UB_20090227>>
  USE gamma_functions_mp_seifert
! UB_20090227<<
  USE konstanten,         ONLY: pi

  IMPLICIT NONE


CONTAINS

  SUBROUTINE ice_nucleation()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Nukleation der Wolkeneispartikel                        *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q,                &
         &                        q_cloud,q_ice,q_rain,q_snow,q_graupel,          &
         &                        n_ice, n_snow, p_g, t_g, rho_g,w,rho,dz3d,dt,   &
         &                        s_i,s_w,dSidz,dT0dz,                            &
         &                        dqdt,speichere_dqdt,                &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: r,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0,dichte
    USE geometrie,          ONLY: x3_x3

    IMPLICIT NONE

    ! Locale Variablen 
    DOUBLE PRECISION            :: q_si            !..Wasserdampfdichte bei Eissaettigung
    DOUBLE PRECISION            :: e_si            !..Wasserpartialdruck bei Eissaettigung
    DOUBLE PRECISION            :: p_a,rho_a,q_d,e_d,dTdz_w,dSdz_w
    DOUBLE PRECISION            :: nuc_n, nuc_q, ndiag
    DOUBLE PRECISION            :: a_ld,a_dl
    INTEGER                     :: i,j,k,nuc_typ,stat_var

    DOUBLE PRECISION, PARAMETER :: a_md = -0.639
    DOUBLE PRECISION, PARAMETER :: b_md = 12.960
    DOUBLE PRECISION, PARAMETER :: N_m0 = 1.0d+3
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    nuc_typ = nuc_i_typ

    IF(isIO() .AND. isdebug)THEN
      IF (nuc_typ == 1) THEN
        WRITE (6, *) "CLOUD ice_nucleation: fletcher-formel"
      ELSEIF (nuc_typ == 2) THEN
        WRITE (6, *) "CLOUD ice_nucleation: meyers-formel" 
      ELSEIF (nuc_typ == 3) THEN
        WRITE (6, *) "CLOUD ice_nucleation: cooper-formel" 
      ELSE
        WRITE (6, *) "CLOUD ice_nucleation: nuc_typ = ",nuc_typ
      ENDIF
      WRITE (6,*) "CLOUD ice_nucleation: max s_i    = ", MAXVAL(s_i(:,0,1))
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
            nuc_n = 0d0
            IF (nuc_typ < 8) THEN
              ndiag = n_ice_diagnostic(T_0(i,j,k),MIN(s_i(i,j,k),0.25d0),s_w(i,j,k),nuc_typ)
              nuc_n = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
            ELSEIF (nuc_typ == 8 .AND. dSidz(i,j,k) > 0.0 ) THEN
              dSdz_w = dSidz(i,j,k)
              nuc_n =  N_m0 * EXP( a_md + b_md * s_i(i,j,k) ) * b_md * dSdz_w * dt
              nuc_n = MAX(nuc_n,0.d0)
            ELSEIF (nuc_typ == 9 .AND. T_0(i,j,k) >= 245.15 .AND. &
                 (s_i(i,j,k) >= 0.08 .OR. (T_0(i,j,k) <= 265.16 .AND. s_w(i,j,k) >= 0.0))) THEN
              ! Zeitableitungsformulierung auf Basis der Cooper-Formel:
              dTdz_w = dT0dz(i,j,k)
              nuc_n = 5.0 * 0.304 * EXP( 0.304 * (T_3 - T_0(i,j,k)) ) * dTdz_w * dt
              nuc_n = MAX(nuc_n,0.d0)
            ENDIF

            nuc_q = MIN(nuc_n * ice%x_min, q(i,j,k))
            !nuc_n = nuc_q / ice%x_min                !AXEL 20040416

            n_ice(i,j,k) = n_ice(i,j,k) + nuc_n
            q_ice(i,j,k) = q_ice(i,j,k) + nuc_q
            q(i,j,k)     = q(i,j,k)     - nuc_q

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,2) = nuc_q
            END IF
#endif
            IF (speichere_precipstat) THEN
              cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + nuc_q
            END IF
            ! ub<<

            !WRF! T(i,j,k) = T(i,j,k) + L_ew / cp * nuc_q / rho_0(i,j,k)
            !WRF! p(i,j,k) = p(i,j,k) + L_ew / cv * nuc_q * r
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !DEALLOCATE(T_a,s_si)

  END SUBROUTINE ice_nucleation

  REAL(KIND=8) ELEMENTAL FUNCTION dep_growth_timescale(b,y,x_0)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(in) :: b,y,x_0
  
  REAL(KIND=8), PARAMETER  :: &
    SQ31 = 0.577350269,       &
    SIX1 = 0.166666666
  
  REAL(KIND=8) ::  X,F1,F10,F2,F20
    
  IF (y.LE.x_0) THEN
    dep_growth_timescale = 0.
  ELSE
    X     = MIN( y, 0.9999999 )
    F1    = SIX1 * LOG( (1.+X +X**2)  / (1.-X)**2 )
    F10   = SIX1 * LOG( (1.+x_0+x_0**2) / (1.-x_0)**2 )
    F2    = SQ31 * ATAN( SQ31*(1.+2.*X) )
    F20   = SQ31 * ATAN( SQ31*(1.+2.*x_0) )
    dep_growth_timescale = (b+1.)*(F1-F10) + (b-1.)*(F2-F20)
  END IF

  END FUNCTION dep_growth_timescale

  SUBROUTINE ice_nucleation_homhet()
  !*******************************************************************************
  !                                                                              *
  ! Berechnung der Nukleation der Wolkeneispartikel                              *
  !                                                                              *
  ! Nucleation scheme is based on the papers:                                    *
  !                                                                              *
  ! "A parametrization of cirrus cloud formation: Homogenous                     *
  ! freezing of supercooled aerosols" by B. Kaercher and                         *
  ! U. Lohmann 2002 (KL02 hereafter)                                             *
  !                                                                              *
  ! "Physically based parameterization of cirrus cloud formation                 *
  ! for use in global atmospheric models" by B. Kaercher, J. Hendricks           *
  ! and U. Lohmann 2006 (KHL06 hereafter)                                        *
  !                                                                              *
  !*******************************************************************************

  USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q,                &
    &                           q_cloud,q_ice,q_rain,q_snow,q_graupel,          &
    &                           n_ice, n_snow, p_g, t_g, rho_g,w,rho,dz3d,dt,   &
    &                           s_i,s_w,dSidz,dT0dz,                            &
    &                           dqdt,speichere_dqdt,                &
    &                           cond_neu_sb, evap_neu_sb, speichere_precipstat
  USE konstanten,         ONLY: pi,r,cv,cp,wolke_typ 
  USE parallele_umgebung, ONLY: isIO,abortparallel
  USE initialisierung,    ONLY: T_0,p_0,rho_0,dichte
  USE geometrie,          ONLY: x3_x3

  IMPLICIT NONE

  ! Locale Variablen 
  DOUBLE PRECISION            :: T_a,p_a,rho_a,q_d,e_d
  DOUBLE PRECISION            :: q_i,n_i,x_i,r_i
  DOUBLE PRECISION            :: ndiag, n_m, n_f
  INTEGER                     :: i,j,k,nuc_typ,stat_var

  DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

  ! switch for version of Phillips et al. scheme 
  ! (but make sure you have the correct INCLUDE file)
  INTEGER                     :: iphillips = 2010

  ! coeffs of Meyer formula
  DOUBLE PRECISION, PARAMETER :: a_md = -0.639
  DOUBLE PRECISION, PARAMETER :: b_md = 12.960
  DOUBLE PRECISION, PARAMETER :: N_m0 = 1.0d+3
  
  ! some constants needed for Kaercher and Lohmann approach
  DOUBLE PRECISION, PARAMETER :: &
    r_0     = 0.25e-6          , &    ! aerosol particle radius prior to freezing
    alpha_d = 0.5              , &    ! deposition coefficient (KL02; Spichtinger & Gierens 2009)
    k_b     = 1.38065e-23      , &    ! Boltzmann constant [J/K]
    M_w     = 18.01528e-3      , &    ! molecular mass of water [kg/mol]
    M_a     = 28.96e-3         , &    ! molecular mass of air [kg/mol]
    N_avo   = 6.022e23         , &    ! Avogadro number [1/mol]
    ma_w    = M_w / N_avo      , &    ! mass of water molecule [kg]
    grav    = 9.81             , &    ! acceleration of gravity [m/s2]
    p0      = 1013.25e2        , &
    svol    = ma_w / rho_ice          ! specific volume of a water molecule in ice
  
  REAL(KIND=8)  :: Si, e_si
  REAL(KIND=8)  :: ni_hom,ri_hom,mi_hom
  REAL(KIND=8)  :: v_th,n_sat,flux,phi,D_v,cool,tau,delta,w_p,scr
  REAL(KIND=8)  :: ctau, tau_g,acoeff(3),bcoeff(2), ri_dot
  REAL(KIND=8)  :: kappa,sqrtkap,ren,R_imfc,R_im,R_ik,ri_0
  REAL(KIND=8)  :: tgrow,ri_max,beta,xj,dxj,xmid,fmid,nimax
  REAL(KIND=8)  :: xt,xs
  INTEGER       :: jj,ss,tt

  LOGICAL :: use_hetnuc, use_homnuc, use_wp

  REAL(KIND=8), DIMENSION(3) :: infrac

  DOUBLE PRECISION, ALLOCATABLE :: nuc_n(:,:,:), nuc_q(:,:,:)

  nuc_typ = nuc_i_typ
    
  SELECT CASE (nuc_typ)
  CASE(0)
    ! Heterogeneous nucleation ONLY 
    use_hetnuc = .TRUE.
    use_homnuc = .FALSE.
    use_wp     = .FALSE.
  CASE(1)
    ! Homogeneous nucleation only (KHL06)"
    use_hetnuc = .FALSE.
    use_homnuc = .TRUE.
    use_wp     = .FALSE.
  CASE(2)
    ! Homog. and het. nucleation, but wp = 0 in KHL06"
    use_hetnuc = .TRUE.
    use_homnuc = .TRUE.
    use_wp     = .FALSE.
  CASE(3:9)
    ! Homog. and het. nucleation, using KHL06 scheme"
    use_hetnuc = .TRUE.
    use_homnuc = .TRUE.
    use_wp     = .TRUE.
  END SELECT


  IF(isIO() .AND. isdebug)THEN
    IF (use_hetnuc .AND. .NOT.use_homnuc) THEN
      WRITE (6, *) "CLOUD ice_nucleation: Heterogeneous nucleation only"
    ELSE IF (.NOT.use_hetnuc .AND. use_homnuc) THEN
      WRITE (6, *) "CLOUD ice_nucleation: Homogeneous nucleation only (KHL06)"
    ELSE IF ( use_hetnuc .AND. use_homnuc .AND. .NOT.use_wp ) THEN
      WRITE (6, *) "CLOUD ice_nucleation: Homog. and het. nucleation, but wp = 0 in KHL06"
    ELSE IF ( use_hetnuc .AND. use_homnuc .AND. use_wp ) THEN
      WRITE (6, *) "CLOUD ice_nucleation: Homog. and het. nucleation, using KHL06 scheme"
    END IF
  END IF

  ALLOCATE (nuc_n(0:loc_ix, 1:loc_iy, 1:loc_iz), stat=stat_var)
  ALLOCATE (nuc_q(0:loc_ix, 1:loc_iy, 1:loc_iz), stat=stat_var)
    
  nuc_n = 0.0d0
  nuc_q = 0.0d0

  ! Heterogenous nucleation using Meyers parameterization
  ! with an upper limit of n_het_max
  IF (use_hetnuc) THEN
    IF (nuc_typ.EQ.3) THEN
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            
            IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0 .AND. dSidz(i,j,k) > 0.0 &
              & .AND. ( n_ice(i,j,k)+n_snow(i,j,k) < ni_het_max ) )THEN
            
              nuc_n(i,j,k) =  N_m0 * EXP( a_md + b_md * s_i(i,j,k) ) * b_md * dSidz(i,j,k) * dt
              
              nuc_n(i,j,k) = MAX(nuc_n(i,j,k), 0.d0)
              nuc_q(i,j,k) = MIN(nuc_n(i,j,k) * ice%x_min, q(i,j,k))
              
              n_ice(i,j,k) = n_ice(i,j,k) + nuc_n(i,j,k)
              q_ice(i,j,k) = q_ice(i,j,k) + nuc_q(i,j,k)
              q(i,j,k)     = q(i,j,k)     - nuc_q(i,j,k)
            
            ENDIF

          ENDDO
        ENDDO
      ENDDO
    ELSE

      ! for sensitivity experiments
      IF (iphillips == 2008) THEN
        IF (nuc_typ.EQ.5) THEN  
          ! reduce soot and organics by factor of 10
          na_dust    = 162.e3 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e5 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e5 ! number density of organics [1/m3], Phillips08 value 177e6
        ELSEIF (nuc_typ.EQ.6) THEN 
          ! increase soot and organics by factor of 10
          na_dust    = 162.e3 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e7 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e7 ! number density of organics [1/m3], Phillips08 value 177e6
        ELSEIF (nuc_typ.EQ.7) THEN  
          ! increase dust by factor of 10
          na_dust    = 162.e4 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e6 ! number density of organics [1/m3], Phillips08 value 177e6
        ELSEIF (nuc_typ.EQ.8) THEN  
          ! increase dust by factor of 100
          na_dust    = 162.e5 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e6 ! number density of organics [1/m3], Phillips08 value 177e6
        ELSEIF (nuc_typ.EQ.9) THEN  
          ! increase dust by factor of 1000
          na_dust    = 162.e6 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e6 ! number density of organics [1/m3], Phillips08 value 177e6
        END IF
      ELSEIF (iphillips == 2010) THEN
        IF (nuc_typ.EQ.4) THEN  
          ! reduce organics by factor of 10, increase dust by 10
          na_dust    = 162.e4 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e5 ! number density of organics [1/m3], Phillips08 value 177e6
        END IF
        IF (nuc_typ.EQ.5) THEN  
          ! reduce soot by factor 10, reduce organics by factor of 10
          na_dust    = 162.e3 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e5 ! number density of organics [1/m3], Phillips08 value 177e6
        END IF
        IF (nuc_typ.EQ.6) THEN  
          ! increase dust by factor 100, reduce organics by factor of 10
          na_dust    = 162.e5 ! number density of dust [1/m³], Phillips08 value 162e3
          na_soot    =  15.e6 ! number density of soot [1/m³], Phillips08 value 15e6
          na_orga    = 177.e5 ! number density of organics [1/m3], Phillips08 value 177e6
        END IF
      END IF
      
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0  &
              & .AND. ( n_ice(i,j,k)+n_snow(i,j,k) < ni_het_max ) )THEN

              IF (q_cloud(i,j,k) > 0.0) THEN

                ! immersion freezing at water saturation
                xt = (274. - T_0(i,j,k))  / ttstep
                xt = MIN(xt,ttmax-1.)
                tt = INT(xt)
                infrac(1) = (tt+1-xt) * afrac_dust(tt,99) + (xt-tt) * afrac_dust(tt+1,99) 
                infrac(2) = (tt+1-xt) * afrac_soot(tt,99) + (xt-tt) * afrac_soot(tt+1,99) 
                infrac(3) = (tt+1-xt) * afrac_orga(tt,99) + (xt-tt) * afrac_orga(tt+1,99) 
                
              ELSE

                ! calculate indices used for look-up tables
                xt = (274. - T_0(i,j,k))  / ttstep
                xs = (100*s_i(i,j,k)) / ssstep    
                xt = MIN(xt,ttmax-1.)
                xs = MIN(xs,ssmax-1.)          
                tt = INT(xt)
                ss = INT(xs)
                
                ! bi-linear interpolation in look-up tables
                infrac(1) = (tt+1-xt) * (ss+1-xs) * afrac_dust(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_dust(tt+1,ss  ) &
                          + (tt+1-xt) * (xs-ss)   * afrac_dust(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_dust(tt+1,ss+1)
                infrac(2) = (tt+1-xt) * (ss+1-xs) * afrac_soot(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_soot(tt+1,ss  ) &
                          + (tt+1-xt) * (xs-ss)   * afrac_soot(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_soot(tt+1,ss+1)
                infrac(3) = (tt+1-xt) * (ss+1-xs) * afrac_orga(tt,ss  ) +  (xt-tt)*(ss+1-xs) * afrac_orga(tt+1,ss  ) &
                          + (tt+1-xt) * (xs-ss)   * afrac_orga(tt,ss+1) +  (xt-tt)*(xs-ss)   * afrac_orga(tt+1,ss+1)
              ENDIF

              ndiag = na_dust * infrac(1) + na_soot * infrac(2) + na_orga * infrac(3)
              ndiag = MIN(ndiag,ni_het_max)

              nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              
              nuc_q(i,j,k) = MIN(nuc_n(i,j,k) * ice%x_min, q(i,j,k))

              n_ice(i,j,k) = n_ice(i,j,k) + nuc_n(i,j,k)
              q_ice(i,j,k) = q_ice(i,j,k) + nuc_q(i,j,k)
              q(i,j,k)     = q(i,j,k)     - nuc_q(i,j,k)

            ENDIF

          END DO
        END DO
      END DO
    END IF

  END IF

  ! Homogeneous nucleation using KHL06 approach
  IF (use_homnuc) THEN
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          p_a  = p_0(i,j,k)
          T_a  = T_0(i,j,k)
          e_si = e_es(T_a)
          Si   = q(i,j,k) * R_d * T_a / e_si

          ! critical supersaturation for homogenous nucleation, cf. Eq. (1) of KB08
          Scr  = 2.349 - T_a / 259.00
          
          IF (Si > Scr  &
            .AND. n_ice(i,j,k) < ni_hom_max ) THEN

            n_i = n_ice(i,j,k) ! + n_snow(i,j,k) 
            q_i = q_ice(i,j,k) ! + q_snow(i,j,k) 
            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max) 
            !r_i = 0.5 * ice%a_geo * x_i**ice%b_geo  
            r_i = (x_i/(4./3.*pi*rho_ice))**(1./3.)

            D_v     = 0.211e-4 * (T_a/T_3)**1.94 * (p0/p_a)
            v_th    = SQRT( 8.*k_b*T_a/(pi*ma_w) ) 
            flux    = alpha_d * v_th/4.
            n_sat   = e_si / (k_b*T_a)  

            ! coeffs of supersaturation equation
            acoeff(1) = (L_ed * grav) / (cp * R_d * T_a**2) - grav/(R_l * T_a)
            acoeff(2) = 1.0/n_sat
            acoeff(3) = (L_ed**2 * M_w * ma_w)/(cp * p_a * T_a * M_a) 
          
            ! coeffs of depositional growth equation
            bcoeff(1) = flux * svol * n_sat * (Si - 1.)
            bcoeff(2) = flux / D_v

            IF (use_wp) THEN
              ! pre-existing ice crystals included as reduced updraft speed
              ri_dot = bcoeff(1) / (1. + bcoeff(2) * r_i)
              R_ik   = (4 * pi) / svol * n_i * r_i**2 * ri_dot
              w_p    = (acoeff(2) + acoeff(3) * Si)/(acoeff(1) * Si) * R_ik  ! KHL06 Eq. 19
              w_p    = MAX(dble(w_p),0.d0)
            ELSE
              w_p    = 0.d0
            END IF

            IF (w(i,j,k) > w_p) THEN   ! homogenous nucleation event

              ! timescales of freezing event (see KL02, RM05, KHL06)
              cool    = grav / cp * w(i,j,k)
              ctau    = T_a * ( 0.004*T_a - 2. ) + 304.4         
              tau     = 1.0 / (ctau * cool)                       ! freezing timescale, eq. (5)
              delta   = (bcoeff(2) * r_0)                         ! dimless aerosol radius, eq.(4)  
              tau_g   = (bcoeff(1) / r_0) / (1 + delta)           ! timescale for initial growth, eq.(4)
              phi     = acoeff(1)*Si / ( acoeff(2) + acoeff(3)*Si) * (w(i,j,k) - w_p) 
     
              ! monodisperse approximation following KHL06
              kappa   = 2. * bcoeff(1) * bcoeff(2) * tau / (1.+ delta)**2  ! kappa, Eq. 8 KHL06
              sqrtkap = SQRT(kappa)                                        ! root of kappa
              ren     = 3. * sqrtkap / ( 2. + SQRT(1.+9.*kappa/pi) )       ! analy. approx. of erfc by RM05
              R_imfc  = 4. * pi * bcoeff(1)/bcoeff(2)**2 / svol  
              R_im    = R_imfc / (1.+ delta) * ( delta**2 - 1. &
                & + (1.+0.5*kappa*(1.+ delta)**2) * ren/sqrtkap)           ! RIM Eq. 6 KHL06
              
              ! number concentration and radius of ice particles
              ni_hom  = phi / R_im                                         ! ni Eq.9 KHL06
              ri_0    = 1. + 0.5 * sqrtkap * ren                           ! for Eq. 3 KHL06
              ri_hom  = (ri_0 * (1. + delta) - 1. ) / bcoeff(2)            ! Eq. 3 KHL06 * REN = Eq.23 KHL06
              mi_hom  = (4./3. * pi * rho_ice) * ni_hom * ri_hom**3
              mi_hom  = MAX(dble(mi_hom),ice%x_min)

              nuc_n(i,j,k) = MAX(MIN(dble(ni_hom), ni_hom_max), 0.d0)
              nuc_q(i,j,k) = MIN(nuc_n(i,j,k) * mi_hom, q(i,j,k))
              
              n_ice(i,j,k) = n_ice(i,j,k) + nuc_n(i,j,k)
              q_ice(i,j,k) = q_ice(i,j,k) + nuc_q(i,j,k)
              q(i,j,k)     = q(i,j,k)     - nuc_q(i,j,k)
              
            END IF
          END IF

        ENDDO
      ENDDO
    ENDDO
  END IF

  DEALLOCATE(nuc_n, nuc_q)
    
  END SUBROUTINE ice_nucleation_homhet

  SUBROUTINE ice_nucleation_vec()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Nukleation der Wolkeneispartikel                        *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q,                &
         &                        q_cloud,q_ice,q_rain,q_snow,q_graupel,          &
         &                        n_ice, n_snow, p_g, t_g, rho_g,w,rho,dz3d,dt,   &
         &                        s_i,s_w,dSidz,dT0dz,                            &
         &                        dqdt,speichere_dqdt,                &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: r,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0,dichte
    USE geometrie,          ONLY: x3_x3

    IMPLICIT NONE

    ! Locale Variablen 
    DOUBLE PRECISION            :: p_a,rho_a,q_d,e_d,dTdz_w,dSdz_w,sitmp,eitmp,ewtmp
    DOUBLE PRECISION            :: ndiag, n_m, n_f
    DOUBLE PRECISION            :: a_ld,a_dl
    INTEGER                     :: i,j,k,nuc_typ,stat_var

    DOUBLE PRECISION, PARAMETER :: a_md = -0.639
    DOUBLE PRECISION, PARAMETER :: b_md = 12.960
    DOUBLE PRECISION, PARAMETER :: N_m0 = 1.0d+3
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION, ALLOCATABLE :: nuc_n(:,:,:), nuc_q(:,:,:)

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    nuc_typ = nuc_i_typ

    IF(isIO())THEN
    !IF(isIO() .AND. isdebug)THEN
      IF (nuc_typ == 1) THEN
        WRITE (6, *) "CLOUD ice_nucleation: fletcher-formel"
      ELSEIF (nuc_typ == 2) THEN
        WRITE (6, *) "CLOUD ice_nucleation: meyers-formel" 
      ELSEIF (nuc_typ == 3) THEN
        WRITE (6, *) "CLOUD ice_nucleation: cooper-formel" 
      ELSE
        WRITE (6, *) "CLOUD ice_nucleation: nuc_typ = ",nuc_typ
      ENDIF
      WRITE (6,*) "CLOUD ice_nucleation: max s_i    = ", MAXVAL(s_i)
    END IF

    ALLOCATE (nuc_n(0:loc_ix, 1:loc_iy, 1:loc_iz), stat=stat_var)
    ALLOCATE (nuc_q(0:loc_ix, 1:loc_iy, 1:loc_iz), stat=stat_var)
    
    nuc_n = 0.0d0
    nuc_q = 0.0d0

    IF (nuc_typ < 8) THEN

      SELECT CASE (nuc_typ)
      CASE(1)

        DO k = 1, loc_iz
          DO j = 1, loc_iy
            DO i = 0, loc_ix
              IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
                ndiag = n_ice_fletcher(T_0(i,j,k))
                nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              END IF
            END DO
          END DO
        END DO

      CASE(2)
        
        DO k = 1, loc_iz
          DO j = 1, loc_iy
            DO i = 0, loc_ix
              IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
                ndiag = n_ice_meyers( MIN(s_i(i,j,k),0.25d0) )
                nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              END IF
            END DO
          END DO
        END DO

      CASE(3)

        DO k = 1, loc_iz
          DO j = 1, loc_iy
            DO i = 0, loc_ix
              IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
                sitmp = MIN(s_i(i,j,k),0.25d0)
                ! Constraints suggested by Hugh Morrison:
                IF (sitmp >= 0.08 .OR. (T_0(i,j,k) <= 265.16 .AND. s_w(i,j,k) >= 0.0)) THEN
                  ndiag = n_ice_cooper( T_0(i,j,k) )
                ELSE
                  ndiag = 0.0d0
                END IF
                nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              END IF
            END DO
          END DO
        END DO


      CASE(4)

        DO k = 1, loc_iz
          DO j = 1, loc_iy
            DO i = 0, loc_ix
              IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
                sitmp = MIN(s_i(i,j,k),0.25d0)
                ewtmp = e_ws( T_0(i,j,k) )
                eitmp = e_es( T_0(i,j,k) )
                ndiag = n_ice_huffmann_vali_vec( T_0(i,j,k), sitmp, ewtmp, eitmp )
                nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              END IF
            END DO
          END DO
        END DO
        
      CASE(6)

        DO k = 1, loc_iz
          DO j = 1, loc_iy
            DO i = 0, loc_ix
              IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
                ndiag = n_ice_fletcher_contact( T_0(i,j,k) )
                nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              END IF
            END DO
          END DO
        END DO

      CASE(7)

        DO k = 1, loc_iz
          DO j = 1, loc_iy
            DO i = 0, loc_ix
              IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN
                N_m = n_ice_meyers( MIN(s_i(i,j,k),0.25d0) )
                N_f = n_ice_fletcher( T_0(i,j,k) )
                ndiag = MAX(MIN(MAX(N_m,0.1*N_f),10.0*N_f),0.01d3)
                nuc_n(i,j,k) = MAX( ndiag - (n_ice(i,j,k)+n_snow(i,j,k)),0.d0)
              END IF
            END DO
          END DO
        END DO
        
      END SELECT
    
    ELSEIF (nuc_typ == 8) THEN

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0 .AND. dSidz(i,j,k) > 0.0) THEN
              sitmp = MIN(s_i(i,j,k),0.25d0)
              dSdz_w = dSidz(i,j,k)
              nuc_n(i,j,k) =  N_m0 * EXP( a_md + b_md * sitmp ) * b_md * dSdz_w * dt
              nuc_n(i,j,k) = MAX(nuc_n(i,j,k),0.d0)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ELSEIF (nuc_typ == 9 .AND. T_0(i,j,k) >= 245.15 .AND. &
         (s_i(i,j,k) >= 0.08 .OR. (T_0(i,j,k) <= 265.16 .AND. s_w(i,j,k) >= 0.0))) THEN

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0 .AND. &
                 T_0(i,j,k) >= 245.15 .AND. &
                 (s_i(i,j,k) >= 0.08 .OR. (T_0(i,j,k) <= 265.16 .AND. s_w(i,j,k) >= 0.0)) &
                 ) THEN
              ! Zeitableitungsformulierung auf Basis der Cooper-Formel:
              dTdz_w = dT0dz(i,j,k)
              nuc_n(i,j,k) = 5.0 * 0.304 * EXP( 0.304 * (T_3 - T_0(i,j,k)) ) * dTdz_w * dt
              nuc_n(i,j,k) = MAX(nuc_n(i,j,k),0.d0)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
    ENDIF

    ! Finalize:
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          IF (T_0(i,j,k) < T_nuc .AND. s_i(i,j,k) > 0.0) THEN

            nuc_q(i,j,k) = MIN(nuc_n(i,j,k) * ice%x_min, q(i,j,k))

            n_ice(i,j,k) = n_ice(i,j,k) + nuc_n(i,j,k)
            q_ice(i,j,k) = q_ice(i,j,k) + nuc_q(i,j,k)
            q(i,j,k)     = q(i,j,k)     - nuc_q(i,j,k)

          ENDIF
        ENDDO
      ENDDO
    ENDDO

#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,2) = nuc_q
    END IF
#endif
    IF (speichere_precipstat) THEN
      cond_neu_sb = cond_neu_sb + nuc_q
    END IF

    DEALLOCATE(nuc_n, nuc_q)

  END SUBROUTINE ice_nucleation_vec

  DOUBLE PRECISION FUNCTION n_ice_diagnostic(T_a,S,S_w,nuc_typ)
    IMPLICIT NONE

    INTEGER           :: nuc_typ
    DOUBLE PRECISION  :: T_a  !..Absolute Temperatur
    DOUBLE PRECISION  :: S    !..Uebersaettigung bzgl. Eis
    DOUBLE PRECISION  :: S_w  !..Uebersaettigung bzgl. Liquid
    DOUBLE PRECISION  :: N_m,N_f

    S = MIN(10.d0,S)

    IF (nuc_typ == 1) THEN
      n_ice_diagnostic = n_ice_fletcher(T_a)
    ELSEIF (nuc_typ == 2) THEN
      n_ice_diagnostic = n_ice_meyers(S)
    ELSEIF (nuc_typ == 3) THEN
      ! Constraints suggested by Hugh Morrison:
      IF (S >= 0.08 .OR. (T_a <= 265.16 .AND. S_w >= 0.0)) THEN
        n_ice_diagnostic = n_ice_cooper(T_a)
      ELSE
        n_ice_diagnostic = 0.0
      END IF
    ELSEIF (nuc_typ == 4) THEN
      n_ice_diagnostic = n_ice_huffmann_vali(T_a,S)
    ELSEIF (nuc_typ == 6) THEN
      n_ice_diagnostic = n_ice_fletcher_contact(T_a)
    ELSEIF (nuc_typ == 7) THEN
      N_m = n_ice_meyers(S)
      N_f = n_ice_fletcher(T_a)
      n_ice_diagnostic = MAX(MIN(MAX(N_m,0.1*N_f),10.0*N_f),0.01d3)
    ENDIF

    RETURN
  END FUNCTION n_ice_diagnostic

  DOUBLE PRECISION FUNCTION n_ice_meyers_contact(T_a,S)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Meyers (1992) 
    IMPLICIT NONE

    DOUBLE PRECISION            :: S,T_a
    DOUBLE PRECISION, PARAMETER :: N_0 = 1.0d+3 
    DOUBLE PRECISION, PARAMETER :: N_m = 1.0d+3 
    DOUBLE PRECISION, PARAMETER :: a_d = -0.639
    DOUBLE PRECISION, PARAMETER :: b_d = 12.960
    DOUBLE PRECISION, PARAMETER :: c_d = -2.8
    DOUBLE PRECISION, PARAMETER :: d_d = 0.262

    n_ice_meyers_contact = N_0 * EXP( a_d + b_d * S )         &
         &               + N_m * EXP( c_d + d_d * (T_a - T_3) )

    RETURN
  END FUNCTION n_ice_meyers_contact

  DOUBLE PRECISION FUNCTION n_ice_fletcher(T_a)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Fletcher (1962) 
    ! mit Temperaturlimit nach Reisner et al. (1998) 
    IMPLICIT NONE

    DOUBLE PRECISION            :: T_a            !..Absolute Temperatur
    DOUBLE PRECISION, PARAMETER :: T_m = 246.00   !..Temperaturlimit
    DOUBLE PRECISION, PARAMETER :: N_0 = 1.0d-2
    DOUBLE PRECISION, PARAMETER :: k_f = 0.6

    n_ice_fletcher = N_0 * EXP( - k_f * (MAX(T_a,T_m) - T_3) )

    RETURN
  END FUNCTION n_ice_fletcher

  DOUBLE PRECISION FUNCTION n_ice_fletcher_contact(T_a)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Fletcher (1962) 
    ! mit Temperaturlimit nach Reisner et al. (1998)
    ! + Contactnukleation Meyers (1992)
    IMPLICIT NONE

    DOUBLE PRECISION            :: T_a            !..Absolute Temperatur
    DOUBLE PRECISION, PARAMETER :: T_m = 246.00   !..Temperaturlimit
    DOUBLE PRECISION, PARAMETER :: N_f = 1.0d-2
    DOUBLE PRECISION, PARAMETER :: k_f = 0.6
    DOUBLE PRECISION, PARAMETER :: c = -2.8
    DOUBLE PRECISION, PARAMETER :: d = 0.262
    DOUBLE PRECISION, PARAMETER :: N_m = 1.0d+3 

    n_ice_fletcher_contact = N_f * EXP( - k_f * (MAX(T_a,T_m) - T_3) )&
         &                 + N_m * EXP( c + d * (T_a - T_3) )

    RETURN
  END FUNCTION n_ice_fletcher_contact

  DOUBLE PRECISION FUNCTION n_ice_fletcher_DT(T_a)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Fletcher (1962) 
    ! -> Ableitung nach T_a (siehe Murakami, 1990)
    IMPLICIT NONE

    DOUBLE PRECISION            :: T_a            !..Absolute Temperatur
    DOUBLE PRECISION, PARAMETER :: N_0 = 1.0d-2
    DOUBLE PRECISION, PARAMETER :: k_f = 0.6

    n_ice_fletcher_DT = k_f * N_0 * EXP( - k_f * (T_a - T_3) )

    RETURN
  END FUNCTION n_ice_fletcher_DT

  DOUBLE PRECISION FUNCTION n_ice_huffmann_vali(T_a,S)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Fletcher (1962) 
    ! und Huffmann und Vali (1973) (siehe auch Murakami, 1990)
    ! mit Temperaturlimit nach Reisner et al. (1998) 
    IMPLICIT NONE

    DOUBLE PRECISION            :: T_a            !..Absolute Temperatur
    DOUBLE PRECISION            :: S              !..Absolute Temperatur
    DOUBLE PRECISION, PARAMETER :: T_m = 246.00   !..Temperaturlimit
    DOUBLE PRECISION, PARAMETER :: N_0 = 1.0d-2
    DOUBLE PRECISION, PARAMETER :: k_f = 0.6
    DOUBLE PRECISION, PARAMETER :: B   = 4.5

    DOUBLE PRECISION            :: e_w,e_i,S_0

    e_w = e_ws(T_a)
    e_i = e_es(T_a)
    S_0 = e_i/e_w*(S+1.0)-1.0

    n_ice_huffmann_vali = N_0 * (S/S_0)**B * EXP( - k_f * (MAX(T_a,T_m) - T_3) )

    RETURN

  END FUNCTION n_ice_huffmann_vali

  DOUBLE PRECISION FUNCTION n_ice_huffmann_vali_vec(T_a,S,e_w,e_i)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Fletcher (1962) 
    ! und Huffmann und Vali (1973) (siehe auch Murakami, 1990)
    ! mit Temperaturlimit nach Reisner et al. (1998) 
    IMPLICIT NONE

    DOUBLE PRECISION            :: T_a            !..Absolute Temperatur
    DOUBLE PRECISION            :: S              !..Absolute Temperatur
    DOUBLE PRECISION, PARAMETER :: T_m = 246.00   !..Temperaturlimit
    DOUBLE PRECISION, PARAMETER :: N_0 = 1.0d-2
    DOUBLE PRECISION, PARAMETER :: k_f = 0.6
    DOUBLE PRECISION, PARAMETER :: B   = 4.5

    DOUBLE PRECISION            :: e_w,e_i,S_0

!    e_w = e_ws(T_a)
!    e_i = e_es(T_a)
    S_0 = e_i/e_w*(S+1.0)-1.0

    n_ice_huffmann_vali_vec = N_0 * (S/S_0)**B * EXP( - k_f * (MAX(T_a,T_m) - T_3) )

    RETURN

  END FUNCTION n_ice_huffmann_vali_vec

  DOUBLE PRECISION FUNCTION n_ice_cooper(T_a)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Cooper (1986) 
    ! (siehe auch Rasmussen et al 2002, JAS; Thompson etal 2004, MWR)
    IMPLICIT NONE

    DOUBLE PRECISION            :: T_a            !..Absolute Temperatur
! ub>> Suggestion by Hugh Morrison and Roy Rasmussen and Greg Thompson:
!    DOUBLE PRECISION, PARAMETER :: T_m = 233.00   !..Temperaturlimit
    DOUBLE PRECISION, PARAMETER :: T_m = 238.15   !..Temperaturlimit -35 Grad C
! ub<<
    DOUBLE PRECISION, PARAMETER :: N_0 = 5.0      
    DOUBLE PRECISION, PARAMETER :: k_f = 0.304

! weiteres limit: hoechstens 500 IN pro Liter
    n_ice_cooper = MIN( 5.0d5, N_0 * EXP( k_f * (T_3 - MAX(T_a,T_m)) ) )

    RETURN
  END FUNCTION n_ice_cooper

  DOUBLE PRECISION FUNCTION n_ice_meyers(S)
    ! Diagnostische Beziehung fuer Anzahldichte der Eisteilchen nach Meyers (1992) 
    IMPLICIT NONE

    DOUBLE PRECISION            :: S              !..Uebersaettigung bzgl. Eis
    DOUBLE PRECISION, PARAMETER :: a_d = -0.639
    DOUBLE PRECISION, PARAMETER :: b_d = 12.960
    DOUBLE PRECISION, PARAMETER :: N_0 = 1.0d+3 

    n_ice_meyers = N_0 * EXP( a_d + b_d * S )

    RETURN
  END FUNCTION n_ice_meyers

  SUBROUTINE ice_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Eispartikel                               *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_ice, n_ice, q_cloud, n_cloud, q_graupel, n_graupel,dt, dqdt, speichere_dqdt, &
         &                        rimeqcrate_ice, rimencrate_ice

    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i
    DOUBLE PRECISION            :: q_c,n_c,x_c,d_c,v_c,e_coll,x_coll_c
    DOUBLE PRECISION            :: rime_n,rime_q
    DOUBLE PRECISION            :: conv_n,conv_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_ii,delta_n_ic,delta_n_cc
    DOUBLE PRECISION, SAVE      :: delta_q_ii,delta_q_ic,delta_q_cc
    DOUBLE PRECISION, SAVE      :: theta_n_ii,theta_n_ic,theta_n_cc
    DOUBLE PRECISION, SAVE      :: theta_q_ii,theta_q_ic,theta_q_cc
    DOUBLE PRECISION            :: const1,const2,const3,const4,const5

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      IF (ice_multiplication) THEN
        WRITE (6, *) "CLOUD ice_cloud_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD ice_cloud_riming ohne ice_multiplication"
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_ii = coll_delta_11(ice,cloud,0)
      delta_n_ic = coll_delta_12(ice,cloud,0)
      delta_n_cc = coll_delta_22(ice,cloud,0)
      delta_q_ii = coll_delta_11(ice,cloud,0) 
      delta_q_ic = coll_delta_12(ice,cloud,1)
      delta_q_cc = coll_delta_22(ice,cloud,1)

      theta_n_ii = coll_theta_11(ice,cloud,0)
      theta_n_ic = coll_theta_12(ice,cloud,0)
      theta_n_cc = coll_theta_22(ice,cloud,0)
      theta_q_ii = coll_theta_11(ice,cloud,0)
      theta_q_ic = coll_theta_12(ice,cloud,1)
      theta_q_cc = coll_theta_22(ice,cloud,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:        a_ice      = ",ice%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_ice      = ",ice%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Eis:            alf_ice    = ",ice%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_ice    = ",ice%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentr.:  a_cloud    = ",cloud%a_geo
        WRITE (6,'(A,D10.3)') "                                        b_cloud    = ",cloud%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Wolkentr.:      alf_cloud  = ",cloud%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_cloud  = ",cloud%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming ice-cloud:         delta_n_ii = ",delta_n_ii
        WRITE (6,'(A,D10.3)') "                                        delta_n_ic = ",delta_n_ic
        WRITE (6,'(A,D10.3)') "                                        delta_n_cc = ",delta_n_cc
        WRITE (6,'(A,D10.3)') "                                        theta_n_ii = ",theta_n_ii
        WRITE (6,'(A,D10.3)') "                                        theta_n_ic = ",theta_n_ic
        WRITE (6,'(A,D10.3)') "                                        theta_n_cc = ",theta_n_cc
        WRITE (6,'(A,D10.3)') "                                        delta_q_ii = ",delta_q_ii
        WRITE (6,'(A,D10.3)') "                                        delta_q_ic = ",delta_q_ic
        WRITE (6,'(A,D10.3)') "                                        delta_q_cc = ",delta_q_cc
        WRITE (6,'(A,D10.3)') "                                        theta_q_ii = ",theta_q_ii
        WRITE (6,'(A,D10.3)') "                                        theta_q_ic = ",theta_q_ic
        WRITE (6,'(A,D10.3)') "                                        theta_q_cc = ",theta_q_cc
      END IF
      firstcall = 1
    ENDIF

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..Mindestmasse fuer collection, begrenzt rime_n

    const1 = e_ic/(D_coll_c - D_krit_c)
    const2 = 1/x_coll_c
    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_c = n_cloud(i,j,k)                                   !..Anzahldichte
          q_c = q_cloud(i,j,k)                                   !..Fluessigwassergehalt
          n_i = n_ice(i,j,k)                                     !..Anzahldichte
          q_i = q_ice(i,j,k)                                     !..Massendichte

          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)  !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                   !..mittlerer Durchmesser

          x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)      !..mittlere Masse
          D_i = ice%a_geo * x_i**ice%b_geo                       !..mittlerer Durchmesser

          T_a = T_0(i,j,k) !WRF! + T_0(i,j,k) + T_g(i,j,k)

          IF (q_c > q_krit_c .AND. q_i > q_krit_ic .AND. D_i > D_krit_ic .AND. D_c > D_krit_c) THEN

            v_c = cloud%a_vel * x_c**cloud%b_vel * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.
            v_i = ice%a_vel   * x_i**ice%b_vel   * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            e_coll = MIN(e_ic, MAX(const1*(D_c - D_krit_c), e_min))

            rime_n = pi/4.0 * e_coll * n_i * n_c * dt & 
                 &   *     (delta_n_ii * D_i*D_i + delta_n_ic * D_i*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_ii * v_i*v_i - theta_n_ic * v_i*v_c + theta_n_cc * v_c*v_c  &
                 &         +ice_s_vel**2)

            rime_q = pi/4.0 * e_coll * n_i * q_c * dt & 
                 &   *     (delta_q_ii * D_i*D_i + delta_q_ic * D_i*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_ii * v_i*v_i - theta_q_ic * v_i*v_c + theta_q_cc * v_c*v_c  &
                 &          +ice_s_vel**2)

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_q = MIN(q_c,rime_q)
!ub              rime_n = MIN(n_c,MIN(rime_n,const2*rime_q))
              rime_n = MIN(n_c,rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n
              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,11) = rime_q
              END IF
#endif
              ! ub<<

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
              ENDIF

              ! Umwandlung ice -> graupel

              IF (D_i > D_conv_ig) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_i**3/x_i - 1.0) )
                ! D_i darf durch CONV nicht kleiner werden als D_conv_ig
                !conv_q = MIN(q_ice(i,j,k)-n_i*(D_conv_ig/ice%a_geo)**(1.0/ice%b_geo),conv_q)
                conv_q = MIN(q_ice(i,j,k),conv_q)
                ! ub >>
                x_i = MIN(MAX((q_ice(i,j,k))/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse incl. Riming
                ! ub <<
                conv_n = conv_q / MAX(x_i,x_conv) 
                conv_n = MIN(n_ice(i,j,k),conv_n)
              ELSE
                conv_q = 0.0
                conv_n = 0.0
              ENDIF

              q_ice(i,j,k)     = q_ice(i,j,k)     - conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q

              n_ice(i,j,k)     = n_ice(i,j,k)     - conv_n
              n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n
              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,12) = conv_q
              END IF
#endif

            ELSE

              rimeqcrate_ice(i,j,k) = rimeqcrate_ice(i,j,k) + rime_q
              rimencrate_ice(i,j,k) = rimencrate_ice(i,j,k) + rime_n

            END IF
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_cloud_riming

  SUBROUTINE snow_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Schneepartikel                            *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, q_ice, n_ice, &
         &                        q_snow, n_snow, q_cloud, n_cloud, q_graupel, n_graupel,dt, dqdt, speichere_dqdt, &
         &                        rimeqcrate_snow, rimencrate_snow

    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s
    DOUBLE PRECISION            :: q_c,n_c,x_c,d_c,v_c,e_coll
    DOUBLE PRECISION            :: rime_n,rime_q
    DOUBLE PRECISION            :: conv_n,conv_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_ss,delta_n_sc,delta_n_cc
    DOUBLE PRECISION, SAVE      :: delta_q_ss,delta_q_sc,delta_q_cc
    DOUBLE PRECISION, SAVE      :: theta_n_ss,theta_n_sc,theta_n_cc
    DOUBLE PRECISION, SAVE      :: theta_q_ss,theta_q_sc,theta_q_cc
    DOUBLE PRECISION            :: const1,const3,const4,const5

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      IF (ice_multiplication) THEN
        WRITE (6, *) "CLOUD snow_cloud_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD snow_cloud_riming "
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,cloud,0)
      delta_n_sc = coll_delta_12(snow,cloud,0)
      delta_n_cc = coll_delta_22(snow,cloud,0)
      delta_q_ss = coll_delta_11(snow,cloud,0) 
      delta_q_sc = coll_delta_12(snow,cloud,1)
      delta_q_cc = coll_delta_22(snow,cloud,1)

      theta_n_ss = coll_theta_11(snow,cloud,0)
      theta_n_sc = coll_theta_12(snow,cloud,0)
      theta_n_cc = coll_theta_22(snow,cloud,0)
      theta_q_ss = coll_theta_11(snow,cloud,0)
      theta_q_sc = coll_theta_12(snow,cloud,1)
      theta_q_cc = coll_theta_22(snow,cloud,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:        a_snow      = ",snow%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_snow      = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Eis:            alf_snow    = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_snow    = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentr.:  a_cloud    = ",cloud%a_geo
        WRITE (6,'(A,D10.3)') "                                        b_cloud    = ",cloud%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Wolkentr.:      alf_cloud  = ",cloud%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_cloud  = ",cloud%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming snow-cloud:        delta_n_ss = ",delta_n_ss
        WRITE (6,'(A,D10.3)') "                                        delta_n_sc = ",delta_n_sc
        WRITE (6,'(A,D10.3)') "                                        delta_n_cc = ",delta_n_cc
        WRITE (6,'(A,D10.3)') "                                        theta_n_ss = ",theta_n_ss
        WRITE (6,'(A,D10.3)') "                                        theta_n_sc = ",theta_n_sc
        WRITE (6,'(A,D10.3)') "                                        theta_n_cc = ",theta_n_cc
        WRITE (6,'(A,D10.3)') "                                        delta_q_ss = ",delta_q_ss
        WRITE (6,'(A,D10.3)') "                                        delta_q_sc = ",delta_q_sc
        WRITE (6,'(A,D10.3)') "                                        delta_q_cc = ",delta_q_cc
        WRITE (6,'(A,D10.3)') "                                        theta_q_ss = ",theta_q_ss
        WRITE (6,'(A,D10.3)') "                                        theta_q_sc = ",theta_q_sc
        WRITE (6,'(A,D10.3)') "                                        theta_q_cc = ",theta_q_cc
      END IF
      firstcall = 1
    ENDIF

! UB_20090316: re-invent constx 
    const1 = e_ic/(D_coll_c - D_krit_c)
    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_c = n_cloud(i,j,k)                                   !..Anzahldichte
          q_c = q_cloud(i,j,k)                                   !..Fluessigwassergehalt
          n_s = n_snow(i,j,k)                                    !..Anzahldichte
          q_s = q_snow(i,j,k)                                    !..Massendichte

          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)  !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                   !..mittlerer Durchmesser

          x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse
          D_s = snow%a_geo * x_s**snow%b_geo                     !..mittlerer Durchmesser

          T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

          IF (q_c > q_krit_c .AND. q_s > q_krit_sc .AND. D_s > D_krit_sc .AND. D_c > D_krit_c) THEN

            v_c = cloud%a_vel * x_c**cloud%b_vel * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.
            v_s = snow%a_vel  * x_s**snow%b_vel  * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            e_coll = MIN(e_ic, MAX(const1*(D_c - D_krit_c), e_min))

            rime_n = pi/4.0 * e_coll * n_s * n_c * dt & 
                 &   * (delta_n_ss * D_s**2 + delta_n_sc * D_s*D_c + delta_n_cc * D_c**2) &
                 &   * (theta_n_ss * v_s**2 - theta_n_sc * v_s*v_c + theta_n_cc * v_c**2  &
                 &     +snow_s_vel**2)**0.5

            rime_q = pi/4.0 * e_coll * n_s * q_c * dt & 
                 &   * (delta_q_ss * D_s**2 + delta_q_sc * D_s*D_c + delta_q_cc * D_c**2) &
                 &   * (theta_q_ss * v_s**2 - theta_q_sc * v_s*v_c + theta_q_cc * v_c**2  &
                 &     +snow_s_vel**2)**0.5

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_q = MIN(q_c,rime_q)
              rime_n = MIN(n_c,rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n
              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,13) = rime_q
              END IF
#endif
              ! ub<<

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,14) = mult_q
              END IF
#endif
              ! ub<<

              ! Umwandlung snow -> graupel

              IF (D_s > D_conv_sg) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_s**3/x_s - 1.0) )
                !conv_q = MIN(q_snow(i,j,k)-n_s*(D_conv_sg/snow%a_geo)**(1.0/snow%b_geo),conv_q)
                !conv_q = MAX(0.d0,conv_q)
                conv_q = MIN(q_snow(i,j,k),conv_q)
                ! ub >>
                x_s = MIN(MAX((q_snow(i,j,k))/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse incl. Riming
                ! ub <<
                conv_n = conv_q / MAX(x_s,x_conv) 
                conv_n = MIN(n_snow(i,j,k),conv_n)
              ELSE
                conv_q = 0.0
                conv_n = 0.0
              ENDIF

              q_snow(i,j,k)    = q_snow(i,j,k)    - conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q

              n_snow(i,j,k)    = n_snow(i,j,k)    - conv_n
              n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n
              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,15) = conv_q
              END IF
#endif

            ELSE

              rimeqcrate_snow(i,j,k) = rimeqcrate_snow(i,j,k) + rime_q
              rimencrate_snow(i,j,k) = rimencrate_snow(i,j,k) + rime_n

            END IF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_cloud_riming

  SUBROUTINE graupel_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Graupelpartikel                           *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g,        &
         &                        q_graupel, n_graupel, q_cloud, n_cloud, q_rain, n_rain,  &
         &                        q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g
    DOUBLE PRECISION            :: q_c,n_c,x_c,d_c,v_c,x_coll_c
    DOUBLE PRECISION            :: rime_n,rime_q,e_coll_n
    DOUBLE PRECISION            :: melt_n,melt_q,e_coll_q
    DOUBLE PRECISION            :: shed_n,shed_q,x_shed
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_gg,delta_n_gc,delta_n_cc
    DOUBLE PRECISION, SAVE      :: delta_q_gg,delta_q_gc,delta_q_cc
    DOUBLE PRECISION, SAVE      :: theta_n_gg,theta_n_gc,theta_n_cc
    DOUBLE PRECISION, SAVE      :: theta_q_gg,theta_q_gc,theta_q_cc
    DOUBLE PRECISION            :: const1,const2,const3,const4

    DOUBLE PRECISION, PARAMETER :: eps = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: riming_n,riming_q,rate_q

    ALLOCATE(riming_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(riming_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD graupel_cloud_riming",4)

    IF (isIO() .AND. isdebug) THEN
      IF (enhanced_melting .AND. ice_multiplication) THEN
        WRITE (6, *) "CLOUD graupel_cloud_riming mit ice_multiplication und enhanced_melting"
      ELSEIF (enhanced_melting .AND. .NOT. ice_multiplication) THEN
        WRITE (6, *) "CLOUD graupel_cloud_riming mit enhanced_melting"
      ELSEIF (ice_multiplication .AND. .NOT. enhanced_melting) THEN
        WRITE (6, *) "CLOUD graupel_cloud_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD graupel_cloud_riming" 
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,cloud,0)
      delta_n_gc = coll_delta_12(graupel,cloud,0)
      delta_n_cc = coll_delta_22(graupel,cloud,0)
      delta_q_gg = coll_delta_11(graupel,cloud,0) 
      delta_q_gc = coll_delta_12(graupel,cloud,1)
      delta_q_cc = coll_delta_22(graupel,cloud,1)

      theta_n_gg = coll_theta_11(graupel,cloud,0)
      theta_n_gc = coll_theta_12(graupel,cloud,0)
      theta_n_cc = coll_theta_22(graupel,cloud,0)
      theta_q_gg = coll_theta_11(graupel,cloud,0)
      theta_q_gc = coll_theta_12(graupel,cloud,1)
      theta_q_cc = coll_theta_22(graupel,cloud,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:   a_graupel   = ",graupel%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_graupel   = ",graupel%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Graupel:       alf_graupel = ",graupel%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_graupel = ",graupel%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentr.: a_cloud     = ",cloud%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_cloud     = ",cloud%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Wolkentr.:     alf_cloud   = ",cloud%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_cloud   = ",cloud%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming graupel-cloud:    delta_n_gg  = ",delta_n_gg
        WRITE (6,'(A,D10.3)') "                                       delta_n_gc  = ",delta_n_gc
        WRITE (6,'(A,D10.3)') "                                       delta_n_cc  = ",delta_n_cc
        WRITE (6,'(A,D10.3)') "                                       theta_n_gg  = ",theta_n_gg
        WRITE (6,'(A,D10.3)') "                                       theta_n_gc  = ",theta_n_gc
        WRITE (6,'(A,D10.3)') "                                       theta_n_cc  = ",theta_n_cc
        WRITE (6,'(A,D10.3)') "                                       delta_q_gg  = ",delta_q_gg
        WRITE (6,'(A,D10.3)') "                                       delta_q_gc  = ",delta_q_gc
        WRITE (6,'(A,D10.3)') "                                       delta_q_cc  = ",delta_q_cc
        WRITE (6,'(A,D10.3)') "                                       theta_q_gg  = ",theta_q_gg
        WRITE (6,'(A,D10.3)') "                                       theta_q_gc  = ",theta_q_gc
        WRITE (6,'(A,D10.3)') "                                       theta_q_cc  = ",theta_q_cc
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3.*pi * rho_w * r_shedding**3     !..Mittlere Masse der Sheddingtropfen

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..Mindestmasse fuer collection, begrenzt rime_n

    const1 = e_gc/(D_coll_c - D_krit_c)    
    const2 = 1.0/(T_mult_opt - T_mult_min)
    const3 = 1.0/(T_mult_opt - T_mult_max)
    const4 = c_w / L_ew
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix             
          q_c = q_cloud(i,j,k)                                      !..Massendichte
          q_g = q_graupel(i,j,k)                                    !..Massendichte
          n_c = n_cloud(i,j,k)                                      !..Anzahldichte
          n_g = n_graupel(i,j,k)                                    !..Anzahldichte
          x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse
          D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)     !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                      !..mittlerer Durchmesser
          T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)                  !..abs. Temperatur

          IF (q_c > q_krit_c .AND. q_g > q_krit_gc .AND. D_g > D_krit_gc .AND. &
               & D_c > D_krit_c) THEN

            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            v_c = cloud%a_vel   * x_c**cloud%b_vel   * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.

            e_coll_n = MIN(e_gc, MAX(const1*(D_c - D_krit_c),e_min))
            e_coll_q = e_coll_n

            rime_n = pi/4.0 * e_coll_n * n_g * n_c * dt & 
                 &   * (delta_n_gg * D_g*D_g + delta_n_gc * D_g*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_gg * v_g*v_g - theta_n_gc * v_g*v_c + theta_n_cc * v_c*v_c)

            rime_q = pi/4.0 * e_coll_q * n_g * q_c * dt & 
                 &   * (delta_q_gg * D_g*D_g + delta_q_gc * D_g*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_gg * v_g*v_g - theta_q_gc * v_g*v_c + theta_q_cc * v_c*v_c)

            rime_q = MIN(q_c,rime_q)
!ub            rime_n = MIN(n_c,MIN(rime_n,rime_q/x_coll_c))
            rime_n = MIN(n_c,rime_n)

            riming_q(i,j,k) = rime_q

            q_graupel(i,j,k) = q_graupel(i,j,k) + rime_q
            q_cloud(i,j,k)   = q_cloud(i,j,k)   - rime_q
            n_cloud(i,j,k)   = n_cloud(i,j,k)   - rime_n
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,16) = rime_q
            END IF
#endif
            ! ub<<

            ! Eismultiplikation nach Hallet und Mossop

            mult_q = 0.0
            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = const2*(T_a - T_mult_min) 
              mult_2 = const3*(T_a - T_mult_max) 
              mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
              mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,j,k)     = n_ice(i,j,k)     + mult_n
              q_ice(i,j,k)     = q_ice(i,j,k)     + mult_q
              q_graupel(i,j,k) = q_graupel(i,j,k) - mult_q

            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,17) = mult_q
            END IF
#endif
            ! ub<<

            ! enhancement of melting of graupel

            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = const4 * (T_a - T_3) * rime_q
              melt_n = melt_q / x_g

              melt_q = MIN(q_graupel(i,j,k),melt_q)
              melt_n = MIN(n_graupel(i,j,k),melt_n)

              rate_q(i,j,k) = melt_q

              q_graupel(i,j,k) = q_graupel(i,j,k) - melt_q
              q_rain(i,j,k)    = q_rain(i,j,k)    + melt_q

              n_graupel(i,j,k) = n_graupel(i,j,k) - melt_n
              n_rain(i,j,k)    = n_rain(i,j,k)    + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,18) = melt_q
            END IF
#endif
            ! ub<<

            ! Shedding

            IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
              !IF (graupel_shedding .AND. T_a > T_shed ) THEN
              q_g = q_graupel(i,j,k)
              n_g = n_graupel(i,j,k)
              x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI

              shed_q = MIN(q_g,rime_q)
              shed_n = shed_q / MIN(x_shed,x_g)

              q_graupel(i,j,k) = q_graupel(i,j,k) - shed_q
              q_rain(i,j,k)    = q_rain(i,j,k)    + shed_q                   
              n_rain(i,j,k)    = n_rain(i,j,k)    + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,19) = shed_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(riming_n,riming_q,rate_q)

  END SUBROUTINE graupel_cloud_riming

  ! hn >>
  SUBROUTINE hail_cloud_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Hagelpartikel                             *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g,        &
         &                        q_hail, n_hail, q_cloud, n_cloud, q_rain, n_rain,  &
         &                        q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h
    DOUBLE PRECISION            :: q_c,n_c,x_c,d_c,v_c,x_coll_c
    DOUBLE PRECISION            :: rime_n,rime_q,e_coll_n
    DOUBLE PRECISION            :: melt_n,melt_q,e_coll_q
    DOUBLE PRECISION            :: shed_n,shed_q,x_shed
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_hh,delta_n_hc,delta_n_cc
    DOUBLE PRECISION, SAVE      :: delta_q_hh,delta_q_hc,delta_q_cc
    DOUBLE PRECISION, SAVE      :: theta_n_hh,theta_n_hc,theta_n_cc
    DOUBLE PRECISION, SAVE      :: theta_q_hh,theta_q_hc,theta_q_cc
    DOUBLE PRECISION            :: const1,const2,const3,const4

    DOUBLE PRECISION, PARAMETER :: eps = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: riming_n,riming_q,rate_q

    ALLOCATE(riming_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(riming_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD hail_cloud_riming",4)

    IF (isIO() .AND. isdebug) THEN
      IF (enhanced_melting .AND. ice_multiplication) THEN
        WRITE (6, *) "CLOUD hail_cloud_riming mit ice_multiplication und enhanced_melting"
      ELSEIF (enhanced_melting .AND. .NOT. ice_multiplication) THEN
        WRITE (6, *) "CLOUD hail_cloud_riming mit enhanced_melting"
      ELSEIF (ice_multiplication .AND. .NOT. enhanced_melting) THEN
        WRITE (6, *) "CLOUD hail_cloud_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD hail_cloud_riming" 
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,cloud,0)
      delta_n_hc = coll_delta_12(hail,cloud,0)
      delta_n_cc = coll_delta_22(hail,cloud,0)
      delta_q_hh = coll_delta_11(hail,cloud,0) 
      delta_q_hc = coll_delta_12(hail,cloud,1)
      delta_q_cc = coll_delta_22(hail,cloud,1)

      theta_n_hh = coll_theta_11(hail,cloud,0)
      theta_n_hc = coll_theta_12(hail,cloud,0)
      theta_n_cc = coll_theta_22(hail,cloud,0)
      theta_q_hh = coll_theta_11(hail,cloud,0)
      theta_q_hc = coll_theta_12(hail,cloud,1)
      theta_q_cc = coll_theta_22(hail,cloud,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hagel  :   a_hail      = ",hail%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_hail      = ",hail%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Hagel:         alf_hail    = ",hail%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_hail    = ",hail%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentr.: a_cloud     = ",cloud%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_cloud     = ",cloud%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Wolkentr.:     alf_cloud   = ",cloud%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_cloud   = ",cloud%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming hail-cloud:       delta_n_hh  = ",delta_n_hh
        WRITE (6,'(A,D10.3)') "                                       delta_n_hc  = ",delta_n_hc
        WRITE (6,'(A,D10.3)') "                                       delta_n_cc  = ",delta_n_cc
        WRITE (6,'(A,D10.3)') "                                       theta_n_hh  = ",theta_n_hh
        WRITE (6,'(A,D10.3)') "                                       theta_n_hc  = ",theta_n_hc
        WRITE (6,'(A,D10.3)') "                                       theta_n_cc  = ",theta_n_cc
        WRITE (6,'(A,D10.3)') "                                       delta_q_hh  = ",delta_q_hh
        WRITE (6,'(A,D10.3)') "                                       delta_q_hc  = ",delta_q_hc
        WRITE (6,'(A,D10.3)') "                                       delta_q_cc  = ",delta_q_cc
        WRITE (6,'(A,D10.3)') "                                       theta_q_hh  = ",theta_q_hh
        WRITE (6,'(A,D10.3)') "                                       theta_q_hc  = ",theta_q_hc
        WRITE (6,'(A,D10.3)') "                                       theta_q_cc  = ",theta_q_cc
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3.*pi * rho_w * r_shedding**3     !..Mittlere Masse der Sheddingtropfen

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..Mindestmasse fuer collection, begrenzt rime_n

    const1 = e_hc/(D_coll_c - D_krit_c)    
    const2 = 1.0/(T_mult_opt - T_mult_min)
    const3 = 1.0/(T_mult_opt - T_mult_max)
    const4 = c_w / L_ew
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix             
          q_c = q_cloud(i,j,k)                                      !..Massendichte
          q_h = q_hail(i,j,k)                                    !..Massendichte
          n_c = n_cloud(i,j,k)                                      !..Anzahldichte
          n_h = n_hail(i,j,k)                                    !..Anzahldichte
          x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse
          D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
          x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)     !..mittlere Masse
          D_c = cloud%a_geo * x_c**cloud%b_geo                      !..mittlerer Durchmesser
          T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)                  !..abs. Temperatur

          IF (q_c > q_krit_c .AND. q_h > q_krit_hc .AND. D_h > D_krit_hc .AND. &
               & D_c > D_krit_c) THEN

            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            v_c = cloud%a_vel   * x_c**cloud%b_vel   * rrho_c(i,j,k)  !..mittlere Sedimentationsgeschw.

            e_coll_n = MIN(e_hc, MAX(const1*(D_c - D_krit_c),e_min))
            e_coll_q = e_coll_n

            rime_n = pi/4.0 * e_coll_n * n_h * n_c * dt & 
                 &   * (delta_n_hh * D_h*D_h + delta_n_hc * D_h*D_c + delta_n_cc * D_c*D_c) &
                 &   * SQRT(theta_n_hh * v_h*v_h - theta_n_hc * v_h*v_c + theta_n_cc * v_c*v_c)

            rime_q = pi/4.0 * e_coll_q * n_h * q_c * dt & 
                 &   * (delta_q_hh * D_h*D_h + delta_q_hc * D_h*D_c + delta_q_cc * D_c*D_c) &
                 &   * SQRT(theta_q_hh * v_h*v_h - theta_q_hc * v_h*v_c + theta_q_cc * v_c*v_c)

            rime_q = MIN(q_c,rime_q)
!ub            rime_n = MIN(n_c,MIN(rime_n,rime_q/x_coll_c))
            rime_n = MIN(n_c,rime_n)

            riming_q(i,j,k) = rime_q

            q_hail(i,j,k)  = q_hail(i,j,k) + rime_q
            q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
            n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,44) = rime_q
            END IF
#endif
            ! ub<<

            ! Eismultiplikation nach Hallet und Mossop

            mult_q = 0.0
            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = const2*(T_a - T_mult_min) 
              mult_2 = const3*(T_a - T_mult_max) 
              mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
              mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
              q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
              q_hail(i,j,k) = q_hail(i,j,k) - mult_q

            ENDIF
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,45) = mult_q
            END IF
#endif
            ! ub<<

            ! enhancement of melting of hail

            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = const4 * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(q_hail(i,j,k),melt_q)
              melt_n = MIN(n_hail(i,j,k),melt_n)

              rate_q(i,j,k) = melt_q

              q_hail(i,j,k) = q_hail(i,j,k) - melt_q
              q_rain(i,j,k) = q_rain(i,j,k) + melt_q

              n_hail(i,j,k) = n_hail(i,j,k) - melt_n
              n_rain(i,j,k) = n_rain(i,j,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,46) = melt_q
            END IF
#endif
            ! ub<<

            ! Shedding

            IF ((D_h > D_shed_h .AND. T_a > T_shed .AND. hail_shedding) .OR. T_a > T_3 ) THEN
              q_h = q_hail(i,j,k)
              n_h = n_hail(i,j,k)
              x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse in SI

              shed_q = MIN(q_h,rime_q)
              shed_n = shed_q / MIN(x_shed,x_h)

              q_hail(i,j,k) = q_hail(i,j,k) - shed_q
              q_rain(i,j,k) = q_rain(i,j,k) + shed_q                   
              n_rain(i,j,k) = n_rain(i,j,k) + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,47) = shed_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(riming_n,riming_q,rate_q)

  END SUBROUTINE hail_cloud_riming

  DOUBLE PRECISION FUNCTION coll_delta(p1,n)
    IMPLICIT NONE
    TYPE(PARTICLE) :: p1
    INTEGER        :: n

    coll_delta = gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)         &
         &                  / gfct((p1%nu+1.0  )/p1%mu)         &
         &        * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_geo+n)   &
         &        / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_geo+n)

    RETURN
  END FUNCTION coll_delta

  DOUBLE PRECISION FUNCTION coll_delta_11(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_delta_11 = coll_delta(p1,n)

    RETURN
  END FUNCTION coll_delta_11

  DOUBLE PRECISION FUNCTION coll_delta_22(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_delta_22 = coll_delta(p2,n)

    RETURN
  END FUNCTION coll_delta_22

  DOUBLE PRECISION FUNCTION coll_delta_12(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_delta_12 = 2.0 * gfct((p1%b_geo+p1%nu+1.0)/p1%mu)               &
         &                       / gfct((p1%nu+1.0)/p1%mu)               &
         &                * gfct((p1%nu+1.0)/p1%mu)**(p1%b_geo)          &
         &                / gfct((p1%nu+2.0)/p1%mu)**(p1%b_geo)          &
         &              * gfct((p2%b_geo+p2%nu+1.0+n)/p2%mu)             &
         &                        /gfct((p2%nu+1.0  )/p2%mu)             &
         &                * gfct((p2%nu+1.0)/p2%mu)**(p2%b_geo+n)        &
         &                / gfct((p2%nu+2.0)/p2%mu)**(p2%b_geo+n)

    RETURN
  END FUNCTION coll_delta_12

  DOUBLE PRECISION FUNCTION coll_theta(p1,n)

    TYPE(PARTICLE) :: p1
    INTEGER        :: n

    coll_theta = gfct((2.0*p1%b_vel+2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &               / gfct((2.0*p1%b_geo+p1%nu+1.0+n)/p1%mu)    &
         &          * gfct((p1%nu+1.0)/p1%mu)**(2.0*p1%b_vel)        &
         &          / gfct((p1%nu+2.0)/p1%mu)**(2.0*p1%b_vel)

    RETURN
  END FUNCTION coll_theta

  DOUBLE PRECISION FUNCTION coll_theta_11(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_theta_11 = coll_theta(p1,n)

    RETURN
  END FUNCTION coll_theta_11

  DOUBLE PRECISION FUNCTION coll_theta_22(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_theta_22 = coll_theta(p2,n)

    RETURN
  END FUNCTION coll_theta_22

  DOUBLE PRECISION FUNCTION coll_theta_12(p1,p2,n)

    TYPE(PARTICLE) :: p1,p2
    INTEGER        :: n

    coll_theta_12 = 2.0 * gfct((p1%b_vel+2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &  
         &                    / gfct((2.0*p1%b_geo+p1%nu+1.0)/p1%mu)       &
         &              * gfct((p1%nu+1.0)/p1%mu)**(p1%b_vel)              &
         &              / gfct((p1%nu+2.0)/p1%mu)**(p1%b_vel)              &
         &           * gfct((p2%b_vel+2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &  
         &                    / gfct((2.0*p2%b_geo+p2%nu+1.0+n)/p2%mu)     &
         &              * gfct((p2%nu+1.0)/p2%mu)**(p2%b_vel)              &
         &              / gfct((p2%nu+2.0)/p2%mu)**(p2%b_vel)

    RETURN
  END FUNCTION coll_theta_12

  ! Faktor zur Berechnung des mittleren Durchmessers von verallg. gammaverteilten Hydrometeoren:
  FUNCTION D_average_factor (parti)

    ! gueltig fuer D = a_geo * x^b_geo
    ! Berechnung des mittleren Durchmessers: D_average = parti%b_geo * D_average_factor * (q/qn)**parti%b_geo
  
    IMPLICIT NONE
    
    DOUBLE PRECISION :: D_average_factor
    TYPE(PARTICLE), INTENT(in) :: parti
    
    D_average_factor = &
         ( gfct( (parti%b_geo+parti%nu+1.0)/parti%mu ) / &
         gfct( (parti%nu+1.0)/parti%mu ) ) * &
         (gfct( (parti%nu+1.0)/parti%mu ) / gfct( (parti%nu+2.0)/parti%mu ) ) ** parti%b_geo
    
  END FUNCTION D_average_factor

  SUBROUTINE ice_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Eispartikel                               *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_ice, n_ice, q_rain, n_rain, q_hail, n_hail, &
         &                        q_graupel, n_graupel, dt, dqdt, speichere_dqdt, &
         &                        rimenrrate_ice, rimeqirate_ice, rimeqrrate_ice, &
         &                        d_id_sp, d_rd_sp_ice
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i,d_id
    DOUBLE PRECISION            :: q_r,n_r,x_r,d_r,v_r,d_rd
    DOUBLE PRECISION            :: rime_n,rime_qi,rime_qr
    DOUBLE PRECISION            :: conv_n,conv_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_ii,delta_n_ir,           delta_n_rr
    DOUBLE PRECISION, SAVE      :: delta_q_ii,delta_q_ir,delta_q_ri,delta_q_rr
    DOUBLE PRECISION, SAVE      :: theta_n_ii,theta_n_ir,           theta_n_rr
    DOUBLE PRECISION, SAVE      :: theta_q_ii,theta_q_ir,theta_q_ri,theta_q_rr,D_av_fakt_i,D_av_fakt_r
    DOUBLE PRECISION            :: const3,const4

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      IF (ice_multiplication) THEN
        WRITE (6, *) "CLOUD ice_rain_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD ice_rain_riming "
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_ii = coll_delta_11(ice,rain,0)
      delta_n_ir = coll_delta_12(ice,rain,0)
      delta_n_rr = coll_delta_22(ice,rain,0)
      delta_q_ii = coll_delta_11(ice,rain,1) 
      delta_q_ir = coll_delta_12(ice,rain,1)
      delta_q_ri = coll_delta_12(rain,ice,1)
      delta_q_rr = coll_delta_22(ice,rain,1)

      theta_n_ii = coll_theta_11(ice,rain,0)
      theta_n_ir = coll_theta_12(ice,rain,0)
      theta_n_rr = coll_theta_22(ice,rain,0)
      theta_q_ii = coll_theta_11(ice,rain,1)
      theta_q_ir = coll_theta_12(ice,rain,1)
      theta_q_ri = coll_theta_12(rain,ice,1)
      theta_q_rr = coll_theta_22(ice,rain,1)

      D_av_fakt_r = D_average_factor(rain)
      D_av_fakt_i = D_average_factor(ice)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:        a_ice      = ",ice%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_ice      = ",ice%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Eis:            alf_ice    = ",ice%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_ice    = ",ice%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Regentr..:  a_rain    = ",rain%a_geo
        WRITE (6,'(A,D10.3)') "                                        b_rain    = ",rain%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Regentr.:       alf_rain  = ",rain%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_rain  = ",rain%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming ice-rain:          delta_n_ii = ",delta_n_ii
        WRITE (6,'(A,D10.3)') "                                        delta_n_ir = ",delta_n_ir
        WRITE (6,'(A,D10.3)') "                                        delta_n_rr = ",delta_n_rr
        WRITE (6,'(A,D10.3)') "                                        theta_n_ii = ",theta_n_ii
        WRITE (6,'(A,D10.3)') "                                        theta_n_ir = ",theta_n_ir
        WRITE (6,'(A,D10.3)') "                                        theta_n_rr = ",theta_n_rr
        WRITE (6,'(A,D10.3)') "                                        delta_q_ii = ",delta_q_ii
        WRITE (6,'(A,D10.3)') "                                        delta_q_ir = ",delta_q_ir
        WRITE (6,'(A,D10.3)') "                                        delta_q_ri = ",delta_q_ri
        WRITE (6,'(A,D10.3)') "                                        delta_q_rr = ",delta_q_rr
        WRITE (6,'(A,D10.3)') "                                        theta_q_ii = ",theta_q_ii
        WRITE (6,'(A,D10.3)') "                                        theta_q_ir = ",theta_q_ir
        WRITE (6,'(A,D10.3)') "                                        theta_q_ri = ",theta_q_ri
        WRITE (6,'(A,D10.3)') "                                        theta_q_rr = ",theta_q_rr
        WRITE (6,'(A,D10.3)') "                                        D_av_fakt_r = ",D_av_fakt_r
        WRITE (6,'(A,D10.3)') "                                        D_av_fakt_i = ",D_av_fakt_i
      END IF
      firstcall = 1
    ENDIF

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                              !..Fluessigwassergehalt in SI
          q_i = q_ice(i,j,k)                               !..Fluessigwassergehalt in SI
          n_r = n_rain(i,j,k)                              !..Anzahldichte in SI
          n_i = n_ice(i,j,k)                               !..Anzahldichte in SI

          T_a = T_0(i,j,k) !WRF! + T(i,j,k) + T_g(i,j,k)

          x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse in SI     
          D_i = ice%a_geo * x_i**ice%b_geo                     !..mittlerer Durchmesser
          D_id = D_av_fakt_i * D_i 

          IF (q_r > q_krit .AND. q_i > q_krit_ir .AND. D_i > D_krit_ir) THEN

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)  !..mittlere Masse in SI     

            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

            D_r = rain%a_geo * x_r**rain%b_geo                   !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.
            D_rd = D_av_fakt_r * D_r

            rime_n  = pi/4.0 * n_i * n_r * dt & 
                 &   * (delta_n_ii * D_i*D_i + delta_n_ir * D_i*D_r + delta_n_rr * D_r*D_r) &
                 &   * SQRT(theta_n_ii * v_i*v_i - theta_n_ir * v_i*v_r + theta_n_rr * v_r*v_r  &
                 &     +ice_s_vel**2)

            rime_qr = pi/4.0 * n_i * q_r * dt & 
                 &   * (delta_n_ii * D_i*D_i + delta_q_ir * D_i*D_r + delta_q_rr * D_r*D_r) &
                 &   * SQRT(theta_n_ii * v_i*v_i - theta_q_ir * v_i*v_r + theta_q_rr * v_r*v_r  &
                 &     +ice_s_vel**2)

            rime_qi = pi/4.0 * n_r * q_i * dt & 
                 &   * (delta_q_ii * D_i*D_i + delta_q_ri * D_i*D_r + delta_n_rr * D_r*D_r) &
                 &   * SQRT(theta_q_ii * v_i*v_i - theta_q_ri * v_i*v_r + theta_n_rr * v_r*v_r  &
                 &     +ice_s_vel**2)

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_n  = MIN(MIN(n_r,n_i),rime_n)
              rime_qr = MIN(q_r,rime_qr)
              rime_qi = MIN(q_i,rime_qi)

              n_ice(i,j,k)  = n_ice(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_ice(i,j,k)  = q_ice(i,j,k)  - rime_qi
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0  ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_id > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_ice(i,j,k) = n_ice(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_ice(i,j,k) = q_ice(i,j,k) + rime_qi     ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Eisteilchen schmelzen instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qi   ! UB_20081120 - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,29) = rime_qi
                  END IF
#endif
                END IF
              ELSE
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n  ! UB_20081120
                q_ice(i,j,k) = q_ice(i,j,k)  + mult_q  ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Eis + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,21) = mult_q
                    dqdt(i,j,k,20) = rime_qi + rime_qr
                  END IF
#endif
                ELSE
                  IF (D_id > D_rd) THEN
                    ! Eis + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,21) = mult_q
                      dqdt(i,j,k,20) = rime_qi + rime_qr
                    END IF
#endif
                  ELSE
                    ! Eis + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,21) = mult_q
                      dqdt(i,j,k,61) = rime_qi + rime_qr
                    END IF
#endif
                  END IF
                END IF
              END IF

            ELSE

              rimenrrate_ice(i,j,k) = rimenrrate_ice(i,j,k) + rime_n
              rimeqirate_ice(i,j,k) = rimeqirate_ice(i,j,k) + rime_qi
              rimeqrrate_ice(i,j,k) = rimeqrrate_ice(i,j,k) + rime_qr
              d_id_sp(i,j,k)        = D_id
              d_rd_sp_ice(i,j,k)    = D_rd

            END IF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_rain_riming

  SUBROUTINE snow_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Schneepartikel                            *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_snow, n_snow, q_rain, n_rain, q_hail, n_hail, &
         &                        q_graupel, n_graupel, n_ice, q_ice, dt, dqdt, speichere_dqdt, &
         &                        rimenrrate_snow, rimeqirate_snow, rimeqrrate_snow, &
         &                        d_sd_sp, d_rd_sp_snow
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s,d_sd
    DOUBLE PRECISION            :: q_r,n_r,x_r,d_r,v_r,d_rd
    DOUBLE PRECISION            :: rime_n,rime_qr,rime_qs
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_ss,delta_n_sr,delta_n_rr
    DOUBLE PRECISION, SAVE      :: delta_q_ss,delta_q_sr,delta_q_rs,delta_q_rr
    DOUBLE PRECISION, SAVE      :: theta_n_ss,theta_n_sr,theta_n_rr
    DOUBLE PRECISION, SAVE      :: theta_q_ss,theta_q_sr,theta_q_rs,theta_q_rr,D_av_fakt_s,D_av_fakt_r
    DOUBLE PRECISION            :: const3,const4

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      IF (ice_multiplication) THEN
        WRITE (6, *) "CLOUD snow_rain_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD snow_rain_riming "
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,rain,0)
      delta_n_sr = coll_delta_12(snow,rain,0)
      delta_n_rr = coll_delta_22(snow,rain,0)
      delta_q_ss = coll_delta_11(snow,rain,1) 
      delta_q_sr = coll_delta_12(snow,rain,1)
      delta_q_rs = coll_delta_12(rain,snow,1)
      delta_q_rr = coll_delta_22(snow,rain,1)

      theta_n_ss = coll_theta_11(snow,rain,0)
      theta_n_sr = coll_theta_12(snow,rain,0)
      theta_n_rr = coll_theta_22(snow,rain,0)
      theta_q_ss = coll_theta_11(snow,rain,1)
      theta_q_sr = coll_theta_12(snow,rain,1)
      theta_q_rs = coll_theta_12(rain,snow,1)
      theta_q_rr = coll_theta_22(snow,rain,1)

      D_av_fakt_r = D_average_factor(rain)
      D_av_fakt_s = D_average_factor(snow)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:        a_snow      = ",snow%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_snow      = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Eis:            alf_snow    = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_snow    = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentr.:  a_rain    = ",rain%a_geo
        WRITE (6,'(A,D10.3)') "                                        b_rain    = ",rain%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Wolkentr.:      alf_rain  = ",rain%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_rain  = ",rain%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming snow-rain:         delta_n_ss = ",delta_n_ss
        WRITE (6,'(A,D10.3)') "                                        delta_n_sr = ",delta_n_sr
        WRITE (6,'(A,D10.3)') "                                        delta_n_rr = ",delta_n_rr
        WRITE (6,'(A,D10.3)') "                                        theta_n_ss = ",theta_n_ss
        WRITE (6,'(A,D10.3)') "                                        theta_n_sr = ",theta_n_sr
        WRITE (6,'(A,D10.3)') "                                        theta_n_rr = ",theta_n_rr
        WRITE (6,'(A,D10.3)') "                                        delta_q_ss = ",delta_q_ss
        WRITE (6,'(A,D10.3)') "                                        delta_q_sr = ",delta_q_sr
        WRITE (6,'(A,D10.3)') "                                        delta_q_rs = ",delta_q_rs
        WRITE (6,'(A,D10.3)') "                                        delta_q_rr = ",delta_q_rr
        WRITE (6,'(A,D10.3)') "                                        theta_q_ss = ",theta_q_ss
        WRITE (6,'(A,D10.3)') "                                        theta_q_sr = ",theta_q_sr
        WRITE (6,'(A,D10.3)') "                                        theta_q_rs = ",theta_q_rs
        WRITE (6,'(A,D10.3)') "                                        theta_q_rr = ",theta_q_rr
        WRITE (6,'(A,D10.3)') "                                        D_av_fakt_r = ",D_av_fakt_r
        WRITE (6,'(A,D10.3)') "                                        D_av_fakt_s = ",D_av_fakt_s
      END IF
      firstcall = 1
    ENDIF

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                              !..Fluessigwassergehalt in SI
          q_s = q_snow(i,j,k)                              !..Fluessigwassergehalt in SI
          n_r = n_rain(i,j,k)                              !..Anzahldichte in SI
          n_s = n_snow(i,j,k)                              !..Anzahldichte in SI

          x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max) !..mittlere Masse in SI     
          D_s = snow%a_geo * x_s**snow%b_geo                  !..mittlerer Durchmesser
          D_sd = D_av_fakt_s * D_s

          T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

          IF (q_r > q_krit .AND. q_s > q_krit_sr .AND. D_s > D_krit_sr) THEN

            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)  !..mittlere Masse in SI     
            D_r = rain%a_geo * x_r**rain%b_geo                   !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.
            D_rd = D_av_fakt_r * D_r

            rime_n = pi/4.0 * n_s * n_r * dt & 
                 &   * (delta_n_ss * D_s**2 + delta_n_sr * D_s*D_r + delta_n_rr * D_r**2) &
                 &   * (theta_n_ss * v_s**2 - theta_n_sr * v_s*v_r + theta_n_rr * v_r**2  &
                 &     +snow_s_vel**2)**0.5

            rime_qr = pi/4.0 * n_s * q_r * dt & 
                 &   * (delta_n_ss * D_s**2 + delta_q_sr * D_s*D_r + delta_q_rr * D_r**2) &
                 &   * (theta_n_ss * v_s**2 - theta_q_sr * v_s*v_r + theta_q_rr * v_r**2  &
                 &     +snow_s_vel**2)**0.5

            rime_qs = pi/4.0 * n_r * q_s * dt & 
                 &   * (delta_q_ss * D_s**2 + delta_q_rs * D_s*D_r + delta_n_rr * D_r**2) &
                 &   * (theta_q_ss * v_s**2 - theta_q_rs * v_s*v_r + theta_n_rr * v_r**2  &
                 &     +snow_s_vel**2)**0.5

            IF (.NOT.use_ice_graupel_conv_uli) THEN

              rime_qr = MIN(q_r,rime_qr)
              rime_qs = MIN(q_s,rime_qs)
              rime_n  = MIN(n_r,rime_n)
              rime_n  = MIN(n_s,rime_n)

              n_snow(i,j,k) = n_snow(i,j,k) - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_snow(i,j,k) = q_snow(i,j,k) - rime_qs
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0 ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_sd > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_snow(i,j,k) = n_snow(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_snow(i,j,k) = q_snow(i,j,k) + rime_qs ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Schneeflocken werden instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qs ! UB_20081120 - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,30) = rime_qs
                  END IF
#endif
                END IF
              ELSE
                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n ! UB_20081120
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,23) = mult_q
                    dqdt(i,j,k,22) = rime_qr + rime_qs
                  END IF
#endif
                ELSE
                  IF (D_sd > D_rd) THEN
                    ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,23) = mult_q
                      dqdt(i,j,k,22) = rime_qr + rime_qs
                    END IF
#endif
                  ELSE
                    ! Schnee + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,23) = mult_q
                      dqdt(i,j,k,57) = rime_qr + rime_qs
                  END IF
#endif
                  END IF
                END IF
              END IF


            ELSE

              rimenrrate_snow(i,j,k)  = rimenrrate_snow(i,j,k)  + rime_n
              rimeqirate_snow(i,j,k)  = rimeqirate_snow(i,j,k) + rime_qs
              rimeqrrate_snow(i,j,k)  = rimeqrrate_snow(i,j,k) + rime_qr
              d_sd_sp(i,j,k) = D_sd
              d_rd_sp_snow(i,j,k) = D_rd

            END IF

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_rain_riming

  SUBROUTINE graupel_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Graupelpartikel                           *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g,     &
         &                        q_graupel, n_graupel, q_rain, n_rain, q_hail, n_hail, &
         &                        q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g,d_h,x_h,d_gd
    DOUBLE PRECISION            :: q_r,n_r,x_r,d_r,v_r,x_shed,q_h,n_h,d_rd
    DOUBLE PRECISION            :: rime_n,rime_q,rime_qg,rime_qr,melt_n,melt_q,shed_n,shed_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_gg,delta_n_gr,delta_n_rr
    DOUBLE PRECISION, SAVE      :: delta_q_gg,delta_q_gr,delta_q_rg,delta_q_rr
    DOUBLE PRECISION, SAVE      :: theta_n_gg,theta_n_gr,theta_n_rr
    DOUBLE PRECISION, SAVE      :: theta_q_gg,theta_q_gr,theta_q_rg,theta_q_rr,D_av_fakt_g,D_av_fakt_r
    DOUBLE PRECISION            :: const3,const4

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF (isIO() .AND. isdebug) THEN
      IF (enhanced_melting .AND. ice_multiplication) THEN
        WRITE (6, *) "CLOUD graupel_rain_riming mit ice_multiplication und enhanced_melting"
      ELSEIF (enhanced_melting .AND. .NOT. ice_multiplication) THEN
        WRITE (6, *) "CLOUD graupel_rain_riming mit enhanced_melting"
      ELSEIF (ice_multiplication .AND. .NOT. enhanced_melting) THEN
        WRITE (6, *) "CLOUD graupel_rain_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD graupel_rain_riming" 
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,rain,0)
      delta_n_gr = coll_delta_12(graupel,rain,0)
      delta_n_rr = coll_delta_22(graupel,rain,0)
!      delta_q_gg = coll_delta_11(graupel,rain,0) ! UB_20081120
      delta_q_gg = coll_delta_11(graupel,rain,1)  ! UB_20081120
      delta_q_gr = coll_delta_12(graupel,rain,1)
      delta_q_rg = coll_delta_12(rain,graupel,1)
      delta_q_rr = coll_delta_22(graupel,rain,1)

      theta_n_gg = coll_theta_11(graupel,rain,0)
      theta_n_gr = coll_theta_12(graupel,rain,0)
      theta_n_rr = coll_theta_22(graupel,rain,0)
!      theta_q_gg = coll_theta_11(graupel,rain,0) ! UB_20081120
      theta_q_gg = coll_theta_11(graupel,rain,1)  ! UB_20081120
      theta_q_gr = coll_theta_12(graupel,rain,1)
      theta_q_rg = coll_theta_12(rain,graupel,1)
      theta_q_rr = coll_theta_22(graupel,rain,1)

      D_av_fakt_r = D_average_factor(rain)
      D_av_fakt_g = D_average_factor(graupel)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:   a_graupel   = ",graupel%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_graupel   = ",graupel%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Graupel:       alf_graupel = ",graupel%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_graupel = ",graupel%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Regentr.:  a_rain   = ",rain%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_rain   = ",rain%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Regentr.:      alf_rain = ",rain%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_rain = ",rain%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming graupel-rain:     delta_n_gg = ",delta_n_gg
        WRITE (6,'(A,D10.3)') "                                       delta_n_gr = ",delta_n_gr
        WRITE (6,'(A,D10.3)') "                                       delta_n_rr = ",delta_n_rr
        WRITE (6,'(A,D10.3)') "                                       theta_n_gg = ",theta_n_gg
        WRITE (6,'(A,D10.3)') "                                       theta_n_gr = ",theta_n_gr
        WRITE (6,'(A,D10.3)') "                                       theta_n_rr = ",theta_n_rr
        WRITE (6,'(A,D10.3)') "                                       delta_q_gg = ",delta_q_gg
        WRITE (6,'(A,D10.3)') "                                       delta_q_gr = ",delta_q_gr
        WRITE (6,'(A,D10.3)') "                                       delta_q_rg = ",delta_q_rg
        WRITE (6,'(A,D10.3)') "                                       delta_q_rr = ",delta_q_rr
        WRITE (6,'(A,D10.3)') "                                       theta_q_gg = ",theta_q_gg
        WRITE (6,'(A,D10.3)') "                                       theta_q_gr = ",theta_q_gr
        WRITE (6,'(A,D10.3)') "                                       theta_q_rg = ",theta_q_rg
        WRITE (6,'(A,D10.3)') "                                       theta_q_rr = ",theta_q_rr
        WRITE (6,'(A,D10.3)') "                                       D_av_fakt_r = ",D_av_fakt_r
        WRITE (6,'(A,D10.3)') "                                       D_av_fakt_g = ",D_av_fakt_g
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3. * pi * rho_w * r_shedding**3

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                                    !..Fluessigwassergehalt in SI
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI

          T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)
          IF (q_r > q_krit .AND. q_g > q_krit) THEN

            n_r = n_rain(i,j,k)                                 !..Anzahldichte in SI
            n_g = n_graupel(i,j,k)                              !..Anzahldichte in SI

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)       !..mittlere Masse in SI     
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI     

            D_r = rain%a_geo * x_r**rain%b_geo                        !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.
            D_rd = D_av_fakt_r * D_r

            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            D_gd = D_av_fakt_g * D_g 

            IF (T_a >= T_3 .OR. ice_typ < 3 .OR. D_gd > D_rd) THEN

              rime_n = pi/4.0 * n_g * n_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_n_gr * D_g*D_r + delta_n_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_n_gr * v_g*v_r + theta_n_rr * v_r**2)**0.5

              rime_q = pi/4.0 * n_g * q_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_q_gr * D_g*D_r + delta_q_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_q_gr * v_g*v_r + theta_q_rr * v_r**2)**0.5

              rime_n = MIN(n_r,rime_n)
              rime_q = MIN(q_r,rime_q)

              q_graupel(i,j,k) = q_graupel(i,j,k) + rime_q
              q_rain(i,j,k)    = q_rain(i,j,k)    - rime_q
              n_rain(i,j,k)    = n_rain(i,j,k)    - rime_n

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,24) = rime_q
              END IF
#endif
              ! ub<<

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)     = n_ice(i,j,k)     + mult_n
                q_ice(i,j,k)     = q_ice(i,j,k)     + mult_q
                q_graupel(i,j,k) = q_graupel(i,j,k) - mult_q
              ENDIF

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,25) = mult_q
              END IF
#endif
              ! ub<<

              ! enhancement of melting of graupel

              IF (T_a > T_3 .AND. enhanced_melting) THEN
                melt_q = c_w / L_ew * (T_a - T_3) * rime_q
                melt_n = melt_q / x_g

                melt_q = MIN(q_graupel(i,j,k),melt_q)
                melt_n = MIN(n_graupel(i,j,k),melt_n)

                q_graupel(i,j,k) = q_graupel(i,j,k) - melt_q
                q_rain(i,j,k)    = q_rain(i,j,k)    + melt_q

                n_graupel(i,j,k) = n_graupel(i,j,k) - melt_n
                n_rain(i,j,k)    = n_rain(i,j,k)    + melt_n
              ELSE
                melt_q = 0.0
              ENDIF

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,26) = melt_q
              END IF
#endif

              ! Shedding

              IF ((graupel_shedding .AND. D_g > D_shed_g .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
                q_g = q_graupel(i,j,k)
                n_g = n_graupel(i,j,k)
                x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI     

                shed_q = MIN(q_g,rime_q)
                IF (T_a <= T_3) THEN
                  shed_n = shed_q / MIN(x_shed,x_g)
                ELSE
                  shed_n = shed_q / MAX(x_r,x_g)
                ENDIF

                q_graupel(i,j,k) = q_graupel(i,j,k) - shed_q
                q_rain(i,j,k)    = q_rain(i,j,k)    + shed_q                   
                n_rain(i,j,k)    = n_rain(i,j,k)    + shed_n
              ELSE
                shed_q = 0.0
              ENDIF

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,27) = shed_q
              END IF
#endif


            ELSE   ! T_a < T_3 .and. ice_typ >= 3 .and. D_g <= D_r

              rime_n = pi/4.0 * n_g * n_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_n_gr * D_g*D_r + delta_n_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_n_gr * v_g*v_r + theta_n_rr * v_r**2)**0.5

              rime_qr = pi/4.0 * n_g * q_r * dt & 
                   &   * (delta_n_gg * D_g**2 + delta_q_gr * D_g*D_r + delta_q_rr * D_r**2) &
                   &   * (theta_n_gg * v_g**2 - theta_q_gr * v_g*v_r + theta_q_rr * v_r**2)**0.5

              rime_qg = pi/4.0 * n_r * q_g * dt & 
                   &   * (delta_q_gg * D_g**2 + delta_q_rg * D_g*D_r + delta_n_rr * D_r**2) &
                   &   * (theta_q_gg * v_g**2 - theta_q_rg * v_g*v_r + theta_n_rr * v_r**2)**0.5

              rime_n = MIN(MIN(n_r,n_g),rime_n)
              rime_qr = MIN(q_r,rime_qr)
              rime_qg = MIN(q_g,rime_qg)

              n_graupel(i,j,k)  = n_graupel(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_graupel(i,j,k)  = q_graupel(i,j,k)  - rime_qg
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr

              n_hail(i,j,k) = n_hail(i,j,k) + rime_n
              q_hail(i,j,k) = q_hail(i,j,k) + rime_qg + rime_qr

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,58) = rime_qg + rime_qr
              END IF
#endif
              ! ub<<

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

                n_ice(i,j,k)     = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)     = q_ice(i,j,k)  + mult_q
                q_hail(i,j,k)    = q_hail(i,j,k) - mult_q
              ENDIF

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,25) = mult_q
              END IF
#endif
              ! ub<<

              ! Shedding

              x_h = MIN(MAX(q_hail(i,j,k)/(n_hail(i,j,k)+eps),hail%x_min),hail%x_max)       !..mittlere Masse in SI     
              D_h = hail%a_geo * x_h**hail%b_geo                        !..mittlerer Durchmesser

              IF (hail_shedding .AND. D_h > D_shed_h .AND. T_a > T_shed) THEN
                ! ub                IF (graupel_shedding .AND. (D_g > D_shed_g .OR. T_a > T_shed) ) THEN
                q_h = q_hail(i,j,k)
                n_h = n_hail(i,j,k)

                ! Vorher wurde Graupel + Regen zu Hagel umgewandelt und n_h und q_h haben sich erhoeht.
                ! Durch Shedding soll der Hagel nun die durch den Regen hinzugekommene Masse wieder verlieren.
                ! graupel_rain_riming wirkt also hier als Katalysator zur Umwandlung von Graupel in Hagel und
                ! zur Umwandlung von Regen in kleine Regentropfen.
                ! Zu Ueberlegen waere, ob nicht in diesem Falle der Graupel Graupel bleiben soll!!!
                shed_q = MIN(q_h,rime_qr)
                shed_n = shed_q / MIN(x_shed,x_h)

                q_hail(i,j,k)    = q_hail(i,j,k)    - shed_q
                q_rain(i,j,k)    = q_rain(i,j,k)    + shed_q                   
                n_rain(i,j,k)    = n_rain(i,j,k)    + shed_n
              ELSE
                shed_q = 0.0
              ENDIF

              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,27) = shed_q
              END IF
#endif
              ! ub<<

            ENDIF  ! T_a >= T_3 .or. ice_typ < 3
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_rain_riming

  ! hn >>
  SUBROUTINE hail_rain_riming()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Bereifens der Hagelpartikel                             *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g,     &
         &                        q_hail, n_hail, q_rain, n_rain, &
         &                        q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h
    DOUBLE PRECISION            :: q_r,n_r,x_r,d_r,v_r,x_shed
    DOUBLE PRECISION            :: rime_n,rime_q,melt_n,melt_q,shed_n,shed_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION, SAVE      :: delta_n_hh,delta_n_hr,delta_n_rr
    DOUBLE PRECISION, SAVE      :: delta_q_hh,delta_q_hr,delta_q_rr
    DOUBLE PRECISION, SAVE      :: theta_n_hh,theta_n_hr,theta_n_rr
    DOUBLE PRECISION, SAVE      :: theta_q_hh,theta_q_hr,theta_q_rr
    DOUBLE PRECISION            :: const3,const4

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF (isIO() .AND. isdebug) THEN
      IF (enhanced_melting .AND. ice_multiplication) THEN
        WRITE (6, *) "CLOUD hail_rain_riming mit ice_multiplication und enhanced_melting"
      ELSEIF (enhanced_melting .AND. .NOT. ice_multiplication) THEN
        WRITE (6, *) "CLOUD hail_rain_riming mit enhanced_melting"
      ELSEIF (ice_multiplication .AND. .NOT. enhanced_melting) THEN
        WRITE (6, *) "CLOUD hail_rain_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD hail_rain_riming" 
      ENDIF
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,rain,0)
      delta_n_hr = coll_delta_12(hail,rain,0)
      delta_n_rr = coll_delta_22(hail,rain,0)
      delta_q_hh = coll_delta_11(hail,rain,0) 
      delta_q_hr = coll_delta_12(hail,rain,1)
      delta_q_rr = coll_delta_22(hail,rain,1)

      theta_n_hh = coll_theta_11(hail,rain,0)
      theta_n_hr = coll_theta_12(hail,rain,0)
      theta_n_rr = coll_theta_22(hail,rain,0)
      theta_q_hh = coll_theta_11(hail,rain,0)
      theta_q_hr = coll_theta_12(hail,rain,1)
      theta_q_rr = coll_theta_22(hail,rain,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hagel:     a_hail     = ",hail%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_hail     = ",hail%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Hagel:         alf_hail   = ",hail%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_hail   = ",hail%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Regentr.:  a_rain     = ",rain%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_rain     = ",rain%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Regentr.:      alf_rain   = ",rain%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_rain   = ",rain%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer riming hail-rain:        delta_n_hh = ",delta_n_hh
        WRITE (6,'(A,D10.3)') "                                       delta_n_hr = ",delta_n_hr
        WRITE (6,'(A,D10.3)') "                                       delta_n_rr = ",delta_n_rr
        WRITE (6,'(A,D10.3)') "                                       theta_n_hh = ",theta_n_hh
        WRITE (6,'(A,D10.3)') "                                       theta_n_hr = ",theta_n_hr
        WRITE (6,'(A,D10.3)') "                                       theta_n_rr = ",theta_n_rr
        WRITE (6,'(A,D10.3)') "                                       delta_q_hh = ",delta_q_hh
        WRITE (6,'(A,D10.3)') "                                       delta_q_hr = ",delta_q_hr
        WRITE (6,'(A,D10.3)') "                                       delta_q_rr = ",delta_q_rr
        WRITE (6,'(A,D10.3)') "                                       theta_q_hh = ",theta_q_hh
        WRITE (6,'(A,D10.3)') "                                       theta_q_hr = ",theta_q_hr
        WRITE (6,'(A,D10.3)') "                                       theta_q_rr = ",theta_q_rr
      END IF
      firstcall = 1
    ENDIF

    x_shed = 4./3. * pi * rho_w * r_shedding**3

    const3 = 1/(T_mult_opt - T_mult_min)
    const4 = 1/(T_mult_opt - T_mult_max)

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)                                    !..Fluessigwassergehalt in SI
          q_h = q_hail(i,j,k)                                    !..Fluessigwassergehalt in SI

          T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)
          IF (q_r > q_krit .AND. q_h > q_krit) THEN
            n_r = n_rain(i,j,k)                                 !..Anzahldichte in SI
            n_h = n_hail(i,j,k)                                 !..Anzahldichte in SI

            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)       !..mittlere Masse in SI     
            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)       !..mittlere Masse in SI     

            D_r = rain%a_geo * x_r**rain%b_geo                        !..mittlerer Durchmesser
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.

            D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.


            rime_n = pi/4.0 * n_h * n_r * dt & 
                 &   * (delta_n_hh * D_h**2 + delta_n_hr * D_h*D_r + delta_n_rr * D_r**2) &
                 &   * (theta_n_hh * v_h**2 - theta_n_hr * v_h*v_r + theta_n_rr * v_r**2)**0.5

            rime_q = pi/4.0 * n_h * q_r * dt & 
                 &   * (delta_q_hh * D_h**2 + delta_q_hr * D_h*D_r + delta_q_rr * D_r**2) &
                 &   * (theta_q_hh * v_h**2 - theta_q_hr * v_h*v_r + theta_q_rr * v_r**2)**0.5

            rime_n = MIN(n_r,rime_n)
            rime_q = MIN(q_r,rime_q)

            q_hail(i,j,k) = q_hail(i,j,k) + rime_q
            q_rain(i,j,k) = q_rain(i,j,k) - rime_q
            n_rain(i,j,k) = n_rain(i,j,k) - rime_n
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,48) = rime_q
            END IF
#endif
            ! ub<<

            ! Eismultiplikation nach Hallet und Mossop

            mult_q = 0.0
            IF (T_a < T_3 .AND. ice_multiplication) THEN
              mult_1 = (T_a - T_mult_min) * const3
              mult_2 = (T_a - T_mult_max) * const4
              mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
              mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
              mult_n = C_mult * mult_1 * mult_2 * rime_q
              mult_q = mult_n * ice%x_min
              mult_q = MIN(rime_q,mult_q)

              n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
              q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
              q_hail(i,j,k) = q_hail(i,j,k) - mult_q
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,49) = mult_q
            END IF
#endif
            ! ub<<

            ! enhancement of melting of hail

            IF (T_a > T_3 .AND. enhanced_melting) THEN
              melt_q = c_w / L_ew * (T_a - T_3) * rime_q
              melt_n = melt_q / x_h

              melt_q = MIN(q_hail(i,j,k),melt_q)
              melt_n = MIN(n_hail(i,j,k),melt_n)

              q_hail(i,j,k) = q_hail(i,j,k) - melt_q
              q_rain(i,j,k) = q_rain(i,j,k) + melt_q

              n_hail(i,j,k) = n_hail(i,j,k) - melt_n
              n_rain(i,j,k) = n_rain(i,j,k) + melt_n
            ELSE
              melt_q = 0.0
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,50) = melt_q
            END IF
#endif
            ! ub<<

            ! Shedding

            IF ((hail_shedding .AND. D_h > D_shed_h .AND. T_a > T_shed) .OR. T_a > T_3 ) THEN
              ! ub                IF (hail_shedding .AND. (D_h > D_shed_h .OR. T_a > T_shed) ) THEN
              q_h = q_hail(i,j,k)
              n_h = n_hail(i,j,k)
              x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse in SI     

              shed_q = MIN(q_h,rime_q)
              IF (T_a <= T_3) THEN
                shed_n = shed_q / MIN(x_shed,x_h)
              ELSE
                !                shed_n = shed_q / x_shed
                shed_n = shed_q / MAX(x_r,x_h)
              ENDIF

              q_hail(i,j,k) = q_hail(i,j,k) - shed_q
              q_rain(i,j,k) = q_rain(i,j,k)    + shed_q                   
              n_rain(i,j,k) = n_rain(i,j,k)    + shed_n
            ELSE
              shed_q = 0.0
            ENDIF

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,51) = shed_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_rain_riming
  ! << hn

  ! ub>>
  SUBROUTINE complete_ice_snow_riming ()
    !*******************************************************************************
    !                                                                              *
    !       Im Falle von use_ice_graupel_conv_uli = .true. :                       *
    !                                                                              *
    !       Summierung der vorher lediglich gespeicherten riming-Raten             *
    !       auf die entsprechenden Felder (ice, snow, cloud, rain)                 *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_ice, n_ice, q_snow, n_snow, q_cloud, n_cloud, q_rain, n_rain, &
         &                        q_graupel, n_graupel, q_hail, n_hail, &
         &                        dt, dqdt, speichere_dqdt, &
         &                        deprate_ice, deprate_snow, &
         &                        rimeqcrate_ice, rimencrate_ice, &
         &                        rimeqirate_ice, rimeqrrate_ice, rimenrrate_ice, &
         &                        rimeqcrate_snow, rimencrate_snow, &
         &                        rimeqirate_snow, rimeqrrate_snow, rimenrrate_snow, &
         &                        d_id_sp, d_rd_sp_ice, d_sd_sp, d_rd_sp_snow

    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen
    INTEGER                     :: i,j,k

    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,D_id
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,D_sd
    DOUBLE PRECISION            :: x_coll_c,x_r
    DOUBLE PRECISION            :: D_rd
    DOUBLE PRECISION            :: rime_n,rime_q, rime_qr, rime_qs, rime_qi
    DOUBLE PRECISION            :: conv_n,conv_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION            :: const2,const3,const4,const5
    DOUBLE PRECISION, SAVE      :: D_av_fakt_i,D_av_fakt_s

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      IF (ice_multiplication) THEN
        WRITE (6, *) "CLOUD complete_ice_snow_riming mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD complete_ice_snow_riming ohne ice_multiplication"
      ENDIF
    END IF

    IF (firstcall .NE. 1) THEN
      D_av_fakt_i = D_average_factor(ice)
      D_av_fakt_s = D_average_factor(snow)
      firstcall = 1
    END IF

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..Mindestmasse fuer collection, begrenzt rime_n bei n_cloud

    const2 = 1.0/x_coll_c
    const3 = 1.0/(T_mult_opt - T_mult_min)
    const4 = 1.0/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a = T_0(i,j,k)

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !.. Vollenden des ice-rimings:
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          IF (deprate_ice(i,j,k) > 0.0d0 .AND. deprate_ice(i,j,k) >= rimeqcrate_ice(i,j,k)+rimeqrrate_ice(i,j,k)) THEN

            !.. Deposition ist groesser als gesamtes riming, deswegen bleibt Eis Eis:

            IF (rimeqcrate_ice(i,j,k) > 0.0d0 .OR. rimencrate_ice(i,j,k) > 0.0d0) THEN

              rime_q = rimeqcrate_ice(i,j,k)
              rime_n = rimencrate_ice(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,11) = rime_q
              END IF
#endif

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n
              ENDIF

            END IF

            IF (rimeqrrate_ice(i,j,k) > 0.0d0 .OR. rimenrrate_ice(i,j,k) > 0.0d0) THEN

              rime_q = rimeqrrate_ice(i,j,k)
              rime_n = rimenrrate_ice(i,j,k)
              rime_q = MIN(q_rain(i,j,k),rime_q)
              rime_n  = MIN(n_rain(i,j,k),rime_n)
              
              q_ice(i,j,k)  = q_ice(i,j,k)  + rime_q
              q_rain(i,j,k) = q_rain(i,j,k) - rime_q
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,62) = rime_q
              END IF
#endif
              
              ! Eismultiplikation nach Hallet und Mossop
              
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n
              ENDIF
              
            END IF

          ELSE

            !.. Depositionsrate ist kleiner als riming. Jetzt findet ice -> graupel bzw. ice -> hail statt:

            !.. Operator-Splitting: Zuerst ice_cloud_riming behandeln:

            IF (rimeqcrate_ice(i,j,k) > 0.0d0 .OR. rimencrate_ice(i,j,k) > 0.0d0) THEN

              n_i = n_ice(i,j,k)
              x_i = MIN(MAX(q_ice(i,j,k)/(n_i+eps),ice%x_min),ice%x_max)
              D_i = ice%a_geo * x_i**ice%b_geo

              rime_q = rimeqcrate_ice(i,j,k)
              rime_n = rimencrate_ice(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              !              rime_n = MIN(n_cloud(i,j,k),MIN(rime_n,const2*rime_q))
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n              
              ENDIF

              ! Umwandlung ice -> graupel

              IF (D_i > D_conv_ig) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_i**3/x_i - 1.0) )
                conv_q = MIN(q_ice(i,j,k),conv_q)
                x_i = MIN(MAX(q_ice(i,j,k)/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse incl. Riming
                conv_n = conv_q / MAX(x_i,x_conv) 
                conv_n = MIN(n_ice(i,j,k),conv_n)
              ELSE
                conv_q = 0.0
                conv_n = 0.0
              ENDIF

              q_ice(i,j,k)     = q_ice(i,j,k)     - conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q

              n_ice(i,j,k)     = n_ice(i,j,k)     - conv_n
              n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,12) = conv_q
                dqdt(i,j,k,11) = rime_q
              END IF
#endif
            END IF

            !.. Operator-Splitting: Dann ice_rain_riming behandeln:

            IF (rimeqirate_ice(i,j,k) > 0.0d0 .OR. rimeqrrate_ice(i,j,k) > 0.0d0 .OR.  rimenrrate_ice(i,j,k) > 0.0d0) THEN
              D_id = d_id_sp(i,j,k)
              D_rd = d_rd_sp_ice(i,j,k)

              x_r = MIN(MAX(q_rain(i,j,k)/(n_rain(i,j,k)+eps),rain%x_min),rain%x_max)

              rime_qi = rimeqirate_ice(i,j,k)
              rime_qr = rimeqrrate_ice(i,j,k)
              rime_n = rimenrrate_ice(i,j,k)
              rime_n  = MIN(MIN(n_rain(i,j,k),n_ice(i,j,k)),rime_n)
              rime_qr = MIN(q_rain(i,j,k),rime_qr)
              rime_qi = MIN(q_ice(i,j,k),rime_qi)

              n_ice(i,j,k)  = n_ice(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_ice(i,j,k)  = q_ice(i,j,k)  - rime_qi
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0  ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF


              IF (T_a >= T_3) THEN
                IF (D_id > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_ice(i,j,k) = n_ice(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_ice(i,j,k) = q_ice(i,j,k) + rime_qi     ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Eisteilchen schmelzen instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qi   ! UB_20081120 - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,29) = rime_qi
                  END IF
#endif
                END IF
              ELSE
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n  ! UB_20081120
                q_ice(i,j,k) = q_ice(i,j,k)  + mult_q  ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Eis + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,21) = mult_q
                    dqdt(i,j,k,20) = rime_qi + rime_qr
                  END IF
#endif
                ELSE
                  IF (D_id > D_rd) THEN
                    ! Eis + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,21) = mult_q
                      dqdt(i,j,k,20) = rime_qi + rime_qr
                    END IF
#endif
                  ELSE
                    ! Eis + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,21) = mult_q
                      dqdt(i,j,k,61) = rime_qi + rime_qr
                    END IF
#endif
                  END IF
                END IF
              END IF

            END IF

          END IF

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Vollenden des snow-rimings:
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          IF (deprate_snow(i,j,k) > 0.0d0 .AND. deprate_snow(i,j,k) >= rimeqcrate_snow(i,j,k)+rimeqrrate_snow(i,j,k)) THEN

            !.. Deposition ist groesser als gesamtes riming, deswegen bleibt Schnee Schnee:

            IF (rimeqcrate_snow(i,j,k) > 0.0d0 .OR. rimencrate_snow(i,j,k) > 0.0d0) THEN

              rime_q = rimeqcrate_snow(i,j,k)
              rime_n = rimencrate_snow(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,14) = mult_q
                dqdt(i,j,k,13) = rime_q
              END IF
#endif

            END IF

            IF (rimeqrrate_snow(i,j,k) > 0.0d0 .OR. rimenrrate_snow(i,j,k) > 0.0d0) THEN

              rime_q = rimeqrrate_snow(i,j,k)
              rime_n = rimenrrate_snow(i,j,k)
              rime_q = MIN(q_rain(i,j,k),rime_q)
              rime_n = MIN(n_rain(i,j,k),rime_n)

              q_snow(i,j,k) = q_snow(i,j,k) + rime_q
              q_rain(i,j,k) = q_rain(i,j,k) - rime_q
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,23) = mult_q
                dqdt(i,j,k,63) = rime_q
              END IF
#endif

            END IF

          ELSE

            !.. Operator-Splitting: Zuerst snow_cloud_riming behandeln:

            IF (rimeqcrate_snow(i,j,k) > 0.0d0 .OR. rimencrate_snow(i,j,k) > 0.0d0) THEN

              n_s = n_snow(i,j,k)
              x_s = MIN(MAX(q_snow(i,j,k)/(n_s+eps),snow%x_min),snow%x_max)
              D_s = snow%a_geo * x_s**snow%b_geo

              rime_q = rimeqcrate_snow(i,j,k)
              rime_n = rimencrate_snow(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

              ! Umwandlung snow -> graupel

              IF (D_s > D_conv_sg) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_s**3/x_s - 1.0) )
                conv_q = MIN(q_snow(i,j,k),conv_q)
                x_s = MIN(MAX(q_snow(i,j,k)/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse incl. Riming
                conv_n = conv_q / MAX(x_s,x_conv) 
                conv_n = MIN(n_snow(i,j,k),conv_n)
              ELSE
                conv_q = 0.0
                conv_n = 0.0
              ENDIF

              q_snow(i,j,k)    = q_snow(i,j,k)    - conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q

              n_snow(i,j,k)    = n_snow(i,j,k)    - conv_n
              n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,14) = mult_q
                dqdt(i,j,k,13) = rime_q
                dqdt(i,j,k,15) = conv_q
              END IF
#endif

            END IF

            !.. Operator-Splitting: Dann snow_rain_riming behandeln:

            IF (rimeqirate_snow(i,j,k) > 0.0d0 .OR. rimeqrrate_snow(i,j,k) > 0.0d0 .OR.  rimenrrate_snow(i,j,k) > 0.0d0) THEN

              D_sd = d_sd_sp(i,j,k)
              D_rd = d_rd_sp_snow(i,j,k)

              x_r = MIN(MAX(q_rain(i,j,k)/(n_rain(i,j,k)+eps),rain%x_min),rain%x_max)

              rime_qs = rimeqirate_snow(i,j,k)
              rime_qr = rimeqrrate_snow(i,j,k)
              rime_n  = rimenrrate_snow(i,j,k)
              rime_qr = MIN(q_rain(i,j,k),rime_qr)
              rime_qs = MIN(q_snow(i,j,k),rime_qs)
              rime_n  = MIN(n_rain(i,j,k),rime_n)
              rime_n  = MIN(n_snow(i,j,k),rime_n)

              n_snow(i,j,k) = n_snow(i,j,k) - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_snow(i,j,k) = q_snow(i,j,k) - rime_qs
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0 ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_sd > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_snow(i,j,k) = n_snow(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_snow(i,j,k) = q_snow(i,j,k) + rime_qs ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Schneeflocken werden instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qs ! UB_20081120 - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,30) = rime_qs
                  END IF
#endif
                END IF
              ELSE
                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n ! UB_20081120
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,23) = mult_q
                    dqdt(i,j,k,22) = rime_qr + rime_qs
                  END IF
#endif
                ELSE
                  IF (D_sd > D_rd) THEN
                    ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,23) = mult_q
                      dqdt(i,j,k,22) = rime_qr + rime_qs
                    END IF
#endif
                  ELSE
                    ! Schnee + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,23) = mult_q
                      dqdt(i,j,k,57) = rime_qr + rime_qs
                  END IF
#endif
                  END IF
                END IF
              END IF

            END IF

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE complete_ice_snow_riming

! Vectorized version (hopefully ...)
  SUBROUTINE complete_ice_snow_rim_vec ()
    !*******************************************************************************
    !                                                                              *
    !       Im Falle von use_ice_graupel_conv_uli = .true. :                       *
    !                                                                              *
    !       Summierung der vorher lediglich gespeicherten riming-Raten             *
    !       auf die entsprechenden Felder (ice, snow, cloud, rain)                 *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_ice, n_ice, q_snow, n_snow, q_cloud, n_cloud, q_rain, n_rain, &
         &                        q_graupel, n_graupel, q_hail, n_hail, &
         &                        dt, dqdt, speichere_dqdt, &
         &                        deprate_ice, deprate_snow, &
         &                        rimeqcrate_ice, rimencrate_ice, &
         &                        rimeqirate_ice, rimeqrrate_ice, rimenrrate_ice, &
         &                        rimeqcrate_snow, rimencrate_snow, &
         &                        rimeqirate_snow, rimeqrrate_snow, rimenrrate_snow, &
         &                        d_id_sp, d_rd_sp_ice, d_sd_sp, d_rd_sp_snow

    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen
    INTEGER                     :: i,j,k

    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,D_id
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,D_sd
    DOUBLE PRECISION            :: x_coll_c,x_r
    DOUBLE PRECISION            :: D_rd
    DOUBLE PRECISION            :: rime_n,rime_q, rime_qr, rime_qs, rime_qi
    DOUBLE PRECISION            :: conv_n,conv_q
    DOUBLE PRECISION            :: mult_n,mult_q,mult_1,mult_2
    DOUBLE PRECISION            :: const2,const3,const4,const5
    DOUBLE PRECISION, SAVE      :: D_av_fakt_i,D_av_fakt_s

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      IF (ice_multiplication) THEN
        WRITE (6, *) "CLOUD complete_ice_snow_rim_vec mit ice_multiplication"
      ELSE
        WRITE (6, *) "CLOUD complete_ice_snow_rim_vec ohne ice_multiplication"
      ENDIF
    END IF

    IF (firstcall .NE. 1) THEN
      D_av_fakt_i = D_average_factor(ice)
      D_av_fakt_s = D_average_factor(snow)
      firstcall = 1
    END IF

    x_coll_c = (D_coll_c/cloud%a_geo)**3          !..Mindestmasse fuer collection, begrenzt rime_n bei n_cloud

    const2 = 1.0/x_coll_c
    const3 = 1.0/(T_mult_opt - T_mult_min)
    const4 = 1.0/(T_mult_opt - T_mult_max)
    const5 = alpha_spacefilling * rho_w/rho_ice

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !.. Vollenden des ice-rimings:
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !.. 1) Depositional growth is stronger than riming growth, therefore ice stays ice:
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a = T_0(i,j,k)

          IF (deprate_ice(i,j,k) > 0.0d0 .AND. deprate_ice(i,j,k) >= rimeqcrate_ice(i,j,k)+rimeqrrate_ice(i,j,k)) THEN

            !.. Operator-Splitting: Zuerst ice_cloud_riming behandeln:

            IF (rimeqcrate_ice(i,j,k) > 0.0d0 .OR. rimencrate_ice(i,j,k) > 0.0d0) THEN

              rime_q = rimeqcrate_ice(i,j,k)
              rime_n = rimencrate_ice(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,11) = rime_q
              END IF
#endif

              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n
              ENDIF

            END IF

            !.. Operator-Splitting: Dann ice_rain_riming behandeln:

            IF (rimeqrrate_ice(i,j,k) > 0.0d0 .OR. rimenrrate_ice(i,j,k) > 0.0d0) THEN

              rime_q = rimeqrrate_ice(i,j,k)
              rime_n = rimenrrate_ice(i,j,k)
              rime_q = MIN(q_rain(i,j,k),rime_q)
              rime_n  = MIN(n_rain(i,j,k),rime_n)
              
              q_ice(i,j,k)  = q_ice(i,j,k)  + rime_q
              q_rain(i,j,k) = q_rain(i,j,k) - rime_q
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,62) = rime_q
              END IF
#endif
              
              ! Eismultiplikation nach Hallet und Mossop
              
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n
              ENDIF
              
            END IF

          END IF
        END DO
      END DO
    END DO

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !.. Vollenden des ice-rimings:
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !.. 2) Depositional growth is smaller than riming growth, therefore ice is 
    !      allowed to convert to graupel and / or hail:
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a = T_0(i,j,k)

          IF (deprate_ice(i,j,k) <= 0.0d0 .OR. deprate_ice(i,j,k) < rimeqcrate_ice(i,j,k)+rimeqrrate_ice(i,j,k)) THEN

            !.. Operator-Splitting: Zuerst ice_cloud_riming behandeln:

            IF (rimeqcrate_ice(i,j,k) > 0.0d0 .OR. rimencrate_ice(i,j,k) > 0.0d0) THEN

              n_i = n_ice(i,j,k)
              x_i = MIN(MAX(q_ice(i,j,k)/(n_i+eps),ice%x_min),ice%x_max)
              D_i = ice%a_geo * x_i**ice%b_geo

              rime_q = rimeqcrate_ice(i,j,k)
              rime_n = rimencrate_ice(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
!              rime_n = MIN(n_cloud(i,j,k),MIN(rime_n,const2*rime_q))
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_ice(i,j,k)   = q_ice(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4 
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q

                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n              
              ENDIF

              ! Umwandlung ice -> graupel

              IF (D_i > D_conv_ig) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_i**3/x_i - 1.0) )
                conv_q = MIN(q_ice(i,j,k),conv_q)
                x_i = MIN(MAX(q_ice(i,j,k)/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse incl. Riming
                conv_n = conv_q / MAX(x_i,x_conv) 
                conv_n = MIN(n_ice(i,j,k),conv_n)
              ELSE
                conv_q = 0.0
                conv_n = 0.0
              ENDIF

              q_ice(i,j,k)     = q_ice(i,j,k)     - conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q

              n_ice(i,j,k)     = n_ice(i,j,k)     - conv_n
              n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,12) = conv_q
                dqdt(i,j,k,11) = rime_q
              END IF
#endif
            END IF

            !.. Operator-Splitting: Dann ice_rain_riming behandeln:

            IF (rimeqirate_ice(i,j,k) > 0.0d0 .OR. rimeqrrate_ice(i,j,k) > 0.0d0 .OR.  rimenrrate_ice(i,j,k) > 0.0d0) THEN
              D_id = d_id_sp(i,j,k)
              D_rd = d_rd_sp_ice(i,j,k)

              x_r = MIN(MAX(q_rain(i,j,k)/(n_rain(i,j,k)+eps),rain%x_min),rain%x_max)

              rime_qi = rimeqirate_ice(i,j,k)
              rime_qr = rimeqrrate_ice(i,j,k)
              rime_n = rimenrrate_ice(i,j,k)
              rime_n  = MIN(MIN(n_rain(i,j,k),n_ice(i,j,k)),rime_n)
              rime_qr = MIN(q_rain(i,j,k),rime_qr)
              rime_qi = MIN(q_ice(i,j,k),rime_qi)

              n_ice(i,j,k)  = n_ice(i,j,k)  - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_ice(i,j,k)  = q_ice(i,j,k)  - rime_qi
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0  ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)
              ENDIF


              IF (T_a >= T_3) THEN
                IF (D_id > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_ice(i,j,k) = n_ice(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_ice(i,j,k) = q_ice(i,j,k) + rime_qi     ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Eisteilchen schmelzen instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qi   ! UB_20081120 - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,29) = rime_qi
                  END IF
#endif
                END IF
              ELSE
                n_ice(i,j,k) = n_ice(i,j,k)  + mult_n  ! UB_20081120
                q_ice(i,j,k) = q_ice(i,j,k)  + mult_q  ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Eis + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,21) = mult_q
                    dqdt(i,j,k,20) = rime_qi + rime_qr
                  END IF
#endif
                ELSE
                  IF (D_id > D_rd) THEN
                    ! Eis + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,21) = mult_q
                      dqdt(i,j,k,20) = rime_qi + rime_qr
                    END IF
#endif
                  ELSE
                    ! Eis + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qi + rime_qr - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,21) = mult_q
                      dqdt(i,j,k,61) = rime_qi + rime_qr
                    END IF
#endif
                  END IF
                END IF
              END IF

            END IF

          END IF
        END DO
      END DO
    END DO

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Vollenden des snow-rimings:
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !.. 1) Depositional growth is stronger than riming growth, therefore snow stays snow:
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a = T_0(i,j,k)

          IF (deprate_snow(i,j,k) > 0.0d0 .AND. deprate_snow(i,j,k) >= rimeqcrate_snow(i,j,k)+rimeqrrate_snow(i,j,k)) THEN

            !.. Operator-Splitting: Zuerst snow_cloud_riming behandeln:

            IF (rimeqcrate_snow(i,j,k) > 0.0d0 .OR. rimencrate_snow(i,j,k) > 0.0d0) THEN

              rime_q = rimeqcrate_snow(i,j,k)
              rime_n = rimencrate_snow(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,14) = mult_q
                dqdt(i,j,k,13) = rime_q
              END IF
#endif

            END IF

            !.. Operator-Splitting: Dann snow_rain_riming behandeln:

            IF (rimeqrrate_snow(i,j,k) > 0.0d0 .OR. rimenrrate_snow(i,j,k) > 0.0d0) THEN

              rime_q = rimeqrrate_snow(i,j,k)
              rime_n = rimenrrate_snow(i,j,k)
              rime_q = MIN(q_rain(i,j,k),rime_q)
              rime_n = MIN(n_rain(i,j,k),rime_n)

              q_snow(i,j,k) = q_snow(i,j,k) + rime_q
              q_rain(i,j,k) = q_rain(i,j,k) - rime_q
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min)*const3
                mult_2 = (T_a - T_mult_max)*const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,23) = mult_q
                dqdt(i,j,k,63) = rime_q
              END IF
#endif

            END IF

          END IF
        END DO
      END DO
    END DO

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !.. Vollenden des snow-rimings:
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    !.. 2) Depositional growth is smaller than riming growth, therefore snow is 
    !      allowed to convert to graupel and / or hail:
    !
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a = T_0(i,j,k)

          IF (deprate_snow(i,j,k) <= 0.0d0 .AND. deprate_snow(i,j,k) < rimeqcrate_snow(i,j,k)+rimeqrrate_snow(i,j,k)) THEN

            !.. Operator-Splitting: Zuerst snow_cloud_riming behandeln:

            IF (rimeqcrate_snow(i,j,k) > 0.0d0 .OR. rimencrate_snow(i,j,k) > 0.0d0) THEN

              n_s = n_snow(i,j,k)
              x_s = MIN(MAX(q_snow(i,j,k)/(n_s+eps),snow%x_min),snow%x_max)
              D_s = snow%a_geo * x_s**snow%b_geo

              rime_q = rimeqcrate_snow(i,j,k)
              rime_n = rimencrate_snow(i,j,k)
              rime_q = MIN(q_cloud(i,j,k),rime_q)
              rime_n = MIN(n_cloud(i,j,k),rime_n)

              q_snow(i,j,k)  = q_snow(i,j,k)  + rime_q
              q_cloud(i,j,k) = q_cloud(i,j,k) - rime_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - rime_n

              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_q
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_q,mult_q)

                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q
                q_snow(i,j,k) = q_snow(i,j,k) - mult_q
              ENDIF

              ! Umwandlung snow -> graupel

              IF (D_s > D_conv_sg) THEN
                conv_q = (rime_q - mult_q) / ( const5 * (pi/6.0*rho_ice*d_s**3/x_s - 1.0) )
                conv_q = MIN(q_snow(i,j,k),conv_q)
                x_s = MIN(MAX(q_snow(i,j,k)/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse incl. Riming
                conv_n = conv_q / MAX(x_s,x_conv) 
                conv_n = MIN(n_snow(i,j,k),conv_n)
              ELSE
                conv_q = 0.0
                conv_n = 0.0
              ENDIF

              q_snow(i,j,k)    = q_snow(i,j,k)    - conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) + conv_q

              n_snow(i,j,k)    = n_snow(i,j,k)    - conv_n
              n_graupel(i,j,k) = n_graupel(i,j,k) + conv_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,14) = mult_q
                dqdt(i,j,k,13) = rime_q
                dqdt(i,j,k,15) = conv_q
              END IF
#endif

            END IF

            !.. Operator-Splitting: Dann snow_rain_riming behandeln:

            IF (rimeqirate_snow(i,j,k) > 0.0d0 .OR. rimeqrrate_snow(i,j,k) > 0.0d0 .OR.  rimenrrate_snow(i,j,k) > 0.0d0) THEN

              D_sd = d_sd_sp(i,j,k)
              D_rd = d_rd_sp_snow(i,j,k)

              x_r = MIN(MAX(q_rain(i,j,k)/(n_rain(i,j,k)+eps),rain%x_min),rain%x_max)

              rime_qs = rimeqirate_snow(i,j,k)
              rime_qr = rimeqrrate_snow(i,j,k)
              rime_n  = rimenrrate_snow(i,j,k)
              rime_qr = MIN(q_rain(i,j,k),rime_qr)
              rime_qs = MIN(q_snow(i,j,k),rime_qs)
              rime_n  = MIN(n_rain(i,j,k),rime_n)
              rime_n  = MIN(n_snow(i,j,k),rime_n)

              n_snow(i,j,k) = n_snow(i,j,k) - rime_n
              n_rain(i,j,k) = n_rain(i,j,k) - rime_n
              q_snow(i,j,k) = q_snow(i,j,k) - rime_qs
              q_rain(i,j,k) = q_rain(i,j,k) - rime_qr


              ! Eismultiplikation nach Hallet und Mossop

              mult_q = 0.0
              mult_n = 0.0 ! UB_20081120
              IF (T_a < T_3 .AND. ice_multiplication) THEN
                mult_1 = (T_a - T_mult_min) * const3
                mult_2 = (T_a - T_mult_max) * const4
                mult_1 = MAX(0.d0,MIN(mult_1,1.d0))
                mult_2 = MAX(0.d0,MIN(mult_2,1.d0))
                mult_n = C_mult * mult_1 * mult_2 * rime_qr
                mult_q = mult_n * ice%x_min
                mult_q = MIN(rime_qr,mult_q)

              ENDIF

              IF (T_a >= T_3) THEN
                IF (D_sd > D_rd) THEN
                  ! Regen wird wieder abgeshedded
                  n_snow(i,j,k) = n_snow(i,j,k) + rime_n
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_qr / x_r
                  q_snow(i,j,k) = q_snow(i,j,k) + rime_qs ! UB_20081120 - mult_q
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr
                ELSE
                  ! Schneeflocken werden instantan zu Regen.
                  n_rain(i,j,k) = n_rain(i,j,k) + rime_n
                  q_rain(i,j,k) = q_rain(i,j,k) + rime_qr + rime_qs ! UB_20081120 - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,30) = rime_qs
                  END IF
#endif
                END IF
              ELSE
                n_ice(i,j,k)  = n_ice(i,j,k)  + mult_n ! UB_20081120
                q_ice(i,j,k)  = q_ice(i,j,k)  + mult_q ! UB_20081120
                IF (ice_typ < 3) THEN
                  ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                  n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                  q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                  IF (speichere_dqdt) THEN
                    dqdt(i,j,k,23) = mult_q
                    dqdt(i,j,k,22) = rime_qr + rime_qs
                  END IF
#endif
                ELSE
                  IF (D_sd > D_rd) THEN
                    ! Schnee + angefrorenes Regenwasser ergibt Graupel:
                    n_graupel(i,j,k) = n_graupel(i,j,k) + rime_n
                    q_graupel(i,j,k) = q_graupel(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,23) = mult_q
                      dqdt(i,j,k,22) = rime_qr + rime_qs
                    END IF
#endif
                  ELSE
                    ! Schnee + angefrorenes Regenwasser ergibt Hagel:
                    n_hail(i,j,k) = n_hail(i,j,k) + rime_n
                    q_hail(i,j,k) = q_hail(i,j,k) + rime_qr + rime_qs - mult_q
#ifdef SAVE_CONVERSIONRATES
                    IF (speichere_dqdt) THEN
                      dqdt(i,j,k,23) = mult_q
                      dqdt(i,j,k,57) = rime_qr + rime_qs
                  END IF
#endif
                  END IF
                END IF
              END IF

            END IF

          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE complete_ice_snow_rim_vec

  SUBROUTINE graupel_melting()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Schmelzens der Graupelpartikel                          *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_graupel, n_graupel, q_rain, n_rain, dt, dqdt,   &
         &                        speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g,T_a,N_re,D_T,e_a
    DOUBLE PRECISION            :: melt,melt_v,melt_h,melt_n,melt_q
    DOUBLE PRECISION            :: fh_n,fv_n
    DOUBLE PRECISION            :: fh_q,fv_q
    DOUBLE PRECISION, SAVE      :: a_melt_n,b_melt_n
    DOUBLE PRECISION, SAVE      :: a_melt_q,b_melt_q

    DOUBLE PRECISION, PARAMETER :: eps   = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD graupel_melting",4)

    IF (firstcall.NE.1) THEN
      a_melt_n = vent_coeff_a(graupel,0)
      b_melt_n = vent_coeff_b(graupel,0)
      a_melt_q = vent_coeff_a(graupel,1)
      b_melt_q = vent_coeff_b(graupel,1)
      firstcall = 1
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD graupel_melting " 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:   a_geo  = ",graupel%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_geo  = ",graupel%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Graupel:       a_vel  = ",graupel%a_vel 
        WRITE (6,'(A,D10.3)') "                                       b_vel  = ",graupel%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Ventilation von Graupel: a_ven  = ",graupel%a_ven 
        WRITE (6,'(A,D10.3)') "                                       b_ven  = ",graupel%b_ven 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Schmelzen:               a_melt_n = ",a_melt_n
        WRITE (6,'(A,D10.3)') "                                       b_melt_n = ",b_melt_n
        WRITE (6,'(A,D10.3)') "                                       a_melt_q = ",a_melt_q
        WRITE (6,'(A,D10.3)') "                                       b_melt_q = ",b_melt_q
      ENDIF
    ELSEIF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD graupel_melting " 
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (T_a > T_3 .AND. q_g > 0.0) THEN
            e_a = e_ws(T_a)                                     !..Saettigungsdampfdruck
            n_g = n_graupel(i,j,k)                              !..Anzahldichte in SI

            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max)  !..mittlere Masse in SI     

            D_g = graupel%a_geo * x_g**graupel%b_geo                   !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_g * D_g / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
  !          fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

! UB_20081125:            D_T  = K_T / (cp * rho_0(i,j,k))!WRF!+rho_g(i,j,k)))
! UB_20081125:            fh_q = D_T / D_v * fv_q  ! <-- gives wrong values!
! UB_20081125: After formulas of Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
!              for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q
            !           fh_n = D_T / D_v * fv_n

            melt   = 2.0*pi / L_ew * D_g * n_g * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_g

! ub>> setzte melt_n so, dass x_h beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q_g) / x_g + n_g, 0.0d0), n_g)

            melt_q = MIN(q_g,melt_q)
            melt_n = MIN(n_g,melt_n)

            melt_q = MAX(0.d0,melt_q)
            melt_n = MAX(0.d0,melt_n)

            rate_q(i,j,k) = melt_q

            q_graupel(i,j,k) = q_graupel(i,j,k) - melt_q
            q_rain(i,j,k)    = q_rain(i,j,k)    + melt_q

            n_graupel(i,j,k) = n_graupel(i,j,k) - melt_n
            n_rain(i,j,k)    = n_rain(i,j,k)    + melt_n
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,31) = melt_q
            END IF
#endif
            ! ub<<

            !WRF!T(i,j,k) = T(i, j, k) - L_ew / cp * melt_q / (rho_0(i,j,k)+rho_g(i,j,k)) 
            !WRF!p(i,j,k) = p(i, j, k) - L_ew / cv * melt_q * R_l
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(rate_q)

  END SUBROUTINE graupel_melting

  SUBROUTINE hail_melting()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Schmelzens der Hagelpartikel                            *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_hail, n_hail, q_rain, n_rain, dt, dqdt,         &
         &                        speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h,T_a,N_re,D_T,e_a
    DOUBLE PRECISION            :: melt,melt_v,melt_h,melt_n,melt_q
    DOUBLE PRECISION            :: fh_n,fv_n
    DOUBLE PRECISION            :: fh_q,fv_q
    DOUBLE PRECISION, SAVE      :: a_melt_n,b_melt_n
    DOUBLE PRECISION, SAVE      :: a_melt_q,b_melt_q

    DOUBLE PRECISION, PARAMETER :: eps   = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD hail_melting",4)

    IF (firstcall.NE.1) THEN
      a_melt_n = vent_coeff_a(hail,0)
      b_melt_n = vent_coeff_b(hail,0)
      a_melt_q = vent_coeff_a(hail,1)
      b_melt_q = vent_coeff_b(hail,1)
      firstcall = 1
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD hail_melting " 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hail:      a_geo  = ",hail%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_geo  = ",hail%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Hail:          a_vel  = ",hail%a_vel 
        WRITE (6,'(A,D10.3)') "                                       b_vel  = ",hail%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Ventilation von Hail:    a_ven  = ",hail%a_ven 
        WRITE (6,'(A,D10.3)') "                                       b_ven  = ",hail%b_ven 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Schmelzen:               a_melt_n = ",a_melt_n
        WRITE (6,'(A,D10.3)') "                                       b_melt_n = ",b_melt_n
        WRITE (6,'(A,D10.3)') "                                       a_melt_q = ",a_melt_q
        WRITE (6,'(A,D10.3)') "                                       b_melt_q = ",b_melt_q
      ENDIF
    ELSEIF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD hail_melting " 
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          q_h = q_hail(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (T_a > T_3 .AND. q_h > 0.0) THEN
            e_a = e_ws(T_a)                                     !..Saettigungsdampfdruck
            n_h = n_hail(i,j,k)                              !..Anzahldichte in SI

            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)  !..mittlere Masse in SI     

            D_h = hail%a_geo * x_h**hail%b_geo                   !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_h * D_h / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
!            fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

! UB_20081125:            D_T  = K_T / (cp * rho_0(i,j,k))!WRF!+rho_g(i,j,k)))
! UB_20081125:            fh_q = D_T / D_v * fv_q  ! <-- gives wrong values!
! UB_20081125: After formulas of Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
!              for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q
!            fh_n = D_T / D_v * fv_n

            melt   = 2.0*pi / L_ew * D_h * n_h * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_h

! ub>> setzte melt_n so, dass x_h beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q_h) / x_h + n_h, 0.0d0), n_h)

            melt_q = MIN(q_h,melt_q)
            melt_n = MIN(n_h,melt_n)

            melt_q = MAX(0.d0,melt_q)
            melt_n = MAX(0.d0,melt_n)

            rate_q(i,j,k) = melt_q

            q_hail(i,j,k) = q_hail(i,j,k) - melt_q
            q_rain(i,j,k) = q_rain(i,j,k)    + melt_q

            n_hail(i,j,k) = n_hail(i,j,k) - melt_n
            n_rain(i,j,k) = n_rain(i,j,k)    + melt_n
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,52) = melt_q
            END IF
#endif
            ! ub<<

            !WRF!T(i,j,k) = T(i, j, k) - L_ew / cp * melt_q / (rho_0(i,j,k)+rho_g(i,j,k)) 
            !WRF!p(i,j,k) = p(i, j, k) - L_ew / cv * melt_q * R_l
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(rate_q)

  END SUBROUTINE hail_melting

  SUBROUTINE cloud_freeze ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Wolkentropfen                      *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_cloud, n_cloud, q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER          :: i, j, k
    DOUBLE PRECISION :: fr_q, fr_n,T_a,q_c,x_c,n_c,j_het,j_hom,j_tot,T_c
    !DOUBLE PRECISION, PARAMETER :: a_HET = 6.5d-1 ! 1/K,      Messung nach Barklie and Gokhale
    !DOUBLE PRECISION, PARAMETER :: b_HET = 2.0d+2 ! 1/(m3 s), Messung nach Barklie and Gokhale
    DOUBLE PRECISION, PARAMETER :: a_HET = 6.6d-1 ! 1/K,      Messung nach  (PK S.350)
    DOUBLE PRECISION, PARAMETER :: b_HET = 1.0d+2 ! 1/(m3 s), Messung nach Bigg (PK S.350) 
    DOUBLE PRECISION, PARAMETER :: eps   = 1.d-20 
    DOUBLE PRECISION            :: facg

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD cloud_freeze" 
    END IF

    facg = moment_gamma(cloud,2)    ! <hn

    !..Test auf Schmelzen oder Gefrieren von Wolkenteilchen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a  = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)
          IF (T_a < T_3) THEN
            q_c = q_cloud(i,j,k)
            n_c = n_cloud(i,j,k)
            T_c = T_a - T_3
            IF (q_c > 0.0) THEN
              IF (T_c < -50.0) THEN            
                fr_q = q_c                                               !..Komplettes hom. Gefrieren
                fr_n = n_c                                               !..unterhalb -50 C
              ELSE                                      
                x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max)    !..mittlere Masse
                
                !..Hom. Gefrieren nach Jeffrey und Austin (1997), siehe auch Cotton und Field (2001)
                IF (T_c > -30.0) THEN            
                  j_hom = 1.0d6 * EXP(-7.63-2.996*(T_c+30.0))           !..J in 1/(m3 s) 
                ELSE
                  j_hom = 1.0d6 * EXP(-243.4-14.75*T_c-0.307*T_c**2-0.00287*T_c**3-0.0000102*T_c**4)
                ENDIF

                !..Het. Gefrieren: stochastisches Modell nach Bigg, Daten nach Barklie and Gokhale
                !j_het = b_HET * ( EXP( - a_HET * T_c) - 1.0 )            !..J in 1/(m3 s)
                j_het = 0.0 ! neglected for cloud droplets
                
                !..Umwandlungsraten fuer Anzahl- und Massendichte
                j_tot = (j_hom + j_het) / rho_w * dt                     !..J*dt in 1/kg
                fr_n  = j_tot * q_c
                fr_q  = j_tot * q_c * x_c * facg
              END IF
              fr_q  = MIN(fr_q,q_c)
              fr_n  = MIN(fr_n,n_c)

              !..Berechnung der H2O-Komponenten
              q_cloud(i,j,k) = q_cloud(i,j,k) - fr_q
              n_cloud(i,j,k) = n_cloud(i,j,k) - fr_n

              fr_n  = MAX(fr_n,fr_q/cloud%x_max)
              
              IF (nuc_c_typ .EQ. 0) THEN
                IF (isdebug) WRITE (*,*) '  ... force upper bound in cloud_freeze'
                fr_n = MAX(MIN(fr_n,qnc_const-n_ice(i,j,k)),0.d0)
              ENDIF
              
              q_ice(i,j,k)   = q_ice(i,j,k)   + fr_q
              n_ice(i,j,k)   = n_ice(i,j,k)   + fr_n

              ENDIF
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,3) = fr_q
            END IF
#endif
            ! ub<<

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE cloud_freeze

  SUBROUTINE rain_freeze_gamlook ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Regentropfen                       *
    ! 
    ! Auftretende unvollst. gamma-Funktionen werden mittels lookup tables
    ! effizienter berechnet als in rain_freeze()
    !
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_rain, n_rain, q_graupel, n_graupel,             &
         &                        q_hail, n_hail, q_ice, n_ice,                     &
         &                        dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0
    USE gamma_functions_mp_seifert, ONLY: incgfct_lower_lookupcreate, incgfct_lower_lookup, &
         &                        gamlookuptable, nlookup, nlookuphr_dummy

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: nuc_typ
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het, &
         &  fr_q_i,fr_n_i,fr_q_g,fr_n_g,fr_q_h,fr_n_h,n_0,lam,xmax_ice,xmax_gr,fr_q_tmp,fr_n_tmp
    DOUBLE PRECISION, PARAMETER :: a_HET = 6.5d-1 ! Messung nach Barklie and Gokhale (PK S.350)
    DOUBLE PRECISION, PARAMETER :: b_HET = 2.0d+2 ! Messung nach Barklie and Gokhale (PK S.350)
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    TYPE(gamlookuptable), SAVE :: ltable1, ltable2, ltable3

    DOUBLE PRECISION, SAVE      :: coeff_z
    DOUBLE PRECISION, SAVE :: nm1, nm2, nm3, g1, g2, g3


    IF (firstcall.NE.1) THEN
      firstcall = 1
      !..Koeff. fuer Reflektivitaet Z (2. Moment)
      coeff_z = moment_gamma(rain,2)
      nm1 = (rain%nu+1.0)/rain%mu
      nm2 = (rain%nu+2.0)/rain%mu
      nm3 = (rain%nu+3.0)/rain%mu
      CALL incgfct_lower_lookupcreate(nm1, ltable1, nlookup, nlookuphr_dummy)
      CALL incgfct_lower_lookupcreate(nm2, ltable2, nlookup, nlookuphr_dummy)
      CALL incgfct_lower_lookupcreate(nm3, ltable3, nlookup, nlookuphr_dummy)
      ! ordinary gamma function of nm1 is the last value in the lookup table 1:
      g1 = ltable1%igf(ltable1%n)
      ! ordinary gamma function of nm2 is the last value in the lookup table 2:
      g2 = ltable2%igf(ltable2%n)
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD rain_freeze_gamlook" 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Z von Regentropfen: coeff_z= ",coeff_z
      ENDIF
    ELSE IF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD rain_freeze_gamlook" 
    ENDIF

    xmax_ice = ( (D_rainfrz_ig/rain%a_geo)**(1.0d0/rain%b_geo) ) ** rain%mu
    xmax_gr  = ( (D_rainfrz_gh/rain%a_geo)**(1.0d0/rain%b_geo) ) ** rain%mu

    nuc_typ = nuc_i_typ

    !..Test auf Schmelzen oder Gefrieren von Regentropfen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k)
          q_r = q_rain(i,j,k)
          n_r = n_rain(i,j,k)

          IF (T_a < T_freeze) THEN
            IF (q_r <= q_krit_fr) THEN
              IF (T_a < T_f) THEN
                fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                fr_n = n_r
                fr_n_i= n_r
                fr_q_i= q_r
                fr_n_g= 0.0
                fr_q_g= 0.0
                fr_n_h= 0.0
                fr_q_h= 0.0
                ! UB_20080220>
                fr_n_tmp = 1.0
                fr_q_tmp = 1.0
                ! <UB_20080220
              ELSE
                fr_q = 0.0
                fr_n = 0.0
                fr_n_i= 0.0
                fr_q_i= 0.0
                fr_n_g= 0.0
                fr_q_g= 0.0
                fr_n_h= 0.0
                fr_q_h= 0.0
                ! UB_20080220>
                fr_n_tmp = 0.0
                fr_q_tmp = 0.0
                ! <UB_20080220
              END IF
            ELSE
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              n_r = q_r / x_r
              IF (T_a < T_f) THEN            !..Nur Eis
!!! Diesen Zweig koennte man auch weglassen. ist zudem zwar quanitativ richtig, aber nicht konsistent zum
!!! Grenzfall fuer komplettes Gefrieren der Rechnung im T_a >= T_f - Zweig weiter unten
                fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                fr_n = n_r

                ! ub>> Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                !      oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                !      bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                !      und von xmax_gr bis unendlich (--> Hagel).

                lam = ( g1 / g2 * x_r)**(-rain%mu)
                n_0 = rain%mu * n_r * lam**(nm1) / g1
                fr_n_i = n_0/(rain%mu*lam**(nm1))* &
                     incgfct_lower_lookup(lam*xmax_ice, ltable1)
                fr_q_i = n_0/(rain%mu*lam**(nm2))* &
                     incgfct_lower_lookup(lam*xmax_ice, ltable2)
                fr_n_g = n_0/(rain%mu*lam**(nm1))* &
                     incgfct_lower_lookup(lam*xmax_gr,  ltable1)
                fr_q_g = n_0/(rain%mu*lam**(nm2))* &
                     incgfct_lower_lookup(lam*xmax_gr,  ltable2)


                fr_n_h = fr_n - fr_n_g
                fr_q_h = fr_q - fr_q_g
                fr_n_g = fr_n_g - fr_n_i
                fr_q_g = fr_q_g - fr_q_i
                fr_n_tmp = n_r/MAX(fr_n,n_r)
                fr_q_tmp = q_r/MAX(fr_q,q_r)

              ELSE                           !..Heterogenes Gefrieren
                j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.d0) / rho_w * dt

                ! ub>> Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                !      oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                !      bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                !      und von xmax_gr bis unendlich (--> Hagel).

                IF (j_het >= 1d-20) THEN
                  fr_n  = j_het * q_r
                  fr_q  = j_het * q_r * x_r * coeff_z

                  lam = ( g1 / g2 * x_r)**(-rain%mu)
                  n_0 = rain%mu * n_r * lam**(nm1) / g1
                  fr_n_i = j_het * n_0/(rain%mu*lam**(nm2))* &
                       incgfct_lower_lookup(lam*xmax_ice, ltable2)
                  fr_q_i = j_het * n_0/(rain%mu*lam**(nm3))* &
                       incgfct_lower_lookup(lam*xmax_ice, ltable3)
                  fr_n_g = j_het * n_0/(rain%mu*lam**(nm2))* &
                       incgfct_lower_lookup(lam*xmax_gr,  ltable2)
                  fr_q_g = j_het * n_0/(rain%mu*lam**(nm3))* &
                       incgfct_lower_lookup(lam*xmax_gr,  ltable3)

                  fr_n_h = fr_n - fr_n_g
                  fr_q_h = fr_q - fr_q_g
                  fr_n_g = fr_n_g - fr_n_i
                  fr_q_g = fr_q_g - fr_q_i
                  fr_n_tmp = n_r/MAX(fr_n,n_r)
                  fr_q_tmp = q_r/MAX(fr_q,q_r)
                ELSE
                  fr_n= 0.0
                  fr_q= 0.0
                  fr_n_i= 0.0
                  fr_q_i= 0.0
                  fr_n_g= 0.0
                  fr_q_g= 0.0
                  ! UB_20080212>
                  fr_n_h= 0.0
                  fr_q_h= 0.0
                  ! <UB_20080212
                  fr_n_tmp = 0.0
                  fr_q_tmp = 0.0
                END IF

              END IF

              fr_n = fr_n * fr_n_tmp
              fr_q = fr_q * fr_q_tmp

              fr_n_i = fr_n_i * fr_n_tmp
              fr_n_g = fr_n_g * fr_n_tmp
              fr_n_h = fr_n_h * fr_n_tmp
              fr_q_i = fr_q_i * fr_q_tmp
              fr_q_g = fr_q_g * fr_q_tmp
              fr_q_h = fr_q_h * fr_q_tmp

            END IF

            !..Berechnung der H2O-Komponenten

            q_rain(i,j,k) = q_rain(i,j,k) - fr_q
            n_rain(i,j,k) = n_r - fr_n

            !          IF (q_rain(i,j,k) < 0.0) THEN
            !            !WRITE (*,*) 'SEIFERT RAIN_FREEZE_GAMLOOK: Qrain < 0.0, ', i,j,k, T_a, q_r, j_het, fr_q, fr_q_tmp
            !            q_rain(i,j,k) = 0.0d0
            !          END IF
            !          IF (n_rain(i,j,k) < 0.0) THEN
            !            !WRITE (*,*) 'SEIFERT RAIN_FREEZE_GAMLOOK: Nrain < 0.0, ', i,j,k, T_a, n_r, j_het, fr_n, fr_n_tmp
            !            n_rain(i,j,k) = 0.0d0
            !          END IF

            IF (ice_typ < 2) THEN 
              ! ohne Hagelklasse,  gefrierender Regen wird Eis oder Graupel
              q_ice(i,j,k) = q_ice(i,j,k)  + fr_q_i
              n_ice(i,j,k) = n_ice(i,j,k)  + fr_n_i
              q_graupel(i,j,k) = q_graupel(i,j,k)  + fr_q_h + fr_q_g
              n_graupel(i,j,k) = n_graupel(i,j,k)  + fr_n_h + fr_n_g
            ELSE
              ! mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel
              q_ice(i,j,k) = q_ice(i,j,k)  + fr_q_i
              n_ice(i,j,k) = n_ice(i,j,k)  + fr_n_i
              q_graupel(i,j,k) = q_graupel(i,j,k)  + fr_q_g
              n_graupel(i,j,k) = n_graupel(i,j,k)  + fr_n_g
              q_hail(i,j,k) = q_hail(i,j,k)  + fr_q_h
              n_hail(i,j,k) = n_hail(i,j,k)  + fr_n_h
            ENDIF

            IF (q_graupel(i,j,k) < 0.0) THEN
              q_graupel(i,j,k) = 0.0d0
            END IF
            IF (n_graupel(i,j,k) < 0.0) THEN
              n_graupel(i,j,k) = 0.0d0
            END IF            
            IF (q_hail(i,j,k) < 0.0) THEN
              q_hail(i,j,k) = 0.0d0
            END IF
            IF (n_hail(i,j,k) < 0.0) THEN
              n_hail(i,j,k) = 0.0d0
            END IF            

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE rain_freeze_gamlook

  SUBROUTINE rain_freeze ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Regentropfen                       *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_rain, n_rain, q_graupel, n_graupel,             &
         &                        q_hail, n_hail, q_ice, n_ice,                     &
         &                        dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0
    USE gamma_functions_mp_seifert,    ONLY: incgfct_lower
    USE wolken_konstanten,  ONLY: gfct

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het,&
         &  fr_q_i,fr_n_i,fr_q_g,fr_n_g,fr_q_h,fr_n_h,n_0,lam,xmax_ice,xmax_gr,fr_q_tmp,fr_n_tmp
    DOUBLE PRECISION, SAVE      :: coeff_z
    DOUBLE PRECISION, PARAMETER :: a_HET = 6.5d-1 ! Messung nach Barklie and Gokhale (PK S.350)
    DOUBLE PRECISION, PARAMETER :: b_HET = 2.0d+2 ! Messung nach Barklie and Gokhale (PK S.350)
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20


    IF (firstcall.NE.1) THEN
      !..Koeff. fuer Reflektivitaet Z (2. Moment)
      coeff_z = moment_gamma(rain,2)
      firstcall = 1       
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD rain_freeze" 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Z von Regentropfen: coeff_z= ",coeff_z
      ENDIF
    ELSE IF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD rain_freeze" 
    ENDIF

    xmax_ice = (D_rainfrz_ig/rain%a_geo)**(1.0d0/rain%b_geo)
    xmax_gr  = (D_rainfrz_gh/rain%a_geo)**(1.0d0/rain%b_geo)

    !..Test auf Schmelzen oder Gefrieren von Regentropfen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          q_r = q_rain(i,j,k)
          n_r = n_rain(i,j,k)

          IF (T_a < T_freeze) THEN
            IF (q_r <= q_krit_fr) THEN
              IF (T_a < T_f) THEN
                fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                fr_n = n_r
                fr_n_i= n_r
                fr_q_i= q_r
                fr_n_g= 0.0
                fr_q_g= 0.0
                fr_n_h= 0.0
                fr_q_h= 0.0
! UB_20080220>
                fr_n_tmp = 1.0
                fr_q_tmp = 1.0
! <UB_20080220
              ELSE
                fr_q = 0.0
                fr_n = 0.0
                fr_n_i= 0.0
                fr_q_i= 0.0
                fr_n_g= 0.0
                fr_q_g= 0.0
                fr_n_h= 0.0
                fr_q_h= 0.0
! UB_20080220>
                fr_n_tmp = 0.0
                fr_q_tmp = 0.0
! <UB_20080220
              END IF
            ELSE
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              n_r = q_r / x_r
              IF (T_a < T_f) THEN            !..Nur Eis
                fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
                fr_n = n_r

                ! ub>> Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                !      oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                !      bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                !      und von xmax_gr bis unendlich (--> Hagel).
                
                lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                fr_n_i = n_0/(rain%mu*lam**((rain%nu+1.0)/rain%mu))* &
                    incgfct_lower((rain%nu+1.0)/rain%mu, lam*xmax_ice**rain%mu)
                fr_q_i = n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                    incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                fr_n_g = n_0/(rain%mu*lam**((rain%nu+1.0)/rain%mu))* &
                    incgfct_lower((rain%nu+1.0)/rain%mu, lam*xmax_gr**rain%mu)
                fr_q_g = n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                    incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_gr**rain%mu)
                
                fr_n_h = fr_n - fr_n_g
                fr_q_h = fr_q - fr_q_g
                fr_n_g = fr_n_g - fr_n_i
                fr_q_g = fr_q_g - fr_q_i
                fr_n_tmp = n_r/MAX(fr_n,n_r)
                fr_q_tmp = q_r/MAX(fr_q,q_r)

              ELSE                           !..Heterogenes Gefrieren
                j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.d0) / rho_w * dt

                ! ub>> Je nach Groesse werden die gefrorenen Regentropfen dem Wolkeneis zugeschlagen
                !      oder dem Graupel oder Hagel. Hierzu erfolgt eine partielle Integration des Spektrums von 0
                !      bis zu einer ersten Trennmasse xmax_ice (--> Eis), von dort bis zu xmax_gr (--> Graupel)
                !      und von xmax_gr bis unendlich (--> Hagel).
                
                IF (j_het >= 1d-20) THEN
                  fr_n  = j_het * q_r
                  fr_q  = j_het * q_r * x_r * coeff_z
                  
                  lam = ( gfct((rain%nu+1.0)/rain%mu) / gfct((rain%nu+2.0)/rain%mu) * x_r)**(-rain%mu)
                  n_0 = rain%mu * n_r * lam**((rain%nu+1.0)/rain%mu) / gfct((rain%nu+1.0)/rain%mu)
                  fr_n_i = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                      incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_ice**rain%mu)
                  fr_q_i = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                      incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_ice**rain%mu)
                  fr_n_g = j_het * n_0/(rain%mu*lam**((rain%nu+2.0)/rain%mu))* &
                      incgfct_lower((rain%nu+2.0)/rain%mu, lam*xmax_gr**rain%mu)
                  fr_q_g = j_het * n_0/(rain%mu*lam**((rain%nu+3.0)/rain%mu))* &
                      incgfct_lower((rain%nu+3.0)/rain%mu, lam*xmax_gr**rain%mu)
                  
                  fr_n_h = fr_n - fr_n_g
                  fr_q_h = fr_q - fr_q_g
                  fr_n_g = fr_n_g - fr_n_i
                  fr_q_g = fr_q_g - fr_q_i
                  fr_n_tmp = n_r/MAX(fr_n,n_r)
                  fr_q_tmp = q_r/MAX(fr_q,q_r)
                ELSE
                  fr_n= 0.0
                  fr_q= 0.0
                  fr_n_i= 0.0
                  fr_q_i= 0.0
                  fr_n_g= 0.0
                  fr_q_g= 0.0
! UB_20080212>
                  fr_n_h= 0.0
                  fr_q_h= 0.0
! <UB_20080212
                  fr_n_tmp = 0.0
                  fr_q_tmp = 0.0
                END IF
                
              END IF

              fr_n = fr_n * fr_n_tmp
              fr_q = fr_q * fr_q_tmp

              fr_n_i = fr_n_i * fr_n_tmp
              fr_n_g = fr_n_g * fr_n_tmp
              fr_n_h = fr_n_h * fr_n_tmp
              fr_q_i = fr_q_i * fr_q_tmp
              fr_q_g = fr_q_g * fr_q_tmp
              fr_q_h = fr_q_h * fr_q_tmp

            END IF

            !..Berechnung der H2O-Komponenten

            q_rain(i,j,k) = q_rain(i,j,k) - fr_q
            n_rain(i,j,k) = n_r - fr_n

            IF (q_rain(i,j,k) < 0.0) THEN
              !WRITE (*,*) 'SEIFERT RAIN_FREEZE: Qrain < 0.0, ', i,j,k, T_a, q_r, j_het, fr_q, fr_q_tmp
              q_rain(i,j,k) = 0.0d0
            END IF
            IF (n_rain(i,j,k) < 0.0) THEN
              !WRITE (*,*) 'SEIFERT RAIN_FREEZE: Nrain < 0.0, ', i,j,k, T_a, n_r, j_het, fr_n, fr_n_tmp
              n_rain(i,j,k) = 0.0d0
            END IF            

            IF (ice_typ < 2) THEN 
              ! ohne Hagelklasse,  gefrierender Regen wird Eis oder Graupel
              q_ice(i,j,k) = q_ice(i,j,k)  + fr_q_i
              n_ice(i,j,k) = n_ice(i,j,k)  + fr_n_i
              q_graupel(i,j,k) = q_graupel(i,j,k)  + fr_q_h + fr_q_g
              n_graupel(i,j,k) = n_graupel(i,j,k)  + fr_n_h + fr_n_g
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,60) = fr_q_i
                dqdt(i,j,k,28) = fr_q_g
              END IF
#endif
            ELSE
              ! mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel
              q_ice(i,j,k) = q_ice(i,j,k)  + fr_q_i
              n_ice(i,j,k) = n_ice(i,j,k)  + fr_n_i
              q_graupel(i,j,k) = q_graupel(i,j,k)  + fr_q_g
              n_graupel(i,j,k) = n_graupel(i,j,k)  + fr_n_g
              q_hail(i,j,k) = q_hail(i,j,k)  + fr_q_h
              n_hail(i,j,k) = n_hail(i,j,k)  + fr_n_h
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,60) = fr_q_i
                dqdt(i,j,k,28) = fr_q_g
                dqdt(i,j,k,41) = fr_q_h
              END IF
#endif
            ENDIF

            IF (q_ice(i,j,k) < 0.0) THEN
!              WRITE (*,*) 'SEIFERT RAIN_FREEZE: Qice < 0.0, ', i,j,k, T_a, q_ice(i,j,k), j_het, fr_q_i, fr_q_tmp
              q_ice(i,j,k) = 0.0d0
            END IF
            IF (n_ice(i,j,k) < 0.0) THEN
!              WRITE (*,*) 'SEIFERT RAIN_FREEZE: Nice < 0.0, ', i,j,k, T_a, n_ice(i,j,k), j_het, fr_n_i, fr_n_tmp
              n_ice(i,j,k) = 0.0d0
            END IF            

            IF (q_graupel(i,j,k) < 0.0) THEN
!              WRITE (*,*) 'SEIFERT RAIN_FREEZE: Qgraupel < 0.0, ', i,j,k, T_a, q_graupel(i,j,k), j_het, fr_q_g, fr_q_tmp
              q_graupel(i,j,k) = 0.0d0
            END IF
            IF (n_graupel(i,j,k) < 0.0) THEN
!              WRITE (*,*) 'SEIFERT RAIN_FREEZE: Ngraupel < 0.0, ', i,j,k, T_a, n_graupel(i,j,k), j_het, fr_n_g, fr_n_tmp
              n_graupel(i,j,k) = 0.0d0
            END IF            

            IF (q_hail(i,j,k) < 0.0) THEN
!              WRITE (*,*) 'SEIFERT RAIN_FREEZE: Qhail < 0.0, ', i,j,k, T_a, q_hail(i,j,k), j_het, fr_q, fr_q_i, fr_q_g, fr_q_h, fr_q_tmp
              q_hail(i,j,k) = 0.0d0
            END IF
            IF (n_hail(i,j,k) < 0.0) THEN
!              WRITE (*,*) 'SEIFERT RAIN_FREEZE: Nhail < 0.0, ', i,j,k, T_a, n_hail(i,j,k), j_het, fr_n, fr_n_i, fr_n_g, fr_n_h, fr_n_tmp
              n_hail(i,j,k) = 0.0d0
            END IF            

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE rain_freeze

  SUBROUTINE rain_freeze_old ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Gefrieren von Regentropfen                       *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_rain, n_rain, q_graupel, n_graupel, q_hail, n_hail, &
         &                        dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: fr_q,fr_n,T_a,q_r,x_r,n_r,j_het
    DOUBLE PRECISION, SAVE      :: coeff_z
    DOUBLE PRECISION, PARAMETER :: a_HET = 6.5d-1 ! Messung nach Barklie and Gokhale (PK S.350)
    DOUBLE PRECISION, PARAMETER :: b_HET = 2.0d+2 ! Messung nach Barklie and Gokhale (PK S.350)
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q, rate_n

    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), &
         rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)  
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD rain_freeze",4)
    rate_q = 0.0d0
    rate_n = 0.0d0

    IF (firstcall.NE.1) THEN
      !..Koeff. fuer Reflektivitaet Z (2. Moment)
      coeff_z = moment_gamma(rain,2)
      firstcall = 1       
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD rain_freeze" 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Z von Regentropfen: coeff_z= ",coeff_z
      ENDIF
    ELSE IF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD rain_freeze" 
    ENDIF

    !..Test auf Schmelzen oder Gefrieren von Regentropfen
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          q_r = q_rain(i,j,k)
          n_r = n_rain(i,j,k)
          IF (T_a < T_freeze .AND. q_r > q_krit_fr) THEN
            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
            IF (T_a < T_f) THEN            !..Nur Eis
              fr_q = q_r                  !  Ausfrieren unterhalb T_f \approx -40 C
              fr_n = n_r
            ELSE                           !..Heterogenes Gefrieren
              j_het = MAX(b_HET * ( EXP( a_HET * (T_3 - T_a)) - 1.0 ),0.d0) / rho_w * dt
              fr_n  = j_het * q_r
              fr_q  = j_het * q_r * x_r * coeff_z
            END IF
            fr_q  = MIN(fr_q,q_r)
            fr_n  = MIN(fr_n,n_r)

            rate_q(i,j,k) = fr_q
            rate_n(i,j,k) = fr_n

          END IF
        END DO
      END DO
    END DO

    !..Berechnung der H2O-Komponenten
    q_rain = q_rain - rate_q
    n_rain = n_rain - rate_n

    IF (ice_typ < 2) THEN 
      ! ohne Hagelklasse,  gefrierender Regen wird Graupel
      q_graupel = q_graupel  + rate_q
      n_graupel = n_graupel  + rate_n
#ifdef SAVE_CONVERSIONRATES
      IF (speichere_dqdt) THEN
        dqdt(:,:,:,28) = rate_q
      END IF
#endif
    ELSE
      ! mit Hagelklasse, gefrierender Regen wird Hagel
      q_hail = q_hail  + rate_q
      n_hail = n_hail  + rate_n
#ifdef SAVE_CONVERSIONRATES
      IF (speichere_dqdt) THEN
        dqdt(:,:,:,41) = rate_q
      END IF
#endif
    ENDIF

    DEALLOCATE(rate_q, rate_n)

  END SUBROUTINE rain_freeze_old

  SUBROUTINE ice_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Selbsteinfangs der Eispartikel                          *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_ice, n_ice, q_snow, n_snow, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..absolute Temperatur
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i,e_coll,x_conv
    DOUBLE PRECISION            :: self_n,self_q,weight
    DOUBLE PRECISION            :: delta_n_11,delta_n_12,delta_n_22
    DOUBLE PRECISION            :: delta_q_11,delta_q_12,delta_q_22
    DOUBLE PRECISION            :: theta_n_11,theta_n_12,theta_n_22
    DOUBLE PRECISION            :: theta_q_11,theta_q_12,theta_q_22
    DOUBLE PRECISION,SAVE       :: delta_n,delta_q
    DOUBLE PRECISION,SAVE       :: theta_n,theta_q

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD ice_selfcollection " 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(ice,ice,0)
      delta_n_12 = coll_delta_12(ice,ice,0)
      delta_n_22 = coll_delta_22(ice,ice,0)
      delta_q_11 = coll_delta_11(ice,ice,0) 
      delta_q_12 = coll_delta_12(ice,ice,1)
      delta_q_22 = coll_delta_22(ice,ice,1)

      theta_n_11 = coll_theta_11(ice,ice,0)
      theta_n_12 = coll_theta_12(ice,ice,0)
      theta_n_22 = coll_theta_22(ice,ice,0)
      theta_q_11 = coll_theta_11(ice,ice,0)
      theta_q_12 = coll_theta_12(ice,ice,1)
      theta_q_22 = coll_theta_22(ice,ice,1)

      delta_n = delta_n_11 + delta_n_12 + delta_n_22
      delta_q = delta_q_11 + delta_q_12 + delta_q_22
      theta_n = theta_n_11 - theta_n_12 + theta_n_22
      theta_q = theta_q_11 - theta_q_12 + theta_q_22

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:        a_ice      = ",ice%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_ice      = ",ice%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Eis:            alf_ice    = ",ice%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_ice    = ",ice%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer selfcollection ice:       delta_n_11 = ",delta_n_11
        WRITE (6,'(A,D10.3)') "                                        delta_n_12 = ",delta_n_12
        WRITE (6,'(A,D10.3)') "                                        delta_n_22 = ",delta_n_22
        WRITE (6,'(A,D10.3)') "                                        delta_n    = ",delta_n
        WRITE (6,'(A,D10.3)') "                                        theta_n_11 = ",theta_n_11
        WRITE (6,'(A,D10.3)') "                                        theta_n_12 = ",theta_n_12
        WRITE (6,'(A,D10.3)') "                                        theta_n_22 = ",theta_n_22
        WRITE (6,'(A,D10.3)') "                                        theta_n    = ",theta_n
        WRITE (6,'(A,D10.3)') "                                        delta_q_11 = ",delta_q_11
        WRITE (6,'(A,D10.3)') "                                        delta_q_12 = ",delta_q_12
        WRITE (6,'(A,D10.3)') "                                        delta_q_22 = ",delta_q_22
        WRITE (6,'(A,D10.3)') "                                        delta_q    = ",delta_q
        WRITE (6,'(A,D10.3)') "                                        theta_q_11 = ",theta_q_11 
        WRITE (6,'(A,D10.3)') "                                        theta_q_12 = ",theta_q_12
        WRITE (6,'(A,D10.3)') "                                        theta_q_22 = ",theta_q_22
        WRITE (6,'(A,D10.3)') "                                        theta_q    = ",theta_q
      END IF
      firstcall = 1
    ENDIF

    x_conv = (D_conv_ii/snow%a_geo)**(1./snow%b_geo)

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)                                   !..Fluessigwassergehalt in SI
          n_i = n_ice(i,j,k)                                   !..Anzahldichte in SI
          x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse in SI     
          D_i = ice%a_geo * x_i**ice%b_geo                     !..mittlerer Durchmesser
          IF ( n_i > 0.0 .AND. q_i > q_krit_ii .AND. D_i > D_krit_ii ) THEN

            T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
              !   (siehe auch Straka, 1989; S. 53)
              e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2d0)
              !.. Temperaturabhaengige Efficiency nach Lin et al. (1983)
              !e_coll = MIN(exp(0.09*(T_a-T_3)),1.0d0)
              !e_coll = MAX(e_ii,MIN(exp(0.09*(T_a-T_3)),1.0d0))
            END IF


            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)      !..mittlere Sedimentationsgeschw.

            self_n = pi * 0.25d0 * e_coll * delta_n * n_i * n_i * D_i * D_i &
                 & * ( theta_n * v_i * v_i + 2.0*ice_s_vel**2 )**0.5 * dt

            self_q = pi * 0.25d0 * e_coll * delta_q * n_i * q_i * D_i * D_i &
                 & * ( theta_q * v_i * v_i + 2.0*ice_s_vel**2 )**0.5 * dt

            self_q = MIN(self_q,q_i)
            self_n = MIN(MIN(self_n,self_q/x_conv),n_i)

            q_ice(i,j,k)  = q_ice(i,j,k)  - self_q
            q_snow(i,j,k) = q_snow(i,j,k) + self_q

            n_ice(i,j,k)  = n_ice(i,j,k)  - self_n
            n_snow(i,j,k) = n_snow(i,j,k) + self_n / 2.0

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,54) = self_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ice_selfcollection
 
  SUBROUTINE snow_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Selbsteinfangs der Eispartikel                          *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_snow, n_snow, dt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s,e_ss,e_coll
    DOUBLE PRECISION            :: self_n,self_q,weight
    DOUBLE PRECISION            :: delta_n_11,delta_n_12
    DOUBLE PRECISION            :: delta_q_11,delta_q_12
    DOUBLE PRECISION            :: theta_n_11,theta_n_12
    DOUBLE PRECISION            :: theta_q_11,theta_q_12
    DOUBLE PRECISION,SAVE       :: delta_n,delta_q
    DOUBLE PRECISION,SAVE       :: theta_n,theta_q

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD snow_selfcollection " 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(snow,snow,0)
      delta_n_12 = coll_delta_12(snow,snow,0)

      theta_n_11 = coll_theta_11(snow,snow,0)
      theta_n_12 = coll_theta_12(snow,snow,0)

      delta_n = (2.0*delta_n_11 + delta_n_12)
      theta_n = (2.0*theta_n_11 - theta_n_12)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:     a_snow      = ",snow%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_snow      = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:         alf_snow    = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_snow    = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer selfcollection snow:      delta_n_11 = ",delta_n_11
        WRITE (6,'(A,D10.3)') "                                        delta_n_12 = ",delta_n_12
        WRITE (6,'(A,D10.3)') "                                        delta_n    = ",delta_n
        WRITE (6,'(A,D10.3)') "                                        theta_n_11 = ",theta_n_11
        WRITE (6,'(A,D10.3)') "                                        theta_n_12 = ",theta_n_12
        WRITE (6,'(A,D10.3)') "                                        theta_n    = ",theta_n
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF ( q_s > q_krit ) THEN
            T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
              !   (siehe auch Straka, 1989; S. 53)
              !e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2d0)
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              !e_coll = MIN(exp(0.09*(T_a-T_3)),1.0d0)
              e_coll = MAX(0.1d0,MIN(EXP(0.09*(T_a-T_3)),1.0d0))
            ENDIF

            n_s = n_snow(i,j,k)                        !..Anzahldichte in SI
            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse in SI     

            D_s = snow%a_geo * x_s**snow%b_geo                     !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

            !self_n = pi/4.0 * e_coll * n_s**2 * delta_n * D_s**2 * &
            !     &                ( theta_n * v_s**2 + 2.0 * snow_s_vel**2 )**0.5 * dt

            ! >>hn
            ! fehlenden Faktor 1/2 ergaenzt                                
            self_n = pi * 0.125d0 * e_coll * n_s * n_s * delta_n * D_s * D_s * &
                 &                ( theta_n * v_s * v_s + 2.0 * snow_s_vel**2 )**0.5 * dt
            ! << hn

            self_n = MIN(self_n,n_s)

            n_snow(i,j,k) = n_snow(i,j,k) - self_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_selfcollection

  SUBROUTINE snow_melting()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Schmelzens der Schneepartikel                           *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_snow, n_snow, q_rain, n_rain, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s,T_a,N_re,D_T,e_a
    DOUBLE PRECISION            :: melt,melt_v,melt_h,melt_n,melt_q
    DOUBLE PRECISION            :: fh_n,fv_n
    DOUBLE PRECISION            :: fh_q,fv_q
    DOUBLE PRECISION, SAVE      :: a_melt_n,b_melt_n
    DOUBLE PRECISION, SAVE      :: a_melt_q,b_melt_q

    DOUBLE PRECISION, PARAMETER :: eps   = 1.d-20

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD snow_melting " 
    END IF

    IF (firstcall.NE.1) THEN
      a_melt_n = vent_coeff_a(snow,0)
      b_melt_n = vent_coeff_b(snow,0)
      a_melt_q = vent_coeff_a(snow,1)
      b_melt_q = vent_coeff_b(snow,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:    a_geo  = ",snow%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_geo  = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:        a_vel  = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                       b_vel  = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Ventilation von Schnee:  a_ven  = ",snow%a_ven 
        WRITE (6,'(A,D10.3)') "                                       b_ven  = ",snow%b_ven 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Schmelzen:               a_melt_n = ",a_melt_n
        WRITE (6,'(A,D10.3)') "                                       b_melt_n = ",b_melt_n
        WRITE (6,'(A,D10.3)') "                                       a_melt_q = ",a_melt_q
        WRITE (6,'(A,D10.3)') "                                       b_melt_q = ",b_melt_q
      ENDIF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)
          q_s = q_snow(i,j,k)                                    !..Fluessigwassergehalt in SI
          IF (T_a > T_3 .AND. q_s > 0.0) THEN
            e_a = e_ws(T_a)                                     !..Saettigungsdampfdruck
            n_s = n_snow(i,j,k)                                 !..Anzahldichte in SI

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max) !..mittlere Masse in SI     

            D_s = snow%a_geo * x_s**snow%b_geo                  !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            N_re = v_s * D_s / nu_l                             !..mittlere Reynoldszahl
            fv_q = a_melt_q + b_melt_q * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf
!            fv_n = a_melt_n + b_melt_n * N_sc**n_f * N_re**m_f  !..mittlerer Vent.Koeff. Wasserdampf

! UB_20081111:           D_T  = K_T/(cp * (rho_0(i,j,k)+rho_g(i,j,k)))
! UB_20081125:            D_T  = K_T / (cp * rho_0(i,j,k))
! UB_20081125:            fh_q = D_T/D_v * fv_q  ! <-- gives wrong values!
! UB_20081125: After formulas of Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q is approx. 1.05
!              for a wide temperature- and pressure range:
            fh_q = 1.05 * fv_q
!            fh_n = D_T/D_v * fv_n

            melt   = 2.0*pi/L_ew * D_s * n_s * dt

            melt_h = melt * K_T * (T_a - T_3)
            melt_v = melt * D_v*L_wd/R_d * (e_a/T_a - e_3/T_3)

            melt_q = (melt_h * fh_q + melt_v * fv_q)
!            melt_n = (melt_h * fh_n + melt_v * fv_n) / x_s

! ub>> setzte melt_n so, dass x_s beim Schmelzvorgang erhalten bleibt:
            melt_n = MIN(MAX( (melt_q - q_s) / x_s + n_s, 0.0d0), n_s)
 
            melt_q = MIN(q_s,MAX(melt_q,0.d0))
            melt_n = MIN(n_s,MAX(melt_n,0.d0))

            IF (T_a - T_3 > 10.0) THEN
              melt_q = q_s
              melt_n = n_s
            ENDIF

            q_snow(i,j,k) = q_snow(i,j,k) - melt_q
            q_rain(i,j,k) = q_rain(i,j,k) + melt_q

            n_snow(i,j,k) = n_snow(i,j,k) - melt_n
            n_rain(i,j,k) = n_rain(i,j,k) + melt_n

            n_snow(i,j,k) = MAX(n_snow(i,j,k), q_snow(i,j,k)/snow%x_max)
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,30) = melt_q
            END IF
#endif
            ! ub<<

            !WRF!T(i,j,k) = T(i, j, k) - L_ew / cp * melt_q / (rho_0(i,j,k)+rho_g(i,j,k)) 
            !WRF!p(i,j,k) = p(i, j, k) - L_ew / cv * melt_q * R_l
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_melting

  SUBROUTINE graupel_snow_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Graupel und Schnee                      *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_graupel, n_graupel, q_snow, n_snow, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s
    DOUBLE PRECISION            :: coll_n,coll_q,e_coll
    DOUBLE PRECISION, SAVE      :: delta_n_gg,delta_n_gs,delta_n_ss
    DOUBLE PRECISION, SAVE      :: delta_q_gg,delta_q_gs,delta_q_ss
    DOUBLE PRECISION, SAVE      :: theta_n_gg,theta_n_gs,theta_n_ss
    DOUBLE PRECISION, SAVE      :: theta_q_gg,theta_q_gs,theta_q_ss

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION, PARAMETER :: e_gs = 0.05 ! Alter Wert 0.001

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD graupel_snow_collection:"
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,snow,0)
      delta_n_gs = coll_delta_12(graupel,snow,0)
      delta_n_ss = coll_delta_22(graupel,snow,0)
      delta_q_gg = coll_delta_11(graupel,snow,0) 
      delta_q_gs = coll_delta_12(graupel,snow,1)
      delta_q_ss = coll_delta_22(graupel,snow,1)

      theta_n_gg = coll_theta_11(graupel,snow,0)
      theta_n_gs = coll_theta_12(graupel,snow,0)
      theta_n_ss = coll_theta_22(graupel,snow,0)
      theta_q_gg = coll_theta_11(graupel,snow,0)
      theta_q_gs = coll_theta_12(graupel,snow,1)
      theta_q_ss = coll_theta_22(graupel,snow,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:   a_graupel   = ",graupel%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_graupel   = ",graupel%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Graupel:       alf_graupel = ",graupel%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_graupel = ",graupel%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee.:   a_snow   = ",snow%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_snow   = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:        alf_snow = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_snow = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Koag. graupel-snow:      delta_n_gg = ",delta_n_gg
        WRITE (6,'(A,D10.3)') "                                       delta_n_gs = ",delta_n_gs
        WRITE (6,'(A,D10.3)') "                                       delta_n_ss = ",delta_n_ss
        WRITE (6,'(A,D10.3)') "                                       theta_n_gg = ",theta_n_gg
        WRITE (6,'(A,D10.3)') "                                       theta_n_gs = ",theta_n_gs
        WRITE (6,'(A,D10.3)') "                                       theta_n_ss = ",theta_n_ss
        WRITE (6,'(A,D10.3)') "                                       delta_q_gg = ",delta_q_gg
        WRITE (6,'(A,D10.3)') "                                       delta_q_gs = ",delta_q_gs
        WRITE (6,'(A,D10.3)') "                                       delta_q_ss = ",delta_q_ss
        WRITE (6,'(A,D10.3)') "                                       theta_q_gg = ",theta_q_gg
        WRITE (6,'(A,D10.3)') "                                       theta_q_gs = ",theta_q_gs
        WRITE (6,'(A,D10.3)') "                                       theta_q_ss = ",theta_q_ss
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)                                    !..Fluessigwassergehalt in SI
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_s > q_krit .AND. q_g > q_krit) THEN
            T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0d0)
            ENDIF

            n_s = n_snow(i,j,k)                                       !..Anzahldichte
            n_g = n_graupel(i,j,k)                                    !..Anzahldichte

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)       !..mittlere Masse
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse

            D_s = snow%a_geo * x_s**snow%b_geo                        !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.

            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_g * n_s * e_coll * dt & 
                 &   * (delta_n_gg * D_g**2 + delta_n_gs * D_g*D_s + delta_n_ss * D_s**2) &
                 &   * (theta_n_gg * v_g**2 - theta_n_gs * v_g*v_s + theta_n_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_q = pi/4.0 * n_g * q_s * e_coll * dt & 
                 &   * (delta_q_gg * D_g**2 + delta_q_gs * D_g*D_s + delta_q_ss * D_s**2) &
                 &   * (theta_q_gg * v_g**2 - theta_q_gs * v_g*v_s + theta_q_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            q_graupel(i,j,k) = q_graupel(i,j,k) + coll_q
            q_snow(i,j,k)    = q_snow(i,j,k)    - coll_q
            n_snow(i,j,k)    = n_snow(i,j,k)    - coll_n

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,10) = coll_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_snow_collection


  ! hn>>
  SUBROUTINE hail_snow_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Hagel und Schnee                      *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_hail, n_hail, q_snow, n_snow, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s
    DOUBLE PRECISION            :: coll_n,coll_q,e_coll
    DOUBLE PRECISION, SAVE      :: delta_n_hh,delta_n_hs,delta_n_ss
    DOUBLE PRECISION, SAVE      :: delta_q_hh,delta_q_hs,delta_q_ss
    DOUBLE PRECISION, SAVE      :: theta_n_hh,theta_n_hs,theta_n_ss
    DOUBLE PRECISION, SAVE      :: theta_q_hh,theta_q_hs,theta_q_ss

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION, PARAMETER :: e_hs = 0.05 ! Alter Wert 0.001

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD hail_snow_collection:"
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,snow,0)
      delta_n_hs = coll_delta_12(hail,snow,0)
      delta_n_ss = coll_delta_22(hail,snow,0)
      delta_q_hh = coll_delta_11(hail,snow,0) 
      delta_q_hs = coll_delta_12(hail,snow,1)
      delta_q_ss = coll_delta_22(hail,snow,1)

      theta_n_hh = coll_theta_11(hail,snow,0)
      theta_n_hs = coll_theta_12(hail,snow,0)
      theta_n_ss = coll_theta_22(hail,snow,0)
      theta_q_hh = coll_theta_11(hail,snow,0)
      theta_q_hs = coll_theta_12(hail,snow,1)
      theta_q_ss = coll_theta_22(hail,snow,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hail:      a_hail   = ",hail%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_hail   = ",hail%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Hail:          alf_hail = ",hail%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_hail = ",hail%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee.:   a_snow   = ",snow%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_snow   = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:        alf_snow = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_snow = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Koag. hail-snow:         delta_n_hh = ",delta_n_hh
        WRITE (6,'(A,D10.3)') "                                       delta_n_hs = ",delta_n_hs
        WRITE (6,'(A,D10.3)') "                                       delta_n_ss = ",delta_n_ss
        WRITE (6,'(A,D10.3)') "                                       theta_n_hh = ",theta_n_hh
        WRITE (6,'(A,D10.3)') "                                       theta_n_hs = ",theta_n_hs
        WRITE (6,'(A,D10.3)') "                                       theta_n_ss = ",theta_n_ss
        WRITE (6,'(A,D10.3)') "                                       delta_q_hh = ",delta_q_hh
        WRITE (6,'(A,D10.3)') "                                       delta_q_hs = ",delta_q_hs
        WRITE (6,'(A,D10.3)') "                                       delta_q_ss = ",delta_q_ss
        WRITE (6,'(A,D10.3)') "                                       theta_q_hh = ",theta_q_hh
        WRITE (6,'(A,D10.3)') "                                       theta_q_hs = ",theta_q_hs
        WRITE (6,'(A,D10.3)') "                                       theta_q_ss = ",theta_q_ss
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)                                 !..Fluessigwassergehalt in SI
          q_h = q_hail(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_s > q_krit .AND. q_h > q_krit) THEN
            T_a = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0d0)
            ENDIF

            n_s = n_snow(i,j,k)                                    !..Anzahldichte
            n_h = n_hail(i,j,k)                                    !..Anzahldichte

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max) !..mittlere Masse
            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max) !..mittlere Masse

            D_s = snow%a_geo * x_s**snow%b_geo                        !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)       !..mittlere Sedimentationsgeschw.

            D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_h * n_s * e_coll * dt & 
                 &   * (delta_n_hh * D_h**2 + delta_n_hs * D_h*D_s + delta_n_ss * D_s**2) &
                 &   * (theta_n_hh * v_h**2 - theta_n_hs * v_h*v_s + theta_n_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_q = pi/4.0 * n_h * q_s * e_coll * dt & 
                 &   * (delta_q_hh * D_h**2 + delta_q_hs * D_h*D_s + delta_q_ss * D_s**2) &
                 &   * (theta_q_hh * v_h**2 - theta_q_hs * v_h*v_s + theta_q_ss * v_s**2  &
                 &     +snow_s_vel**2)**0.5

            coll_n = MIN(n_s,coll_n)
            coll_q = MIN(q_s,coll_q)

            q_hail(i,j,k) = q_hail(i,j,k) + coll_q
            q_snow(i,j,k) = q_snow(i,j,k) - coll_q
            n_snow(i,j,k) = n_snow(i,j,k) - coll_n

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,43) = coll_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_snow_collection
  ! <<hn


  SUBROUTINE graupel_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Graupel und Schnee                      *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_graupel, n_graupel, q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i
    DOUBLE PRECISION            :: coll_n,coll_q,e_coll
    DOUBLE PRECISION, SAVE      :: delta_n_gg,delta_n_gi,delta_n_ii
    DOUBLE PRECISION, SAVE      :: delta_q_gg,delta_q_gi,delta_q_ii
    DOUBLE PRECISION, SAVE      :: theta_n_gg,theta_n_gi,theta_n_ii
    DOUBLE PRECISION, SAVE      :: theta_q_gg,theta_q_gi,theta_q_ii

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION, PARAMETER :: e_gi = 0.05

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD graupel_ice_collection:"
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_gg = coll_delta_11(graupel,ice,0)
      delta_n_gi = coll_delta_12(graupel,ice,0)
      delta_n_ii = coll_delta_22(graupel,ice,0)
      delta_q_gg = coll_delta_11(graupel,ice,0) 
      delta_q_gi = coll_delta_12(graupel,ice,1)
      delta_q_ii = coll_delta_22(graupel,ice,1)

      theta_n_gg = coll_theta_11(graupel,ice,0)
      theta_n_gi = coll_theta_12(graupel,ice,0)
      theta_n_ii = coll_theta_22(graupel,ice,0)
      theta_q_gg = coll_theta_11(graupel,ice,0)
      theta_q_gi = coll_theta_12(graupel,ice,1)
      theta_q_ii = coll_theta_22(graupel,ice,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:   a_graupel   = ",graupel%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_graupel   = ",graupel%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Graupel:       alf_graupel = ",graupel%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_graupel = ",graupel%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee.:   a_ice   = ",ice%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_ice   = ",ice%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:        alf_ice = ",ice%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_ice = ",ice%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Koag. graupel-ice:      delta_n_gg = ",delta_n_gg
        WRITE (6,'(A,D10.3)') "                                       delta_n_gi = ",delta_n_gi
        WRITE (6,'(A,D10.3)') "                                       delta_n_ii = ",delta_n_ii
        WRITE (6,'(A,D10.3)') "                                       theta_n_gg = ",theta_n_gg
        WRITE (6,'(A,D10.3)') "                                       theta_n_gi = ",theta_n_gi
        WRITE (6,'(A,D10.3)') "                                       theta_n_ii = ",theta_n_ii
        WRITE (6,'(A,D10.3)') "                                       delta_q_gg = ",delta_q_gg
        WRITE (6,'(A,D10.3)') "                                       delta_q_gi = ",delta_q_gi
        WRITE (6,'(A,D10.3)') "                                       delta_q_ii = ",delta_q_ii
        WRITE (6,'(A,D10.3)') "                                       theta_q_gg = ",theta_q_gg
        WRITE (6,'(A,D10.3)') "                                       theta_q_gi = ",theta_q_gi
        WRITE (6,'(A,D10.3)') "                                       theta_q_ii = ",theta_q_ii
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)                                     !..Fluessigwassergehalt in SI
          q_g = q_graupel(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_i > q_krit .AND. q_g > q_krit) THEN
            T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

            !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0d0)
            END IF

            n_i = n_ice(i,j,k)                                        !..Anzahldichte
            n_g = n_graupel(i,j,k)                                    !..Anzahldichte

            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)         !..mittlere Masse
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse

            D_i = ice%a_geo * x_i**ice%b_geo                          !..mittlerer Durchmesser
            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)         !..mittlere Sedimentationsgeschw.

            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_g * n_i * e_coll * dt & 
                 &   * (delta_n_gg * D_g**2 + delta_n_gi * D_g*D_i + delta_n_ii * D_i**2) &
                 &   * (theta_n_gg * v_g**2 - theta_n_gi * v_g*v_i + theta_n_ii * v_i**2)**0.5

            coll_q = pi/4.0 * n_g * q_i * e_coll * dt & 
                 &   * (delta_q_gg * D_g**2 + delta_q_gi * D_g*D_i + delta_q_ii * D_i**2) &
                 &   * (theta_q_gg * v_g**2 - theta_q_gi * v_g*v_i + theta_q_ii * v_i**2)**0.5

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_graupel(i,j,k) = q_graupel(i,j,k) + coll_q
            q_ice(i,j,k)     = q_ice(i,j,k)     - coll_q
            n_ice(i,j,k)     = n_ice(i,j,k)     - coll_n

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,9) = coll_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_ice_collection

  ! hn>>
  SUBROUTINE hail_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Hagel und Schnee                        *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_hail, n_hail, q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i
    DOUBLE PRECISION            :: coll_n,coll_q,e_coll
    DOUBLE PRECISION, SAVE      :: delta_n_hh,delta_n_hi,delta_n_ii
    DOUBLE PRECISION, SAVE      :: delta_q_hh,delta_q_hi,delta_q_ii
    DOUBLE PRECISION, SAVE      :: theta_n_hh,theta_n_hi,theta_n_ii
    DOUBLE PRECISION, SAVE      :: theta_q_hh,theta_q_hi,theta_q_ii

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION, PARAMETER :: e_hi = 0.05

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD hail_ice_collection:"
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_hh = coll_delta_11(hail,ice,0)
      delta_n_hi = coll_delta_12(hail,ice,0)
      delta_n_ii = coll_delta_22(hail,ice,0)
      delta_q_hh = coll_delta_11(hail,ice,0) 
      delta_q_hi = coll_delta_12(hail,ice,1)
      delta_q_ii = coll_delta_22(hail,ice,1)

      theta_n_hh = coll_theta_11(hail,ice,0)
      theta_n_hi = coll_theta_12(hail,ice,0)
      theta_n_ii = coll_theta_22(hail,ice,0)
      theta_q_hh = coll_theta_11(hail,ice,0)
      theta_q_hi = coll_theta_12(hail,ice,1)
      theta_q_ii = coll_theta_22(hail,ice,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hagel:     a_hail   = ",hail%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_hail   = ",hail%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Hagel:         alf_hail = ",hail%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_hail = ",hail%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee.:   a_ice    = ",ice%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_ice    = ",ice%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:        alf_ice  = ",ice%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_ice  = ",ice%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Koag. hail-ice:          delta_n_hh = ",delta_n_hh
        WRITE (6,'(A,D10.3)') "                                       delta_n_hi = ",delta_n_hi
        WRITE (6,'(A,D10.3)') "                                       delta_n_ii = ",delta_n_ii
        WRITE (6,'(A,D10.3)') "                                       theta_n_hh = ",theta_n_hh
        WRITE (6,'(A,D10.3)') "                                       theta_n_hi = ",theta_n_hi
        WRITE (6,'(A,D10.3)') "                                       theta_n_ii = ",theta_n_ii
        WRITE (6,'(A,D10.3)') "                                       delta_q_hh = ",delta_q_hh
        WRITE (6,'(A,D10.3)') "                                       delta_q_hi = ",delta_q_hi
        WRITE (6,'(A,D10.3)') "                                       delta_q_ii = ",delta_q_ii
        WRITE (6,'(A,D10.3)') "                                       theta_q_hh = ",theta_q_hh
        WRITE (6,'(A,D10.3)') "                                       theta_q_hi = ",theta_q_hi
        WRITE (6,'(A,D10.3)') "                                       theta_q_ii = ",theta_q_ii
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)                                  !..Fluessigwassergehalt in SI
          q_h = q_hail(i,j,k)                                 !..Fluessigwassergehalt in SI
          IF (q_i > q_krit .AND. q_h > q_krit) THEN
            T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MIN(EXP(0.09*(T_a-T_3)),1.0d0)
            ENDIF

            n_i = n_ice(i,j,k)                                     !..Anzahldichte
            n_h = n_hail(i,j,k)                                    !..Anzahldichte

            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)      !..mittlere Masse
            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)    !..mittlere Masse

            D_i = ice%a_geo * x_i**ice%b_geo                       !..mittlerer Durchmesser
            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)      !..mittlere Sedimentationsgeschw.

            D_h = hail%a_geo * x_h**hail%b_geo                  !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_h * n_i * e_coll * dt & 
                 &   * (delta_n_hh * D_h**2 + delta_n_hi * D_h*D_i + delta_n_ii * D_i**2) &
                 &   * (theta_n_hh * v_h**2 - theta_n_hi * v_h*v_i + theta_n_ii * v_i**2)**0.5

            coll_q = pi/4.0 * n_h * q_i * e_coll * dt & 
                 &   * (delta_q_hh * D_h**2 + delta_q_hi * D_h*D_i + delta_q_ii * D_i**2) &
                 &   * (theta_q_hh * v_h**2 - theta_q_hi * v_h*v_i + theta_q_ii * v_i**2)**0.5

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_hail(i,j,k) = q_hail(i,j,k) + coll_q
            q_ice(i,j,k)  = q_ice(i,j,k)  - coll_q
            n_ice(i,j,k)  = n_ice(i,j,k)  - coll_n

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,42) = coll_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE hail_ice_collection
  ! <<hn

  SUBROUTINE snow_ice_collection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Koagulation von Schnee und Eis                          *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_snow, n_snow, q_ice, n_ice, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i
    DOUBLE PRECISION            :: coll_n,coll_q,e_coll
    DOUBLE PRECISION, SAVE      :: delta_n_ss,delta_n_si,delta_n_ii
    DOUBLE PRECISION, SAVE      :: delta_q_ss,delta_q_si,delta_q_ii
    DOUBLE PRECISION, SAVE      :: theta_n_ss,theta_n_si,theta_n_ii
    DOUBLE PRECISION, SAVE      :: theta_q_ss,theta_q_si,theta_q_ii

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION            :: e_si

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD snow_ice_collection:"
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_ss = coll_delta_11(snow,ice,0)
      delta_n_si = coll_delta_12(snow,ice,0)
      delta_n_ii = coll_delta_22(snow,ice,0)
      delta_q_ss = coll_delta_11(snow,ice,0) 
      delta_q_si = coll_delta_12(snow,ice,1)
      delta_q_ii = coll_delta_22(snow,ice,1)

      theta_n_ss = coll_theta_11(snow,ice,0)
      theta_n_si = coll_theta_12(snow,ice,0)
      theta_n_ii = coll_theta_22(snow,ice,0)
      theta_q_ss = coll_theta_11(snow,ice,0)
      theta_q_si = coll_theta_12(snow,ice,1)
      theta_q_ii = coll_theta_22(snow,ice,1)

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:    a_snow   = ",snow%a_geo 
        WRITE (6,'(A,D10.3)') "                                       b_snow   = ",snow%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:        alf_snow = ",snow%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_snow = ",snow%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:       a_ice   = ",ice%a_geo
        WRITE (6,'(A,D10.3)') "                                       b_ice   = ",ice%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Eis:           alf_ice = ",ice%a_vel 
        WRITE (6,'(A,D10.3)') "                                       bet_ice = ",ice%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Koag. snow-ice:          delta_n_ss = ",delta_n_ss
        WRITE (6,'(A,D10.3)') "                                       delta_n_si = ",delta_n_si
        WRITE (6,'(A,D10.3)') "                                       delta_n_ii = ",delta_n_ii
        WRITE (6,'(A,D10.3)') "                                       theta_n_ss = ",theta_n_ss
        WRITE (6,'(A,D10.3)') "                                       theta_n_si = ",theta_n_si
        WRITE (6,'(A,D10.3)') "                                       theta_n_ii = ",theta_n_ii
        WRITE (6,'(A,D10.3)') "                                       delta_q_ss = ",delta_q_ss
        WRITE (6,'(A,D10.3)') "                                       delta_q_si = ",delta_q_si
        WRITE (6,'(A,D10.3)') "                                       delta_q_ii = ",delta_q_ii
        WRITE (6,'(A,D10.3)') "                                       theta_q_ss = ",theta_q_ss
        WRITE (6,'(A,D10.3)') "                                       theta_q_si = ",theta_q_si
        WRITE (6,'(A,D10.3)') "                                       theta_q_ii = ",theta_q_ii
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)
          q_g = q_snow(i,j,k)
          IF (q_i > q_krit .AND. q_g > q_krit) THEN
            T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

            IF (T_a > T_3) THEN
              e_coll = 1.0
            ELSE
              !.. Temperaturabhaengige sticking efficiency nach Lin (1983)
              e_coll = MAX(0.1d0,MIN(EXP(0.09*(T_a-T_3)),1.0d0))
              !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
              !   (siehe auch Straka, 1989; S. 53)
              !e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2d0)
            ENDIF

            n_i = n_ice(i,j,k)                       !..Anzahldichte in SI
            n_g = n_snow(i,j,k)                      !..Anzahldichte in SI

            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)     !..mittlere Masse in SI     
            x_g = MIN(MAX(q_g/(n_g+eps),snow%x_min),snow%x_max)   !..mittlere Masse in SI     

            D_i = ice%a_geo * x_i**ice%b_geo                      !..mittlerer Durchmesser
            v_i = ice%a_vel * x_i**ice%b_vel * rrho_04(i,j,k)     !..mittlere Sedimentationsgeschw.

            D_g = snow%a_geo * x_g**snow%b_geo                    !..mittlerer Durchmesser
            v_g = snow%a_vel * x_g**snow%b_vel * rrho_04(i,j,k)   !..mittlere Sedimentationsgeschw.

            coll_n = pi/4.0 * n_g * n_i * e_coll * dt & 
                 &   * (delta_n_ss * D_g**2 + delta_n_si * D_g*D_i + delta_n_ii * D_i**2) &
                 &   * (theta_n_ss * v_g**2 - theta_n_si * v_g*v_i + theta_n_ii * v_i**2  &
                 &     +snow_s_vel**2 + ice_s_vel**2)**0.5

            coll_q = pi/4.0 * n_g * q_i * e_coll * dt & 
                 &   * (delta_q_ss * D_g**2 + delta_q_si * D_g*D_i + delta_q_ii * D_i**2) &
                 &   * (theta_q_ss * v_g**2 - theta_q_si * v_g*v_i + theta_q_ii * v_i**2  &
                 &     +snow_s_vel**2 + ice_s_vel**2)**0.5

            coll_n = MIN(n_i,coll_n)
            coll_q = MIN(q_i,coll_q)

            q_snow(i,j,k) = q_snow(i,j,k) + coll_q
            q_ice(i,j,k)  = q_ice(i,j,k)  - coll_q
            n_ice(i,j,k)  = n_ice(i,j,k)  - coll_n

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,8) = coll_q
            END IF
#endif
            ! ub<<

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE snow_ice_collection

  SUBROUTINE graupel_selfcollection()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Selbsteinfangs der Eispartikel                          *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_graupel, n_graupel, dt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall

    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g,eff_g
    DOUBLE PRECISION            :: self_n,self_q,weight
    DOUBLE PRECISION            :: delta_n_11,delta_n_12
    DOUBLE PRECISION            :: delta_q_11,delta_q_12
    DOUBLE PRECISION            :: theta_n_11,theta_n_12
    DOUBLE PRECISION            :: theta_q_11,theta_q_12
    DOUBLE PRECISION            :: delta_n,delta_q
    DOUBLE PRECISION            :: theta_n,theta_q
    DOUBLE PRECISION,SAVE       :: coll_n,coll_q

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20
    DOUBLE PRECISION, PARAMETER :: e_gg = 0.1
    DOUBLE PRECISION, PARAMETER :: e_gg_wet = 0.4

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD graupel_selfcollection " 
    END IF

    IF (firstcall.NE.1) THEN
      delta_n_11 = coll_delta_11(graupel,graupel,0)
      delta_n_12 = coll_delta_12(graupel,graupel,0)

      theta_n_11 = coll_theta_11(graupel,graupel,0)
      theta_n_12 = coll_theta_12(graupel,graupel,0)

      delta_n = (2.0*delta_n_11 + delta_n_12)
      theta_n = (2.0*theta_n_11 - theta_n_12)**0.5

      coll_n  = delta_n * theta_n

      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:     a_graupel      = ",graupel%a_geo 
        WRITE (6,'(A,D10.3)') "                                        b_graupel      = ",graupel%b_geo
        WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedi. von Schnee:         alf_graupel    = ",graupel%a_vel 
        WRITE (6,'(A,D10.3)') "                                        bet_graupel    = ",graupel%b_vel 
        WRITE (6,'(A,D10.3)') "  Koeff. fuer selfcollection graupel:   delta_n_11 = ",delta_n_11
        WRITE (6,'(A,D10.3)') "                                        delta_n_12 = ",delta_n_12
        WRITE (6,'(A,D10.3)') "                                        delta_n    = ",delta_n
        WRITE (6,'(A,D10.3)') "                                        theta_n_11 = ",theta_n_11
        WRITE (6,'(A,D10.3)') "                                        theta_n_12 = ",theta_n_12
        WRITE (6,'(A,D10.3)') "                                        theta_n    = ",theta_n
        WRITE (6,'(A,D10.3)') "                                        coll_n     = ",coll_n
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_g = q_graupel(i,j,k)
          n_g = n_graupel(i,j,k)                                       !..Anzahldichte in SI
          IF ( n_g > 0.0 ) THEN
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse in SI
            D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k) !..mittlere Sedimentationsgeschw.
            !self_n = pi/4.0 * e_gg * coll_n * n_g**2 * D_g**2 * v_g * dt
            ! hn>> Korrektur:
            ! Faktor 1/2 ergaenzt
            ! ub>> efficiency von nassem Graupel etwas erhoeht:
            IF (T_0(i,j,k) > T_3) THEN
              self_n = pi/8.0 * e_gg_wet * coll_n * n_g**2 * D_g**2 * v_g * dt
            ELSE
              self_n = pi/8.0 * e_gg * coll_n * n_g**2 * D_g**2 * v_g * dt
            END IF
            ! <<hn

            self_n = MIN(self_n,n_g)

            n_graupel(i,j,k) = n_graupel(i,j,k) - self_n

          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_selfcollection

  SUBROUTINE ice_melting ()
    !*******************************************************************************
    !                                                                              *
    ! Diese Routine behandelt das Schmelzen der Eispartikel                        *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_ice, n_ice, q_cloud, n_cloud, q_rain, n_rain, dt, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ 
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    DOUBLE PRECISION            :: T_a             !..absolute Temperatur
    DOUBLE PRECISION            :: q_c,q_i,x_c,x_i,n_i,weight
    DOUBLE PRECISION            :: melt_q, melt_n
    DOUBLE PRECISION, PARAMETER :: eps = 1.00d-20

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD ice_melting" 
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_i = q_ice(i,j,k)
          T_a = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)
          IF (T_a > T_3 .AND. q_i > 0.0) THEN            !..Nur Fluessigwasser
            n_i = n_ice(i,j,k)
            x_i = MIN(MAX(q_i/(n_i+eps),ice%x_min),ice%x_max)    !..mittlere Masse in SI     

            melt_q = q_i                !  spontanes Schmelzen oberhalb des Tripelpunkts
            melt_n = n_i
            IF (x_i > cloud%x_max) THEN
              weight = 0.0
            ELSE
              weight = 1.0
            ENDIF

            q_ice(i, j, k)   = q_ice(i, j, k) - melt_q
            n_ice(i, j, k)   = n_ice(i, j, k) - melt_n

            q_cloud(i,j,k) = q_cloud(i,j,k) + melt_q * weight
            n_cloud(i,j,k) = n_cloud(i,j,k) + melt_n * weight
            q_rain(i,j,k)  = q_rain(i,j,k)  + melt_q * (1.0 - weight)
            n_rain(i,j,k)  = n_rain(i,j,k)  + melt_n * (1.0 - weight)
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,29) = melt_q
            END IF
#endif
            ! ub<<

            !WRF!T(i,j,k) = T(i,j,k) - L_ew / cp * melt_q / (rho_0(i,j,k)+rho_g(i,j,k))
            !WRF!p(i,j,k) = p(i,j,k) - L_ew / cv * melt_q * R_l

          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ice_melting

  SUBROUTINE rain_evaporation_uli ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstungung von Regentropfen
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, q_cloud, &
         &                        q_rain, n_rain, dt, S_w, dqdt, speichere_dqdt, &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0
    USE gamma_functions_mp_seifert,    ONLY: incgfct_lower
    USE wolken_konstanten,  ONLY: gfct

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a,p_a,e_sw,q_sw,s_sw,R_f,g_d,eva,eva_q,eva_n,a_ld,a_dl
    DOUBLE PRECISION            :: q_d,q_r,n_r,x_r,d_r,v_r,f_v,N_re,x_d,e_d,f_q,f_n
    DOUBLE PRECISION            :: n,n_0,mue,d_m,lam,xmax,xmin,G1,G4
    DOUBLE PRECISION, SAVE      :: c_r               !..Koeff. fuer mittlere Kapazitaet
    DOUBLE PRECISION, SAVE      :: a_n,b_n,a_q,b_q   !..Koeff. fuer mittleren Ventilationkoeff.
    TYPE(PARTICLE) :: rai

    ! ub: neue Parametrisierung der Verdunstung von QNR, nur in Kombination 
    !mit use_mu_Dm_rain_evap = .true. wirksam:
    ! Achtung: Z.Zt. wirkungslos, da auskommentiert!!!!
    LOGICAL, PARAMETER          :: use_evap_qn_uli = .TRUE.

    ! In case of exponential rain DSD (nu=-2/3) in the below rain particle definition,
    ! set this switch to true (in init_seifert()) if intended to use exponential DSD 
    ! also for evaporation:
    LOGICAL :: use_exponentialdsd_rain_evap = .FALSE.

    DOUBLE PRECISION, PARAMETER :: eps = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD rain_evaporation",4)

    rate_q = 0.0
    rate_n = 0.0

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF (cloud_typ > 1 .AND. use_mu_Dm_rain_evap) THEN
      WRITE (*,*) 'ERROR: This option is currently not supported here'
      STOP
    END IF


    IF (firstcall.NE.1) THEN
!      a_q = vent_coeff_a(rain,1)
!      b_q = vent_coeff_b(rain,1)
!      a_n = vent_coeff_a(rain,0)
!      b_n = vent_coeff_b(rain,0)
      a_q = vent_coeff_a(rain,1)
      b_q = vent_coeff_b(rain,1)
      IF (.NOT.use_exponentialdsd_rain_evap) THEN
        a_n = vent_coeff_a(rain,0)
        b_n = vent_coeff_b(rain,0)
      ELSE
        ! For any distribution with nue <= -1/3, 
        ! a_n and b_n are not defined. An exponential 
        ! distribution w.r. to D has nue = -2/3.
        ! This is another manifestation, that the
        ! original SB parameterization of the evaporative
        ! tendency for QNRAIN is wrong!
        ! Therefore, we cheat here and set nu to the smallest
        ! possible value. Remember that QNRAIN-tendendy is wrong anyways...
        rai = rain
        rai%nu = -0.3
        a_n = vent_coeff_a(rai,0)
        b_n = vent_coeff_b(rai,0)
      END IF
      c_r = 1.0 / rain%cap
      firstcall = 1
      IF (isIO()) THEN
        !IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "CLOUD rain_evaporation:"
        WRITE (6,'(A,D10.3)') "     nu    = ",rain%nu
        WRITE (6,'(A,D10.3)') "     mu    = ",rain%mu
        WRITE (6,'(A,D10.3)') "     a_geo = ",rain%a_geo
        WRITE (6,'(A,D10.3)') "     b_geo = ",rain%b_geo
        WRITE (6,'(A,D10.3)') "     a_n   = ",a_n
        WRITE (6,'(A,D10.3)') "     b_n   = ",b_n
        WRITE (6,'(A,D10.3)') "     a_q   = ",a_q
        WRITE (6,'(A,D10.3)') "     b_q   = ",b_q
        WRITE (6,'(A,D10.3)') "     c_r   = ",c_r
      END IF
    ELSEIF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD rain_evaporation " 
    END IF

    G1 = 1              ! gfct(n+1) 
    G4 = pi * rho_w     ! gfct(n+4) * pi/6. * rho_w
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)
          p_a  = p_0(i,j,k) !WRF!+ p(i,j,k)
          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          q_d  = q(i,j,k)
          x_d  = q(i,j,k) / (rho_0(i,j,k))!WRF!+rho_g(i,j,k))
          e_d  = q(i,j,k) * R_d * T_a
          e_sw = e_ws(T_a)
          s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser
!          IF(s_sw < 0.0 .AND. q_r > q_krit .AND. q_cloud(i,j,k) < q_krit)THEN
          IF(s_sw < 0.0 .AND. q_r > 0.0d0 .AND. q_cloud(i,j,k) < q_krit)THEN


            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_v * e_sw) )

            n_r = n_rain(i,j,k)
            IF (cloud_typ <= 1) THEN
              n_0 = 1.00d+07
              n_r = n_0 * (pi*rho_w*n_0/q_r)**(-1./4.)     !..Nur q_r  (n_0 = const) 
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              x_r = MAX(1e3/2e7*q_r,MIN(x_r,1e4/2e5*q_r))
            ELSEIF (use_mu_Dm_rain_evap) THEN

              x_r = q_r / (n_r+eps)
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
              n_r = q_r / x_r

              D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.)
              mue = rain_cmu1*TANH((rain_cmu2*(D_m-rain_cmu3))**rain_cmu5) + rain_cmu4
              IF (D_m < rain_cmu3) THEN
                mue = MAX( (rain%nu+1.0d0)*3.0d0 -1.0d0 , rain_cmu4)
              ELSE
                mue = MAX(MIN( mue, 10.0d0), MAX((rain%nu+1.0d0)*3.0d0-1.0d0, rain_cmu4))
              END IF
              rai = rain
              rai%nu = (mue-2.0)/3.0
              rai%mu = 0.3333333
              a_q = vent_coeff_a(rai,1)
              b_q = vent_coeff_b(rai,1)
              a_n = vent_coeff_a(rai,0)
              b_n = vent_coeff_b(rai,0)
            ELSEIF (use_mu_orig_rain_evap) THEN
              x_r = q_r / (n_r+eps)
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)
              n_r = q_r / x_r
              rai = rain
              rai%mu = 0.3333333
              a_q = vent_coeff_a(rai,1)
              b_q = vent_coeff_b(rai,1)
              a_n = vent_coeff_a(rai,0)
              b_n = vent_coeff_b(rai,0)
            ELSEIF (use_exponentialdsd_rain_evap) THEN
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              n_0 = n_r/G1 * (G4/G1/x_r)**(1./3.)          !..n0
              n_0 = MIN(MAX(250.0d+03,n_0),20000.0d+03)
              n_r = n_0 * (pi*rho_w*n_0/q_r)**(-0.25)     !..Nur q_r  (n_0 = const) 
            ELSE
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              x_r = MAX(1e3/2e7*q_r,MIN(x_r,1e4/2e5*q_r))
              n_r = q_r / x_r
            ENDIF
            d_r = rain%a_geo * x_r**rain%b_geo                 
            v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)

            N_re = v_r * d_r / nu_l                            
            f_v  = N_sc**n_f * N_re**m_f           

            f_n  = a_n + b_n * f_v
            f_q  = a_q + b_q * f_v

            eva   = g_d * n_r * c_r * d_r * s_sw * dt 
!            eva_q = 0.5 * f_q * eva
            eva_q = f_q * eva

            ! UB: Neue Parametrisierung der Verdunstung von Regentropfen-Anzahldichte im Falle
            ! der Anwendung der mue-D-Beziehung von Axel:
            ! eva_n wird bestimmt als Anzahl der Regentropfen (D > 80 mue-m), die waehrend des Zeitschrittes
            ! dt auf einen Durchmesser kleiner als 80 mue-m schrumpfen koennen oder die sowieso
            ! schon kleiner als 80 mue-m sind. Von diesen Tropfen wird angenommen, dass
            ! sie im folgenden sehr schnell verdunsten und dass deren Masse schon in der
            ! vorher berechneten eva_q steckt, was so nicht ganz richtig ist: nur ein Teil ihrer Masse steckt in eva_q ...
            IF (use_mu_Dm_rain_evap .AND. use_evap_qn_uli) THEN
!!!!              xmin = (2.67304d-6/rai%a_geo)**(1.0d0/rai%b_geo) ! xmin von D=80 mue-m
!              xmax = (dmax_evap_von_dt_80(dt,s_sw,T_a,p_a*1e-2)*1.0d-6/rai%a_geo)**(1.0d0/rai%b_geo)
!              lam = ( gfct((rai%nu+1.0)/rai%mu) / gfct((rai%nu+2.0)/rai%mu) * x_r)**(-rai%mu)
!              n_0 = rai%mu * n_r * lam**((rai%nu+1.0)/rai%mu) / gfct((rai%nu+1.0)/rai%mu)
!              eva_n = n_0/(rai%mu*lam**((rai%nu+1.0)/rai%mu))* &
!                   incgfct_lower((rai%nu+1.0)/rai%mu, lam*xmax**rai%mu)

! ub>> setzte melt_n so, dass x_r beim Schmelzvorgang erhalten bleibt:
              eva_n = MAX(MIN( (eva_q + q_r) / x_r - n_r, 0.0d0), -n_r)

            ELSE
!              eva_n = f_n * eva / x_r
              eva_n = MAX(MIN( (eva_q + q_r) / x_r - n_r, 0.0d0), -n_r)
            END IF

            eva_q = MAX(-eva_q,0.d0) 
            eva_n = MAX(-eva_n,0.d0) 

            eva_q = MIN(eva_q,q_rain(i,j,k)) 
            eva_n = MIN(eva_n,n_rain(i,j,k)) 

            rate_q(i,j,k) = eva_q

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)      = q(i,j,k)      + eva_q
            q_rain(i,j,k) = q_rain(i,j,k) - eva_q
            IF (cloud_typ > 1) THEN
              n_rain(i,j,k) = n_rain(i,j,k) - eva_n
            END IF
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,36) = eva_q
            END IF
#endif
            IF (speichere_precipstat) THEN
              evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) + eva_q
            END IF
            ! ub<<

          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(rate_n,rate_q)

  CONTAINS

    FUNCTION dmax_evap_von_dt_80(t0, ssw, temp, pres)
      IMPLICIT NONE
      DOUBLE PRECISION :: dmax_evap_von_dt_80
      DOUBLE PRECISION, INTENT(in) :: t0, ssw, temp, pres

      ! t0 Verdunstungszeit in s
      ! ssw Uebersaettigung (dimensionslos), <= 0.0!!!
      ! temp Temperatur in K
      ! pres Luftdruck in hPa

      ! dmax_evap_von_dt_fit ist max. Tropfendurchmesser, der waehrend t0 zu einem
      ! Durchmesser kleiner 80 mue-m schrumpfen kann, in mue-m

      ! Fit an exakte Loesung der Verdunstungsgleichung fuer einen Wassertropfen, UB 28.12.2006

      IF (ssw >= 0.0d0 .OR. t0 <= 0.0d0) THEN
        dmax_evap_von_dt_80 = 80.0d0
      ELSE
        IF (pres >= 350.0d0 .AND. temp >= 255.0) THEN
          dmax_evap_von_dt_80 = dmax_evap_von_dt_80_fit(t0, ssw, temp, pres)
        ELSEIF (pres < 350.0d0 .AND. temp >= 255.0) THEN
          dmax_evap_von_dt_80 = dmax_evap_von_dt_80_fit(t0, ssw, temp, 350d0)
        ELSEIF (pres >= 350.0d0 .AND. temp < 255.0) THEN
          dmax_evap_von_dt_80 = dmax_evap_von_dt_80_fit(t0, ssw, 255d0, pres)
        ELSE
          dmax_evap_von_dt_80 = dmax_evap_von_dt_80_fit(t0, ssw, 255d0, 350d0)
        END IF
      END IF

      RETURN

    END FUNCTION dmax_evap_von_dt_80

    FUNCTION dmax_evap_von_dt_80_fit(t0, ssw, temp, pres)
      IMPLICIT NONE
      DOUBLE PRECISION :: dmax_evap_von_dt_80_fit
      DOUBLE PRECISION, INTENT(in) :: t0, ssw, temp, pres

      ! t0 Verdunstungszeit in s
      ! ssw Uebersaettigung (dimensionslos), <= 0.0!!!
      ! temp Temperatur in K
      ! pres Luftdruck in hPa

      ! dmax_evap_von_dt_fit ist max. Tropfendurchmesser, der waehrend t0 zu einem
      ! Durchmesser kleiner 80 mue-m schrumpfen kann, in mue-m

      ! Fit an exakte Loesung der Verdunstungsgleichung fuer einen Wassertropfen, UB 28.12.2006

      dmax_evap_von_dt_80_fit = &
           80.0d0 + &
           (-ssw) ** (8.46314e-01 - 7.25603e-02*(-ssw)) * &
           t0 ** (8.45971e-01 - 3.05889e-05*t0) * &
           EXP( -4.57221e-01*pres + 3.38631e-04*pres**2 - &
	   7.52532e-02*temp + 5.74819e-04*temp**2 - 1.00609e-06*temp**3 + &
	   4.76251e-03*temp*pres - 1.66550e-05*temp**2*pres - &
	   3.62118e-06*temp*pres**2 + 1.29988e-08*temp**2*pres**2 + &
	   1.95564e-08*temp**3*pres + 3.00090e-11*temp*pres**3 - &
	   1.56593e-11*temp**3*pres**2 - 2.57761e-13*temp**2*pres**3 + &
	   5.34057e-16*temp**3*pres**3)

      RETURN

    END FUNCTION dmax_evap_von_dt_80_fit

    FUNCTION dmax_evap_von_dt(t0, ssw, temp, pres)
      IMPLICIT NONE
      DOUBLE PRECISION :: dmax_evap_von_dt
      DOUBLE PRECISION, INTENT(in) :: t0, ssw, temp, pres

      ! t0 Verdunstungszeit in s
      ! ssw Uebersaettigung (dimensionslos), <= 0.0!!!
      ! temp Temperatur in K
      ! pres Luftdruck in hPa

      ! dmax_evap_von_dt_fit ist max. Tropfendurchmesser, der waehrend t0 zu einem
      ! Durchmesser kleiner 3 mue-m schrumpfen kann, in mue-m

      ! Fit an exakte Loesung der Verdunstungsgleichung fuer einen Wassertropfen, UB 28.12.2006

      IF (ssw >= 0.0d0 .OR. t0 <= 0.0d0) THEN
        dmax_evap_von_dt = 2.67304d0
      ELSE
        IF (pres >= 350.0d0 .AND. temp >= 255.0) THEN
          dmax_evap_von_dt = dmax_evap_von_dt_fit(t0, ssw, temp, pres)
        ELSEIF (pres < 350.0d0 .AND. temp >= 255.0) THEN
          dmax_evap_von_dt = dmax_evap_von_dt_fit(t0, ssw, temp, 350d0)
        ELSEIF (pres >= 350.0d0 .AND. temp < 255.0) THEN
          dmax_evap_von_dt = dmax_evap_von_dt_fit(t0, ssw, 255d0, pres)
        ELSE
          dmax_evap_von_dt = dmax_evap_von_dt_fit(t0, ssw, 255d0, 350d0)
        END IF
      END IF

      RETURN

    END FUNCTION dmax_evap_von_dt

    FUNCTION dmax_evap_von_dt_fit(t0, ssw, temp, pres)
      IMPLICIT NONE
      DOUBLE PRECISION :: dmax_evap_von_dt_fit
      DOUBLE PRECISION, INTENT(in) :: t0, ssw, temp, pres
      DOUBLE PRECISION :: k(34), temp0, pres0, ssw0, t00

      ! t0 Verdunstungszeit in s
      ! ssw Uebersaettigung (dimensionslos), <= 0.0!!!
      ! temp Temperatur in K
      ! pres Luftdruck in hPa

      ! dmax_evap_von_dt_fit ist max. Tropfendurchmesser, der waehrend t0 zu einem
      ! Durchmesser kleiner 3 mue-m schrumpfen kann, in mue-m

      ! Fit an exakte Loesung der Verdunstungsgleichung fuer einen Wassertropfen, UB 28.12.2006

      k(1) =  3.0723494d+00
      k(2) =  3.4323082d+00
      k(3) = -3.1757109d+02
      k(4) =  2.2300284d+02
      k(5) = -3.1900499d+01
      k(6) =  3.2519116d+02
      k(7) = -5.9528492d+02
      k(8) = -1.5190419d+00
      k(9) = -5.6175505d+00
      k(10) =  3.2666038d+03
      k(11) = -1.1267302d+04
      k(12) = -2.3444542d+03
      k(13) =  8.2523606d+03
      k(14) =  1.3018925d+04
      k(15) =  5.2177278d+00
      k(16) = -9.7173927d+03
      k(17) = -5.3470122d+01
      k(18) =  1.2083440d+02
      k(19) =  4.4115234d+00
      k(20) =  1.9539047d-01
      k(21) = -3.6778135d-03
      k(22) =  3.0932550d-01
      k(23) = -1.6607613d+01
      k(24) =  2.7951320d+01
      k(25) = -1.9397343d+01
      k(26) =  3.2484952d+01
      k(27) =  3.6264186d+01
      k(28) =  9.9492881d+00
      k(29) = -1.5765954d+01
      k(30) = -1.9654622d+00
      k(31) =  3.0199854d+00
      k(32) = -6.0361783d+01
      k(33) = -3.1046165d+01
      k(34) =  5.7600929d+01

      temp0 = temp * 1d-3
      pres0 = pres * 1d-3
      ssw0 = -ssw
      t00 = t0 * 1d-2

      dmax_evap_von_dt_fit = 2.67304d0 + &
           (ssw0)**( k(1) + k(9)*ssw0 + k(19)*ssw**2 + &
                         k(20)*ssw0*(t00) + & 
			 k(23)*(temp0) + k(24)*(temp0)**2 + &
			 k(27)*(temp0)*(ssw0)+ &
			 k(32)*(temp0)**2*(ssw0) + &
			 k(33)*(temp0)*(ssw0)**2 + &
			 k(34)*(temp0)**2*(ssw0)**2 ) * &
           (t0)**( k(2) + k(8)*t00 + k(21)*ssw0*(t00) + &
                         k(22)*(t00)**2 + & 
                         k(25)*(temp0) + k(26)*(temp0)**2 + &
                         k(28)*(temp0)*(t00) + &
                         k(29)*(temp0)**2*(t00) + &
                         k(30)*(temp0)*(t00)**2 + &
                         k(31)*(temp0)**2*(t00)**2 ) * &
           EXP( k(3)*pres0 + k(4)*(pres0)**2 + k(5)*temp0 + &
                         k(6)*(temp0)**2 + k(7)*(temp0)**3 + &
                         k(10)*(pres0)*(temp0) + &
                         k(11)*(pres0)*(temp0)**2 + &
                         k(12)*(pres0)**2*(temp0) + &
                         k(13)*(pres0)**2*(temp0)**2 + &
                         k(14)*(pres0)*(temp0)**3 + &
                         k(15)*(pres0)**3*(temp0) + &
                         k(16)*(pres0)**2*(temp0)**3 + &
                         k(17)*(pres0)**3*(temp0)**2 + &
                         k(18)*(pres0)**3*(temp0)**3 )

      RETURN

    END FUNCTION dmax_evap_von_dt_fit

  END SUBROUTINE rain_evaporation_uli

  SUBROUTINE rain_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstungung von Regentropfen
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, q_cloud, &
         &                        q_rain, n_rain, dt, S_w, dqdt, speichere_dqdt, &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a,p_a,e_sw,q_sw,s_sw,R_f,g_d,eva,eva_q,eva_n,a_ld,a_dl
    DOUBLE PRECISION            :: q_d,q_r,n_r,x_r,d_r,v_r,f_v,N_re,x_d,e_d,f_q,f_n
    DOUBLE PRECISION            :: n,mue,d_m,gamma_eva,lam,n0r,d_vtp,gfak
    DOUBLE PRECISION, SAVE      :: c_r               !..Koeff. fuer mittlere Kapazitaet
    DOUBLE PRECISION, SAVE      :: a_n,b_n,a_q,b_q   !..Koeff. fuer mittleren Ventilationkoeff.
    TYPE(PARTICLE) :: rai

    DOUBLE PRECISION, PARAMETER :: aa = 9.65e+00     ! in SI [m/s]
    DOUBLE PRECISION, PARAMETER :: bb = 1.03e+01     ! in SI [m/s]
    DOUBLE PRECISION, PARAMETER :: cc = 6.00e+02     ! in SI [1/m]

    DOUBLE PRECISION, PARAMETER :: eps = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD rain_evaporation",4)

    rate_q = 0.0
    rate_n = 0.0

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF (firstcall.NE.1) THEN
      a_q = vent_coeff_a(rain,1)
      b_q = vent_coeff_b(rain,1)
      c_r = 1.0 / rain%cap
      firstcall = 1
      IF (isIO()) THEN
        !IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "CLOUD rain_evaporation:"
        WRITE (6,'(A,D10.3)') "     nu    = ",rain%nu
        WRITE (6,'(A,D10.3)') "     mu    = ",rain%mu
        WRITE (6,'(A,D10.3)') "     a_geo = ",rain%a_geo
        WRITE (6,'(A,D10.3)') "     b_geo = ",rain%b_geo
        WRITE (6,'(A,D10.3)') "     a_q   = ",a_q
        WRITE (6,'(A,D10.3)') "     b_q   = ",b_q
        WRITE (6,'(A,D10.3)') "     c_r   = ",c_r
        WRITE (6,'(A,D10.3)') "     cmu0  = ",rain_cmu0
        WRITE (6,'(A,D10.3)') "     cmu1  = ",rain_cmu1
        WRITE (6,'(A,D10.3)') "     cmu2  = ",rain_cmu2
        WRITE (6,'(A,D10.3)') "     cmu3  = ",rain_cmu3
        WRITE (6,'(A,D10.3)') "     cmu4  = ",rain_cmu4
        WRITE (6,'(A,D10.3)') "     g_fak = ",rain_gfak
      END IF
    ELSEIF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD rain_evaporation " 
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_r = q_rain(i,j,k)
          n_r  = n_rain(i,j,k)
          q_d  = q(i,j,k)
          p_a  = p_0(i,j,k)
          T_a  = T_0(i,j,k)
          x_d  = q(i,j,k) / rho_0(i,j,k)
          e_d  = q(i,j,k) * R_d * T_a
          e_sw = e_ws(T_a)
          s_sw = e_d / e_sw - 1.0  !..Uebersaettigung bzgl. Wasser

          IF(s_sw < 0.0 .AND. q_r > 0.0d0 .AND. q_cloud(i,j,k) < q_krit)THEN

            D_vtp = 8.7602e-5 * T_a**(1.81) / p_a
            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )

            IF (use_mu_Dm_rain_evap) THEN
              x_r = q_r / (n_r+eps)
              x_r = MIN(MAX(x_r,rain%x_min),rain%x_max)

              D_m = rain%a_geo * x_r**rain%b_geo                 

              IF (D_m.LE.rain_cmu3) THEN ! see Seifert (2008)            
                mue = rain_cmu0*TANH((4.*rain_cmu2*(D_m-rain_cmu3))**rain_cmu5) &
                  & + rain_cmu4
              ELSE
                mue = rain_cmu1*TANH((1.*rain_cmu2*(D_m-rain_cmu3))**rain_cmu5) &
                  & + rain_cmu4
              ENDIF

              lam = (pi/6.*rho_w &
                &      * (mue+3.0)*(mue+2.0)*(mue+1.0) / x_r)**(1./3.)

!AS20081207 optimierte Variante

              ! chebyshev approximation of gamma_fct(zrmue+5/2)/gamma_fct(zrmue+2)
              ! (for mue in [0,31], error smaller than 2%)
!             gfak =  0.1357940435E+01              &
!                  &   + mue * 0.3033273220E+00 * ( &
!                  &   - mue * 0.1299313363E-01 * ( &
!                  &   + mue * 0.4002257774E-03 * ( &
!                  &   - mue * 0.4856703981E-05 ) ) )

! Uli              
              gfak =  0.1357940435E+01 &
                &  + mue * ( +0.3033273220E+00  &
                &  + mue * ( -0.1299313363E-01  &
                &  + mue * ( +0.4002257774E-03  &
                &  - mue * 0.4856703981E-05 ) ) )

!              gfak = gfct(mue+5/2)/gfct(mue+2)

!               gfak = 0.1357940435E+01          &
!                  & + 0.3033273220E+00 * mue    &
!                  & - 0.1299313363E-01 * mue**2 &
!                  & + 0.4002257774E-03 * mue**3 & 
!                  & - 0.4856703981E-05 * mue**4

              f_q  = rain%a_ven                                          &
              &    + rain%b_ven * N_sc**n_f                              &
              &                 * (aa/nu_l*rrho_04(i,j,k))**m_f          &
              &    * gfak / SQRT(lam)                                    &
              &    * (1.0                                                &
              &      - 1./2.  * (bb/aa)**1 * (lam/(1.*cc+lam))**(mue+5.0/2.0) &
              &      - 1./8.  * (bb/aa)**2 * (lam/(2.*cc+lam))**(mue+5.0/2.0) &
              &      - 1./16. * (bb/aa)**3 * (lam/(3.*cc+lam))**(mue+5.0/2.0) &
              &      - 5./127.* (bb/aa)**4 * (lam/(4.*cc+lam))**(mue+5.0/2.0) &
              &      )

              
              IF (rain_gfak.GT.0) THEN
                gamma_eva = rain_gfak * (1.1e-3/D_m) * EXP(-0.2*mue)
              ELSE
                gamma_eva = 1.0
              END IF

              eva_q = g_d * c_r * n_r * (mue+1.0) / lam * f_q * s_sw * dt
              eva_n = gamma_eva * eva_q / x_r

            ELSE ! old evaporation

              IF (cloud_typ <= 1) THEN
                n0r = 1.00d+07
                n_r = n0r * (pi*rho_w*n0r/q_r)**(-1./4.)     !..Nur q_r  (n0 = const) 
              ENDIF
              x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
              x_r = MAX(1e3/2e7*q_r,MIN(x_r,1e4/2e5*q_r))

              d_r = rain%a_geo * x_r**rain%b_geo                 
              v_r = rain%a_vel * x_r**rain%b_vel * rrho_04(i,j,k)
              
              N_re = v_r * d_r / nu_l                            
              f_v  = N_sc**n_f * N_re**m_f           
              
              f_q  = a_q + b_q * f_v
              
              eva   = g_d * n_r * c_r * d_r * s_sw * dt 
              eva_q = f_q * eva
              !eva_n = f_q * eva / x_r
              !eva_n = MAX(MIN( (eva_q + q_r) / x_r - n_r, 0.0d0), -n_r)
              eva_n = MAX(MIN( eva_q / x_r, 0.0d0), -n_r)
            ENDIF

            eva_q = MAX(-eva_q,0.d0) 
            eva_n = MAX(-eva_n,0.d0) 

            eva_q = MIN(eva_q,q_r) 
            eva_n = MIN(eva_n,n_r) 

            rate_q(i,j,k) = eva_q

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)      = q(i,j,k)      + eva_q
            q_rain(i,j,k) = q_rain(i,j,k) - eva_q
            IF (cloud_typ > 1) THEN
              n_rain(i,j,k) = n_rain(i,j,k) - eva_n
            END IF
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,36) = eva_q
            END IF
#endif
            IF (speichere_precipstat) THEN
              evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) + eva_q
            END IF
            ! ub<<

          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE rain_evaporation

  SUBROUTINE graupel_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstung von schmelzendem Graupel                        *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_graupel, n_graupel, dt, dqdt, speichere_dqdt, &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a,p_a,e_sw,q_sw,s_sw,R_f,g_d,eva,a_ld,a_dl
    DOUBLE PRECISION            :: q_d,q_g,n_g,x_g,d_g,v_g,f_v,N_re,x_d,e_d
    DOUBLE PRECISION, SAVE      :: c_g             !..Koeff. fuer mittlere Kapazitaet
    DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD graupel_evaporation",4)

    rate_q = 0.0
    rate_n = 0.0

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(graupel,1)
      b_f = vent_coeff_b(graupel,1)
      c_g = 1.0 / graupel%cap
      firstcall = 1
      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "CLOUD graupel_evaporation:"
        WRITE (6,'(A,D10.3)') "     a_f = ",a_f
        WRITE (6,'(A,D10.3)') "     b_f = ",b_f
        WRITE (6,'(A,D10.3)') "     c_g = ",c_g
      END IF
    ELSEIF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD graupel_evaporation " 
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_g = q_graupel(i,j,k)
          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          IF(q_g > 0.d0 .AND. T_a > T_3)THEN
            q_d = q(i,j,k)
            p_a  = p(i,j,k) + p_0(i,j,k)
            x_d  = q(i,j,k) / (rho_0(i,j,k)+rho_g(i,j,k))
            e_d  = q(i,j,k) * R_d * T_a 
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser

            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_g = n_graupel(i,j,k)                          !..Anzahldichte in SI
            q_g = q_graupel(i,j,k)                          !..Fluessigwassergehalt in SI

            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max)   !..mittlere Masse in SI

            d_g = graupel%a_geo * x_g**graupel%b_geo                   !..mittlerer Durchmesser
            v_g = graupel%a_vel * x_g**graupel%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_g * d_g / nu_l                         !..mittlere Reynoldszahl
            f_v  = a_f + b_f * N_sc**n_f * N_re**m_f        !..mittlerer Vent.Koeff.

            eva = g_d * n_g * c_g * d_g * f_v * s_sw * dt

            eva = MAX(-eva,0.d0) 
            eva = MIN(eva,q_g) 

            rate_q(i,j,k) = eva

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)         = q(i,j,k)         + eva
            q_graupel(i,j,k) = q_graupel(i,j,k) - eva
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,33) = eva
            END IF
#endif
            IF (speichere_precipstat) THEN
              evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) + eva
            END IF
            ! ub<<

            !..Berechnung der Verdampfungswaerme
            !WRF!T(i,j,k) = T(i,j,k) - L_wd / cp * eva / (rho_0(i,j,k)+rho_g(i,j,k))
            !WRF!p(i,j,k) = p(i,j,k) - L_wd / cv * eva * R_l
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE graupel_evaporation

  SUBROUTINE hail_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstung von schmelzendem Hagel                          *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_hail, n_hail, dt, dqdt, speichere_dqdt, &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a,p_a,e_sw,q_sw,s_sw,R_f,g_d,eva,a_ld,a_dl
    DOUBLE PRECISION            :: q_d,q_h,n_h,x_h,d_h,v_h,f_v,N_re,x_d,e_d
    DOUBLE PRECISION, SAVE      :: c_h             !..Koeff. fuer mittlere Kapazitaet
    DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD hail_evaporation",4)

    rate_q = 0.0
    rate_n = 0.0

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(hail,1)
      b_f = vent_coeff_b(hail,1)
      c_h = 1.0 / hail%cap
      firstcall = 1
      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "CLOUD hail_evaporation:"
        WRITE (6,'(A,D10.3)') "     a_f = ",a_f
        WRITE (6,'(A,D10.3)') "     b_f = ",b_f
        WRITE (6,'(A,D10.3)') "     c_h = ",c_h
      END IF
    ELSEIF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD hail_evaporation " 
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_h = q_hail(i,j,k)
          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          IF(q_h > 0.d0 .AND. T_a > T_3)THEN
            q_d = q(i,j,k)
            p_a  = p(i,j,k) + p_0(i,j,k)
            x_d  = q(i,j,k) / (rho_0(i,j,k)+rho_g(i,j,k))
            e_d  = q(i,j,k) * R_d * T_a 
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser

            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_h = n_hail(i,j,k)                          !..Anzahldichte in SI
            q_h = q_hail(i,j,k)                          !..Fluessigwassergehalt in SI

            x_h = MIN(MAX(q_h/(n_h+eps),hail%x_min),hail%x_max)   !..mittlere Masse in SI

            d_h = hail%a_geo * x_h**hail%b_geo                   !..mittlerer Durchmesser
            v_h = hail%a_vel * x_h**hail%b_vel * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

            N_re = v_h * d_h / nu_l                         !..mittlere Reynoldszahl
            f_v  = a_f + b_f * N_sc**n_f * N_re**m_f        !..mittlerer Vent.Koeff.

            eva = g_d * n_h * c_h * d_h * f_v * s_sw * dt

            eva = MAX(-eva,0.d0) 
            eva = MIN(eva,q_h) 

            rate_q(i,j,k) = eva

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)      = q(i,j,k)      + eva
            q_hail(i,j,k) = q_hail(i,j,k) - eva
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,53) = eva
            END IF
#endif
            IF (speichere_precipstat) THEN
              evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) + eva
            END IF
            ! ub<<

            !..Berechnung der Verdampfungswaerme
            !WRF!T(i,j,k) = T(i,j,k) - L_wd / cp * eva / (rho_0(i,j,k)+rho_g(i,j,k))
            !WRF!p(i,j,k) = p(i,j,k) - L_wd / cv * eva * R_l
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE hail_evaporation

  SUBROUTINE snow_evaporation ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Verdunstung von schmelzendem Schnee                         *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g, &
         &                        q_snow, n_snow, dt, dqdt, speichere_dqdt, &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    ! .. Local Variables ..
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a,p_a,e_sw,q_sw,s_sw,x_d,R_f,g_d,eva,a_ld,a_dl
    DOUBLE PRECISION            :: q_d,q_s,n_s,x_s,d_s,v_s,f_v,N_re,e_d
    DOUBLE PRECISION, SAVE      :: c_s             !..Koeff. fuer mittlere Kapazitaet
    DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.

    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD snow_evaporation",4)

    rate_q = 0.0
    rate_n = 0.0

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF (firstcall.NE.1) THEN
      a_f = vent_coeff_a(snow,1)
      b_f = vent_coeff_b(snow,1)
      c_s = 1.0/snow%cap
      firstcall = 1
      IF (isIO() .AND. isdebug) THEN
        WRITE (6,'(A,D10.3)') "CLOUD snow_evaporation:"
        WRITE (6,'(A,D10.3)') "     a_f = ",a_f
        WRITE (6,'(A,D10.3)') "     b_f = ",b_f
        WRITE (6,'(A,D10.3)') "     c_s = ",c_s
      END IF
    ELSEIF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD snow_evaporation " 
    END IF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_s = q_snow(i,j,k)
          T_a  = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)
          IF(q_s > 0.d0 .AND. T_a > T_3)THEN
            q_d = q(i,j,k)
            p_a  = p(i,j,k) + p_0(i,j,k)
            x_d  = q(i,j,k) / (rho_0(i,j,k)+rho_g(i,j,k))
            e_d  = q(i,j,k) * R_d * T_a 
            e_sw = e_ws(T_a)
            s_sw = e_d / e_sw - 1.0                   !..Uebersaettigung bzgl. Wasser


            g_d = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_3**2) + R_d * T_3 / (D_v * e_sw) )

            n_s = n_snow(i,j,k)                             !..Anzahldichte in SI
            q_s = q_snow(i,j,k)                             !..Fluessigwassergehalt in SI

            x_s = MIN(MAX(q_s/(n_s+eps),snow%x_min),snow%x_max)    !..mittlere Masse in SI     

            d_s = snow%a_geo * x_s**snow%b_geo                     !..mittlerer Durchmesser
            v_s = snow%a_vel * x_s**snow%b_vel * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

            N_re = v_s * d_s / nu_l                          !..mittlere Reynoldszahl
            f_v  = a_f + b_f * N_sc**n_f * N_re**m_f   !..mittlerer Vent.Koeff.

            eva = g_d * n_s * c_s * d_s * f_v * s_sw * dt

            eva = MAX(-eva,0.d0) 
            eva = MIN(eva,q_s) 

            rate_q(i,j,k) = eva

            !..Berechnung der H_2O-Komponenten
            q(i,j,k)         = q(i,j,k)         + eva
            q_snow(i,j,k) = q_snow(i,j,k) - eva
            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              dqdt(i,j,k,32) = eva
            END IF
#endif
            IF (speichere_precipstat) THEN
              evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) + eva
            END IF
            ! ub<<

          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE snow_evaporation

  !===========================================================================
  ! Subroutinen fuer die Wet Growth Parametrisierung:
  ! Initialisierung: Einlesen der Lookup-table aus einer Textdatei.
  ! Diese Subroutine muss von der Interface-Routine des 2-M-Schemas 
  ! aufgerufen werden.
  ! Eventuelle Verteilung der Table auf alle Knoten bei Parallelbetrieb
  ! muss ebenfalls von der Interface-Routine besorgt werden.

  SUBROUTINE init_dmin_wetgrowth(dateiname, unitnr)

    USE globale_variablen, ONLY : dmin_wg_g, pvec_wg_g, Tvec_wg_g, &
         & qwvec_wg_g, qivec_wg_g, anzp_wg, anzT_wg, anzi_wg, anzw_wg

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: dateiname
    INTEGER, INTENT(in) :: unitnr
    INTEGER :: i, j, k, l, error

    OPEN(unitnr, file=TRIM(dateiname), status='old', form='formatted', iostat=error)
    IF (error /= 0) THEN
      WRITE (*,*) 'init_dmin_wetgrowth: lookup-table ' &
        // TRIM(dateiname) // ' not found'
      STOP
    END IF

    READ (unitnr,*) anzp_wg, anzT_wg, anzi_wg, anzw_wg

    ALLOCATE(pvec_wg_g(anzp_wg))
    ALLOCATE(Tvec_wg_g(anzT_wg))
    ALLOCATE(qwvec_wg_g(anzw_wg))
    ALLOCATE(qivec_wg_g(anzi_wg))
    ALLOCATE(dmin_wg_g(anzp_wg,anzT_wg,anzw_wg,anzi_wg))

    READ (unitnr,*,iostat=error) pvec_wg_g(1:anzp_wg)
    IF (error /= 0) THEN
      WRITE (*,*) 'init_dmin_wetgrowth: Error reading pvec from ' &
        // TRIM(dateiname)
      STOP
    END IF
    READ (unitnr,*,iostat=error) Tvec_wg_g(1:anzT_wg)
    IF (error /= 0) THEN
      WRITE (*,*) 'init_dmin_wetgrowth: Error reading Tvec from ' &
        // TRIM(dateiname)
      STOP
    END IF
    READ (unitnr,*,iostat=error) qwvec_wg_g(1:anzw_wg)
    IF (error /= 0) THEN
      WRITE (*,*) 'init_dmin_wetgrowth: Error reading qwvec from ' &
        // TRIM(dateiname)
      STOP
    END IF
    READ (unitnr,*,iostat=error) qivec_wg_g(1:anzi_wg)
    IF (error /= 0) THEN
      WRITE (*,*) 'init_dmin_wetgrowth: Error reading qivec from ' &
        // TRIM(dateiname)
      STOP
    END IF

    READ (unitnr,*,iostat=error) dmin_wg_g(1:anzp_wg,1:anzT_wg,1:anzw_wg,1:anzi_wg)
    IF (error /= 0) THEN
      WRITE (*,*) 'init_dmin_wetgrowth: Error reading dmin from ' &
        // TRIM(dateiname)
    END IF
    
    CLOSE(unitnr)

    RETURN

  END SUBROUTINE init_dmin_wetgrowth


  ! wet growth Grenzdurchmesser fuer graupelhail2test in m:
  FUNCTION dmin_wetgrowth_graupel(p_a,T_a,qw_a,qi_a)

    USE globale_variablen, ONLY : dmin_wg_g, pvec_wg_g, Tvec_wg_g, &
         & qwvec_wg_g, qivec_wg_g, anzp_wg, anzT_wg, anzi_wg, anzw_wg

    IMPLICIT NONE

    DOUBLE PRECISION :: dmin_wetgrowth_graupel
    DOUBLE PRECISION, INTENT(in) :: p_a,T_a,qw_a,qi_a
    DOUBLE PRECISION :: p_lok,T_lok,qw_lok,qi_lok
    
    INTEGER :: i
    INTEGER :: iu, io, ju, jo, ku, ko, lu, lo

    DOUBLE PRECISION :: hilf1(2,2,2,2), hilf2(2,2,2), hilf3(2,2), hilf4(2)

    LOGICAL :: found_p, found_T, found_w, found_i

    found_p = .FALSE.
    found_T = .FALSE.
    found_w = .FALSE.
    found_i = .FALSE.
    dmin_wetgrowth_graupel = 999.99

    p_lok = MIN(MAX(p_a,pvec_wg_g(1)),pvec_wg_g(anzp_wg))
    IF (p_a <= pvec_wg_g(1)) THEN
      found_p = .TRUE.
      iu = 1
      io = 2
    ELSE IF (p_a >= pvec_wg_g(anzp_wg)) THEN
      found_p = .TRUE.
      iu = anzp_wg - 1
      io = anzp_wg
    ELSE
      iu = 1
      DO i=1, anzp_wg-1
        IF (p_a >= pvec_wg_g(i) .AND. p_a < pvec_wg_g(i+1)) THEN
          iu = i
          found_p = .TRUE.
          EXIT
        END IF
      END DO
      io = iu + 1
    END IF

    T_lok = MIN(MAX(T_a,Tvec_wg_g(1)),Tvec_wg_g(anzT_wg))
    IF (T_a <= Tvec_wg_g(1)) THEN
      found_T = .TRUE.
      ju = 1
      jo = 2
    ELSE IF (T_a >= Tvec_wg_g(anzT_wg)) THEN
      found_T = .TRUE.
      dmin_wetgrowth_graupel = 0.0
      RETURN
    ELSE
      ju = 1
      DO i=1, anzT_wg-1
        IF (T_a >= Tvec_wg_g(i) .AND. T_a < Tvec_wg_g(i+1)) THEN
          ju = i
          found_T = .TRUE.
          EXIT
        END IF
      END DO
      jo = ju + 1
    END IF

    qw_lok = MIN(MAX(qw_a,qwvec_wg_g(1)),qwvec_wg_g(anzw_wg))
    IF (qw_a <= qwvec_wg_g(1)) THEN
      found_w = .TRUE.
      dmin_wetgrowth_graupel = 999.99
      RETURN
    ELSE IF (qw_a >= qwvec_wg_g(anzw_wg)) THEN
      found_w = .TRUE.
      ku = anzw_wg - 1
      ko = anzw_wg
    ELSE
      ku = 1
      DO i=1, anzw_wg-1
        IF (qw_a >= qwvec_wg_g(i) .AND. qw_a < qwvec_wg_g(i+1)) THEN
          ku = i
          found_w = .TRUE.
          EXIT
        END IF
      END DO
      ko = ku + 1
    END IF

    qi_lok = MIN(MAX(qi_a,qivec_wg_g(1)),qivec_wg_g(anzi_wg))
    IF (qi_a <= qivec_wg_g(1)) THEN
      found_i = .TRUE.
      lu = 1
      lo = 2
    ELSE IF (qi_a >= qivec_wg_g(anzi_wg)) THEN
      found_i = .TRUE.
      lu = anzi_wg - 1
      lo = anzi_wg
    ELSE
      lu = 1
      DO i=1, anzi_wg-1
        IF (qi_a >= qivec_wg_g(i) .AND. qi_a < qivec_wg_g(i+1)) THEN
          lu = i
          found_i = .TRUE.
          EXIT
        END IF
      END DO
      lo = lu + 1
    END IF

    IF (.NOT.found_p .OR. .NOT.found_T .OR. .NOT.found_w .OR. .NOT. found_i) THEN
      WRITE (*,*) 'Seifert dmin_wetgrowth_graupel: interpolation point not found in lookup table'
      dmin_wetgrowth_graupel = 999.99
    ELSE

      ! Tetra-lineare Interpolation von Dmin:
      hilf1 = dmin_wg_g(iu:io,ju:jo,ku:ko,lu:lo)
      hilf2 = hilf1(1,:,:,:) + &
           (hilf1(2,:,:,:)-hilf1(1,:,:,:)) / &
           (pvec_wg_g(io)-pvec_wg_g(iu)) * (p_lok-pvec_wg_g(iu))
      hilf3 = hilf2(1,:,:) + &
           (hilf2(2,:,:)-hilf2(1,:,:)) / (Tvec_wg_g(jo)-Tvec_wg_g(ju)) * (T_lok-Tvec_wg_g(ju))
      
      hilf4 = hilf3(1,:) + &
           (hilf3(2,:)-hilf3(1,:)) / (qwvec_wg_g(ko)-qwvec_wg_g(ku)) * (qw_lok-qwvec_wg_g(ku))
      
      dmin_wetgrowth_graupel = hilf4(1) + &
           (hilf4(2)-hilf4(1)) / (qivec_wg_g(lo)-qivec_wg_g(lu)) * (qi_lok-qivec_wg_g(lu))
      
    END IF

    RETURN

  END FUNCTION dmin_wetgrowth_graupel

  !===========================================================================
  !
  ! Subroutinen fuer die Wet Growth Parametrisierung mit aequidistantem table lookup
  ! fuer eine bessere Vektorisierung.
  !
  ! - Initialisierung: Einlesen der Lookup-table aus einer Textdatei.
  !   Diese Subroutine muss von der Interface-Routine des 2-M-Schemas 
  !   aufgerufen werden.
  !   Eventuelle Verteilung der Table auf alle Knoten bei Parallelbetrieb
  !   muss ebenfalls von der Interface-Routine besorgt werden.
  !
  ! - Die eingelesene lookup table muss bereits in p, qw und qi aequidistant sein. Nur bzgl. T kann eine 
  !   nicht-aequidistanter grid-Vektor vorliegen.
  !
  ! Fuer die eigentliche Table wird das Struct lookupt_4D verwendet:
  !
  !===========================================================================

  SUBROUTINE init_dmin_wg_gr_ltab_equi(dateiname, unitnr, ndT, ltab)

    USE wolken_konstanten, ONLY : lookupt_4D

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: dateiname
    INTEGER, INTENT(in) :: unitnr
    ! Desired number of elements for the fine equidistant grid vector for T:
    INTEGER, INTENT(in) :: ndT
    TYPE(lookupt_4D), INTENT(out) :: ltab

    ! grid spacings of the desired fine grid vectors:
    DOUBLE PRECISION :: minT, maxT
    INTEGER :: i, j, k, l, error, ii

    DOUBLE PRECISION, ALLOCATABLE :: Tvec_wg_g_loc(:), dmin_wg_g_loc(:,:,:,:)
    INTEGER :: anzT_wg_loc
    INTEGER :: ju, jo


    !!! 1) Read the original lookup table from a file. This table may be made of a nonequidistant grid vector for T.
    !!!    The grid vectors for p, qw and qi have to be equidistant.

    OPEN(unitnr, file=TRIM(dateiname), status='old', form='formatted', iostat=error)
    IF (error /= 0) THEN
      WRITE (*,*) 'dmin_wg_gr_ltab_equi: lookup-table ' &
           // TRIM(dateiname) // ' not found'
      STOP
    END IF

    READ (unitnr,*) ltab%n1, anzT_wg_loc, ltab%n3, ltab%n4

    ALLOCATE( Tvec_wg_g_loc(anzT_wg_loc) )
    ALLOCATE( ltab%x1(ltab%n1) )
    ALLOCATE( ltab%x3(ltab%n3) )
    ALLOCATE( ltab%x4(ltab%n4) )
    ALLOCATE( dmin_wg_g_loc(ltab%n1,anzT_wg_loc,ltab%n3,ltab%n4) )

    READ (unitnr,*,iostat=error) ltab%x1(1:ltab%n1)
    IF (error /= 0) THEN
      WRITE (*,*) 'dmin_wg_gr_ltab_equi: Error reading pvec from ' &
           // TRIM(dateiname)
      STOP
    END IF
    READ (unitnr,*,iostat=error) Tvec_wg_g_loc(1:anzT_wg_loc)
    IF (error /= 0) THEN
      WRITE (*,*) 'min_wg_gr_ltab_equi: Error reading Tvec from ' &
           // TRIM(dateiname)
      STOP
    END IF
    READ (unitnr,*,iostat=error) ltab%x3(1:ltab%n3)
    IF (error /= 0) THEN
      WRITE (*,*) 'dmin_wg_gr_ltab_equi: Error reading qwvec from ' &
           // TRIM(dateiname)
      STOP
    END IF
    READ (unitnr,*,iostat=error) ltab%x4(1:ltab%n4)
    IF (error /= 0) THEN
      WRITE (*,*) 'dmin_wg_gr_ltab_equi: Error reading qivec from ' &
           // TRIM(dateiname)
      STOP
    END IF

    DO l=1, ltab%n4
      DO k=1, ltab%n3
        DO j=1, anzT_wg_loc
          DO i=1,ltab%n1
            READ (unitnr,*,iostat=error) dmin_wg_g_loc(i,j,k,l)
            IF (error /= 0) THEN
              WRITE (*,*) l,k,j,i
              WRITE (*,*) 'dmin_wg_gr_ltab_equi: Error reading dmin from ' &
                   // TRIM(dateiname)
              STOP
            END IF
          END DO
        END DO
      END DO
    END DO

    CLOSE(unitnr)


    !!! 2) Generate equidistant table vectors and construct the 
    !!!    equidistant Dmin-lookuptable by linear oversampling:
    ltab%n2 = ndT

    ALLOCATE( ltab%x2(ltab%n2) )
    ALLOCATE( ltab%ltable(ltab%n1,ltab%n2,ltab%n3,ltab%n4) )

    minT  = Tvec_wg_g_loc (1)
    maxT  = Tvec_wg_g_loc (anzT_wg_loc)
    
    ltab%dx1      = ltab%x1(2) - ltab%x1(1)
    ltab%odx1     = 1.0d0 / ltab%dx1
    ltab%dx2      = (maxT - minT) / (ndT - 1.0d0)
    ltab%odx2     = 1.0d0 / ltab%dx2
    ltab%dx3      = ltab%x3(2) - ltab%x3(1)
    ltab%odx3     = 1.0d0 / ltab%dx3
    ltab%dx4      = ltab%x4(2) - ltab%x4(1)
    ltab%odx4     = 1.0d0 / ltab%dx4

    ! Equidistant grid vectors for T:
    DO j=1, ltab%n2
      ltab%x2(j) = minT + (j-1) * ltab%dx2
    END DO

    ! Linear interpolation w.r.t. T of the equidistant Dmin-lookuptable from
    ! the original table in the datafile, which may be non-equidistant
    ! w.r.t. T:
    DO j=1, ltab%n2

      ju = 1
      DO ii=1, anzT_wg_loc-1
        IF (ltab%x2(j) >= Tvec_wg_g_loc(ii) .AND. ltab%x2(j) <= Tvec_wg_g_loc(ii+1)) THEN
          ju = ii
          EXIT
        END IF
      END DO
      jo = ju + 1

      ! Linear interplation of Dmin with respect to T:
      ltab%ltable(:,j,:,:) = dmin_wg_g_loc(:,ju,:,:) + &
           (dmin_wg_g_loc(:,jo,:,:) - dmin_wg_g_loc(:,ju,:,:)) / &
           (Tvec_wg_g_loc(jo)-Tvec_wg_g_loc(ju))  * (ltab%x2(j)-Tvec_wg_g_loc(ju))

    END DO

    ! clean up memory:
    DEALLOCATE(Tvec_wg_g_loc,dmin_wg_g_loc)

    RETURN

  END SUBROUTINE init_dmin_wg_gr_ltab_equi

  ! wet growth Grenzdurchmesser fuer graupelhail2test3 in m:
  FUNCTION dmin_wg_gr_ltab_equi(p_a,T_a,qw_a,qi_a,ltab) RESULT (dmin_loc)

    USE wolken_konstanten, ONLY : lookupt_4D

    IMPLICIT NONE

    DOUBLE PRECISION :: dmin_loc
    DOUBLE PRECISION, INTENT(in) :: p_a,T_a,qw_a,qi_a
    TYPE(lookupt_4D), INTENT(in) :: ltab
    DOUBLE PRECISION :: p_lok,T_lok,qw_lok,qi_lok

    INTEGER :: iu, io, ju, jo, ku, ko, lu, lo
    DOUBLE PRECISION :: hilf1(2,2,2,2), hilf2(2,2,2), hilf3(2,2), hilf4(2)

    IF (T_a >= ltab%x2(ltab%n2)) THEN

      dmin_loc = 0.0d0

    ELSE IF (T_a < ltab%x2(1)) THEN

      dmin_loc = 999.99d0

    ELSE

      p_lok = MIN(MAX(p_a,ltab%x1(1)),ltab%x1(ltab%n1))
      iu = MIN(FLOOR((p_lok - ltab%x1(1)) * ltab%odx1 ) + 1, ltab%n1-1)
      io = iu + 1
      T_lok = MIN(MAX(T_a,ltab%x2(1)),ltab%x2(ltab%n2))
      ju = MIN(FLOOR((T_lok - ltab%x2(1)) * ltab%odx2 ) + 1, ltab%n2-1)
      jo = ju + 1
      qw_lok = MIN(MAX(qw_a,ltab%x3(1)),ltab%x3(ltab%n3))
      ku = MIN(FLOOR((qw_lok - ltab%x3(1)) * ltab%odx3 ) + 1, ltab%n3-1)
      ko = ku + 1
      qi_lok = MIN(MAX(qi_a,ltab%x4(1)),ltab%x4(ltab%n4))
      lu = MIN(FLOOR((qi_lok - ltab%x4(1)) * ltab%odx4 ) + 1, ltab%n4-1)
      lo = lu + 1


      ! Tetra-lineare Interpolation von Dmin:
      hilf1 = ltab%ltable(iu:io,ju:jo,ku:ko,lu:lo)
      hilf2 = hilf1(1,:,:,:) + (hilf1(2,:,:,:) - hilf1(1,:,:,:)) * ltab%odx1 * (p_lok-ltab%x1(iu) )
      hilf3 = hilf2(1,:,:)   + (hilf2(2,:,:)   - hilf2(1,:,:)  ) * ltab%odx2 * (T_lok-ltab%x2(ju) )
      hilf4 = hilf3(1,:)     + (hilf3(2,:)     - hilf3(1,:)    ) * ltab%odx3 * (qw_lok-ltab%x3(ku))


      dmin_loc = hilf4(1) + (hilf4(2) - hilf4(1))  * ltab%odx4 * (qi_lok-ltab%x4(lu))


    END IF

    RETURN

  END FUNCTION dmin_wg_gr_ltab_equi

  SUBROUTINE graupel_hail_conv_wet()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Wet Growth Umwandlungsprozesses von Graupel zu Hagel    *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g,        &
         &                        q_graupel, n_graupel, q_cloud, n_cloud, q_rain, n_rain,  &
         &                        q_ice, n_ice, q_snow, n_snow, q_hail, n_hail, dt, dqdt,  &
         &                        speichere_dqdt, ltabdminwgg
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0
    USE gamma_functions_mp_seifert,    ONLY: incgfct_lower

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a, p_a, d_trenn, qw_a, qi_a, N_0, lam, xmin
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g
    DOUBLE PRECISION            :: q_c,n_c,x_c,d_c
    DOUBLE PRECISION            :: q_r,n_r,x_r,d_r
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s
    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h
    DOUBLE PRECISION            :: conv_n, conv_q

    DOUBLE PRECISION, PARAMETER :: eps = 1.d-20
    DOUBLE PRECISION, PARAMETER :: eps2 = 1.d-12  ! UB_20081202


    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix             
          q_c = q_cloud(i,j,k)                                      !..Massendichte
          q_r = q_rain(i,j,k)                                       !..Massendichte
          q_g = q_graupel(i,j,k)                                    !..Massendichte
          n_g = n_graupel(i,j,k)                                    !..Anzahldichte
          x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse
          D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
          n_g = q_g / x_g

          T_a = T_0(i,j,k)
          p_a = p_0(i,j,k)

          !.. Umgebungsgehalt unterkuehltes Wasser:
          qw_a = q_r + q_c

          IF (T_a < T_3 .AND. q_g > q_krit_gc .AND. qw_a > 1e-3) THEN

            !.. Umgebungsgehalt Eispartikel (vernachl. werden Graupel und Hagel wg. geringer Kollisionseff.)
            !.. koennte problematisch sein, weil in konvekt. Wolken viel mehr Graupel und Hagel enthalten ist!!!
            qi_a = q_ice(i,j,k) + q_snow(i,j,k)
            d_trenn = dmin_wetgrowth_graupel(p_a,T_a,qw_a,qi_a)

            !.. Bereich im Graupelspektrum mit D > d_trenn wird zu Hagel:
            xmin = (d_trenn/graupel%a_geo)**(1.0d0/graupel%b_geo)

! UB_20081202>>            IF (xmin > graupel%x_min .AND. d_trenn < 10.0d0 * D_g) THEN
            IF (xmin > eps2 .AND. d_trenn < 10.0d0 * D_g) THEN

              lam = ( gfct((graupel%nu+1.0)/graupel%mu) / gfct((graupel%nu+2.0)/graupel%mu) * x_g)**(-graupel%mu)
              n_0 = graupel%mu * n_g * lam**((graupel%nu+1.0)/graupel%mu) / gfct((graupel%nu+1.0)/graupel%mu)
              conv_n = n_0/(graupel%mu*lam**((graupel%nu+1.0)/graupel%mu))* &
                   incgfct_upper((graupel%nu+1.0)/graupel%mu, lam*xmin**graupel%mu)
              conv_q = n_0/(graupel%mu*lam**((graupel%nu+2.0)/graupel%mu))* &
                   incgfct_upper((graupel%nu+2.0)/graupel%mu, lam*xmin**graupel%mu)
              
              conv_n = MIN(conv_n,n_g)
              conv_q = MIN(conv_q,q_g)

              q_graupel(i,j,k) = q_graupel(i,j,k) - conv_q
              n_graupel(i,j,k) = n_g - conv_n

              q_hail(i,j,k) = q_hail(i,j,k) + conv_q
              n_hail(i,j,k) = n_hail(i,j,k) + conv_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,59) = conv_q
              END IF
#endif

            END IF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_hail_conv_wet

  SUBROUTINE graupel_hail_conv_wet_gamlook()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Wet Growth Umwandlungsprozesses von Graupel zu Hagel    *
    !                                                                              *
    !       Berechnung der incomplete gamma functions mit table lookup
    !
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, p_g, t_g, rho_g,        &
         &                        q_graupel, n_graupel, q_cloud, n_cloud, q_rain, n_rain,  &
         &                        q_ice, n_ice, q_snow, n_snow, q_hail, n_hail, dt, dqdt,  &
         &                        speichere_dqdt, ltabdminwgg
    USE konstanten,         ONLY: pi,cv,cp,wolke_typ
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0
    USE gamma_functions_mp_seifert, ONLY: incgfct_lower_lookupcreate, incgfct_lower_lookup, &
         &                        nlookup, nlookuphr_dummy

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER                     :: i,j,k
    INTEGER                     :: stat_var = 0
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: T_a, p_a, d_trenn, qw_a, qi_a, N_0, lam, xmin
    DOUBLE PRECISION            :: q_g,n_g,x_g,d_g
    DOUBLE PRECISION            :: q_c,n_c,x_c,d_c
    DOUBLE PRECISION            :: q_r,n_r,x_r,d_r
    DOUBLE PRECISION            :: q_i,n_i,x_i,d_i
    DOUBLE PRECISION            :: q_s,n_s,x_s,d_s
    DOUBLE PRECISION            :: q_h,n_h,x_h,d_h
    DOUBLE PRECISION            :: conv_n, conv_q

    TYPE(gamlookuptable), SAVE  :: ltable1, ltable2
    DOUBLE PRECISION, SAVE      :: nm1, nm2, g1, g2

    DOUBLE PRECISION, PARAMETER :: eps = 1.d-20
    DOUBLE PRECISION, PARAMETER :: eps2 = 1.d-12

    IF (firstcall.NE.1) THEN
      firstcall = 1
      nm1 = (graupel%nu+1.0)/graupel%mu
      nm2 = (graupel%nu+2.0)/graupel%mu
      CALL incgfct_lower_lookupcreate(nm1, ltable1, nlookup, nlookuphr_dummy)
      CALL incgfct_lower_lookupcreate(nm2, ltable2, nlookup, nlookuphr_dummy)
      ! ordinary gamma function of nm1 is the last value in the lookup table 1:
      g1 = ltable1%igf(ltable1%n)
      ! ordinary gamma function of nm2 is the last value in the lookup table 2:
      g2 = ltable2%igf(ltable2%n)
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD graupel_hail_conv_wet_gamlook" 
        WRITE (6,'(A,ES12.5)') "(nug+1.0)/mug :    ", nm1
        WRITE (6,'(A,ES12.5)') "(nug+2.0)/mug :    ", nm2
      ENDIF
    ELSE IF (isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD graupel_hail_conv_wet_gamlook" 
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix             
          q_c = q_cloud(i,j,k)                                      !..Massendichte
          q_r = q_rain(i,j,k)                                       !..Massendichte
          q_g = q_graupel(i,j,k)                                    !..Massendichte
          n_g = n_graupel(i,j,k)                                    !..Anzahldichte
          x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max) !..mittlere Masse
          D_g = graupel%a_geo * x_g**graupel%b_geo                  !..mittlerer Durchmesser
          n_g = q_g / x_g

          T_a = T_0(i,j,k)
          p_a = p_0(i,j,k)

          !.. Umgebungsgehalt unterkuehltes Wasser:
          qw_a = q_r + q_c

          IF (T_a < T_3 .AND. q_g > q_krit_gc .AND. qw_a > 1e-3) THEN

            !.. Umgebungsgehalt Eispartikel (vernachl. werden Graupel und Hagel wg. geringer Kollisionseff.)
            !.. koennte problematisch sein, weil in konvekt. Wolken viel mehr Graupel und Hagel enthalten ist!!!
            qi_a = q_ice(i,j,k) + q_snow(i,j,k)
            d_trenn = dmin_wg_gr_ltab_equi(p_a,T_a,qw_a,qi_a,ltabdminwgg)

            !.. Bereich im Graupelspektrum mit D > d_trenn wird zu Hagel:
            xmin = (d_trenn/graupel%a_geo)**(1.0d0/graupel%b_geo)

! UB_20081202            IF (xmin > graupel%x_min .AND. d_trenn < 10.0d0 * D_g) THEN
            IF (xmin > eps2 .AND. d_trenn < 10.0d0 * D_g) THEN

              lam = ( g2 / (g1 * x_g) )**(graupel%mu)
              n_0 = graupel%mu * n_g * lam**(nm1) / g1

              xmin = xmin**graupel%mu
              conv_n = n_0/(graupel%mu*lam**nm1)* &
                   incgfct_upper_lookup(lam*xmin, ltable1)
              conv_q = n_0/(graupel%mu*lam**nm2)* &
                   incgfct_upper_lookup(lam*xmin, ltable2)

              conv_n = MIN(conv_n,n_g)
              conv_q = MIN(conv_q,q_g)

              q_graupel(i,j,k) = q_graupel(i,j,k) - conv_q
              n_graupel(i,j,k) = n_g - conv_n

              q_hail(i,j,k) = q_hail(i,j,k) + conv_q
              n_hail(i,j,k) = n_hail(i,j,k) + conv_n

#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                dqdt(i,j,k,59) = conv_q
              END IF
#endif

            END IF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE graupel_hail_conv_wet_gamlook



END MODULE wolken_eis

!==============================================================================

MODULE wolken

  USE wolken_konstanten
  USE wolken_eis
  ! use isnan_int
  

  IMPLICIT NONE 

  DOUBLE PRECISION, PRIVATE, ALLOCATABLE, DIMENSION (:,:,:) :: ps,  ts,  qs
  DOUBLE PRECISION, PRIVATE, ALLOCATABLE, DIMENSION (:,:,:) :: q_clouds, q_ices, q_rains


CONTAINS 

  SUBROUTINE clouds ()
    !*******************************************************************************
    !                                                                              *
    !     Berechnung der Wolkenphysik mit Eisphase!                                *
    !                                                                              *
    !*******************************************************************************
    !                                                                              *
    ! Eingabe:                                                                     *
    ! dt ............. Mikrophysik-Zeitschritt in s                                *
    ! prec_s ......... Niederschlagssumme in m am Boden                            *
    ! q_cloud ........ Wolkenwasserkonzentration in kg/m**3                        *
    ! q .............. Dampfkonzentration in kg/m**3                               *
    ! q_ice .......... Wolkeneiskonzentration in kg/m**3                           *
    ! q_rain ......... Regenwasserkonzentration in kg/m**3                         *
    ! rho_0 .......... Luftdichte im Grundzustand                                  *
    ! t .............. Stoerung der Temperatur                                     *
    !                                                                              *
    !*******************************************************************************
    !                                                                              *
    ! Abgeleitete Groessen:                                                        *
    ! ts, ps ......... Abs. Temperatur und Gesamtdruck                             *
    ! qs ............. Spezifische Feuchte                                         *
    ! q_clouds, q_ices, q_rains .. Partialmassen der Hydrometeore                  *
    !                                                                              *
    !*******************************************************************************
    !                                                                              *
    ! Ausgabe:                                                                     *
    ! q_cloud, q_ice, q_rain ..... Partialdichten der Hydrometeore                 *
    ! p, q, t ........ Gleichgewichtswerte bei Saettigung                          *
    ! prec ........... Mittlere Niederschlagsrate in kg / m**2 s am Boden          *
    ! prec_s ......... Niederschlagssumme in m am Boden                            *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen, ONLY: loc_ix, loc_iy, loc_iz, dt,                &
         &                       p, t, rho, q, q_cloud, n_cloud,            &
         &                       q_ice, q_rain, q_snow, q_graupel, q_hail,  &
         &                       n_ice, n_rain, n_snow, n_graupel, n_hail,  &
         &                       prec_cloud,prec_ice, prec_rain, prec_snow, &
         &                       prec_graupel, prec_hail,                   &
         &                       prec,prec_s,T_g,p_g,rho_g,w,               &
         &                       speichere_dqdt,                &
         &                       dqdt, nrates, nmicrorates
    USE initialisierung,   ONLY: p_0, t_0, rho_0, dichte
    USE konstanten
    USE parallele_umgebung, ONLY: abortparallel,isIO,global_maxval
    USE wolken_konstanten, ONLY: satad_nach_mikrophysik, &
                                 use_rain_freeze_uli, use_ice_graupel_conv_uli

    ! .. Local Variables ..
    INTEGER          :: i, j, k, stat_var, n_dt
    DOUBLE PRECISION :: dtrho, e_sat, rdrl, rdrlm1, rlrd, rlrdm1, t_
    DOUBLE PRECISION :: dt_local
    REAL             :: wmax, qvmax,qcmax,qrmax,qimax,qsmax,qgmax,precmax,ncmax, nimax

    ! .. Skalar-Initialisierung ..
    dtrho = dt / rho_w

    rdrl  = r1 / r ; rdrlm1 = rdrl - 1.0d0
    rlrd  = r / r1 ; rlrdm1 = rlrd - 1.0d0

    IF (isdebug) WRITE (*,*) 'CLOUDS: Anfang: loc_ix, loc_iy, loc_iz = ',&
      &                               loc_ix, loc_iy, loc_iz
    IF(.NOT.ASSOCIATED(rrho_04))THEN
      ALLOCATE(rrho_04(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      ALLOCATE(rrho_c(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      IF (stat_var/=0) CALL abortparallel("Allokierungsfehler SUB(clouds(1))",4)
      !WRF! rrho_04 = (rho0/(rho_0+rho_g))**rho_vel
      !WRF! rrho_c  = (rho0/(rho_0+rho_g))**rho_vel_c
      rrho_04 = (rho0/rho_0)**rho_vel
      rrho_c  = (rho0/rho_0)**rho_vel_c
      IF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "CLOUD Wolken: initialize rrho" 
      END IF
    END IF

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD Wolken: START" 
    END IF

    WHERE (q         < 0.0d0) q         = 0.0d0
    WHERE (n_ice     < 0.0d0) n_ice     = 0.0d0
    WHERE (q_ice     < 0.0d0) q_ice     = 0.0d0
    WHERE (n_rain    < 0.0d0) n_rain    = 0.0d0
    WHERE (q_rain    < 0.0d0) q_rain    = 0.0d0
    WHERE (n_snow    < 0.0d0) n_snow    = 0.0d0
    WHERE (q_snow    < 0.0d0) q_snow    = 0.0d0
    WHERE (n_cloud   < 0.0d0) n_cloud   = 0.0d0
    WHERE (q_cloud   < 0.0d0) q_cloud   = 0.0d0
    WHERE (q_cloud  == 0.0d0) n_cloud   = 0.0d0
    WHERE (n_graupel < 0.0d0) n_graupel = 0.0d0
    WHERE (q_graupel < 0.0d0) q_graupel = 0.0d0
    WHERE (n_hail    < 0.0d0) n_hail    = 0.0d0
    WHERE (q_hail    < 0.0d0) q_hail    = 0.0d0

    !ncmax = global_maxval(n_cloud)
    !IF (isIO()) WRITE (*,'(a,2e11.4)') 'CLOUDS start: max n_cloud = ',ncmax

    IF (isdebug) WRITE (*,*) 'CLOUDS: cloud_nucleation'
    IF (cloud_typ > 1 ) THEN
      IF (nuc_c_typ .EQ. 0) THEN
        IF (isdebug) WRITE (*,*) '  ... force constant cloud droplet number conc.'
        n_cloud = qnc_const
! UB_20090316      ELSEIF (nuc_c_typ < 9) THEN
      ELSEIF (nuc_c_typ < 6) THEN
        IF (isdebug) WRITE (*,*) '  ... according to SB2006'
        CALL cloud_nucleation()
      ELSE
        IF (isdebug) WRITE (*,*) '  ... look-up tables according to Segal& Khain'
! UB_20090227>> 
!        CALL cloud_nucleation_SK()
! Use equidistant table lookup for better vectorization instead:
        CALL clnuc_sk_4D()
! UB_20090227<<
      END IF
    END IF

    IF (nuc_c_typ.ne.0) THEN
    n_cloud = MAX(n_cloud, q_cloud / cloud%x_max)
    n_cloud = MIN(n_cloud, q_cloud / cloud%x_min)
    end if

    IF (ice_typ .NE. 0) THEN

      ! Eisnukleation
      IF (isdebug) WRITE (*,*) 'CLOUDS: ice_nucleation'
      IF (nuc_i_typ > 1000) THEN
!      CALL ice_nucleation()
! UB_20090316>> vectorized version:
        CALL ice_nucleation_vec()
! UB_20090316<<
      ELSE
! AS_20090609>> new version with hom. nucleation
        CALL ice_nucleation_homhet()
      END IF

      n_ice     = MIN(n_ice, q_ice/ice%x_min)
      n_ice     = MAX(n_ice, q_ice/ice%x_max)
! AS_20090609<<

      ! Gefrieren der Wolkentropfen:
      IF (isdebug) WRITE (*,*) 'CLOUDS: cloud_freeze'
      CALL cloud_freeze ()

      ! Depositionswachstum mit dt <= 10 s
      n_dt = MAX(CEILING(dt/10.0), 1)
      dt_local = dt/n_dt
      IF (isdebug) WRITE (*,*) 'CLOUDS: vapor_deposition_growth' 
      DO i=1,n_dt
! UB_20090227>>
! Use simpler version without the fancy limiting of evaporative tendencies
! to enable vectorization:
!        CALL vapor_deposition_growth(dt_local,i)
        CALL vapor_dep_simple(dt_local,i)
! UB_20090227<<
        !CALL saturation_adjust_h2o ()
      ENDDO

      IF (isdebug) WRITE (*,*) 'CLOUDS: ice_selfcollection'
      CALL ice_selfcollection ()

      IF (isdebug) WRITE (*,*) 'CLOUDS: snow_selfcollection'
      CALL snow_selfcollection ()

      IF (isdebug) WRITE (*,*) 'CLOUDS: snow_ice_selfcollection'
      CALL snow_ice_collection ()

      IF (isdebug) WRITE (*,*) 'CLOUDS: graupel_selfcollection'
      CALL graupel_selfcollection ()

      IF (isdebug) WRITE (*,*) 'CLOUDS: graupel_ice_collection'
      CALL graupel_ice_collection ()

      IF (isdebug) WRITE (*,*) 'CLOUDS: graupel_snow_collection'
      CALL graupel_snow_collection ()

      IF (ice_typ > 1) THEN

        IF (isdebug) WRITE (*,*) 'CLOUDS: graupel_hail_conversion_wetgrowth'
! UB_20090227>> use of lookup tables for the inc. gamma-fct. 
!        CALL graupel_hail_conv_wet ()
        CALL graupel_hail_conv_wet_gamlook ()
! UB_20090227<<

        IF (isdebug) WRITE (*,*) 'CLOUDS: hail_ice_collection'
        CALL hail_ice_collection ()

        IF (isdebug) WRITE (*,*) 'CLOUDS: hail_snow_collection'
        CALL hail_snow_collection ()

      ENDIF

      IF (.NOT. use_ice_graupel_conv_uli) THEN 

        ! Old SB2006 scheme

        ! if (isdebug) write (*,*) 'CLOUDS: ice_cloud_riming'
        CALL ice_cloud_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: snow_cloud_riming'
        CALL snow_cloud_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: graupel_cloud_riming'
        CALL graupel_cloud_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: hail_cloud_riming'
        IF (ice_typ > 1) CALL hail_cloud_riming ()
        
        ! Bereifen mit Regentropfen
        ! if (isdebug) write (*,*) 'CLOUDS: ice_rain_riming'
        CALL ice_rain_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: snow_rain_riming'
        CALL snow_rain_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: graupel_rain_riming'
        CALL graupel_rain_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: hail_rain_riming'
        IF (ice_typ > 1) CALL hail_rain_riming ()

      ELSE

        ! new scheme

        ! in diesem Falle ist eine Umstellung der Reihenfolge sinnvoll, um moeglichst
        ! einfach die positive Definitheit sicherzustellen:
        
        ! Bereifen mit Wolkentropfen
        ! if (isdebug) write (*,*) 'CLOUDS: ice_cloud_riming'
        CALL ice_cloud_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: snow_cloud_riming'
        CALL snow_cloud_riming ()
        ! Bereifen mit Regentropfen
        ! if (isdebug) write (*,*) 'CLOUDS: ice_rain_riming'
        CALL ice_rain_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: snow_rain_riming'
        CALL snow_rain_riming ()
        ! an dieser Stelle werden nun die vorher lediglich gespeicherten 
        ! riming-Raten summiert:
!        CALL complete_ice_snow_riming ()
        ! vectorized version (hopefully ...) 
        CALL complete_ice_snow_rim_vec ()

        CALL graupel_cloud_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: hail_cloud_riming'
        IF (ice_typ > 1) CALL hail_cloud_riming ()
        CALL graupel_rain_riming ()
        ! if (isdebug) write (*,*) 'CLOUDS: hail_rain_riming'
        IF (ice_typ > 1) CALL hail_rain_riming ()

      ENDIF

      ! Gefrieren der Regentropfen:
      IF (isdebug) WRITE (*,*) 'CLOUDS: rain_freeze'
! UB_20080212>
!      IF (use_rain_freeze_uli .and. ice_typ > 1) THEN
      IF (use_rain_freeze_uli) THEN
! <UB_20080212
        ! new freezing scheme
        ! with a partitioning into ice/graupel/hail
! UB_20090227>> use of lookup tables for the inc. gamma-fct. 
!        CALL rain_freeze ()
        CALL rain_freeze_gamlook ()
! UB_20090227<<
      ELSE
        ! simpler SB2006 rain-to-graupel freezing scheme
        CALL rain_freeze_old ()
      END IF

      ! Schmelzen der Eispartikel
      ! if (isdebug) write (*,*) 'CLOUDS: ice_melting'
      CALL ice_melting ()
      ! if (isdebug) write (*,*) 'CLOUDS: snow_melting'
      CALL snow_melting ()
      ! if (isdebug) write (*,*) 'CLOUDS: graupel_melting'
      CALL graupel_melting ()
      ! if (isdebug) write (*,*) 'CLOUDS: hail_melting'
      IF (ice_typ > 1) CALL hail_melting ()

      ! Verdunstung von schmelzenden Eispartikeln
      ! if (isdebug) write (*,*) 'CLOUDS: graupel_evaporation'
      CALL graupel_evaporation ()
      ! if (isdebug) write (*,*) 'CLOUDS: hail_evaporation'
      IF (ice_typ > 1) CALL hail_evaporation ()
      ! if (isdebug) write (*,*) 'CLOUDS: snow_evaporation'
      CALL snow_evaporation ()

      ! Groesse der Partikel beschraenken
      n_snow    = MAX(n_snow, q_snow / snow%x_max)
      n_graupel = MAX(n_graupel, q_graupel / graupel%x_max)
      n_hail    = MAX(n_hail, q_hail / hail%x_max)
      n_ice     = MIN(n_ice, q_ice/ice%x_min)
      n_ice     = MAX(n_ice, q_ice/ice%x_max)

      ! Limit ice number conc. to 1000 per liter in any case:
      ! n_ice = MIN(n_ice, 1000d3)

    ENDIF
    IF (isdebug) CALL process_evaluation (1) ! Optional

    !ncmax = global_maxval(n_cloud)
    !IF (isIO()) WRITE (*,'(a,2e11.4)') 'CLOUDS before warm rain: max n_cloud = ',ncmax

    ! if (isdebug) write (*,*) 'CLOUDS: warm rain processes'
    ! Koagulation der Wolken- und Regentropfen
    IF (cloud_typ == 0) THEN
      !CALL autoconversionKS ()   ! Kessler (1-moment-scheme)
      !CALL accretionKS ()
    ELSE IF (cloud_typ == 1) THEN 
      CALL autoconversionSB ()   ! Seifert and Beheng (2000) (1-moment-scheme)
      CALL accretionSB ()
    ELSE IF (cloud_typ == 2) THEN 
      CALL autoconversionKS ()   ! Kessler (1969) (2-moment-scheme)
      CALL accretionKS ()
      CALL rain_selfcollectionSB ()
    ELSE IF (cloud_typ == 3 .OR. &
         &   cloud_typ == 6 .OR. &
         &   cloud_typ == 7 .OR. &
         &   cloud_typ == 8) THEN
      CALL autoconversionSB ()   ! Seifert and Beheng (2000) (2-moment-scheme)
      CALL accretionSB ()
      CALL rain_selfcollectionSB ()
    ELSE IF (cloud_typ == 4) THEN
      CALL autoconversionKK ()   ! Khai.. and Kogan (2000)
      CALL accretionKK ()
      CALL rain_selfcollectionSB ()
    ELSE IF (cloud_typ == 5) THEN
      CALL autoconversionKB ()   ! Beheng (1994)
      CALL accretionKB ()
      CALL rain_selfcollectionSB ()
    ELSE
      IF(isIO())THEN
        WRITE (6, *) "FEHLER in clouds: cloud_typ == ", cloud_typ 
      END IF
    ENDIF

    ! Verdunstung von Regentropfen
    IF (isdebug) WRITE (*,*) 'CLOUDS: rain_evaporation'
    CALL rain_evaporation ()

    ! Saettigungsadjustierung 
    ! if (isdebug) write (*,*) 'CLOUDS: saturation_adjust_h2o'
    IF (.NOT. satad_nach_mikrophysik) THEN
      CALL saturation_adjust_h2o ()
    END IF

    IF (nuc_c_typ > 0) THEN
    n_cloud = MIN(n_cloud, q_cloud/cloud%x_min)
    n_cloud = MAX(n_cloud, q_cloud/cloud%x_max)

!AS20080408>
    ! Hard upper limit for cloud number conc.
    n_cloud = MIN(n_cloud, 5000d6)
!<AS20080408
    end if

    ! Hydrometeore in ihrer Groesse beschraenken:
    n_rain = MIN(n_rain, q_rain/rain%x_min)
    n_rain = MAX(n_rain, q_rain/rain%x_max)
    IF (ice_typ > 0) THEN
    n_ice = MIN(n_ice, q_ice/ice%x_min)
    n_ice = MAX(n_ice, q_ice/ice%x_max)
    n_snow = MIN(n_snow, q_snow/snow%x_min)
    n_snow = MAX(n_snow, q_snow/snow%x_max)
    n_graupel = MIN(n_graupel, q_graupel/graupel%x_min)
    n_graupel = MAX(n_graupel, q_graupel/graupel%x_max)
    END IF
    IF (ice_typ > 1) THEN
      n_hail = MIN(n_hail, q_hail/hail%x_min)
      n_hail = MAX(n_hail, q_hail/hail%x_max)
    END IF

    ! Analyse der Prozesse
    !CALL process_evaluation (0) ! Optional

    ! diagnostische Beziehung fuer die Dichte
    ! rho = dichte(p,t,q,q_cloud,q_ice,q_rain,q_snow,q_graupel)

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD Wolken: END" 
    END IF

    if (.false.) then
    WHERE (q         < 0.0d0) q         = 0.0d0
    WHERE (n_ice     < 0.0d0) n_ice     = 0.0d0
    WHERE (q_ice     < 0.0d0) q_ice     = 0.0d0
    WHERE (n_rain    < 0.0d0) n_rain    = 0.0d0
    WHERE (q_rain    < 0.0d0) q_rain    = 0.0d0
    WHERE (n_snow    < 0.0d0) n_snow    = 0.0d0
    WHERE (q_snow    < 0.0d0) q_snow    = 0.0d0
    WHERE (n_cloud   < 0.0d0) n_cloud   = 0.0d0
    WHERE (q_cloud   < 0.0d0) q_cloud   = 0.0d0
    WHERE (n_cloud   < 0.0d0) n_cloud   = 0.0d0
    WHERE (n_graupel < 0.0d0) n_graupel = 0.0d0
    WHERE (q_graupel < 0.0d0) q_graupel = 0.0d0
    WHERE (n_hail    < 0.0d0) n_hail    = 0.0d0
    WHERE (q_hail    < 0.0d0) q_hail    = 0.0d0
    end if

  CONTAINS

    SUBROUTINE autoconversionKS ()
      !*******************************************************************************
      !                                                                              *
      !   Berechnung der Autokonversion von Wolkenwasser zu Niederschlagswasser      *
      !                                                                              *
      !*******************************************************************************

      USE globale_variablen,  ONLY: dqdt, speichere_dqdt


      IMPLICIT NONE

      !..Parameter fuer Kessler-Ansatz
      DOUBLE PRECISION, PARAMETER :: C_au = 1.00d-3
      DOUBLE PRECISION, PARAMETER :: C_qc = 5.00d-4
      DOUBLE PRECISION, PARAMETER :: x_s  = 2.60d-10     ! Trennmasse Wolken-Regentropfen

      !..Locale Variablen
      INTEGER          :: i, j, k
      DOUBLE PRECISION :: au, dt_cau, q_c

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

      ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD autoconversionKS",4)

      !..Skalar-Initialisierung
      dt_cau = dt * C_au

      rate_q = 0.0
      rate_n = 0.0
      !..Autoconversionrate nach Kessler
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            q_c = q_cloud(i,j,k)
            IF (q_c >= C_qc) THEN

              au = dt_cau * MAX(q_c - C_qc,0.d0)
              au = MIN(q_c,au)

              !..Berechnung der H_2O-Komponenten
              rate_n(i,j,k)  = au / x_s
              rate_q(i,j,k)  = au
            ENDIF
          END DO
        END DO
      END DO

      !..Zeitintegration
      n_rain  = n_rain  + rate_n 
      q_rain  = q_rain  + rate_q
      n_cloud = n_cloud - rate_n * 2.0
      q_cloud = q_cloud - rate_q
      ! ub>>
#ifdef SAVE_CONVERSIONRATES
      IF (speichere_dqdt) THEN
        dqdt(:,:,:,34) = rate_q
      END IF
#endif
      ! ub<<

      IF(isIO() .AND. isdebug)THEN
        WRITE (6, *) "CLOUD autoconversion KS: cloud_typ == ", cloud_typ 
      END IF

      DEALLOCATE(rate_n,rate_q)

    END SUBROUTINE autoconversionKS

    SUBROUTINE accretionKS ()
      !*******************************************************************************
      !                                                                              *
      !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
      !                                                                              *
      !*******************************************************************************

      USE globale_variablen,  ONLY: dqdt, speichere_dqdt

      IMPLICIT NONE

      !..Parameter fuer Kessler-Ansatz
      DOUBLE PRECISION, PARAMETER :: C_ac  = 0.94d0
      DOUBLE PRECISION, PARAMETER :: eps  = 1.00d-25

      !..Lokale Variablen
      INTEGER          :: i, j, k
      DOUBLE PRECISION :: ac,rho_q_rain,rho_vr,n_c,L_c,x_c,L_r

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q,rate_n

      ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD accretionKS",4)

      rate_n = 0.0
      rate_q = 0.0
      !..Akkreszenzrate nach Kessler (siehe Dotzek, 1999, p. 39)
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            L_c = q_cloud(i,j,k)
            N_c = n_cloud(i,j,k)
            L_r = q_rain(i,j,k)
            IF (L_c > 0.d0 .AND. L_r > 0.d0) THEN

              ac = C_ac * L_c * L_r**0.875 * rrho_04(i,j,k) * dt

              ac = MIN(L_c,ac)

              x_c = MIN(MAX(L_c/(N_c+eps),cloud%x_min),cloud%x_max) !..mittlere Masse in SI
              rate_q(i,j,k)  = ac
              rate_n(i,j,k)  = ac / x_c
            ENDIF
          END DO
        END DO
      END DO

      !..Zeitintegration
      q_rain  = q_rain  + rate_q
      q_cloud = q_cloud - rate_q
      n_cloud = n_cloud - rate_n
      ! ub>>
#ifdef SAVE_CONVERSIONRATES
      IF (speichere_dqdt) THEN
        dqdt(:,:,:,35) = rate_q
      END IF
#endif
      ! ub<<

      DEALLOCATE(rate_n,rate_q)

    END SUBROUTINE accretionKS

    SUBROUTINE process_evaluation (num)
      !*******************************************************************************
      !                                                                              *
      !        Ueberpruefung der Ergebnisse der mikrophysikalischen Prozesse,        *
      !                                                                              *
      !*******************************************************************************

      USE parallele_umgebung, ONLY:reduce_1_int,isIO,mpi_sum,abortparallel

      IMPLICIT NONE
      INTEGER :: num

      ! .. Local Variables ..
      INTEGER :: i, j, k
      INTEGER :: nq_cloud0, nq_d0, nq_ice0, nq_rain0, nq_graupel0, nq_snow0, nq_hail0 
      INTEGER :: nn_cloud0,        nn_ice0, nn_rain0, nn_graupel0, nn_snow0, nn_hail0 
      INTEGER :: nq_cloud1, nq_d1, nq_ice1, nq_rain1, nq_graupel1, nq_snow1, nq_hail1, nsum1
      INTEGER :: nn_cloud1,        nn_ice1, nn_rain1, nn_graupel1, nn_snow1, nn_hail1, nsum2
      INTEGER :: nq_cloud2, nq_d2, nq_ice2, nq_rain2, nq_graupel2, nq_snow2, nq_hail2, nsum3

      ! TEST: auf negative Werte fuer die Massendichten q_k

      ! Skalar-Initialisierung
      nq_d0       = 0
      nq_ice0     = 0
      nq_rain0    = 0
      nq_snow0    = 0
      nq_cloud0   = 0
      nq_graupel0 = 0
      nq_hail0    = 0

      nq_d1       = 0
      nq_ice1     = 0
      nq_rain1    = 0
      nq_snow1    = 0
      nq_cloud1   = 0
      nq_graupel1 = 0
      nq_hail1    = 0

      nq_d2       = 0
      nq_ice2     = 0
      nq_rain2    = 0
      nq_snow2    = 0
      nq_cloud2   = 0
      nq_graupel2 = 0
      nq_hail2    = 0

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            IF (q(i,j,k)         < 0.0d0) nq_d0       = nq_d0 + 1
            IF (q_ice(i,j,k)     < 0.0d0) nq_ice0     = nq_ice0 + 1
            IF (q_rain(i,j,k)    < 0.0d0) nq_rain0    = nq_rain0 + 1
            IF (q_snow(i,j,k)    < 0.0d0) nq_snow0    = nq_snow0 + 1
            IF (q_cloud(i,j,k)   < 0.0d0) nq_cloud0   = nq_cloud0 + 1
            IF (q_graupel(i,j,k) < 0.0d0) nq_graupel0 = nq_graupel0 + 1
            IF (q_hail(i,j,k)    < 0.0d0) nq_hail0    = nq_hail0 + 1       

          END DO
        END DO
      END DO
      CALL reduce_1_int(nq_d0,      nq_d1,      mpi_sum)
      CALL reduce_1_int(nq_ice0,    nq_ice1,    mpi_sum)
      CALL reduce_1_int(nq_rain0,   nq_rain1,   mpi_sum)
      CALL reduce_1_int(nq_snow0,   nq_snow1,   mpi_sum)
      CALL reduce_1_int(nq_cloud0,  nq_cloud1,  mpi_sum)
      CALL reduce_1_int(nq_graupel0,nq_graupel1,mpi_sum)
      CALL reduce_1_int(nq_hail0,   nq_hail1   ,mpi_sum)

      ! TEST: auf negative Werte fuer die Anzahldichten n_k

      ! Skalar-Initialisierung
      nn_ice0     = 0
      nn_rain0    = 0
      nn_snow0    = 0
      nn_cloud0   = 0
      nn_graupel0 = 0
      nn_hail0    = 0

      nn_ice1     = 0
      nn_rain1    = 0
      nn_snow1    = 0
      nn_cloud1   = 0
      nn_graupel1 = 0
      nn_hail1    = 0

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            IF (n_ice(i,j,k)     < 0.0d0) nn_ice0     = nn_ice0 + 1
            IF (n_rain(i,j,k)    < 0.0d0) nn_rain0    = nn_rain0 + 1
            IF (n_snow(i,j,k)    < 0.0d0) nn_snow0    = nn_snow0 + 1
            IF (n_cloud(i,j,k)   < 0.0d0) nn_cloud0   = nn_cloud0 + 1
            IF (n_graupel(i,j,k) < 0.0d0) nn_graupel0 = nn_graupel0 + 1
            IF (n_hail(i,j,k)    < 0.0d0) nn_hail0    = nn_hail0 + 1     

          END DO
        END DO
      END DO
      CALL reduce_1_int(nn_ice0,    nn_ice1,    mpi_sum)
      CALL reduce_1_int(nn_rain0,   nn_rain1,   mpi_sum)
      CALL reduce_1_int(nn_snow0,   nn_snow1,   mpi_sum)
      CALL reduce_1_int(nn_cloud0,  nn_cloud1,  mpi_sum)
      CALL reduce_1_int(nn_graupel0,nn_graupel1,mpi_sum)
      CALL reduce_1_int(nn_hail0,   nn_hail1,   mpi_sum)

      nsum1 = nq_d1 + nq_cloud1 + nq_ice1 + nq_rain1 + nq_snow1 + nq_graupel1 + nq_hail1
      nsum2 =         nn_cloud1 + nn_ice1 + nn_rain1 + nn_snow1 + nn_graupel1 + nn_hail1
      nsum3 = nq_d2 + nq_cloud2 + nq_ice2 + nq_rain2 + nq_snow2 + nq_graupel2 + nq_hail2

      !IF(isIO())THEN
      IF (nsum1 /= 0 .OR. nsum2 /= 0) THEN
        IF (nsum1 /= 0) THEN
          WRITE (6, *) "CLOUD proc_eval",num,": q ", nq_d1, "mal, q_cloud: ", nq_cloud1, &
               &    "mal, q_ice: ", nq_ice1,                                   &
               &    "mal, q_rain: ", nq_rain1,                                 &
               &    "mal, q_snow: ", nq_snow1,                                 &
               &    "mal, q_hail: ", nq_hail1,                                 &
               &    "mal, q_graupel: ", nq_graupel1, "mal < 0!"
        ENDIF
        IF (nsum2 /= 0) THEN
          WRITE (6, *) "CLOUD proc_eval",num,": n_ice: ", nn_ice1,                 &
               &    "mal, n_cloud: ", nn_cloud1,                               &
               &    "mal, n_rain: ", nn_rain1,                                 &
               &    "mal, n_snow: ", nn_snow1,                                 &
               &    "mal, n_hail: ", nn_hail1,                                 &
               &    "mal, n_graupel: ", nn_graupel1, "mal < 0!"
        ENDIF
      ELSE
        !WRITE (6, *) "CLOUD proc_eval:  ok"
      ENDIF
      !END IF

    END SUBROUTINE process_evaluation

    SUBROUTINE saturation_adjust_h2o ()
      !*******************************************************************************
      !                                                                              *
      ! Diese Routine fuehrt in JEDEM Fall eine Adjustierung der H2O-Komponenten und *
      ! der Temperatur auf ihre Gleichgewichtswerte hin durch. Dabei wird angenommen,*
      ! dass sofort bei 100% Saettigung Kondensation eintritt und dass die Phasenum- *
      ! wandlungen vollstaendig reversibel ablaufende thermodynamische Prozesse sind *
      !                                                                              *
      ! Die Saettigungsadjustierung ist immer der LETZTE Schritt eines Wolkenmodells *
      ! Alle anderen Prozesse, wie Advektion, Diffusion und Mikrophysik werden davor *
      ! berechnet. Die Saettigungsadjustierung bestimmt dann nur noch die Nukleation *
      ! die max. zum Erreichen von 100% rel. Feuchte im Wolkenraum erforderlich ist. *
      !                                                                              *
      ! Die tatsaechliche rel. Feuchte in Wolken kann mit der Subroutine saturation_ *
      ! evaluation ermittelt werden. Da das direkte Adjustierungsverfahren eine Nae- *
      ! herung darstellt, wird nicht immer exakt 100% Feuchte erreicht.  Erfahrungs- *
      ! gemaess liegen die Werte aber ungefaehr im Intervall                         *
      !                                                                              *
      !                        [ 99.85 % <= RH <= 100.03 % ]                         *
      !                                                                              *
      ! oder sogar noch naeher bei 100%. Da auch in der Realitaet ganz leichte Unter *
      ! und Uebersaettigungen in der Wolke beobachtet werden, ist dieses Resultat OK *
      !                                                                              *
      ! Die hier programmierte Saettigungsadjustierung setzt die Dampfdruckformel in *
      ! der Form nach Tetens (1930) voraus. Siehe UP clouds sowie die Referenz Te30! *
      !                                                                              *
      !*******************************************************************************

      USE globale_variablen,  ONLY: dqdt, speichere_dqdt, &
        & cond_neu_sb, evap_neu_sb, speichere_precipstat

      ! .. Local Scalars ..
      INTEGER          :: stat_var = 0
      INTEGER          :: i, j, k
      DOUBLE PRECISION :: a_, a1, b_, b1, dq_d, lrcv, q_sat, rcv, lrcp, T_a, p_a, rho_a, hlp

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q
      ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
      IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD saturation_adj",4)

      IF (isIO() .AND. isdebug) WRITE (6, *) "CLOUD saturation_adj " 

      rate_q = 0.0

      rcv  = 1.0d0 /  cv
      lrcv = L_wd  * rcv
      lrcp = L_wd  / cp

      ! diagnostische Beziehung fuer die Dichte
      !rho = dichte(p,t,q,q_cloud,q_ice,q_rain,q_snow,q_graupel)

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            T_a   = T_0(i,j,k) !WRF!+ T_0(i,j,k) + T_g(i,j,k)

            IF (T_a > T_f .OR. ice_typ == 0) THEN
              p_a   = p_0(i,j,k) !WRF!+ p_0(i,j,k) + p_g(i,j,k)
              rho_a = rho_0(i,j,k) !WRF!+ rho_g(i,j,k)

              a_    = A_w
              b_    = B_w
              e_sat = e_ws (T_a) * 1.00

              !...Berechnung der spezifischen Saettigungsfeuchte, Dotzek (4.33)
              q_sat = rlrd / (p_a / e_sat + rlrdm1)

              !...Berechnung der Adjustierungs-Terme, Dotzek, S. 35
              a1 = a_ * (T_3 - b_) / (T_a - b_)**2
              b1 = a1 * lrcp * q_sat

              !...Berechnung des Phasenuebergangs, Dotzek, S. 36
              dq_d = -(q(i,j,k) - q_sat*rho_a) / (1.0d0 + b1)

              !...Verdunstung bis zur Saettigung oder bis kein Kondensat mehr da ist
              dq_d = MAX( MIN (dq_d, q_cloud(i,j,k)), -q(i,j,k))

              rate_q(i,j,k) = dq_d

              !...Berechnung der Verdampfungswaerme
              !WRF!T(i,j,k) = T(i,j,k) - L_wd/cp * dq_d / rho_a
              !WRF!p(i,j,k) = p(i,j,k) - L_wd/cv * dq_d * R_l

              ! ub>>
              !...Berechnung der Uebersaettigung vor der Adjustierung:
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                hlp = q(i,j,k) * p_a / (e_sat*(rlrd-rlrdm1*q(i,j,k)/rho_a)) / rho_a - 1.0
                IF (hlp < 0.0d0) THEN
                  dqdt(i,j,k,39) = hlp
                ELSE
                  dqdt(i,j,k,40) = hlp
                END IF
              END IF
#endif
              ! ub<<

              !...Berechnung der neuen H_2O-Komponenten
              q(i,j,k)       = q(i,j,k)       + dq_d
              q_cloud(i,j,k) = q_cloud(i,j,k) - dq_d
              ! ub>>
#ifdef SAVE_CONVERSIONRATES
              IF (speichere_dqdt) THEN
                IF (dq_d >= 0.0) THEN
                  dqdt(i,j,k,37) = dq_d
                ELSE
                  dqdt(i,j,k,38) = -dq_d
                END IF
              END IF
#endif
              IF (speichere_precipstat) THEN
                IF (dq_d >= 0.0) THEN
                  evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) + dq_d
                ELSE
                  cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) - dq_d
                END IF
              END IF
              ! ub<<
            END IF

          END DO
        END DO
      END DO

      DEALLOCATE(rate_q)

    END SUBROUTINE saturation_adjust_h2o

    SUBROUTINE saturation_adjust_ice ()
      !*******************************************************************************
      !                                                                              *
      ! Diese Routine fuehrt in JEDEM Fall eine Adjustierung der H2O-Komponenten und *
      ! der Temperatur auf ihre Gleichgewichtswerte hin durch. Dabei wird angenommen,*
      ! dass sofort bei 100% Saettigung Kondensation eintritt und dass die Phasenum- *
      ! wandlungen vollstaendig reversibel ablaufende thermodynamische Prozesse sind *
      !                                                                              *
      ! Im Gegensatz zum rein massengewichteten Verfahren aus Ta89 wird hier ein ei- *
      ! gen-entwickeltes Zweischrittverfahren eingefuehrt, bei dem im ersten Schritt *
      ! die Saettigungsadjustierung ueber Wasser gerechnet wird, bei der auch allein *
      ! Fluessigwasser entstehen kann. (Solche Uebersaettigungen werden unter 0 Grad *
      ! Celsius aber kaum je erreicht...). Im zweiten Schritt folgt die Adjustierung *
      ! bzgl. Eis. Liegt hier Untersaettigung vor, sublimiert Eis (Wasser kann hier- *
      ! bei *NICHT* mehr im System sein!), liegt Uebersaettigung vor, wird unterhalb *
      ! von T_f = - 40 Grad C nur Eis, zwischen T_f und T_3 Eis und  Wasser in einem *
      ! temperaturabhaengigen Verhaeltnis (z.Z. bilinear wie auch in Ta89) gebildet. *
      !                                                                              *
      ! Vorteile gegenueber Ta89 sind hierbei:                                       *
      !                                                                              *
      !  - es bleiben keine Uebersaettigungen bzgl. Eis uebrig                       *
      !  - das Verfahren funktioniert immer, d. h. auch bei der Bildung              *
      !    von Eis- und Wasserteilchen in einer entstehenden Mischwolke              *
      !                                                                              *
      ! Die Saettigungsadjustierung ist immer der LETZTE Schritt eines Wolkenmodells *
      ! Alle anderen Prozesse, wie Advektion, Diffusion und Mikrophysik werden davor *
      ! berechnet. Die Saettigungsadjustierung bestimmt dann nur noch die Nukleation *
      ! die max. zum Erreichen von 100% rel. Feuchte im Wolkenraum erforderlich ist. *
      !                                                                              *
      ! Die tatsaechliche rel. Feuchte in Wolken kann mit der Subroutine saturation_ *
      ! evaluation ermittelt werden. Da das direkte Adjustierungsverfahren eine Nae- *
      ! herung darstellt, wird nicht immer exakt 100% Feuchte erreicht.  Erfahrungs- *
      ! gemaess liegen die Werte aber ungefaehr im Intervall                         *
      !                                                                              *
      !                        [ 99.85 % <= RH <= 100.15 % ]                         *
      !                                                                              *
      ! oder sogar noch naeher bei 100%. Da auch in der Realitaet ganz leichte Unter *
      ! und Uebersaettigungen in der Wolke beobachtet werden, ist dieses Resultat OK *
      !                                                                              *
      ! Die hier programmierte Saettigungsadjustierung setzt die Dampfdruckformel in *
      ! der Form nach Tetens (1930) voraus. Siehe UP clouds sowie die Referenz Te30! *
      !                                                                              *
      !*******************************************************************************

      ! .. Local Variables ..
      INTEGER          :: i, j, k
      DOUBLE PRECISION :: a1, a2, a3, b1, b2, cond, depo, dq_d, e_sat, lrcp, lrcv
      DOUBLE PRECISION :: qes, qws, rcp, rcv, rt3mtf

      ! Skalar-Initialisierung
      rcp    = 1.0d0 / cp ; lrcp = L_wd * rcp
      rcv    = 1.0d0 / cv ; lrcv = L_wd * rcv
      rt3mtf = 1.0d0 / (T_3 - T_f)

      DO k = 1, loc_iz

        !   Schichtweise Adjustierung
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            !       1. Teilschritt:
            !       ===============
            IF (ts(i, j, k) > T_f) THEN ! Saettigung ueber Wasser

              e_sat = e_ws (ts(i, j, k))

              !         Berechnung der spezifischen Saettigungsfeuchte
              qws = rlrd / (ps(i, j, k) / e_sat + rlrdm1)

              !         Berechnung der Terme
              a1 = A_w * (T_3 - B_w) / (ts(i, j, k) - B_w) ** 2
              b1 = a1  * lrcp * qws

              !         Berechnung der Art und des Betrags des Phasenuebergangs
              dq_d = -(qs(i, j, k) - qws) / (1.0d0 + b1)

              !         Verdunstung bis zur Saettigung oder bis kein Kondensat mehr da ist
              dq_d = MIN (dq_d, q_clouds(i, j, k))

              !         Berechnung der Verdampfungswaerme
              t(i, j, k) = t(i, j, k) -                      lrcp * dq_d
              p(i, j, k) = p(i, j, k) - r * rho_0(i, j, k) * lrcv * dq_d

              !         Berechnung der neuen H_2O-Komponenten
              qs(i, j, k)       = qs(i, j, k)       + dq_d
              q_clouds(i, j, k) = q_clouds(i, j, k) - dq_d

              !         Neuberechnung abs. Temp. und ges. Druck
              ts(i, j, k) = t(i, j, k) + t_0(i, j, k) + &
                   p(i, j, k) / (cp * rho_0(i, j, k))
              !         ps(i, j, k) = p(i, j, k) + p_0(i, j, k)

            END IF
          END DO
        END DO

      END DO

      !       2. Teilschritt:
      !       ===============
      DO k = 1, loc_iz

        !   Schichtweise Adjustierung
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (ts(i, j, k) < T_3) THEN ! Saettigung ueber Eis

              e_sat = e_es (ts(i, j, k))

              !..Berechnung der spezifischen Saettigungsfeuchten
              qes = rlrd / (ps(i, j, k) / e_sat + rlrdm1)

              !..Berechnung der Terme (Welcher Terme?)
              a2 = A_e * (T_3 - B_e) / (ts(i, j, k) - B_e) ** 2

              b1 = qs(i, j, k) - qes
              b2 = a2 * qes

              IF (b1 > 0.0d0) THEN
                cond = MAX (ts(i, j, k) - T_f, 0.0d0) * rt3mtf ! Wichtungsfaktoren
                depo = 1.0d0 - cond                            ! Eis <-> Wasser...

                a3 = (L_wd * cond + (L_ew + L_wd) * depo)
              ELSE
                a3 =              (L_ew + L_wd)
              END IF

              !..Berechnung der Art und des Betrags des Phasenuebergangs
              dq_d = -b1 / (1.0d0 + b2 * a3 * rcp)

              !..Verdunstung bis zur Saettigung oder bis kein Kondensat mehr da ist
              dq_d = MIN (dq_d, q_ices(i, j, k))

              !..Berechnung der Verdampfungswaerme
              t(i, j, k) = t(i, j, k) -                        a3 * rcp * dq_d
              p(i, j, k) = p(i, j, k) - r * rho_0(i, j, k) *   a3 * rcv * dq_d

              !..Berechnung der neuen H_2O-Komponenten
              qs(i, j, k) = qs(i, j, k) + dq_d

              IF (dq_d < 0.0d0) THEN
                q_clouds(i, j, k) = q_clouds(i, j, k) - cond * dq_d
                q_ices(i, j, k)   = q_ices(i, j, k)   - depo * dq_d
              ELSE
                q_ices(i, j, k) = q_ices(i, j, k) - dq_d
              END IF

              !         Neuberechnung abs. Temp. und ges. Druck
              ts(i, j, k) = t(i, j, k) + t_0(i, j, k) + &
                   p(i, j, k) / (cp * rho_0(i, j, k))
              ps(i, j, k) = p(i, j, k) + p_0(i, j, k)

            END IF

          END DO
        END DO

      END DO

    END SUBROUTINE saturation_adjust_ice

    SUBROUTINE saturation_evaluation ()
      !*******************************************************************************
      !                                                                              *
      !           Ueberpruefung der Ergebnisse der Saettigungsadjustierung           *
      !                                                                              *
      !*******************************************************************************

      USE parallele_umgebung,ONLY: reduce_1_int,reduce_1, &
           mpi_sum,mpi_max,mpi_min,isIO
      ! .. Local Scalars ..
      INTEGER          :: i, j, k
      DOUBLE PRECISION :: rh, rh_max, rh_min,rh_max0,rh_min0,p_a,T_a,rho_a,q_a

      rh_min = 1.0d1
      rh_max = 0.0d0

      ! Treten auffaellige Werte von rh im Wolkenraum auf?
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (q_cloud(i,j,k) > 0.0d0) THEN
              p_a  = p(i,j,k) + p_0(i,j,k) + p_g(i,j,k)
              T_a  = T(i,j,k) + T_0(i,j,k) + T_g(i,j,k)
              rho_a = rho_0(i,j,k) + rho_g(i,j,k)
              q_a  = q(i,j,k) / rho_a

              !         Finden des Saettigungsdampfdrucks
              e_sat = e_ws (T_a)

              !         Rel. Feuchte
              rh =  p_a * q_a / (e_sat * (rlrd - rlrdm1 * q_a))

              rh_max = MAX (rh, rh_max)
              rh_min = MIN (rh, rh_min)
            END IF

          END DO
        END DO
      END DO

      ! Rel. Feuchte in %
      rh_max = 1.0d2 * rh_max
      rh_min = 1.0d2 * rh_min
      CALL reduce_1(rh_max,rh_max0,mpi_max)
      CALL reduce_1(rh_min,rh_min0,mpi_min)

      ! Ergebnis der 2. Fehlerdiagnose
      IF(isIO())THEN
        !IF (rh_max0 - rh_min0 > 1.0d0) THEN
        WRITE (6,'(A, 2(F8.3, A))') " CLOUD saturation_eval: Wolken-Luftfeuchte: [", &
             &                              rh_min0, " : ", rh_max0, " ] %"
        !END IF
      END IF

    END SUBROUTINE saturation_evaluation

  END SUBROUTINE clouds

  SUBROUTINE alloc_wolken ()
    !*******************************************************************************
    !           Allokierung der Wolken- und Niederschlagsfelder                    *
    !*******************************************************************************

    USE globale_variablen,  ONLY: q_cloud, q_rain, q_ice, q_snow, q_graupel,  &
         &                        n_cloud, n_rain, n_ice, n_snow, n_graupel,  &
         &                        prec_cloud,prec_rain,prec_ice,              &
         &                        prec_snow,prec_graupel,prec,prec_s,         &
         &                        q_hail, n_hail, prec_hail,                  & ! <hn
         &                        loc_ix, loc_iy, loc_iz,                     &
                                ! ub>>
         &                        dqdt, nmicrorates, speichere_dqdt, &
         &                        cond_neu_sb, evap_neu_sb, speichere_precipstat
    ! ub<<
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: rho_0    
    USE konstanten,         ONLY: wolke_typ

    IMPLICIT NONE

    INTEGER       :: stat_var = 0
    INTEGER, SAVE :: firstcall

    ice_typ   = wolke_typ/1000
    nuc_i_typ = MOD(wolke_typ/100,10)
    nuc_c_typ = MOD(wolke_typ/10,10)
    cloud_typ = MOD(wolke_typ,10)
    !cloud_typ = 3

    IF(isIO() .AND. firstcall.NE.1) THEN
      WRITE (6,*) "CLOUD alloc_wolken: wolke_typ = ", wolke_typ
      WRITE (6,*) "                    ice_typ   = ", ice_typ
      WRITE (6,*) "                    nuc_i_typ = ", nuc_i_typ
      WRITE (6,*) "                    nuc_c_typ = ", nuc_c_typ
      WRITE (6,*) "                    cloud_typ = ", cloud_typ
      firstcall = 1
    ENDIF

    ALLOCATE (n_cloud(0:loc_ix, 1:loc_iy, 1:loc_iz), &
         q_cloud(0:loc_ix, 1:loc_iy, 1:loc_iz), &
         n_rain(0:loc_ix, 1:loc_iy, 1:loc_iz),  &
         q_rain(0:loc_ix, 1:loc_iy, 1:loc_iz),  &
         n_ice(0:loc_ix, 1:loc_iy, 1:loc_iz),   &
         q_ice(0:loc_ix, 1:loc_iy, 1:loc_iz),   &             
         n_snow(0:loc_ix, 1:loc_iy, 1:loc_iz),   &
         q_snow(0:loc_ix, 1:loc_iy, 1:loc_iz),   &             
         n_hail(0:loc_ix, 1:loc_iy, 1:loc_iz),   &  ! <hn
         q_hail(0:loc_ix, 1:loc_iy, 1:loc_iz),   &  ! <hn
         n_graupel(0:loc_ix, 1:loc_iy, 1:loc_iz),   &
         q_graupel(0:loc_ix, 1:loc_iy, 1:loc_iz), STAT = stat_var)


    IF (stat_var /= 0) THEN
      CALL abortparallel ("Fehler bei Allokierung Wolkenvariablen", stat_var)
    END IF

    q_cloud   = 0.0
    q_rain    = 0.0
    q_ice     = 0.0
    q_snow    = 0.0
    q_graupel = 0.0
    q_hail    = 0.0   ! <hn
    n_cloud   = 0.0
    n_rain    = 0.0
    n_ice     = 0.0
    n_snow    = 0.0
    n_hail    = 0.0   ! <hn
    n_graupel = 0.0

    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN

      ALLOCATE (dqdt(0:loc_ix, 1:loc_iy, 1:loc_iz, nmicrorates), STAT = stat_var)
      IF (stat_var /= 0) THEN
        CALL abortparallel ("Fehler bei Allokierung der Umwandlungsraten", stat_var)
      END IF

      dqdt = 0.0
    END IF
#endif
    IF (speichere_precipstat) THEN
      ALLOCATE(cond_neu_sb(0:loc_ix, 1:loc_iy, 1:loc_iz), evap_neu_sb(0:loc_ix, 1:loc_iy, 1:loc_iz), STAT = stat_var)
      IF (stat_var /= 0) THEN
        CALL abortparallel ("Fehler bei Allokierung cond_neu_sb, evap_neu_sb", stat_var)
      END IF
      cond_neu_sb = 0.0
      evap_neu_sb = 0.0
    END IF
    ! ub<<       

    !Allokierung normalerweise in alloc_boden
    ALLOCATE (prec(0:loc_ix, 1:loc_iy), &
         prec_s(0:loc_ix, 1:loc_iy), &
         prec_cloud(0:loc_ix, 1:loc_iy), &
         prec_rain(0:loc_ix, 1:loc_iy),  &
         prec_ice(0:loc_ix, 1:loc_iy),   &
         prec_snow(0:loc_ix, 1:loc_iy),  &
         prec_hail(0:loc_ix, 1:loc_iy),  &   ! <hn
         prec_graupel(0:loc_ix, 1:loc_iy), STAT = stat_var)
    IF (stat_var /= 0) THEN
      CALL abortparallel ("Fehler bei Allokierung Regenraten", stat_var)
    END IF




    prec   = 0.0
    prec_s = 0.0
    prec_cloud = 0.0
    prec_rain  = 0.0
    prec_ice   = 0.0
    prec_snow  = 0.0
    prec_hail  = 0.0   ! <hn
    prec_graupel = 0.0

    stat_var = 0

  END SUBROUTINE alloc_wolken

  SUBROUTINE dealloc_wolken ()
    !****************************************************************************
    !           Allokierung der Wolken- und Niederschlagsfelder                 *
    !****************************************************************************

    USE globale_variablen,  ONLY: q_cloud, q_rain, q_ice, q_snow, q_graupel,  &
         &                        n_cloud, n_rain, n_ice, n_snow, n_graupel,  &
         &                        prec_cloud,prec_rain,prec_ice,              &
         &                        prec_snow,prec_graupel,                     &
         &                        q_hail, n_hail, prec_hail,                  &  ! <hn
         &                        loc_ix, loc_iy, loc_iz, prec, prec_s, rho_g, &
                                ! ub>>        
         &                        dqdt, speichere_dqdt, cond_neu_sb, evap_neu_sb, speichere_precipstat
    ! ub<<
    USE parallele_umgebung, ONLY: isIO
    IMPLICIT NONE

    DEALLOCATE (n_cloud,q_cloud,n_rain,q_rain,n_ice,q_ice, n_snow, &
         &  q_snow,n_graupel,q_graupel, q_hail, n_hail)

    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) DEALLOCATE (dqdt)
#endif
    IF (speichere_precipstat) DEALLOCATE(cond_neu_sb, evap_neu_sb)
    ! ub<<

    DEALLOCATE (prec,prec_s,prec_cloud, &
         &  prec_rain,prec_ice,prec_snow,prec_graupel, prec_hail)

  END SUBROUTINE dealloc_wolken

  SUBROUTINE autoconversionSB ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Autokonversion und Selfcollection der Wolkentropfen         *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, dt, &
         &                        q_cloud, n_cloud,q_rain, n_rain, dqdt, speichere_dqdt
    USE konstanten,         ONLY: pi
    USE parallele_umgebung, ONLY: isIO,abortparallel

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: eps  = 1.00d-25

    !..Locale Variablen
    INTEGER          :: stat_var = 0
    INTEGER          :: i,j,k
    INTEGER, SAVE    :: firstcall
    DOUBLE PRECISION :: q_c, q_r, q_g, n_c, x_c, nu, mu, &
         & tau, phi, k_au, k_sc, x_s, au, sc, k_c, k_1, k_2

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_nc,rate_nr,rate_q

    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_nr(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_nc(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD autoconversionSB",4)

    !..Skalar-Initialisierung

    IF (cloud_typ <= 3) THEN
      !..Parameter fuer Seifert & Beheng (2001)
      k_c  = 9.44d+9   !..Long-Kernel
      k_1  = 6.00d+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.68d+0   !..Parameter fuer Phi-Fkt.
    ELSEIF (cloud_typ == 6) THEN
      !..Parameter fuer Pinsky et al (2000) Kernel
      k_c  = 4.44d+9   !..CC-Kernel
      k_1  = 4.00d+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.70d+0   !..Parameter fuer Phi-Fkt.
    ELSEIF (cloud_typ == 7) THEN
      !..Parameter fuer Pinsky et al (2000) Kernel
      k_c  = 10.58d+9  !..CC-Kernel
      k_1  = 4.00d+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.70d+0   !..Parameter fuer Phi-Fkt.
    ELSEIF (cloud_typ == 8) THEN
      !..Parameter fuer HUCM Kernel
      k_c  = 30.00d+9  !..CC-Kernel
      k_1  = 4.00d+2   !..Parameter fuer Phi-Fkt.
      k_2  = 0.70d+0   !..Parameter fuer Phi-Fkt.
    ENDIF

    nu    = cloud%nu
    mu    = cloud%mu
    x_s   = cloud%x_max                     !..Trennmasse

    IF (.true.) THEN 
    !IF (mu == 1.0) THEN 
      k_au  = k_c / (2.0d1*x_s) * (nu+2.0d0)*(nu+4.0d0)/(nu+1.0d0)**2
      k_sc  = k_c * (nu+2.0d0)/(nu+1.0d0)
    ELSE
      k_au = k_c / (2.0d1*x_s)                                       &
        & * ( 2.0 * gfct((nu+4.0)/mu)**1                          &
        &         * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2   &
        &   - 1.0 * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 ) &
        &   / gfct((nu+2.0)/mu)**4
      k_sc = k_c * moment_gamma(cloud,2)
    ENDIF

    IF (.false. .and. isIO()) THEN
       WRITE (*,'(a,2e11.4)') 'CLOUDS autocon: nu  = ',nu
       WRITE (*,'(a,2e11.4)') 'CLOUDS autocon: mu  = ',mu
       WRITE (*,'(a,2e11.4)') 'CLOUDS autocon: kau = ',k_au
       WRITE (*,'(a,2e11.4)') 'CLOUDS autocon: xmin = ',cloud%x_min
       WRITE (*,'(a,2e11.4)') 'CLOUDS autocon: xmax = ',cloud%x_max
    END IF

    rate_q  = 0.0
    rate_nc = 0.0
    rate_nr = 0.0

    !..Parametrisierung nach Seifert und Beheng (2000)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          IF (q_c > eps) THEN
            n_c = n_cloud(i,j,k)                                  !..Anzahldichte
            q_r = q_rain(i,j,k)                                   !..Fluessigwassergehalt
            IF (nuc_c_typ .EQ. 0) THEN
              x_c = q_c/qnc_const
            ELSE
              x_c = MIN(MAX(q_c/(n_c+eps),cloud%x_min),cloud%x_max) !..mittlere Masse in SI
            end if

            !..Berechnung der Autokonversionsrate nach SB2000
            au  =  k_au * q_c**2 * x_c**2 * dt 
            !au  = k_au * q_c**2 * x_c**2 * dt * rrho_c(i,j,k)
            IF (q_c > 1.0d-6) THEN
              tau = MIN(MAX(1.0-q_c/(q_c+q_r+eps),eps),0.9d0)
              phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
              au  = au * (1.0 + phi/(1.0 - tau)**2)
            ENDIF
            au = MAX(MIN(q_c,au),0d0)

            sc = k_sc * q_c**2 * dt 
            !sc = k_sc * q_c**2 * dt * rrho_c(i,j,k)

            rate_q(i,j,k)  = au
            rate_nr(i,j,k) = au / x_s
            rate_nc(i,j,k) = MIN(n_c,sc)     !..Selfcollection und Autokonversion

          ENDIF
        END DO
      END DO
    END DO

    !..Zeitintegration
    n_rain  = n_rain  + rate_nr 
    q_rain  = q_rain  + rate_q
    n_cloud = n_cloud - rate_nc
    q_cloud = q_cloud - rate_q
    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,34) = rate_q
    END IF
#endif
    ! ub<<

    !IF(isIO() .AND. isdebug) THEN
    IF(isIO()) THEN
      IF (firstcall.NE.1) THEN
        WRITE (6,*) "CLOUD autoconversion SB:"
        WRITE (6,'(A,D10.3)') "    Trennmasse Wolken-Regentropfen:   x_s  = ", x_s
        WRITE (6,'(A,D10.3)') "    Koeff. fuer Autokonversionsrate:  k_au = ", k_au
        firstcall = 1
      ENDIF
    END IF

    DEALLOCATE(rate_nr,rate_nc,rate_q)

  END SUBROUTINE autoconversionSB

  SUBROUTINE accretionSB ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, dt,              &
         &                        q_cloud, n_cloud, q_rain, n_rain, dqdt,  &
         &                        speichere_dqdt
    USE konstanten,         ONLY: pi
    USE parallele_umgebung, ONLY: isIO,abortparallel

    IMPLICIT NONE

    !..Parameter fuer Seifert & Beheng (2001)
    DOUBLE PRECISION, PARAMETER :: k_r = 5.78d+0   ! Parameter Kernel
    DOUBLE PRECISION, PARAMETER :: k_1 = 5.00d-4   ! Parameter fuer Phi-Fkt.

    !..Parameter fuer Seifert (2002), not recommended
    !DOUBLE PRECISION, PARAMETER :: k_r = 5.25d+0   ! Parameter Kernel
    !DOUBLE PRECISION, PARAMETER :: k_1 = 5.00d-5   ! Parameter fuer Phi-Fkt.
    DOUBLE PRECISION, PARAMETER :: eps = 1.00d-25

    !..Lokale Variablen
    INTEGER          :: stat_var = 0
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: ac
    DOUBLE PRECISION :: L_c, L_r, L_g, tau, phi, n_c, x_c

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q,rate_n

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD accretionnSB",4)

    x_c  = 4./3. * pi * rho_w * r_c**3     !..Mittlere Masse der Wolkentropfen
    rate_n = 0.0
    rate_q = 0.0

    !..Parametrisierung nach Seifert und Beheng (2001)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_c = n_cloud(i,j,k) !..Anzahldichte
          L_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          L_r = q_rain(i,j,k)  !..Fluessigwassergehalt

          IF (L_c > 0.0d0.AND.L_r > 0.0d0) THEN

            !..Berechnung der Akkreszenzrate nach SB2001
            tau = MIN(MAX(1.0-L_c/(L_c+L_r+eps),eps),1.0d0)
            phi = (tau/(tau+k_1))**4
            ac  = k_r *  L_c * L_r * phi * dt
            !ac  = k_r *  L_c * L_r * phi  * rrho_04(i,j,k) * dt

            ac = MIN(L_c,ac)

            x_c = MIN(MAX(L_c/(n_c+eps),cloud%x_min),cloud%x_max)
            rate_q(i,j,k)  = ac
            rate_n(i,j,k)  = MIN(n_c,ac/x_c)
          ENDIF
        END DO
      END DO
    END DO

    !..Zeitintegration
    q_rain  = q_rain  + rate_q
    q_cloud = q_cloud - rate_q
    n_cloud = n_cloud - rate_n
    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,35) = rate_q
    END IF
#endif
    ! ub<<

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE accretionSB

  SUBROUTINE rain_selfcollectionSB ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Selfcollection von Regentropfen                         *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, n_rain, q_rain, dt
    USE konstanten,         ONLY: pi 
    USE parallele_umgebung, ONLY: isIO

    IMPLICIT NONE

    !..Parameter fuer Seifert & Beheng (2001)
    DOUBLE PRECISION, PARAMETER :: k_r = 5.78d+00

    !..Parameter fuer Seifert (2002), DO NOT USE!
    !DOUBLE PRECISION, PARAMETER :: k_sc = 2.87d+00
    !DOUBLE PRECISION, PARAMETER :: k_rr = 7.12d+00
    !DOUBLE PRECISION, PARAMETER :: k_br = 1.00d+03
    !DOUBLE PRECISION, PARAMETER :: D_br = 0.90d-03
    !DOUBLE PRECISION, PARAMETER :: kap  = 6.07d+01

    !..Parameter fuer Seifert (2007)
    DOUBLE PRECISION, PARAMETER :: D_br = 1.10d-03
    DOUBLE PRECISION, PARAMETER :: k_rr = 4.33d+00
    DOUBLE PRECISION, PARAMETER :: k_br = 1.00d+03

    DOUBLE PRECISION, PARAMETER :: eps = 1.00d-25
    DOUBLE PRECISION, PARAMETER :: rho = 1.00d+03   ! in SI [kg/m^3]     

    DOUBLE PRECISION, PARAMETER :: d_min = 80.0d-6  ! in SI [m]     
    DOUBLE PRECISION, PARAMETER :: d_max = 1.00d-2  ! in SI [m]    

    !..Lokale Variablen
    INTEGER          :: i, j, k
    DOUBLE PRECISION :: sc, br, kapexp
    DOUBLE PRECISION :: q_r, n_r, x_r, d_r, lam, phi1, phi2

    !..Parametrisierung nach Seifert und Beheng

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          n_r = n_rain(i,j,k)    !..Anzahldichte
          q_r = q_rain(i,j,k)    !..Fluessigwassergehalt

          IF (q_r > 0.0d0) THEN
            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max)
            D_r = rain%a_geo * x_r**rain%b_geo
            !lam = lambda_gamma(rain,x_r)

            !..Berechnung der Selfcollection nach S2002
            !sc = k_r *  n_r * q_r * rrho_04(i,j,k) * dt
            !sc = k_sc *  n_r * q_r * rrho_04(i,j,k) * dt
            !sc = k_rr *  n_r * q_r * (1.0 + kap/lam)**(-9) * rrho_04(i,j,k) * dt
            sc = k_rr *  n_r * q_r * rrho_04(i,j,k) * dt

            !..Berechnung Breakup nach S2002
            br = 0.0
            IF (D_r.GT.0.30e-3) THEN
              phi1 = (k_br * (D_r - D_br) + 1.0)
              !phi2 = 2.0 * exp(2.3e3 * (D_r - D_br)) - 1.0
              !br = max(phi1,phi2) * sc                   
              br = phi1 * sc                   
            ENDIF
            sc = MIN(n_r,sc-br)

            n_r  = n_r  - sc

            ! Untere und obere Schranke fuer d_rain bzw. n_rain
            x_r = MIN(MAX(q_r/(n_r+eps),rain%x_min),rain%x_max) !..mittlere Masse in SI
            n_rain(i,j,k) = q_rain(i,j,k) / x_r
          ENDIF
        END DO
      END DO
    END DO
    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD rain_selfcollection SB"
    END IF

  END SUBROUTINE rain_selfcollectionSB

  SUBROUTINE autoconversionKB ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Autokonversion von Wolkenwasser zu Niederschlagswasser      *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, q_cloud, q_rain, n_rain, dt
    USE konstanten,         ONLY: pi
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE globale_variablen,  ONLY: dqdt, speichere_dqdt

    IMPLICIT NONE

    !..Parameter fuer Beheng (1994)
    DOUBLE PRECISION, PARAMETER :: eps  = 1.00d-25

    !..Locale Variablen
    INTEGER          :: stat_var = 0
    INTEGER          :: i,j,k
    INTEGER, SAVE    :: firstcall
    DOUBLE PRECISION :: q_c, q_r, q_g, x_c, nu_c, tau, phi, n_c, k_a, x_s, au, nu_r

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD autoconversionSB",4)


    !..Skalar-Initialisierung
    nu_r = 9.59
    x_s  = cloud%x_max                     !..Trennmasse
    x_c  = 4./3. * pi * rho_w * r_c**3     !..Mittlere Masse der Wolkentropfen
    k_a  = 6.0d+25 * nu_r**(-1.7)

    rate_q = 0.0
    rate_n = 0.0
    !..Parametrisierung nach Beheng (1994)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_c = q_cloud(i,j,k)!..Fluessigwassergehalt
          IF (q_c > eps) THEN

            !..Berechnung der Autokonversionsrate nach Beheng (1994)
            au = k_a * (x_c*1e3)**(3.3) * (q_c*1e-3)**(1.4) * dt * 1e3
            au = MIN(q_c,au)

            rate_q(i,j,k) = au
            rate_n(i,j,k) = au / x_s * 2.0

          ENDIF
        END DO
      END DO
    END DO

    !..Zeitintegration
    n_rain  = n_rain  + rate_n
    q_rain  = q_rain  + rate_q
    q_cloud = q_cloud - rate_q
    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,34) = rate_q
    END IF
#endif
    ! ub<<

    IF(isIO() .AND. isdebug) THEN
      IF (firstcall.NE.1) THEN
        WRITE (6,*) "CLOUD autoconversion KB:"
        WRITE (6,'(A,D10.3)') "    Trennmasse Wolken-Regentropfen:   x_s = ", x_s
        WRITE (6,'(A,D10.3)') "    mittlere Masse der Wolkentropfen: x_c = ", x_c
        WRITE (6,'(A,D10.3)') "    Koeff. fuer Autokonversionsrate:  k_a = ", k_a
        firstcall = 1
      ENDIF
      WRITE (6, *) "CLOUD autoconversion KB: cloud_typ == ", cloud_typ 
    END IF

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE autoconversionKB

  SUBROUTINE accretionKB ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, q_cloud, q_rain, n_rain, dt
    USE konstanten,         ONLY: pi
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE globale_variablen,  ONLY: dqdt, speichere_dqdt

    IMPLICIT NONE

    !..Parameter fuer Beheng (1994)
    DOUBLE PRECISION, PARAMETER :: k_r = 6.00d+00   ! Parameter Kernel
    DOUBLE PRECISION, PARAMETER :: eps = 1.00d-25

    !..Lokale Variablen
    INTEGER          :: stat_var = 0
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: ac
    DOUBLE PRECISION :: L_c, L_r, L_g, tau, phi

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q

    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD accretionnSB",4)

    rate_q = 0.0
    !..Parametrisierung nach Seifert und Beheng (2000)
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          L_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          L_r = q_rain(i,j,k)  !..Fluessigwassergehalt

          IF (L_c > 0.0d0.AND.L_r > 0.0d0) THEN
            !..Berechnung der Akkreszenzrate nach Beheng 1994

            ac = k_r *  L_c * L_r * dt

            ac = MIN(L_c,ac)

            rate_q(i,j,k)  = ac
          ENDIF
        END DO
      END DO
    END DO

    !..Zeitintegration
    q_rain  = q_rain  + rate_q
    q_cloud = q_cloud - rate_q
    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,35) = rate_q
    END IF
#endif
    ! ub<<

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD accretion KB"
    END IF

    DEALLOCATE(rate_q)

  END SUBROUTINE accretionKB

  SUBROUTINE autoconversionKK ()
    !*******************************************************************************
    !                                                                              *
    !   Berechnung der Autokonversion von Wolkenwasser zu Niederschlagswasser      *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, q_cloud, q_rain, n_rain, dt
    USE konstanten,         ONLY: pi
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE globale_variablen,  ONLY: dqdt, speichere_dqdt

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: eps  = 1.00d-25
    DOUBLE PRECISION, PARAMETER :: k_a  = 3.47d+07

    !..Locale Variablen
    INTEGER          :: stat_var = 0
    INTEGER          :: i,j,k
    INTEGER, SAVE    :: firstcall
    DOUBLE PRECISION :: q_c, q_r, q_g, x_c, x_s, au

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_n,rate_q

    ALLOCATE(rate_n(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD autoconversionSB",4)


    !..Skalar-Initialisierung
    x_s  = cloud%x_max                     !..Trennmasse
    x_c  = 4./3. * pi * rho_w * r_c**3     !..Mittlere Masse der Wolkentropfen

    rate_q = 0.0
    rate_n = 0.0
    !..Parametrisierung nach  Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          q_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          IF (q_c > eps) THEN

            !..Berechnung der Autokonversionsrate nach KK2000
            au = k_a * (q_c*1e-3)**(0.68) * (x_c*1e3)**(1.79) * dt *1e3
            au = MIN(q_c,au)

            rate_q(i,j,k) = au
            rate_n(i,j,k) = au / x_s * 2.0

          ENDIF
        END DO
      END DO
    END DO

    !..Zeitintegration
    n_rain  = n_rain  + rate_n
    q_rain  = q_rain  + rate_q
    q_cloud = q_cloud - rate_q
    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,34) = rate_q
    END IF
#endif
    ! ub<<

    IF(isIO() .AND. isdebug) THEN
      IF (firstcall.NE.1) THEN
        WRITE (6,*) "CLOUD autoconversion KK:"
        WRITE (6,'(A,D10.3)') "    Trennmasse Wolken-Regentropfen:   x_s = ", x_s
        WRITE (6,'(A,D10.3)') "    mittlere Masse der Wolkentropfen: x_c = ", x_c
        WRITE (6,'(A,D10.3)') "    Koeff. fuer Autokonversionsrate:  k_a = ", k_a
        firstcall = 1
      ENDIF
      WRITE (6, *) "CLOUD autoconversion SB: cloud_typ == ", cloud_typ 
    END IF

    DEALLOCATE(rate_n,rate_q)

  END SUBROUTINE autoconversionKK

  SUBROUTINE accretionKK ()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Akkretion von Wolkentroepfchen an Regentroepfchen       *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, q_cloud, q_rain, n_rain, dt
    USE konstanten,         ONLY: pi
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE globale_variablen,  ONLY: dqdt, speichere_dqdt

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: eps = 1.00d-25
    DOUBLE PRECISION, PARAMETER :: k_a = 5.32d+05
    !DOUBLE PRECISION, PARAMETER :: k_a = 6.70d+01

    !..Lokale Variablen
    INTEGER          :: stat_var = 0
    INTEGER          :: i,j,k
    DOUBLE PRECISION :: ac
    DOUBLE PRECISION :: L_c, L_r, L_g, tau, phi

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rate_q

    ALLOCATE(rate_q(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    IF (stat_var/=0) CALL abortparallel("Allokierungsfehler CLOUD accretionnSB",4)

    rate_q = 0.0
    !..Parametrisierung nach Khairoutdinov and Kogan (2000), MWR 128, 229-243
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          L_c = q_cloud(i,j,k) !..Fluessigwassergehalt
          L_r = q_rain(i,j,k)  !..Fluessigwassergehalt

          IF (L_c > 0.0d0 .AND. L_r > 0.0d0) THEN
            ac  = k_a *  (L_c * L_r * 1e-6)**1.15 * dt * 1e3

            ac = MIN(L_c,ac)

            rate_q(i,j,k)  = ac
          ENDIF
        END DO
      END DO
    END DO

    !..Zeitintegration
    q_rain  = q_rain  + rate_q
    q_cloud = q_cloud - rate_q
    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      dqdt(:,:,:,35) = rate_q
    END IF
#endif
    ! ub<<

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD accretion SB"
    END IF

    DEALLOCATE(rate_q)

  END SUBROUTINE accretionKK

  SUBROUTINE vapor_deposition_growth(dt_local,i_local)
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Wachstums der Eispartikel durch Wasserdampfdiffusion    *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q,T_g,p_g, rho_g,       &
         &                        q_cloud, n_cloud, q_ice, n_ice, q_graupel, n_graupel, &
         &                        q_hail, n_hail,                                       & 
         &                        q_snow, n_snow, q_rain, s_i, s_w, &
         &                        dqdt, speichere_dqdt, cond_neu_sb, &
         &                        evap_neu_sb, speichere_precipstat, &
         &                        deprate_ice, deprate_snow
    USE konstanten,         ONLY: pi,cv,cp
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0,dichte

    IMPLICIT NONE

    INTEGER                     :: i_local
    DOUBLE PRECISION            :: dt_local
    DOUBLE PRECISION            :: a_dl,a_ld,D_vtp
    DOUBLE PRECISION            :: q_g,n_g,x_g,D_g,q_s,n_s,x_s,conv_q,conv_n,x_conv
    ! ub>> 
    DOUBLE PRECISION            :: necessary_d,available_d,weight,depitmp,depstmp,depgtmp,dephtmp,&
         tv_i,tv_s,tv_g,tv_h,tv(4),tvsort(4),depvec(4),depvec2(4),qxvec(4)
    INTEGER :: indsort(4)
    ! ub<<
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: s_si,g_i
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: s_sw,g_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_ice
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_snow
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_cloud
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_graupel, dep_graupel_n
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_hail, dep_hail_n

    ! Locale Variablen 
    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_si            !..Wasserdampfdichte bei Eissaettigung
    DOUBLE PRECISION            :: e_si            !..Wasserpartialdruck bei Eissaettigung
    DOUBLE PRECISION            :: e_sw            !..Wasserpartialdruck bei saettigung
    DOUBLE PRECISION            :: q_d,x_d,e_d,R_f,p_a,rho_a,dep_sum
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: qimax, qsmax, qgmax, qhmax        ! <hn

    ALLOCATE(s_si(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(s_sw(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(g_i(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(g_w(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_ice(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_snow(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_cloud(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_graupel(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_graupel_n(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_hail(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_hail_n(0:loc_ix,1:loc_iy,1:loc_iz))

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD vapor_deposition_growth: " 
    END IF

    IF (firstcall.NE.1) THEN
      IF(isIO() .AND. isdebug)THEN
        WRITE (6,'(A,D10.3)') "  Diffusivitaet von Wasserdampf in Luft:    D_v   = ",D_v
        WRITE (6,'(A,D10.3)') "  Waermeleitfaehigkeit von Luft:            K_T   = ",K_T
        WRITE (6,'(A,D10.3)') "  Sublimationswaerme:                       L_ed  = ",L_ed
        WRITE (6,'(A,D10.3)') "  Kinematische Viskositaet von Luft:        nu_l  = ",nu_l 
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          p_a  = p_0(i,j,k) !WRF!+ p(i,j,k) + p_g(i,j,k)
          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          x_d  = q(i,j,k) / (rho_0(i,j,k))!WRF!+rho_g(i,j,k))
          e_d  = q(i,j,k) * R_d * T_a
          e_si = e_es(T_a)
          e_sw = e_ws(T_a)
          s_si(i,j,k) = e_d / e_si - 1.0                    !..Uebersaettigung bzgl. Eis
          !s_sw(i,j,k) = e_d / e_sw - 1.0                    !..Uebersaettigung bzgl. Fluessigwasser
! UB_20090316:          D_vtp = D_v
          D_vtp = 8.7602e-5 * T_a**(1.81) / p_a
          IF (T_a < T_3) THEN
            g_i(i,j,k) = 4.0*pi / ( L_ed**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_si) )
            ! g_w(i,j,k) = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )
          ELSE
            ! g_w(i,j,k)  = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )
            g_i(i,j,k)  = 0.0
            s_si(i,j,k) = 0.0
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !CALL vapor_deposition_cloud()
    !CALL saturation_adjust_limiter() ! Limitieren der Kondensation durch Saettigungsadjustierung

    dep_ice     = 0.0
    dep_snow    = 0.0
    dep_cloud   = 0.0
    dep_graupel = 0.0
    dep_hail    = 0.0

    ! ub>>
    ! Die Routinen vapor_deposition_ice(), ..., wurden zeitsparender programmiert, indem
    ! Ausdruecke wie "a ** b" durch exp(b*log(a)) ersetzt wurden und einige Berechnungen
    ! mit Konstanten vor die Schleifen gezogen wurde. Insbesondere durch die andere
    ! Potenzierungstechnik ergeben sich aber numerisch leicht unterschiedliche Ergebnisse.
    ! Im Rahmen der hier verwendeten Depositions-Approximation spielt das aber keine
    ! grosse Rolle, weil die Approximation an sich schon recht ungenau ist.
    ! ub<<

    CALL vapor_deposition_ice()
    CALL vapor_deposition_snow()
    CALL vapor_deposition_graupel()
    IF (ice_typ > 1) CALL vapor_deposition_hail()

    !q_cloud = q_cloud + dep_cloud
    !q       = q       - dep_cloud
    !T       = T       + dep_cloud * L_ed / cp / rho_0(i,j,k) 

    x_conv = (250e-6/snow%a_geo)**(1./snow%b_geo)

    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt .AND. i_local == 1) THEN
      ! Da die Deposition u.U. mit einem kleineren Zeitschritt (dt_local) gerechnet wird,
      ! muessen weiter unten die Umwandlungsraten waehrend eines grossen Modellzeitschrittes
      ! (dt) aufsummiert werden. Wenn i_local == 1, dann steht die Zeit auf dem Beginn eines
      ! grossen Zeitschrittes, also dem Anfangszeitpunkt der Summierung. Deshalb:
      ! Zur Sicherheit die Abschnitte des Umwandlungsratenspeicherfeldes dqdt nullen,
      ! die fuer die Depositionsraten gebraucht werden.
      dqdt(:,:,:,4:7) = 0.0
    END IF
#endif
    ! ub<<

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k) 

          ! Nur bei Temperaturen unterhalb des Tripelpunktes wird Deposition/Sublimation parametrisiert.
          ! Bei T > T_3 kann Eis nur dann verdunsten, wenn es vorher anschmilzt.

          IF (T_a < T_3) THEN

            indsort = 0

            IF (s_si(i,j,k) >= 0.0) THEN
              ! uebersaettigt: dep_ice, dep_snow, dep_graupel sind >= 0. Beschraenken auf 
              ! maximal moeglichen uebersaettigten Wasserdampfanteil:
              q_d  = q(i,j,k)

!!! wurde rausgenommen, weil fuer die nachfolgende Reduzierung mittels weight weiter unten
!!! das urspruengliche Verhaeltnis der Depositionsfluesse massgeblich sein soll.
              !                  dep_ice(i,j,k)     = MIN(dep_ice(i,j,k),    q_d)
              !                  dep_snow(i,j,k)    = MIN(dep_snow(i,j,k),   q_d) 
              !                  dep_graupel(i,j,k) = MIN(dep_graupel(i,j,k),q_d)
              !                  dep_hail(i,j,k)    = MIN(dep_hail(i,j,k),q_d) 

              ! Saettigungsueberschuss (zu kondensierendes Wasser, kg/m^3): (Naeherung fuer kleine Dampfdruecke)
              available_d = q_d * s_si(i,j,k) / (s_si(i,j,k) + 1.0)
              ! von der Parametrisierung bestimmte Kondensationsmenge waehrend des Zeitschritts:
              necessary_d = dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k)

              IF (ABS(necessary_d) < eps) THEN
                dep_ice(i,j,k)     = 0.0
                dep_snow(i,j,k)    = 0.0
                dep_graupel(i,j,k) = 0.0               
                dep_hail(i,j,k)    = 0.0
              ELSEIF (necessary_d > available_d) THEN
                weight = available_d / necessary_d
                dep_ice(i,j,k)     = weight * dep_ice(i,j,k)
                dep_snow(i,j,k)    = weight * dep_snow(i,j,k)
                dep_graupel(i,j,k) = weight * dep_graupel(i,j,k)
                dep_hail(i,j,k)    = weight * dep_hail(i,j,k)
              ENDIF

            ELSE
              ! untersaettigt: dep_ice, dep_snow, dep_graupel sind < 0. Beschraenken auf 
              ! maximal moeglichen untersaettigten Wasserdampfanteil sowie auf maximal vorhandenes Eis.
              ! Die Notwendigkeit einer Begrenzung entsteht durch einen i.a. viel zu grossen
              ! Depositionszeitschritt, so dass das hier verwendete 
              ! Euler-Vorwaerts-Zeitintegrationsschema hoffnungslos instabil ist.

              ! Hier wird im Gegensatz zum uebersaettigten Fall zuerst jeder einzelne Depositionsfluss
              ! auf das maximal vorhandene Reservoir begrenzt, bevor durch weight auf den
              ! maximal moeglichen Dampfinput fuer die untersaettigte Luft begrenzt wird.
              depitmp    = MAX(dep_ice(i,j,k),    -q_ice(i,j,k))
              depstmp    = MAX(dep_snow(i,j,k),   -q_snow(i,j,k))
              dephtmp    = MAX(dep_hail(i,j,k),   -q_hail(i,j,k))  
              depgtmp    = MAX(dep_graupel(i,j,k),-q_graupel(i,j,k)) 

              ! Saettigungsdefizit (Wasserdampf-Aufnahmereservoir der Luft, kg/m^3):
              available_d = s_si(i,j,k) * e_es(T_0(i,j,k)) / (R_d * T_0(i,j,k))
              ! von der Parametrisierung bestimmte Verdunstungsmenge waehrend des Zeitschritts:
              necessary_d = depitmp + depstmp + depgtmp + dephtmp

              ! available_d und necessary_d sind hier < 0!
              IF (necessary_d < -eps) THEN
                IF (necessary_d < available_d) THEN
                  ! In diesem Falle kann waehrend des Zeitschritts nicht soviel verdunstet werden 
                  ! wie von der Parametrisierung berechnet, weil das Saettigungsdefizit nicht gross genug ist. 
                  ! Die Verdunstung wird auf das Saettigungsdefizit begrenzt, und die Begrenzung
                  ! wird gemaess der einzelnen Verdunstungsraten "gerecht" auf die Niederschlagsarten aufgeteilt.
                  ! Dabei wird beruecksichtigt, dass waehrend des Zeitschritts von mindestens einer Niederschlagsart
                  ! ggf. alles Wasser verdunsten kann. Dann wird die Begrenzung auf die verbleibende(n)
                  ! Niederschlagsart(en) aufgeteilt, unter Beruecksichtigung, dass evtl. dafuer nicht
                  ! genuegend Reservoir (Saettigungsdefizit) in der Umgebungsluft vorhanden ist.

                  ! Dazu werden benoetigt: Verdunstungszeitskalen der einzelnen Hydrometeorarten:

                  tv_i = dt_local / eps ! initialis. mit seeeehr langem Zeitscale ....
                  tv_s = tv_i
                  tv_g = tv_i  
                  tv_h = tv_i 

                  IF (ABS(dep_ice(i,j,k))     > eps) tv_i = -q_ice(i,j,k) / dep_ice(i,j,k) * dt_local
                  IF (ABS(dep_snow(i,j,k))    > eps) tv_s = -q_snow(i,j,k) / dep_snow(i,j,k) * dt_local
                  IF (ABS(dep_graupel(i,j,k)) > eps) tv_g = -q_graupel(i,j,k) / dep_graupel(i,j,k) * dt_local
                  IF (ABS(dep_hail(i,j,k))    > eps) tv_h = -q_hail(i,j,k) / dep_hail(i,j,k) * dt_local 


                  IF (tv_i >= dt_local .AND. tv_s >= dt_local .AND. tv_g >= dt_local .AND. tv_h >= dt_local) THEN
                    ! Keine Niederschlagsart wuerde waehrend des Zeitschritts vollstaendig aufgezehrt.
                    ! Der Begrenzungsfaktor ist folglich bei allen gleich.
                    weight = available_d / necessary_d
                    dep_ice(i,j,k)     = weight * depitmp
                    dep_snow(i,j,k)    = weight * depstmp
                    dep_graupel(i,j,k) = weight * depgtmp
                    dep_hail(i,j,k)    = weight * dephtmp

                  ELSE 
                    ! Mindestens eine Niederschlagsart wird waehrend des Zeitschrittes vollstaendig aufgezehrt.
                    ! Zusammen mit der Tatsache, dass die Verdunstung aufgrund zu kleinen
                    ! Saettigungsdefizits begrenzt werden muss, ergibt sich das Problem der "sinnvollen"
                    ! Gewichtung dieser Begrenzung. Hier wird der Ansatz verfolgt, dass Hydrometeorarten,
                    ! die waehrend des Zeitschritts aufgrund der Parametrisierung vollstaendig verdunsten
                    ! wuerden, dies auch tun, wenn fuer den in der entsprechenden Zeit insgesamt freigesetzten 
                    ! Wasserdampf genuegend Reservoir
                    ! in der Umgebungsluft vorhanden ist. Die Begrenzung wird dann auf die verbleibende(n) Hydrometeorart(en)
                    ! "gerecht" anhand der Verdunstungsraten aufgeteilt.

                    ! Man erreicht dies z.B. durch eine Aufteilung in Zeitbereiche anhand der
                    ! Verdunstungszeiten der einzelnen Hydrometeorarten, welche man vorher nach ebendieser
                    ! Zeit sortiert. 

                    tv = (/tv_i, tv_s, tv_g, tv_h/)
                    depvec = (/dep_ice(i,j,k), dep_snow(i,j,k), dep_graupel(i,j,k), dep_hail(i,j,k)/)
                    depvec2 = (/depitmp, depstmp, depgtmp, dephtmp/)
                    qxvec(1) = -q_ice(i,j,k)
                    qxvec(2) = -q_snow(i,j,k)
                    qxvec(3) = -q_graupel(i,j,k)
                    qxvec(4) = -q_hail(i,j,k)

                    IF (ice_typ > 1) THEN
                      ! Rechnung mit Hagel,
                      ! aus Rechenzeitgruenden hier Fallunterscheidung

                      tvsort = sortiere4(tv, indsort)

                      IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) >= dt_local .AND. &
                           tvsort(3) >= dt_local .AND. &
                           tvsort(4) >= dt_local) THEN
                        ! Die am schnellsten verdunstende Niederschlagsart wuerde waehrend dt_local  aufgezehrt, die beiden anderen nicht.
                        IF (available_d < (depvec2(indsort(1)) + &
                             (depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4)))/dt_local*tv(indsort(1)))) THEN
                          ! Es ist genuegend Feuchtedefizit in der Luft vorhanden, damit die am schnellsten verdunstende
                          ! Niederschlagsart vollstaendig verdunsten kann. Begrenzt werden die Raten der anderen beiden Arten.
                          depvec(indsort(1))     = depvec2(indsort(1))
                          IF (ABS(depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4))) >= eps) THEN
                            weight = (available_d-depvec2(indsort(1))) &
                                   & / (depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4)))
                          ELSE
                            weight = 0.0
                          END IF
                          depvec(indsort(2)) = weight * depvec(indsort(2))
                          depvec(indsort(3)) = weight * depvec(indsort(3))
                          depvec(indsort(4)) = weight * depvec(indsort(4))

                        ELSE

                          ! Es ist nicht genuegend Feuchtedefizit in der Luft vorhanden, damit die am schnellsten verdunstende
                          ! Niederschlagsart vollstaendig verdunsten kann. Begrenzt werden alle 4 Arten:
                          weight = available_d &
                            & / (depvec(indsort(1))+depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4)))
                          depvec(indsort(1)) = weight * depvec(indsort(1))
                          depvec(indsort(2)) = weight * depvec(indsort(2))
                          depvec(indsort(3)) = weight * depvec(indsort(3))
                          depvec(indsort(4)) = weight * depvec(indsort(4))

                        END IF

                      ELSE IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) < dt_local .AND. &
                           tvsort(3) >= dt_local .AND. &
                           tvsort(4) >= dt_local) THEN
                        ! Die zwei am schnellsten verdunstenden Niederschlagsarten wuerden waehrend dt_local 
                        ! aufgezehrt, die dritte und vierte nicht.
                        IF (available_d < depvec2(indsort(1))+depvec2(indsort(2)) + &
                             (depvec(indsort(3))+depvec(indsort(4)))/dt_local*tv(indsort(2))) THEN
                          ! Genuegend Feuchte-Reservoir in der Umgebungsluft, deshalb koennen die am schnellsten und zweitschnellsten
                          ! verdunstenden Niederschlagsarten auch wirklich verdunsten.
                          depvec(indsort(1))     = depvec2(indsort(1))
                          depvec(indsort(2))     = depvec2(indsort(2))
                          IF (ABS(depvec(indsort(3))+depvec(indsort(4))) >= eps) THEN
                            weight = (available_d-depvec2(indsort(1))-depvec2(indsort(2))) &
                                   & / (depvec(indsort(3))+depvec(indsort(4)))
                          ELSE
                            weight = 0.0
                          END IF
                          depvec(indsort(3)) = weight * depvec(indsort(3))
                          depvec(indsort(4)) = weight * depvec(indsort(4))

                        ELSE

                          IF (available_d < depvec2(indsort(1))+(depvec(indsort(2))+depvec(indsort(3))+ &
                               depvec(indsort(4)))/dt_local*tv(indsort(1))) THEN
                            ! Nur die am schnellsten verdunstende Niederschlagsart wird aufgezehrt; 
                            ! fuer die zweite ist die Zeit bis zum Fuellen des Reservoirs nicht mehr lang genug. 
                            ! Die zweite, dritte und vierte werden nicht aufgezehrt.
                            depvec(indsort(1)) = depvec2(indsort(1))
                            IF (ABS(depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4))) >= eps) THEN
                              weight = (available_d-depvec2(indsort(1))) &
                                     & / (depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4)))
                            ELSE
                              weight = 0.0
                            END IF
                            depvec(indsort(2)) = weight * depvec(indsort(2))
                            depvec(indsort(3)) = weight * depvec(indsort(3))
                            depvec(indsort(4)) = weight * depvec(indsort(4))

                          ELSE
                            ! Aufgrund zu kleinen Reservoirs kann gar keine Niederschlagsart vollstaendig verdunsten;
                            ! Die Zeit bis zum Auffuellen des Reservoirs ist hierfuer zu kurz.
                            weight = available_d / (SUM(depvec))
                            depvec(indsort(1)) = weight * depvec(indsort(1))
                            depvec(indsort(2)) = weight * depvec(indsort(2))
                            depvec(indsort(3)) = weight * depvec(indsort(3))
                            depvec(indsort(4)) = weight * depvec(indsort(4))

                          END IF
                        END IF

                      ELSE IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) < dt_local .AND. &
                           tvsort(3) < dt_local .AND. &
                           tvsort(4) >= dt_local) THEN
                        ! Die drei am schnellsten verdunstenden Niederschlagsarten wuerden waehrend 
                        ! dt_local aufgezehrt werden, die vierte nicht.
                        IF (available_d < depvec2(indsort(1))+depvec2(indsort(2))+depvec2(indsort(3)) + &
                             depvec(indsort(4))/dt_local*tv(indsort(3))) THEN
                          ! Genuegend Feuchte-Reservoir in der Umgebungsluft, deshalb koennen die am schnellsten, zweit- und drittschnellsten
                          ! verdunstenden Niederschlagsarten auch wirklich verdunsten.
                          depvec(indsort(1)) = depvec2(indsort(1))
                          depvec(indsort(2)) = depvec2(indsort(2))
                          depvec(indsort(3)) = depvec2(indsort(3))
                          depvec(indsort(4)) = available_d - (depvec2(indsort(1))+depvec2(indsort(2))+depvec2(indsort(3)))
                        ELSE

                          IF (available_d < depvec2(indsort(1))+depvec2(indsort(2))+ &
                               (depvec(indsort(3))+depvec(indsort(4)))/dt_local*tv(indsort(2))) THEN
                            ! Nur die zwei am schnellsten verdunstenden Niederschlagsarten werden aufgezehrt; 
                            ! fuer die dritte ist die Zeit bis zum Fuellen des Reservoirs nicht mehr lang genug. 
                            ! Die dritte und vierte wird nicht aufgezehrt.
                            depvec(indsort(1)) = depvec2(indsort(1))
                            depvec(indsort(2)) = depvec2(indsort(2))
                            IF (ABS(depvec(indsort(3))+depvec(indsort(4))) >= eps) THEN
                              weight = (available_d-depvec2(indsort(1))-depvec2(indsort(2)))  &
                                     & / (depvec(indsort(3))+depvec(indsort(4)))
                            ELSE
                              weight = 0.0
                            END IF
                            depvec(indsort(3)) = weight * depvec(indsort(3))
                            depvec(indsort(4)) = weight * depvec(indsort(4))

                          ELSE

                            IF (available_d < depvec2(indsort(1))+(depvec(indsort(2))+depvec(indsort(3))+ &
                                 depvec(indsort(4)))/dt_local*tv(indsort(1))) THEN
                              ! Nur die am schnellsten verdunstende Niederschlagsart wird aufgezehrt; 
                              ! fuer die zweite und dritte ist die Zeit bis zum Fuellen des Reservoirs nicht mehr lang genug. 
                              ! Die zweite, dritte und vierte werden nicht aufgezehrt.
                              depvec(indsort(1)) = depvec2(indsort(1))
                              IF (ABS(depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4))) >= eps) THEN
                                weight = (available_d-depvec2(indsort(1))) &
                                       & / (depvec(indsort(2))+depvec(indsort(3))+depvec(indsort(4)))
                              ELSE
                                weight = 0.0
                              END IF
                              depvec(indsort(2)) = weight * depvec(indsort(2))
                              depvec(indsort(3)) = weight * depvec(indsort(3))
                              depvec(indsort(4)) = weight * depvec(indsort(4))

                            ELSE

                              ! Aufgrund zu kleinen Reservoirs kann gar keine Niederschlagsart vollstaendig verdunsten;
                              ! Die Zeit bis zum Auffuellen des Reservoirs ist hierfuer zu kurz.
                              weight = available_d / (SUM(depvec))
                              depvec(indsort(1)) = weight * depvec(indsort(1))
                              depvec(indsort(2)) = weight * depvec(indsort(2))
                              depvec(indsort(3)) = weight * depvec(indsort(3))
                              depvec(indsort(4)) = weight * depvec(indsort(4))

                            END IF
                          END IF
                        END IF

                      ELSE IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) < dt_local .AND. &
                           tvsort(3) < dt_local .AND. &
                           tvsort(4) < dt_local) THEN
                        ! der Fall dass necessary < available und trotzdem alle tv_x < dt_local
                        ! kann eintreten, da necessary aus depxtmp berechnet wird, tv_x aber aus dep_x.
                        ! Alle Eisteilchen verdunsten vollstaendig
                        depvec(indsort(1)) = depvec2(indsort(1))   
                        depvec(indsort(2)) = depvec2(indsort(2))   
                        depvec(indsort(3)) = depvec2(indsort(3))   
                        depvec(indsort(4)) = depvec2(indsort(4))                                         

                      ELSE 
                        WRITE (*,*) 'VAPOR_DEPOSITION: Da ist was faul im untersaettigten Fall excl 1'
                        WRITE (*,*) 'tvsort =', tvsort
                        WRITE (*,*) 'tv     = ', tv
                        WRITE (*,*) 'depvec = ', depvec
                        WRITE (*,*) 'depvec2= ', depvec2
                        WRITE (*,*) 'qxvec  = ', qxvec
                        WRITE (*,*) 'dt_local= ', dt_local, '  available_d=', available_d, ' necesary_d=',necessary_d
                      END IF

                    ELSE ! ice_typ <= 1, kein Hagel, vereinfachter Entscheidungsbaum:

                      tvsort(1:3) = sortiere3(tv(1:3), indsort(1:3))

                      IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) >= dt_local .AND. &
                           tvsort(3) >= dt_local) THEN
                        ! Die am schnellsten verdunstende Niederschlagsart wuerde waehrend dt_local 
                        ! aufgezehrt, die beiden anderen nicht.
                        IF (available_d < (depvec2(indsort(1)) + &
                             (depvec(indsort(2))+depvec(indsort(3)))/dt_local*tv(indsort(1)))) THEN
                          ! Es ist genuegend Feuchtedefizit in der Luft vorhanden, damit die am schnellsten verdunstende
                          ! Niederschlagsart vollstaendig verdunsten kann. Begrenzt werden die Raten der anderen beiden Arten.
                          depvec(indsort(1))     = depvec2(indsort(1))
                          IF (ABS(depvec(indsort(2))+depvec(indsort(3))) >= eps) THEN
                            weight = (available_d-depvec2(indsort(1))) / (depvec(indsort(2))+depvec(indsort(3)))
                          ELSE
                            weight = 0.0
                          END IF
                          depvec(indsort(2)) = weight * depvec(indsort(2))
                          depvec(indsort(3)) = weight * depvec(indsort(3))

                        ELSE

                          ! Es ist nicht genuegend Feuchtedefizit in der Luft vorhanden, damit die am schnellsten verdunstende
                          ! Niederschlagsart vollstaendig verdunsten kann. Begrenzt werden alle 3 Arten:
                          weight = available_d / (depvec(indsort(1))+depvec(indsort(2))+depvec(indsort(3)))
                          depvec(indsort(1)) = weight * depvec(indsort(1))
                          depvec(indsort(2)) = weight * depvec(indsort(2))
                          depvec(indsort(3)) = weight * depvec(indsort(3))

                        END IF

                      ELSE IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) < dt_local .AND. &
                           tvsort(3) >= dt_local) THEN
                        ! Die zwei am schnellsten verdunstende Niederschlagsarten wuerden waehrend dt_local 
                        ! aufgezehrt, die dritte nicht.

                        IF (available_d < depvec2(indsort(1))+depvec2(indsort(2)) + &
                             (depvec(indsort(3)))/dt_local*tv(indsort(2))) THEN
                          ! Es ist genuegend Feuchtedefizit in der Luft vorhanden, damit die zwei am schnellsten verdunstenden
                          ! Niederschlagsarten vollstaendig verdunsten koennen. Begrenzt wird die Rate der dritten NS-Art.
                          depvec(indsort(1))     = depvec2(indsort(1))
                          depvec(indsort(2))     = depvec2(indsort(2))
                          depvec(indsort(3)) = available_d - depvec2(indsort(1)) - depvec2(indsort(2))

                        ELSE

                          IF (available_d < depvec2(indsort(1))+(depvec(indsort(2))+depvec(indsort(3))+ &
                               depvec(indsort(4)))/dt_local*tv(indsort(1))) THEN
                            ! Nur die am schnellsten verdunstende Niederschlagsart wird aufgezehrt; 
                            ! fuer die zweite ist die Zeit bis zum Fuellen des Reservoirs nicht mehr lang genug. 
                            ! Die zweite und dritte werden nicht aufgezehrt.
                            depvec(indsort(1)) = depvec2(indsort(1))
                            IF (ABS(depvec(indsort(2))+depvec(indsort(3))) >= eps) THEN
                              weight = (available_d-depvec2(indsort(1))) / (depvec(indsort(2))+depvec(indsort(3)))
                            ELSE
                              weight = 0.0
                            END IF
                            depvec(indsort(2)) = weight * depvec(indsort(2))
                            depvec(indsort(3)) = weight * depvec(indsort(3))

                          ELSE
                            ! Aufgrund zu kleinen Reservoirs kann gar keine Niederschlagsart vollstaendig verdunsten;
                            ! Die Zeit bis zum Auffuellen des Reservoirs ist hierfuer zu kurz.
                            weight = available_d / (SUM(depvec(1:3)))
                            depvec(indsort(1)) = weight * depvec(indsort(1))
                            depvec(indsort(2)) = weight * depvec(indsort(2))
                            depvec(indsort(3)) = weight * depvec(indsort(3))

                          END IF
                        END IF

                      ELSE IF (tvsort(1) < dt_local .AND. &
                           tvsort(2) < dt_local .AND. &
                           tvsort(3) < dt_local) THEN
                        ! der Fall dass necessary < available und trotzdem alle tv_x < dt_local
                        ! kann eintreten, da necessary aus depxtmp berechnet wird, tv_x aber aus dep_x.
                        ! Alle Eisteilchen verdunsten vollstaendig
                        depvec(indsort(1)) = depvec2(indsort(1))   
                        depvec(indsort(2)) = depvec2(indsort(2))   
                        depvec(indsort(3)) = depvec2(indsort(3))   

                      ELSE 
                        WRITE (*,*) 'VAPOR_DEPOSITION: Da ist was faul im untersaettigten Fall excl 2'
                        WRITE (*,*) 'tvsort =', tvsort
                        WRITE (*,*) 'tv     = ', tv
                        WRITE (*,*) 'depvec = ', depvec
                        WRITE (*,*) 'depvec2= ', depvec2
                        WRITE (*,*) 'qxvec  = ', qxvec
                        WRITE (*,*) 'dt_local= ', dt_local, '  available_d=', available_d, ' necesary_d=',necessary_d
                      END IF

                    END IF

                    ! Zurueckverteilen:
                    dep_ice(i,j,k) = depvec(1)
                    dep_snow(i,j,k) = depvec(2)
                    dep_graupel(i,j,k) = depvec(3)
                    dep_hail(i,j,k) = depvec(4)

                  END IF


                ELSE

                  ! In diesem Falle kann die vollstaendige von der Parametrisierung berechnete
                  ! Verdunstungsmenge aller Hydrometeorarten verdunstet werden. Keine Begrenzung ist noetig 
                  ! (bis auf die Begrenzung auf die vorhandene Hydrometeormasse).
                  dep_ice(i,j,k)     = depitmp
                  dep_snow(i,j,k)    = depstmp
                  dep_graupel(i,j,k) = depgtmp
                  dep_hail(i,j,k)    = dephtmp

                END IF

              ELSE

                ! Deposition so gering, dass man sie getrost zu 0 setzen kann:
                dep_ice(i,j,k)     = 0.0
                dep_snow(i,j,k)    = 0.0
                dep_graupel(i,j,k) = 0.0    
                dep_hail(i,j,k)    = 0.0                    

              END IF

            END IF

            ! Zur Sicherheit: Pruefen auf Qx < 0.0:
            ! An dieser Stelle kann es vorkommen, dass ein |dep_xxx| > q_xxx ist, und zwar im Falle sehr kleiner absoluter Werte.
            IF (q_ice(i,j,k)+dep_ice(i,j,k) < 0.0 .OR. &
                 q_snow(i,j,k)+dep_snow(i,j,k) < 0.0 .OR. &
                 q_graupel(i,j,k)+dep_graupel(i,j,k) < 0.0 .OR. &
                 q_hail(i,j,k)+dep_hail(i,j,k) < 0.0) THEN

              IF (q_ice(i,j,k)+dep_ice(i,j,k) < -1d-20 .OR. &
                   q_snow(i,j,k)+dep_snow(i,j,k) < -1d-20 .OR. &
                   q_graupel(i,j,k)+dep_graupel(i,j,k) < -1d-20 .OR. &
                   q_hail(i,j,k)+dep_hail(i,j,k) < -1d-20) THEN

                WRITE (*,*) 'VAPOR_DEPOSITION: Warnung --- Qx < 0.0 aufgetreten (wird spaeter begrenzt):'
                WRITE (*,*) '     ----- Kann vorkommen bei sehr kleinen Werten von Qx und dep_x -----'
                WRITE (*,*) 'Punkt ', i, ',', j, ',', k
                WRITE (*,*) 'q_ice     = ', q_ice(i,j,k)+dep_ice(i,j,k), '  dep_ice = ', dep_ice(i,j,k), '  tv_i = ', tv_i
                WRITE (*,*) 'q_snow    = ', q_snow(i,j,k)+dep_snow(i,j,k), '  dep_snow = ', dep_snow(i,j,k), '  tv_s = ', tv_s
                WRITE (*,*) 'q_graupel = ', q_graupel(i,j,k)+dep_graupel(i,j,k), '  dep_graupel = ', dep_graupel(i,j,k), '  tv_g = ', tv_g
                WRITE (*,*) 'q_hail    = ', q_hail(i,j,k)+dep_hail(i,j,k), '  dep_hail = ', dep_hail(i,j,k), '  tv_h = ', tv_h
                WRITE (*,*) 'tv_x_max = dt_local / eps = ', dt_local / eps
              END IF
              !AS 20060220>
              ! Deshalb hier nochmals eine letzte Begrenzung:
              dep_hail(i,j,k)    = MAX(dep_hail(i,j,k),   -q_hail(i,j,k)   )
              dep_graupel(i,j,k) = MAX(dep_graupel(i,j,k),-q_graupel(i,j,k))
              dep_snow(i,j,k)    = MAX(dep_snow(i,j,k),   -q_snow(i,j,k)   )
              dep_ice(i,j,k)     = MAX(dep_ice(i,j,k),    -q_ice(i,j,k)    )
              !<end AS
            END IF

            dep_sum = dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k)
            IF (dep_sum > q(i,j,k) .AND. dep_sum > 0.0) THEN
              IF (dep_sum > 1d-20) THEN
                WRITE (*,*) 'VAPOR_DEPOSITION: Warnung --- dep_sum > q aufgetreten (wird spaeter begrenzt):'
                WRITE (*,*) '     ----- Sollte eigentlich nicht vorkommen -----'
                WRITE (*,*) 'Punkt ', i, ',', j, ',', k
                WRITE (*,*) 'q_ice     = ', q_ice(i,j,k)+dep_ice(i,j,k), '  dep_ice = ', dep_ice(i,j,k)
                WRITE (*,*) 'q_snow    = ', q_snow(i,j,k)+dep_snow(i,j,k), '  dep_snow = ', dep_snow(i,j,k)
                WRITE (*,*) 'q_graupel = ', q_graupel(i,j,k)+dep_graupel(i,j,k), '  dep_graupel = ', dep_graupel(i,j,k)
                WRITE (*,*) 'q_hail    = ', q_hail(i,j,k)+dep_hail(i,j,k), '  dep_hail = ', dep_hail(i,j,k)
                WRITE (*,*) 'q = ', q(i,j,k)+dep_sum, '   dep_sum = ', dep_sum
              END IF
              ! Kondensation begrenzen auf max. vorhandene Feuchtigkeit 
              ! (sollte aber nach dem vorhergegangenen Begrenzen wirklich nicht noetig sein ...):
              weight = q(i,j,k) / dep_sum
              dep_ice(i,j,k)     = weight * dep_ice(i,j,k)
              dep_snow(i,j,k)    = weight * dep_snow(i,j,k)
              dep_graupel(i,j,k) = weight * dep_graupel(i,j,k)
              dep_hail(i,j,k)    = weight * dep_hail(i,j,k)
              dep_sum = q(i,j,k)
            END IF

            ! Moeglichkeit der hail-to-snow conversion momentan nicht vorgesehen
            n_g = n_graupel(i,j,k)
            q_g = q_graupel(i,j,k)                  
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max)
            D_g = graupel%a_geo * x_g**graupel%b_geo   
            IF (.FALSE. .AND. dep_graupel(i,j,k) > 0 .AND. D_g < 800d-6) THEN
              !if (dep_graupel(i,j,k) > 0 .AND. D_g < 800d-6 .AND. q_cloud(i,j,k) < 1e-4) then
              !..Graupel to snow conversion
              conv_q = MIN(q_g, 3.0 * dep_graupel(i,j,k))
              !conv_n = conv_q/x_g
              !conv_n = MIN(n_g, 0.25 * dep_graupel_n(i,j,k))
              q_snow(i,j,k)    = q_snow(i,j,k)    + conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) - conv_q + dep_graupel(i,j,k)
              n_snow(i,j,k)    = n_snow(i,j,k)    + conv_q / MIN(x_g,x_conv)
              n_graupel(i,j,k) = n_graupel(i,j,k) - conv_q / x_g      
            ELSE
              q_graupel(i,j,k) = q_graupel(i,j,k) + dep_graupel(i,j,k)
              conv_q = 0.0
            ENDIF

            q_ice(i,j,k)     = q_ice(i,j,k)     + dep_ice(i,j,k)
            q_snow(i,j,k)    = q_snow(i,j,k)    + dep_snow(i,j,k)
            IF (ice_typ > 1) q_hail(i,j,k)    = q_hail(i,j,k)    + dep_hail(i,j,k)    ! <hn

            q(i,j,k) = q(i,j,k) - dep_sum


            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              ! Wegen kleinerem Zeitschritt dt_local: Aufsummieren, nicht einfach den Wert setzen.
              IF (dep_sum >= 0.0) THEN
                dqdt(i,j,k,4) = dqdt(i,j,k,4) + dep_ice(i,j,k)
                dqdt(i,j,k,5) = dqdt(i,j,k,5) + dep_snow(i,j,k)
                dqdt(i,j,k,6) = dqdt(i,j,k,6) + dep_graupel(i,j,k)
                dqdt(i,j,k,55) = dqdt(i,j,k,55) + dep_hail(i,j,k)
                dqdt(i,j,k,7) = dqdt(i,j,k,7) + conv_q
              ELSE
                dqdt(i,j,k,56) = dqdt(i,j,k,56) - dep_ice(i,j,k)
                dqdt(i,j,k,32) = dqdt(i,j,k,32) - dep_snow(i,j,k)
                dqdt(i,j,k,33) = dqdt(i,j,k,33) - dep_graupel(i,j,k)
                dqdt(i,j,k,53) = dqdt(i,j,k,53) - dep_hail(i,j,k)
              END IF
            END IF
#endif

            IF (speichere_precipstat) THEN
              IF (dep_sum >= 0.0) THEN
                cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + dep_sum
              ELSE
                evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) - dep_sum
              END IF
            END IF
            ! ub<<

            IF (use_ice_graupel_conv_uli) THEN
              deprate_ice(i,j,k) = deprate_ice(i,j,k) + dep_ice(i,j,k)
              deprate_snow(i,j,k) = deprate_snow(i,j,k) + dep_snow(i,j,k)
            END IF

            !WRF!T(i,j,k) = T(i,j,k) + dep_sum * L_ed / cp / rho_0(i,j,k) 
            !WRF!p(i,j,k) = p(i,j,k) + dep_sum * L_ed / cv * R_l


          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(s_si,s_sw,g_i,g_w,dep_ice,dep_cloud,dep_snow,dep_graupel,dep_graupel_n, &
         &  dep_hail,dep_hail_n)

  CONTAINS

    ! ub>>
    FUNCTION sortiere4(x, indsort)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x(4)
      DOUBLE PRECISION :: sortiere4(4)
      INTEGER, INTENT(out) :: indsort(4)
      INTEGER :: i, j, itmp
      DOUBLE PRECISION :: tmp

      sortiere4 = x

      DO i=1,4
        indsort(i) = i
      END DO

      ! Bubble-Sort, Ergebnis in aufsteigender Reihenfolge
      DO i=1,3
        DO j=1,3
          IF (sortiere4(j) > sortiere4(j+1)) THEN
            tmp = sortiere4(j)
            sortiere4(j) = sortiere4(j+1)
            sortiere4(j+1) = tmp
            itmp = indsort(j)
            indsort(j) = indsort(j+1)
            indsort(j+1) = itmp
          END IF
        END DO
      END DO

      !      WRITE(*,*) 'VAPOR_DEPOSITION_GROWTH: Sortiere4 ...'

    END FUNCTION sortiere4

    FUNCTION sortiere3(x, indsort)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x(3)
      DOUBLE PRECISION :: sortiere3(3)
      INTEGER, INTENT(out) :: indsort(3)
      INTEGER :: i, j, itmp
      DOUBLE PRECISION :: tmp

      sortiere3 = x

      DO i=1,3
        indsort(i) = i
      END DO

      ! Bubble-Sort, Ergebnis in aufsteigender Reihenfolge
      DO i=1,2
        DO j=1,2
          IF (sortiere3(j) > sortiere3(j+1)) THEN
            tmp = sortiere3(j)
            sortiere3(j) = sortiere3(j+1)
            sortiere3(j+1) = tmp
            itmp = indsort(j)
            indsort(j) = indsort(j+1)
            indsort(j+1) = itmp
          END IF
        END DO
      END DO

      !      WRITE(*,*) 'VAPOR_DEPOSITION_GROWTH: Sortiere3 ...'

    END FUNCTION sortiere3
    ! ub<<

    SUBROUTINE vapor_deposition_cloud()
      IMPLICIT NONE

      DOUBLE PRECISION            :: q_c,n_c,x_c,d_c,v_c,f_v,N_re,f_v_fakt
      DOUBLE PRECISION, SAVE      :: c_c             !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall

      IF (firstcall.NE.1) THEN
        c_c = 1.0 / cloud%cap
        a_f = vent_coeff_a(cloud,1)
        b_f = vent_coeff_b(cloud,1)
        IF(isIO())THEN
          WRITE (6, *) "  CLOUD vapor_deposition_cloud: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentropfen:  a_cloud   = ",cloud%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentropfen:  b_cloud   = ",cloud%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation:                alf_cloud = ",cloud%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation:                bet_cloud = ",cloud%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_c       = ",c_c 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f       = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f       = ",b_f 
        END IF
        firstcall = 1
      ELSEIF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "  CLOUD vapor_deposition_cloud " 
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            n_c = n_cloud(i,j,k)                                   !..Anzahldichte
            q_c = q_cloud(i,j,k)                                   !..Massendichte

            x_c = MIN(MAX(q_c/(n_c+eps),MAX(cloud%x_min,eps)),cloud%x_max)  !..mittlere Masse
            ! x_c**b_geo durch EXP(b_geo*LOG(x_c)) ersetzen, da x_c == 0 nicht mehr vorkommen kann; 20 % schneller:
            d_c = cloud%a_geo * EXP(cloud%b_geo*LOG(x_c))          !..mittlerer Durchmesser
            v_c = cloud%a_vel * EXP(cloud%b_vel*LOG(x_c)) * rrho_c(i,j,k)   !..mittlere Sedimentationsgeschw.

            N_re = v_c * d_c / nu_l                                !..mittlere Reynoldszahl
            ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
            f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))             !..mittlerer Vent.Koeff.
            f_v  = MAX(f_v,1.d0) !unnoetig??

            dep_cloud(i,j,k) = g_w(i,j,k) * n_c * c_c * d_c * f_v * s_sw(i,j,k) * dt_local
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_cloud

    SUBROUTINE vapor_deposition_ice()
      IMPLICIT NONE

      DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i,f_v,N_re,f_v_fakt
      DOUBLE PRECISION, SAVE      :: c_i             !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall

      IF (firstcall.NE.1) THEN
        c_i = 1.0 / ice%cap
        a_f = vent_coeff_a(ice,1)
        b_f = vent_coeff_b(ice,1)
        IF(isIO() .AND. isdebug)THEN
          WRITE (6, *) "  CLOUD vapor_deposition_ice: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:            a_ice   = ",ice%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:            b_ice   = ",ice%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation von Eis:        alf_ice = ",ice%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation von Eis:        bet_ice = ",ice%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_i     = ",c_i 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f     = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f     = ",b_f 
        END IF
        firstcall = 1
      ELSEIF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "  CLOUD vapor_deposition_ice " 
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_ice=0, dep_ice has to be zero too

            IF (q_ice(i,j,k) == 0.0d0) THEN
              dep_ice(i,j,k)   = 0.0d0
            ELSE  
              n_i = n_ice(i,j,k)                                   !..Anzahldichte
              q_i = q_ice(i,j,k)                                   !..Massendichte

              x_i = MIN(MAX(q_i/(n_i+eps),MAX(ice%x_min,eps)),ice%x_max)    !..mittlere Masse
              ! x_i**b_geo durch EXP(b_geo*LOG(x_i)) ersetzen, da x_i == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_i = ice%a_geo * EXP(ice%b_geo*LOG(x_i))            !..mittlerer Durchmesser
              v_i = ice%a_vel * EXP(ice%b_vel*LOG(x_i)) * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

              N_re = v_i * d_i / nu_l                              !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))     !..mittlerer Vent.Koeff.
              f_v  = MAX(f_v,1.d0) !unnoetig??

              dep_ice(i,j,k) = g_i(i,j,k) * n_i * c_i * d_i * f_v * s_si(i,j,k) * dt_local
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_ice

    SUBROUTINE vapor_deposition_graupel()
      IMPLICIT NONE

      ! Locale Variablen 
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall
      DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g,f_v,f_n,N_re,f_v_fakt,vent_fakt
      DOUBLE PRECISION, SAVE      :: c_g                 !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f,a_n,b_n     !..Koeff. fuer mittleren Ventilationkoeff.

      IF (firstcall.NE.1) THEN
        c_g = 1.0 / graupel%cap
        a_n = vent_coeff_a(graupel,0)
        b_n = vent_coeff_b(graupel,0)
        a_f = vent_coeff_a(graupel,1)
        b_f = vent_coeff_b(graupel,1)
        IF (isIO() .AND. isdebug) THEN
          WRITE (6, *) "  CLOUD vapor_deposition_graupel: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:        a_geo = ",graupel%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:        b_geo = ",graupel%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Graupel:        a_vel = ",graupel%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Graupel:        b_vel = ",graupel%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Gaskonstante feuchter Luft:   a_dl  = ",a_dl 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_g   = ",c_g 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f   = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f   = ",b_f 
        END IF
        firstcall = 1
      ELSEIF(isIO() .AND. isdebug)THEN
        WRITE (6, *) "  CLOUD vapor_deposition_graupel" 
      ENDIF

      f_v_fakt = N_sc**n_f
      vent_fakt = b_n / b_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_garupel=0, dep_graupel has to be zero too

            IF (q_graupel(i,j,k) == 0.0d0) THEN
              dep_graupel(i,j,k)   = 0.0d0
              dep_graupel_n(i,j,k) = 0.0d0
            ELSE

              n_g = n_graupel(i,j,k)                                     !..Anzahldichte
              q_g = q_graupel(i,j,k)                                     !..Massendichte

              x_g = MIN(MAX(q_g/(n_g+eps),MAX(graupel%x_min,eps)),graupel%x_max)  !..mittlere Masse
              ! x_g**b_geo durch EXP(b_geo*LOG(x_g)) ersetzen, da x_g == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_g = graupel%a_geo * EXP(graupel%b_geo*LOG(x_g))          !..mittlerer Durchmesser
              v_g = graupel%a_vel * EXP(graupel%b_vel*LOG(x_g)) * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

              N_re = v_g * d_g / nu_l                                    !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))           !..mittlerer Vent.Koeff.
              ! ub>> Schnellere Berechnung wg. Verzicht auf "**":
              !              f_n  = a_n + b_n * f_v_fakt * N_re**m_f                    !..mittlerer Vent.Koeff.
              f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer Vent.Koeff.
              ! ub<<
              f_v  = MAX(f_v,1.d0) !unnoetig??
              f_n  = MAX(f_n,1.d0) !unnoetig??

              dep_graupel(i,j,k) = g_i(i,j,k) * n_g * c_g * d_g * f_v * s_si(i,j,k) * dt_local
              dep_graupel_n(i,j,k) = dep_graupel(i,j,k) * f_n/f_v / x_g
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_graupel

    SUBROUTINE vapor_deposition_hail()
      IMPLICIT NONE

      ! Locale Variablen 
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall
      DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h,f_v,f_n,N_re,f_v_fakt,vent_fakt
      DOUBLE PRECISION, SAVE      :: c_h                 !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f,a_n,b_n     !..Koeff. fuer mittleren Ventilationkoeff.

      IF (firstcall.NE.1) THEN
        c_h = 1.0 / hail%cap
        a_n = vent_coeff_a(hail,0)
        b_n = vent_coeff_b(hail,0)
        a_f = vent_coeff_a(hail,1)
        b_f = vent_coeff_b(hail,1)
        IF (isIO() .AND. isdebug) THEN
          WRITE (6, *) "  CLOUD vapor_deposition_hail: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hail:        a_geo = ",hail%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hail:        b_geo = ",hail%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Hail:        a_vel = ",hail%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Hail:        b_vel = ",hail%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Gaskonstante feuchter Luft:   a_dl  = ",a_dl 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_h   = ",c_h 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f   = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f   = ",b_f 
        END IF
        firstcall = 1
      ELSEIF(isIO() .AND. isdebug)THEN
        WRITE (6, *) "  CLOUD vapor_deposition_hail" 
      ENDIF

      f_v_fakt = N_sc**n_f
      vent_fakt = b_n / b_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_hail=0, dep_hail has to be zero too

            IF (q_hail(i,j,k) == 0.0d0) THEN
              dep_hail(i,j,k)   = 0.0d0
              dep_hail_n(i,j,k) = 0.0d0
            ELSE
              n_h = n_hail(i,j,k)                                     !..Anzahldichte
              q_h = q_hail(i,j,k)                                     !..Massendichte

              x_h = MIN(MAX(q_h/(n_h+eps),MAX(hail%x_min,eps)),hail%x_max)  !..mittlere Masse
              ! x_h**b_geo durch EXP(b_geo*LOG(x_h)) ersetzen, da x_h == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_h = hail%a_geo * EXP(hail%b_geo*LOG(x_h))          !..mittlerer Durchmesser
              v_h = hail%a_vel * EXP(hail%b_vel*LOG(x_h)) * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

              N_re = v_h * d_h / nu_l                                    !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))
              ! ub>> Schnellere Berechnung wg. Verzicht auf "**":
              !              f_n  = a_n + b_n * f_v_fakt * N_re**m_f                   !..mittlerer Vent.Koeff.
              f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer Vent.Koeff.
              ! ub<<
              f_v  = MAX(f_v,1.d0) !unnoetig??
              f_n  = MAX(f_n,1.d0) !unnoetig??

              dep_hail(i,j,k) = g_i(i,j,k) * n_h * c_h * d_h * f_v * s_si(i,j,k) * dt_local
              dep_hail_n(i,j,k) = dep_hail(i,j,k) * f_n/f_v / x_h
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_hail

    SUBROUTINE vapor_deposition_snow()
      IMPLICIT NONE

      ! Locale Variablen 
      DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s,f_v,N_re,f_v_fakt
      DOUBLE PRECISION, SAVE      :: c_s             !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall

      IF (firstcall.NE.1) THEN
        c_s = 1.0 / snow%cap
        a_f = vent_coeff_a(snow,1)
        b_f = vent_coeff_b(snow,1)
        IF(isIO() .AND. isdebug)THEN
          WRITE (6, *) "  CLOUD vapor_deposition_snow: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:        a_geo = ",snow%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:        b_geo = ",snow%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Schnee:        a_vel = ",snow%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Schnee:        b_vel = ",snow%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Gaskonstante feuchter Luft:   a_dl = ",a_dl 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittlerer Kapazitaet:         c_s  = ",c_s 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f  = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f  = ",b_f 
        END IF
        firstcall = 1
      ELSEIF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "  CLOUD vapor_deposition_snow: " 
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (q_snow(i,j,k) == 0.0d0) THEN
              dep_snow(i,j,k) = 0.d0
            ELSE
              n_s = n_snow(i,j,k)                        !..Anzahldichte in SI
              q_s = q_snow(i,j,k)                        !..Fluessigwassergehalt in SI

              x_s = MIN(MAX(q_s/(n_s+eps),MAX(snow%x_min,eps)),snow%x_max)   !..mittlere Masse in SI     
              ! x_s**b_geo durch EXP(b_geo*LOG(x_s)) ersetzen, da x_s == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_s = snow%a_geo * EXP(snow%b_geo*LOG(x_s))           !..mittlerer Durchmesser
              v_s = snow%a_vel * EXP(snow%b_vel*LOG(x_s)) * rrho_04(i,j,k)   !..mittlere Sedimentationsgeschw.

              N_re = v_s * d_s / nu_l                         !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))        !..mittlerer Vent.Koeff.
              f_v  = MAX(f_v,1.d0) !unnoetig??

              dep_snow(i,j,k) = g_i(i,j,k) * n_s * c_s * d_s * f_v * s_si(i,j,k) * dt_local
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_snow

    SUBROUTINE saturation_adjust_limiter ()

      ! .. Local Scalars ..
      DOUBLE PRECISION :: e_sat, rdrl, rdrlm1, rlrd, rlrdm1, t_
      DOUBLE PRECISION :: a_, a1, b_, b1, dq_d, lrcv, q_sat, rcv, lrcp, rho_a, q_c

      rcv  = 1.0d0 /  cv
      lrcv = L_wd  * rcv
      lrcp = L_wd  / cp

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            q_d   = q(i,j,k)
            q_c   = q_cloud(i,j,k)
            T_a   = T_0(i,j,k) ! WRF + T(i,j,k) + T_g(i,j,k)
            p_a   = p_0(i,j,k) ! WRF + p(i,j,k) + p_g(i,j,k)
            rho_a = rho_0(i,j,k) ! WRF + rho_g(i,j,k)

            a_    = A_w
            b_    = B_w
            e_sat = e_ws (T_a)

            !...Berechnung der spezifischen Saettigungsfeuchte, Dotzek (4.33)
            q_sat = rlrd / (p_a / e_sat + rlrdm1)

            !...Berechnung der Adjustierungs-Terme, Dotzek, S. 35
            a1 = a_ * (T_3 - b_) / (T_a - b_)**2
            b1 = a1 * lrcp * q_sat

            !...Berechnung des Phasenuebergangs, Dotzek, S. 36
            dq_d = (q_d - q_sat*rho_a) / (1.0d0 + b1)

            !...Limitieren von dep_cloud
            dep_cloud(i,j,k) = MIN (dq_d, dep_cloud(i,j,k))
            dep_cloud(i,j,k) = MAX (q_c,  dep_cloud(i,j,k))

          END DO
        END DO

      END DO

    END SUBROUTINE saturation_adjust_limiter

  END SUBROUTINE vapor_deposition_growth

  SUBROUTINE vapor_dep_simple(dt_local,i_local)
    !*******************************************************************************
    !                                                                              *
    !       Berechnung des Wachstums der Eispartikel durch Wasserdampfdiffusion    *
    !                                                                              *
    !*******************************************************************************


    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, q,        &
         &                        q_cloud, n_cloud, q_ice, n_ice, q_graupel, n_graupel, &
         &                        q_hail, n_hail,                                       & 
         &                        q_snow, n_snow, q_rain, &
         &                        dqdt, speichere_dqdt, cond_neu_sb, &
         &                        evap_neu_sb, speichere_precipstat, &
         &                        deprate_ice, deprate_snow
    USE konstanten,         ONLY: pi,cv,cp
    USE parallele_umgebung, ONLY: isIO
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    INTEGER                     :: i_local
    DOUBLE PRECISION            :: dt_local
    DOUBLE PRECISION            :: a_dl,a_ld,D_vtp
    DOUBLE PRECISION            :: q_g,n_g,x_g,D_g,q_s,n_s,x_s,conv_q,conv_n,x_conv
    ! ub>> 
    DOUBLE PRECISION            :: necessary_d,available_d,weight,depitmp,depstmp,depgtmp,dephtmp
    ! ub<<
    DOUBLE PRECISION, PARAMETER :: eps  = 1.d-20

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: s_si,g_i
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: s_sw,g_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_ice
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_snow
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_cloud
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_graupel, dep_graupel_n
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dep_hail, dep_hail_n

    ! Locale Variablen 
    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: q_si            !..Wasserdampfdichte bei Eissaettigung
    DOUBLE PRECISION            :: e_si            !..Wasserpartialdruck bei Eissaettigung
    DOUBLE PRECISION            :: e_sw            !..Wasserpartialdruck bei saettigung
    DOUBLE PRECISION            :: q_d,x_d,e_d,R_f,p_a,rho_a,dep_sum
    INTEGER                     :: i,j,k
    INTEGER, SAVE               :: firstcall
    DOUBLE PRECISION            :: qimax, qsmax, qgmax, qhmax        ! <hn

    ALLOCATE(s_si(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(s_sw(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(g_i(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(g_w(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_ice(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_snow(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_cloud(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_graupel(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_graupel_n(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_hail(0:loc_ix,1:loc_iy,1:loc_iz))
    ALLOCATE(dep_hail_n(0:loc_ix,1:loc_iy,1:loc_iz))

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF(isIO() .AND. isdebug)THEN
      WRITE (6, *) "CLOUD vapor_deposition_growth: " 
    END IF

    IF (firstcall.NE.1) THEN
      IF(isIO() .AND. isdebug)THEN
        WRITE (6,'(A,D10.3)') "  Diffusivitaet von Wasserdampf in Luft:    D_v   = ",D_v
        WRITE (6,'(A,D10.3)') "  Waermeleitfaehigkeit von Luft:            K_T   = ",K_T
        WRITE (6,'(A,D10.3)') "  Sublimationswaerme:                       L_ed  = ",L_ed
        WRITE (6,'(A,D10.3)') "  Kinematische Viskositaet von Luft:        nu_l  = ",nu_l 
      END IF
      firstcall = 1
    ENDIF

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix
          p_a  = p_0(i,j,k) !WRF!+ p(i,j,k) + p_g(i,j,k)
          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)
          x_d  = q(i,j,k) / (rho_0(i,j,k))!WRF!+rho_g(i,j,k))
          e_d  = q(i,j,k) * R_d * T_a
          e_si = e_es(T_a)
          e_sw = e_ws(T_a)
          s_si(i,j,k) = e_d / e_si - 1.0                    !..Uebersaettigung bzgl. Eis
          !s_sw(i,j,k) = e_d / e_sw - 1.0                    !..Uebersaettigung bzgl. Fluessigwasser
! UB_20090316:          D_vtp = D_v
          D_vtp = 8.7602e-5 * T_a**(1.81) / p_a
          IF (T_a < T_3) THEN
            g_i(i,j,k) = 4.0*pi / ( L_ed**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_si) )
            ! g_w(i,j,k) = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )
          ELSE
            ! g_w(i,j,k)  = 4.0*pi / ( L_wd**2 / (K_T * R_d * T_a**2) + R_d * T_a / (D_vtp * e_sw) )
            g_i(i,j,k)  = 0.0
            s_si(i,j,k) = 0.0
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !CALL vapor_deposition_cloud()
    !CALL saturation_adjust_limiter() ! Limitieren der Kondensation durch Saettigungsadjustierung

    dep_ice     = 0.0
    dep_snow    = 0.0
    dep_cloud   = 0.0
    dep_graupel = 0.0
    dep_hail    = 0.0

    ! ub>>
    ! Die Routinen vapor_deposition_ice(), ..., wurden zeitsparender programmiert, indem
    ! Ausdruecke wie "a ** b" durch exp(b*log(a)) ersetzt wurden und einige Berechnungen
    ! mit Konstanten vor die Schleifen gezogen wurde. Insbesondere durch die andere
    ! Potenzierungstechnik ergeben sich aber numerisch leicht unterschiedliche Ergebnisse.
    ! Im Rahmen der hier verwendeten Depositions-Approximation spielt das aber keine
    ! grosse Rolle, weil die Approximation an sich schon recht ungenau ist.
    ! ub<<

    CALL vapor_deposition_ice()
    CALL vapor_deposition_snow()
    CALL vapor_deposition_graupel()
    IF (ice_typ > 1) CALL vapor_deposition_hail()

    !q_cloud = q_cloud + dep_cloud
    !q       = q       - dep_cloud
    !T       = T       + dep_cloud * L_ed / cp / rho_0(i,j,k) 

    x_conv = (250e-6/snow%a_geo)**(1./snow%b_geo)

    ! ub>>
#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt .AND. i_local == 1) THEN
      ! Da die Deposition u.U. mit einem kleineren Zeitschritt (dt_local) gerechnet wird,
      ! muessen weiter unten die Umwandlungsraten waehrend eines grossen Modellzeitschrittes
      ! (dt) aufsummiert werden. Wenn i_local == 1, dann steht die Zeit auf dem Beginn eines
      ! grossen Zeitschrittes, also dem Anfangszeitpunkt der Summierung. Deshalb:
      ! Zur Sicherheit die Abschnitte des Umwandlungsratenspeicherfeldes dqdt nullen,
      ! die fuer die Depositionsraten gebraucht werden.
      dqdt(:,:,:,4:7) = 0.0
    END IF
#endif
    ! ub<<

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k) 

          ! Nur bei Temperaturen unterhalb des Tripelpunktes wird Deposition/Sublimation parametrisiert.
          ! Bei T > T_3 kann Eis nur dann verdunsten, wenn es vorher anschmilzt.

          IF (T_a < T_3) THEN

            IF (s_si(i,j,k) >= 0.0) THEN
              ! uebersaettigt: dep_ice, dep_snow, dep_graupel sind >= 0. Beschraenken auf 
              ! maximal moeglichen uebersaettigten Wasserdampfanteil:
              q_d  = q(i,j,k)

              ! Saettigungsueberschuss (zu kondensierendes Wasser, kg/m^3): (Naeherung fuer kleine Dampfdruecke)
              available_d = q_d * s_si(i,j,k) / (s_si(i,j,k) + 1.0)
              ! von der Parametrisierung bestimmte Kondensationsmenge waehrend des Zeitschritts:
              necessary_d = dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k)

              IF (ABS(necessary_d) < eps) THEN
                dep_ice(i,j,k)     = 0.0
                dep_snow(i,j,k)    = 0.0
                dep_graupel(i,j,k) = 0.0               
                dep_hail(i,j,k)    = 0.0
              ELSEIF (necessary_d > available_d) THEN
                weight = available_d / necessary_d
                dep_ice(i,j,k)     = weight * dep_ice(i,j,k)
                dep_snow(i,j,k)    = weight * dep_snow(i,j,k)
                dep_graupel(i,j,k) = weight * dep_graupel(i,j,k)
                dep_hail(i,j,k)    = weight * dep_hail(i,j,k)
              ENDIF

            ELSE
              ! untersaettigt: dep_ice, dep_snow, dep_graupel sind < 0. Beschraenken auf 
              ! maximal moeglichen untersaettigten Wasserdampfanteil sowie auf maximal vorhandenes Eis.
              ! Die Notwendigkeit einer Begrenzung entsteht durch einen i.a. viel zu grossen
              ! Depositionszeitschritt, so dass das hier verwendete 
              ! Euler-Vorwaerts-Zeitintegrationsschema hoffnungslos instabil ist.

              ! Hier wird im Gegensatz zum uebersaettigten Fall zuerst jeder einzelne Depositionsfluss
              ! auf das maximal vorhandene Reservoir begrenzt, bevor durch weight auf den
              ! maximal moeglichen Dampfinput fuer die untersaettigte Luft begrenzt wird.
              depitmp    = MAX(dep_ice(i,j,k),    -q_ice(i,j,k))
              depstmp    = MAX(dep_snow(i,j,k),   -q_snow(i,j,k))
              dephtmp    = MAX(dep_hail(i,j,k),   -q_hail(i,j,k))  
              depgtmp    = MAX(dep_graupel(i,j,k),-q_graupel(i,j,k)) 

              ! Saettigungsdefizit (Wasserdampf-Aufnahmereservoir der Luft, kg/m^3):
              available_d = s_si(i,j,k) * e_es(T_0(i,j,k)) / (R_d * T_0(i,j,k))
              ! von der Parametrisierung bestimmte Verdunstungsmenge waehrend des Zeitschritts:
              necessary_d = depitmp + depstmp + depgtmp + dephtmp

              ! available_d und necessary_d sind hier < 0!
              IF (necessary_d < -eps) THEN
                IF (necessary_d < available_d) THEN
                  ! In diesem Falle kann waehrend des Zeitschritts nicht soviel verdunstet werden 
                  ! wie von der Parametrisierung berechnet, weil das Saettigungsdefizit nicht gross genug ist. 
                  ! Die Verdunstung wird auf das Saettigungsdefizit begrenzt

                  weight = available_d / necessary_d
                  dep_ice(i,j,k)     = weight * depitmp
                  dep_snow(i,j,k)    = weight * depstmp
                  dep_graupel(i,j,k) = weight * depgtmp
                  dep_hail(i,j,k)    = weight * dephtmp


                ELSE

                  ! In diesem Falle kann die vollstaendige von der Parametrisierung berechnete
                  ! Verdunstungsmenge aller Hydrometeorarten verdunstet werden. Keine Begrenzung ist noetig 
                  ! (bis auf die Begrenzung auf die vorhandene Hydrometeormasse).
                  dep_ice(i,j,k)     = depitmp
                  dep_snow(i,j,k)    = depstmp
                  dep_graupel(i,j,k) = depgtmp
                  dep_hail(i,j,k)    = dephtmp

                END IF

              ELSE

                ! Deposition so gering, dass man sie getrost zu 0 setzen kann:
                dep_ice(i,j,k)     = 0.0
                dep_snow(i,j,k)    = 0.0
                dep_graupel(i,j,k) = 0.0    
                dep_hail(i,j,k)    = 0.0                    

              END IF

            END IF

            ! Zur Sicherheit: Pruefen auf Qx < 0.0:
            ! An dieser Stelle kann es vorkommen, dass ein |dep_xxx| > q_xxx ist, und zwar im Falle sehr kleiner absoluter Werte.
            IF (q_ice(i,j,k)+dep_ice(i,j,k) < 0.0 .OR. &
                 q_snow(i,j,k)+dep_snow(i,j,k) < 0.0 .OR. &
                 q_graupel(i,j,k)+dep_graupel(i,j,k) < 0.0 .OR. &
                 q_hail(i,j,k)+dep_hail(i,j,k) < 0.0) THEN

!!$              IF (q_ice(i,j,k)+dep_ice(i,j,k) < -1d-20 .OR. &
!!$                   q_snow(i,j,k)+dep_snow(i,j,k) < -1d-20 .OR. &
!!$                   q_graupel(i,j,k)+dep_graupel(i,j,k) < -1d-20 .OR. &
!!$                   q_hail(i,j,k)+dep_hail(i,j,k) < -1d-20) THEN
!!$
!!$                WRITE (*,*) 'VAPOR_DEPOSITION: Warnung --- Qx < 0.0 aufgetreten (wird spaeter begrenzt):'
!!$                WRITE (*,*) '     ----- Kann vorkommen bei sehr kleinen Werten von Qx und dep_x -----'
!!$                WRITE (*,*) 'Punkt ', i, ',', j, ',', k
!!$                WRITE (*,*) 'q_ice     = ', q_ice(i,j,k)+dep_ice(i,j,k), '  dep_ice = ', dep_ice(i,j,k)
!!$                WRITE (*,*) 'q_snow    = ', q_snow(i,j,k)+dep_snow(i,j,k), '  dep_snow = ', dep_snow(i,j,k)
!!$                WRITE (*,*) 'q_graupel = ', q_graupel(i,j,k)+dep_graupel(i,j,k), '  dep_graupel = ', dep_graupel(i,j,k)
!!$                WRITE (*,*) 'q_hail    = ', q_hail(i,j,k)+dep_hail(i,j,k), '  dep_hail = ', dep_hail(i,j,k)
!!$                WRITE (*,*) 'tv_x_max = dt_local / eps = ', dt_local / eps
!!$              END IF
              !AS 20060220>
              ! Deshalb hier nochmals eine letzte Begrenzung:
              dep_hail(i,j,k)    = MAX(dep_hail(i,j,k),   -q_hail(i,j,k)   )
              dep_graupel(i,j,k) = MAX(dep_graupel(i,j,k),-q_graupel(i,j,k))
              dep_snow(i,j,k)    = MAX(dep_snow(i,j,k),   -q_snow(i,j,k)   )
              dep_ice(i,j,k)     = MAX(dep_ice(i,j,k),    -q_ice(i,j,k)    )
              !<end AS
            END IF

            dep_sum = dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k)
            IF (dep_sum > q(i,j,k) .AND. dep_sum > 0.0) THEN
!!$              IF (dep_sum > 1d-20) THEN
!!$                WRITE (*,*) 'VAPOR_DEPOSITION: Warnung --- dep_sum > q aufgetreten (wird spaeter begrenzt):'
!!$                WRITE (*,*) '     ----- Sollte eigentlich nicht vorkommen -----'
!!$                WRITE (*,*) 'Punkt ', i, ',', j, ',', k
!!$                WRITE (*,*) 'q_ice     = ', q_ice(i,j,k)+dep_ice(i,j,k), '  dep_ice = ', dep_ice(i,j,k)
!!$                WRITE (*,*) 'q_snow    = ', q_snow(i,j,k)+dep_snow(i,j,k), '  dep_snow = ', dep_snow(i,j,k)
!!$                WRITE (*,*) 'q_graupel = ', q_graupel(i,j,k)+dep_graupel(i,j,k), '  dep_graupel = ', dep_graupel(i,j,k)
!!$                WRITE (*,*) 'q_hail    = ', q_hail(i,j,k)+dep_hail(i,j,k), '  dep_hail = ', dep_hail(i,j,k)
!!$                WRITE (*,*) 'q = ', q(i,j,k)+dep_sum, '   dep_sum = ', dep_sum
!!$              END IF
              ! Kondensation begrenzen auf max. vorhandene Feuchtigkeit 
              ! (sollte aber nach dem vorhergegangenen Begrenzen wirklich nicht noetig sein ...):
              weight = q(i,j,k) / dep_sum
              dep_ice(i,j,k)     = weight * dep_ice(i,j,k)
              dep_snow(i,j,k)    = weight * dep_snow(i,j,k)
              dep_graupel(i,j,k) = weight * dep_graupel(i,j,k)
              dep_hail(i,j,k)    = weight * dep_hail(i,j,k)
              dep_sum = q(i,j,k)
            END IF

            ! Moeglichkeit der hail-to-snow conversion momentan nicht vorgesehen
            n_g = n_graupel(i,j,k)
            q_g = q_graupel(i,j,k)                  
            x_g = MIN(MAX(q_g/(n_g+eps),graupel%x_min),graupel%x_max)
            D_g = graupel%a_geo * x_g**graupel%b_geo   
            IF (.FALSE. .AND. dep_graupel(i,j,k) > 0 .AND. D_g < 800d-6) THEN
              !if (dep_graupel(i,j,k) > 0 .AND. D_g < 800d-6 .AND. q_cloud(i,j,k) < 1e-4) then
              !..Graupel to snow conversion
              conv_q = MIN(q_g, 3.0 * dep_graupel(i,j,k))
              !conv_n = conv_q/x_g
              !conv_n = MIN(n_g, 0.25 * dep_graupel_n(i,j,k))
              q_snow(i,j,k)    = q_snow(i,j,k)    + conv_q
              q_graupel(i,j,k) = q_graupel(i,j,k) - conv_q + dep_graupel(i,j,k)
              n_snow(i,j,k)    = n_snow(i,j,k)    + conv_q / MIN(x_g,x_conv)
              n_graupel(i,j,k) = n_graupel(i,j,k) - conv_q / x_g      
            ELSE
              q_graupel(i,j,k) = q_graupel(i,j,k) + dep_graupel(i,j,k)
              conv_q = 0.0
            ENDIF

            q_ice(i,j,k)     = q_ice(i,j,k)     + dep_ice(i,j,k)
            q_snow(i,j,k)    = q_snow(i,j,k)    + dep_snow(i,j,k)
            IF (ice_typ > 1) q_hail(i,j,k)    = q_hail(i,j,k)    + dep_hail(i,j,k)    ! <hn

            q(i,j,k) = q(i,j,k) - dep_sum

            ! ub>>
#ifdef SAVE_CONVERSIONRATES
            IF (speichere_dqdt) THEN
              ! Wegen kleinerem Zeitschritt dt_local: Aufsummieren, nicht einfach den Wert setzen.
              IF (dep_sum >= 0.0) THEN
                dqdt(i,j,k,4) = dqdt(i,j,k,4) + dep_ice(i,j,k)
                dqdt(i,j,k,5) = dqdt(i,j,k,5) + dep_snow(i,j,k)
                dqdt(i,j,k,6) = dqdt(i,j,k,6) + dep_graupel(i,j,k)
                dqdt(i,j,k,55) = dqdt(i,j,k,55) + dep_hail(i,j,k)
                dqdt(i,j,k,7) = dqdt(i,j,k,7) + conv_q
              ELSE
                dqdt(i,j,k,56) = dqdt(i,j,k,56) - dep_ice(i,j,k)
                dqdt(i,j,k,32) = dqdt(i,j,k,32) - dep_snow(i,j,k)
                dqdt(i,j,k,33) = dqdt(i,j,k,33) - dep_graupel(i,j,k)
                dqdt(i,j,k,53) = dqdt(i,j,k,53) - dep_hail(i,j,k)
              END IF
            END IF
#endif

            IF (speichere_precipstat) THEN
              IF (dep_sum >= 0.0) THEN
                cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + dep_sum
              ELSE
                evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) - dep_sum
              END IF
            END IF
            ! ub<<

 

            IF (use_ice_graupel_conv_uli) THEN
              deprate_ice(i,j,k) = deprate_ice(i,j,k) + dep_ice(i,j,k)
              deprate_snow(i,j,k) = deprate_snow(i,j,k) + dep_snow(i,j,k)
            END IF

            !WRF!T(i,j,k) = T(i,j,k) + dep_sum * L_ed / cp / rho_0(i,j,k) 
            !WRF!p(i,j,k) = p(i,j,k) + dep_sum * L_ed / cv * R_l


          ENDIF
        ENDDO
      ENDDO
    ENDDO

!!$#ifdef SAVE_CONVERSIONRATES
    IF (speichere_dqdt) THEN
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            IF (dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k) >= 0.0) THEN
              dqdt(i,j,k,4) = dqdt(i,j,k,4) + dep_ice(i,j,k)
              dqdt(i,j,k,5) = dqdt(i,j,k,5) + dep_snow(i,j,k)
              dqdt(i,j,k,6) = dqdt(i,j,k,6) + dep_graupel(i,j,k)
              dqdt(i,j,k,55) = dqdt(i,j,k,55) + dep_hail(i,j,k)
              dqdt(i,j,k,7) = dqdt(i,j,k,7) + conv_q
            ELSE
              dqdt(i,j,k,56) = dqdt(i,j,k,56) - dep_ice(i,j,k)
              dqdt(i,j,k,32) = dqdt(i,j,k,32) - dep_snow(i,j,k)
              dqdt(i,j,k,33) = dqdt(i,j,k,33) - dep_graupel(i,j,k)
              dqdt(i,j,k,53) = dqdt(i,j,k,53) - dep_hail(i,j,k)
            END IF

          END DO
        END DO
      END DO
    END IF
!!$#endif

    IF (speichere_precipstat) THEN
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            dep_sum = dep_ice(i,j,k) + dep_graupel(i,j,k) + dep_snow(i,j,k) + dep_hail(i,j,k)
            IF (dep_sum >= 0.0) THEN
              cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + dep_sum
            ELSE
              evap_neu_sb(i,j,k) = evap_neu_sb(i,j,k) - dep_sum
            END IF

          END DO
        END DO
      END DO
    END IF


    DEALLOCATE(s_si,s_sw,g_i,g_w,dep_ice,dep_cloud,dep_snow,dep_graupel,dep_graupel_n, &
         &  dep_hail,dep_hail_n)

  CONTAINS


    SUBROUTINE vapor_deposition_cloud()
      IMPLICIT NONE

      DOUBLE PRECISION            :: q_c,n_c,x_c,d_c,v_c,f_v,N_re,f_v_fakt
      DOUBLE PRECISION, SAVE      :: c_c             !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall

      IF (firstcall.NE.1) THEN
        c_c = 1.0 / cloud%cap
        a_f = vent_coeff_a(cloud,1)
        b_f = vent_coeff_b(cloud,1)
        IF(isIO())THEN
          WRITE (6, *) "  CLOUD vapor_deposition_cloud: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentropfen:  a_cloud   = ",cloud%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Wolkentropfen:  b_cloud   = ",cloud%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation:                alf_cloud = ",cloud%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation:                bet_cloud = ",cloud%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_c       = ",c_c 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f       = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f       = ",b_f 
        END IF
        firstcall = 1
      ELSEIF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "  CLOUD vapor_deposition_cloud " 
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            n_c = n_cloud(i,j,k)                                   !..Anzahldichte
            q_c = q_cloud(i,j,k)                                   !..Massendichte

            x_c = MIN(MAX(q_c/(n_c+eps),MAX(cloud%x_min,eps)),cloud%x_max)  !..mittlere Masse
            ! x_c**b_geo durch EXP(b_geo*LOG(x_c)) ersetzen, da x_c == 0 nicht mehr vorkommen kann; 20 % schneller:
            d_c = cloud%a_geo * EXP(cloud%b_geo*LOG(x_c))          !..mittlerer Durchmesser
            v_c = cloud%a_vel * EXP(cloud%b_vel*LOG(x_c)) * rrho_c(i,j,k)   !..mittlere Sedimentationsgeschw.

            N_re = v_c * d_c / nu_l                                !..mittlere Reynoldszahl
            ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
            f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))             !..mittlerer Vent.Koeff.
            f_v  = MAX(f_v,1.d0) !unnoetig??

            dep_cloud(i,j,k) = g_w(i,j,k) * n_c * c_c * d_c * f_v * s_sw(i,j,k) * dt_local
          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_cloud

    SUBROUTINE vapor_deposition_ice()
      IMPLICIT NONE

      DOUBLE PRECISION            :: q_i,n_i,x_i,d_i,v_i,f_v,N_re,f_v_fakt
      DOUBLE PRECISION, SAVE      :: c_i             !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall

      IF (firstcall.NE.1) THEN
        c_i = 1.0 / ice%cap
        a_f = vent_coeff_a(ice,1)
        b_f = vent_coeff_b(ice,1)
        IF(isIO() .AND. isdebug)THEN
          WRITE (6, *) "  CLOUD vapor_deposition_ice: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:            a_ice   = ",ice%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Eis:            b_ice   = ",ice%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation von Eis:        alf_ice = ",ice%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sedimentation von Eis:        bet_ice = ",ice%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_i     = ",c_i 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f     = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f     = ",b_f 
        END IF
        firstcall = 1
      ELSEIF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "  CLOUD vapor_deposition_ice " 
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_ice=0, dep_ice has to be zero too

            IF (q_ice(i,j,k) == 0.0d0) THEN
              dep_ice(i,j,k)   = 0.0d0
            ELSE  
              n_i = n_ice(i,j,k)                                   !..Anzahldichte
              q_i = q_ice(i,j,k)                                   !..Massendichte

              x_i = MIN(MAX(q_i/(n_i+eps),MAX(ice%x_min,eps)),ice%x_max)    !..mittlere Masse
              ! x_i**b_geo durch EXP(b_geo*LOG(x_i)) ersetzen, da x_i == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_i = ice%a_geo * EXP(ice%b_geo*LOG(x_i))            !..mittlerer Durchmesser
              v_i = ice%a_vel * EXP(ice%b_vel*LOG(x_i)) * rrho_04(i,j,k)    !..mittlere Sedimentationsgeschw.

              N_re = v_i * d_i / nu_l                              !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))     !..mittlerer Vent.Koeff.
              !f_v  = MAX(f_v,1.d0) !unnoetig??

              dep_ice(i,j,k) = g_i(i,j,k) * n_i * c_i * d_i * f_v * s_si(i,j,k) * dt_local
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_ice

    SUBROUTINE vapor_deposition_graupel()
      IMPLICIT NONE

      ! Locale Variablen 
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall
      DOUBLE PRECISION            :: q_g,n_g,x_g,d_g,v_g,f_v,f_n,N_re,f_v_fakt,vent_fakt
      DOUBLE PRECISION, SAVE      :: c_g                 !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f,a_n,b_n     !..Koeff. fuer mittleren Ventilationkoeff.

      IF (firstcall.NE.1) THEN
        c_g = 1.0 / graupel%cap
        a_n = vent_coeff_a(graupel,0)
        b_n = vent_coeff_b(graupel,0)
        a_f = vent_coeff_a(graupel,1)
        b_f = vent_coeff_b(graupel,1)
        IF (isIO() .AND. isdebug) THEN
          WRITE (6, *) "  CLOUD vapor_deposition_graupel: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:        a_geo = ",graupel%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Graupel:        b_geo = ",graupel%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Graupel:        a_vel = ",graupel%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Graupel:        b_vel = ",graupel%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Gaskonstante feuchter Luft:   a_dl  = ",a_dl 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_g   = ",c_g 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f   = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f   = ",b_f 
        END IF
        firstcall = 1
      ELSEIF(isIO() .AND. isdebug)THEN
        WRITE (6, *) "  CLOUD vapor_deposition_graupel" 
      ENDIF

      f_v_fakt = N_sc**n_f
      vent_fakt = b_n / b_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_garupel=0, dep_graupel has to be zero too

            IF (q_graupel(i,j,k) == 0.0d0) THEN
              dep_graupel(i,j,k)   = 0.0d0
              dep_graupel_n(i,j,k) = 0.0d0
            ELSE

              n_g = n_graupel(i,j,k)                                     !..Anzahldichte
              q_g = q_graupel(i,j,k)                                     !..Massendichte

              x_g = MIN(MAX(q_g/(n_g+eps),MAX(graupel%x_min,eps)),graupel%x_max)  !..mittlere Masse
              ! x_g**b_geo durch EXP(b_geo*LOG(x_g)) ersetzen, da x_g == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_g = graupel%a_geo * EXP(graupel%b_geo*LOG(x_g))          !..mittlerer Durchmesser
              v_g = graupel%a_vel * EXP(graupel%b_vel*LOG(x_g)) * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

              N_re = v_g * d_g / nu_l                                    !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))           !..mittlerer Vent.Koeff.
              ! ub>> Schnellere Berechnung wg. Verzicht auf "**":
              !              f_n  = a_n + b_n * f_v_fakt * N_re**m_f                    !..mittlerer Vent.Koeff.
              f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer Vent.Koeff.
              ! ub<<
              f_v  = MAX(f_v,1.d0) !unnoetig??
              f_n  = MAX(f_n,1.d0) !unnoetig??

              dep_graupel(i,j,k) = g_i(i,j,k) * n_g * c_g * d_g * f_v * s_si(i,j,k) * dt_local
              dep_graupel_n(i,j,k) = dep_graupel(i,j,k) * f_n/f_v / x_g
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_graupel

    SUBROUTINE vapor_deposition_hail()
      IMPLICIT NONE

      ! Locale Variablen 
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall
      DOUBLE PRECISION            :: q_h,n_h,x_h,d_h,v_h,f_v,f_n,N_re,f_v_fakt,vent_fakt
      DOUBLE PRECISION, SAVE      :: c_h                 !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f,a_n,b_n     !..Koeff. fuer mittleren Ventilationkoeff.

      IF (firstcall.NE.1) THEN
        c_h = 1.0 / hail%cap
        a_n = vent_coeff_a(hail,0)
        b_n = vent_coeff_b(hail,0)
        a_f = vent_coeff_a(hail,1)
        b_f = vent_coeff_b(hail,1)
        IF (isIO() .AND. isdebug) THEN
          WRITE (6, *) "  CLOUD vapor_deposition_hail: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hail:        a_geo = ",hail%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Hail:        b_geo = ",hail%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Hail:        a_vel = ",hail%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Hail:        b_vel = ",hail%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Gaskonstante feuchter Luft:   a_dl  = ",a_dl 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Kapazitaet:                   c_h   = ",c_h 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f   = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f   = ",b_f 
        END IF
        firstcall = 1
      ELSEIF(isIO() .AND. isdebug)THEN
        WRITE (6, *) "  CLOUD vapor_deposition_hail" 
      ENDIF

      f_v_fakt = N_sc**n_f
      vent_fakt = b_n / b_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix

            ! hn: in case q_hail=0, dep_hail has to be zero too

            IF (q_hail(i,j,k) == 0.0d0) THEN
              dep_hail(i,j,k)   = 0.0d0
              dep_hail_n(i,j,k) = 0.0d0
            ELSE
              n_h = n_hail(i,j,k)                                     !..Anzahldichte
              q_h = q_hail(i,j,k)                                     !..Massendichte

              x_h = MIN(MAX(q_h/(n_h+eps),MAX(hail%x_min,eps)),hail%x_max)  !..mittlere Masse
              ! x_h**b_geo durch EXP(b_geo*LOG(x_h)) ersetzen, da x_h == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_h = hail%a_geo * EXP(hail%b_geo*LOG(x_h))          !..mittlerer Durchmesser
              v_h = hail%a_vel * EXP(hail%b_vel*LOG(x_h)) * rrho_04(i,j,k)  !..mittlere Sedimentationsgeschw.

              N_re = v_h * d_h / nu_l                                    !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))
              ! ub>> Schnellere Berechnung wg. Verzicht auf "**":
              !              f_n  = a_n + b_n * f_v_fakt * N_re**m_f                   !..mittlerer Vent.Koeff.
              f_n = a_n + vent_fakt * (f_v - a_f)                        !..mittlerer Vent.Koeff.
              ! ub<<
              f_v  = MAX(f_v,1.d0) !unnoetig??
              f_n  = MAX(f_n,1.d0) !unnoetig??

              dep_hail(i,j,k) = g_i(i,j,k) * n_h * c_h * d_h * f_v * s_si(i,j,k) * dt_local
              dep_hail_n(i,j,k) = dep_hail(i,j,k) * f_n/f_v / x_h
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_hail

    SUBROUTINE vapor_deposition_snow()
      IMPLICIT NONE

      ! Locale Variablen 
      DOUBLE PRECISION            :: q_s,n_s,x_s,d_s,v_s,f_v,N_re,f_v_fakt
      DOUBLE PRECISION, SAVE      :: c_s             !..Koeff. fuer mittlere Kapazitaet
      DOUBLE PRECISION, SAVE      :: a_f,b_f         !..Koeff. fuer mittleren Ventilationkoeff.
      INTEGER                     :: i,j,k
      INTEGER, SAVE               :: firstcall

      IF (firstcall.NE.1) THEN
        c_s = 1.0 / snow%cap
        a_f = vent_coeff_a(snow,1)
        b_f = vent_coeff_b(snow,1)
        IF(isIO() .AND. isdebug)THEN
          WRITE (6, *) "  CLOUD vapor_deposition_snow: " 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:        a_geo = ",snow%a_geo  
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Geometrie von Schnee:        b_geo = ",snow%b_geo
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Schnee:        a_vel = ",snow%a_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Sediment. von Schnee:        b_vel = ",snow%b_vel 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer Gaskonstante feuchter Luft:   a_dl = ",a_dl 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittlerer Kapazitaet:         c_s  = ",c_s 
          WRITE (6,'(A,D10.3)') "  Koeff. fuer mittleren Ventilationskoeff.: a_f  = ",a_f 
          WRITE (6,'(A,D10.3)') "                                            b_f  = ",b_f 
        END IF
        firstcall = 1
      ELSEIF (isIO() .AND. isdebug) THEN
        WRITE (6, *) "  CLOUD vapor_deposition_snow: " 
      ENDIF

      f_v_fakt = N_sc**n_f
      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            IF (q_snow(i,j,k) == 0.0d0) THEN
              dep_snow(i,j,k) = 0.d0
            ELSE
              n_s = n_snow(i,j,k)                        !..Anzahldichte in SI
              q_s = q_snow(i,j,k)                        !..Fluessigwassergehalt in SI

              x_s = MIN(MAX(q_s/(n_s+eps),MAX(snow%x_min,eps)),snow%x_max)   !..mittlere Masse in SI     
              ! x_s**b_geo durch EXP(b_geo*LOG(x_s)) ersetzen, da x_s == 0 nicht mehr vorkommen kann; 20 % schneller:
              d_s = snow%a_geo * EXP(snow%b_geo*LOG(x_s))           !..mittlerer Durchmesser
              v_s = snow%a_vel * EXP(snow%b_vel*LOG(x_s)) * rrho_04(i,j,k)   !..mittlere Sedimentationsgeschw.

              N_re = v_s * d_s / nu_l                         !..mittlere Reynoldszahl
              ! N_re**m_f durch EXP(m_f*LOG(N_re)) ersetzen, da N_re == 0 nicht mehr vorkommen kann:
              f_v  = a_f + b_f * f_v_fakt * EXP(m_f*LOG(N_re))        !..mittlerer Vent.Koeff.
              f_v  = MAX(f_v,1.d0) !unnoetig??

              dep_snow(i,j,k) = g_i(i,j,k) * n_s * c_s * d_s * f_v * s_si(i,j,k) * dt_local
            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE vapor_deposition_snow

    SUBROUTINE saturation_adjust_limiter ()

      ! .. Local Scalars ..
      DOUBLE PRECISION :: e_sat, rdrl, rdrlm1, rlrd, rlrdm1, t_
      DOUBLE PRECISION :: a_, a1, b_, b1, dq_d, lrcv, q_sat, rcv, lrcp, rho_a, q_c

      rcv  = 1.0d0 /  cv
      lrcv = L_wd  * rcv
      lrcp = L_wd  / cp

      DO k = 1, loc_iz
        DO j = 1, loc_iy
          DO i = 0, loc_ix
            q_d   = q(i,j,k)
            q_c   = q_cloud(i,j,k)
            T_a   = T_0(i,j,k) ! WRF + T(i,j,k) + T_g(i,j,k)
            p_a   = p_0(i,j,k) ! WRF + p(i,j,k) + p_g(i,j,k)
            rho_a = rho_0(i,j,k) ! WRF + rho_g(i,j,k)

            a_    = A_w
            b_    = B_w
            e_sat = e_ws (T_a)

            !...Berechnung der spezifischen Saettigungsfeuchte, Dotzek (4.33)
            q_sat = rlrd / (p_a / e_sat + rlrdm1)

            !...Berechnung der Adjustierungs-Terme, Dotzek, S. 35
            a1 = a_ * (T_3 - b_) / (T_a - b_)**2
            b1 = a1 * lrcp * q_sat

            !...Berechnung des Phasenuebergangs, Dotzek, S. 36
            dq_d = (q_d - q_sat*rho_a) / (1.0d0 + b1)

            !...Limitieren von dep_cloud
            dep_cloud(i,j,k) = MIN (dq_d, dep_cloud(i,j,k))
            dep_cloud(i,j,k) = MAX (q_c,  dep_cloud(i,j,k))

          END DO
        END DO

      END DO

    END SUBROUTINE saturation_adjust_limiter

  END SUBROUTINE vapor_dep_simple

  SUBROUTINE cloud_nucleation()
    !*******************************************************************************
    !                                                                              *
    !       Berechnung der Nukleation der Wolkentropfen                            *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, q_cloud, n_cloud, &
         &                        T_g, p_g, rho_g, w, dt,dz3d, S_w, dSwdz, zml_k,    &
         &                        dqdt, speichere_dqdt, cond_neu_sb,     &
         &                        evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: r,cv,cp,dz 
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    DOUBLE PRECISION            :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION            :: p_a,q_d,x_d,e_d,n_c,rho_a
    DOUBLE PRECISION            :: nuc_n, nuc_q, S, dS, dSdz_w
    DOUBLE PRECISION            :: a_ld,a_dl,e_sw,q_sw
    DOUBLE PRECISION            :: N_ccn,k_ccn,N_max,N_min,S_max
    DOUBLE PRECISION            :: z0_nccn, z1e_nccn
    INTEGER                     :: i,j,k,nuc_typ,stat_var

    !DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: s_sw,dSdz
    !ALLOCATE(s_sw(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    !ALLOCATE(dSdz(0:loc_ix,1:loc_iy,1:loc_iz), STAT=stat_var)
    !IF (stat_var/=0) CALL abortparallel("Allokierungsfehler cloud_nucleation",4)

    nuc_typ = nuc_c_typ

    IF(isIO() .AND. isdebug) THEN
      WRITE (6, *) "CLOUD cloud_nucleation: nuc_typ = ",nuc_typ
    ENDIF

    a_dl = R_d/R_l
    a_ld = 1.0/a_dl

    IF (nuc_typ == 1) THEN
      !...HUCM maritime case
      N_ccn = 100.0d06
      N_max = 150.0d06
      N_min =  10.0d06
      S_max = 1.100
      k_ccn = 0.462
    ELSEIF (nuc_typ == 3) THEN
      !... maritime case
      N_ccn = 300.0d06
      N_max = 600.0d06
      N_min =  50.0d06
      S_max = 2.000
      k_ccn = 0.462
    ELSEIF (nuc_typ == 2) THEN
      !...HUCM continental case (Texas CCN)
      N_ccn = 1260.0d06
      N_max = 3000.0d06
      N_min =  300.0d06
      S_max = 20.000
      k_ccn = 0.308
    ELSEIF (nuc_typ == 4) THEN
      !...conti / 2
      N_ccn = 1260.0d06 / 2.0
      N_max = 3000.0d06 / 2.0
      N_min =  300.0d06 / 2.0 
      S_max = 20.000
      k_ccn = 0.308
    ELSEIF (nuc_typ == 5) THEN
      !...conti / 2
      N_ccn = 1260.0d06 / 3.0
      N_max = 3000.0d06 / 3.0
      N_min =  300.0d06 / 3.0 
      S_max = 20.000
      k_ccn = 0.308
    ELSEIF (nuc_typ == 6) THEN
      !...conti / 2
      N_ccn = 1260.0d06 / 3.0 * 2.
      N_max = 3000.0d06 / 3.0 * 2.
      N_min =  300.0d06 / 3.0 * 2.
      S_max = 20.000
      k_ccn = 0.308
    ELSEIF (nuc_typ == 7) THEN
      !...conti / 2
      N_ccn = 1260.0d06 / 3.0 * 5.
      N_max = 3000.0d06 / 3.0 * 5.
      N_min =  300.0d06 / 3.0 * 5.
      S_max = 20.000
      k_ccn = 0.308
    ELSE
      S_max = 1d3
    ENDIF

    !..parameter for exponential decrease of N_ccn with height:
    !  1) up to this height (m) constant unchanged value:
    z0_nccn  = 99999.0d0
    !  2) above, height interval at which N_ccn decreses by factor 1/e:
    z1e_nccn = 99999.0d0

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

          S  = 1.0d2 * s_w(i,j,k)

          dSdz_w = 1.0d2 * dSwdz(i,j,k)

          T_a  = T_0(i,j,k) !WRF!+ T(i,j,k) + T_g(i,j,k)

          IF ( S <= 0.0 .OR. dSdz_w < 0.0 .OR. S > S_max .OR. T_a < T_f) CYCLE

          nuc_n = 0.0
          nuc_q = 0.0
          n_c   = n_cloud(i,j,k)
          rho_a = rho_0(i,j,k)!WRF!+rho_g(i,j,k)

          IF (nuc_typ > 9) THEN
!AS20081207> alte Variante auskommentiert wegen Vektorisierung
!            nuc_n = n_ccn_diagnostic()
!            nuc_n = nuc_n * MIN(EXP((z0_nccn-zml_k(i,j,k))/z1e_nccn), 1.0d0)
!            nuc_n = MAX(nuc_n - n_c,0.d0)
!            nuc_q = nuc_n * cloud%x_min
          ELSE
            IF (zml_k(i,j,k) <= z0_nccn) THEN
              nuc_n = N_ccn * k_ccn * S**(k_ccn-1.0) * dSdz_w * dt
            ELSE
              ! exponential decrease of N_ccn with height
              nuc_n = ( N_ccn * k_ccn * S**(k_ccn-1.0) * dSdz_w - &
                   N_ccn / z1e_nccn * S**k_ccn ) * EXP((z0_nccn-zml_k(i,j,k))/z1e_nccn) * dt
            END IF
            nuc_n = MAX(MIN(nuc_n,N_max-n_c),0.d0)
            nuc_q = nuc_n * cloud%x_min
          ENDIF
          IF (nuc_q > q(i,j,k)) THEN
            nuc_q = q(i,j,k)
            nuc_n = nuc_q / cloud%x_min
          ENDIF

          n_cloud(i,j,k) = n_cloud(i,j,k) + nuc_n
          q_cloud(i,j,k) = q_cloud(i,j,k) + nuc_q
          q(i,j,k)       = q(i,j,k)       - nuc_q             

!AS20080426> hard upper and lower limit for number conc that
!            eliminates also unrealistic high value
!            that would come from the dynamical core

          n_cloud(i,j,k) = MIN(n_cloud(i,j,k),N_max)
          n_cloud(i,j,k) = MAX(n_cloud(i,j,k),N_min)

!<AS20080426

          ! ub>>
#ifdef SAVE_CONVERSIONRATES
          IF (speichere_dqdt) THEN
            dqdt(i,j,k,1) = nuc_q
          END IF
#endif
          IF (speichere_precipstat) THEN
            cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + nuc_q
          END IF
          ! ub<<

          !WRF!T(i, j, k) = T(i, j, k) + L_wd / cp * nuc_q / rho_a
          !WRF!p(i, j, k) = p(i, j, k) + L_wd / cv * nuc_q * r

        ENDDO
      ENDDO
    ENDDO

    WHERE (q_cloud > 1e-9) 
      n_cloud = MAX(n_cloud,N_min)
    ENDWHERE

  CONTAINS

    DOUBLE PRECISION FUNCTION n_ccn_diagnostic()
      IMPLICIT NONE

      IF (nuc_typ == 3) THEN
        n_ccn_diagnostic = n_ccn_maritim()
      ELSEIF (nuc_typ == 4) THEN
        n_ccn_diagnostic = n_ccn_texas()
      ELSEIF (nuc_typ == 5) THEN
        n_ccn_diagnostic = n_ccn_synthetic_maritime()
      ELSEIF (nuc_typ == 6) THEN
        n_ccn_diagnostic = n_ccn_synthetic_continental()
      ENDIF

      RETURN
    END FUNCTION n_ccn_diagnostic

    DOUBLE PRECISION FUNCTION n_ccn_maritim()

      IMPLICIT NONE

      DOUBLE PRECISION :: S

      DOUBLE PRECISION, PARAMETER :: N_a   = 100.0d06
      DOUBLE PRECISION, PARAMETER :: N_max = 150.0d06
      DOUBLE PRECISION, PARAMETER :: k_a    = 0.462

      !..CCN-Spektrum nach HUCM
      !  (maritime case)

      S = MIN(1.1d0,0.5*1.0d2*s_w(i,j,k))              !..Uebersaettigung in Prozent

      n_ccn_maritim = N_a * S**k_a                      !..Power-Law

      n_ccn_maritim = MIN(N_max,n_ccn_maritim)          !..Limit                        

      RETURN

    END FUNCTION n_ccn_maritim

    DOUBLE PRECISION FUNCTION n_ccn_texas()

      IMPLICIT NONE

      DOUBLE PRECISION :: S

      DOUBLE PRECISION, PARAMETER :: N_a = 1260.0d06
      DOUBLE PRECISION, PARAMETER :: k_a = 0.308

      !..CCN-Spektrum nach Khain et al. (2001), Rosenfeld & Woodley (2000)
      !  (continental case, Midland/Texas)

      S = MIN(1.0d1,1.0d2*s_w(i,j,k))                   !..Uebersaettigung in Prozent

      n_ccn_texas = N_a * S**k_a                         !..Power-Law

      RETURN

    END FUNCTION n_ccn_texas

    DOUBLE PRECISION FUNCTION n_ccn_polluted()

      IMPLICIT NONE

      DOUBLE PRECISION :: S

      DOUBLE PRECISION, PARAMETER :: N_a = 986.0d06

      DOUBLE PRECISION, PARAMETER :: C   = 1865.0d6
      DOUBLE PRECISION, PARAMETER :: k_a  = 0.86
      DOUBLE PRECISION, PARAMETER :: A   = 6.947


      !..CCN-Spektrum nach Cohard et al. (1998)
      !  (polluted case, experimental Hudson-Li, siehe Tabelle 2)

      S = MIN(1.0d1,1.0d2*s_w(i,j,k))                   !..Uebersaettigung in Prozent

      n_ccn_polluted = C/A**k_a * LOG( A**k_a * S**k_a + 1.0)  !..Approximation der Hypergeo-Funk.

      ! n_ccn_polluted = N_a * s_w**0.64                !..Power-Law

      RETURN

    END FUNCTION n_ccn_polluted

    DOUBLE PRECISION FUNCTION n_ccn_synthetic_maritime()

      IMPLICIT NONE

      DOUBLE PRECISION :: S

      DOUBLE PRECISION, PARAMETER :: C   = 1.93d14
      DOUBLE PRECISION, PARAMETER :: k_a = 4.16
      DOUBLE PRECISION, PARAMETER :: A   = 62.82


      !..CCN-Spektrum nach Cohard et al. (1998)
      !  (synthetic maritime case)

      S = MIN(1.0d1,1.0d2*s_w(i,j,k))                !..Uebersaettigung in Prozent

      !..Approximation der Hypergeo-Funk.
      n_ccn_synthetic_maritime = C/A**k_a * LOG( A**k_a * S**k_a + 1.0)

      RETURN

    END FUNCTION n_ccn_synthetic_maritime

    DOUBLE PRECISION FUNCTION n_ccn_synthetic_continental()

      IMPLICIT NONE

      DOUBLE PRECISION :: S

      DOUBLE PRECISION, PARAMETER :: C   = 3.27d9
      DOUBLE PRECISION, PARAMETER :: k_a = 1.56
      DOUBLE PRECISION, PARAMETER :: A   = 8.904


      !..CCN-Spektrum nach Cohard et al. (1998)
      !  (synthetic continental case)

      S = MIN(1.0d1,1.0d2*s_w(i,j,k))                !..Uebersaettigung in Prozent

      !..Approximation der Hypergeo-Funk.
      n_ccn_synthetic_continental = C/A**k_a * LOG( A**k_a * S**k_a + 1.0)

      RETURN

    END FUNCTION n_ccn_synthetic_continental

  END SUBROUTINE cloud_nucleation

  SUBROUTINE cloud_nucleation_SK()
    !*******************************************************************************
    !                                                                              *
    !       Calculation of cloud droplet nucleation                                *
    !       using the look-up tables by Segal and Khain 2006 (JGR, vol.11)         *
    !       (implemnented by Heike Noppel, IMK)                                    *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, q_cloud, n_cloud, &
         &                        T_g, p_g, rho_g, dt,dz3d, w_cb, dqdt, zml_k,       &
                                  speichere_dqdt, cond_neu_sb,           &
                                  &evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: r,cv,cp,dz 
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0

    IMPLICIT NONE

    ! Locale Variablen 
    INTEGER, PARAMETER :: n_ncn=8, n_r2=3, n_lsigs=5, n_wcb=4
    INTEGER            :: i_lsigs, i_R2
    DOUBLE PRECISION   :: T_a             !..Absolute Temperatur
    DOUBLE PRECISION   :: p_a,q_d,x_d,e_d,n_c,q_c,rho_a
    DOUBLE PRECISION   :: nuc_n, nuc_q, S, dS, dSdz_w, nucn_max
    DOUBLE PRECISION   :: a_ld,a_dl,e_sw,q_sw
    DOUBLE PRECISION   :: N_cn, N_cn0, lsigs, Ndrop, R2, wcb, etas
    DOUBLE PRECISION   :: z0_nccn, z1e_nccn
    DOUBLE PRECISION   :: tab_Ncn(n_ncn),         & ! look-up-table for Ncn 
                          tab_R2(n_r2),           & ! look-up_tbale for R2
                          tab_lsigs(n_lsigs),     & ! look-up-table for log(sigma_s)
                          tab_wcb(n_wcb),         & ! look-up-table for w at cloud base
                          tab_Ndrop(n_wcb,n_ncn), & ! number of cloud droplets in look_up-Table
                          tab_Ndrop_i(n_wcb)        ! number of cloud droplets in look_up-Table, 
                                                    !           interpolated with respect to Ncn

    INTEGER                     :: i,j,k,nuc_typ,stat_var

    nuc_typ = nuc_c_typ

    tab_R2  = (/0.02d0, 0.03d0, 0.04d0/)     ! in 10^(-6) m
    tab_wcb = (/0.5d0, 1.0d0, 2.5d0, 5.0d0/)
    tab_Ncn = (/50.d06, 100.d06, 200.d06, 400.d06, 800.d06, 1600.d06, 3200.d06, 6400.d06/) 
    tab_lsigs = (/0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0/)

    ! ATTENTION: At the moment only the values given above can be chosen for R2, 
    ! and lsigs (see below).
    ! Only wcb and N_cn0 can be interpolated. For the others this possibility is still missing.
    ! For high N_cn0 and small values of wcb a kind of "saturation" for cloud droplet 
    ! nucleation was assumed, because the original look-up tables don't give these values. 
    ! This assumption might be wrong!!! (comment by Heike Noppel)
         
    IF(isIO() .AND. isdebug) THEN
       WRITE (6, *) "CLOUD cloud_nucleation_SK: nuc_typ = ",nuc_typ
    ENDIF

    !..parameter for exponential decrease of N_ccn with height:
    !  1) up to this height (m) constant unchanged value:
    z0_nccn = 4000.0d0
    !  2)  height interval at which N_ccn decreses by factor 1/e above z0_nccn:
    z1e_nccn = 2000.0d0
!AS10080428>
    !z1e_nccn = 5000.0d0
!<AS10080428

    ! characteristics of different kinds of CN

    SELECT CASE(nuc_typ)
    CASE(6) 
       !... maritime case
       N_cn0 = 100.0d06   ! CN concentration at ground
       lsigs = 0.4d0      ! log(sigma_s)
       R2    = 0.03d0     ! in mum
       etas  = 0.9        ! soluble fraction
    CASE(7) 
       !... intermediate case
       N_cn0 = 500.0d06
       lsigs = 0.4d0
       R2    = 0.03d0       ! in mum
       etas  = 0.8          ! soluble fraction
    CASE(8)
       !... continental case
       N_cn0 = 1700.0d06
       lsigs = 0.2d0
       R2    = 0.03d0       ! in mum
       etas  = 0.7          ! soluble fraction
    CASE(9) 
       !... "polluted" continental 
       N_cn0 = 3200.0d06
       lsigs = 0.2d0
       R2    = 0.03d0       ! in mum
       etas  = 0.7          ! soluble fraction 
    CASE DEFAULT
      IF (isIO()) THEN
        WRITE(*,*) "ERROR: only the values 6,7,8,9 are allowed in cloud_nucleation_SK"
        WRITE(*,*) "       Actual nuc_type = " ,nuc_typ
      END IF
      CALL abortparallel ("Invalid value for nuc_typ in cloud_nucleation_sk", 1)
    END SELECT  
    !CASE(4)
    !   !... "seeded" continental (seeding with large hygroscopic particles)
    !   N_cn0 = 1800.0d06
    !   lsigs = 0.3d0
    !   R2    = 0.04d0       ! in mum
    !   etas  = 0.7          ! soluble fraction    
    !CASE(5) 
    !   !... "seeded" maritime (seeding with small hygroscopic particles) 
    !   N_cn0 = 200.0d06
    !   lsigs = 0.3d0
    !   R2    = 0.02d0       ! in mum
    !   etas  = 0.7          ! soluble fraction     
    !END SELECT

    i_lsigs = 0
    i_R2   = 0
    DO k=1,n_lsigs
      IF (lsigs==tab_lsigs(k)) THEN
         i_lsigs=k
         EXIT
      ENDIF
    END DO
    DO k=1,n_r2
      IF (R2==tab_R2(k)) THEN
         i_R2=k
         EXIT
      ENDIF
    END DO
    IF (i_lsigs==0) CALL abortparallel ("!!! INVALID VALUE FOR LSIGS !!! ", 1)
    IF (i_R2==0) CALL abortparallel ("!!! INVALID VALUE FOR R2 !!! ", 1)
  
    CALL lookuptable(tab_Ndrop,i_lsigs,i_R2)

!   Nucleation    
    nucn_max=0.d0
 
    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

!          IF (w_cb(i,j,k) <= 0.00d0 .or. T_0(i,j,k) < 263.0d0 ) CYCLE

          nuc_q = 0.0d0
          nuc_n = 0.d0
          n_c   = n_cloud(i,j,k)
          q_c   = q_cloud(i,j,k)
          rho_a = rho_0(i,j,k)
          wcb   = w_cb(i,j,k)

!AS20080419> hard upper limit for number conc that
!            eliminates also unrealistic high value
!            that would come from the dynamical core

          n_cloud(i,j,k) = MIN(n_cloud(i,j,k),N_cn0)

          IF (w_cb(i,j,k) <= 0.00d0) CYCLE
!<AS20080419

          ! --- N_cn depends on height (to avoid strong in-cloud nucleation)--
          IF(zml_k(i,j,k) > z0_nccn) THEN
             N_cn = N_cn0*EXP((z0_nccn-zml_k(i,j,k))/z1e_nccn)  ! exponential decrease with height
          ELSE
             N_cn = N_cn0
          ENDIF 
! UB_20090316:       N_cn = N_cn0

          ! Interpolation of the look-up tables with respect to Ncn
          tab_Ndrop_i = ip_ndrop_ncn(tab_Ndrop,tab_Ncn,N_cn)
          ! interpol. with respect to wcb = ip_ndrop_wcb          
          nuc_n = etas * ip_ndrop_wcb(tab_Ndrop_i,tab_wcb,wcb) - n_c

          nuc_n = MAX(nuc_n,0.0d0)

          nuc_q = MIN(nuc_n * cloud%x_min, q(i,j,k))
          nuc_n = nuc_q / cloud%x_min

!AS20080409> hard upper limit for nucleation rate
          nuc_n = MIN(nuc_n,N_cn-n_c)
!<AS20080409

          n_cloud(i,j,k) = n_cloud(i,j,k) + nuc_n
          q_cloud(i,j,k) = q_cloud(i,j,k) + nuc_q
          q(i,j,k)       = q(i,j,k)       - nuc_q


          ! nucn_max=max(nuc_n,nucn_max)     

#ifdef SAVE_CONVERSIONRATES
          IF (speichere_dqdt) THEN
            dqdt(i,j,k,1) = nuc_q
          END IF
#endif
          IF (speichere_precipstat) THEN
            cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + nuc_q
          END IF

       END DO
    END DO
  END DO

  CONTAINS
     
    SUBROUTINE lookuptable(tab_ndrop,i_lsigs,i_R2)

      IMPLICIT NONE

      INTEGER, PARAMETER ::  n_wcb = 4, n_ncn = 8
      INTEGER            ::  i_lsigs, i_R2
      DOUBLE PRECISION   ::  &
        ndrop1_11(n_ncn), ndrop1_12(n_ncn), ndrop1_13(n_ncn), ndrop1_14(n_ncn), ndrop1_15(n_ncn), &
        ndrop1_21(n_ncn), ndrop1_22(n_ncn), ndrop1_23(n_ncn), ndrop1_24(n_ncn), ndrop1_25(n_ncn), & 
        ndrop1_31(n_ncn), ndrop1_32(n_ncn), ndrop1_33(n_ncn), ndrop1_34(n_ncn), ndrop1_35(n_ncn), & 
        ndrop1_41(n_ncn), ndrop1_42(n_ncn), ndrop1_43(n_ncn), ndrop1_44(n_ncn), ndrop1_45(n_ncn), &
        ndrop2_11(n_ncn), ndrop2_12(n_ncn), ndrop2_13(n_ncn), ndrop2_14(n_ncn), ndrop2_15(n_ncn), &
        ndrop2_21(n_ncn), ndrop2_22(n_ncn), ndrop2_23(n_ncn), ndrop2_24(n_ncn), ndrop2_25(n_ncn), & 
        ndrop2_31(n_ncn), ndrop2_32(n_ncn), ndrop2_33(n_ncn), ndrop2_34(n_ncn), ndrop2_35(n_ncn), & 
        ndrop2_41(n_ncn), ndrop2_42(n_ncn), ndrop2_43(n_ncn), ndrop2_44(n_ncn), ndrop2_45(n_ncn), &
        ndrop3_11(n_ncn), ndrop3_12(n_ncn), ndrop3_13(n_ncn), ndrop3_14(n_ncn), ndrop3_15(n_ncn), &
        ndrop3_21(n_ncn), ndrop3_22(n_ncn), ndrop3_23(n_ncn), ndrop3_24(n_ncn), ndrop3_25(n_ncn), & 
        ndrop3_31(n_ncn), ndrop3_32(n_ncn), ndrop3_33(n_ncn), ndrop3_34(n_ncn), ndrop3_35(n_ncn), & 
        ndrop3_41(n_ncn), ndrop3_42(n_ncn), ndrop3_43(n_ncn), ndrop3_44(n_ncn), ndrop3_45(n_ncn), &
        tab_ndrop(n_wcb,n_ncn)
 
      ! look up tables
      ! Ncn              50       100       200       400       800       1600      3200      6400
      ! table4a (R2=0.02mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 "extrapolated")
      ndrop1_11 =  (/  42.2d06,  70.2d06, 112.2d06, 173.1d06, 263.7d06, 397.5d06, 397.5d06, 397.5d06/)
      ndrop1_12 =  (/  35.5d06,  60.1d06, 100.0d06, 163.9d06, 264.5d06, 418.4d06, 418.4d06, 418.4d06/)
      ndrop1_13 =  (/  32.6d06,  56.3d06,  96.7d06, 163.9d06, 272.0d06, 438.5d06, 438.5d06, 438.5d06/)
      ndrop1_14 =  (/  30.9d06,  54.4d06,  94.6d06, 162.4d06, 271.9d06, 433.5d06, 433.5d06, 433.5d06/)
      ndrop1_15 =  (/  29.4d06,  51.9d06,  89.9d06, 150.6d06, 236.5d06, 364.4d06, 364.4d06, 364.4d06/)
      ! table4b (R2=0.02mum, wcb=1.0m/s) (for Ncn=50 "interpolted" and Ncn=6400 extrapolated)
      ndrop1_21 =  (/  45.3d06,  91.5d06, 158.7d06, 264.4d06, 423.1d06, 672.5d06, 397.5d06, 397.5d06/)
      ndrop1_22 =  (/  38.5d06,  77.1d06, 133.0d06, 224.9d06, 376.5d06, 615.7d06, 418.4d06, 418.4d06/)
      ndrop1_23 =  (/  35.0d06,  70.0d06, 122.5d06, 212.0d06, 362.1d06, 605.3d06, 438.5d06, 438.5d06/)
      ndrop1_24 =  (/  32.4d06,  65.8d06, 116.4d06, 204.0d06, 350.6d06, 584.4d06, 433.5d06, 433.5d06/)
      ndrop1_25 =  (/  31.2d06,  62.3d06, 110.1d06, 191.3d06, 320.6d06, 501.3d06, 364.4d06, 364.4d06/)
      ! table4c (R2=0.02mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 "interpolated")
      ndrop1_31 =  (/  50.3d06, 100.5d06, 201.1d06, 373.1d06, 664.7d06,1132.8d06,1876.8d06,2973.7d06/)
      ndrop1_32 =  (/  44.1d06,  88.1d06, 176.2d06, 314.0d06, 546.9d06, 941.4d06,1579.2d06,2542.2d06/)
      ndrop1_33 =  (/  39.7d06,  79.5d06, 158.9d06, 283.4d06, 498.9d06, 865.9d06,1462.6d06,2355.8d06/)
      ndrop1_34 =  (/  37.0d06,  74.0d06, 148.0d06, 264.6d06, 468.3d06, 813.3d06,1371.3d06,2137.2d06/)
      ndrop1_35 =  (/  34.7d06,  69.4d06, 138.8d06, 246.9d06, 432.9d06, 737.8d06,1176.7d06,1733.0d06/)
      ! table4d (R2=0.02mum, wcb=5.0m/s) (for Ncn=50,100,200 "interpolated")
      ndrop1_41 =  (/  51.5d06, 103.1d06, 206.1d06, 412.2d06, 788.1d06,1453.1d06,2585.1d06,4382.5d06/)
      ndrop1_42 =  (/  46.6d06,  93.2d06, 186.3d06, 372.6d06, 657.2d06,1202.8d06,2098.0d06,3556.9d06/)
      ndrop1_43 =  (/  70.0d06,  70.0d06, 168.8d06, 337.6d06, 606.7d06,1078.5d06,1889.0d06,3206.9d06/)
      ndrop1_44 =  (/  42.2d06,  84.4d06, 166.4d06, 312.7d06, 562.2d06,1000.3d06,1741.1d06,2910.1d06/)
      ndrop1_45 =  (/  36.5d06,  72.9d06, 145.8d06, 291.6d06, 521.0d06, 961.1d06,1551.1d06,2444.6d06/)
      ! table5a (R2=0.03mum, wcb=0.5m/s) 
      ndrop2_11 =  (/  50.0d06,  95.8d06, 176.2d06, 321.6d06, 562.3d06, 835.5d06, 835.5d06, 835.5d06/)
      ndrop2_12 =  (/  44.7d06,  81.4d06, 144.5d06, 251.5d06, 422.7d06, 677.8d06, 677.8d06, 677.8d06/)
      ndrop2_13 =  (/  40.2d06,  72.8d06, 129.3d06, 225.9d06, 379.9d06, 606.5d06, 606.5d06, 606.5d06/)
      ndrop2_14 =  (/  37.2d06,  67.1d06, 119.5d06, 206.7d06, 340.5d06, 549.4d06, 549.4d06, 549.4d06/)
      ndrop2_15 =  (/  33.6d06,  59.0d06,  99.4d06, 150.3d06, 251.8d06, 466.0d06, 466.0d06, 466.0d06/)
      ! table5b (R2=0.03mum, wcb=1.0m/s) (Ncn=50 "interpolated", Ncn=6400 "extrapolated)
      ndrop2_21 =  (/  50.7d06, 101.4d06, 197.6d06, 357.2d06, 686.6d06,1186.4d06,1892.2d06,1892.2d06/)
      ndrop2_22 =  (/  46.6d06,  93.3d06, 172.2d06, 312.1d06, 550.7d06, 931.6d06,1476.6d06,1476.6d06/)
      ndrop2_23 =  (/  42.2d06,  84.4d06, 154.0d06, 276.3d06, 485.6d06, 811.2d06,1271.7d06,1271.7d06/)
      ndrop2_24 =  (/  39.0d06,  77.9d06, 141.2d06, 251.8d06, 436.7d06, 708.7d06,1117.7d06,1117.7d06/)
      ndrop2_25 =  (/  35.0d06,  70.1d06, 123.9d06, 210.2d06, 329.9d06, 511.9d06, 933.4d06, 933.4d06/)
      ! table5c (R2=0.03mum, wcb=2.5m/s)
      ndrop2_31 =  (/  51.5d06, 103.0d06, 205.9d06, 406.3d06, 796.4d06,1524.0d06,2781.4d06,4609.3d06/)
      ndrop2_32 =  (/  49.6d06,  99.1d06, 198.2d06, 375.5d06, 698.3d06,1264.1d06,2202.8d06,3503.6d06/)
      ndrop2_33 =  (/  45.8d06,  91.6d06, 183.2d06, 339.5d06, 618.9d06,1105.2d06,1881.8d06,2930.9d06/)
      ndrop2_34 =  (/  42.3d06,  84.7d06, 169.3d06, 310.3d06, 559.5d06, 981.7d06,1611.6d06,2455.6d06/)
      ndrop2_35 =  (/  38.2d06,  76.4d06, 152.8d06, 237.3d06, 473.3d06, 773.1d06,1167.9d06,1935.0d06/)
      ! table5d (R2=0.03mum, wcb=5.0m/s)
      ndrop2_41 =  (/  51.9d06, 103.8d06, 207.6d06, 415.1d06, 819.6d06,1616.4d06,3148.2d06,5787.9d06/)
      ndrop2_42 =  (/  50.7d06, 101.5d06, 203.0d06, 405.9d06, 777.0d06,1463.8d06,2682.6d06,4683.0d06/)
      ndrop2_43 =  (/  47.4d06,  94.9d06, 189.7d06, 379.4d06, 708.7d06,1301.3d06,2334.3d06,3951.8d06/)
      ndrop2_44 =  (/  44.0d06,  88.1d06, 176.2d06, 352.3d06, 647.8d06,1173.0d06,2049.7d06,3315.6d06/)
      ndrop2_45 =  (/  39.7d06,  79.4d06, 158.8d06, 317.6d06, 569.5d06, 988.5d06,1615.6d06,2430.3d06/)
      ! table6a (R2=0.04mum, wcb=0.5m/s)
      ndrop3_11 =  (/  50.6d06, 100.3d06, 196.5d06, 374.7d06, 677.3d06,1138.9d06,1138.9d06,1138.9d06/)
      ndrop3_12 =  (/  48.4d06,  91.9d06, 170.6d06, 306.9d06, 529.2d06, 862.4d06, 862.4d06, 862.4d06/)
      ndrop3_13 =  (/  44.4d06,  82.5d06, 150.3d06, 266.4d06, 448.0d06, 740.7d06, 740.7d06, 740.7d06/)
      ndrop3_14 =  (/  40.9d06,  75.0d06, 134.7d06, 231.9d06, 382.1d06, 657.6d06, 657.6d06, 657.6d06/)
      ndrop3_15 =  (/  34.7d06,  59.3d06,  93.5d06, 156.8d06, 301.9d06, 603.8d06, 603.8d06, 603.8d06/)
      ! table6b (R2=0.04mum, wcb=1.0m/s)
      ndrop3_21 =  (/  50.9d06, 101.7d06, 201.8d06, 398.8d06, 773.7d06,1420.8d06,2411.8d06,2411.8d06/)
      ndrop3_22 =  (/  49.4d06,  98.9d06, 189.7d06, 356.2d06, 649.5d06,1117.9d06,1805.2d06,1805.2d06/)
      ndrop3_23 =  (/  45.6d06,  91.8d06, 171.5d06, 214.9d06, 559.0d06, 932.8d06,1501.6d06,1501.6d06/)
      ndrop3_24 =  (/  42.4d06,  84.7d06, 155.8d06, 280.5d06, 481.9d06, 779.0d06,1321.9d06,1321.9d06/)
      ndrop3_25 =  (/  36.1d06,  72.1d06, 124.4d06, 198.4d06, 319.1d06, 603.8d06,1207.6d06,1207.6d06/)
      ! table6c (R2=0.04mum, wcb=2.5m/s)
      ndrop3_31 =  (/  51.4d06, 102.8d06, 205.7d06, 406.9d06, 807.6d06,1597.5d06,3072.2d06,5393.9d06/)
      ndrop3_32 =  (/  50.8d06, 101.8d06, 203.6d06, 396.0d06, 760.4d06,1422.1d06,2517.4d06,4062.8d06/)
      ndrop3_33 =  (/  48.2d06,  96.4d06, 193.8d06, 367.3d06, 684.0d06,1238.3d06,2087.3d06,3287.1d06/)
      ndrop3_34 =  (/  45.2d06,  90.4d06, 180.8d06, 335.7d06, 611.2d06,1066.3d06,1713.4d06,2780.3d06/)
      ndrop3_35 =  (/  38.9d06,  77.8d06, 155.5d06, 273.7d06, 455.2d06, 702.2d06,1230.7d06,2453.7d06/)
      ! table6d (R2=0.04mum, wcb=5.0m/s)
      ndrop3_41 =  (/  53.1d06, 106.2d06, 212.3d06, 414.6d06, 818.3d06,1622.2d06,3216.8d06,6243.9d06/)
      ndrop3_42 =  (/  51.6d06, 103.2d06, 206.3d06, 412.5d06, 805.3d06,1557.4d06,2940.4d06,5210.1d06/)
      ndrop3_43 =  (/  49.6d06,  99.2d06, 198.4d06, 396.7d06, 755.5d06,1414.5d06,2565.3d06,4288.1d06/)
      ndrop3_44 =  (/  46.5d06,  93.0d06, 186.0d06, 371.9d06, 692.9d06,1262.0d06,2188.3d06,3461.2d06/)
      ndrop3_45 =  (/  39.9d06,  79.9d06, 159.7d06, 319.4d06, 561.7d06, 953.9d06,1493.9d06,2464.7d06/)

      SELECT CASE (i_lsigs)
      CASE(1)
          SELECT CASE(i_R2)
          CASE(1)
            tab_Ndrop(1,:)=ndrop1_11
            tab_Ndrop(2,:)=ndrop1_21
            tab_Ndrop(3,:)=ndrop1_31
            tab_Ndrop(4,:)=ndrop1_41
          CASE(2)
            tab_Ndrop(1,:)=ndrop2_11
            tab_Ndrop(2,:)=ndrop2_21
            tab_Ndrop(3,:)=ndrop2_31
            tab_Ndrop(4,:)=ndrop2_41
          CASE(3)
            tab_Ndrop(1,:)=ndrop3_11
            tab_Ndrop(2,:)=ndrop3_21
            tab_Ndrop(3,:)=ndrop3_31
            tab_Ndrop(4,:)=ndrop3_41
          CASE DEFAULT
            WRITE(*,*) "!!!! wrong value for R2 in cloud_nucleation_SK !!!!!"
          END SELECT
      CASE(2)
          SELECT CASE(i_R2)
          CASE(1)
            tab_Ndrop(1,:)=ndrop1_12
            tab_Ndrop(2,:)=ndrop1_22
            tab_Ndrop(3,:)=ndrop1_32
            tab_Ndrop(4,:)=ndrop1_42
          CASE(2)
            tab_Ndrop(1,:)=ndrop2_12
            tab_Ndrop(2,:)=ndrop2_22
            tab_Ndrop(3,:)=ndrop2_32
            tab_Ndrop(4,:)=ndrop2_42
          CASE(3)
            tab_Ndrop(1,:)=ndrop3_12
            tab_Ndrop(2,:)=ndrop3_22
            tab_Ndrop(3,:)=ndrop3_32
            tab_Ndrop(4,:)=ndrop3_42
          CASE DEFAULT
            WRITE(*,*) "!!!! wrong value for R2 in cloud_nucleation_SK !!!!!"
          END SELECT
      CASE(3)
          SELECT CASE(i_R2)
          CASE(1)
            tab_Ndrop(1,:)=ndrop1_13
            tab_Ndrop(2,:)=ndrop1_23
            tab_Ndrop(3,:)=ndrop1_33
            tab_Ndrop(4,:)=ndrop1_43
          CASE(2)
            tab_Ndrop(1,:)=ndrop2_13
            tab_Ndrop(2,:)=ndrop2_23
            tab_Ndrop(3,:)=ndrop2_33
            tab_Ndrop(4,:)=ndrop2_43
          CASE(3)
            tab_Ndrop(1,:)=ndrop3_13
            tab_Ndrop(2,:)=ndrop3_23
            tab_Ndrop(3,:)=ndrop3_33
            tab_Ndrop(4,:)=ndrop3_43
          CASE DEFAULT
            WRITE(*,*) "!!!! wrong value for R2 in cloud_nucleation_SK !!!!!"
          END SELECT
      CASE(4)
          SELECT CASE(i_R2)
          CASE(1)
            tab_Ndrop(1,:)=ndrop1_14
            tab_Ndrop(2,:)=ndrop1_24
            tab_Ndrop(3,:)=ndrop1_34
            tab_Ndrop(4,:)=ndrop1_44
          CASE(2)
            tab_Ndrop(1,:)=ndrop2_14
            tab_Ndrop(2,:)=ndrop2_24
            tab_Ndrop(3,:)=ndrop2_34
            tab_Ndrop(4,:)=ndrop2_44
          CASE(3)
            tab_Ndrop(1,:)=ndrop3_14
            tab_Ndrop(2,:)=ndrop3_24
            tab_Ndrop(3,:)=ndrop3_34
            tab_Ndrop(4,:)=ndrop3_44
          CASE DEFAULT
            WRITE(*,*) "!!!! wrong value for R2 in cloud_nucleation_SK !!!!!"
          END SELECT
      CASE(5)
          SELECT CASE(i_R2)
          CASE(1)
            tab_Ndrop(1,:)=ndrop1_15
            tab_Ndrop(2,:)=ndrop1_25
            tab_Ndrop(3,:)=ndrop1_35
            tab_Ndrop(4,:)=ndrop1_45
          CASE(2)
            tab_Ndrop(1,:)=ndrop2_15
            tab_Ndrop(2,:)=ndrop2_25
            tab_Ndrop(3,:)=ndrop2_35
            tab_Ndrop(4,:)=ndrop2_45
          CASE(3)
            tab_Ndrop(1,:)=ndrop3_15
            tab_Ndrop(2,:)=ndrop3_25
            tab_Ndrop(3,:)=ndrop3_35
            tab_Ndrop(4,:)=ndrop3_45
          CASE DEFAULT
            WRITE(*,*) "!!!! wrong value for R2 in cloud_nucleation_SK !!!!!"
          END SELECT
      CASE DEFAULT
        WRITE(*,*)  "!!!! wrong value for lsigs in cloud_nucleation_SK !!!!!"
      END SELECT
            
      RETURN

    END SUBROUTINE lookuptable
     
    !--------------

    FUNCTION ip_ndrop_ncn(tab_ndrop,tab_ncn,Ncn)

      ! Interpolation of the look-up table with respect to aerosol concentration Ncn

      IMPLICIT NONE

      INTEGER, PARAMETER ::  n_wcb = 4, n_ncn=8
      INTEGER            ::  ki
      DOUBLE PRECISION   ::  ip_ndrop_ncn(n_wcb)
      DOUBLE PRECISION   ::  tab_ncn(n_ncn), tab_ndrop(n_wcb,n_ncn), Ncn
      LOGICAL :: found
  
      !... Interpolation of Ndrop from the values of Ncn in the look-up-table to given Ncn
            
      found = .FALSE.
      IF (Ncn <= tab_ncn(1)) THEN
        ip_ndrop_ncn = tab_ndrop(:,1) 
        found = .TRUE.
      ELSE IF (Ncn >= tab_ncn(n_ncn)) THEN
        ip_ndrop_ncn = tab_ndrop(:,n_ncn)      
        found = .TRUE.
      ELSE
        DO ki = 1,n_ncn-1
          IF (Ncn >= tab_ncn(ki) .AND.  Ncn <= tab_ncn(ki+1)) THEN
            ip_ndrop_ncn(:) = tab_ndrop(:,ki) + &
                 (tab_ndrop(:,ki+1)-tab_ndrop(:,ki)) / (tab_ncn(ki+1)-tab_ncn(ki)) * (Ncn-tab_ncn(ki))
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF (.NOT.found) THEN
        WRITE (*,*) "CLOUD_NUCLEATION ip_ndrop_ncn: Lookup table interpolation failed!" 
        STOP
      END IF

      RETURN

    END FUNCTION  ip_ndrop_ncn

    !------------------

    DOUBLE PRECISION FUNCTION ip_ndrop_wcb(tab_Ndrop_i,tab_wcb,wcb)

      ! Interpolation of the interpolated look-up table with respect to w_cb

      IMPLICIT NONE

      INTEGER, PARAMETER ::  n_wcb = 4
      INTEGER            ::  ki
      DOUBLE PRECISION   ::  tab_wcb(1:n_wcb), tab_Ndrop_i(1:n_wcb), wcb
      LOGICAL :: found

      !... Interpolation for Ndrop from the values of wcb in the look-up-tabel to detected wcb
            
      found = .FALSE.
      IF (wcb <= tab_wcb(1)) THEN
        ip_ndrop_wcb = tab_ndrop_i(1) / tab_wcb(1) * wcb
        found = .TRUE.
      ELSE IF (wcb  >= tab_wcb(n_wcb)) THEN
        ip_ndrop_wcb = tab_ndrop_i(n_wcb)      
        found = .TRUE.
      ELSE
        DO ki = 1,n_wcb-1
          IF (wcb >= tab_wcb(ki) .AND.  wcb <= tab_wcb(ki+1)) THEN
            ip_ndrop_wcb = tab_ndrop_i(ki) + &
              (tab_ndrop_i(ki+1)-tab_ndrop_i(ki)) / (tab_wcb(ki+1)-tab_wcb(ki)) * (wcb-tab_wcb(ki))
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF (.NOT.found) THEN
        WRITE (*,*) "CLOUD_NUCLEATION ip_ndrop_wcb: Lookup table interpolation failed!" 
        STOP
      END IF

      RETURN

    END FUNCTION  ip_ndrop_wcb

  END SUBROUTINE cloud_nucleation_SK

  SUBROUTINE clnuc_sk_4d()
    !*******************************************************************************
    !                                                                              *
    !       Calculation of cloud droplet nucleation                                *
    !       using the look-up tables by Segal and Khain 2006 (JGR, vol.11)         *
    !                                                                              *
    !       Difference to Heikes routine cloud_nucleation_SK()                     *
    !       Equidistant lookup table is used to enable better vectorization        *
    !       properties.                                                            *
    !                                                                              *
    !*******************************************************************************

    USE globale_variablen,  ONLY: loc_ix, loc_iy, loc_iz, T, p, q, q_cloud, n_cloud, &
         &                        T_g, p_g, rho_g, dt,dz3d, w_cb, dqdt, zml_k,       &
                                  speichere_dqdt, cond_neu_sb,           &
                                  &evap_neu_sb, speichere_precipstat
    USE konstanten,         ONLY: r,cv,cp,dz 
    USE parallele_umgebung, ONLY: isIO,abortparallel
    USE initialisierung,    ONLY: T_0,p_0,rho_0


    IMPLICIT NONE

    ! Locale Variablen

    ! grid sizes of the original table:
    INTEGER, PARAMETER       :: n_r2 = 3, n_lsigs = 5, n_ncn = 8 , n_wcb = 4
    ! desired grid sizes of the new equidistant table:
    INTEGER, PARAMETER       :: nr2  = 3, nlsigs  = 5, nncn  = 129, nwcb  = 11
    DOUBLE PRECISION         :: n_c,q_c
    DOUBLE PRECISION         :: nuc_n, nuc_q, nucn_max
    DOUBLE PRECISION         :: n_cn, n_cn0, lsigs, nccn, r2, wcb, etas
    DOUBLE PRECISION         :: r2_loc, lsigs_loc, ncn_loc, wcb_loc
    DOUBLE PRECISION         :: z0_nccn, z1e_nccn
    INTEGER                  :: i, j, k, nuc_typ
    INTEGER, SAVE            :: firstcall
    INTEGER                  :: iu, ju, ku, lu
    TYPE(lookupt_4D), SAVE   :: otab
    TYPE(lookupt_4D), SAVE   :: tab
    DOUBLE PRECISION         :: hilf1(2,2,2,2), hilf2(2,2,2), hilf3(2,2), hilf4(2)


    nuc_typ = nuc_c_typ


    IF (firstcall .NE. 1) THEN
      !========================================================
      ! original look-up-table from Segal and Khain
      ! get the original 4D table values:
      CALL get_otab()
      !========================================================
      ! construct the new equidistant table tab:
      CALL equi_table()
      firstcall = 1
      IF(isIO() .AND. isdebug) THEN
        WRITE (*,*) "CLNUC_SK_4D: Equidistant lookup table created "
      ENDIF
    END IF

    IF(isIO() .AND. isdebug) THEN
      WRITE (*,*) "CLNUC_SK_4D: nuc_typ = ",nuc_typ
    ENDIF

    !..parameter for exponential decrease of N_ccn with height:
    !  1) up to this height (m) constant unchanged value:
    z0_nccn = 4000.0d0
    !  2)  height interval at which N_ccn decreses by factor 1/e above z0_nccn:
    z1e_nccn = 2000.0d0


    !======================================================
    ! characteristics of different kinds of prototype CN:

    SELECT CASE(nuc_typ)
    CASE(6) 
      !... maritime case
      n_cn0 = 100.0d06   ! CN concentration at ground
      lsigs = 0.4d0      ! log(sigma_s)
      r2    = 0.03d0     ! in mum
      etas  = 0.9        ! soluble fraction
    CASE(7) 
      !... intermediate case
      n_cn0 = 500.0d06
      lsigs = 0.4d0
      r2    = 0.03d0       ! in mum
      etas  = 0.8          ! soluble fraction
    CASE(8)
      !... continental case
      n_cn0 = 1700.0d06
      lsigs = 0.2d0
      r2    = 0.03d0       ! in mum
      etas  = 0.7          ! soluble fraction
    CASE(9) 
      !... "polluted" continental 
      n_cn0 = 3200.0d06
      lsigs = 0.2d0
      r2    = 0.03d0       ! in mum
      etas  = 0.7          ! soluble fraction 
    CASE DEFAULT
      !      IF (isIO()) THEN
      WRITE(*,*) "ERROR: only the values 6,7,8,9 are allowed in cloud_nucleation_SK"
      WRITE(*,*) "       Actual nuc_type = " ,nuc_typ
      !      END IF
      CALL abortparallel ("Invalid value for nuc_typ in cloud_nucleation_sk", 1)
    END SELECT


!!!! Probably do the following in a future program version:
    ! Interpolation of the original 4D table to a 2D table. Interpolation is
    ! done bilinearily w.r.t. r2 and lsigs. In this way, it is assumed that those parameters
    ! do not change during the simulation. 


    ! Interpolate linearily w.r.t. all 4 parameters:
    nucn_max=0.d0

    DO k = 1, loc_iz
      DO j = 1, loc_iy
        DO i = 0, loc_ix

!          IF (w_cb(i,j,k) > 0.00d0 .AND. T_0(i,j,k) >= 263.0d0 ) THEN
          IF (w_cb(i,j,k) > 0.00d0) THEN

            nuc_q = 0.0d0
            nuc_n = 0.d0
            n_c   = n_cloud(i,j,k)
            q_c   = q_cloud(i,j,k)
            wcb   = w_cb(i,j,k)

!AS20080419> hard upper limit for number conc that
!            eliminates also unrealistic high value
!            that would come from the dynamical core

            n_cloud(i,j,k) = MIN(n_cloud(i,j,k),N_cn0)

            ! write(*,*) "cloud base detected", zml_k(i,j,k), w_cb(i,j,k)
            ! --- N_cn depends on height (to avoid strong in-cloud nucleation)--
            n_cn = n_cn0 * MIN( EXP((z0_nccn-zml_k(i,j,k))/z1e_nccn), 1.0d0)  ! exponential decrease with height

            ! Interpolation of the look-up tables with respect to all 4 parameters:
            ! (clip values outside range to the marginal values)
            r2_loc    = MIN(MAX(r2,     tab%x1(1)), tab%x1(tab%n1))
            iu = MIN(FLOOR((r2_loc -    tab%x1(1)) * tab%odx1 ) + 1, tab%n1-1)
            lsigs_loc = MIN(MAX(lsigs,  tab%x2(1)), tab%x2(tab%n2))
            ju = MIN(FLOOR((lsigs_loc - tab%x2(1)) * tab%odx2 ) + 1, tab%n2-1)
            ncn_loc   = MIN(MAX(n_cn,   tab%x3(1)), tab%x3(tab%n3))
            ku = MIN(FLOOR((ncn_loc -   tab%x3(1)) * tab%odx3 ) + 1, tab%n3-1)
            wcb_loc   = MIN(MAX(wcb,    tab%x4(1)), tab%x4(tab%n4))
            lu = MIN(FLOOR((wcb_loc -   tab%x4(1)) * tab%odx4 ) + 1, tab%n4-1)

            hilf1 = tab%ltable( iu:iu+1, ju:ju+1, ku:ku+1, lu:lu+1)
            hilf2 = hilf1(1,:,:,:) + (hilf1(2,:,:,:) - hilf1(1,:,:,:)) * tab%odx1 * ( r2_loc    - tab%x1(iu) )
            hilf3 = hilf2(1,:,:)   + (hilf2(2,:,:)   - hilf2(1,:,:)  ) * tab%odx2 * ( lsigs_loc - tab%x2(ju) )
            hilf4 = hilf3(1,:)     + (hilf3(2,:)     - hilf3(1,:)    ) * tab%odx3 * ( ncn_loc   - tab%x3(ku) )
            nccn  = hilf4(1)       + (hilf4(2)       - hilf4(1)      ) * tab%odx4 * ( wcb_loc   - tab%x4(lu) )

! HN_20080401>
            ! If n_cn is outside the range of the lookup table values, resulting 
            ! NCCN are clipped to the margin values. For the case of these margin values
            ! beeing larger than n_cn (which happens sometimes, unfortunately), limit NCCN by n_cn:
            nccn = MIN(nccn, n_cn)
! HN_20080401<

            nuc_n = etas * nccn - n_c
            nuc_n = MAX(nuc_n,0.0d0)

            nuc_q = MIN(nuc_n * cloud%x_min, q(i,j,k))
            nuc_n = nuc_q / cloud%x_min

            n_cloud(i,j,k) = n_cloud(i,j,k) + nuc_n
            q_cloud(i,j,k) = q_cloud(i,j,k) + nuc_q
            q(i,j,k)       = q(i,j,k)       - nuc_q
            nucn_max=MAX(nuc_n,nucn_max)     

#ifdef SAVE_CONVERSIONRATES
          IF (speichere_dqdt) THEN
            dqdt(i,j,k,1) = nuc_q
          END IF
#endif
          IF (speichere_precipstat) THEN
            cond_neu_sb(i,j,k) = cond_neu_sb(i,j,k) + nuc_q
          END IF

          END IF

       END DO
    END DO
  END DO

  CONTAINS

    SUBROUTINE get_otab()
      IMPLICIT NONE

      otab%n1 = n_r2
      otab%n2 = n_lsigs
      otab%n3 = n_ncn + 1
      otab%n4 = n_wcb + 1

      ALLOCATE( otab%x1(otab%n1) )
      ALLOCATE( otab%x2(otab%n2) )
      ALLOCATE( otab%x3(otab%n3) )
      ALLOCATE( otab%x4(otab%n4) )
      ALLOCATE( otab%ltable(otab%n1,otab%n2,otab%n3,otab%n4) )
 
      ! original (non-)equidistant table vectors:
      ! r2:
      otab%x1  = (/0.02d0, 0.03d0, 0.04d0/)     ! in 10^(-6) m
      ! lsigs:
      otab%x2  = (/0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0/)
      ! n_cn: (UB: um 0.0 m**-3 ergaenzt zur linearen Interpolation zw. 0.0 und 50e6 m**-3)
      otab%x3  = (/0.0d6, 50.d06, 100.d06, 200.d06, 400.d06, 800.d06, 1600.d06, 3200.d06, 6400.d06/) ! in m**-3
      ! wcb: (UB: um 0.0 m/s ergaenzt zur linearen Interpolation zw. 0.0 und 0.5 m/s)
      otab%x4  = (/0.0d0, 0.5d0, 1.0d0, 2.5d0, 5.0d0/)

      ! look up table for NCCN activated at given R2, lsigs, Ncn and wcb:

      ! Ncn              50       100       200       400       800       1600      3200      6400
      ! table4a (R2=0.02mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 "extrapolated")
      otab%ltable(1,1,2:otab%n3,2) =  (/  42.2d06,  70.2d06, 112.2d06, 173.1d06, 263.7d06, 397.5d06, 397.5d06, 397.5d06/)
      otab%ltable(1,2,2:otab%n3,2) =  (/  35.5d06,  60.1d06, 100.0d06, 163.9d06, 264.5d06, 418.4d06, 418.4d06, 418.4d06/)
      otab%ltable(1,3,2:otab%n3,2) =  (/  32.6d06,  56.3d06,  96.7d06, 163.9d06, 272.0d06, 438.5d06, 438.5d06, 438.5d06/)
      otab%ltable(1,4,2:otab%n3,2) =  (/  30.9d06,  54.4d06,  94.6d06, 162.4d06, 271.9d06, 433.5d06, 433.5d06, 433.5d06/)
      otab%ltable(1,5,2:otab%n3,2) =  (/  29.4d06,  51.9d06,  89.9d06, 150.6d06, 236.5d06, 364.4d06, 364.4d06, 364.4d06/)
      ! table4b (R2=0.02mum, wcb=1.0m/s) (for Ncn=50 "interpolted" and Ncn=6400 extrapolated)
      otab%ltable(1,1,2:otab%n3,3) =  (/  45.3d06,  91.5d06, 158.7d06, 264.4d06, 423.1d06, 672.5d06, 397.5d06, 397.5d06/)
      otab%ltable(1,2,2:otab%n3,3) =  (/  38.5d06,  77.1d06, 133.0d06, 224.9d06, 376.5d06, 615.7d06, 418.4d06, 418.4d06/)
      otab%ltable(1,3,2:otab%n3,3) =  (/  35.0d06,  70.0d06, 122.5d06, 212.0d06, 362.1d06, 605.3d06, 438.5d06, 438.5d06/)
      otab%ltable(1,4,2:otab%n3,3) =  (/  32.4d06,  65.8d06, 116.4d06, 204.0d06, 350.6d06, 584.4d06, 433.5d06, 433.5d06/)
      otab%ltable(1,5,2:otab%n3,3) =  (/  31.2d06,  62.3d06, 110.1d06, 191.3d06, 320.6d06, 501.3d06, 364.4d06, 364.4d06/)
      ! table4c (R2=0.02mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 "interpolated")
      otab%ltable(1,1,2:otab%n3,4) =  (/  50.3d06, 100.5d06, 201.1d06, 373.1d06, 664.7d06,1132.8d06,1876.8d06,2973.7d06/)
      otab%ltable(1,2,2:otab%n3,4) =  (/  44.1d06,  88.1d06, 176.2d06, 314.0d06, 546.9d06, 941.4d06,1579.2d06,2542.2d06/)
      otab%ltable(1,3,2:otab%n3,4) =  (/  39.7d06,  79.5d06, 158.9d06, 283.4d06, 498.9d06, 865.9d06,1462.6d06,2355.8d06/)
      otab%ltable(1,4,2:otab%n3,4) =  (/  37.0d06,  74.0d06, 148.0d06, 264.6d06, 468.3d06, 813.3d06,1371.3d06,2137.2d06/)
      otab%ltable(1,5,2:otab%n3,4) =  (/  34.7d06,  69.4d06, 138.8d06, 246.9d06, 432.9d06, 737.8d06,1176.7d06,1733.0d06/)
      ! table4d (R2=0.02mum, wcb=5.0m/s) (for Ncn=50,100,200 "interpolated")
      otab%ltable(1,1,2:otab%n3,5) =  (/  51.5d06, 103.1d06, 206.1d06, 412.2d06, 788.1d06,1453.1d06,2585.1d06,4382.5d06/)
      otab%ltable(1,2,2:otab%n3,5) =  (/  46.6d06,  93.2d06, 186.3d06, 372.6d06, 657.2d06,1202.8d06,2098.0d06,3556.9d06/)
      otab%ltable(1,3,2:otab%n3,5) =  (/  70.0d06,  70.0d06, 168.8d06, 337.6d06, 606.7d06,1078.5d06,1889.0d06,3206.9d06/)
      otab%ltable(1,4,2:otab%n3,5) =  (/  42.2d06,  84.4d06, 166.4d06, 312.7d06, 562.2d06,1000.3d06,1741.1d06,2910.1d06/)
      otab%ltable(1,5,2:otab%n3,5) =  (/  36.5d06,  72.9d06, 145.8d06, 291.6d06, 521.0d06, 961.1d06,1551.1d06,2444.6d06/)
      ! table5a (R2=0.03mum, wcb=0.5m/s) 
      otab%ltable(2,1,2:otab%n3,2) =  (/  50.0d06,  95.8d06, 176.2d06, 321.6d06, 562.3d06, 835.5d06, 835.5d06, 835.5d06/)
      otab%ltable(2,2,2:otab%n3,2) =  (/  44.7d06,  81.4d06, 144.5d06, 251.5d06, 422.7d06, 677.8d06, 677.8d06, 677.8d06/)
      otab%ltable(2,3,2:otab%n3,2) =  (/  40.2d06,  72.8d06, 129.3d06, 225.9d06, 379.9d06, 606.5d06, 606.5d06, 606.5d06/)
      otab%ltable(2,4,2:otab%n3,2) =  (/  37.2d06,  67.1d06, 119.5d06, 206.7d06, 340.5d06, 549.4d06, 549.4d06, 549.4d06/)
      otab%ltable(2,5,2:otab%n3,2) =  (/  33.6d06,  59.0d06,  99.4d06, 150.3d06, 251.8d06, 466.0d06, 466.0d06, 466.0d06/)
      ! table5b (R2=0.03mum, wcb=1.0m/s) (Ncn=50 "interpolated", Ncn=6400 "extrapolated)
      otab%ltable(2,1,2:otab%n3,3) =  (/  50.7d06, 101.4d06, 197.6d06, 357.2d06, 686.6d06,1186.4d06,1892.2d06,1892.2d06/)
      otab%ltable(2,2,2:otab%n3,3) =  (/  46.6d06,  93.3d06, 172.2d06, 312.1d06, 550.7d06, 931.6d06,1476.6d06,1476.6d06/)
      otab%ltable(2,3,2:otab%n3,3) =  (/  42.2d06,  84.4d06, 154.0d06, 276.3d06, 485.6d06, 811.2d06,1271.7d06,1271.7d06/)
      otab%ltable(2,4,2:otab%n3,3) =  (/  39.0d06,  77.9d06, 141.2d06, 251.8d06, 436.7d06, 708.7d06,1117.7d06,1117.7d06/)
      otab%ltable(2,5,2:otab%n3,3) =  (/  35.0d06,  70.1d06, 123.9d06, 210.2d06, 329.9d06, 511.9d06, 933.4d06, 933.4d06/)
      ! table5c (R2=0.03mum, wcb=2.5m/s)
      otab%ltable(2,1,2:otab%n3,4) =  (/  51.5d06, 103.0d06, 205.9d06, 406.3d06, 796.4d06,1524.0d06,2781.4d06,4609.3d06/)
      otab%ltable(2,2,2:otab%n3,4) =  (/  49.6d06,  99.1d06, 198.2d06, 375.5d06, 698.3d06,1264.1d06,2202.8d06,3503.6d06/)
      otab%ltable(2,3,2:otab%n3,4) =  (/  45.8d06,  91.6d06, 183.2d06, 339.5d06, 618.9d06,1105.2d06,1881.8d06,2930.9d06/)
      otab%ltable(2,4,2:otab%n3,4) =  (/  42.3d06,  84.7d06, 169.3d06, 310.3d06, 559.5d06, 981.7d06,1611.6d06,2455.6d06/)
      otab%ltable(2,5,2:otab%n3,4) =  (/  38.2d06,  76.4d06, 152.8d06, 237.3d06, 473.3d06, 773.1d06,1167.9d06,1935.0d06/)
      ! table5d (R2=0.03mum, wcb=5.0m/s)
      otab%ltable(2,1,2:otab%n3,5) =  (/  51.9d06, 103.8d06, 207.6d06, 415.1d06, 819.6d06,1616.4d06,3148.2d06,5787.9d06/)
      otab%ltable(2,2,2:otab%n3,5) =  (/  50.7d06, 101.5d06, 203.0d06, 405.9d06, 777.0d06,1463.8d06,2682.6d06,4683.0d06/)
      otab%ltable(2,3,2:otab%n3,5) =  (/  47.4d06,  94.9d06, 189.7d06, 379.4d06, 708.7d06,1301.3d06,2334.3d06,3951.8d06/)
      otab%ltable(2,4,2:otab%n3,5) =  (/  44.0d06,  88.1d06, 176.2d06, 352.3d06, 647.8d06,1173.0d06,2049.7d06,3315.6d06/)
      otab%ltable(2,5,2:otab%n3,5) =  (/  39.7d06,  79.4d06, 158.8d06, 317.6d06, 569.5d06, 988.5d06,1615.6d06,2430.3d06/)
      ! table6a (R2=0.04mum, wcb=0.5m/s)
      otab%ltable(3,1,2:otab%n3,2) =  (/  50.6d06, 100.3d06, 196.5d06, 374.7d06, 677.3d06,1138.9d06,1138.9d06,1138.9d06/)
      otab%ltable(3,2,2:otab%n3,2) =  (/  48.4d06,  91.9d06, 170.6d06, 306.9d06, 529.2d06, 862.4d06, 862.4d06, 862.4d06/)
      otab%ltable(3,3,2:otab%n3,2) =  (/  44.4d06,  82.5d06, 150.3d06, 266.4d06, 448.0d06, 740.7d06, 740.7d06, 740.7d06/)
      otab%ltable(3,4,2:otab%n3,2) =  (/  40.9d06,  75.0d06, 134.7d06, 231.9d06, 382.1d06, 657.6d06, 657.6d06, 657.6d06/)
      otab%ltable(3,5,2:otab%n3,2) =  (/  34.7d06,  59.3d06,  93.5d06, 156.8d06, 301.9d06, 603.8d06, 603.8d06, 603.8d06/)
      ! table6b (R2=0.04mum, wcb=1.0m/s)
      otab%ltable(3,1,2:otab%n3,3) =  (/  50.9d06, 101.7d06, 201.8d06, 398.8d06, 773.7d06,1420.8d06,2411.8d06,2411.8d06/)
      otab%ltable(3,2,2:otab%n3,3) =  (/  49.4d06,  98.9d06, 189.7d06, 356.2d06, 649.5d06,1117.9d06,1805.2d06,1805.2d06/)
      otab%ltable(3,3,2:otab%n3,3) =  (/  45.6d06,  91.8d06, 171.5d06, 214.9d06, 559.0d06, 932.8d06,1501.6d06,1501.6d06/)
      otab%ltable(3,4,2:otab%n3,3) =  (/  42.4d06,  84.7d06, 155.8d06, 280.5d06, 481.9d06, 779.0d06,1321.9d06,1321.9d06/)
      otab%ltable(3,5,2:otab%n3,3) =  (/  36.1d06,  72.1d06, 124.4d06, 198.4d06, 319.1d06, 603.8d06,1207.6d06,1207.6d06/)
      ! table6c (R2=0.04mum, wcb=2.5m/s)
      otab%ltable(3,1,2:otab%n3,4) =  (/  51.4d06, 102.8d06, 205.7d06, 406.9d06, 807.6d06,1597.5d06,3072.2d06,5393.9d06/)
      otab%ltable(3,2,2:otab%n3,4) =  (/  50.8d06, 101.8d06, 203.6d06, 396.0d06, 760.4d06,1422.1d06,2517.4d06,4062.8d06/)
      otab%ltable(3,3,2:otab%n3,4) =  (/  48.2d06,  96.4d06, 193.8d06, 367.3d06, 684.0d06,1238.3d06,2087.3d06,3287.1d06/)
      otab%ltable(3,4,2:otab%n3,4) =  (/  45.2d06,  90.4d06, 180.8d06, 335.7d06, 611.2d06,1066.3d06,1713.4d06,2780.3d06/)
      otab%ltable(3,5,2:otab%n3,4) =  (/  38.9d06,  77.8d06, 155.5d06, 273.7d06, 455.2d06, 702.2d06,1230.7d06,2453.7d06/)
      ! table6d (R2=0.04mum, wcb=5.0m/s)
      otab%ltable(3,1,2:otab%n3,5) =  (/  53.1d06, 106.2d06, 212.3d06, 414.6d06, 818.3d06,1622.2d06,3216.8d06,6243.9d06/)
      otab%ltable(3,2,2:otab%n3,5) =  (/  51.6d06, 103.2d06, 206.3d06, 412.5d06, 805.3d06,1557.4d06,2940.4d06,5210.1d06/)
      otab%ltable(3,3,2:otab%n3,5) =  (/  49.6d06,  99.2d06, 198.4d06, 396.7d06, 755.5d06,1414.5d06,2565.3d06,4288.1d06/)
      otab%ltable(3,4,2:otab%n3,5) =  (/  46.5d06,  93.0d06, 186.0d06, 371.9d06, 692.9d06,1262.0d06,2188.3d06,3461.2d06/)
      otab%ltable(3,5,2:otab%n3,5) =  (/  39.9d06,  79.9d06, 159.7d06, 319.4d06, 561.7d06, 953.9d06,1493.9d06,2464.7d06/)

      ! Additional values for wcb = 0.0 m/s, which are used for linear interpolation between
      ! wcb = 0.0 and 0.5 m/s. Values of 0.0 are reasonable here, because if no
      ! updraft is present, no new nucleation will take place:
      otab%ltable(:,:,:,1) = 0.0d0
      ! Additional values for n_cn = 0.0 m**-3, which are used for linear interpolation between
      ! n_cn = 0.0 and 50 m**-3. Values of 0.0 are reasonable, because if no aerosol
      ! particles are present, no nucleation will take place:
      otab%ltable(:,:,1,:) = 0.0d0

      !!! otab%dx1 ... otab%odx4 remain empty because this is a non-equidistant table.

      RETURN
    END SUBROUTINE get_otab
    
    SUBROUTINE equi_table()
      IMPLICIT NONE

      INTEGER :: i, j, k, l, ii, iu, ju,ku, lu
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iuv, juv, kuv, luv
      DOUBLE PRECISION :: odx1, odx2, odx3, odx4
      DOUBLE PRECISION :: hilf1(2,2,2,2), hilf2(2,2,2), hilf3(2,2), hilf4(2)

      tab%n1 = nr2
      tab%n2 = nlsigs
      tab%n3 = nncn
      tab%n4 = nwcb

      ALLOCATE( tab%x1(tab%n1) )
      ALLOCATE( tab%x2(tab%n2) )
      ALLOCATE( tab%x3(tab%n3) )
      ALLOCATE( tab%x4(tab%n4) )
      ALLOCATE( tab%ltable(tab%n1,tab%n2,tab%n3,tab%n4) )

      !===========================================================
      ! construct equidistant table:
      !===========================================================

      ! grid distances (also inverse):
      tab%dx1  = (otab%x1(otab%n1) - otab%x1(1)) / (tab%n1 - 1.0d0)  ! dr2
      tab%odx1 = 1.0d0 / tab%dx1
      tab%dx2  = (otab%x2(otab%n2) - otab%x2(1)) / (tab%n2 - 1.0d0)  ! dlsigs
      tab%odx2 = 1.0d0 / tab%dx2
      tab%dx3  = (otab%x3(otab%n3) - otab%x3(1)) / (tab%n3 - 1.0d0)  ! dncn
      tab%odx3 = 1.0d0 / tab%dx3
      tab%dx4  = (otab%x4(otab%n4) - otab%x4(1)) / (tab%n4 - 1.0d0)  ! dwcb
      tab%odx4 = 1.0d0 / tab%dx4

      ! grid vectors:
      DO i=1, tab%n1
        tab%x1(i) = otab%x1(1) + (i-1) * tab%dx1
      END DO
      DO i=1, tab%n2
        tab%x2(i) = otab%x2(1) + (i-1) * tab%dx2
      END DO
      DO i=1, tab%n3
        tab%x3(i) = otab%x3(1) + (i-1) * tab%dx3
      END DO
      DO i=1, tab%n4
        tab%x4(i) = otab%x4(1) + (i-1) * tab%dx4
      END DO

      
      ! Tetra-linear interpolation of the new equidistant lookuptable from
      ! the original non-equidistant table:

      ALLOCATE(iuv(tab%n1))
      ALLOCATE(juv(tab%n2))
      ALLOCATE(kuv(tab%n3))
      ALLOCATE(luv(tab%n4))


      DO l=1, tab%n1

        iuv(l) = 1
        DO ii=1, otab%n1 - 1
          IF (tab%x1(l) >= otab%x1(ii) .AND. tab%x1(l) <= otab%x1(ii+1)) THEN
            iuv(l) = ii
            EXIT
          END IF
        END DO

      END DO

      DO l=1, tab%n2

        juv(l) = 1
        DO ii=1, otab%n2 - 1
          IF (tab%x2(l) >= otab%x2(ii) .AND. tab%x2(l) <= otab%x2(ii+1)) THEN
            juv(l) = ii
            EXIT
          END IF
        END DO

      END DO

      DO l=1, tab%n3

        kuv(l) = 1
        DO ii=1, otab%n3 - 1
          IF (tab%x3(l) >= otab%x3(ii) .AND. tab%x3(l) <= otab%x3(ii+1)) THEN
            kuv(l) = ii
            EXIT
          END IF
        END DO

      END DO

      DO l=1, tab%n4

        luv(l) = 1
        DO ii=1, otab%n4 - 1
          IF (tab%x4(l) >= otab%x4(ii) .AND. tab%x4(l) <= otab%x4(ii+1)) THEN
            luv(l) = ii
            EXIT
          END IF
        END DO

      END DO

      ! Tetra-linear interpolation:
      DO i=1, tab%n1
        iu = iuv(i)
        odx1 = 1.0d0 / ( otab%x1(iu+1) - otab%x1(iu) )
        DO j=1, tab%n2
          ju = juv(j)
          odx2 = 1.0d0 / ( otab%x2(ju+1) - otab%x2(ju) )
          DO l=1, tab%n4
            lu = luv(l)
            odx4 = 1.0d0 / ( otab%x4(lu+1) - otab%x4(lu) )
!cdir nodep(hilf1)
            DO k=1, tab%n3
              ku = kuv(k)
              odx3 = 1.0d0 / ( otab%x3(ku+1) - otab%x3(ku) )
              hilf1 = otab%ltable( iu:iu+1, ju:ju+1, ku:ku+1, lu:lu+1)
              hilf2 = hilf1(1,1:2,1:2,1:2) + (hilf1(2,1:2,1:2,1:2) - hilf1(1,1:2,1:2,1:2)) * odx1 * ( tab%x1(i) - otab%x1(iu) )
              hilf3 = hilf2(1,1:2,1:2)     + (hilf2(2,1:2,1:2)     - hilf2(1,1:2,1:2)  )   * odx2 * ( tab%x2(j) - otab%x2(ju) )
              hilf4 = hilf3(1,1:2)         + (hilf3(2,1:2)         - hilf3(1,1:2)    )     * odx3 * ( tab%x3(k) - otab%x3(ku) )
              tab%ltable(i,j,k,l) = hilf4(1) +  ( hilf4(2) - hilf4(1) ) * odx4 * ( tab%x4(l) - otab%x4(lu) )
            END DO
          END DO
        END DO
      END DO

      ! clean up memory:
      DEALLOCATE(iuv,juv,kuv,luv)

      RETURN
    END SUBROUTINE equi_table

  END SUBROUTINE clnuc_sk_4d

END MODULE wolken


