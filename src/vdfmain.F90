!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFMAIN    ( CDCONF, &
 & KIDIA  , KFDIA  , KLON   , KLEV   , KLEVS  , KSTEP  , KTILES , &
 & KTRAC  , KLEVSN , KLEVI  , KDHVTLS, KDHFTLS, KDHVTSS, KDHFTSS, &
 & KDHVTTS, KDHFTTS, KDHVTIS, KDHFTIS, &
 & PTSPHY , KTVL   , KTVH   , KCNT   , PCVL   , PCVH   , PSIGFLT, &
 & PUM1   , PVM1   , PTM1   , PQM1   , PLM1   , PIM1   , PAM1   , PCM1   , &
 & PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , PTSKM1M, PTSAM1M, PWSAM1M, &
 & PSSRFL , PSLRFL , PEMIS  , PHRLW  , PHRSW  , &
 & PTSNOW , PTICE  , &
 & PSST   , PFRTI  , PALBTI , PWLMX  , &
 & PCHAR  , PUCURR , PVCURR , PTSKRAD, PCFLX  , &
 & PSOTEU , PSOTEV , PSOBETA, PVERVEL, &
 & PZ0M   , PZ0H   , &
 & PVDIS  , PVDISG , PAHFLEV, PAHFLSB, PFWSB  , PBIR   , PVAR   , &
 & PU10M  , PV10M  , PT2M   , PD2M   , PQ2M   , PZINV  , PBLH   , KHPBLN , KVARTOP, &
 & PSSRFLTI,PEVAPSNW,PGUST  , PZIDLWV, PWUAVG , LDNODECP,KPBLTYPE, PLDIFF, &
 & PFPLVL , PFPLVN , PFHPVL , PFHPVN , &
 & PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , &
 & PTE    , PQE    , PLE    , PIE    , PAE    , PVOM   , PVOL   , &
 & PTENC  , PTSKE1 , &
 & PUSTRTI, PVSTRTI, PAHFSTI, PEVAPTI, PTSKTI , &
 & PDIFTS , PDIFTQ , PDIFTL , PDIFTI , PSTRTU , PSTRTV , PTOFDU , PTOFDV, &
 & PSTRSOU, PSTRSOV,   PKH  , &
 & PDHTLS , PDHTSS , PDHTTS , PDHTIS)

!***

!**   *VDFMAIN* - DOES THE VERTICAL EXCHANGE OF U,V,SLG,QT BY TURBULENCE.

!     J.F.GELEYN       20/04/82   Original  
!     C.A.BLONDIN      18/12/86
!     A.C.M. BELJAARS  20/10/89   IFS-VERSION (TECHNICAL REVISION OF CY34)
!     A.C.M. BELJAARS  26/03/90   OBUKHOV-L UPDATE 
!     A.C.M. BELJAARS  30/09/98   SURFACE TILES 
!     P. Viterbo       17/05/2000 Surface DDH for TILES
!     D. Salmond       15/10/2001 FULLIMP mods
!     S. Abdalla       27/11/2001 Passing Zi/L to waves
!     A. Beljaars       2/05/2003 New tile coupling     
!     P.Viterbo        24/05/2004 Change surface units
!     M. Ko"hler        3/12/2004 Moist Advection-Diffusion
!     A. Beljaars       4/04/2005 Turb. orogr. drag
!     A. Beljaars      30/09/2005 Include Subgr. Oro. in solver
!     A.Beljaars       31/03/2005 Introduction of ocean current b.c.
!     P. Viterbo       17/06/2005 surf external library
!     M. Ko"hler        6/06/2006 Single Column Model option (LSCMEC) 
!     R. Neggers       24/05/2006 PBL dualM scheme + bimodal cloud scheme
 
!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE FOUR
!     PROGNOSTIC VARIABLES U,V,T AND Q DUE TO THE VERTICAL EXCHANGE BY
!     TURBULENT (= NON-MOIST CONVECTIVE) PROCESSES. THESE TENDENCIES ARE
!     OBTAINED AS THE DIFFERENCE BETWEEN THE RESULT OF AN IMPLICIT
!     TIME-STEP STARTING FROM VALUES AT T-1 AND THESE T-1 VALUES. ALL
!     THE DIAGNOSTIC COMPUTATIONS (EXCHANGE COEFFICIENTS, ...) ARE DONE
!      FROM THE T-1 VALUES. AS A BY-PRODUCT THE ROUGHNESS LENGTH OVER SEA
!     IS UPDATED ACCORDINGLY TO THE *CHARNOCK FORMULA. HEAT AND MOISTURE
!     SURFACE FLUXES AND THEIR DERIVATIVES AGAINST TS, WS AND WL
!     (THE LATTER WILL BE LATER WEIGHTED WITH THE SNOW FACTOR IN
!     *VDIFF*), LATER TO BE USED FOR SOIL PROCESSES TREATMENT, ARE ALSO
!     COMPUTED AS WELL AS A STABILITY VALUE TO BE USED AS A DIAGNOSTIC
!     OF THE DEPTH OF THE WELL MIXED LAYER IN CONVECTIVE COMPUTATIONS.

!     INTERFACE.
!     ----------
!          *VDIFF* TAKES THE MODEL VARIABLES AT T-1 AND RETURNS THE VALUES
!     FOR THE PROGNOSTIC TIME T+1 DUE TO VERTICAL DIFFUSION.
!     THE MODEL VARIABLES, THE MODEL DIMENSIONS AND THE DIAGNOSTICS DATA
!     ARE PASSED AS SUBROUTINE ARGUMENTS. CONSTANTS THAT DO NOT CHANGE
!     DURING A MODEL RUN (E.G. PHYSICAL CONSTANTS, SWITCHES ETC.) ARE
!     STORED IN A SINGLE COMMON BLOCK *YOMVDF*, WHICH IS INITIALIZED
!     BY SET-UP ROUTINE *SUVDF*.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*        NUMBER OF SOIL LAYERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                   OF SURFACE BOUNDARY CONDITION)
!    *KTRAC*        Number of tracers
!    *KLEVSN*       Number of snow layers (diagnostics) 
!    *KLEVI*        Number of sea ice layers (diagnostics)
!    *KDHVTLS*      Number of variables for individual tiles
!    *KDHFTLS*      Number of fluxes for individual tiles
!    *KDHVTSS*      Number of variables for snow energy budget
!    *KDHFTSS*      Number of fluxes for snow energy budget
!    *KDHVTTS*      Number of variables for soil energy budget
!    *KDHFTTS*      Number of fluxes for soil energy budget
!    *KDHVTIS*      Number of variables for sea ice energy budget
!    *KDHFTIS*      Number of fluxes for sea ice energy budget

!    *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION

!    *KCNT*         Index of vdf sub steps.

!     INPUT PARAMETERS (LOGICAL)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):

!    *PCVL*         LOW VEGETATION COVER                          -  
!    *PCVH*         HIGH VEGETATION COVER                         -  
!    *PSIGFLT*      STANDARD DEVIATION OF FILTERED OROGRAPHY      M
!    *PUM1*         X-VELOCITY COMPONENT                          M/S
!    *PVM1*         Y-VELOCITY COMPONENT                          M/S
!    *PTM1*         TEMPERATURE                                   K
!    *PQM1*         SPECIFIC HUMIDITY                             KG/KG
!    *PLM1*         SPECIFIC CLOUD LIQUID WATER                   KG/KG
!    *PIM1*         SPECIFIC CLOUD ICE                            KG/KG
!    *PAM1*         CLOUD FRACTION                                1
!    *PCM1*         TRACER CONCENTRATION                          KG/KG
!    *PAPM1*        PRESSURE ON FULL LEVELS                       PA
!    *PAPHM1*       PRESSURE ON HALF LEVELS                       PA
!    *PGEOM1*       GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL AT HALF LEVELS                   M2/S2
!    *PTSKM1M*      SKIN TEMPERATURE                              K
!    *PTSAM1M*      SURFACE TEMPERATURE                           K
!    *PWSAM1M*      SOIL MOISTURE ALL LAYERS                      M**3/M**3
!    *PSSRFL*       NET SHORTWAVE RADIATION FLUX AT SURFACE       W/M2
!    *PSLRFL*       NET LONGWAVE RADIATION FLUX AT SURFACE        W/M2
!    *PEMIS*        MODEL SURFACE LONGWAVE EMISSIVITY
!    *PHRLW*        LONGWAVE HEATING RATE                         K/s
!    *PHRSW*        SHORTWAVE HEATING RATE                        K/s
!    *PTSNOW*       SNOW TEMPERATURE                              K
!    *PTICE*        ICE TEMPERATURE (TOP SLAB)                    K
!    *PSST*         (OPEN) SEA SURFACE TEMPERATURE                K
!    *PFRTI*        TILE FRACTIONS                                (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PALBTI*       BROADBAND ALBEDO FOR TILE FRACTIONS
!    *PWLMX*        MAXIMUM SKIN RESERVOIR CAPACITY               kg/m**2
!    *PCHAR*        "EQUIVALENT" CHARNOCK PARAMETER
!    *PUCURR*       OCEAN CURRENT X_COMPONENT      
!    *PVCURR*       OCEAN CURRENT Y_COMPONENT      
!    *PTSKRAD*      SKIN TEMPERATURE OF LATEST FULL RADIATION
!                      TIMESTEP                                   K
!    *PCFLX*        TRACER SURFACE FLUX                           kg/(m2 s)
!    *PSOTEU*       Explicit part of U-tendency from subgrid orography scheme    
!    *PSOTEV*       Explicit part of V-tendency from subgrid orography scheme     
!    *PSOBETA*      Implicit part of subgrid orography 

!    *PVERVEL*      VERTICAL VELOCITY

!     INPUT PARAMETERS (LOGICAL):

!     CONTRIBUTIONS TO BUDGETS (OUTPUT,REAL):

!    *PVDIS*        TURBULENT DISSIPATION                         W/M2
!    *PVDISG*        SUBGRID OROGRAPHY DISSIPATION                 W/M2
!    *PAHFLEV*      LATENT HEAT FLUX  (SNOW/ICE FREE PART)        W/M2
!    *PAHFLSB*      LATENT HEAT FLUX  (SNOW/ICE COVERED PART)     W/M2

!     UPDATED PARAMETERS (REAL):

!    *PTE*          TEMPERATURE TENDENCY                          K/S
!    *PQE*          MOISTURE TENDENCY                             KG/(KG S)
!    *PLE*          LIQUID WATER TENDENCY                         KG/(KG S)
!    *PIE*          ICE WATER TENDENCY                            KG/(KG S)
!    *PAE*          CLOUD FRACTION TENDENCY                       1/S)
!    *PVOM*         MERIODINAL VELOCITY TENDENCY (DU/DT)          M/S2
!    *PVOL*         LATITUDE TENDENCY            (DV/DT)          M/S2
!    *PTENC*        TRACER TENDENCY                               KG/(KG S)
!    *PTSKE1*       SKIN TEMPERATURE TENDENCY                     K/S
!    *PZ0M*         AERODYNAMIC ROUGHNESS LENGTH                  M
!    *PZ0H*         ROUGHNESS LENGTH FOR HEAT                     M

!     UPDATED PARAMETERS FOR TILES (REAL): 

!    *PUSTRTI*      SURFACE U-STRESS                              N/M2 
!    *PVSTRTI*      SURFACE V-STRESS                              N/M2
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2
!    *PEVAPTI*      SURFACE MOISTURE FLUX                         KG/M2/S
!    *PTSKTI*       SKIN TEMPERATURE                              K

!     OUTPUT PARAMETERS (REAL):

!    *PFPLVL*       PBL PRECIPITATION FLUX AS RAIN                KG/(M**2*S)
!    *PFPLVN*       PBL PRECIPITATION FLUX AS SNOW                KG/(M**2*S)
!    *PFHPVL*       ENTHALPY FLUX OF PBL PRECIPITATION AS RAIN    J/(M**2*S)
!    *PFHPVN*       ENTHALPY FLUX OF PBL PRECIPITATION AS SNOW    J/(M**2*S)

!    *PLDIFF*       CONTRIB TO PBL CONDENSATE BY PASSIVE CLOUDS   KG/KG

!    *PFWSB*        EVAPORATION OF SNOW                           KG/(M**2*S)
!    *PU10M*        U-COMPONENT WIND AT 10 M                      M/S
!    *PV10M*        V-COMPONENT WIND AT 10 M                      M/S
!    *PT2M*         TEMPERATURE AT 2M                             K
!    *PD2M*         DEW POINT TEMPERATURE AT 2M                   K
!    *PQ2M*         SPECIFIC HUMIDITY AT 2M                       KG/KG
!    *PGUST*        GUST AT 10 M                                  M/S
!    *PZIDLWV*      Zi/L used for gustiness in wave model         M/M
!                   (NOTE: Positive values of Zi/L are set to ZERO)
!    *PBLH*         PBL HEIGHT (dry diagnostic based on Ri#)      M
!    *PZINV*        PBL HEIGHT (moist parcel, not for stable PBL) M
!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                      EACH TILE                                  W/M2
!    *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST            KG/(M2*S)
!    *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM            KG*(M/S)/(M2*S)
!    *PTOFDU*       TOFD COMP. OF TURBULENT FLUX OF U-MOMEMTUM   KG*(M/S)/(M2*S)
!    *PTOFDV*       TOFD COMP. OF TURBULENT FLUX OF V-MOMEMTUM   KG*(M/S)/(M2*S)
!    *PDIFTS*       TURBULENT FLUX OF HEAT                         J/(M2*S)
!    *PDIFTQ*       TURBULENT FLUX OF SPECIFIC HUMIDITY           KG/(M2*S)
!    *PDIFTL*       TURBULENT FLUX OF LIQUID WATER                KG/(M2*S)
!    *PDIFTI*       TURBULENT FLUX OF ICE WATER                   KG/(M2*S)
!    *PSTRSOU*      SUBGRID OROGRAPHY FLUX OF U-MOMEMTUM    KG*(M/S)/(M2*S)
!    *PSTRSOV*      SUBGRID OROGRAPHY FLUX OF V-MOMEMTUM    KG*(M/S)/(M2*S)

!    *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!                   IN SURFACE LAYER: CH*U                        (M/S)
!    *PDHTLS*       Diagnostic array for tiles (see module yomcdh)
!                      (Wm-2 for energy fluxes, kg/(m2s) for water fluxes)
!    *PDHTSS*       Diagnostic array for snow T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTTS*       Diagnostic array for soil T (see module yomcdh)
!                      (Wm-2 for fluxes)
!    *PDHTIS*       Diagnostic array for ice T (see module yomcdh)
!                      (Wm-2 for fluxes)

!     Additional parameters for flux boundary condition (in SCM model):

!    *LSFCFLX*      If .TRUE. flux boundary condtion is used 
!    *REXTSHF*      Specified sensible heat flux [W/m2]
!    *REXTLHF*      Specified latent heat flux [W/m2]

!     METHOD.
!     -------

!          FIRST AN AUXIALIARY VARIABLE CP(Q)T+GZ IS CREATED ON WHICH
!     THE VERTICAL DIFFUSION PROCESS WILL WORK LIKE ON U,V AND Q. THEN
!     ALONG THE VERTICAL AND AT THE SURFACE, EXCHANGE COEFFICIENTS (WITH
!     THE DIMENSION OF A PRESSURE THICKNESS) ARE COMPUTED FOR MOMENTUM
!     AND FOR HEAT (SENSIBLE PLUS LATENT). THE LETTERS M AND H ARE USED
!     TO DISTINGUISH THEM AND THE COMPUTATION IS THE RESULT OF A
!     CONDITIONAL MERGE BETWEEN THE STABLE AND THE UNSTABLE CASE
!     (DEPENDING ON THE SIGN OF THE *RICHARDSON BULK NUMBER).
!          IN THE SECOND PART OF THE ROUTINE THE IMPLICIT LINEAR
!     SYSTEMS FOR U,V FIRST AND T,Q SECOND ARE SOLVED BY A *GAUSSIAN
!     ELIMINATION BACK-SUBSTITUTION METHOD. FOR T AND Q THE LOWER
!     BOUNDARY CONDITION DEPENDS ON THE SURFACE STATE.
!     OVER LAND, TWO DIFFERENT REGIMES OF EVAPORATION PREVAIL:
!     A STOMATAL RESISTANCE DEPENDENT ONE OVER THE VEGETATED PART
!     AND A SOIL RELATIVE HUMIDITY DEPENDENT ONE OVER THE
!     BARE SOIL PART OF THE GRID MESH.
!     POTENTIAL EVAPORATION TAKES PLACE OVER THE SEA, THE SNOW
!     COVERED PART AND THE LIQUID WATER COVERED PART OF THE
!     GRID MESH AS WELL AS IN CASE OF DEW DEPOSITION.
!          FINALLY ONE RETURNS TO THE VARIABLE TEMPERATURE TO COMPUTE
!     ITS TENDENCY AND THE LATER IS MODIFIED BY THE DISSIPATION'S EFFECT
!     (ONE ASSUMES NO STORAGE IN THE TURBULENT KINETIC ENERGY RANGE) AND
!     THE EFFECT OF MOISTURE DIFFUSION ON CP. Z0 IS UPDATED AND THE
!     SURFACE FLUXES OF T AND Q AND THEIR DERIVATIVES ARE PREPARED AND
!     STORED LIKE THE DIFFERENCE BETWEEN THE IMPLICITELY OBTAINED
!     CP(Q)T+GZ AND CP(Q)T AT THE SURFACE.

!     EXTERNALS.
!     ----------

!     *VDFMAIN* CALLS SUCESSIVELY:
!         *SURFEXCDRIVER*
!         *VDFEXCU*
!         *VDFTOFDC*
!         *VDFDIFM*
!         *VDFDIFH*
!         *VDFDIFC*
!         *VDFINCR*
!         *VDFSDRV*
!         *VDFPPCFL*
!         *VDFUPDZ0*

!     REFERENCE.
!     ----------

!          SEE VERTICAL DIFFUSION'S PART OF THE MODEL'S DOCUMENTATION
!     FOR DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     ------------------------------------------------------------------
use garbage, only : foealfa, foeewm, foeldcpm, surfpp, surfexcdriver, vdfdpbl, vdftofdc, vdffblend
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

USE YOMCT0   , ONLY : LSCMEC   ,LSFCFLX  ,REXTSHF  ,REXTLHF
USE YOEVDF   , ONLY : RVDIFTS  ,LLDIAG
USE YOMCST   , ONLY : RG       ,RD       ,&
                    & RCPD     ,RETV     ,RLVTT    ,RLSTT    ,RTT 
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
                    & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
                    & R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
                    & RTICECU  ,RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOMJFH   , ONLY : N_VMASS
USE YOEPHY   , ONLY : LVDFTRAC 

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTLS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTSS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTTS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHVTIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDHFTIS 
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCNT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOTEV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOBETA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDISG(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAHFLSB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFWSB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBIR(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVAR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBLH(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KHPBLN(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSRFLTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZIDLWV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUAVG(KLON)
LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLDIFF(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KVARTOP(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PIE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQ(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOFDU(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOFDV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRSOV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTLS(KLON,KTILES,KDHVTLS+KDHFTLS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTSS(KLON,KLEVSN,KDHVTSS+KDHFTSS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTTS(KLON,KLEVS,KDHVTTS+KDHFTTS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHTIS(KLON,KLEVI,KDHVTIS+KDHFTIS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)     :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZDIFTQT(KLON,0:KLEV), ZDIFTSLG(KLON,0:KLEV) 

REAL(KIND=JPRB) ::    ZCPTGZ(KLON,KLEV) , ZCFM(KLON,KLEV)   , ZCFH(KLON,KLEV)   ,&
                    & ZUDIF(KLON,KLEV)  , ZVDIF(KLON,KLEV)  ,&
                    & ZQTDIF(KLON,KLEV) , ZSLGDIF(KLON,KLEV),&
                    & ZSLGM1(KLON,KLEV) , ZQTM1(KLON,KLEV)  , ZQTE(KLON,KLEV)   ,&
                    & ZSLGE(KLON,KLEV)  , ZTOFDC(KLON,KLEV) , ZSOC(KLON,KLEV)
REAL(KIND=JPRB) ::    ZKHFL(KLON)       , ZKQFL(KLON)       , ZKMFL(KLON)  
REAL(KIND=JPRB) ::    ZQEA(KLON,KLEV)   , ZLEA(KLON,KLEV)   , ZIEA(KLON,KLEV)   ,&
                    & ZQTEA(KLON,KLEV)  , ZSLGEA(KLON,KLEV) , ZAEA(KLON,KLEV)   ,&
                    & ZTEA(KLON,KLEV)   , ZUEA(KLON,KLEV)   , ZVEA(KLON,KLEV)   ,&
                    & ZSLGEWODIS(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZQTEP(KLON,KLEV)  , ZSLGEP(KLON,KLEV), ZDZRHOI

!RN --- variables associated with dualM scheme -------------------------

INTEGER(KIND=JPIM), PARAMETER   ::   IDRAFT = 3    ! nr of updrafts

   !  thermodynamic transport
REAL(KIND=JPRB) ::  ZMFLX(KLON,0:KLEV,IDRAFT)  ,  ZQTUH(KLON,0:KLEV,IDRAFT),&   
                  & ZSLGUH(KLON,0:KLEV,IDRAFT) ,  ZWUH(KLON,0:KLEV,IDRAFT)
   !  momentum transport
REAL(KIND=JPRB) ::    ZMFLXM(KLON,0:KLEV,IDRAFT)            , ZUUH(KLON,0:KLEV,IDRAFT),&
                    & ZVUH(KLON,0:KLEV,IDRAFT) ,&
                    & ZUCURR(KLON)      , ZVCURR(KLON)      , ZTAUX(KLON)       ,&
                    & ZTAUY(KLON)

REAL(KIND=JPRB) ::   ZFRACB(KLON,IDRAFT) , ZZPTOP(KLON,IDRAFT) , ZZPLCL(KLON,IDRAFT)
INTEGER(KIND=JPIM) :: IPTOP(KLON,IDRAFT) , IPLCL(KLON,IDRAFT)  , IPLZB(KLON,IDRAFT)  

!RN --------------------------------------------------------------------

REAL(KIND=JPRB) ::    ZZ0MW(KLON)       , ZZ0HW(KLON)       , ZZ0QW(KLON)       ,&
                    & ZBLEND(KLON)      , ZFBLEND(KLON)
REAL(KIND=JPRB) ::    ZZCPTS(KLON)      , ZZQSA(KLON)       , ZZBUOM(KLON)      ,&
                    & ZZZDL(KLON)
REAL(KIND=JPRB) ::    ZTUPD(KLON,KLEV)  , ZQUPD(KLON,KLEV)  , ZLUPD(KLON,KLEV)  ,&
                    & ZIUPD(KLON,KLEV)  , ZQTUPD(KLON,KLEV) , ZLIUPD(KLON,KLEV) ,&
                    & ZSLGUPD(KLON,KLEV), ZAUPD(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZTAUXCG(KLON,KLEV), ZTAUYCG(KLON,KLEV)
REAL(KIND=JPRB) ::    ZQSVAR(KLON,KLEV) , ZANEW(KLON,KLEV)  , ZLNEW(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZSVFLUXCLD(KLON,0:KLEV)               , ZSVFLUXSUB(KLON,0:KLEV),&
                    & ZSVFLUX(KLON,0:KLEV),ZBUOYPOS(KLON)   , ZBUOYNEG(KLON)    ,&
                    & ZDZH(KLON,0:KLEV)
REAL(KIND=JPRB) ::    ZALFA1            , ZALFA2            , ZDELQ             ,&
                    & ZCORQS(KLON,KLEV) , ZDQSDTEMP(KLON,KLEV)

REAL(KIND=JPRB) ::    ZCLDBASE(KLON)    , ZCLDTOP(KLON)     , &
                    & ZRICUI(KLON)      , ZMCU(KLON)        , ZDTHV(KLON)

REAL(KIND=JPRB) ::    ZPFLXUSUM(KLON,0:KLEV)

!RN --- VDF qt variance budget & bimodal cloud scheme variables ---------

REAL(KIND=JPRB) ::    ZVARGEN           , ZVARTRANS         , ZVARDISS          ,&
                    & ZTAU(KLON,KLEV)   , ZTAUNEW           , &
                    & ZDQTDZ            , ZWQTF             , &
                    & ZWQT2(KLON,0:KLEV), & 
                    & ZQTUP(KLON,KLEV)  , ZSLGUP(KLON,KLEV) , &               
                    & ZQTTEST(KLON,KLEV)   
                    
REAL(KIND=JPRB) ::    ZCLDFRAC(KLON,KLEV),ZQLAV(KLON,KLEV)

!RN ---------------------------------------------------------------------

REAL(KIND=JPRB) ::    ZTSKINTIOLD(KLON,KTILES)

REAL(KIND=JPRB) ::    ZEXTSHF(KLON)     , ZEXTLHF(KLON)

REAL(KIND=JPRB) ::    ZCPTSTI(KLON,KTILES), ZQSTI(KLON,KTILES)  ,&
                    & ZDQSTI(KLON,KTILES) , ZCSATTI(KLON,KTILES),&
                    & ZCAIRTI(KLON,KTILES), ZCFHTI(KLON,KTILES) ,&
                    & ZCFQTI(KLON,KTILES) , ZAHFLTI(KLON,KTILES),&
                    & ZTSKTIP1(KLON,KTILES)
REAL(KIND=JPRB) ::    ZSTR(KLON,KTILES)   , ZG0(KLON,KTILES)

LOGICAL ::            LLRUNDRY(KLON), LLPBL(KLON,KLEV), &
                      LLTROPDAMP, LLPERTURB

INTEGER(KIND=JPIM) :: ITOP, JK, JL, JD, KCAP,JKK, KPERTURB

REAL(KIND=JPRB) ::    ZGDPH, ZRHO, ZTMST, ZRG, ZRTMST
LOGICAL ::            LLSFCFLX
REAL(KIND=JPRB) ::    ZHU1
REAL(KIND=JPRB) ::    ZALFAW(KLON,KLEV), ZFACW, ZFACI, ZFAC, ZESDP, ZCOR

REAL(KIND=JPRB) ::    ZDUMMY2, ZDUMMY3, ZCFNC1, ZMGEOM

REAL(KIND=JPRB) ::    ZDETR(KLON,KLEV), ZAWAKE(KLON,KLEV)

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

INTERFACE
#include "surfexcdriver.h"
#include "surfpp.h"
END INTERFACE

#include "vdfdifh.intfb.h"
#include "vdfdifh5.intfb.h"
#include "vdfdifm.intfb.h"
#include "vdfdifm2.intfb.h"
#include "vdfdifc.intfb.h"
#include "vdfdpbl.intfb.h"
#include "vdfexcu.intfb.h"
#include "vdfhghtn.intfb.h"
#include "vdfincr.intfb.h"
#include "vdffblend.intfb.h"
#include "vdfcloud.intfb.h"
#include "vdftofdc.intfb.h"

#include "fcttre.h"



!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFMAIN',0,ZHOOK_HANDLE)

ZTMST       = PTSPHY
ZRTMST      = 1.0_JPRB/PTSPHY    ! optimization
ZRG         = 1.0_JPRB/RG        !     -"-
LLRUNDRY(:) = .FALSE.  ! option to run dry updrafts with no condensation
!DO JL=KIDIA,KFDIA
! IF ( .NOT. LDNODECP(JL) )  LLRUNDRY(JL) = .TRUE. ! run dry for real vdfmain's
! IF ( .NOT. LDNODECP(JL) )  PBIR(JL) = 1.0        ! always decouple for real vdfmain's
!ENDDO

!-- Switch for eliminating LS tendencies in troposphere --
!LLTROPDAMP = .TRUE.
LLTROPDAMP = .FALSE.

!-- Switch for ML humidity perturbation --
!LLPERTURB = .TRUE.
LLPERTURB = .FALSE.
KPERTURB = 2*96


!*         1.0  SCM: Fixed fluxes for flux boundary condition ([W/m^2] downward)

IF (LSCMEC) THEN
  LLSFCFLX   = LSFCFLX   ! scm namelist parameters
  ZEXTSHF(:) = REXTSHF
  ZEXTLHF(:) = REXTLHF
ELSE
  LLSFCFLX   = .FALSE.
  ZEXTSHF(:) = 0.0_JPRB
  ZEXTLHF(:) = 0.0_JPRB
ENDIF


!*         1.1  Store initial tendencies for flux calculation
!*              and initialize variable.

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQEA(JL,JK)=PQE(JL,JK)
    ZLEA(JL,JK)=PLE(JL,JK)
    ZIEA(JL,JK)=PIE(JL,JK)
    ZAEA(JL,JK)=PAE(JL,JK)
    ZTEA(JL,JK)=PTE(JL,JK)
    ZUEA(JL,JK)=PVOM(JL,JK)
    ZVEA(JL,JK)=PVOL(JL,JK)
    ZQTEP(JL,JK)  = 0._JPRB
    ZSLGEP(JL,JK) = 0._JPRB
    PLDIFF(JL,JK) = 0._JPRB
  ENDDO
ENDDO

DO JD=1,IDRAFT
  DO JK=0,KLEV
    DO JL=KIDIA,KFDIA
      ZMFLX(JL,JK,JD)  = 0.0_JPRB
      ZSLGUH(JL,JK,JD) = 0.0_JPRB
      ZQTUH(JL,JK,JD)  = 0.0_JPRB
      ZWUH(JL,JK,JD)   = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    PFPLVL(JL,JK) = 0.0_JPRB
    PFPLVN(JL,JK) = 0.0_JPRB
    PFHPVL(JL,JK) = 0.0_JPRB
    PFHPVN(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO



!     ------------------------------------------------------------------

!*         2.     NEW VARIABLES S, SLG, QT
!*                (at initial time level)
!                 ------------------------

!RN ####### !temporary for testing purposes! #######
!RN - disable any tendency above PBL
!IF (LLTROPDAMP) THEN
!
!  KCAP = 45  !3000m in L60 
!  IF (KLEV>70) THEN 
!    KCAP = 71  !3000m in L91
!  ENDIF
!  
!  DO JK=1,KCAP
!    DO JL=KIDIA,KFDIA
!      PQE(JL,JK) = 0._JPRB
!      PTE(JL,JK) = 0._JPRB
!    ENDDO
!  ENDDO
!  
!ENDIF  


DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

!*         2.1  dry static energy cp(q)*T + gz

    ZCPTGZ(JL,JK)  =PGEOM1(JL,JK) + PTM1(JL,JK)*RCPD * (1.0_JPRB + RVTMP2*PQM1(JL,JK))

!*         2.2  total water and generalized liquid water static energy 
!*              slg = cp*T + gz - Lcond*ql - Ldep*qi

    ZSLGM1(JL,JK) = ZCPTGZ(JL,JK) - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)
    ZSLGE(JL,JK)  = RCPD * ( ( 1.0_JPRB + RVTMP2 * PQM1(JL,JK) ) * PTE(JL,JK)   &  !dcpT/dt
                &                       + RVTMP2 * PTM1(JL,JK)   * PQE(JL,JK) ) &  
                & - RLVTT * PLE(JL,JK) - RLSTT * PIE(JL,JK)                        !dLqli/dt  

    ZQTM1(JL,JK ) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)
    ZQTE(JL,JK)   = PQE(JL,JK)  + PLE(JL,JK)  + PIE(JL,JK)             !rad+dyn. qt tendency
    
    ZSLGEA(JL,JK) = ZSLGE(JL,JK)
    ZQTEA(JL,JK)  = ZQTE(JL,JK)
    
  ENDDO
ENDDO



!     ------------------------------------------------------------------

!*         3.  Compute all surface related quantities
!          ------------------------------------------

DO JL=KIDIA,KFDIA
  DO JK = 1,KTILES
    ZTSKINTIOLD(JL,JK) = PTSKTI(JL,JK)
  ENDDO   
ENDDO   

CALL SURFEXCDRIVER(CDCONF=CDCONF, &
 & KIDIA=KIDIA, KFDIA=KFDIA, KLON=KLON, KLEVS=KLEVS, KTILES=KTILES, KSTEP=KSTEP, &
 & KLEVSN=KLEVSN, KLEVI=KLEVI, KDHVTLS=KDHVTLS, KDHFTLS=KDHFTLS, &
 & KDHVTSS=KDHVTSS, KDHFTSS=KDHFTSS, KDHVTTS=KDHVTTS, KDHFTTS=KDHFTTS, &
 & KDHVTIS=KDHVTIS, KDHFTIS=KDHFTIS, K_VMASS=N_VMASS, &
 & PTSTEP=PTSPHY, PRVDIFTS=RVDIFTS, &
! input data, non-tiled
 & KTVL=KTVL, KTVH=KTVH, PCVL=PCVL, PCVH=PCVH, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PTMLEV=PTM1(:,KLEV), &
 & PQMLEV=PQM1(:,KLEV), PAPHMS=PAPHM1(:,KLEV), PGEOMLEV=PGEOM1(:,KLEV), &
 & PCPTGZLEV=ZCPTGZ(:,KLEV), PSST=PSST, PTSKM1M=PTSKM1M, PCHAR=PCHAR, &
 & PSSRFL=PSSRFL, PSLRFL=PSLRFL, PEMIS=PEMIS, PTICE=PTICE, PTSNOW=PTSNOW, &
 & PWLMX=PWLMX, PUCURR=PUCURR, PVCURR=PVCURR, &
! input data, soil
 & PTSAM1M=PTSAM1M, PWSAM1M=PWSAM1M, &
! input data, tiled
 & PFRTI=PFRTI, PALBTI=PALBTI, &
! updated data, tiled
 & PUSTRTI=PUSTRTI, PVSTRTI=PVSTRTI, PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, &
 & PTSKTI=PTSKTI, &
! updated data, non-tiled
 & PZ0M=PZ0M, PZ0H=PZ0H, &
! output data, tiled
 & PSSRFLTI=PSSRFLTI, PQSTI=ZQSTI, PDQSTI=ZDQSTI, PCPTSTI=ZCPTSTI, &
 & PCFHTI=ZCFHTI, PCFQTI=ZCFQTI, PCSATTI=ZCSATTI, PCAIRTI=ZCAIRTI, &
! output data, non-tiled
 & PKHLEV=PKH(:,KLEV), PCFMLEV=ZCFM(:,KLEV), PKMFL=ZKMFL, PKHFL=ZKHFL, &
 & PKQFL=ZKQFL, PEVAPSNW=PEVAPSNW, PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
 & PBLENDPP=ZBLEND, PCPTSPP=ZZCPTS, PQSAPP=ZZQSA, PBUOMPP=ZZBUOM, &
 & PZDLPP=ZZZDL, &
! output data, diagnostics
 & PDHTLS=PDHTLS, PDHTSS=PDHTSS, PDHTTS=PDHTTS, PDHTIS=PDHTIS &
 & )


                       
!RN --- surface heat flux check ---
DO JL=KIDIA,KFDIA
!  If (.false.) THEN
  IF (ZKHFL(JL).lt.-10._JPRB.or.ZKQFL(JL).lt.-10._JPRB) then
        
!    write(0,'(a,i,f10.6)') 'vdfmain:',JL,PTSKM1M(JL)
    write(0,'(a,2f10.6)')   '  vdfmain: sfluxes:',ZKHFL(JL),ZKQFL(JL)
    
!    DO JK = 1,KTILES
!     write(0,'(a,i,2f10.6)')'  vdfmain: Tskin:',JK,PFRTI(JL,JK),PTSKTI(JL,JK)
!    ENDDO   
    
  ENDIF
!  ENDIF
ENDDO


!RN --- surface heat flux limiters ---
!DO JL=KIDIA,KFDIA
!  ZKHFL(JL) = MAX( ZKHFL(JL), -1000.0_JPRB / RCPD)
!  ZKQFL(JL) = MAX( ZKQFL(JL), -2000.0_JPRB / RLVTT )
!ENDDO



!     ------------------------------------------------------------------

!*         4.     EXCHANGE COEFFICIENTS
!                 ---------------------

!*         4.4  COMPUTATION OF THE PBL EXTENSION


!          SET PBL HEIGHT-INDEX TO 1

DO JL=KIDIA,KFDIA
  KHPBLN(JL)=1
  KVARTOP(JL)  = 0
ENDDO
ITOP=1  !ITOP is used in some solvers. ITOP=1 means: always integrate over whole atmosphere depth


!          FLUX BOUNDARY CONDITION

IF (LLSFCFLX) THEN
  DO JL=KIDIA,KFDIA
    !fixed prescribed surface fluxes
    ZRHO = PAPHM1(JL,KLEV)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    ZKHFL(JL) = ZEXTSHF(JL) / ( RCPD*(1.0_JPRB+RVTMP2*PQM1(JL,KLEV)) ) / ZRHO
    ZKQFL(JL) = ZEXTLHF(JL) / RLVTT / ZRHO
  ENDDO
ENDIF


!*         4.5  BOUNDARY LAYER HEIGHT FOR DIANOSTICS ONLY

CALL VDFDPBL(KIDIA,KFDIA,KLON,KLEV,&
 & PUM1,PVM1,PTM1,PQM1,PGEOM1,&
 & ZKMFL,ZKHFL,ZKQFL,PBLH)  


!*         4.6  ORGANIZED UPDRAFTS

CALL VDFHGHTN (KIDIA   , KFDIA   , KLON    , KLEV    , IDRAFT   , ZTMST , KSTEP, &
             & PUM1    , PVM1    , PTM1    , PQM1    , PLM1    , PIM1    , PAM1     , &
             & PAPHM1  , PAPM1   , PGEOM1  , PGEOH   , PVERVEL , PQE     , PTE      , &
             & ZKMFL   , ZKHFL   , ZKQFL   , ZMFLX   , &
! DIAGNOSTIC OUTPUT
             & PEXTR2  , KFLDX2  , PEXTRA  , KLEVX   , KFLDX    , &
!              
             & ZUUH    , ZVUH    , ZSLGUH  , ZQTUH   , ZFRACB   , ZWUH  , &
             & ZZPTOP  , IPTOP   , ZZPLCL  , IPLCL   , IPLZB    , &
             & PWUAVG  , ZRICUI  , ZMCU    , ZDTHV   , &
             & PFPLVL  , PFPLVN  , ZDETR   , &
             & PBIR    , LDNODECP, LLRUNDRY, KPBLTYPE, ZWQT2 )


!--- set some PBL heights based on i) PBL type and ii) various updraft properties ---
DO JL=KIDIA,KFDIA
    
  SELECT CASE (KPBLTYPE(JL))
    
      CASE(0)
          !Stable PBL
          PZINV(JL)    = 0.0_JPRB
          KHPBLN(JL)   = KLEV
          ZCLDBASE(JL) = -100._JPRB
          ZCLDTOP(JL)  = -100._JPRB
          KVARTOP(JL)  = 0
     
      CASE(1)
          !Dry convective PBL
          PZINV(JL)    = ZZPTOP(JL,2)
          KHPBLN(JL)   = IPTOP(JL,2)
          ZCLDBASE(JL) = -100._JPRB
          ZCLDTOP(JL)  = -100._JPRB
          KVARTOP(JL)  = IPTOP(JL,2)

      CASE(2)
          !Stratocumulus
          PZINV(JL)    = ZZPTOP(JL,3)
          KHPBLN(JL)   = IPTOP(JL,3)
          ZCLDBASE(JL) = ZZPLCL(JL,3)
          ZCLDTOP(JL)  = ZZPTOP(JL,3)
          KVARTOP(JL)  = IPTOP(JL,3)
	  
      CASE(3)
          !Shallow cumulus
          PZINV(JL)    = ZZPLCL(JL,3)
          KHPBLN(JL)   = IPLCL(JL,3)+1   !equivalent to definition of top level: use first half level *BELOW* boundary
          ZCLDBASE(JL) = ZZPLCL(JL,3)
          ZCLDTOP(JL)  = ZZPTOP(JL,3)
          KVARTOP(JL)  = IPTOP(JL,1)   !PBL includes the inversion between moist updraft top and test updraft top
		
      CASE(4)
          !Deep cumulus - only do a dry subcloud ML
          PZINV(JL)    = ZZPTOP(JL,2)
          KHPBLN(JL)   = IPTOP(JL,2)
          ZCLDBASE(JL) = ZZPLCL(JL,1)
          ZCLDTOP(JL)  = ZZPTOP(JL,1)
          KVARTOP(JL)  = IPTOP(JL,2)
    
  END SELECT 
  
ENDDO !JL


!--- set PBL indicator ---   
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    LLPBL(JL,JK) = JK >= KVARTOP(JL) .AND. KVARTOP(JL)>0
  ENDDO
ENDDO


!*         4.7  EXCHANGE COEFFICIENTS ABOVE THE SURFACE LAYER

CALL VDFEXCU(KIDIA  , KFDIA  , KLON   , KLEV    , IDRAFT  , ZTMST  , PZ0M   , &
           & PHRLW  , PHRSW  , &
           & PUM1   , PVM1   , PTM1   , PQM1    , PLM1    , PIM1   , &
           & PAPHM1 , PAPM1  , PGEOM1 , PGEOH   , ZCPTGZ  , &
! DIAGNOSTIC OUTPUT
            &PEXTR2 , KFLDX2 , PEXTRA , KLEVX   , KFLDX   , &
!              
           & ZKMFL  , ZKHFL  , ZKQFL  , ZCFM    , ZCFH    , ZTAUXCG, ZTAUYCG, &
           & ZRICUI , ZMCU   , ZDTHV  , ZMFLX   , KVARTOP , &
           & PZINV  , KHPBLN , PKH    , ZCLDBASE, ZCLDTOP , KPBLTYPE)  


!*         4.8  MASS FLUX MODIFICATIONS (FOR SOLVER STABILITY)

DO JD=2,IDRAFT

  !-- Remove single and/or double massflux layers --
  DO JL=KIDIA,KFDIA
    IF ( ZMFLX(JL,KLEV-2,JD) < 1.E-40_JPRB ) THEN
    !IF ( ZMFLX(JL,KLEV-3,JD) < 1.E-40_JPRB ) THEN
      ZMFLX(JL,KLEV-1,JD) = 0._JPRB
      !ZMFLX(JL,KLEV-2,JD) = 0._JPRB
    ENDIF
  ENDDO
  
  !-- Prune massflux-spikes in top PBL layer --
  DO JK=2,KLEV-2
  DO JL=KIDIA,KFDIA
      IF ( ZMFLX(JL,JK-1,JD) < 1.E-40_JPRB .AND. ZMFLX(JL,JK,JD) > 1.E-40_JPRB ) THEN
        ZMFLX(JL,JK,JD) = MIN( ZMFLX(JL,JK,JD) , 1.5_JPRB*ZMFLX(JL,JK+1,JD) )
      ENDIF
  ENDDO
  ENDDO
  
ENDDO


!*         4.9     TURBULENT OROGRAPHIC DRAG COEFFICIENTS 

CALL VDFTOFDC(KIDIA,KFDIA,KLON,KLEV,ZTMST,&
 & PUM1,PVM1,PGEOM1,PSIGFLT,&
 & ZTOFDC)  



!     ------------------------------------------------------------------

!*         5.     SOLVE ADVECTION-DIFFUSION EQUATION
!                 ----------------------------------

!*         5.1  MOMENTUM

ZMFLXM(:,:,:)=ZMFLX(:,:,:)
ZUCURR(:)    =0._JPRB   ! ocean currents not yet active
ZVCURR(:)    =0._JPRB   ! ...
ZHU1=RVDIFTS*ZTMST
ZSOC(KIDIA:KFDIA,1:KLEV)=PSOBETA(KIDIA:KFDIA,1:KLEV)*ZHU1

!...counter-gradient momentum transport (Brown et al 2007)
!CALL VDFDIFM2 (KIDIA, KFDIA, KLON , KLEV  , ITOP, &
!            & ZTMST, PUM1 , PVM1 , PAPHM1, ZCFM, ZTOFDC, &
!            & ZTAUXCG, ZTAUYCG, PSOTEU,PSOTEV,ZSOC ,&
!            & PVOM,PVOL,PUCURR,PVCURR,ZUDIF,ZVDIF)

CALL VDFDIFM (KIDIA, KFDIA, KLON , KLEV  , IDRAFT , ITOP  , &
            & ZTMST, PUM1 , PVM1 , PAPHM1, ZCFM   , ZMFLXM, ZUUH , ZVUH ,&
            & ZTOFDC,PSOTEU,PSOTEV,ZSOC  , &
            & PVOM , PVOL , ZUCURR,ZVCURR, ZUDIF  , ZVDIF , ZTAUX, ZTAUY)  


!*         5.2  GENERALIZED LIQUID WATER STATIC ENERGY AND TOTAL WATER

CALL VDFDIFH (KIDIA  , KFDIA  , KLON   , KLEV   , IDRAFT , ITOP   , KTILES, &
            & ZTMST  , ZEXTSHF, ZEXTLHF, LLSFCFLX, &
            & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
            & ZSLGM1 , PTM1   , PQM1   , ZQTM1  , PAPHM1 , &
            & ZCFH   , ZCFHTI , ZCFQTI , ZMFLX  , ZSLGUH , ZQTUH  , &
            & ZSLGDIF, ZQTDIF , ZCPTSTI, ZQSTI  , ZCAIRTI, ZCSATTI, &
            & ZDQSTI , PTSKTI , PTSKRAD, PTSAM1M(1,1)    , PTSNOW , PTICE  , PSST, &
            & ZTSKTIP1,ZSLGE  , PTE    , ZQTE, &
            & PEVAPTI, PAHFSTI, ZAHFLTI, ZSTR   , ZG0)

!...fully upstream M (phi_up - phi_bar) term
!CALL VDFDIFH5(KIDIA  , KFDIA  , KLON   , KLEV   , IDRAFT , ITOP   , KTILES, &
!            & ZTMST  , ZEXTSHF, ZEXTLHF, LLSFCFLX, &
!            & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
!            & ZSLGM1 , PTM1   , PQM1   , ZQTM1  , PAPHM1 , &
!            & ZCFH   , ZCFHTI , ZCFQTI , ZMFLX  , ZSLGUH , ZQTUH  , &
!            & ZSLGDIF, ZQTDIF , ZCPTSTI, ZQSTI  , ZCAIRTI, ZCSATTI, &
!            & ZDQSTI , PTSKTI , PTSKRAD, PTSAM1M(1,1)    , PTSNOW , PTICE  , PSST, &
!            & ZTSKTIP1,ZSLGE  , PTE    , ZQTE, &
!            & PEVAPTI, PAHFSTI, ZAHFLTI, ZSTR   , ZG0)


!*         5.3  INCREMENTATION OF U AND V TENDENCIES, STORAGE OF
!*              THE DISSIPATION, COMPUTATION OF MULTILEVEL FLUXES.

CALL VDFINCR (KIDIA  , KFDIA  , KLON   , KLEV   , ITOP   , ZTMST  , &
            & PUM1   , PVM1   , ZSLGM1 , PTM1   , ZQTM1  , PAPHM1 , PGEOM1 , &
            & ZCFM   , ZTOFDC , PSOTEU , PSOTEV , ZSOC   ,&
            & ZUDIF  , ZVDIF  , PUCURR , PVCURR ,ZSLGDIF, ZQTDIF , &
            & PVOM   , PVOL   , ZSLGE  , ZQTE   , ZSLGEWODIS, &
            & PVDIS  , PVDISG , PSTRTU , PSTRTV , PSTRSOU, PSTRSOV , PTOFDU , PTOFDV)  


!          5.4  Solve for tracers

IF (LVDFTRAC .AND. KTRAC > 0) THEN 
  CALL VDFDIFC(KIDIA,KFDIA,KLON,KLEV,ITOP,KTRAC,&
             & ZTMST,PCM1,PTENC,PAPHM1,ZCFH,PCFLX)
ENDIF



!     ------------------------------------------------------------------
!
!*        6.     TIME INTEGRATION OF CONSERVED STATE VARIABLES QT AND SLG
!                 --------------------------------------------------------
!
!         Compute the conserved state after rad+dyn *AND* pbl conv+diff.
!         This will be used later to obtain the tendencies of the non-conserved 
!         prognostic variables T and QV.
!
  
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      
      !calculate updraft precipitation flux divergence (tendencies)
      ZDZRHOI = RG/( PAPHM1(JL,JK)-PAPHM1(JL,JK-1) )
      ZQTEP(JL,JK)  = -( PFPLVL(JL,JK) - PFPLVL(JL,JK-1) ) * ZDZRHOI &
                    & -( PFPLVN(JL,JK) - PFPLVN(JL,JK-1) ) * ZDZRHOI
      ZSLGEP(JL,JK) =   RLVTT * ( PFPLVL(JL,JK) - PFPLVL(JL,JK-1) ) * ZDZRHOI &
                    & + RLSTT * ( PFPLVN(JL,JK) - PFPLVN(JL,JK-1) ) * ZDZRHOI
    ENDDO
  ENDDO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      
      !account for updraft precipitation tendencies
      ZQTE(JL,JK)  = ZQTE(JL,JK)  + ZQTEP(JL,JK)
      ZSLGE(JL,JK) = ZSLGE(JL,JK) + ZSLGEP(JL,JK)
         
      !--wipe out any tendency above PBL--
      IF (LLTROPDAMP .AND. .NOT.LLPBL(JL,JK)) THEN
        ZQTE(JL,JK)  = ZQTE(JL,JK) - ZQTEA(JL,JK)
        ZSLGE(JL,JK) = ZSLGE(JL,JK) - ZSLGEA(JL,JK)
      ENDIF
      
      !integrate in time
      ZQTUPD(JL,JK)  = ZQTM1(JL,JK)  + ZQTE(JL,JK)  * ZTMST
      ZSLGUPD(JL,JK) = ZSLGM1(JL,JK) + ZSLGE(JL,JK) * ZTMST
      
      !total specific humidity limiter
      ZQTUPD(JL,JK) = MAX( 0._JPRB, ZQTUPD(JL,JK))
      ZQTE(JL,JK)   = (ZQTUPD(JL,JK) - ZQTM1(JL,JK) ) * ZRTMST
      
      IF (LLDIAG) THEN
        !--output of rad+dyn and conv+diff tendencies of QT and THL --
        PEXTRA(JL,JK,71) = ZSLGEA(JL,JK) * 3600._JPRB * 24._JPRB / RCPD
        PEXTRA(JL,JK,72) = (ZSLGE(JL,JK)-ZSLGEA(JL,JK)) * 3600._JPRB * 24._JPRB / RCPD

        PEXTRA(JL,JK,73) = ZQTEA(JL,JK) * 3600._JPRB * 24._JPRB * 1000._JPRB
        PEXTRA(JL,JK,74) = (ZQTE(JL,JK)-ZQTEA(JL,JK)) * 3600._JPRB * 24._JPRB * 1000._JPRB
      ENDIF  

    ENDDO
  ENDDO



!     ------------------------------------------------------------------

!*         7.     SURFACE FLUXES - TILES
!                 ----------------------
!*         AND    COMPUTE 2M TEMPERATURE AND HUMIDITY, 10M WIND,
!*                  and gustiness

!  Compute wind speed at blending height

CALL VDFFBLEND(KIDIA,KFDIA,KLON,KLEV, &
 & PUM1, PVM1, PGEOM1, PUCURR, PVCURR, ZBLEND, &
 & ZFBLEND)

! Wrap-up computations for the surface and 2T/2D/10U/10V/gustiness computation

CALL SURFPP( KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,KTILES=KTILES, &
 & KDHVTLS=KDHVTLS,KDHFTLS=KDHFTLS, &
 & PTSTEP=PTSPHY, &
! input
 & PFRTI=PFRTI, PAHFLTI=ZAHFLTI, PG0TI=ZG0, &
 & PSTRTULEV=PSTRTU(:,KLEV), PSTRTVLEV=PSTRTV(:,KLEV), PTSKM1M=PTSKM1M, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PQMLEV=PQM1(:,KLEV), &
 & PGEOMLEV=PGEOM1(:,KLEV), PCPTSPP=ZZCPTS, PCPTGZLEV=ZCPTGZ(:,KLEV), &
 & PAPHMS=PAPHM1(:,KLEV), PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
 & PZDL=ZZZDL, PQSAPP=ZZQSA, PBLEND=ZBLEND, PFBLEND=ZFBLEND, PBUOM=ZZBUOM, &
 & PZ0M=PZ0M, PEVAPSNW=PEVAPSNW,PSSRFLTI=PSSRFLTI, PSLRFL=PSLRFL, PSST=PSST, &
 & PUCURR=PUCURR, PVCURR=PVCURR, &
! updated
 & PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, PTSKE1=PTSKE1,PTSKTIP1=ZTSKTIP1, &
! output
 & PDIFTSLEV=PDIFTS(:,KLEV), PDIFTQLEV=PDIFTQ(:,KLEV), PUSTRTI=PUSTRTI, &
 & PVSTRTI=PVSTRTI,  PTSKTI=PTSKTI, PAHFLEV=PAHFLEV, PAHFLSB=PAHFLSB, &
 & PFWSB=PFWSB, PU10M=PU10M, PV10M=PV10M, PT2M=PT2M, PD2M=PD2M, PQ2M=PQ2M, &
 & PGUST=PGUST, PZIDLWV=PZIDLWV, &
! output DDH
 & PDHTLS=PDHTLS &
 & )

PDIFTL  (KIDIA:KFDIA,KLEV) = 0.0_JPRB
PDIFTI  (KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTQT (KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTSLG(KIDIA:KFDIA,KLEV) = 0.0_JPRB
ZDIFTQT (KIDIA:KFDIA,KLEV) = PDIFTQ(KIDIA:KFDIA,KLEV)
ZDIFTSLG(KIDIA:KFDIA,KLEV) = PDIFTS(KIDIA:KFDIA,KLEV)



!     ------------------------------------------------------------------

!*         8.     SLG, QT, U, V FLUX COMPUTATIONS AND T,SKIN TENDENCY
!                 ---------------------------------------------------

DO JL=KIDIA,KFDIA
  ZDIFTQT (JL,0) = 0.0_JPRB
  PDIFTQ  (JL,0) = 0.0_JPRB
  PDIFTL  (JL,0) = 0.0_JPRB
  PDIFTI  (JL,0) = 0.0_JPRB
  PDIFTS  (JL,0) = 0.0_JPRB
  ZDIFTSLG(JL,0) = 0.0_JPRB
  PSTRTU  (JL,0) = 0.0_JPRB
  PSTRTV  (JL,0) = 0.0_JPRB
ENDDO

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
!...change in slg,qt,u,v tendencies are converted to fluxes
    ZDIFTSLG(JL,JK) = ( ZSLGEWODIS(JL,JK+1) - ZSLGEA(JL,JK+1) ) * ZGDPH &
                    & + ZDIFTSLG(JL,JK+1)  
    ZDIFTQT(JL,JK)  = (ZQTE (JL,JK+1)-ZQTEA(JL,JK+1))*ZGDPH + ZDIFTQT(JL,JK+1)
  ENDDO
ENDDO



!     ------------------------------------------------------------------
! 
!*         9.     OLD SATURATION SPECIFIC HUMIDITY
!                 --------------------------------
!
!          This part will be redundant soon.
!          Only ZALFAW is still used to separate ice and liquid water.
!


  !-- state after dynamics and radiation --
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      
      ZAUPD(JL,JK)  = PAM1(JL,JK) + PAE(JL,JK) * ZTMST
      ZLUPD(JL,JK)  = PLM1(JL,JK) + PLE(JL,JK) * ZTMST
      ZIUPD(JL,JK)  = PIM1(JL,JK) + PIE(JL,JK) * ZTMST
      ZQUPD(JL,JK)  = PQM1(JL,JK) + PQE(JL,JK) * ZTMST
      ZTUPD(JL,JK)  = PTM1(JL,JK) + PTE(JL,JK) * ZTMST
!          total condensate (water + ice)
      ZLIUPD(JL,JK) = ZLUPD(JL,JK) + ZIUPD(JL,JK)
      
    ENDDO
  ENDDO


  !-- calculate saturation specific humidity --
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZQSVAR(JL,JK) = FOEEWM(ZTUPD(JL,JK))/PAPM1(JL,JK)
      ZQSVAR(JL,JK) = MIN(0.5_JPRB,ZQSVAR(JL,JK))
      ZQSVAR(JL,JK) = ZQSVAR(JL,JK)/(1.0_JPRB-RETV*ZQSVAR(JL,JK))

      !PEXTRA(JL,JK,75) = 1000._JPRB * ZQSVAR(JL,JK)

      ZQSVAR(JL,JK) = MAX(ZQSVAR(JL,JK),ZQUPD(JL,JK)) !don't allow input supersat.
      !PEXTRA(JL,JK,76) = 1000._JPRB * ZQSVAR(JL,JK)

!          dqsat/dT correction factor (1+L/cp*dqsat/dT) & alfa

      ZALFAW(JL,JK)=FOEALFA(ZTUPD(JL,JK))
      ZFACW=R5LES/((ZTUPD(JL,JK)-R4LES)**2)
      ZFACI=R5IES/((ZTUPD(JL,JK)-R4IES)**2)
      ZFAC=ZALFAW(JL,JK)*ZFACW+(1.0_JPRB-ZALFAW(JL,JK))*ZFACI
      ZESDP=FOEEWM(ZTUPD(JL,JK))/PAPM1(JL,JK)
      ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
      ZDQSDTEMP(JL,JK)=ZFAC*ZCOR*ZQSVAR(JL,JK)  !dqsat/dT
      ZCORQS(JL,JK)=MAX(1.0_JPRB,1.0_JPRB+FOELDCPM(ZTUPD(JL,JK))*ZDQSDTEMP(JL,JK))
    ENDDO
  ENDDO
  
  
  
!     ------------------------------------------------------------------

!*         10.    QT VARIANCE BUDGET
!                 ------------------
!
!          Flux-gradient production, transport and dissipation
!

  !-- recall variance from previous timestep ----
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      !PVAR(JL,JK) = PEXTRA(JL,JK,43)**2
      PVAR(JL,JK) = 0._JPRB
    ENDDO
  ENDDO
      
      
      
  DO JK=KLEV-1,2,-1
    DO JL=KIDIA,KFDIA
      
      
      !-- recall and update the variance dissipation timescale ZTAU --
!      ZTAU(JL,JK)    = PEXTRA(JL,JK,65)
!      IF (LLPBL(JL,JK)) THEN
!        ZTAU(JL,JK) =  PZINV(JL) / MAX(PWUAVG(JL),0.01_JPRB)
!        !ZTAU(JL,JK) =  2._JPRB * PZINV(JL) / MAX(PWUAVG(JL),0.01_JPRB)
!        
!        !  implicit
!        !ZTAUNEW = 2._JPRB * PZINV(JL) / MAX(PWUAVG(JL),0.01_JPRB)
!        !ZTAU(JL,JK) =  ZTAUNEW +  ( ZTAU(JL,JK) - ZTAUNEW ) * EXP( - ZTMST / ZTAU(JL,JK) )
!        
!        ZTAU(JL,JK) = MAX(ZTAU(JL,JK), 100.0_JPRB)
!      ELSE
!        ZTAU(JL,JK) = ZTAU(JL,JK) + ZTMST
!        !ZTAU(JL,JK) = 2.0E3_JPRB
!      ENDIF
!      ZTAU(JL,JK) = MIN( 1.0E3_JPRB, ZTAU(JL,JK) )
!      !ZTAU(JL,JK) = MIN( 2.0E3_JPRB, ZTAU(JL,JK) )
      
      ZTAU(JL,JK) =  PZINV(JL) / MAX(PWUAVG(JL),0.01_JPRB)

      
      !-- do the individual variance budget terms --
      ZVARDISS   = 0._JPRB
      ZVARGEN    = 0._JPRB
      ZVARTRANS  = 0._JPRB
      ZDQTDZ  = 0._JPRB
      ZWQTF   = 0._JPRB
      
      IF (LLPBL(JL,JK)) THEN
      
         !--- I  flux-gradient variance production at full level ---
        IF (KPBLTYPE(JL)==2 .AND. JK<=KVARTOP(JL)+1 ) THEN
          !for stratocumulus, protect variance production in top PBL layer against strong capping gradient
          ZDQTDZ = (ZQTUPD(JL,JK)-ZQTUPD(JL,JK+1)) * RG / (PGEOM1(JL,JK)-PGEOM1(JL,JK+1))
        ELSE
          !otherwise, do it truely centered
          ZDQTDZ = (ZQTUPD(JL,JK-1)-ZQTUPD(JL,JK+1)) * RG / (PGEOM1(JL,JK-1)-PGEOM1(JL,JK+1))
        ENDIF
          
        !IF ( JK==KVARTOP(JL) ) THEN
        !  ZWQTF = 0._JPRB
        !ELSE
          ZWQTF = (ZDIFTQT(JL,JK-1) + ZDIFTQT(JL,JK))/2._JPRB  
        !ENDIF
        ZWQTF = -(RD * PTM1(JL,JK) / PAPM1(JL,JK)) * ZWQTF
           
        ZVARGEN = -2._JPRB * ZWQTF * ZDQTDZ
        ZVARGEN = MAX(ZVARGEN,0._JPRB)   ! exclude countergradient flow
          
          
        !--- II   variance transport at full level ---
        ZVARTRANS = - (ZWQT2(JL,JK-1) - ZWQT2(JL,JK) ) * RG / (PGEOH(JL,JK-1)-PGEOH(JL,JK))


        !--- III  variance dissipation at full level (for output only) ---
        ZVARDISS = MIN(0._JPRB, -ZVARGEN -ZVARTRANS)       
           
!      ELSE
!      
!        !--- III  variance dissipation (for output only) ---
!        ZVARDISS = MIN(0._JPRB, -PVAR(JL,JK)/ZTAU(JL,JK) )
             
             
       !IF ( JK==KVARTOP(JL)+1 ) THEN
       !   ZVARGEN = ZVARGEN * 4._JPRB
       !ENDIF
       
       
      ENDIF
      
      
      !--- update the variance ---
!      !IF (PVAR(JL,JK)>0._JPRB) THEN
!      IF (PVAR(JL,JK)>1.0E-10_JPRB) THEN
!        !--- if variance exists at t-1, do it prognostically (implicit) ---
!        PVAR(JL,JK) = (ZVARGEN+ZVARTRANS) * ZTAU(JL,JK) + &
!                    & (PVAR(JL,JK) - (ZVARGEN+ZVARTRANS) * ZTAU(JL,JK) ) * &
!                    &     EXP( - ZTMST / ZTAU(JL,JK) )
!      ELSE
        !--- if no variance exists at t-1, do it diagnostically ---
        PVAR(JL,JK) = (ZVARGEN+ZVARTRANS) * ZTAU(JL,JK) 
!      ENDIF              

      PVAR(JL,JK) = MAX(PVAR(JL,JK),0._JPRB)
      
      IF (LLDIAG) THEN
        !PEXTRA(JL,JK,57) = 1000000._JPRB * (ZVARGEN+ZVARTRANS) * ZTAU(JL,JK)
        PEXTRA(JL,JK,58) = 1000000._JPRB * PVAR(JL,JK) 
        PEXTRA(JL,JK,59) = 1000000._JPRB * ZVARTRANS 
        PEXTRA(JL,JK,60) = 1000000._JPRB * ZVARGEN           
        PEXTRA(JL,JK,61) = 1000000._JPRB * ZVARDISS
        PEXTRA(JL,JK,62) = 1000._JPRB * ZDQTDZ
        PEXTRA(JL,JK,63) = 1000._JPRB * ZWQTF    
        PEXTRA(JL,JK,64) = 1000000._JPRB * ZWQT2(JL,JK) 
        PEXTRA(JL,JK,65) = ZTAU(JL,JK)
      ENDIF  
      
    ENDDO
  ENDDO
  
  
  !-- update the prognostic variance and variance timescale --
!!  IF (KCNT>0) THEN
!    DO JK=2,KLEV-1
!    DO JL=KIDIA,KFDIA
!!      PEXTRA(JL,JK,65) = ZTAU(JL,JK)
!      PEXTRA(JL,JK,43) = SQRT( PVAR(JL,JK) )
!    ENDDO
!    ENDDO
!!  ENDIF


      
!     ------------------------------------------------------------------
!
!*         11.    VECTORIZED BIMODAL CLOUD SCHEME
!                 -------------------------------
!          
!          The EDMF decomposition is extended into the cloud scheme, by doing a bimodal
!          PDF: one diffusive, one updraft. Each PDF is Gaussian, their 1st and 2nd moments are
!          parameterized. The moments of all PDFs are related, see Lewellen and Yoh (JAS, 1993).
!          The scheme is formulated in {thl,qt} space, using vector calculus.
! 
!          Variance closures:
!            * The overall variance is done through the full budget (see section 10).
!            * The updraft PDF variance is done using properties of multiple updrafts.
!
!          See Neggers (JAS, 2009)
!

  !  use new {QT,SLG} state
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      
      ZQTTEST(JL,JK) = 0._JPRB

      SELECT CASE (KPBLTYPE(JL))
      
        CASE (0,1,4)
          ZSLGUP(JL,JK)  = ZSLGUPD(JL,JK)
          ZQTUP(JL,JK)   = ZQTUPD(JL,JK)
          ZQTTEST(JL,JK) = ZQTUPD(JL,JK)
          
        CASE (2,3)
          IF (LLPBL(JL,JK) .AND. ZSLGUH(JL,JK,3)>0._JPRB .AND. ZSLGUH(JL,JK-1,3)>0._JPRB ) THEN
            !-- interpolate the updraft fields to full levels --
            ZSLGUP(JL,JK)  = ( ZSLGUH(JL,JK,3)+ ZSLGUH(JL,JK-1,3) ) / 2._JPRB + ZSLGE(JL,JK)
            ZQTUP(JL,JK)   = ( ZQTUH(JL,JK,3) + ZQTUH(JL,JK-1,3)  ) / 2._JPRB + ZQTE(JL,JK)
            ZQTTEST(JL,JK) = ( ZQTUH(JL,JK,1) + ZQTUH(JL,JK-1,1)  ) / 2._JPRB + ZQTE(JL,JK)
          ELSE
            ZSLGUP(JL,JK)  = ZSLGUPD(JL,JK)
            ZQTUP(JL,JK)   = ZQTUPD(JL,JK)
            ZQTTEST(JL,JK) = ZQTUPD(JL,JK)
          ENDIF
          
      END SELECT
      
      ZRHO = PAPHM1(JL,JK)/( RD*PTM1(JL,JK)*(1.0_JPRB+RETV*PQM1(JL,JK)) )
      ZAWAKE(JL,JK) = MIN(1._JPRB, MAX(0._JPRB, ZDETR(JL,JK) * ZTAU(JL,JK) / ZRHO ) )
        
    ENDDO
  ENDDO
  
  
  CALL VDFCLOUD ( KIDIA     , KFDIA   , KLON    , KLEV   , IDRAFT , &
                & PAPM1     , PGEOM1  , PGEOH   , &
                & ZQTUPD    , ZSLGUPD , &
                & ZFRACB    , KVARTOP , &
                & ZQTUP     , ZSLGUP  , &
                & ZQTTEST   , &
                & PVAR      , ZAWAKE  , &
! DIAGNOSTIC OUTPUT
                & PEXTR2    , KFLDX2  , PEXTRA  , KLEVX  , KFLDX  , &
!              
                & ZCLDFRAC  , ZQLAV   , PLDIFF)



!     ------------------------------------------------------------------
!
!*         12.    NET TENDENCIES
!                 --------------
!
!          Calculate net tendencies of the prognostic model variables QL, QI, A, QV and T


  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      

      !RN --- temporarily switch off cloudiness, for testing ---
      !ZCLDFRAC(JL,JK) = 0._JPRB
      !ZQLAV(JL,JK)    = 0._JPRB
      
      
      !--- Convert back from conserved variables QT and SLG to non-conserved QV and T ---
      !
      !  New cloud variables (liquid, ice, fraction) for tendency calculation:
      !
      !    Within PBL: use new values from VDFCLOUD
      !    Above PBL:  use t-1 input values + tendencies
      !
      !    Result: * above PBL, VDFMAIN does not change cloud variables:
      !                tendencies only contain contributions by rad+dyn.
      !
      IF (LLPBL(JL,JK)) THEN    !reset cloudiness within PBL
      !IF (.TRUE.) THEN          !reset cloudiness everywhere! (also above PBL)
      
        !-- safety: total condensate can not be larger than total specific humidity ---
        ZQLAV(JL,JK) = MIN( ZQTUPD(JL,JK) ,ZQLAV(JL,JK) )
        ZAUPD(JL,JK) = ZCLDFRAC(JL,JK)
        
        !-- decomposition of total condensate into ice and liquid ---
        !ZALFAW(JL,JK) = FOEALFA(PTM1(JL,JK))   !new alpha?
        ZLUPD(JL,JK) = ZQLAV(JL,JK) * ZALFAW(JL,JK) 
        ZIUPD(JL,JK) = ZQLAV(JL,JK) * ( 1.0_JPRB - ZALFAW(JL,JK))
        
      ELSE
      
        !-- outside PBL, maintain tendencies from rad+dyn ---
        ZAUPD(JL,JK) = PAM1(JL,JK) + PAE(JL,JK)*ZTMST
        ZLUPD(JL,JK) = PLM1(JL,JK) + PLE(JL,JK)*ZTMST
        ZIUPD(JL,JK) = PIM1(JL,JK) + PIE(JL,JK)*ZTMST
        
      ENDIF  
      
      
      !--- Derive non-conserved properties QV and T ---
      ZQUPD(JL,JK)  = ZQTUPD(JL,JK) - ZLUPD(JL,JK) - ZIUPD(JL,JK)
      ZTUPD(JL,JK)  = ( ZSLGUPD(JL,JK) - PGEOM1(JL,JK) &
        &     + RLVTT * ZLUPD(JL,JK) + RLSTT * ZIUPD(JL,JK) &
        &   ) / ( RCPD * ( 1.0_JPRB + RVTMP2 * ZQUPD(JL,JK) ) )   !compare to T->SLG conversion in section 2.2


      !--- Calculate the final tendencies between state at t-1 and state after rad + dyn + pbl ---
      PQE(JL,JK) = ( ZQUPD(JL,JK) - PQM1(JL,JK) ) * ZRTMST
      PTE(JL,JK) = ( ZTUPD(JL,JK) - PTM1(JL,JK) ) * ZRTMST
      PLE(JL,JK) = ( ZLUPD(JL,JK) - PLM1(JL,JK) ) * ZRTMST
      PIE(JL,JK) = ( ZIUPD(JL,JK) - PIM1(JL,JK) ) * ZRTMST
      PAE(JL,JK) = ( ZAUPD(JL,JK) - PAM1(JL,JK) ) * ZRTMST


      !--- T-check ---
      IF ( ZTUPD(JL,JK).gt.400._JPRB.or.ZTUPD(JL,JK).lt.100._JPRB) THEN
      
        write(0,'(a,3i5)')    'vdfmain T alarm:',KCNT,JL,JK
        write(0,'(a,3i5)')    '               : ', KPBLTYPE(JL),KVARTOP(JL),KHPBLN(JL)
        write(0,'(a,2f10.6)')    '        sfluxes: ', ZKHFL(JL),ZKQFL(JL)
        write(0,'(a,2f10.6)')    '    cld heights: ', ZCLDBASE(JL),ZCLDTOP(JL)
        write(0,'(a,4f10.6)')    '             T: ',&
           & PTM1(JL,JK),ZTUPD(JL,JK),PTE(JL,JK)*ZTMST,ZTEA(JL,JK)*ZTMST
        write(0,'(a,5f10.6)')    '           SLG: ',&
           & ZSLGM1(JL,JK)/RCPD       , ZSLGUPD(JL,JK)/RCPD      , &
           & ZSLGE(JL,JK)*ZTMST/RCPD  , ZSLGEA(JL,JK)*ZTMST/RCPD , &
           & ZSLGEP(JL,JK)*ZTMST/RCPD
        write(0,'(a,5f10.6)')    '            QT: ',&
           & ZQTM1(JL,JK)*1000._JPRB        , ZQTUPD(JL,JK)*1000._JPRB      , &
           & ZQTE(JL,JK)*ZTMST*1000._JPRB   , ZQTEA(JL,JK)*ZTMST*1000._JPRB , &
           & ZQTEP(JL,JK)*ZTMST*1000._JPRB
           
        IF (JK<=KLEV-1) THEN
          ZMGEOM  = PGEOM1(JL,JK)-PGEOM1(JL,JK+1)      
          ZCFNC1  = RVDIFTS * ZTMST * RG**2 * PAPHM1(JL,JK) &
                   & /( ZMGEOM * RD * 0.5_JPRB &
                   & *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  )) &
                   & +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1))))  
          ZDUMMY2 = ZMFLX(JL,JK,2) / ( ZCFNC1 * ZMGEOM * ZRG )
          ZDUMMY3 = ZMFLX(JL,JK,3) / ( ZCFNC1 * ZMGEOM * ZRG )
        ELSE
          ZDUMMY2 = 0.
          ZDUMMY3 = 0.
	ENDIF
	
	write(0,'(a,6f10.6)')    '          UPDR: ',&
!           & ZMFLX(JL,JK,2), ZSLGUH(JL,JK,2)/RCPD, ZQTUH(JL,JK,2)*1000._JPRB, &
!           & ZMFLX(JL,JK,3), ZSLGUH(JL,JK,3)/RCPD, ZQTUH(JL,JK,3)*1000._JPRB
           & ZDUMMY2, ZSLGUH(JL,JK,2)/RCPD, ZQTUH(JL,JK,2)*1000._JPRB, &
           & ZDUMMY3, ZSLGUH(JL,JK,3)/RCPD, ZQTUH(JL,JK,3)*1000._JPRB

        DO JKK= KLEV-5,KLEV
          write(0,'(a,i5,4f10.6)')    '      M-STRUCT: ',JKK, ZMFLX(JL,JKK,2), ZWUH(JL,JKK,2), ZMFLX(JL,JKK,3), ZWUH(JL,JKK,3)
        ENDDO
	
      ENDIF
      
      
      !--- q-check ----
      IF (ZQTUPD(JL,JK).lt.0._JPRB.OR.ZQTUPD(JL,JK).gt.0.05_JPRB) THEN
        write(0,'(3i5,a,4f10.6,a,3i5,2f10.6)') KCNT,JL,JK,'  terror! q=',&
           & ZQTUPD(JL,JK),ZQTM1(JL,JK),ZLUPD(JL,JK),PLM1(JL,JK),' - ',&
           & KVARTOP(JL),KHPBLN(JL),KPBLTYPE(JL),ZCLDBASE(JL),ZCLDTOP(JL)
      ENDIF

    ENDDO
  ENDDO



!RN     ---------- add a perturbation to ML ---------------

  IF (LLPERTURB .AND. KSTEP==KPERTURB .AND. KCNT==1 ) THEN
  
    write(0,*) '-------- ML perturbation -------- ',KSTEP,KCNT
       
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
      
!        IF ( JK>=KHPBLN(JL) ) THEN
        IF ( JK>=IPLCL(JL,3) ) THEN
          !humidity perturbation
          PQE(JL,JK) = PQE(JL,JK) + 0.001_JPRB * ZRTMST
        
          !temperature perturbation: preserve Bowen ratio
          !PTE(JL,JK) = PTE(JL,JK) + PQE(JL,JK) * ZKHFL(JL) / ZKQFL(JL)
          
          !write(0,*) '      ', 0.001, ZKHFL(JL) / ZKQFL(JL), 0.001 * ZKHFL(JL) / ZKQFL(JL)
          
        ENDIF  

      ENDDO
    ENDDO
     
  ENDIF
  
  
  
!     ----------- output --------------------------------------

!  IF (LLDIAG) THEN
!    DO JL=KIDIA,KFDIA
!      PEXTRA(JL,1:KLEV,67) = 100._JPRB * ZCLDFRAC(JL,:)
!      PEXTRA(JL,1:KLEV,69) = 1000._JPRB * ZQLAV(JL,:)
!      PEXTRA(JL,:,54) = ZDIFTQT(JL,:)
!      PEXTRA(JL,:,55) = ZDIFTSLG(JL,:)
!    ENDDO
!  ENDIF
  
      

!     ------------------------------------------------------------------

!*         13.    FLUX COMPUTATIONS
!                 -----------------
!


!*         13.1    Q, QL, QI AND S FLUX

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    ZGDPH = - (PAPHM1(JL,JK)-PAPHM1(JL,JK+1)) * ZRG
!...changes in q,l,i tendencies are converted to fluxes
    PDIFTQ(JL,JK) = (PQE (JL,JK+1) - ZQEA(JL,JK+1)) * ZGDPH + PDIFTQ(JL,JK+1)
    PDIFTL(JL,JK) = (PLE (JL,JK+1) - ZLEA(JL,JK+1)) * ZGDPH + PDIFTL(JL,JK+1)
    PDIFTI(JL,JK) = (PIE (JL,JK+1) - ZIEA(JL,JK+1)) * ZGDPH + PDIFTI(JL,JK+1)
!...slg=s-Lc*ql-Ld*qi (same for fluxes)
    PDIFTS(JL,JK) = ZDIFTSLG(JL,JK) &
                & + RLVTT * PDIFTL(JL,JK) + RLSTT * PDIFTI(JL,JK)  
  ENDDO
ENDDO

!*         13.2  PBL PRECIPITATION ENTHALPY FLUXES

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PFHPVL(JL,JK) = -RLVTT*PFPLVL(JL,JK)
    PFHPVN(JL,JK) = -RLSTT*PFPLVN(JL,JK)
  ENDDO    
ENDDO    


IF (LHOOK) CALL DR_HOOK('VDFMAIN',1,ZHOOK_HANDLE)
END SUBROUTINE VDFMAIN
