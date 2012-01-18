!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFHGHTN (KIDIA   , KFDIA   , KLON    , KLEV   , KDRAFT  , PTMST  , KSTEP    , &
                   & PUM1    , PVM1    , PTM1    , PQM1   , PLM1    , PIM1   , PAM1     ,&
                   & PAPHM1  , PAPM1   , PGEOM1  , PGEOH  , PVERVEL , PQE    , PTE      , &
                   & PKMFL   , PKHFL   , PKQFL   , PMFLX  , &
                   & PEXTR2  , KFLDX2  , PEXTRA  , KLEVX  , KFLDX , &
                   & PUUH    , PVUH    , PSLGUH  , PQTUH  , PFRACB  , PWUH  , &
                   & PZPTOP  , KPTOP   , PZPLCL  , KPLCL  , KPLZB   , &
                   & PWUAVG  , PRICUI  , PMCU    , PDTHV  , &
                   & PFPLVL  , PFPLVN  , PDETR   , &
                   & PBIR    , LDNODECP, LDRUNDRY, KPBLTYPE, PWQT2)
                   
                     
!     ------------------------------------------------------------------

!**   *VDFHGHTN* - DETERMINES THE PBL-HEIGHT AND STRONG UPDRAFT FIELDS
!                  USING A ENTRAINING PARCEL ASCENT METHOD.

!     A.P. SIEBESMA    30/06/99    Original (dry)
!     M. Ko"hler        3/12/2004  Moist Version
!     Roel Neggers     12/04/2005  Multiple updraft extension


!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT AND UPDRAFT FIELDS

!     INTERFACE
!     ---------

!     *VDFHGHTN* IS CALLED BY *VDFMAIN*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KDRAFT*       NUMBER OF EXPLICITLY MODELED DRAFTS - CURRENTLY 3:
!                    1: test parcel
!                    2: rising dry thermals which stop at cloud base or inversion
!                    3: rising dry thermals which become cloudy
!                    (4: downdrafts .. to be done?)

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)        S
!     *PUM1*         X-VELOCITY COMPONENT AT T-1                  M/S
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1                  M/S
!     *PTM1*         TEMPERATURE AT T-1                           K
!     *PQM1*         SPECIFIC HUMUDITY AT T-1                     KG/KG
!     *PLM1*         SPECIFIC CLOUD LIQUID WATER AT T-1           KG/KG
!     *PIM1*         SPECIFIC CLOUD ICE AT T-1                    KG/KG
!     *PAM1*         CLOUD FRACTION AT T-1                        KG/KG
!     *PAPHM1*       PRESSURE AT HALF LEVEL AT T-1                PA
!     *PAPM1*        PRESSURE AT FULL LEVEL AT T-1                PA
!     *PGEOM1*       GEOPOTENTIAL AT T-1                          M2/S2
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                   M2/S2
!     *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2  
!     *PKHFL*        SURFACE KINEMATIC HEAT FLUX                  K*M/S
!     *PKQFL*        SURFACE KINEMATIC MOISTURE FLUX              M/S
!     *PBIR*         BUOYANCY-FLUX INTEGRAL RATIO (-N/P)
!                    USED FOR DECOUPLING CRITERIA

!     *PVERVEL*      VERTICAL VELOCITY

!     INPUT PARAMETERS (LOGICAL):

!     *LDNODECP*     TRUE:  NEVER DECOUPLE
!                    FALSE: MAYBE DECOUPLE
!     *LDRUNDRY*     TRUE:  RUN PARCEL WITHOUT CONDENSATION
!                    FALSE: RUN PARCEL WITH CONDENSATION

!     OUTPUT PARAMETERS (REAL):

!     *PFPLVL*       PBL PRECIPITATION FLUX AS RAIN                KG/(M**2*S)
!     *PFPLVN*       PBL PRECIPITATION FLUX AS SNOW                KG/(M**2*S)

!     *PUUH*         UPDRAFT X-MOMENTUM
!     *PVUH*         UPDRAFT Y-MOMENTUM
!     *PSLGUH*       UPDRAFT GENERALIZED LIQUID STATIC ENERGY (SLG)
!                    AT HALF LEVEL                                M2/S2
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL   KG/KG
!     *PMFLX*        PBL MASS FLUX                                M/S
!     *PZPLCL*       HEIGHT OF LIFTING CONDENSATION LEVEL OF UPDRAFT          M
!     *PZPTOP*       HEIGHT OF LEVEL OF ZERO KINETIC ENERGY (W=0) OF UPDRAFT  M

!     OUTPUT PARAMETERS (INTEGER):

!     *KPLCL*         FIRST HALF LEVEL ABOVE REAL HEIGHT OF UPRAFT LCL
!     *KPTOP*         HIGHEST HALF LEVEL BELOW PZTOP, AND
!                       UPDRAFT TOP FULL LEVEL (PZTOP IS WITHIN THAT LAYER)
!     *KPLZB*         LEVEL OF UPRAFT ZERO BUOYANCY (LAST FULL LEVEL THAT IS POS. BUOYANT)
!     *KPBLTYPE*    -1: not defined yet
!                    0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: stratocumulus
!                    3: shallow cumulus
!                    4: deep cumulus

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

use garbage, only : foeewm, surf_inq
USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RG       ,RD       ,RCPD     ,RETV     ,RLVTT    ,&
                    & RLSTT    ,RATM     ,RTT      ,RLMLT

USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
                    & R4IES    ,R5LES    ,R5IES    ,RVTMP2   ,R5ALVCP  ,&
                    & R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,&
                    & RTICECU  ,RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
	      
USE YOEVDF   , ONLY : RKAP     ,RVDIFTS  ,LLDIAG
USE YOECUMF  , ONLY : RTAUMEL

!USE YOMGF1C  , ONLY : NC
!USE YOMLOG1C , ONLY : LCCN

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KPLCL(KLON,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KPTOP(KLON,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KPLZB(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PMFLX(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)     :: PUUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)     :: PVUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PSLGUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PQTUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PWUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(OUT)     :: PFRACB(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PZPLCL(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PZPTOP(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PWUAVG(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDETR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBIR(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRICUI(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDTHV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMCU(KLON) 
LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON) 
!ldrundry not used now
LOGICAL           ,INTENT(IN)    :: LDRUNDRY(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT)   :: KPBLTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWQT2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQE(KLON,KLEV) 

!diagnostic output
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)



!*         0.2    LOCAL VARIABLES

!--- mean & environmental properties ---
REAL(KIND=JPRB) ::    ZUSTAR  (KLON)     , ZWSTAR(KLON)       , ZKHVFL(KLON)       , &
                    & ZUSIGMA(KLON)      , ZWSIGMA(KLON)      , &
                    & ZSLGENH(KLON,0:KLEV),ZQLENH(KLON,0:KLEV), ZQIENH(KLON,0:KLEV), &
                    & ZQTENH(KLON,0:KLEV), ZUENH(KLON,0:KLEV) , ZVENH(KLON,0:KLEV) , &
                    & ZTVEN(KLON,KLEV)   , ZQTM1 (KLON,KLEV)  , &
                    & ZSLGM1(KLON,KLEV)  , ZMGEOM(KLON,0:KLEV), &
                    & ZTENH(KLON,0:KLEV) , ZRHOH (KLON,0:KLEV), ZTHVEN(KLON,KLEV)

!--- updraft parameters ---
REAL(KIND=JPRB) ::    ZWU2H (KLON,0:KLEV,KDRAFT), ZWUH, &
                    & ZQCUH (KLON,0:KLEV,KDRAFT), ZQUH  (KLON,0:KLEV,KDRAFT), &
		    & ZTUH  (KLON,0:KLEV,KDRAFT), ZEPS  (KLON,0:KLEV,KDRAFT), &
		    & ZFRAC (KLON,0:KLEV,KDRAFT), &
		    & ZBUOF (KLON,KLEV,KDRAFT)  , ZMCLD (KLON)              , &
                    & ZABULK(KLON,0:KLEV) , ZWBULK(KLON,0:KLEV)  , &
                    & ZQTBULK(KLON,0:KLEV), ZSLGBULK(KLON,0:KLEV), &
                    & ZUBULK(KLON,0:KLEV) , ZVBULK(KLON,0:KLEV)  , &
                    & ZZPLZB(KLON,KDRAFT) , ZCAPE1(KLON)
                    
REAL(KIND=JPRB) ::    ZQSATM, ZSATDEF, &
                    & ZUPFLXL(KLON,0:KLEV,KDRAFT), ZUPFLXN(KLON,0:KLEV,KDRAFT), &
                    & ZUPGENL(KLON,KLEV,KDRAFT), ZUPGENN(KLON,KLEV,KDRAFT), &
                    & ZALFAW, ZDZRHO, ZPFLXTOT, ZPEVAPUP, ZFAC, ZUPMELT, &
                    & ZAPRECEVAP, ZQTEP

REAL(KIND=JPRB) ::    ZFRACB(KLON,KDRAFT), ZMFLXB(KLON,KDRAFT)
      
REAL(KIND=JPRB) ::    ZFRACMAX , ZFACMAXEXC , ZFRACTEST , ZFACTESTEXC , &
                    & ZFACEXC(KLON,KDRAFT), ZDUMFRAC, ZDUMR, ZMASSCAPDEPTH, &
                    & ZPDFFACPHI(KLON), ZPDFFACW(KLON), ZLOBUKHOV, &
                    & ZSTABILITY(KLON), ZCOUPLING(KLON)

LOGICAL ::            LLDONE(KLON,KDRAFT), LLMASSCAP, LLMCIND(KLON), & 
                    & LLWIPE, LLSTCU


INTEGER(KIND=JPIM) :: IS, JK, JL, JD, IBASE, JKH, ITOP

REAL(KIND=JPRB) ::    ZQEXC   , ZTEXC   , ZDZ     ,  ZDB   , &
                    & ZSPEEDENV         , ZSPEEDUP, ZWINDIR , &
                    & ZZ      , ZTOTW2(KLON)      , ZTOTP(KLON) , &
                    & ZCONS10 , ZCFNC1(KLON,0:KLEV)      , ZTVMEAN     , &
                    & ZRG     , ZMFMAX  , ZMFS(KLON,KDRAFT)

!          REMAINING MODEL PARAMETERS

REAL(KIND=JPRB) ::    ZTAUEPS , ZCLDDEPTH     , &
                    & ZW2THRESH               , ZSTABTHRESH    , ZBIRTHRESH , &
                    & ZTVLIM  , ZCLDDEPTHDP   , ZDZCLOUD(KLON) , ZW2H       , &
                    & ZZFUNC  , ZZFUNC3(KLON) , ZZFUNC4(KLON)  , ZDHRI(KLON)    , ZDHCL      , &
                    & ZDTHVDZ , ZREPUST       , ZGHM1
		    
REAL(KIND=JPRB) ::    ZZI(KLON)

INTEGER(KIND=JPIM) :: IZI(KLON,KDRAFT)

REAL(KIND=JPRB) ::    ZBUOYCU, ZDTHVCUTOP, ZRICUINORM, ZDMDZ

REAL(KIND=JPRB) ::    ZFRACBPLUS(KLON,0:KLEV), ZWUHBPLUS(KLON,0:KLEV), &
                    & ZQTUHBPLUS(KLON,0:KLEV), ZSLGUHBPLUS(KLON,0:KLEV)

LOGICAL ::            LLCAPE, LLCAPETEST

REAL(KIND=JPRB) ::    ZRHS(KLON), ZMASSCAPE(KLON), ZTAUBM, &
                    & ZFRACBCONG(KLON), ZMSCALE(KLON)
		    
REAL(KIND=JPRB) :: ZHOOK_HANDLE


INTERFACE
#include "surf_inq.h"
END INTERFACE


#include "vdfparcel.intfb.h"
#include "vdfstcucrit.intfb.h"
#include "vdfpdftable.intfb.h"
#include "vdfbuoysort.intfb.h"

#include "fcttre.h"



!     ------------------------------------------------------------------

!*         1.     INITIALIZATION
!                 --------------

IF (LHOOK) CALL DR_HOOK('VDFHGHTN',0,ZHOOK_HANDLE)

!if (LCCN) then
!  write(0,*) 'vdfhghtn: nc=',nc
!endif


!-- top % of the PDF associated with the test parcel --
!ZFRACTEST   = 0.0002_JPRB    
!ZFRACTEST   = 0.001_JPRB    
ZFRACTEST   = 0.002_JPRB      !cy32r3 
!ZFRACTEST   = 0.005_JPRB    
CALL VDFPDFTABLE (ZFRACTEST, ZFACTESTEXC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor

!-- total convective area fraction that is done with mass flux --
!ZFRACMAX    = 0.05_JPRB     
ZFRACMAX    = 0.075_JPRB     !cy32r3
!ZFRACMAX    = 0.1_JPRB     
!ZFRACMAX    = 0.15_JPRB     
CALL VDFPDFTABLE (ZFRACMAX, ZFACMAXEXC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor

!-- eddy turnover time scale used in parcel entrainment [s]  (Neggers, Siebesma & Jonker, JAS 2002) --
!ZTAUEPS     = 300._JPRB    
ZTAUEPS     = 400._JPRB    ! cy32 
!ZTAUEPS     = 500._JPRB    

!-- threshold parcel vertical velocity squared [m2/s2] --
!ZW2THRESH  = -1._JPRB     
ZW2THRESH   = 0._JPRB      

!-- threshold cloud thickness for stcu/cu transition [m] --
ZCLDDEPTH   = 2000._JPRB   

!-- threshold cloud thickness used in shallow/deep decision [m] --
!ZCLDDEPTHDP = 2000._JPRB   
!ZCLDDEPTHDP = 3000._JPRB   !cy32
ZCLDDEPTHDP = 100000._JPRB   

ZSTABTHRESH = 20._JPRB     ! threshold stability (Klein & Hartmann criteria) [K]
ZBIRTHRESH  = 0.1_JPRB     ! threshold BIR (TKE decoupling criteria) [1]
ZTVLIM      = 0.1_JPRB     ! cloud fraction limit in Tv,env calculation

CALL SURF_INQ(PREPUST=ZREPUST)
                 
!-- switch for moist mass flux depth limiter --      
!LLMASSCAP     = .TRUE.
LLMASSCAP     = .FALSE.    
ZMASSCAPDEPTH = 3000._JPRB 
!ZMASSCAPDEPTH = 50000._JPRB 

!-- switch for applying Klein-Hartmann criterion for stratocumulus --
!LLSTCU = .TRUE.
LLSTCU = .FALSE.

!-- factor used in updraft initialization --
DO JL=KIDIA,KFDIA
  ZPDFFACW(JL)   = 0.2_JPRB
  ZPDFFACPHI(JL) = 0.2_JPRB
  !ZPDFFACW(JL)   = 1.0_JPRB
  !ZPDFFACPHI(JL) = 1.0_JPRB
ENDDO  

!-- updraft precip evaporation constant --
ZAPRECEVAP = 0.001_JPRB
!ZAPRECEVAP = 0.000544_JPRB
!ZAPRECEVAP = 0.0001_JPRB

!-- settings for Fritsch-Chappell CAPE-removal --
!LLCAPE     = .TRUE.                  !activation switch
LLCAPE     = .FALSE.
ZTAUBM     = 3600._JPRB * 1._JPRB    !the associated CAPE adjustment timescale
!ZTAUBM     = 3600._JPRB * 2._JPRB
LLCAPETEST = .TRUE.                  !Option I: test updraft carries the required additional transport
!LLCAPETEST = .FALSE.                 !Option II: moist updraft carries the required additional transport
  
!-- optimization --
ZRG    = 1.0_JPRB/RG


! set some stuff to zero
DO JL=KIDIA,KFDIA
  
  PWUAVG(JL)     = 0.0_JPRB
  KPBLTYPE(JL)   = -1          ! -1 means: yet unknown
  
  ZZI(JL)        = 0._JPRB      !mixed layer scalings
  ZWSTAR(JL)     = 0._JPRB        
  
  PRICUI(JL)  = 1._JPRB       ! 1 / cumulus inversion Richardson number
  PDTHV(JL)   = 0._JPRB
   
  ZCAPE1(JL)     = 0._JPRB

  ZMCLD(JL)       = 0._JPRB
  PMCU(JL)        = 0._JPRB       !cloud-depth average moist updraft mass flux
  
  ZCOUPLING(JL) = 0._JPRB
  
  ZRHS(JL)      = 0._JPRB 
  
ENDDO


DO JD=1,KDRAFT
  DO JL=KIDIA,KFDIA
    PZPLCL(JL,JD)  = -100._JPRB  ! default value: -100 (no LCL)
    PZPTOP(JL,JD)  = 0._JPRB     
    KPLCL(JL,JD)   = 0           ! default value: 0 (no PBL cloud)
    KPTOP(JL,JD)   = 0          
    KPLZB(JL,JD)   = 0          
    LLDONE(JL,JD)  = .TRUE.       ! default: TRUE (don't launch the parcel)
    ZFRACB(JL,JD)  = 0._JPRB 
    PFRACB(JL,JD)  = 0._JPRB 
    ZFACEXC(JL,JD) = 0._JPRB 
    ZFACEXC(JL,JD) = 0._JPRB 
    ZMFLXB(JL,JD)  = 0._JPRB 
    IZI(JL,JD)     = 0._JPRB        
  ENDDO
ENDDO


DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    PWQT2(JL,JK) = 0._JPRB  
  ENDDO
ENDDO


DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PDETR(JL,JK) = 0._JPRB
  ENDDO
ENDDO    


!--- parcel half level parameters ---
DO JD=1,KDRAFT
  DO JK=0,KLEV
    DO JL=KIDIA,KFDIA
    PUUH(JL,JK,JD)  = 0.0_JPRB
    PVUH(JL,JK,JD)  = 0.0_JPRB
    PSLGUH(JL,JK,JD)= 0.0_JPRB
    PQTUH(JL,JK,JD) = 0.0_JPRB
    PMFLX(JL,JK,JD) = 0.0_JPRB
    ZTUH(JL,JK,JD)  = 0.0_JPRB
    ZQUH(JL,JK,JD)  = 0.0_JPRB
    ZQCUH(JL,JK,JD) = 0.0_JPRB
    ZEPS(JL,JK,JD)  = 0.0_JPRB
    ZWU2H(JL,JK,JD) = 0.0_JPRB
    ZFRAC(JL,JK,JD) = 0.0_JPRB
    ZUPFLXL(JL,JK,JD)  = 0.0_JPRB
    ZUPFLXN(JL,JK,JD)  = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO


!--- parcel full level parameters ---
DO JD=1,KDRAFT
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZBUOF(JL,JK,JD)    = 0.0_JPRB
      ZUPGENL(JL,JK,JD)  = 0.0_JPRB
      ZUPGENN(JL,JK,JD)  = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO


!--- reset output stuff ---
IF (LLDIAG) THEN

  DO JL=KIDIA,KFDIA
    PEXTR2(JL,1:49) = 0._JPRB
  ENDDO

  DO JK=1,KLEVX
  DO JL=KIDIA,KFDIA
    PEXTRA(JL,JK,49) = 0._JPRB
    PEXTRA(JL,JK,50) = 0._JPRB
  ENDDO
  ENDDO

ENDIF


!     -----------------------------------------------------------------
!
!*         2.     PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!*                OF CONSERVED VARIABLES
!                 -----------------------------------------------------

!*         2.1  FULL LEVEL CPM, SLG, QT AND TV
!*

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZSLGM1(JL,JK) = RCPD * PTM1(JL,JK) + PGEOM1(JL,JK) &
                  & - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)  
      ZQTM1 (JL,JK) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)

!          parcel goes through cloud portion of environment
!          (added ql loading; ql,cld=ql,mean/fc; qv = qsat) 
!          safety: fc>0.1; linear interpolation between overcast 
!                  and cloudy portion for 0<fc<0.1
!                  guaranteed to be < tv from mean conditions

!          grid box mean virtual effect
      ZTVMEAN       = PTM1(JL,JK) * ( 1.0_JPRB + RETV * PQM1(JL,JK) &
                  & - PLM1(JL,JK) - PIM1(JL,JK) )       !qli loading  
      ZTVEN(JL,JK) = ZTVMEAN
      ZTHVEN(JL,JK) = ( PAPM1(JL,JK)/RATM )**(-RD/RCPD) * ZTVEN(JL,JK)
    ENDDO
  ENDDO


!*         2.2  HALF-LEVEL ENVIRONMENT INTERPOLATION (QT, QL, QI, SLG)
!*

  DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA

    IF (JK==1) THEN
      ZGHM1 = PGEOH(JL,JK) + 50000._JPRB*RG   !avoid using top half level (=inf)
    ELSE
      ZGHM1 = PGEOH(JL,JK-1)
    ENDIF  
    
    ZQTENH(JL,JK) = ( ZQTM1(JL,JK+1) *(ZGHM1-PGEOH(JL,JK  )) &
                & +   ZQTM1(JL,JK)   *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZQLENH(JL,JK) = ( PLM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PLM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZQIENH(JL,JK) = ( PIM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PIM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZSLGENH(JL,JK)= ( ZSLGM1(JL,JK+1)*(ZGHM1-PGEOH(JL,JK  )) &
                & +   ZSLGM1(JL,JK)  *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZUENH(JL,JK)  = ( PUM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PUM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))
    ZVENH(JL,JK)  = ( PVM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                & +   PVM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(ZGHM1-PGEOH(JL,JK+1))

    !calculate T at half levels from sl, for later use in density calculations		
    !ZTENH        = ( ZSLGENH (JL,JK) - PGEOH(JL,JK) &
    !                   & + RLVTT*ZQLENH(JL,JK) + RLSTT*ZQIENH(JL,JK) &
    !                   & ) / RCPD
    
    ZTENH(JL,JK)  =  ( PTM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                 & +   PTM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                 &   )                /(ZGHM1-PGEOH(JL,JK+1))
                
    !air density at half levels
    ZRHOH(JL,JK) = PAPHM1(JL,JK)/(RD*ZTENH(JL,JK))
      
  ENDDO
  ENDDO



!     -----------------------------------------------------------------

!*         3.     RELEASE THE FIRST (TEST) UPDRAFT TO GET PBL HEIGHTS
  
  
  !* set updraft index to 1
  JD = 1   

  DO JL=KIDIA,KFDIA
    
 
    PFRACB(JL,JD) = ZFRACTEST
 
 
    !* 3.1    DETERMINE STABILITY OF BL USING THE SURFACE BUOYANCY FLUX
    !*
    ZKHVFL(JL)  = ( 1.0_JPRB + RETV *  ZQTM1(JL,KLEV) ) * PKHFL(JL) + &
                & ( RETV * ZSLGM1(JL,KLEV) / RCPD )     * PKQFL(JL) 


    IF ( ZKHVFL(JL) >= 0.0_JPRB ) THEN
      
      ! stable BL (no updrafts expected/needed)
      KPBLTYPE(JL)  = 0

    ELSE

      LLDONE(JL,JD) = .FALSE.  !confirm launch
     
     
      !* 3.2    SURFACE LAYER SCALING
      !*
      ZUSTAR(JL)  = MAX( SQRT(PKMFL(JL)), ZREPUST )               !u* (repust=10e-4)
      ZWSTAR(JL)  = (- ZKHVFL(JL) * RG / PTM1(JL,KLEV) * 1000._JPRB ) &   !zi=1000m
                       & ** ( 1._JPRB/3._JPRB) 
      ZWSIGMA(JL)      = 1.2_JPRB &
       & * ( ZUSTAR(JL)**3 &
       & - 1.5_JPRB * RKAP * ZKHVFL(JL) * PGEOH(JL,KLEV-1) / PTM1(JL,KLEV-1) &
       & ) ** ( 1.0_JPRB/3._JPRB )                         ! Kolmogorov 1/3-power
      ZUSIGMA(JL)      = 2.29_JPRB  &
       & * ( ZUSTAR(JL)**3 &
       & + 0.5_JPRB / 12.0_JPRB * RKAP * ZWSTAR(JL)**3 &
       &   ) ** ( 1._JPRB/3._JPRB)
      
      !  scaling factors between phi* and initial updraft phi excesses
      !ZLOBUKHOV = -ZUSTAR(JL)**3 * PTM1(JL,KLEV) / (RG * RKAP * ZKHVFL(JL))
      !ZPDFFACPHI(JL) = ZRG*PGEOH(JL,KLEV-1)/ZLOBUKHOV
      
      
      !* 3.3    INITIALIZE TEST UPDRAFT
      !*
      
      !get the constant associated with the top ZFRACTEST % of the PDF
      ZFACEXC(JL,1) = ZFACTESTEXC
      
      !calculate the initial excess values
      ZWU2H(JL,KLEV-1,JD) = ( ZPDFFACW(JL) * ZFACEXC(JL,1) * ZWSIGMA(JL) )**2         
      ZTEXC            = - ZPDFFACPHI(JL) * ZFACEXC(JL,1) * PKHFL(JL) / ZWSIGMA(JL) 
      ZQEXC            = - ZPDFFACPHI(JL) * ZFACEXC(JL,1) * PKQFL(JL) / ZWSIGMA(JL) 
      ZTEXC            = MAX(ZTEXC, 0.0_JPRB)
      ZQEXC            = MAX(ZQEXC, 0.0_JPRB)
      PQTUH(JL,KLEV-1,JD) = ZQTENH(JL,KLEV-1) + ZQEXC
      ZQCUH(JL,KLEV-1,JD) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
      ZQUH (JL,KLEV-1,JD) = PQTUH(JL,KLEV-1,JD)  - ZQCUH(JL,KLEV-1,JD)
      PSLGUH(JL,KLEV-1,JD)= ZSLGENH(JL,KLEV-1) + RCPD * ZTEXC
      ZTUH (JL,KLEV-1,JD) = ( PSLGUH (JL,KLEV-1,JD) - PGEOH(JL,KLEV-1) &
                       & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1) &
                       & ) / RCPD

! ... u & v: (wind speed assumed to be negatively correlated with T and q excesses)
      ZSPEEDENV        = SQRT( ZUENH(JL,KLEV-1)**2 + ZVENH(JL,KLEV-1)**2 )
      ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC(JL,1) * ZUSIGMA(JL), 0._JPRB )
!     ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC(JL,1) * ZUSTAR(JL)**2/ ZWSIGMA(JL) , 0._JPRB )

!      PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
!      PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
      PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1)  !mean wind at this half level
      PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1)
!      PUUH(JL,KLEV-1,JD)= PUM1(JL,KLEV)     !10m wind instead of 20m
!      PVUH(JL,KLEV-1,JD)= PVM1(JL,KLEV)
!      PUUH(JL,KLEV-1,JD)= 0._JPRB           !using this is really bad for the scores!
!      PVUH(JL,KLEV-1,JD)= 0._JPRB
!      PUUH(JL,KLEV-1,JD)= PUM1(JL,KLEV-1)   !full level wind
!      PVUH(JL,KLEV-1,JD)= PVM1(JL,KLEV-1)
      
    ENDIF

  ENDDO !JL


  !* 3.4   RELEASE THE TEST UPDRAFT #1
  !*          - USED TO MAKE A FIRST GUESS OF THE HEIGHTS OF CLOUD BASE & INVERSION,
  !*            AND TO DETERMINE PBL TYPE.
  !*
  CALL VDFPARCEL (KIDIA   , KFDIA   , KLON    , KLEV    , KDRAFT  , &
                & PGEOH   , PGEOM1  , PAPHM1  , &
		& PUM1    , PVM1    , ZQTM1   , ZSLGM1  , ZTVEN   , &
		& PUUH    , PVUH    , PSLGUH  , PQTUH   , ZWU2H   , ZQCUH  , ZBUOF , & 
		& ZQUH    , ZTUH    , ZEPS    , ZFACEXC , &
		& PZPLCL  , KPLCL   , PZPTOP  , KPTOP   , KPLZB   , &
		& JD      , ZUPGENL , ZUPGENN , &
		& ZTAUEPS , ZW2THRESH, LLDONE , KPBLTYPE)  



!     -----------------------------------------------------------------

!*         4.     CLASSIFICATION OF THE CONVECTIVE PBL
!                 ------------------------------------


  !* 4.1    CLASSIFY THE CONVECTIVE PBL
  !*
  DO JL=KIDIA,KFDIA
    IF ( KPBLTYPE(JL)/=0 ) THEN
 
      IF ( PZPLCL(JL,1) > PZPTOP(JL,1) .OR. KPLCL(JL,1) == 0 ) THEN
      
        !dry convective PBL
        KPBLTYPE(JL)  = 1                   !dry convective PBL
        ZDZCLOUD(JL)  = 0.0_JPRB            !cloud thickness

      ELSE

        !moist convective PBL
        ZDZCLOUD(JL)  = PZPTOP(JL,1) - PZPLCL(JL,1) !cloud thickness
            
        IF (ZDZCLOUD(JL)>ZCLDDEPTHDP .AND. .NOT.LLMASSCAP ) THEN
        
	  !deep convection
	  KPBLTYPE(JL) = 4
          
	ELSE
        
          IF (LLSTCU) THEN
            KPBLTYPE(JL) = 2   !set the type to stratocumulus for the moment
          ELSE
            KPBLTYPE(JL) = 3   !RN run without Klein-Hartmann criterion!
          ENDIF  
          
        ENDIF
	
      ENDIF

    ENDIF !KPBLTYPE /=0
  ENDDO !JL
  
  
  
  !* 4.2    CHECK THE STRATOCUMULUS/SHALLOW CUMULUS CRITERION (TRIGGER FUNCTION)
  !*        IF SHALLOW CUMULUS IS DIAGNOSED, KPBLTYPE WILL BE SET TO 3
  !*
  CALL VDFSTCUCRIT ( KIDIA , KFDIA  , KLON  , KLEV , KDRAFT , &
		  &    PTM1  , ZSLGM1 , ZQTM1 , PAPM1 , &
		  &    ZSTABTHRESH, ZCLDDEPTH, ZBIRTHRESH, ZDZCLOUD, &
		  &    KPTOP , KPBLTYPE, LDNODECP, &
                  &    ZSTABILITY )
  
  !-- formulate decoupling constraint on moist updraft area fraction --                    
  DO JL=KIDIA,KFDIA
    
    ZCOUPLING(JL) = MAX( 0._JPRB, (ZSTABILITY(JL) - 20._JPRB) / 2._JPRB )
    ZCOUPLING(JL) = MIN( ZCOUPLING(JL), 1._JPRB )
    
    ZZFUNC4(JL) = ZFRACMAX * ZCOUPLING(JL)
    
    !ZZFUNC4(JL) = ZFRACMAX
    !ZZFUNC4(JL) = 0._JPRB
  
!    IF (LLDIAG) THEN
!      PEXTR2(JL,10) = ZSTABILITY(JL)
!    ENDIF 

  ENDDO

  
  
!     -----------------------------------------------------------------

!*         5.     CLOSURE FOR ORGANIZED UPDRAFTS (JD=2,3)
!                 ---------------------------------------


  !* 5.1    MIXED LAYER SCALINGS
  !*

  DO JL=KIDIA,KFDIA
    
    IF ( KPBLTYPE(JL)/=0 ) THEN    !don't do this for stable PBL
      
      !--- Mixed layer scaling depth ---
      SELECT CASE (KPBLTYPE(JL))
    
        CASE(1)
          !Dry convective PBL - Inversion height
          ZZI(JL)   = PZPTOP(JL,1)
          IZI(JL,1) = KPTOP(JL,1)
		
        CASE(2)
          !Stratocumulus - Inversion height
          !CAUTION: During decoupling in the intermediate regime (e.g. ASTEX/ATEX) the
	  !   relevant ML scaling height changes from PBL inversion to level of minimum
          !   buoyancy flux. In the current setup this is not modelled yet!
          ZZI(JL)   = PZPTOP(JL,1)
          IZI(JL,1) = KPTOP(JL,1)
	  
        CASE(3)
          !Shallow cumulus - Level of minimum buoyancy flux
	  !Assume that the moist updraft LCL is very close to this level
          ZZI(JL)   = PZPLCL(JL,1)
          IZI(JL,1) = KPLCL(JL,1)
	  
        CASE(4)
          !Deep cumulus - Only do a dry parcel up to cloud base
          ZZI(JL)   = PZPLCL(JL,1)
          IZI(JL,1) = KPLCL(JL,1)
		
      END SELECT 

      !--- Mixed layer convective velocity scale ---
      ZWSTAR(JL) = ( -ZKHVFL(JL) * RG * ZZI(JL) / ZTHVEN(JL,KLEV)  ) ** (1._JPRB/3._JPRB)
      
    ENDIF
    
  ENDDO  

    
    
  !*  5.2    RI NUMBER OF CUMULUS INVERSION
  !* 
  
  !-- Test-updraft cloudy CAPE --
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA
      IF ( ZQCUH(JL,JK,1)>0._JPRB .AND. ZBUOF(JL,JK,1)>0._JPRB .AND. JK<=KPLCL(JL,1) ) THEN  
        ZDZ = ZRG*( PGEOH(JL,JK-1) - PGEOH(JL,JK) )
        ZCAPE1(JL) = ZCAPE1(JL) + ZDZ * ZBUOF(JL,JK,1) 
      ENDIF
    ENDDO
  ENDDO
  
  DO JL=KIDIA,KFDIA
    IF ( KPBLTYPE(JL)==2 .OR. KPBLTYPE(JL)==3 ) THEN
    
      !-- Interpolate LZB height --
      ZZPLZB(JL,1) = PZPTOP(JL,1)
      IF (KPLZB(JL,1)>2) THEN
        ZDZ = (PGEOH(JL,KPLZB(JL,1)-1) - PGEOH(JL,KPLZB(JL,1)))*ZRG
        ZDB = ZBUOF(JL,KPLZB(JL,1),1)  - ZBUOF(JL,KPLZB(JL,1)-1,1)
        IF (ZDB>0._JPRB) THEN
          ZZPLZB(JL,1) = PGEOH(JL,KPLZB(JL,1)) * ZRG + &
                       & ZDZ * ZBUOF(JL,KPLZB(JL,1),1) / ZDB 
          ZZPLZB(JL,1) = MIN( ZZPLZB(JL,1), PGEOH(JL,KPLZB(JL,1)-1)*ZRG )
        ENDIF               
      ENDIF
      
      !-- Cloud layer average test-updraft positive buoyancy --
      ZBUOYCU = 0._JPRB
      IF ( ZZPLZB(JL,1)-PZPLCL(JL,1)>0._JPRB ) THEN  
	ZBUOYCU = ZCAPE1(JL) / ( ZZPLZB(JL,1) - PZPLCL(JL,1) )
      ENDIF  
  
      !-- Inversion theta_v jump --
      !JK = KPLZB(JL,1)    !use level of zero buoyancy (LZB) of test-updraft
      JK = KPTOP(JL,1)    !use top level of test-updraft
      IF (JK>2) THEN
        ZDTHVCUTOP = MAX( ZTHVEN(JL,JK-1)-ZTHVEN(JL,JK), ZTHVEN(JL,JK)-ZTHVEN(JL,JK+1) )
      ENDIF  
      
      !-- Cumulus Ri number - used again in VDFEXCU --
      IF ( ZDTHVCUTOP > 0._JPRB ) THEN 
        PRICUI(JL) = ZBUOYCU * ZRG * ZTHVEN(JL,KLEV) / ZDTHVCUTOP   
      ENDIF  
      
      !-- RN testing: no top entrainment for stcu (yikes) --
      !IF ( ZSTABILITY(JL) > ZSTABTHRESH ) THEN
      !  PRICUI(JL) = 0._JPRB
      !ENDIF
      !PRICUI(JL) = ( 1._JPRB - ZCOUPLING(JL) ) * PRICUI(JL)
      
!      IF (LLDIAG) THEN
!        PEXTR2(JL,14) = ZBUOYCU
!        PEXTR2(JL,15) = ZDTHVCUTOP
!        PEXTR2(JL,16) = RG * ZDTHVCUTOP / ZTHVEN(JL,KLEV)
!        PEXTR2(JL,17) = MIN(PRICUI(JL),1._JPRB)
!      ENDIF  

    ENDIF
  ENDDO
  
  

  !* 5.3    CLOSURE OF UPDRAFT AREA FRACTIONS (JD=2,3)
  !*
    
  DO JL=KIDIA,KFDIA
    
    IF ( KPBLTYPE(JL)/=0 ) THEN    !don't do this for stable PBL


      IF (KPBLTYPE(JL)>1) THEN
      
        !--- Transition layer depth, Scale I: Dry entrainment layer depth ---
        !       using thv gradient averaged over a number of layers above h
        !
        ZDTHVDZ  = MAX( 0.01_JPRB,ZTHVEN(JL,IZI(JL,1)-2)-ZTHVEN(JL,IZI(JL,1)) ) * RG / &
               & ( PGEOM1(JL,IZI(JL,1)-2) - PGEOM1(JL,IZI(JL,1)) )

        ZW2H      = 0.5_JPRB * ZWSTAR(JL)**2.
        ZDHRI(JL) = (   ZW2H * ZRG * ZTHVEN(JL,KLEV) / (0.5_JPRB * ZDTHVDZ)  )**0.5
        PDTHV(JL) = ZDHRI(JL) * ZDTHVDZ   !used again in VDFEXCU
      
        ZZFUNC3(JL) = MAX( 0._JPRB, 0.2_JPRB * ZDHRI(JL) / PZPLCL(JL,1) )   !cy32r3
      
!        IF (LLDIAG) THEN
!          PEXTR2(JL,50) = PZPLCL(JL,1)
!          PEXTR2(JL,51) = 100._JPRB * ZZFUNC3(JL)
!          PEXTR2(JL,52) = ZDHRI(JL)
!          PEXTR2(JL,48) = ZDTHVDZ
!          PEXTR2(JL,49) = PDTHV(JL)
!          PEXTR2(JL,55) = KPLCL(JL,1)-2 - KPTOP(JL,1)
!        ENDIF      

      ENDIF
      
      
      !--- Calculation of moist updraft area fraction ---
      SELECT CASE (KPBLTYPE(JL))
    
        CASE(1)
	
          !Dry convective PBL
	  !Set area fraction of moist group to zero
	  ZFRACB(JL,3) = 0._JPRB

        CASE(2)
	
          !Stratocumulus
	  !Set area fraction of moist group to ZFRACMAX
	  ZFRACB(JL,3) = ZFRACMAX

        CASE(3)
	  
          !Shallow cumulus
          !Flexible updraft area fractions

	  !-- Transition layer depth, Scale II: Cumulus condensation depth-scale --
          !
          ZDHCL = MIN(200._JPRB, 0.1_JPRB * ZDZCLOUD(JL))
          ZDHCL = MAX(ZDHCL,0._JPRB)
          
          !Use M/w* ~ Dh / h  (Neggers et al., QJ, 2007)
          IF (PZPLCL(JL,1).GT.0._JPRB) THEN
            ZZFUNC =  0.15_JPRB * ( ZDHCL/PZPLCL(JL,1) )     !CY32R3
            !ZZFUNC =  0.2_JPRB * ( ZDHCL/PZPLCL(JL,1) )   
          ELSE
            ZZFUNC = 0._JPRB
          ENDIF
          
!          IF (LLDIAG) THEN
!            PEXTR2(JL,53) = ZDHCL
!            PEXTR2(JL,44) = 100._JPRB*ZZFUNC
!          ENDIF  
          
          !-- Choose the minimum of scales I and II --
          ZFRACB(JL,3)  = MIN( ZFRACMAX, ZZFUNC, ZZFUNC3(JL) )
          ZFRACB(JL,3)  = MAX( ZFRACB(JL,3), ZZFUNC4(JL) )  !superimpose decoupling criterion
              
          !-- 1st guess for the LCL mass flux --
          ZMFLXB(JL,3)  = ZWSTAR(JL) * ZFRACB(JL,3) * ZRHOH(JL,KPLCL(JL,1))
          
          !-- Switch KPBLTYPE to dry convective if moist updraft is not launched --
          IF (ZFRACB(JL,3).EQ.0._JPRB) THEN
            KPBLTYPE(JL)=1
          ENDIF  
	  
	  
        CASE(4)
	
          !Deep cumulus
	  !Set area fraction of moist group to zero (only allow updraft transport in dry mixed layer)
	  ZFRACB(JL,3) = 0._JPRB


      END SELECT !KPBLTYPE
      
      
      !--- Dry updraft area fraction (JD=2) ---
      ZFRACB(JL,2) = MAX( 0._JPRB, ZFRACMAX - ZFRACB(JL,3) )
      
      PFRACB(JL,2) = ZFRACB(JL,2)
      PFRACB(JL,3) = ZFRACB(JL,3)


    ENDIF !KPBLTYPE /=0
    
  ENDDO !JL


      
!     -----------------------------------------------------------------

!*         6.     CALCULATE VERTICAL PROFILES OF ALL UPDRAFTS (JD=2,3)
!                 ----------------------------------------------------


  !*       6.1    CALCULATE THE SCALING FACTORS OF THE UPDRAFT EXCESS WITH THE SURFACE JOINT PDFS
  !*
  DO JD = 2,KDRAFT
    DO JL=KIDIA,KFDIA
      
      IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB ) THEN
        
        !-- Get the PDF scaling factor --
	SELECT CASE (JD)
  	  
	  CASE(2)
	    !lower part of top ZFRACMAX %
	    ZDUMFRAC = ZFRACMAX - ZFRACB(JL,2)
            CALL VDFPDFTABLE(ZDUMFRAC , ZFACEXC(JL,2), ZDUMR, ZDUMR, 0)
	    ZFACEXC(JL,2) = ( ZFRACMAX * ZFACMAXEXC - ZDUMFRAC * ZFACEXC(JL,2) ) / ZFRACB(JL,2)
		    
	  CASE(3)
	    !upper part of top ZFRACMAX %
	    ZDUMFRAC = ZFRACB(JL,JD)
            CALL VDFPDFTABLE(ZDUMFRAC , ZFACEXC(JL,3), ZDUMR, ZDUMR, 0)
	    
	END SELECT
	
      ENDIF !KPBLTYPE & ZFRACB

    ENDDO !JL
  ENDDO !JD
    
          
    
  !*       6.2    VERTICAL INTEGRATION OF DRY & MOIST UPDRAFT BUDGETS (JD=2,3)
  !*
  DO JD = 2,KDRAFT
  
    !-- Initialize updraft --
    DO JL=KIDIA,KFDIA
      
      IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB ) THEN
        
        LLDONE(JL,JD) = .FALSE. !confirm launch
          
        ZWU2H(JL,KLEV-1,JD) = ( ZPDFFACW(JL) * ZFACEXC(JL,JD) * ZWSIGMA(JL) )**2 
        ZTEXC            = - ZPDFFACPHI(JL) * ZFACEXC(JL,JD) * PKHFL(JL) / ZWSIGMA(JL) 
        ZQEXC            = - ZPDFFACPHI(JL) * ZFACEXC(JL,JD) * PKQFL(JL) / ZWSIGMA(JL) 
        ZTEXC            = MAX(ZTEXC, 0.0_JPRB)
        ZQEXC            = MAX(ZQEXC, 0.0_JPRB)
        PQTUH(JL,KLEV-1,JD) = ZQTENH(JL,KLEV-1) + ZQEXC
        ZQCUH(JL,KLEV-1,JD) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
        ZQUH (JL,KLEV-1,JD) = PQTUH(JL,KLEV-1,JD)  - ZQCUH(JL,KLEV-1,JD)
        PSLGUH(JL,KLEV-1,JD)= ZSLGENH(JL,KLEV-1) + RCPD * ZTEXC
        ZTUH (JL,KLEV-1,JD) = ( PSLGUH (JL,KLEV-1,JD) - PGEOH(JL,KLEV-1) &
                       & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1) &
                       & ) / RCPD

!   ... u & v: (wind speed assumed to be negatively correlated with T and q excesses)
        ZSPEEDENV        = SQRT( ZUENH(JL,KLEV-1)**2 + ZVENH(JL,KLEV-1)**2 )
        ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC(JL,JD) * ZUSIGMA(JL), 0._JPRB )
!       ZSPEEDUP         = MAX( ZSPEEDENV - ZFACEXC(JL,JD) * ZUSTAR(JL)**2/ ZWSIGMA(JL) , 0._JPRB )

!        PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
!        PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1) * ZSPEEDUP/ZSPEEDENV
        PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1)  !mean wind at this half level
        PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1)
!        PUUH(JL,KLEV-1,JD)= PUM1(JL,KLEV)     !10m wind instead of 20m
!        PVUH(JL,KLEV-1,JD)= PVM1(JL,KLEV)
!        PUUH(JL,KLEV-1,JD)= 0._JPRB           !using this is really bad for the scores!
!        PVUH(JL,KLEV-1,JD)= 0._JPRB
!        PUUH(JL,KLEV-1,JD)= PUM1(JL,KLEV-1)   !full level wind
!        PVUH(JL,KLEV-1,JD)= PVM1(JL,KLEV-1)
 
      ENDIF !KPBLTYPE & ZFRACB
      
    ENDDO !JL
    
    
    
    
    !-- Release the updraft --
    CALL VDFPARCEL (KIDIA   , KFDIA   , KLON    , KLEV    , KDRAFT  , &
                  & PGEOH   , PGEOM1  , PAPHM1  , &
		  & PUM1    , PVM1    , ZQTM1   , ZSLGM1  , ZTVEN   , &
		  & PUUH    , PVUH    , PSLGUH  , PQTUH   , ZWU2H   , ZQCUH  , ZBUOF , & 
		  & ZQUH    , ZTUH    , ZEPS    , ZFACEXC , &
		  & PZPLCL  , KPLCL   , PZPTOP  , KPTOP   , KPLZB   , &
		  & JD      , ZUPGENL , ZUPGENN , &
		  & ZTAUEPS , ZW2THRESH, LLDONE , KPBLTYPE)  
    
  ENDDO !JD 



  !*        6.3 SOME PROTECTIONAL MEASURES AGAINST UPDRAFT #3 FAILING (ONLY JUST) TO REACH CONDENSATION.
  !*
  DO JL=KIDIA,KFDIA

    IF ( KPBLTYPE(JL)==2 .OR. KPBLTYPE(JL)==3) THEN
      IF ( ZFRACB(JL,3)>0._JPRB .AND. PZPLCL(JL,3)<0._JPRB) THEN
      
        !set cloud base height to updraft #3 top (cloud layer exists but has zero depth)
        KPLCL(JL,3)  = KPTOP(JL,3)
        PZPLCL(JL,3) = PZPTOP(JL,3)
        
      ENDIF
      
    ENDIF
     
  ENDDO !JL 


    
  !*        6.4  LIMITER FOR CLOUDY DEPTH OF MOIST UPDRAFT  *TESTING*
  !* 
  IF (LLMASSCAP) THEN
  
    DO JL=KIDIA,KFDIA
      LLMCIND(JL) = .FALSE.
      IF (KPBLTYPE(JL)==3) THEN
        IF ( PZPLCL(JL,3)>0._JPRB .AND. (PZPTOP(JL,3) - PZPLCL(JL,3))>ZMASSCAPDEPTH ) THEN
          PZPTOP(JL,3) = PZPLCL(JL,3) + ZMASSCAPDEPTH
          LLMCIND(JL) = .TRUE.
        ENDIF
      ENDIF
    ENDDO
    
    DO JK=KLEV-1,1,-1
      DO JL=KIDIA,KFDIA
      
        IF ( LLMCIND(JL) ) THEN
          IF ( PGEOH(JL,JK+1)*ZRG<=PZPTOP(JL,3) .AND. PGEOH(JL,JK)*ZRG>PZPTOP(JL,3) ) THEN
            KPTOP(JL,3) = JK+1
          ENDIF
        ENDIF
        
      ENDDO
    ENDDO
    
  ENDIF


  !*        6.5  UPDRAFT PRECIPITATION FLUXES (RAIN AND SNOW)
  !* 
  DO JD = 3,KDRAFT  !moist updrafts only
    
    DO JK=2,KLEV
      DO JL=KIDIA,KFDIA
      
        ZDZRHO = ZRG * ( PAPHM1(JL,JK)-PAPHM1(JL,JK-1) ) 
        
        !-- Add precip generation to flux [kg /m2 /s: tendency * layer depth * air density] --
        ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK-1,JD) + ZUPGENL(JL,JK,JD) * ZDZRHO
        ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK-1,JD) + ZUPGENN(JL,JK,JD) * ZDZRHO
        
        !-- Do some melting at freezing level (snow->rain) --
        IF (ZUPFLXN(JL,JK,JD)>0._JPRB .AND. PTM1(JL,JK) > RTT) THEN
          ZUPMELT = (1.0_JPRB+0.5_JPRB*(PTM1(JL,JK)-RTT)) * &
                  & (PTM1(JL,JK)-RTT) * RCPD/(RLMLT*RTAUMEL) * ZDZRHO
          ZUPMELT = MIN(ZUPFLXN(JL,JK,JD),ZUPMELT)
          ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK,JD) + ZUPMELT
          ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK,JD) - ZUPMELT
        ENDIF  
        
        !-- Saturation deficit of mean state T --
        ZQSATM = FOEEWM(PTM1(JL,JK))/PAPM1(JL,JK)
        ZQSATM = MIN(0.5_JPRB,ZQSATM)
        ZQSATM = ZQSATM/(1.0_JPRB-RETV*ZQSATM)
        ZSATDEF = MAX( 0._JPRB, ZQSATM-PQM1(JL,JK) )

        ZPFLXTOT = ZUPFLXL(JL,JK,JD) + ZUPFLXN(JL,JK,JD)
        IF (ZPFLXTOT>0._JPRB) THEN
         
          !-- Precip evaporation tendency [kg/kg /s] (Kessler 1969, Tiedtke 1993) --
          ZPEVAPUP = ZAPRECEVAP * ZSATDEF * ( &
              & ( ZPFLXTOT / 0.00509_JPRB ) * &
              & ( PAPM1(JL,JK)/PAPHM1(JL,KLEV) )**0.5_JPRB &
              & )**0.5777_JPRB

          !RN testing: instantaneous evaporation
          !ZPEVAPUP = ( ZUPFLXL(JL,JK,JD) + ZUPFLXN(JL,JK,JD) ) / ZDZRHO

          !RN testing: no evap in subcloud layer
          !IF (JK>=KPLCL(JL,3)) then
          !  ZPEVAPUP = 0._JPRB
          !ENDIF
        
          !-- Back-partition evaporation and substract from fluxes --
          ZFAC = ZUPFLXL(JL,JK,JD) / ZPFLXTOT
          ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK,JD) - ZPEVAPUP * ZDZRHO * ZFAC
          ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK,JD) - ZPEVAPUP * ZDZRHO * (1._JPRB - ZFAC)
          ZUPFLXL(JL,JK,JD) = MAX(0._JPRB,ZUPFLXL(JL,JK,JD))
          ZUPFLXN(JL,JK,JD) = MAX(0._JPRB,ZUPFLXN(JL,JK,JD))
        ENDIF
        
      ENDDO
    ENDDO
    
    !Add contribution to total flux - weight by updraft area fraction
    DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
        PFPLVL(JL,JK) = PFPLVL(JL,JK) + ZFRACB(JL,JD) * ZUPFLXL(JL,JK,JD)
        PFPLVN(JL,JK) = PFPLVN(JL,JK) + ZFRACB(JL,JD) * ZUPFLXN(JL,JK,JD)
      ENDDO
    ENDDO
    
  ENDDO !JD
        
        

!     -----------------------------------------------------------------

!*         7.     CONSTRUCT MASS FLUX PROFILES (JD=2,3)
!                 -------------------------------------


  !*         7.1  DETERMINE THE MIXED LAYER SCALING HEIGHT FOR JD=2,3
  !*
  DO JL=KIDIA,KFDIA
    
    SELECT CASE (KPBLTYPE(JL))
    
      CASE(1)
        !Dry convective PBL - no moist parcel
        IZI(JL,2) = KPTOP(JL,2)      !half level below level of zero kinetic energy
                
      CASE(2)
        !Stratocumulus - no dry parcel
        IZI(JL,3) = KPLCL(JL,3)+1      !half level below lcl
	  
      CASE(3)
        !Shallow cumulus - both dry and moist
        IZI(JL,3) = KPLCL(JL,3)+1    !half level below lcl
        IZI(JL,2) = KPTOP(JL,2)      !half level below level of zero kinetic energy

      CASE(4)
        !Deep cumulus - no moist parcel
        IZI(JL,2) = KPTOP(JL,2)  

    END SELECT 

  ENDDO !JL 



  !*         7.2  CONSTRUCT MIXED LAYER MASS FLUXES
  !*               - USE CONSTANT AREA FRACTION, AND MULTIPLY BY PARCEL W
  !*
  DO JD = 2,KDRAFT

    DO JK=KLEV-1,1,-1

      DO JL=KIDIA,KFDIA

        IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB) THEN
	 IF (JK>=IZI(JL,JD) ) THEN
	  
          ZWUH = MAX( ZWU2H(JL,JK,JD),0._JPRB )
          ZWUH = ZWUH**0.5_JPRB
          PMFLX(JL,JK,JD)  = ZFRACB(JL,JD) * ZWUH * ZRHOH(JL,JK)
          
	  IF (ZWU2H(JL,JK,JD)>0._JPRB) THEN
	    ZFRAC(JL,JK,JD) = ZFRACB(JL,JD)		   
	  ELSE
	    ZFRAC(JL,JK,JD) = 0._JPRB
	  ENDIF
          
	 ENDIF
	ENDIF
	
      ENDDO !JL
    
    ENDDO !JK
    
  ENDDO !JD
  
    
  
  !*         7.3    BUOYANCY SORTING ON MOIST UPDRAFT IN CLOUD LAYER
  !*
  
  CALL VDFBUOYSORT( KIDIA     , KFDIA   , KLON    , KLEV   , KDRAFT , &
                  & PAPM1     , PGEOM1  , PGEOH   , &
                  & ZQTM1     , ZSLGM1  , &
                  & PFRACB    , KPLCL   , KPTOP   , KPLZB  , &
                  & PQTUH     , PSLGUH  , ZWU2H   , PUUH   , PVUH   , &
                ! DIAGNOSTIC OUTPUT
                  & PEXTR2    , KFLDX2  , PEXTRA  , KLEVX  , KFLDX  , &
                !              
                  & ZABULK    , ZWBULK  , ZQTBULK , ZSLGBULK , ZUBULK , ZVBULK )
  


  !*         7.4    CONSTRUCT CLOUDY MASS FLUX PROFILE (JD=3 ONLY)
  !*
  DO JK=KLEV-2,1,-1

    DO JL=KIDIA,KFDIA
   
      IF ( KPBLTYPE(JL)/=0 .AND. KPBLTYPE(JL)/=4 .AND. ZFRACB(JL,3)>0._JPRB ) THEN
        
        IF (JK>=KPTOP(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN

        ZWUH = MAX(0._JPRB,ZWU2H(JL,JK,3))**0.5_JPRB
        
        IF( JK==KPLCL(JL,3) ) THEN 

          !  Special treatment for cloud base level

          ! -- moist updraft --
!          PMFLX(JL,JK,3) = ZFRACB(JL,3)*ZWUH*ZRHOH(JL,JK) 
          PMFLX(JL,JK,3) = MIN( ZFRACB(JL,3)*ZWUH*ZRHOH(JL,JK), 2._JPRB*PMFLX(JL,JK+1,3))   !cap acceleration through cloud base (can cause instability)
          ZFRAC(JL,JK,3) = ZFRACB(JL,3)
          
          ! -- dry updraft --
          !    tie dry updraft flux to its value at layer below, to prevent too sharp gradients
          !PMFLX(JL,JK,2) = 0.3_JPRB*PMFLX(JL,JK+1,2)
          !ZFRAC(JL,JK,2) = 0.3_JPRB*ZFRAC(JL,JK+1,2)
          PMFLX(JL,JK,2) = 1.0_JPRB*PMFLX(JL,JK+1,2)
          ZFRAC(JL,JK,2) = 1.0_JPRB*ZFRAC(JL,JK+1,2)
          !PMFLX(JL,JK,2) = 0.0_JPRB
          !ZFRAC(JL,JK,2) = 0.0_JPRB
          
          PQTUH(JL,JK,2)  = 0.3_JPRB*( PQTUH(JL,JK+1,2)  - ZQTENH(JL,JK+1)  ) + ZQTENH(JL,JK)
          PSLGUH(JL,JK,2) = 0.3_JPRB*( PSLGUH(JL,JK+1,2) - ZSLGENH(JL,JK+1) ) + ZSLGENH(JL,JK)
          PUUH(JL,JK,2)   = 0.3_JPRB*( PUUH(JL,JK+1,2)   - ZUENH(JL,JK+1)   ) + ZUENH(JL,JK)
          PVUH(JL,JK,2)   = 0.3_JPRB*( PVUH(JL,JK+1,2)   - ZVENH(JL,JK+1)   ) + ZVENH(JL,JK)
          
	ELSEIF ( PMFLX(JL,JK+1,3) > 0._JPRB ) THEN
          
          IF ( JK>=KPTOP(JL,3) ) THEN
            
            !-- convective cloud layer --
            !   buoyancy sorting
            
            ZFRAC(JL,JK,3)  = ZABULK(JL,JK)
            ZWU2H(JL,JK,3)  = ZWBULK(JL,JK) ** 2._JPRB
            PMFLX(JL,JK,3)  = ZABULK(JL,JK) * ZWBULK(JL,JK) * ZRHOH(JL,JK)
            PQTUH(JL,JK,3)  = ZQTBULK(JL,JK)
            PSLGUH(JL,JK,3) = ZSLGBULK(JL,JK)
            PUUH(JL,JK,3)   = ZUBULK(JL,JK)
            PVUH(JL,JK,3)   = ZVBULK(JL,JK)
            
          ELSE 
            
            !-- inversion (layer between tops of test parcel and moist parcel) --
            !   prescribed linear decay of flux
            
            !ZFRAC(JL,JK,3)  = ZFRAC(JL,JK+1,3) * 0.25_JPRB
            !PMFLX(JL,JK,3)  = PMFLX(JL,JK+1,3) * 0.25_JPRB
            zfac = ( PZPTOP(JL,1) - PGEOH(JL,JK)*ZRG  ) / ( PZPTOP(JL,1) - PZPTOP(JL,3) )
            zfac = MAX( 0._JPRB, MIN(1._JPRB,zfac) )
            ZFRAC(JL,JK,3)  = ZFRAC(JL,KPTOP(JL,3),3) * zfac
            PMFLX(JL,JK,3)  = PMFLX(JL,KPTOP(JL,3),3) * zfac
            
            PQTUH(JL,JK,3)  = ( PQTUH(JL,JK+1,3)  - ZQTENH(JL,JK+1)  ) + ZQTENH(JL,JK)
            PSLGUH(JL,JK,3) = ( PSLGUH(JL,JK+1,3) - ZSLGENH(JL,JK+1) ) + ZSLGENH(JL,JK)
            PUUH(JL,JK,3)   = ( PUUH(JL,JK+1,3)   - ZUENH(JL,JK+1)   ) + ZUENH(JL,JK)
            PVUH(JL,JK,3)   = ( PVUH(JL,JK+1,3)   - ZVENH(JL,JK+1)   ) + ZVENH(JL,JK)
            
            !ZWU2H(JL,JK,3)  = ZWU2H(JL,KPTOP(JL,3),3)
            ZWU2H(JL,JK,3)  = ZWU2H(JL,KPTOP(JL,3),3) * (ZFAC**2._JPRB)
            !ZWU2H(JL,JK,3)  = ZWU2H(JL,JK,1)
            
          ENDIF
	  
          !make sure that updraft #2 does not do any flux here
          PMFLX(JL,JK,2)  = 0._JPRB
          ZFRAC(JL,JK,2)  = 0._JPRB
          PQTUH(JL,JK,2)  = 0._JPRB
          PSLGUH(JL,JK,2) = 0._JPRB
          PUUH(JL,JK,2)   = 0._JPRB
          PVUH(JL,JK,2)   = 0._JPRB
          
	ENDIF 
	
	
        !--- limit mass flux covering 50% area (M<rho*w,up*0.5) ---
        !    (detrainment is initiated if strong w,up slowdown)
        !inv: comment out
!inv        PMFLX(JL,JK,3) =   MIN( PMFLX(JL,JK,3) , 0.5_JPRB * ZWUH * ZRHOH(JL,JK) )   

	ENDIF 
        
      ENDIF
     
    ENDDO !JL
     
  ENDDO !JK

  

  !*        7.5  COMPUTE I) VARIANCE TRANSPORT FLUX AND II) BULK UPDRAFT DETRAINMENT, AS USED IN VDFMAIN
  !*
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA
!      IF (JK >= KPTOP(JL,1) ) THEN  
      IF (JK > KPTOP(JL,1) ) THEN  
      
        !-- calculate w'qt'qt' at half levels --
        PWQT2(JL,JK) =  PMFLX(JL,JK,3) * (PQTUH(JL,JK,3) - ZQTENH(JL,JK))**2._JPRB + &
                     &  PMFLX(JL,JK,2) * (PQTUH(JL,JK,2) - ZQTENH(JL,JK))**2._JPRB

      ENDIF 
    ENDDO !JL
  ENDDO !JK
  
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF ( KPBLTYPE(JL)/=0 .AND. KPBLTYPE(JL)/=4 .AND. ZFRACB(JL,3)>0._JPRB ) THEN
        
!      IF (JK >= KPTOP(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN  
      IF (JK > KPTOP(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN  
	
        ZDZ   =  (   PGEOH(JL,JK-1)   - PGEOH(JL,JK)   ) * ZRG 
        
	ZDMDZ = -(   PMFLX(JL,JK-1,3) - PMFLX(JL,JK,3) &
!	          & + PMFLX(JL,JK-1,2) - PMFLX(JL,JK,2) &  
		 & ) / ZDZ
                 
        IF (JK==KPTOP(JL,1)+1) THEN
	  ZDMDZ = PMFLX(JL,JK,3) / ZDZ
        ENDIF
        
        ZDMDZ = MAX( 0._JPRB, ZDMDZ )   
           
        PDETR(JL,JK) = ZDMDZ + PMFLX(JL,JK,3) / ( ZTAUEPS * MAX(0.0001_JPRB,ZWU2H(JL,JK,3))**0.5_JPRB )  
        
!        IF (LLDIAG) THEN
!          PEXTRA(JL,JK,98) = ZDMDZ
!          PEXTRA(JL,JK,99) = PDETR(JL,JK) - ZDMDZ
!        ENDIF  
	
      ENDIF
      	
      ENDIF	
    ENDDO !JL
  ENDDO !JK
          
  

  !*        7.6  ESTIMATE MASS FLUX NEEDED FOR CAPE REMOVAL (FRITSCH & CHAPPELL)
  !*
  
  IF ( LLCAPE ) THEN
  
  !-- calculate right-hand-side --
  DO JK=KLEV-1,2,-1
    DO JL=KIDIA,KFDIA
      IF (KPBLTYPE(JL) ==3 .AND. JK >= KPTOP(JL,3) .AND. JK <= KPLCL(JL,3) .AND. PMFLX(JL,KPLCL(JL,3),3)>0._JPRB ) THEN 
        ZRHS(JL) = ZRHS(JL) + &
      & ( PMFLX(JL,JK,3) / PMFLX(JL,KPLCL(JL,3),3) ) * &
      & ( RG/ZTHVEN(JL,KLEV) ) * &
      & ( ZTHVEN(JL,JK-1)-ZTHVEN(JL,JK) ) 
      ENDIF
    ENDDO !JL
  ENDDO !JK
  
  !-- calculate M_CAPE --
  DO JL=KIDIA,KFDIA
    ZFRACBCONG(JL) = 0._JPRB
    ZMSCALE(JL)    = 1._JPRB
    ZMASSCAPE(JL)  = 0._JPRB
    IF (ZRHS(JL)>0._JPRB) THEN
      ZMASSCAPE(JL)  = ZCAPE1(JL) / (ZTAUBM * ZRHS(JL))
      IF (KPLCL(JL,3).GT.0) THEN
        JK = KPLCL(JL,3)
        ZWUH = SQRT( MAX( 0._JPRB, ZWU2H(JL,JK,1) ) )
        IF (ZWUH > 0._JPRB) THEN
          ZFRACBCONG(JL) = MAX( 0._JPRB, ZMASSCAPE(JL)-PMFLX(JL,JK,3) ) / ZWUH
        ENDIF  
        IF (PMFLX(JL,JK,3) > 0._JPRB) THEN
          !ZMSCALE(JL) = ZMASSCAPE(JL)/PMFLX(JL,JK,3)                  !pure CAPE closure
          ZMSCALE(JL) = MAX(1._JPRB, ZMASSCAPE(JL)/PMFLX(JL,JK,3) )   !take the maximum of the two
        ENDIF  
      ENDIF  
    ENDIF  
  ENDDO !JL

  DO JK=KLEV-1,2,-1
    DO JL=KIDIA,KFDIA
      IF (LLCAPETEST) THEN
        ! Option I: assign mass to the test updraft
        PMFLX(JL,JK,1) = ZFRACBCONG(JL) * SQRT( MAX( 0._JPRB, ZWU2H(JL,JK,1) ) )  
      ELSE
        ! Option II: rescale moist mass flux
        PMFLX(JL,JK,3) = ZMSCALE(JL) * PMFLX(JL,JK,3)
      ENDIF  
    ENDDO !JL
  ENDDO !JK
  
  !-- diagnostics --
!  IF ( LLDIAG ) THEN
!    DO JL=KIDIA,KFDIA
!      PEXTR2(JL,5) = ZCAPE1(JL)
!      PEXTR2(JL,6) = ZFRACBCONG(JL)
!      PEXTR2(JL,7) = ZMASSCAPE(JL)
!      PEXTR2(JL,8) = PMFLX(JL,KPLCL(JL,3),3)
!    ENDDO
!    DO JK=1,KLEV
!      DO JL=KIDIA,KFDIA
!        PEXTRA(JL,JK,98) = SQRT( MAX( 0._JPRB, ZWU2H(JL,JK,1) ) )
!        PEXTRA(JL,JK,99) = PMFLX(JL,JK,1)
!      ENDDO
!    ENDDO
!  ENDIF

  ENDIF



  !*        7.6  MASS FLUX LIMIT ACCORDING TO CFL CRITERION
  !*
 
  ZCONS10 = 1.0_JPRB/(RG*PTMST)

!  DO JD = 1,KDRAFT
!    DO JL=KIDIA,KFDIA
!      ZMFS(JL,JD) = 1.0_JPRB    ! default reduction factor
!    ENDDO
!  ENDDO
  
  DO JD = 1,KDRAFT
    DO JK=1,KLEV-1
      DO JL=KIDIA,KFDIA
        IF ( JK >= KPTOP(JL,JD) .AND. KPTOP(JL,JD)>0) THEN
          ZMFMAX = (PAPM1(JL,JK+1)-PAPM1(JL,JK)) * ZCONS10

!          PEXTRA(JL,JK,6) = ZMFMAX

          IF (JD==1) THEN
            ZMFMAX = 1.0_JPRB * ZMFMAX
          ELSE
            ZMFMAX = 3.0_JPRB * ZMFMAX
          ENDIF  
          
!          ZMFMAX = 4.0_JPRB * ZMFMAX
!          ZMFMAX = 2.0_JPRB * ZMFMAX
!          ZMFMAX = 1.0_JPRB * ZMFMAX
!          ZMFMAX = 0.3_JPRB * ZMFMAX

!          PEXTRA(JL,JK,7) = ZMFMAX

          !Option I: preserve vertical structure
!          IF ( PMFLX(JL,JK,JD) > ZMFMAX ) THEN
!            ZMFS(JL,JD) = MIN(ZMFS(JL,JD),ZMFMAX/PMFLX(JL,JK,JD))
!          ENDIF
          
          !Option II: correct level by level
	  PMFLX(JL,JK,JD) = MIN( ZMFMAX, PMFLX(JL,JK,JD) )
	  
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  
!  DO JD = 1,KDRAFT
!    DO JK=1,KLEV
!      DO JL=KIDIA,KFDIA
!        PMFLX(JL,JK,JD) = PMFLX(JL,JK,JD)*ZMFS(JL,JD)
!      ENDDO
!    ENDDO
!  ENDDO
          

  !RN sensitivity test: no dry updraft
!  DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
!        PMFLX(JL,JK,2) = 0.
!    ENDDO
!  ENDDO
          


  !*        7.7  CLOUD-LAYER AVERAGE MOIST UPDRAFT MASS FLUX.
  !*               (USED IN VDFEXCU FOR ESTIMATING ENTRAINMENT-K AT CUMULUS PBL TOP)  
  !*
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA
      !-- cloudy mass flux: use moist updraft M (depth average) --
      IF ( ZQCUH(JL,JK,3)>0._JPRB .AND. PMFLX(JL,JK,3)>0._JPRB ) THEN  
        ZMCLD(JL) = ZMCLD(JL) + PMFLX(JL,JK,3) * &
                &     ZRG*( PGEOH(JL,JK-1) - PGEOH(JL,JK) )
      ENDIF
    ENDDO
  ENDDO
  
  DO JL=KIDIA,KFDIA
    IF (PZPTOP(JL,3)-PZPLCL(JL,3)>0._JPRB ) THEN  
!      PMCU(JL)    = ZMCLD(JL) / (PZPTOP(JL,3)-PZPLCL(JL,3))
      PMCU(JL) = ZMFLXB(JL,3)
    ENDIF
  ENDDO


    
  !*        7.8  SCALING (SEE VDFEXCU)
  !* 
  DO JD = 1,KDRAFT
    DO JK=1,KLEV-1
      DO JL=KIDIA,KFDIA
      
        ZMGEOM(JL,JK)  = PGEOM1(JL,JK)-PGEOM1(JL,JK+1)      
        ZCFNC1(JL,JK)  = RVDIFTS * PTMST * RG**2 * PAPHM1(JL,JK) &
                 & /( ZMGEOM(JL,JK) * RD * 0.5_JPRB &
                 & *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  )) &
                 & +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1))))  
        PMFLX(JL,JK,JD) = ZCFNC1(JL,JK) * ZMGEOM(JL,JK) * ZRG * PMFLX(JL,JK,JD)
        
      ENDDO
    ENDDO
  ENDDO



!     -----------------------------------------------------------------

!*         8.     SET SOME W-SCALES FOR USE IN VDFMAIN
!                 ------------------------------------

  DO JL=KIDIA,KFDIA
    PWUAVG(JL) = ZWSTAR(JL)
    !PWUAVG(JL) = 2.0_JPRB * ZWSTAR(JL)
  ENDDO

  !for use in buoyancy sorting scheme
  DO JD = 1,KDRAFT
    DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
        PWUH(JL,JK,JD) = SQRT( MAX( 0._JPRB, ZWU2H(JL,JK,JD) ) )
      ENDDO
    ENDDO
  ENDDO
  
  

!     -----------------------------------------------------------------

!*         9.     ADVECTIVE FLUX ADJUSTMENTS AT UPDRAFT TOP-LEVEL
!                 -----------------------------------------------
!
!                 Vertical mixing at the top-level is prescribed
!                 and controlled in VDFEXCU.
!

  DO JD = 2,KDRAFT
    DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
        
        !  Remove all mass flux at and above top layer of updraft
	ITOP = MAX( KPTOP(JL,JD), KPTOP(JL,3) )
	
        IF ( JK <= ITOP ) THEN
          
          LLWIPE=.TRUE.    
          
          ! Protect cumulus cloud top: this is treated later in VDFEXCU, dept. on LLRICU=T
          !   Note: do wipe in case of single layer moist convection
          IF ( KPBLTYPE(JL)==3 .AND. KPTOP(JL,3)<KPLCL(JL,3)+1 ) THEN   
            LLWIPE=.FALSE. 
          ENDIF
          
          IF (LLWIPE) THEN
            ZFRAC(JL,JK,JD)  = 0.0_JPRB
            PMFLX(JL,JK,JD)  = 0.0_JPRB
            PWQT2(JL,JK)     = 0.0_JPRB
          ENDIF
          
        ENDIF
        
      ENDDO
    ENDDO
  ENDDO
  
  
  
  !---------------- some output ------------------------
  IF (LLDIAG) THEN
  
  DO JL=KIDIA,KFDIA
  
    !  boundary layer classification
    PEXTR2(JL,30) = KPBLTYPE(JL)

    !  updraft heights
    PEXTR2(JL,31) = PZPTOP(JL,1)
    PEXTR2(JL,32) = PZPTOP(JL,2)
    PEXTR2(JL,33) = PZPTOP(JL,3)
    PEXTR2(JL,34) = PZPLCL(JL,1)
    PEXTR2(JL,35) = PZPLCL(JL,2)
    PEXTR2(JL,36) = PZPLCL(JL,3)
    
    !PEXTR2(JL,37) = PWUAVG(JL)
    !PEXTR2(JL,38) = SQRT ( ZWU2H(JL,KLEV-1,1) )
    
    !  updraft levels
    PEXTR2(JL,24) = KPLCL(JL,1)
    PEXTR2(JL,25) = KPLCL(JL,2)
    PEXTR2(JL,26) = KPLCL(JL,3)
    PEXTR2(JL,27) = KPTOP(JL,1)
    PEXTR2(JL,28) = KPTOP(JL,2)
    PEXTR2(JL,29) = KPTOP(JL,3)
    
    !  various scalings
    !PEXTR2(JL,45) = 100._JPRB*ZFRACB(JL,3)
      
    !  test updraft properties   
    ZFRAC(JL,0:KPTOP(JL,1)-1,1) = 0._JPRB
    ZFRAC(JL,KPTOP(JL,1):KLEV-1,1) = 0.0001_JPRB
    PEXTRA(JL,:,15) = 1000._JPRB * CEILING(ZFRAC(JL,:,1)) * ( PQTUH(JL,:,1)  - ZQTENH(JL,:) )
    PEXTRA(JL,:,16) = CEILING(ZFRAC(JL,:,1)) * ( PSLGUH(JL,:,1) - ZSLGENH(JL,:) )/RCPD 
    PEXTRA(JL,:,17) = ZWU2H(JL,:,1)
    PEXTRA(JL,:,18) = ZEPS(JL,:,1)
    PEXTRA(JL,:,19) = 1000._JPRB * ZQCUH(JL,:,1)
    PEXTRA(JL,1:KLEV,20) = ZBUOF(JL,:,1)
    
    !PEXTR2(JL,13) = ZCAPE1(JL)
    !PEXTR2(JL,18) = PMCU(JL)
    
    !  updraft precip generation
    PEXTRA(JL,:,21) = PFPLVL(JL,:)
    PEXTRA(JL,:,22) = PFPLVN(JL,:)

    !  updraft buoyancy
    PEXTRA(JL,1:KLEV,28) = ZBUOF(JL,:,2)
    PEXTRA(JL,1:KLEV,29) = ZBUOF(JL,:,3)

    !  updraft mass flux
    PEXTRA(JL,:,30) = PMFLX(JL,:,2) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )  
    PEXTRA(JL,:,31) = PMFLX(JL,:,3) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )
    
    !  updraft excesses
    PEXTRA(JL,:,32) = 1000._JPRB * CEILING(ZFRAC(JL,:,2)) * ( PQTUH(JL,:,2)  - ZQTENH(JL,:)  )
    PEXTRA(JL,:,33) = 1000._JPRB * CEILING(ZFRAC(JL,:,3)) * ( PQTUH(JL,:,3)  - ZQTENH(JL,:)  )
    PEXTRA(JL,:,34) = CEILING(ZFRAC(JL,:,2)) * ( PSLGUH(JL,:,2) - ZSLGENH(JL,:) )/RCPD 
    PEXTRA(JL,:,35) = CEILING(ZFRAC(JL,:,3)) * ( PSLGUH(JL,:,3) - ZSLGENH(JL,:) )/RCPD 
    PEXTRA(JL,:,36) = ZWU2H(JL,:,2)
    PEXTRA(JL,:,37) = ZWU2H(JL,:,3)
    PEXTRA(JL,:,26) = MAX(0._JPRB,ZWU2H(JL,:,2))**0.5 / MAX(0.01_JPRB,ZWSTAR(JL))
    PEXTRA(JL,:,27) = MAX(0._JPRB,ZWU2H(JL,:,1))**0.5 / MAX(0.01_JPRB,ZWSTAR(JL))

    !  updraft fractions
    PEXTRA(JL,:,38) = 100._JPRB * ZFRAC(JL,:,2)
    PEXTRA(JL,:,39) = 100._JPRB * ZFRAC(JL,:,3)
    
    !  updraft entrainment / detrainment
    PEXTRA(JL,:,40) = ZEPS(JL,:,2)
    PEXTRA(JL,:,41) = ZEPS(JL,:,3)
    
    !  updraft qt flux
    PEXTRA(JL,:,44) =  RLVTT * ( PMFLX(JL,:,2) * (PQTUH(JL,:,2) - ZQTENH(JL,:)) ) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )
    PEXTRA(JL,:,45) =  RLVTT * ( PMFLX(JL,:,3) * (PQTUH(JL,:,3) - ZQTENH(JL,:)) ) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )
    PEXTRA(JL,:,46) =  RLVTT * ( PMFLX(JL,:,2) * (PQTUH(JL,:,2) - ZQTENH(JL,:)) + &
                    &  PMFLX(JL,:,3) * (PQTUH(JL,:,3) - ZQTENH(JL,:)) ) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG ) 

    !  updraft condensate
    PEXTRA(JL,:,47) = 1000._JPRB * ZQCUH(JL,:,2)
    PEXTRA(JL,:,48) = 1000._JPRB * ZQCUH(JL,:,3)
            
    !  updraft thl flux
    PEXTRA(JL,:,51) =  ( PMFLX(JL,:,2) * (PSLGUH(JL,:,2) - ZSLGENH(JL,:)) ) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )
    PEXTRA(JL,:,52) =  ( PMFLX(JL,:,3) * (PSLGUH(JL,:,3) - ZSLGENH(JL,:)) ) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )
    PEXTRA(JL,:,53) =  ( PMFLX(JL,:,2) * (PSLGUH(JL,:,2) - ZSLGENH(JL,:)) + &
                    &    PMFLX(JL,:,3) * (PSLGUH(JL,:,3) - ZSLGENH(JL,:)) ) / (ZCFNC1(JL,:) * ZMGEOM(JL,:) * ZRG )
    
  ENDDO !JL
  
  ENDIF !LLDIAG


IF (LHOOK) CALL DR_HOOK('VDFHGHTN',1,ZHOOK_HANDLE)

END SUBROUTINE VDFHGHTN
