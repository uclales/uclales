!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFSTCUCRIT (KIDIA   , KFDIA   , KLON    , KLEV    , KDRAFT  , &
                   &    PTM1    , PSLGM1  , PQTM1   , PAPM1   , &
                   &    PSTABTHRESH, PCLDDEPTH, PBIRTHRESH, PDZCLOUD, &
                   &    KPTOP   , KPBLTYPE , LDNODECP, &
                   &    PSTABILITY  )  
!     ------------------------------------------------------------------

!**   *VDFSTCUCRIT* - CRITERIA FOR STRATOCUMULUS OCCURRENCE
!
!     based on original VDFHGHTN.F90 (CY29R1 and earlier) by
!             A.P. SIEBESMA    30/06/99  
!             M. Ko"hler       3/12/2004 
!     put into separate file by   
!             Roel Neggers     12/04/2005 


!     PURPOSE
!     -------

!     DETERMINE STRATOCUMULUS OCCURRENCE

!     INTERFACE
!     ---------

!     *VDFSTCUCRIT* IS CALLED BY *VDFHGHTN*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KPTOP*        HIGHEST HALF LEVEL BELOW PBL HEIGHT, AND
!                    PBL TOP FULL LEVEL (PZPTOP IS WITHIN THAT LAYER)
!     *KDRAFT*       NUMBER OF EXPLICITLY MODELED DRAFTS - CURRENTLY 3:
!                    1: test parcel
!                    2: rising dry thermals which stop at cloud base or inversion
!                    3: rising dry thermals which become cloudy
!                    (4: downdrafts .. to be done)


!     INPUT PARAMETERS (REAL):

!     *PTM1*         TEMPERATURE AT T-1                               K
!     *PSLGM1*       LIQUID STATIC ENERGY (SLG) AT T-1                K
!     *PQTM1*        TOTAL SPECIFIC HUMIDITY AT T-1                   KG/KG
!     *PAPM1*        PRESSURE AT FULL LEVEL AT T-1                    PA
!     *PSTABTHRESH*  STABILITY CRITERION (Klein & Hartmann criteria)  K                      
!     *PCLDDEPTH*    THRESHOLD CLOUD THICKNESS FOR STCU/CU TRANSITION M
!     *PBIRTHTHRASH* THRESHOLD BIR (TKE DECOUPLING CRITERIA)          1
!     *PDZCLOUD*     CLOUD THICKNESS                                  M


!     INPUT PARAMETERS (LOGICAL):

!     *LDNODECP*     TRUE:  NEVER DECOUPLE
!                    FALSE: MAYBE DECOUPLE



!     OUTPUT PARAMETERS (INTEGER):

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

!#include "tsmbkind.h"
use garbage   ,only : foeewm, foealfa
USE PARKIND1  ,ONLY : JPIM     , JPRB

! USE YOMHOOK   ,ONLY : LHOOK    , DR_HOOK

USE YOMCST   , ONLY : RD       , RG      , RCPD     , RETV     , RLVTT, &
                     &RLSTT    ,RTT

USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,R4IES   ,R5LES    ,R5IES, &
                     &R5ALVCP  ,R5ALSCP  ,RALVDCP  ,RALSDCP  ,RTWAT   ,RTICE    ,RTICECU, &
                     &RTWAT_RTICE_R      ,RTWAT_RTICECU_R

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTOP(KLON,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPBLTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTABTHRESH
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDDEPTH
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBIRTHRESH
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDZCLOUD(KLON) 
LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTABILITY(KLON)


!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::   ZRG,  ZDCTEI(KLON)

!          VARIABLES FOR CLOUD BASE ESTIMATION

REAL(KIND=JPRB) ::    ZQS(KLON,0:KLEV)  , ZALFAW  , ZFACW   , ZFACI   , ZFAC        , &
                    & ZESDP   , ZCOR    , ZDQSDTEMP(KLON)   , ZBETA   , ZDSL(KLEV)  , &
		    & ZDQT(KLEV)
		      
INTEGER(KIND=JPIM) :: IS, JK, JL, JD

INTEGER(KIND=JPIM) :: I700(KLON)

! REAL(KIND=JPRB) ::    ZHOOK_HANDLE


!#include "cuadjtq.intfb.h"

!DIR$ VFUNCTION EXPHF
! #include "fcttre.h"



!     -----------------------------------------------------------------

!*         1.     SET SOME CONSTANTS
!                 --------------------

! IF (LHOOK) CALL DR_HOOK('VDFSTCUCRIT',0,ZHOOK_HANDLE)

! optimization
ZRG         = 1.0_JPRB/RG




!     -----------------------------------------------------------------

!*         2.     PREPARE VARIABLES APPEARING IN SOME CRITERIA 
!                 --------------------------------------------

!*         2.1  stability criteria == theta(700hPa) - theta(sfc)

!          find index I700 of pressure closest to 700hPa

  DO JL=KIDIA,KFDIA
    I700(JL) = 0
  ENDDO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF ( I700(JL) == 0  .AND.  PAPM1(JL,JK) > 70000.0_JPRB ) THEN
        I700(JL) = JK
      ENDIF
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    IF ( I700(JL) > 1 ) THEN
      IF ( ABS(PAPM1(JL,I700(JL)-1)-70000.0_JPRB)  <  &
         & ABS(PAPM1(JL,I700(JL))  -70000.0_JPRB)    ) THEN  
        I700(JL) = I700(JL)-1
      ENDIF

      PSTABILITY(JL) = PTM1(JL,I700(JL)) * ( 1.0e5_JPRB/PAPM1(JL,I700(JL)) ) ** (RD/RCPD) &
                   & - PTM1(JL,KLEV)     * ( 1.0e5_JPRB/PAPM1(JL,KLEV) )     ** (RD/RCPD)  
    ELSE
      PSTABILITY(JL) = 0.0_JPRB
    ENDIF

  ENDDO


!*         2.2  cloud top entrainment instability (CTEI) criteria

  DO JL=KIDIA,KFDIA

   IF ( .FALSE. ) THEN
    IF ( KPBLTYPE(JL) == 2) THEN

      JK            = KPTOP(JL,1)  ! PBL top full level taken as K+1
                                   ! full level above PBL taken as K-1

!          qsat (full level)
      ZQS(JL,JK+1)  = FOEEWM(PTM1(JL,JK+1))/PAPM1(JL,JK+1)
      ZQS(JL,JK+1)  = MIN(0.5_JPRB,ZQS(JL,JK+1))
      ZQS(JL,JK+1)  = ZQS(JL,JK+1)/(1.0_JPRB-RETV*ZQS(JL,JK+1))

!          calculate dqs/dT correction factor (full level)
      ZALFAW        = FOEALFA(PTM1(JL,JK+1))
      ZFACW         = R5LES/((PTM1(JL,JK+1)-R4LES)**2)
      ZFACI         = R5IES/((PTM1(JL,JK+1)-R4IES)**2)
      ZFAC          = ZALFAW*ZFACW+(1.0_JPRB-ZALFAW)*ZFACI
      ZESDP         = FOEEWM(PTM1(JL,JK+1))/PAPM1(JL,JK+1)
      ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
      ZDQSDTEMP(JL) = ZFAC*ZCOR*ZQS(JL,JK+1)

!          CTEI
      ZBETA = ( 1.0_JPRB + (1.0_JPRB+RETV) * PTM1(JL,JK+1) * ZDQSDTEMP(JL) ) &
          & / ( 1.0_JPRB + RLVTT/RCPD * ZDQSDTEMP(JL) )  
      ZDSL(JL)  = PSLGM1(JL,JK-1) - PSLGM1(JL,JK+1)
      ZDQT(JL)  = PQTM1 (JL,JK-1) - PQTM1 (JL,JK+1)

      ZDCTEI(JL) = ZBETA * ZDSL(JL) &
               & + ( ZBETA - RCPD/RLVTT * PTM1(JL,JK+1) ) * RLVTT * ZDQT(JL)  

    ELSE

      ZDCTEI(JL) = 1.0

    ENDIF
   ENDIF




!     -----------------------------------------------------------------

!*         3      STRATOCUMULUS - SHALLOW CUMULUS CRITERIA: 
!*                * CLOUD THICKNESS = 1000M
!*                * STABILITY = 15K
!*                * TKE DECOUPLING, BIR = 0.1
!                 -----------------------------------------


    IF ( .NOT. LDNODECP(JL) .AND.  KPBLTYPE(JL) == 2 ) THEN

!..........stability criteria (Klein & Hartmann 1993)
      IF ( PSTABILITY(JL) < PSTABTHRESH ) THEN 

!..........cloud thickness criteria
!     IF ( PDZCLOUD(JL) > PCLDDEPTH ) THEN

!..........CTEI...
!     IF ( ZDCTEI(JL) < 0 ) THEN

!..........TKE decoupling (or cloud thickness criteria for safety)
!     IF ( PBIR(JL) > PBIRTHRESH .OR. PDZCLOUD(JL) > PCLDDEPTH ) THEN

!..........always decouple
!     IF ( .TRUE. ) THEN

        KPBLTYPE(JL) = 3   !decouple: PBL type 3 (shallow cumulus)

      ENDIF
    ENDIF


  ENDDO !JL




! IF (LHOOK) CALL DR_HOOK('VDFSTCUCRIT',1,ZHOOK_HANDLE)
END SUBROUTINE VDFSTCUCRIT
