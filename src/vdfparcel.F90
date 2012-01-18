!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFPARCEL (KIDIA   , KFDIA   , KLON    , KLEV    , KDRAFT  , &
                    & PGEOH   , PGEOM1  , PAPHM1  , &
                    & PUM1    , PVM1    , PQTM1   , PSLGM1  , PTVEN   , &
                    & PUUH    , PVUH    , PSLGUH  , PQTUH   , PWU2H   , PQCUH  , PBUOF , & 
                    & PQUH    , PTUH    , PEPS    , PFACEXC , &
                    & PZPLCL  , KPLCL   , PZPTOP  , KPTOP   , KPLZB   , &
                    & KD      , PUPGENL , PUPGENN , &
                    & PTAUEPS , PW2THRESH, LDDONE , KPBLTYPE)  
!     ------------------------------------------------------------------

!**   *VDFPARCEL* - VERTICAL INTEGRATION FOR PARCEL ASCENT
!
!     based on original VDFHGHTN.F90 (CY29R1 and earlier) by
!             A.P. SIEBESMA    30/06/99  
!             M. Ko"hler       3/12/2004 
!             Roel Neggers     12/04/2005     separated from VDFHGHTN
!                              15/10/2005     pressure term added
!                              12/04/2006     debugged interpolation of LCL height
!                              30/11/2006     updraft precipitation generation added
!                                             level of zero buoyancy determination


!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT AND UPDRAFT FIELDS

!     INTERFACE
!     ---------

!     *VDFPARCEL* IS CALLED BY *VDFHGHTN*

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
!                    (4: downdrafts .. to be done)
!     *KD*           Draft index
!     *KPBLTYPE*    -1: not defined yet
!                    0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: stratocumulus
!                    3: shallow cumulus
!                    4: deep cumulus



!     INPUT PARAMETERS (REAL):

!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                   M2/S2
!     *PGEOM1*       GEOPOTENTIAL AT T-1                          M2/S2
!     *PAPHM1*       PRESSURE AT HALF LEVEL AT T-1                PA
!     *PUM1*         X-VELOCITY COMPONENT AT T-1                  M/S
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1                  M/S
!     *PQTM1*        MEAN SPECIFIC TOTAL WATER AT FULL LEVEL      KG/KG
!     *PSLGM1*       MEAN LIQUID STATIC ENERGY AT FULL LEVEL      M2/S2
!     *PTVEN*        ENVIRONMENTAL VIRTUAL TEMPERATURE            K
!     *PTAUEPS*      UPDRAFT ENTRAINMENT TIMESCALE                S
!     *PW2THRESH*    THRESHOLD UPDRAFT VELOCITY SQUARED           M2/S2


!     INPUT PARAMETERS (LOGICAL):
!     *LDDONE*       PARCEL LIFTOFF CONFIRMATION (.TRUE. means 'don't launch parcel')


!     OUTPUT PARAMETERS (REAL):

!     *PUUH*         UPDRAFT X-MOMENTUM
!     *PVUH*         UPDRAFT Y-MOMENTUM
!     *PSLGUH*       UPDRAFT GENERALIZED LIQUID STATIC ENERGY (SLG)
!                    AT HALF LEVEL                                   M2/S2
!     *PQTUH*        UPDRAFT TOTAL SPECIFIC HUMIDITY AT HALF LEVEL   KG/KG
!     *PWU2H*        UPDRAFT VERTICAL VELOCITY SQUARE AT HALF LEVEL  M2/S2
!     *PQCUH*        UPDRAFT LIQUID WATER AT HALF LEVEL              KG/KG
!     *PQUH*         UPDRAFT SPECIFIC HUMIDITY AT HALF LEVEL         KG/KG
!     *PTUH*         UPDRAFT TEMPERATURE AT HALF LEVEL               K
!     *PBUOF*        UPDRAFT BUOYANCY AT FULL LEVEL                  M/S2
!     *PEPS*         UPDRAFT ENTRAINMENT RATE                        1/M
!     *PZPLCL*       HEIGHT OF LIFTING CONDENSATION LEVEL OF UPDRAFT          M
!     *PZPTOP*       HEIGHT OF LEVEL OF ZERO KINETIC ENERGY (W=0) OF UPDRAFT  M
!     *PUPGENL*      UPDRAFT RAIN GENERATION                         KG/KG /S
!     *PUPGENN*      UPDRAFT SNOW GENERATION                         KG/KG /S
!
!     OUTPUT PARAMETERS (INTEGER):

!     *KPLCL*        FIRST HALF LEVEL ABOVE REAL HEIGHT OF UPRAFT LCL
!     *KPTOP*        HIGHEST HALF LEVEL BELOW PZTOP, AND
!                    UPDRAFT TOP FULL LEVEL (PZTOP IS WITHIN THAT LAYER)
!     *KPLZB*        LEVEL OF UPRAFT ZERO BUOYANCY (HIGHEST FULL LEVEL THAT IS POS. BUOYANT)

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

!#include "tsmbkind.h"
use garbage, only : foealfa, cuadjtq
USE PARKIND1  ,ONLY : JPIM     , JPRB

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOETHF   , ONLY : R2ES     , R3LES    , R3IES    , &
                     &R4LES    , R4IES    , R5LES    , R5IES     , RVTMP2  , & 
                     &RALVDCP  , RALSDCP  , &
		     &RTWAT    , RTICE    , RTICECU  , R5ALVCP   , R5ALSCP , & 
                     &RTWAT_RTICE_R       , RTWAT_RTICECU_R
                     
USE YOMCST   , ONLY : RG  , RLSTT , RCPD , RLVTT , RETV ,RTT         



IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KD
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPLCL(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPTOP(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPLZB(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPBLTYPE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM1 (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTVEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PUUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PVUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PSLGUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PQTUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PWU2H(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PQCUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PQUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PEPS(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PFACEXC(KLON,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PZPTOP(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PZPLCL(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)      :: PTAUEPS
REAL(KIND=JPRB)   ,INTENT(IN)      :: PW2THRESH
LOGICAL           ,INTENT(INOUT)   :: LDDONE(KLON,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PBUOF(KLON,KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PUPGENL(KLON,KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PUPGENN(KLON,KLEV,KDRAFT)

!*         0.2    LOCAL VARIABLES

LOGICAL ::            LLDOIT(KLON), LLCLOUD(KLON), LLPREC

REAL(KIND=JPRB) ::    ZRG   , ZDZ   , ZTVUF , ZQUF    , ZQCUF   , ZMU, ZB, &
                    & ZQLUF , ZQIUF , ZQTUF , ZTUF  , ZSLGUF  , &
		    & ZMIX  (KLON,0:KLEV), ZMIXW (KLON,0:KLEV) !,ZBUOF(KLON,KLEV)

REAL(KIND=JPRB) ::    ZWUH  (KLON,0:KLEV)   , ZQUH  (KLON,0:KLEV)    , ZPH  (KLON), &
                    & ZTTEMP(KLON,KLEV)     , ZQTEMP(KLON,KLEV)      
		    
REAL(KIND=JPRB) ::    ZALFAW  , ZFACW   , ZFACI   , ZFAC    , ZTEMP ,&
                    & ZESDP   , ZCOR    , ZDQSDTEMP(KLON)   , ZQS(KLON,0:KLEV), &
                    & ZPGENUP , ZWUHTEMP(KLON), ZTAUEPS(KLON), ZCEPSZ(KLON)

INTEGER(KIND=JPIM) :: IS, JK, JL, JKMAX, JKMIN

REAL(KIND=JPRB) ::    ZDQSUDZ, ZDQTUDZ, ZLCLFAC(KLON), ZEPSCFLFAC, ZQLWORK, ZLCRIT

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

!DIR$ VFUNCTION EXPHF
#include "fcttre.h"
#include "cuadjtq.intfb.h"


!     -----------------------------------------------------------------

!*         1.     SET SOME CONSTANTS
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFPARCEL',0,ZHOOK_HANDLE)

!  Constant of proportionality between updraft induced pressure term and
!     updraft vertical acceleration (Siebesma, Soares and Teixeira, JAS 2007)
!ZMU = 0._JPRB  
ZMU = 0.15_JPRB   !cy32r3

!  Constant of proportionality between kinematic and thermodynamic entrainment
!ZB = 1.0_JPRB
ZB = 0.5_JPRB   !cy32r3

!  CFL criterion factor for updraft lateral entrainment
ZEPSCFLFAC = 0.6_JPRB   
!ZEPSCFLFAC = 1.0_JPRB

!  Optimization
ZRG         = 1.0_JPRB/RG 

!  Updraft precipitation switch
LLPREC=.TRUE.
!LLPREC=.FALSE.

!  Critical updraft condensate [kg/kg] in Sundqvist precipitation generation
!ZLCRIT = 0.0005_JPRB  
ZLCRIT = 0.001_JPRB     !cy32r3
!ZLCRIT = 0.0015_JPRB   



!     -----------------------------------------------------------------

!*         2.     SOME FINAL INITIALIZATION
!                 ------------------------

  DO JL=KIDIA,KFDIA
    LLCLOUD(JL)        = .FALSE.
    PBUOF (JL,KLEV,KD) = 0.0_JPRB
    ZLCLFAC(JL)        = 0.0_JPRB
    KPLZB(JL,KD)       = KLEV-1
    ZTAUEPS(JL)        = PTAUEPS
    
    !  1/z scaling factor in eps
    ZCEPSZ(JL) = 0._JPRB
    !IF (KD==2 .AND. KPBLTYPE(JL)==1) THEN
    !  ZCEPSZ(JL) = 0.4_JPRB
    !ENDIF
    !IF (.NOT.LDDONE(JL,KD)) THEN
    !  ZCEPSZ(JL) = ( 2._JPRB * PFACEXC(JL,KD) )**(-1._JPRB)
    !ENDIF  
    !write(0,'(a,i,2f)') '    ',KD,PFACEXC(JL,KD),ZCEPSZ(JL)
    
  ENDDO
    
  !integration over total depth for now: 
  !    Note that limiting JKMIN to KPLCL/KPTOP for KD=2 could speed things up a bit..
  JKMAX = KLEV-2
  JKMIN = 1


!     -----------------------------------------------------------------

!*         3.     VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
!                 -----------------------------------------------
  
  DO JK=JKMAX,JKMIN,-1
  
    IS=0
    DO JL=KIDIA,KFDIA
      IF (.NOT.LDDONE(JL,KD)) THEN
        IS            = IS+1


!*         3.1  Updraft entrainment
        
        ZWUH(JL,JK+1) = SQRT( MAX( PWU2H(JL,JK+1,KD), 0.01_JPRB) ) ! w,up > 0.1 m/s (safety)
        PEPS(JL,JK+1,KD) = 1.0_JPRB / ( ZWUH(JL,JK+1) * ZTAUEPS(JL) ) !& ! eps=1/(w,up*tau)
        
        !RN numerical entrainment limiter: maximally c_e/Dz
        ZDZ           = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
        PEPS(JL,JK+1,KD) = MIN( ZEPSCFLFAC/ZDZ, PEPS(JL,JK+1,KD) )
        PEPS(JL,JK+1,KD) = MAX( PEPS(JL,JK+1,KD), ZCEPSZ(JL) * RG/PGEOH(JL,JK+1) )
        

!*         3.2  Ascent of slg and qt (exact)
 
        ZMIX(JL,JK+1) = exp( - ZDZ * PEPS(JL,JK+1,KD) )
        PQTUH(JL,JK,KD)  = ( PQTUH (JL,JK+1,KD) - PQTM1 (JL,JK+1) ) * ZMIX(JL,JK+1) &
                    & + PQTM1 (JL,JK+1)
        PSLGUH(JL,JK,KD) = ( PSLGUH(JL,JK+1,KD) - PSLGM1(JL,JK+1) ) * ZMIX(JL,JK+1) &
                    & + PSLGM1(JL,JK+1)
        PUUH(JL,JK,KD)   = ( PUUH (JL,JK+1,KD)  - PUM1  (JL,JK+1) ) * ZMIX(JL,JK+1) &
                    & + PUM1 (JL,JK+1)
        PVUH(JL,JK,KD)   = ( PVUH (JL,JK+1,KD)  - PVM1  (JL,JK+1) ) * ZMIX(JL,JK+1) &
                    & + PVM1 (JL,JK+1)
                    

!*         3.3  Condensation - diagnose T, qv, ql
!*         (Taylor, becomes more inaccurate for states far from q*
!*          -> double condensation step)

!          cuadjtq initialization (assume qv=qt)

        PQUH(JL,JK,KD)   = PQTUH(JL,JK,KD)
        PQCUH(JL,JK,KD)  = 0.0_JPRB

!          cuadjtq initialization (assume qv=qv(jk+1) for speed)

!       IF ( ZQUH(JL,JK+1) < PQTUH(JL,JK) ) THEN
!         ZQUH(JL,JK) = ZQUH(JL,JK+1)
!         ZQCUH(JL,JK)= PQTUH(JL,JK) - ZQUH(JL,JK+1)
!       ENDIF

        PTUH(JL,JK,KD)   = ( PSLGUH(JL,JK,KD) - PGEOH(JL,JK) + RLVTT*PQCUH(JL,JK,KD) ) &
                    & / RCPD      ! assume liquid phase!
        ZPH(JL)       = PAPHM1(JL,JK)
        ZQTEMP(JL,JK) = PQUH(JL,JK,KD)
        ZTTEMP(JL,JK) = PTUH(JL,JK,KD)
	
      ENDIF
      
      LLDOIT(JL)      = .NOT. LDDONE(JL,KD)
      
      IF ( KD==2 ) THEN    ! condensation not done for dry subcloud thermal
        LLDOIT(JL)    = .FALSE.
      ENDIF
      
    ENDDO


    CALL CUADJTQ &
     & ( KIDIA,    KFDIA,    KLON,     0,       KLEV,&
     &   JK,&
     &   ZPH,      ZTTEMP,   ZQTEMP,   LLDOIT,  4)  


    DO JL=KIDIA,KFDIA
      IF ( LLDOIT(JL) ) THEN
        IF ( ZQTEMP(JL,JK) < PQTUH(JL,JK,KD) ) THEN !allow evaporation up to qt
          PQUH(JL,JK,KD) = ZQTEMP(JL,JK)
          PQCUH(JL,JK,KD)= PQTUH(JL,JK,KD) - PQUH(JL,JK,KD)
          PTUH(JL,JK,KD) = ZTTEMP(JL,JK)
        ELSE                          !case where qv(initial)<qt but qv(final)>qt
          PQUH(JL,JK,KD) = PQTUH(JL,JK,KD)  !(unusual!)
          PQCUH(JL,JK,KD)= 0.0_JPRB
          PTUH(JL,JK,KD) = ( PSLGUH(JL,JK,KD) - PGEOH(JL,JK) + RLVTT*PQCUH(JL,JK,KD) ) &
                    & / RCPD
        ENDIF
      ENDIF
    ENDDO


    DO JL=KIDIA,KFDIA


!*         3.4  Updraft microphysics                *experimental* RN 
!*
      !Precip generation tendency (Sundqvist,1978) 
      !    [kg/kg /s]    in full level below current half level
      IF (LLPREC .AND. LLCLOUD(JL)) THEN

        ZDZ = (PGEOH(JL,JK) - PGEOH(JL,JK+1)) * ZRG
        ZQLWORK = PQCUH(JL,JK,KD)
        !ZQLWORK = PQCUH(JL,JK+1,KD)
        !ZQLWORK = ( PQCUH(JL,JK,KD) + PQCUH(JL,JK+1,KD) ) / 2._JPRB
        ZPGENUP = 0.0015_JPRB * ZQLWORK * (1._JPRB - EXP(-(ZQLWORK/ZLCRIT)**2) )
        ZPGENUP = MIN( ZPGENUP, PQCUH(JL,JK,KD)/ZDZ )              

        ZALFAW = FOEALFA(PTUH(JL,JK,KD))
        PUPGENL(JL,JK+1,KD) = ZPGENUP * ZALFAW
        PUPGENN(JL,JK+1,KD) = ZPGENUP * (1._JPRB - ZALFAW)
        
        !adjust the associated updraft state variables (integrate tendency over layer)      
        PQCUH(JL,JK,KD)  = PQCUH(JL,JK,KD)  - ZPGENUP * ZDZ
        PQTUH(JL,JK,KD)  = PQTUH(JL,JK,KD)  - ZPGENUP * ZDZ
        PSLGUH(JL,JK,KD) = PSLGUH(JL,JK,KD) + ZDZ * &
           & (RLVTT * PUPGENL(JL,JK+1,KD) + RLSTT * PUPGENN(JL,JK+1,KD) )

      ENDIF


!*         3.5  Interpolation of updraft LCL

      IF ( PQCUH(JL,JK,KD) > 0.0_JPRB  .AND.  .NOT. LLCLOUD(JL) ) THEN

        LLCLOUD(JL)   = .TRUE.
	
	!cloud base level is first level with ql>0
        KPLCL(JL,KD)  = JK

        !RN --- new interpolation method incorporating dqt/dz *AND* dqsat/dz ---
        ZDQSUDZ = ( ZQTEMP(JL,JK)   - ZQTEMP(JL,JK+1)   ) * RG / ( PGEOH(JL,JK) - PGEOH(JL,JK+1) )
        ZDQTUDZ = ( PQTUH(JL,JK,KD) - PQTUH(JL,JK+1,KD) ) * RG / ( PGEOH(JL,JK) - PGEOH(JL,JK+1) )
        
        PZPLCL(JL,KD) = PGEOH(JL,JK)*ZRG 
        IF (ZDQSUDZ-ZDQTUDZ.LT.0._JPRB) THEN
          PZPLCL(JL,KD) = PGEOH(JL,JK)*ZRG + PQCUH(JL,JK,KD)/(ZDQSUDZ-ZDQTUDZ)
        ENDIF
        
        IF ( KPLCL(JL,KD) < KLEV ) THEN
          PZPLCL(JL,KD) = MAX( PZPLCL(JL,KD), PGEOH(JL,KPLCL(JL,KD)+1)*ZRG )
        ELSE
          PZPLCL(JL,KD) = MAX( PZPLCL(JL,KD), 0.0_JPRB )
        ENDIF
        
        ZLCLFAC(JL) = ( PGEOH(JL,JK) - PZPLCL(JL,KD)*RG  ) / ( PGEOH(JL,JK) - PGEOH(JL,JK+1) )
        ZLCLFAC(JL) = MAX(0._JPRB, MIN(1._JPRB, ZLCLFAC(JL)) )

      ENDIF


!*         3.6  Updraft buoyancy
!*         (at full level k+1 from interpolation of slg, q, qc)

      IF ( .NOT. LDDONE(JL,KD) ) THEN

        ZSLGUF        = 0.5_JPRB * ( PSLGUH(JL,JK,KD) + PSLGUH(JL,JK+1,KD) )
        ZQTUF         = 0.5_JPRB * ( PQTUH (JL,JK,KD) + PQTUH (JL,JK+1,KD) )
        ZQCUF         = 0.5_JPRB * ( PQCUH (JL,JK,KD) + PQCUH (JL,JK+1,KD) )
	
	! At first full level above real cloud base, interpolate ql between 
	! height of real cloud base and the half level of KPLCL
        IF ( JK == KPLCL(JL,KD) ) THEN
          ZQCUF = ZQCUF * ZLCLFAC(JL)
        ENDIF
        
	ZQUF          = ZQTUF - ZQCUF
        ZTUF          = ( ZSLGUF - PGEOM1(JL,JK+1) & ! preliminary estimate:
                    & + RLVTT * ZQCUF ) / RCPD       ! all liquid 
        ZALFAW        = FOEALFA( ZTUF )
        ZQLUF         = ZALFAW            * ZQCUF
        ZQIUF         = (1.0_JPRB-ZALFAW) * ZQCUF
        ZTUF          = ( ZSLGUF - PGEOM1(JL,JK+1) &
                    & + RLVTT * ZQLUF + RLSTT * ZQIUF ) / RCPD  
        ZTVUF         = ZTUF * ( 1.0_JPRB + RETV * ZQUF - ZQCUF )  
        PBUOF(JL,JK+1,KD) = RG * ( ZTVUF - PTVEN(JL,JK+1) ) / PTVEN(JL,JK+1)
        
        !update level of zero buoyancy
        IF ( PBUOF(JL,JK+1,KD)>0._JPRB ) THEN
          KPLZB(JL,KD) = JK+1
        ENDIF


!*         3.7  Kinetic energy equation (exact)
        
        ZDZ           = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
        ZTEMP         = PBUOF(JL,JK+1,KD) / (ZB*PEPS(JL,JK+1,KD))
!        PWU2H(JL,JK,KD)  = ( PWU2H(JL,JK+1,KD) - ZTEMP ) * ZMIX(JL,JK+1)**2 + ZTEMP
        ZMIXW(JL,JK+1) = exp( - ZDZ * (ZB*PEPS(JL,JK+1,KD)) / (1._JPRB - 2._JPRB*ZMU) )
        PWU2H(JL,JK,KD)  = ( PWU2H(JL,JK+1,KD) - ZTEMP ) * ZMIXW(JL,JK+1)**2 + ZTEMP
        

!*         3.8  Inversion height at w=0  (lin. interpolation in w^2)
        
        IF ( PWU2H(JL,JK,KD) < 0.0_JPRB  .AND.  PWU2H(JL,JK+1,KD) > 0.0_JPRB ) THEN 
	
	  !set top level to last level with positive w
          KPTOP(JL,KD)   = JK+1  
	  
          PZPTOP(JL,KD)   = PGEOH(JL,JK+1) * ZRG &
                    & + ZDZ * PWU2H(JL,JK+1,KD) / ( PWU2H(JL,JK+1,KD) - PWU2H(JL,JK,KD) ) 
                     
        ENDIF

        IF ( PWU2H(JL,JK,KD) < PW2THRESH ) THEN   !allow parcel to overcome layers 
          LDDONE(JL,KD)  = .TRUE.                 !with small negative kin. energy
        ENDIF                                  !but remember last w(z)=0 (z=PZPTOP)


      ENDIF
    ENDDO
    IF (IS == 0) EXIT
    
        
  ENDDO !JK

  
  !protect for updrafts that have reached top level (let's hope this is unnecessary!)
  DO JL=KIDIA,KFDIA
    IF ( .NOT. LDDONE(JL,KD) ) THEN
      LDDONE(JL,KD) = .TRUE.
      KPTOP(JL,KD)  = JKMIN+1
      KPLZB(JL,KD)  = JKMIN+1
      PZPTOP(JL,KD) = PGEOH(JL,JKMIN+1) * ZRG
    ENDIF
    
    !RN testing RICO: 
    !KPLCL(JL,KD) = 81
    
  ENDDO


IF (LHOOK) CALL DR_HOOK('VDFPARCEL',1,ZHOOK_HANDLE)
END SUBROUTINE VDFPARCEL
