!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFBUOYSORT ( KIDIA    , KFDIA   , KLON    , KLEV     , KDRAFT   , &
                       & PAPM1    , PGEOM1  , PGEOH   , &
                       & PQTM     , PSLGM   , &
                       & PFRACB   , KPLCL   , KPTOP   , KPLZB    , &
                       & PQTUH    , PSLGUH  , PWU2H   , PUUH     , PVUH     , &
! DIAGNOSTIC OUTPUT
                       & PEXTR2   , KFLDX2  , PEXTRA  , KLEVX    , KFLDX    , &
!                   
                       & PABULK   , PWBULK  , PQTBULK , PSLGBULK , PUBULK   , PVBULK )

!     ------------------------------------------------------------------
!
!**   *VDFBUOYSORT* - BUOYANCY SORTING ROUTINE
!          
!          Roel Neggers, KNMI, 29/4/2008
!
!     PURPOSE     
!     -------     
!
!          
!     INTERFACE
!     ---------
!
!     *VDFBUOYSORT IS CALLED BY *VDFHGHTN*
!
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!
!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KDRAFT*       NUMBER OF UPDRAFTS
!
!     INPUT PARAMETERS (REAL):
!
!     *PAPM1*        PRESSURE ON FULL LEVELS AT T-1                   PA
!     *PGEOM1*       GEOPOTENTIAL AT T-1                              M2/S2
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                       M2/S2
!     *PQTM*         MEAN STATE TOTAL SPECIFIC HUMIDITY       (FL)    KG/KG
!     *PSLGM*        MEAN STATE LIQUID STATIC ENERGY (SLG)    (FL)    M2/S2
!     *PFRACB*       FRACTIONS OF UPDRAFTS                            0..1
!     *PQTUH*        UPDRAFT TOTAL SPECIFIC HUMIDITY          (HL)    KG/KG
!     *PSLGUH*       UPDRAFT LIQUID STATIC ENERGY (SLG)       (HL)    M2/S2
!     *PWU2H*        UPDRAFT VERTICAL KINETIC ENERGY          (HL)    M2/S2
!     *PUUH*         UPDRAFT X-MOMENTUM                       (HL)    M/S
!     *PVUH*         UPDRAFT Y-MOMENTUM                       (HL)    M/S
!
!     INPUT PARAMETERS (LOGICAL):
!
!
!     OUTPUT PARAMETERS (REAL):
!
!     *PABULK*       BULK UPDRAFT FRACTION
!     *PWBULK*       BULK UPDRAFT VERTICAL VELOCITY
!     *PQTBULK*      BULK UPDRAFT TOTAL SPECIFIC HUMIDITY
!     *PSLGBULK*     BULK UPDRAFT LIQUID STATIC ENERGY (SLG)
!     *PUBULK*       BULK UPDRAFT X-MOMENTUM
!     *PVBULK*       BULK UPDRAFT Y-MOMENTUM
!
!     OUTPUT PARAMETERS (INTEGER):
!
!
!     METHOD
!     ------
!
!

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM   , JPRB

! USE YOMHOOK   ,ONLY : LHOOK  , DR_HOOK

USE YOMCST   , ONLY : RG     , RCPD    , RLVTT

USE YOEVDF   , ONLY : LLDIAG

IMPLICIT NONE



!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KPLCL(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTOP(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPLZB(KLON,KDRAFT)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGM(KLON,KLEV)
 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRACB(KLON,KDRAFT) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWU2H(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUUH(KLON,0:KLEV,KDRAFT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVUH(KLON,0:KLEV,KDRAFT) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PABULK(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWBULK(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQTBULK(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGBULK(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUBULK(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVBULK(KLON,0:KLEV)

!diagnostic output
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)



!*         0.2    LOCAL VARIABLES
      
INTEGER(KIND=JPIM), PARAMETER ::  IPDF = 2
      
REAL(KIND=JPRB) ::    ZDUMMY       , ZFRC        , ZQSAT         , ZTHVM         ,ZQL
     
REAL(KIND=JPRB) ::    ZQTXDEF(KLON,0:KLEV) , ZDQTXDEFDZ(KLON,0:KLEV), &
                 &    ZQTX(KLON,0:KLEV)    , &
                 &    ZBULKDQTXDEFDZBOT(KLON) , IBULKBOT(KLON) , &
                 &    ZBULKDQTXDEFDZTOP(KLON) , IBULKTOP(KLON)
                 
REAL(KIND=JPRB) ::    ZS1(3)       , ZS0(3)      , ZNS(3)        , &
                 &    ZSZB(3)      , ZMZB0(3)    , ZMZB1(3)      , &
                 &    ZMZBSAT0(3)  , ZMZBSAT1(3) , &
                 &    ZXTHV0(3)    , ZXTHV1(3)   , ZXTHVFLAT0(3) , ZXTHVFLAT1(3) , XDIR    , &
                 &    ZNORM        , ZTHV0       , ZTHV1         , &
                 &    ZX0(3)       , ZX1(3)      , ZNX(3)        , ZXMEAN(3)     , ZXMZB(3), &
                 &    ZSIGX        , ZXDIST      , ZXMIN         , ZXBAR

REAL(KIND=JPRB) ::    ZAGRADLCL(KLON), ZAGRADTOP(KLON), &
                 &    ZAGRADLCLMIN, ZAGRADLCLMAX, &
                 &    ZAGRADTOPMIN, ZAGRADTOPMAX, &
                 &    ZSIGW(KLON,0:KLEV)

REAL(KIND=JPRB) ::    ZRG,ZCL,ZASTAR,ZDZ,ZDELTAMINEPS, &
                 &    ZFRAC, ZFAC, ZDUMR
                 
INTEGER(KIND=JPIM) :: JK, JL, JP, JD

LOGICAL ::            LLPBL(KLON,KLEV)


! REAL(KIND=JPRB) ::    ZHOOK_HANDLE

! 
! #include "vdfpdftable.intfb.h"
! #include "vdfsat.intfb.h"
! #include "vdfthermo.intfb.h"



!     ------------------------------------------------------------------
!
!*         1.     INITIALIZE VARIABLES
!                 --------------------
!

! IF (LHOOK) CALL DR_HOOK('VDFBUOYSORT',0,ZHOOK_HANDLE)

ZRG = 1._JPRB/RG

!-- allowed ranges for the dimensionless gradients --
ZAGRADLCLMIN = -10._JPRB
!ZAGRADLCLMAX =   0._JPRB
ZAGRADLCLMAX =  -2._JPRB

ZAGRADTOPMIN = -2.6_JPRB   
!ZAGRADTOPMAX = 0.41_JPRB    
ZAGRADTOPMAX = 0._JPRB    


!-- reset --
ZAGRADLCL(:)  = 0._JPRB    
ZAGRADTOP(:)  = 0._JPRB    


ZDUMMY = 0._JPRB
ZNORM  = 0._JPRB


DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZQTXDEF(JL,JK)      = 0._JPRB
    ZDQTXDEFDZ(JL,JK)   = 0._JPRB
    ZQTX(JL,JK)         = 0._JPRB
  ENDDO
ENDDO
 
      
DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    ZSIGW(JL,JK) = 0._JPRB
    
    PABULK(JL,JK)     = 0._JPRB
    PWBULK(JL,JK)     = 0._JPRB
    PQTBULK(JL,JK)    = 0._JPRB
    PSLGBULK(JL,JK)   = 0._JPRB
    PUBULK(JL,JK)     = 0._JPRB
    PVBULK(JL,JK)     = 0._JPRB
    
  ENDDO
ENDDO
   
    
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    LLPBL(JL,JK) = KPTOP(JL,1)/=0 .AND. JK>=KPTOP(JL,1)
  ENDDO
ENDDO


!     ------------------------------------------------------------------
!
!*         2.     MOIST ZERO BUOYANCY POINT ON MIXING LINE
!                 ----------------------------------------
!      
      
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
  
    IF ( LLPBL(JL,JK) ) THEN                            !only within PBL
      
      !-- mean state thv --
      CALL VDFTHERMO( PSLGM(JL,JK)/RCPD,1000._JPRB*PQTM(JL,JK), ZQSAT, ZQL, PAPM1(JL,JK), PGEOM1(JL,JK) )
      ZQL = ZQL / 1000._JPRB
      ZTHVM = ( PSLGM(JL,JK) + RLVTT*ZQL) * (1._JPRB + 0.61_JPRB*PQTM(JL,JK) - 1.61_JPRB * ZQL) / RCPD
      !write(0,'(6f7.2)') 1000._JPRB*PQTM(JL,JK),ZQSAT,1000._JPRB*ZQL,ZTHVM,PSLGM(JL,JK)/1000._JPRB


      !-- vector of linearized dry saturation curve --
      ZS0   = (/ PSLGM(JL,JK)/RCPD , 0._JPRB, 0._JPRB /)
      ZS1   = ZS0 + (/-1._JPRB,0._JPRB, 0._JPRB/)  !-- minus 1K --
      CALL VDFTHERMO( ZS1(1), 0._JPRB, ZS1(2), ZDUMMY, PAPM1(JL,JK), PGEOM1(JL,JK) )
      CALL VDFTHERMO( ZS0(1), 0._JPRB, ZS0(2), ZDUMMY, PAPM1(JL,JK), PGEOM1(JL,JK) )
      ZNORM = SQRT(DOT_PRODUCT(ZS1(:) - ZS0(:),ZS1(:) - ZS0(:)))
      ZNS(:) = (ZS1(:) - ZS0(:))  / ZNORM
      
      
      !-- zero buoyancy point on saturation curve --
      ZXTHVFLAT0  = (/ ZS0(2), 0._JPRB, 0._JPRB/)
      ZXTHVFLAT1  = (/ ZS1(2), 0._JPRB, 0._JPRB/)
      ZTHV0       = ZS0(1) * (1._JPRB + 0.00061_JPRB*ZS0(2) )
      ZTHV1       = ZS1(1) * (1._JPRB + 0.00061_JPRB*ZS1(2) )
      ZXTHV0      = (/ ZS0(2), ZTHV0-ZTHVM, 0._JPRB/)
      ZXTHV1      = (/ ZS1(2), ZTHV1-ZTHVM, 0._JPRB/)
      !write(0,'(3f7.2,a,f7.2)') ZTHVM, ZTHV1, ZTHV0, "  -  ",ZTHV1-ZTHVM
      
      ZSZB(:) = (/ 0._JPRB, 0._JPRB, 0._JPRB/)
      CALL VDFSAT(ZXTHV1,ZXTHV0,ZXTHVFLAT1,ZXTHVFLAT0,ZSZB,XDIR)  
      !write(0,'(a,5f7.2)') "    ",ZSZB(1),ZS1(2),ZS0(2)
      ZSZB(:) = ZS0(:) + ZNS(:) * ( ZSZB(1) - ZS0(2) ) / ZNS(2)
      !write(0,'(a,3f7.2)') "    ",ZSZB(:)
      
      
      !-- retrieve moist zero buoyancy line --      
      ZMZB0(:)    = ZSZB(:)
      
      ZMZBSAT0(:) = ZMZB0(:) + 2._JPRB * ZNS(:)    !minus 2 unit lengths along the sat curve
      ZTHV0       = ZMZBSAT0(1) * (1._JPRB + 0.00061_JPRB*ZMZBSAT0(2) )
      
      ZMZBSAT1(:) = ZMZBSAT0(:) + (/0._JPRB, 4._JPRB, 0._JPRB /)
      CALL VDFTHERMO( ZMZBSAT1(1), ZMZBSAT1(2), ZQSAT, ZQL, PAPM1(JL,JK), PGEOM1(JL,JK) )  !this point is cloudy per definition
      ZQL = ZQL / 1000._JPRB
      ZTHV1       = (ZMZBSAT1(1) + RLVTT*ZQL/RCPD) * (1._JPRB + 0.00061_JPRB*ZMZBSAT1(2) - 1.61_JPRB*ZQL )
      
      ZXTHVFLAT0  = (/ ZMZBSAT0(2), 0._JPRB, 0._JPRB/)
      ZXTHVFLAT1  = (/ ZMZBSAT1(2), 0._JPRB, 0._JPRB/)
      ZXTHV0      = (/ ZMZBSAT0(2), ZTHV0-ZTHVM, 0._JPRB/)
      ZXTHV1      = (/ ZMZBSAT1(2), ZTHV1-ZTHVM, 0._JPRB/)
      !write(0,'(3f7.2,a,f7.2)') ZTHVM, ZTHV0, ZTHV1
      
      ZMZB1(:) = (/ 0._JPRB, 0._JPRB, 0._JPRB/)
      CALL VDFSAT(ZXTHV1,ZXTHV0,ZXTHVFLAT1,ZXTHVFLAT0,ZMZB1,XDIR)  
      !write(0,'(a,3f7.2)') "    ",ZMZB1
      !write(0,'(a,3f7.2)') "    ",ZMZB1(1),ZMZBSAT1(2),ZMZBSAT0(2)
      ZMZB1(:) = ZMZBSAT0(:) + (/ 0._JPRB, ZMZB1(1)-ZMZBSAT0(2), 0._JPRB /)
      !write(0,'(a,3f7.2,a,3f7.2,a)') "(",ZMZB0,")-(",ZMZB1,")"
      
      !this part is for checking zero bouyancy only
      CALL VDFTHERMO( ZMZB1(1), ZMZB1(2), ZQSAT, ZQL, PAPM1(JL,JK), PGEOM1(JL,JK) )  !this point is cloudy per definition
      ZQL = ZQL / 1000._JPRB
      ZTHV1       = (ZMZB1(1) + RLVTT*ZQL/RCPD) * (1._JPRB + 0.00061_JPRB*ZMZB1(2) - 1.61_JPRB*ZQL )
      !write(0,'(a,2f7.2)') "    buoy check:",ZTHVM,ZTHV1


      !-- updraft PDF axis: lateral mixing line --
!      ZX1(:)    = (/ PSLGUH(JL,JK,3)/RCPD, 1000._JPRB * PQTUH(JL,JK,3), 0._JPRB /)
      ZX1(:)    = (/ PSLGUH(JL,JK,1)/RCPD, 1000._JPRB * PQTUH(JL,JK,1), 0._JPRB /)
      ZX0(:)    = (/ PSLGM(JL,JK)/RCPD   , 1000._JPRB * PQTM(JL,JK)   , 0._JPRB /)
      ZXMEAN(:) = ZX1(:)
      ZNORM     = SQRT(DOT_PRODUCT(ZX1(:) - ZX0(:),ZX1(:) - ZX0(:)))
      IF (ZNORM.EQ.0.) THEN
        !-- Protection against zero length: assume mixing line = dry zero buoyancy line --
        print '(a)',"zero length mixing line vector - terror!"
      ENDIF  
      ZNX(:) = (ZX1(:) - ZX0(:))  / ZNORM
      
      
      !-- intersection point of moist zero buoy line and updraft PDF axis
      CALL VDFSAT(ZMZB1,ZMZB0,ZX1,ZX0,ZXMZB,XDIR)  
      !write(0,'(a,3f7.2,a,3f7.2,a,3f7.2)') "    (",ZMZB0,")-(",ZMZB1,")   ",ZXMZB
      
      
      !-- difference between mean and moist zero buoyancy point --
      ZQTXDEF(JL,JK)   = ZXMZB(2) - ZX0(2)
      !ZSLGMZB(JL,JK)   = ZXMZB(1) * RCPD
      ZQTX(JL,JK)      = ZXMZB(2) / 1000._JPRB

      
    ENDIF
    
  ENDDO
ENDDO
   
      

!     --------------------------------------------------------------------------
!
!*         3.     BULK UPDRAFT FRACTION & BUOYANCY SORTING
!                 ----------------------------------------
!      

DO JL=KIDIA,KFDIA
  ZBULKDQTXDEFDZTOP(JL) = 0._JPRB  
  ZBULKDQTXDEFDZBOT(JL) = 0._JPRB  
  IBULKTOP(JL)          = 0
  IBULKBOT(JL)          = 0
ENDDO


DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
  
    IF ( LLPBL(JL,JK) ) THEN   
    
      ZDQTXDEFDZ(JL,JK) = ( ZQTXDEF(JL,JK) - ZQTXDEF(JL,JK+1)  ) / &
                        & ( PGEOH(JL,JK)   - PGEOH(JL,JK+1)    )
                      
      !--- Vertical gradients of (qtx-qtm) at layer boundaries ---
      !IF ( JK>=KPLZB(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN   
      IF ( JK>=KPLZB(JL,1) .AND. JK>=KPLCL(JL,3)-3 .AND. JK<=KPLCL(JL,3)  ) THEN   
        ZBULKDQTXDEFDZBOT(JL) = ZBULKDQTXDEFDZBOT(JL) + ZDQTXDEFDZ(JL,JK)
        IBULKBOT(JL)          = IBULKBOT(JL) + 1
      ENDIF

      !IF ( JK>=KPLZB(JL,1) .AND. JK<=KPLCL(JL,3) ) THEN   
      IF ( JK>=KPLZB(JL,1) .AND. JK<=KPLZB(JL,1)+3 .AND. JK<=KPLCL(JL,3)  ) THEN   
        ZBULKDQTXDEFDZTOP(JL) = ZBULKDQTXDEFDZTOP(JL) + ZDQTXDEFDZ(JL,JK)
        IBULKTOP(JL)          = IBULKTOP(JL) + 1
      ENDIF

      IF ( JK>=KPTOP(JL,3) .AND. JK<=KPLCL(JL,3)  ) THEN   
        ZSIGW(JL,JK) = PWU2H(JL,JK,1)**0.5 - PWU2H(JL,JK,3)**0.5
      ENDIF  
      
    ENDIF
    
!    IF (LLDIAG) THEN
!      PEXTRA(JL,JK,8)  = 1000._JPRB *ZQTX(JL,JK)
!      PEXTRA(JL,JK,9)  = ZQTXDEF(JL,JK)
!      PEXTRA(JL,JK,10) = ZDQTXDEFDZ(JL,JK)  
!      PEXTRA(JL,JK,11) = ZSIGW(JL,JK)
!    ENDIF
          
  ENDDO
ENDDO


DO JL=KIDIA,KFDIA
  
  
  !--- Bulk cloud-layer gradients ---
  IF (IBULKBOT(JL)>0) THEN
    ZBULKDQTXDEFDZBOT(JL) = ZBULKDQTXDEFDZBOT(JL) / IBULKBOT(JL)
    ZBULKDQTXDEFDZBOT(JL) = ZBULKDQTXDEFDZBOT(JL) * ( PGEOH(JL,KPLZB(JL,1)) - PGEOH(JL,KPLCL(JL,3)) )
  ENDIF
  IF (IBULKTOP(JL)>0) THEN
    ZBULKDQTXDEFDZTOP(JL) = ZBULKDQTXDEFDZTOP(JL) / IBULKTOP(JL)
    ZBULKDQTXDEFDZTOP(JL) = ZBULKDQTXDEFDZTOP(JL) * ( PGEOH(JL,KPLZB(JL,1)) - PGEOH(JL,KPLCL(JL,3)) )
  ENDIF
  
  
  !--- Normalized gradients of a at layer boundaries ---
  ZAGRADLCL(JL) = -1.8_JPRB * ZBULKDQTXDEFDZBOT(JL) / MAX( 1._JPRB, ZQTXDEF(JL,KPLCL(JL,3)) )
  ZAGRADLCL(JL) = MAX( ZAGRADLCLMIN, ZAGRADLCL(JL) )
  ZAGRADLCL(JL) = MIN( ZAGRADLCLMAX, ZAGRADLCL(JL) )
  
  ZAGRADTOP(JL) = -1.8_JPRB * ZBULKDQTXDEFDZTOP(JL) / MAX( 1._JPRB, ZQTXDEF(JL,KPLZB(JL,1)) )
  ZAGRADTOP(JL) = MAX( ZAGRADTOPMIN, ZAGRADTOP(JL) )
  ZAGRADTOP(JL) = MIN( ZAGRADTOPMAX, ZAGRADTOP(JL))
  
  !RN testing: specified mass flux decay
  !ZAGRADLCL(JL) = -4._JPRB
  !ZAGRADTOP(JL) = -2._JPRB
  
  
  !--- Diagnostics ---
!  IF (LLDIAG) THEN
!    PEXTR2(JL,1) = ZBULKDQTXDEFDZBOT(JL)
!    PEXTR2(JL,2) = ZBULKDQTXDEFDZTOP(JL)
!    PEXTR2(JL,3) = ZAGRADLCL(JL)
!    PEXTR2(JL,4) = ZAGRADTOP(JL)
!  ENDIF
  
  
ENDDO



DO JL=KIDIA,KFDIA
  PABULK(JL,KLEV) = PFRACB(JL,3)
ENDDO

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    
    ZFAC   = 0._JPRB
    
    IF ( LLPBL(JL,JK) ) THEN   
    
      IF ( JK>=KPLCL(JL,3)) THEN
        
        PABULK(JL,JK)     = PABULK(JL,JK+1)
      
      ELSE IF ( JK>=KPTOP(JL,3)) THEN
            
        ZCL = ( PGEOH(JL,JK)          - PGEOH(JL,KPLCL(JL,3)) ) / &
            & ( PGEOH(JL,KPTOP(JL,3)) - PGEOH(JL,KPLCL(JL,3)) )
        ZASTAR = (1._JPRB - ZCL) * ZAGRADLCL(JL) + ZCL * ZAGRADTOP(JL)
        ZASTAR = EXP(ZASTAR)
                 
        ZDZ           = ( PGEOH(JL,JK) - PGEOH(JL,JK+1) ) * ZRG
        ZDELTAMINEPS  = LOG(ZASTAR) / MAX( 1._JPRB , ( PGEOH(JL,KPTOP(JL,3))-PGEOH(JL,KPLCL(JL,3)) )*ZRG )
        PABULK(JL,JK) = PABULK(JL,JK+1) * EXP( ZDZ * ZDELTAMINEPS ) 
            
        ZFRAC = PABULK(JL,JK) / PFRACB(JL,3)
        CALL VDFPDFTABLE (ZFRAC, ZFAC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor
        
        ZFAC = ZFAC / 3.9581_JPRB
        
        PWBULK(JL,JK)   = PWU2H(JL,JK,3)**0.5 + ZFAC * ZSIGW(JL,JK)
        PQTBULK(JL,JK)  = (1._JPRB - ZFAC) * PQTUH(JL,JK,3)  + ZFAC * PQTUH(JL,JK,1)
        PSLGBULK(JL,JK) = (1._JPRB - ZFAC) * PSLGUH(JL,JK,3) + ZFAC * PSLGUH(JL,JK,1)
        PUBULK(JL,JK)   = (1._JPRB - ZFAC) * PUUH(JL,JK,3)   + ZFAC * PUUH(JL,JK,1)
        PVBULK(JL,JK)   = (1._JPRB - ZFAC) * PVUH(JL,JK,3)   + ZFAC * PVUH(JL,JK,1)
        
      ENDIF
      
      
    ENDIF
          
!    IF (LLDIAG) THEN
!      PEXTRA(JL,JK,12) = PABULK(JL,JK)
!      PEXTRA(JL,JK,13) = PWBULK(JL,JK)
!      PEXTRA(JL,JK,14) = PABULK(JL,JK) * PWBULK(JL,JK)
!      PEXTRA(JL,JK,7)  = ZFAC
!    ENDIF
  
  ENDDO
ENDDO
  
        
        
! IF (LHOOK) CALL DR_HOOK('VDFBUOYSORT',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE VDFBUOYSORT
      
