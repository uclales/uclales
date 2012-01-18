!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFCLOUD ( KIDIA    , KFDIA   , KLON    , KLEV     , KDRAFT   , &
                    & PAPM1    , PGEOM1  , PGEOH   , &
                    & PQTM     , PSLGM   , &
                    & PFRACB   , KVARTOP , &
                    & PQTUP    , PSLGUP  , &
                    & PQTTEST  , &
                    & PSIGQT2  , PAWAKE , &
                    & PEXTR2   , KFLDX2  , PEXTRA  , KLEVX  , KFLDX , &
                    & PCLDFRAC , PQLAV   , PLDIFF)

!     ------------------------------------------------------------------
!
!**   *VDFCLOUD* - VECTORIZED BIMODAL CLOUD SCHEME
!          
!          As formulated by Neggers et al. (JAS, 2009). 
!
!     PURPOSE     
!     -------     
!
!     DETERMINES THE PBL CLOUD FRACTION AND TOTAL CONDENSATE AS A FUNCTION OF 
!     THE QT VARIANCES, THE MEAN STATE AND THE MULTIPLE UPDRAFT STATES.
!          
!     INTERFACE
!     ---------
!
!     *VDFCLOUD* IS CALLED BY *VDFMAIN*
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
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVELS                      M2/S2
!     *PQTM*         MEAN STATE TOTAL SPECIFIC HUMIDITY       (FL)    KG/KG
!     *PSLGM*        MEAN STATE LIQUID STATIC ENERGY (SLG)    (FL)    M2/S2
!     *PFRACB*       FRACTIONS OF UPDRAFTS                            0..1
!     *PQTUP*        MOIST UPDRAFT TOTAL SPECIFIC HUMIDITY    (FL)    KG/KG
!     *PSLGUP*       MOIST UPDRAFT LIQUID STATIC ENERGY (SLG) (FL)    M2/S2
!     *PQTTEST*      TEST PARCEL TOTAL SPECIFIC HUMIDITY      (FL)    KG/KG
!     *PSIGQT2*      GRIDBOX AVERAGE QT VARIANCE              (FL)    KG^2/KG^2
!
!     INPUT PARAMETERS (LOGICAL):
!
!
!     OUTPUT PARAMETERS (REAL):
!
!     *PCLDFRAC*      CLOUD FRACTION AT FULL LEVEL               0..1
!     *PQLAV*         TOTAL CONDENSATE AT FULL LEVEL             KG/KG
!     *PLDIFF*        CONTRIB TO PBL CONDENSATE BY PASSIVE CLOUDS   KG/KG
!
!     OUTPUT PARAMETERS (INTEGER):
!
!
!     METHOD
!     ------
!
!          The scheme is formulated in {thl,qt} space. The same ED-MF decomposition as
!          done for convective transport is now made in the cloud scheme. The total PDF is
!          split  up into two separate PDFs: a diffusive and an updraft PDF. Both joint-PDFs
!          are linearized: the diffusive PDF is oriented along the zero buoyancy line, 
!          and the updraft PDF along the lateral mixing line between moist updraft state 
!          and the gridbox-mean state.
!
!          For each PDF:
!            - The intersection point of the dry saturation line and the PDF axis is calculated,
!              using a vector calculus method (VDFSAT).
!            - Its width is estimated, using the relations between the moments of a bimodal PDF 
!              (Lewellen & Yoh, JAS 1993),
!
!                 mu = \sum_{i=1}^{2} a_i mu_i
!
!                 sigma^2 + mu^2 = \sum_{i=1}^{2} a_i ( sigma_i^2 + mu_i^2 )
!
!              where mu_i and sigma_i are the 1st and 2nd moments of the PDF i. The updraft
!              PDF variance is derived from the difference in thermodynamic state of two updrafts:
!              the test updraft (used to diagnose PBL height) and the moist updraft.
!            - Finally, the saturated fraction and total condensate of each PDF is estimated,
!              using the  Gaussian PDF lookup table (VDFPDFTABLE)
! 
!          Finally all saturated fractions and condensate
!          For qsat the triple-point routine of Sommeria & Deardorff (1977) is used (VDFTHERMO).
!

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RG       ,RCPD

USE YOEVDF   , ONLY : LLDIAG

IMPLICIT NONE



!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGM(KLON,KLEV)
 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRACB(KLON,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVARTOP(KLON)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTUP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGUP(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTTEST(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGQT2(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLDFRAC(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQLAV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLDIFF(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAWAKE(KLON,KLEV)

!diagnostic output
INTEGER(KIND=JPIM),INTENT(IN) :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)



!*         0.2    LOCAL VARIABLES
      
INTEGER(KIND=JPIM), PARAMETER ::  IPDF = 2
      
REAL(KIND=JPRB) ::    ZQBAR        , ZQMIN       , ZIQBAR(KLON)   
REAL(KIND=JPRB) ::    ZDUMMY, RCPDTMP, ZRG
REAL(KIND=JPRB) ::    ZDZ, ZSGSEFOLD
     
REAL(KIND=JPRB) ::    ZATEST(KLON) , ZAK(KLON)   , ZAUP(KLON)    , &
                 &    ZQTK(KLON,KLEV)      , ZSLGK(KLON,KLEV)    , &
                 &    ZSIGQT2UP(KLON,KLEV) , ZSIGQT2K(KLON,KLEV)

     
REAL(KIND=JPRB) ::    ZX1(IPDF,3)      , ZX0(IPDF,3)    , ZS1(3)           , ZS0(3) , &
                 &    ZNX(IPDF,3)      , ZSMIN(IPDF,3)  , ZXMEAN(IPDF,3)   , &
                 &    ZAPDF(IPDF)      , ZCCPDF(IPDF)   , ZSIGQT2PDF(IPDF) , &
                 &    ZCCTOT           , ZQLAVTOT       , &
                 &    ZDQTDX(IPDF)     , ZDQSDX(IPDF)   , &
                 &    ZQLCCPDF(IPDF)   , ZQLAVPDF(IPDF) , &
                 &    ZXBAR(IPDF)      , ZXMIN(IPDF)    , ZXDIST(IPDF)     , &
                 &    ZSIGX(IPDF)      , ZQSX0          , ZQSX1            , ZNORM, &
                 &    XDIR             , ZFRC

INTEGER(KIND=JPIM) :: JK, JL, JP

LOGICAL ::            LLPBL(KLON,KLEV), LLSGSOVERLAP, LLCLOUDWAKE

REAL(KIND=JPRB) ::    ZQLWAKEFAC

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

! 
! #include "vdfpdftable.intfb.h"
! #include "vdfsat.intfb.h"
! #include "vdfthermo.intfb.h"



!     ------------------------------------------------------------------
!
!*         1.     INITIALIZE VARIABLES
!                 --------------------
!

! IF (LHOOK) CALL DR_HOOK('VDFCLOUD',0,ZHOOK_HANDLE)


!-- Cloud wake parameterization --
LLCLOUDWAKE = .TRUE.      !activation switch
!LLCLOUDWAKE = .FALSE.
!ZQLWAKEFAC = 0.1_JPRB    !proportionality factor with updraft PDF condensate
ZQLWAKEFAC = 0._JPRB

!-- SGS cloud overlap function --
LLSGSOVERLAP = .TRUE.      !activation switch
!LLSGSOVERLAP = .FALSE.
ZSGSEFOLD = 400._JPRB      !efolding depth  [m]


ZRG = 1._JPRB / RG

ZQMIN  = 0._JPRB
ZDUMMY = 0._JPRB
ZQBAR  = 0._JPRB
ZNORM  = 0._JPRB

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZSIGQT2UP(JL,JK) = 0._JPRB
    ZSIGQT2K(JL,JK)  = PSIGQT2(JL,JK)
  ENDDO
ENDDO
      
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    LLPBL(JL,JK) = KVARTOP(JL)/=0 .AND. JK>=KVARTOP(JL)
  ENDDO
ENDDO



!     ------------------------------------------------------------------
!
!*         2.     CLOSURE OF PDF MOMENTS
!                 ----------------------
!      
      
      
DO JL=KIDIA,KFDIA

  !----- fractions ---------
  ZATEST(JL) = PFRACB(JL,1)
  ZAUP(JL)   = PFRACB(JL,3)
  ZAK(JL)    = 1._JPRB - ZAUP(JL)
  
  !----- determine scaling factor between updraft variance and the
  !         difference in updraft qt excess of moist- and test updraft. -----
  IF ( ZAUP(JL)==0._JPRB .OR. ZAUP(JL).LT.1.5_JPRB*ZATEST(JL) ) THEN
    ZIQBAR(JL) = 0._JPRB  !change updraft pdf into a delta function
  ELSE
    ZFRC = ZATEST(JL)/ZAUP(JL)
    CALL VDFPDFTABLE(ZFRC, ZQBAR, ZQMIN, ZDUMMY, 0)
    ZIQBAR(JL) = 1._JPRB / ZQBAR    !iqbar=1/0.54 max, corresponding to ZATEST=1.5*ZAUP)
  ENDIF  
  
ENDDO
      
      
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

  IF (LLPBL(JL,JK) ) THEN

  !----- means ---------
  ZQTK(JL,JK)  = ( PQTM(JL,JK)  - ZAUP(JL) * PQTUP(JL,JK)  ) / ZAK(JL)
  ZSLGK(JL,JK) = ( PSLGM(JL,JK) - ZAUP(JL) * PSLGUP(JL,JK) ) / ZAK(JL)
        
  !----- variances -----
  ZSIGQT2UP(JL,JK) = (  ZIQBAR(JL) * MAX(0._JPRB,PQTTEST(JL,JK)-PQTUP(JL,JK))  )**2._JPRB
      
  ZSIGQT2K(JL,JK) = -(ZQTK(JL,JK)**2._JPRB) +&
                   & (               PQTM(JL,JK)**2._JPRB  + PSIGQT2(JL,JK)        &
                   &   - ZAUP(JL) * (PQTUP(JL,JK)**2._JPRB + ZSIGQT2UP(JL,JK))   ) &
                   & / ZAK(JL)
                   
  ZSIGQT2K(JL,JK) = MAX(0._JPRB,ZSIGQT2K(JL,JK))        
  
        
!  IF (LLDIAG) THEN
!    PEXTRA(JL,JK,81) = 1000._JPRB*PQTTEST(JL,JK)
!    PEXTRA(JL,JK,82) = 1000._JPRB*PQTUP(JL,JK)
!    PEXTRA(JL,JK,83) = 1000._JPRB*( PQTTEST(JL,JK)-PQTUP(JL,JK) )
!  
!    PEXTRA(JL,JK,84) = 1000000._JPRB*ZSIGQT2UP(JL,JK)
!    PEXTRA(JL,JK,85) = 1000000._JPRB*ZSIGQT2K(JL,JK)
!    PEXTRA(JL,JK,86) = 1000000._JPRB*PSIGQT2(JL,JK)
!  
!    PEXTRA(JL,JK,86:92) = 0._JPRB      
!    PEXTRA(JL,JK,79:80) = 0._JPRB      
!  ENDIF
  
  ENDIF

  ENDDO
ENDDO
       
      
      
      
!     ----------------------------------------------------------------
!
!*         3.     SUMMATION OF PDF-INTEGRALS
!                 --------------------------
!      
      
      
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
        
  !-----------------------------------------------------------------
  !   
  !  3.1  Reset
  !
  !-----------------------------------------------------------------
        
  ZCCPDF(:)   = (/0._JPRB, 0._JPRB/)
  ZQLCCPDF(:) = (/0._JPRB, 0._JPRB/)
  ZQLAVPDF(:) = (/0._JPRB, 0._JPRB/)
  ZXDIST(:)   = (/0._JPRB, 0._JPRB/)
        
  PCLDFRAC(JL,JK) = 0._JPRB
  PQLAV(JL,JK)    = 0._JPRB
        
  ZSMIN(1,:) = (/ 0._JPRB, 0._JPRB, 0._JPRB/)
  ZSMIN(2,:) = (/ 0._JPRB, 0._JPRB, 0._JPRB/)
        

  IF (LLPBL(JL,JK) ) THEN       !only within PBL
  
  !-----------------------------------------------------------------
  !   
  !  3.2  Initialization and vectorization
  !
  !       Formulated in {THL,QT} system
  !
  !       PDF index: 1=updraft, 2=diffusive
  !
  !-----------------------------------------------------------------
  
  !RCPDTMP = RCPD*(1._JPRB+RVTMP2*PQTM(JL,JK))
   
  ZAPDF      = (/ ZAUP(JL) , ZAK(JL) /)
  !ZAPDF      = (/ ZAUP(JL) , MAX(0._JPRB,ZAK(JL)-PAWAKE(JL,JK)) /)
  
  ZSIGQT2PDF = (/ MAX( 1000000._JPRB*ZSIGQT2UP(JL,JK),0._JPRB ), &
                & MAX( 1000000._JPRB*ZSIGQT2K(JL,JK) ,0._JPRB )  /)
        
  !-- PDF 2: K on dry zero buoyancy line --
  ZX0(2,:)    = (/ ZSLGK(JL,JK)/RCPD, 1000._JPRB * ZQTK(JL,JK)  , 0._JPRB/)
  ZX1(2,:)    = (/ 0._JPRB     , ZX0(2,2)+2.               , 0._JPRB/)  !-- plus 2 g/kg --
  ZX1(2,1)    = ZX0(2,1) * (1._JPRB + 0.00061_JPRB*ZX0(2,2)) / (1._JPRB + 0.00061_JPRB*ZX1(2,2))
  ZXMEAN(2,:) = ZX0(2,:)
        
  !-- PDF 1: updraft on mixing line --
  ZX1(1,:)    = (/ PSLGUP(JL,JK)/RCPD, 1000._JPRB * PQTUP(JL,JK), 0._JPRB /)
  ZX0(1,:)    = (/ PSLGM(JL,JK)/RCPD , 1000._JPRB * PQTM(JL,JK) , 0._JPRB /)
  ZXMEAN(1,:) = ZX1(1,:)
  ZNORM       = SQRT(DOT_PRODUCT(ZX1(1,:) - ZX0(1,:),ZX1(1,:) - ZX0(1,:)))
  IF (ZNORM.EQ.0.) THEN
    !-- Protection against zero length: assume mixing line = dry zero buoyancy line --
    ZX0(1,:)    = ZX0(2,:)
    ZX1(1,:)    = ZX1(2,:)
    ZXMEAN(1,:) = ZX0(2,:)
    !print '(a)',"zero length mixing line vector - reverting to zero buoy line"
  ENDIF  
        
  DO JP=1,IPDF
        
    !-- Unit vectors along PDF axes --
    ZNORM = SQRT(DOT_PRODUCT(ZX1(JP,:) - ZX0(JP,:),ZX1(JP,:) - ZX0(JP,:)))
    
    IF (ZNORM>0._JPRB) THEN
    
      ZNX(JP,:) = (ZX1(JP,:) - ZX0(JP,:))  / ZNORM
          
      !-- Project qt-variances onto PDF axes -- 
      IF (ABS(ZNX(JP,2)).GT.0.01_JPRB) THEN 
        ZSIGX(JP) = SQRT( (ZSIGQT2PDF(JP)) * (1._JPRB + (ZNX(JP,1)/ZNX(JP,2))**2._JPRB ) )
      ELSE
        ZSIGX(JP) = 0._JPRB
      ENDIF
    
    ELSE
    
      ZNX(JP,:) = (/0._JPRB, 1._JPRB, 0._JPRB /)
      ZSIGX(JP) = 0._JPRB
    
    ENDIF
          
  ENDDO
        
  !-- vector of linearized dry saturation curve --
  ZS0   = (/ PSLGM(JL,JK)/RCPD , 0._JPRB, 0._JPRB /)
  ZS1   = ZS0 + (/-1._JPRB,0._JPRB, 0._JPRB/)  !-- minus 1K --
  CALL VDFTHERMO( ZS1(1), 0._JPRB, ZS1(2), ZDUMMY, PAPM1(JL,JK), PGEOM1(JL,JK) )
  CALL VDFTHERMO( ZS0(1), 0._JPRB, ZS0(2), ZDUMMY, PAPM1(JL,JK), PGEOM1(JL,JK) )



  DO JP=1,IPDF
    
    !-----------------------------------------------------------------
    !   
    !  3.3  Retrieve the intersection point smin of the 
    !       PDF axis and the dry sat curve
    !
    !-----------------------------------------------------------------
          
    CALL VDFSAT(ZS1,ZS0,ZX1(JP,:),ZX0(JP,:),ZSMIN(JP,:),XDIR)  
    


    !-----------------------------------------------------------------
    !   
    !  3.4  Calculate the saturated fraction of the PDF
    !
    !-----------------------------------------------------------------
    
    IF (ZSMIN(JP,3)>-900._JPRB) THEN
      !-- Normalized saturation deficit --
      ZXDIST(JP) = SIGN(1._JPRB,XDIR) * DOT_PRODUCT(ZSMIN(JP,:)-ZXMEAN(JP,:),ZNX(JP,:))
      IF (ZSIGX(JP).GT.0._JPRB) THEN
        ZXMIN(JP) = ZXDIST(JP) / ZSIGX(JP)
      ELSE
        !-- Heaviside function --
        ZXMIN(JP) = SIGN(100._JPRB,ZXDIST(JP))  
      ENDIF  
    ELSE
      !-- Heaviside function --
      ZXDIST(JP) = ZSMIN(JP,2)-ZXMEAN(JP,2)
      ZXMIN(JP)  = SIGN(100._JPRB,ZXDIST(JP))  
    ENDIF  
    
    CALL VDFPDFTABLE(ZCCPDF(JP), ZXBAR(JP), ZXMIN(JP), ZDUMMY, 1)



    !---------------------------------------------------------------------
    !   
    !  3.5  Calculate total condensate of the saturated part of the PDF
    !
    !---------------------------------------------------------------------
    
    IF (ZCCPDF(JP).GT.0._JPRB) THEN
          
      IF (ZCCPDF(JP).LT.1._JPRB) THEN
      
        !-----------------------------------------------------------------
        !   
        !  3.5.1  Partial saturation:
        !
        !         Calculate PDF-average ql using gradients of qt and qsat 
        !         along the PDF axis (x)
        !
        !-----------------------------------------------------------------
      
        ZDQTDX(JP) = ZNX(JP,2)  ! Delta x = 1.
      
        !-- evaluate qs-gradient at SMIN --
        CALL VDFTHERMO(ZSMIN(JP,1),0._JPRB,ZQSX0,ZDUMMY,PAPM1(JL,JK), PGEOM1(JL,JK))
        CALL VDFTHERMO(ZSMIN(JP,1)+ZNX(JP,1),ZQSX0+ZNX(JP,2),ZQSX1,ZDUMMY,PAPM1(JL,JK), PGEOM1(JL,JK))
      
        ZDQSDX(JP) = ( ZQSX1-ZQSX0 ) ! Delta x = 1.
                
       
        !-----------------------------------------------------------------
        !   
        !  3.5.1.1  Distance of mean of moist part of PDF to smin,
        !           along the PDF axis
        !
        !-----------------------------------------------------------------
          
        ZXDIST(JP) = -ZXDIST(JP) + ZXBAR(JP)*ZSIGX(JP)
        ZXDIST(JP) = MAX(ZXDIST(JP),0._JPRB)
          
          
        !-----------------------------------------------------------------
        !   
        !  3.5.1.2  Integration of (qt-qs) over saturated part
        !
        !-----------------------------------------------------------------
          
        ZQLCCPDF(JP) = ZXDIST(JP)  * MAX(ZDQTDX(JP) - ZDQSDX(JP),0._JPRB)
      
      
      ELSE
          
        !-----------------------------------------------------------------
        !   
        !  3.5.2  Total saturation:
        !
        !         Don't do gradients along x: smin can be very far away
        !         from xmean, so that local linearization of qs no longer holds.
        !
        !         Instead, simply calculate qs of the PDF mean.
        !
        !-----------------------------------------------------------------
      
        CALL VDFTHERMO(ZXMEAN(JP,1),ZXMEAN(JP,2),ZQSX0,ZDUMMY,PAPM1(JL,JK),PGEOM1(JL,JK))
        
        ZQLCCPDF(JP) = MAX( ZXMEAN(JP,2)-ZQSX0, 0._JPRB )
        
        ZDQTDX(JP) = 0._JPRB
        ZDQSDX(JP) = 0._JPRB
        ZQSX1 = 0._JPRB
        
        
      ENDIF
      
      
      ZQLAVPDF(JP) = ZAPDF(JP) * ZCCPDF(JP) * ZQLCCPDF(JP)
          
          
    ENDIF
  
  
    !-----------------------------------------------------------------
    !   
    !  3.6  Cumulatives
    !
    !-----------------------------------------------------------------

    PCLDFRAC(JL,JK)  = PCLDFRAC(JL,JK)  + ZAPDF(JP)*ZCCPDF(JP)
    PQLAV(JL,JK)     = PQLAV(JL,JK)     + ZQLAVPDF(JP)

    !output: contribution of each PDF to..
    IF (LLDIAG) THEN
      PEXTRA(JL,JK,86+JP) = 100._JPRB*ZAPDF(JP)*ZCCPDF(JP)  !cloud fraction
      PEXTRA(JL,JK,89+JP) = ZQLAVPDF(JP)               !liquid water
    ENDIF  
        
    
  ENDDO !JP

  PQLAV(JL,JK) = MAX( 0._JPRB, PQLAV(JL,JK) / 1000._JPRB )    !from g/kg back to kg/kg
  
  !store contrib by diffusive PDF to qlav for use in LS precipitation generation routine (CLOUDSC)
  PLDIFF(JL,JK) = MAX( 0._JPRB, ZQLAVPDF(2) ) /1000._JPRB
  
  
    
  !-----------------------------------------------------------------
  !   
  !  3.6  Add cloud wake PDF
  !
  !-----------------------------------------------------------------
  
  IF (LLCLOUDWAKE) THEN
  
    !PCLDFRAC(JL,JK)  = MAX(0._JPRB, MIN( 1._JPRB,PCLDFRAC(JL,JK)+PAWAKE(JL,JK)) )
    PAWAKE(JL,JK)    = MAX(0._JPRB, MIN( PAWAKE(JL,JK), 1._JPRB-PCLDFRAC(JL,JK) ) )
    PCLDFRAC(JL,JK)  = PCLDFRAC(JL,JK)+PAWAKE(JL,JK)

    PQLAV(JL,JK)     = PQLAV(JL,JK) + ( ZQLCCPDF(1) * 0.001_JPRB * ZQLWAKEFAC ) * PAWAKE(JL,JK)
  
!    IF (LLDIAG) THEN
!      PEXTRA(JL,JK,89) = 100._JPRB * PAWAKE(JL,JK)
!      PEXTRA(JL,JK,92) = ( ZQLCCPDF(1) * 0.5_JPRB ) * PAWAKE(JL,JK)
!    ENDIF  
  
    PLDIFF(JL,JK) = PLDIFF(JL,JK) + ( ZQLCCPDF(1) * 0.001_JPRB * ZQLWAKEFAC ) * PAWAKE(JL,JK)
  
  ENDIF
  
 
  
  !-----------------------------------------------------------------
  !   
  !  4.0  Subgrid-scale cloud overlap function
  !
  !-----------------------------------------------------------------
  
  IF (LLDIAG) THEN
    PEXTRA(JL,JK,95) = PCLDFRAC(JL,JK)  !store cloud fraction before SGS overlap
  ENDIF  
    
  IF (LLSGSOVERLAP) THEN
    
    ZDZ = ( PGEOH(JL,JK-1) - PGEOH(JL,JK) ) * ZRG
    PCLDFRAC(JL,JK) = PCLDFRAC(JL,JK) * EXP( ZDZ / ZSGSEFOLD )
    PCLDFRAC(JL,JK) = MIN( 1._JPRB, PCLDFRAC(JL,JK) )
  
!    IF (LLDIAG) THEN
!      PEXTRA(JL,JK,93) = EXP( -ZDZ / ZSGSEFOLD )
!      PEXTRA(JL,JK,94) = ZDZ
!    ENDIF  
  
  ENDIF
  
  

  ENDIF !LLPBL


  ENDDO
ENDDO
         
        
! IF (LHOOK) CALL DR_HOOK('VDFCLOUD',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE VDFCLOUD
      
