!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFPDFTABLE (PFRAC    , PFACEXC    , PFACX     , PQLFAC    , IFLAG)  
!     ------------------------------------------------------------------

!**   *VDFPDFTABLE* - CALCULATES GAUSSIAN PDF SCALING FACTORS
!
!     Roel Neggers   ECMWF     21/04/2005 


!     PURPOSE
!     -------

!     CALCULATES THE SCALING FACTOR BETWEEN THE MEAN VALUE OF A CERTAIN 
!     TOP-PERCENTAGE OF A GAUSSIAN PDF (WITH ZERO MEAN) AND ITS STANDARD DEVIATION:
!
!         PHI(PFRAC)    = PFACEXC(PFRAC) * SIGMA_PHI,
!         PHIMIN(PFRAC) = PFACX(PFRAC)   * SIGMA_PHI,
!
!     WHERE PHI IS AN ARBITRARY VARIABLE, PFRAC IS THE PERCENTAGE, AND
!     SIGMA_PHI IS THE STANDARD DEVIATION. PHIMIN IS THE MINIMUM PHI VALUE
!     OF THE BIN.
!
!     THREE SEARCH OPTIONS (IFLAG):                                        
!                                                                
!       IFLAG=0: PFACEXC and PFACX are returned as a function of PFRAC.
!       IFLAG=1: PFACEXC and PFRAC are returned as a function of PFACX.
!       IFLAG=2: PFACEXC, PFACX and PFRAC are returned as a function of PQLFAC.
!

!     INTERFACE
!     ---------

!     *VDFPDFTABLE* IS CALLED BY *VDFHGHTN*
!                                *VDFCLOUD*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (REAL):

!     *PFRAC*       THE TOP PERCENTAGE OF THE PDF                 %

!     OUTPUT PARAMETERS (REAL):

!     *PFACEXC*     THE FACTOR ASSOC WITH THE MEAN                 -
!     *PFACX*       THE FACTOR ASSOC WITH THE PDF CUT-OFF POINT X  -
!     *PQLFAC*      THE RATIO  (QSAT-QTAV)/QLAV  =  
!                              PFACX / (PFRAC *(PFAXEXC-PFACX) )


!     METHOD
!     ------

!     INTEGRATION OF PART OF THE GAUSSIAN PDF: SEE DOCUMENTATION

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRAC
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFACEXC
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFACX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLFAC
INTEGER(KIND=JPIM),INTENT(IN)    :: IFLAG


!*         0.2    LOCAL VARIABLES

INTEGER(KIND=JPIM) :: IDIM, JI, IA
PARAMETER (IDIM=40)

REAL(KIND=JPRB) ::    ZA(IDIM), ZBAR(IDIM), ZX(IDIM), ZQLFAC(IDIM), ZFAC

REAL(KIND=JPRB) ::    ZHOOK_HANDLE

LOGICAL ::  LLSTOP


!     -----------------------------------------------------------------

!*         1.     SET SOME CONSTANTS
!                 --------------------

! IF (LHOOK) CALL DR_HOOK('VDFPDFTABLE',0,ZHOOK_HANDLE)

!array of percentages
ZA(1) = 0.0001_JPRB
ZA(2) = 0.001_JPRB
!--- 0 to 10% step 0.5% ---
DO JI=3,22
  ZA(JI) = (JI-2)*0.005_JPRB
ENDDO
!--- 10 to 20% step 1% ---
DO JI=23,32
  ZA(JI) = 0.1_JPRB + (JI-22)*0.01_JPRB
ENDDO
!--- 20 to 100% step 10% ---
DO JI=33,IDIM-1
  ZA(JI) = (JI-30)*0.1
ENDDO
ZA(IDIM) = 0.999_JPRB

!array of associated factors
ZBAR = (/ 3.9581_JPRB , &   !0.01
       &  3.3678_JPRB , &   !0.1
       &  2.8982_JPRB , &
       &  2.6733_JPRB , &   !1
       &  2.5253_JPRB , &
       &  2.4247_JPRB , &   !2
       &  2.3376_JPRB , &
       &  2.2672_JPRB , &   !3
       &  2.2074_JPRB , &
       &  2.1532_JPRB , &   !4
       &  2.1078_JPRB , &
       &  2.0622_JPRB , &   !5
       &  2.0232_JPRB , &
       &  1.9847_JPRB , &   !6
       &  1.9510_JPRB , &
       &  1.9175_JPRB , &   !7
       &  1.8861_JPRB , &
       &  1.8585_JPRB , &   !8
       &  1.8312_JPRB , &
       &  1.8037_JPRB , &   !9
       &  1.7782_JPRB , &
       &  1.7540_JPRB , &   !10
       &  1.7091795_JPRB , &
       &  1.6667867_JPRB , &
       &  1.6270676_JPRB , &
       &  1.5896436_JPRB , &
       &  1.5542130_JPRB , &
       &  1.5205319_JPRB , &
       &  1.4883975_JPRB , &
       &  1.4576342_JPRB , &
       &  1.4281014_JPRB , &
       &  1.3996768_JPRB , &   !20
       &  1.1588855_JPRB , &   !30
       &  0.9657842_JPRB , &   !40
       &  0.7978202_JPRB , &   !50
       &  0.6438460_JPRB , &   !60
       &  0.4966645_JPRB , &   !70
       &  0.3499109_JPRB , &   !80
       &  0.1949553_JPRB , &   !90
       &  0.0033274_JPRB /)    !top fraction &    99.9% : whole pdf

ZX = (/ 3.7190497_JPRB , &   !0.01
      & 3.0902362_JPRB , &   !0.1
      & 2.5758295_JPRB , &
      & 2.3263483_JPRB , &   !1
      & 2.1700907_JPRB , &
      & 2.0537486_JPRB , &   !2
      & 1.9599638_JPRB , &
      & 1.8807937_JPRB , &   !3
      & 1.8119106_JPRB , &
      & 1.7506862_JPRB , &   !4
      & 1.6953980_JPRB , &
      & 1.6448534_JPRB , &   !5
      & 1.5981932_JPRB , &
      & 1.5547736_JPRB , &   !6
      & 1.5141019_JPRB , &
      & 1.4757911_JPRB , &   !7
      & 1.4396021_JPRB , &
      & 1.4050716_JPRB , &   !8
      & 1.3722037_JPRB , &
      & 1.3407550_JPRB , &   !9
      & 1.3105791_JPRB , &
      & 1.2815515_JPRB , &   !10
      & 1.2265279_JPRB , &
      & 1.1749867_JPRB , &
      & 1.1263912_JPRB , &
      & 1.0803194_JPRB , &
      & 1.0364335_JPRB , &
      & 0.9944579_JPRB , &
      & 0.9541652_JPRB , &
      & 0.9153650_JPRB , &
      & 0.8778962_JPRB , &
      & 0.8416213_JPRB , &   !20
      & 0.5244006_JPRB , &   !30
      & 0.2533472_JPRB , &   !40
      & 0.0000000_JPRB , &   !50
      & -0.2533471_JPRB , &   !60
      & -0.5244005_JPRB , &   !70
      & -0.8416210_JPRB , &   !80
      & -1.2815515_JPRB , &   !90
      & -3.0902183_JPRB /)    !99.9
      
ZQLFAC = (/ 201451.4062500_JPRB, &
      & 11606.1923828_JPRB, &
      & 1645.2435303_JPRB, &
      & 689.9847412_JPRB, &
      & 409.3978882_JPRB, &
      & 280.4224854_JPRB, &
      & 207.9357910_JPRB, &
      & 162.1739960_JPRB, &
      & 130.9946594_JPRB, &
      & 108.5721817_JPRB, &
      & 91.7832108_JPRB, &
      & 78.8131485_JPRB, &
      & 68.5385056_JPRB, &
      & 60.2317429_JPRB, &
      & 53.4003906_JPRB, &
      & 47.7006798_JPRB, &
      & 42.8894539_JPRB, &
      & 38.7755508_JPRB, &
      & 35.2328568_JPRB, &
      & 32.1543732_JPRB, &
      & 29.4590988_JPRB, &
      & 27.0839195_JPRB, &
      & 23.1020756_JPRB, &
      & 19.9096317_JPRB, &
      & 17.3056850_JPRB, &
      & 15.1506004_JPRB, &
      & 13.3445892_JPRB, &
      & 11.8146152_JPRB, &
      & 10.5061712_JPRB, &
      & 9.3779278_JPRB, &
      & 8.3977880_JPRB, &
      & 7.5406590_JPRB, &
      & 2.7549937_JPRB, &
      & 0.8890161_JPRB, &
      & 0.0000000_JPRB, &
      & -0.4706290_JPRB, &
      & -0.7336884_JPRB, &
      & -0.8829191_JPRB, &
      & -0.9644020_JPRB, &
      & -0.9999243_JPRB/) 
      
!     -----------------------------------------------------------------

!*         2.     INTERPOLATE
!                 -----------

SELECT CASE (IFLAG)


  CASE(0) 
      
      !--- PFRAC is input variable ---
      IF ( ZA(1)>=PFRAC ) THEN

        !PDF is totally unsaturated
        PFACEXC = ZBAR(1)   !set an upper limit to the pdf
        PFACX   = ZX(1)

      ELSE

        LLSTOP = .FALSE.
  
        DO IA = 1,IDIM-1
          IF ( ZA(IA)<PFRAC .AND. ZA(IA+1)>=PFRAC .AND. .NOT. LLSTOP ) THEN
    
          ZFAC = ( PFRAC - ZA(IA) ) / ( ZA(IA+1) - ZA(IA) )
      
          PFACEXC = ZBAR(IA) + ZFAC * ( ZBAR(IA+1) - ZBAR(IA) )
          PFACX   = ZX(IA)   + ZFAC * ( ZX(IA+1) - ZX(IA) )
      
          LLSTOP = .TRUE.
    
          ENDIF
        ENDDO
  
        IF (.NOT. LLSTOP) THEN
          PFACEXC = 0.0_JPRB    !out of bounds: Fraction is larger than 99.9%, so set PFACEXC to 0.0
          PFACX   = -99.9_JPRB    !lower limit can now be anything
        ENDIF

      ENDIF


  CASE(1) 
      
      !--- PFACX is input variable ---
      IF ( ZX(1)<=PFACX ) THEN
            
        !PDF is totally unsaturated
        PFACEXC = PFACX      !=ZBAR(1)
        PFRAC   = 0._JPRB    !=ZA(1)   

      ELSE

        LLSTOP = .FALSE.
  
        DO IA = 1,IDIM-1
          IF ( ZX(IA)>PFACX .AND. ZX(IA+1)<=PFACX .AND. .NOT. LLSTOP ) THEN
    
          ZFAC = ( PFACX - ZX(IA) ) / ( ZX(IA+1) - ZX(IA) )
      
          PFACEXC = ZBAR(IA) + ZFAC * ( ZBAR(IA+1) - ZBAR(IA) )
          PFRAC   = ZA(IA)   + ZFAC * ( ZA(IA+1) - ZA(IA) )
      
          LLSTOP = .TRUE.
          
          ENDIF
        ENDDO
  
        IF (.NOT. LLSTOP) THEN
          !PDF is totally saturated
          PFACEXC = 0.0_JPRB    !set mean to 0.0
          PFRAC   = 1.0_JPRB    
        ENDIF

      ENDIF
      
      
  CASE(2) 
      
      !--- PQLFAC is input variable ---
      IF ( ZQLFAC(1)<=PQLFAC ) THEN
          
        !PDF is totally unsaturated
        PFACEXC = ZBAR(1)
        PFACX   = PFACEXC
        PFRAC   = 0._JPRB   

      ELSE

        LLSTOP = .FALSE.
  
        DO IA = 1,IDIM-1
          IF ( ZQLFAC(IA)>PQLFAC .AND. ZQLFAC(IA+1)<=PQLFAC .AND. .NOT. LLSTOP ) THEN
    
          ZFAC = ( PQLFAC - ZQLFAC(IA) ) / ( ZQLFAC(IA+1) - ZQLFAC(IA) )
      
          PFACEXC = ZBAR(IA) + ZFAC * ( ZBAR(IA+1) - ZBAR(IA) )
          PFACX   = ZX(IA)   + ZFAC * ( ZX(IA+1) - ZX(IA) )
          PFRAC   = ZA(IA)   + ZFAC * ( ZA(IA+1) - ZA(IA) )
      
          LLSTOP = .TRUE.
          
          ENDIF
        ENDDO
  
        IF (.NOT. LLSTOP) THEN
          !PDF is totally saturated
          PFACEXC = 0.0_JPRB    !set mean to 0.0
          PFACX   = -99.9_JPRB  !lower limit can now be anything
          PFRAC   = 1.0_JPRB    
        ENDIF

    ENDIF
      
      
END SELECT !IFLAG



! IF (LHOOK) CALL DR_HOOK('VDFPDFTABLE',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE VDFPDFTABLE
