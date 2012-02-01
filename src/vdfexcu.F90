SUBROUTINE VDFEXCU(KIDIA  , KFDIA  , KLON   , KLEV   , KDRAFT , PTMST  , PZ0MM  , &
                  &PHRLW  , PHRSW  , &
                  &PUM1   , PVM1   , PTM1   , PQM1   , PLM1   , PIM1   , &
                  &PAPHM1 , PAPM1  , PGEOM1 , PGEOH  , PCPTGZ , &
 ! DIAGNOSTIC OUTPUT
                  &PEXTR2 , KFLDX2 , PEXTRA , KLEVX  , KFLDX  , &
 !
                  &PKMFL  , PKHFL  , PKQFL  , &
                  &PCFM   , PCFH   , PTAUXCG, PTAUYCG, &
                  &PRICUI , PMCU   , PDTHV  , PMFLX  , KVARTOP, &
                  &PZINV  , KHPBL  , PKH    , PZCLDBASE , PZCLDTOP     , KPBLTYPE)
!     ------------------------------------------------------------------

!**   *VDFEXCU* - DETERMINES THE EXCHANGE COEFFICIENTS BETWEEN THE
!                 UPPER MODEL LEVELS WITH STABILITY AS A FUNCTION OF
!                 OBUKHOV-L

!     A.C.M. BELJAARS  26/03/90.  Original
!     A.C.M. BELJAARS  26/03/99   Tiling of the land surface.
!     J.HAGUE          13/01/2003 MASS Vector Functions
!     M. Ko"hler        3/12/2004 Moist Advection-Diffusion incl.
!                                 K,cloud and cloud top entrainment
!     P. Lopez         02/06/2005 Removed option for linearized
!                                 physics (now called separately)
!     R. Neggers       01/06/2006 Reorganization (into internal & interface K-modes)
!                                 Entrainment efficiency closure at cumulus PBL top
!     A. Beljaars      29/03/2006 Counter gradient stresses (Brown and Grant, 1997)
!
!     PURPOSE
!     -------

!     DETERMINE EXCHANGE COEFFICIENTS BETWEEN THE UPPER MODEL LEVELS

!     INTERFACE
!     ---------

!     *VDFEXCU* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1
!     *PTM1*         TEMPERATURE AT T-1
!     *PQM1*         SPECIFIC HUMUDITY AT T-1
!     *PAPHM1*       PRESSURE AT HALF LEVELS AT T-1
!     *PAPM1*        PRESSURE AT FULL LEVELS AT T-1
!     *PGEOM1*       GEOPOTENTIAL AT T-1
!     *PCPTGZ*       DRY STATIC ENERGY
!     *PKMFL*        KINEMATIC MOMENTUM FLUX                [#]
!     *PKHFL*        KINEMATIC HEAT FLUX                    [#]
!     *PKQFL*        KINEMATIC MOISTURE FLUX                [#]
!     *PZINV*        INVERSION HEIGHT              	    [M]
!     *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)

!     OUTPUT PARAMETERS (REAL):

!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT     (C-STAR IN DOC.)
!                    (ONLY PCFM(*,1:KLEV-1) AND
!                          PCFH(*,1:KLEV-1) ARE COMPUTED)
!     *PTAUXCG*      COUNTER GRADIENT STRESS X-COMPONENT    (N/m2)
!     *PTAUYCG*      COUNTER GRADIENT STRESS Y-COMPONENT    (N/m2)

!     REMARK: [#] UNUSED PARAMETERS IN TANGENT LINEAR AND ADJOINT VERSIONS
!     ------

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------
use garbage, only : phims, phihs, phimu, phihu
USE PARKIND1  ,ONLY : JPIM     ,JPRB
! ! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE yos_cst   , ONLY : RG       ,RD       ,RCPD     ,RETV     ,RATM
USE YOETHF   , ONLY : RVTMP2
USE YOEVDF   , ONLY : RLAM     ,RKAP     ,RVDIFTS  ,REPDU2   ,LLDIAG   
USE YOEVDFS  , ONLY : JPRITBL  ,RITBL    ,ARITBL   ,RCHBA    ,&
                    & RCHBB    ,RCHBD    ,RCHB23A  ,RCHBBCD  ,RCHBCD   ,&
                    & RCHETA   ,RCHETB   ,RCDHALF  ,RCDHPI2  ,RIMAX    ,&
                    & DRITBL   ,DRI26  
USE YOEPHLI  , ONLY : RLPMIXL  ,RLPBETA
USE YOMJFH   , ONLY : N_VMASS
use yos_exc, only : repust

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)                  :: PAPM1(KLON,KLEV) ! Argument NOT used
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTAUXCG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTAUYCG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZINV(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KHPBL(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZCLDBASE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZCLDTOP(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPBLTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRICUI(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTHV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMCU(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFLX(KLON,0:KLEV,KDRAFT) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVARTOP(KLON)
!          DIAGNOSTIC OUTPUT
INTEGER(KIND=JPIM),INTENT(IN)     :: KFLDX2, KLEVX, KFLDX
REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PEXTR2(KLON,KFLDX2), PEXTRA(KLON,KLEVX,KFLDX)

!*         0.2    LOCAL VARIABLES

REAL(KIND=JPRB) ::    ZRI(KLON),ZMGEOM(KLON),ZUST(KLON),&
                    & ZDTV(KLON),ZKHVFL(KLON),ZL(KLON),ZPHIM(KLON),&
                    & ZPHIH(KLON),ZTAUXCG(KLON),ZTAUYCG(KLON)    
REAL(KIND=JPRB) ::    ZDU2(KLON+N_VMASS)

INTEGER(KIND=JPIM) :: IRIB, JK, JL, JLEN
REAL(KIND=JPRB) ::    ZENTRSFC, ZENTRRAD, ZENTRTOP, &
                    & Z2GEOMF, ZA, ZALH2, ZALM2, ZB, ZCB, &
                    & ZCD, ZCFNC1, ZRHO, ZUABS, &
                    & ZCONS13, ZCONS1, ZHU1, ZHU2, ZWST3, ZRG, &
                    & ZDH, ZDL, ZDRORO, ZEPS, ZETA, ZHLM2, &
                    & ZLIM, ZLIM2, ZPHIKH, ZPHIKM, ZSCF, &
                    & ZX2, ZZ, ZWTVENTR, ZKH, ZCFHNEW, &
                    & ZML, ZBASE, ZVSC, ZKCLD, &
                    & ZKFACEDMF, ZTAUX, ZTAUY

REAL(KIND=JPRB) ::    ZZH, ZIFLTGM, ZIFLTGH, ZIFMOM, ZIFMOH, ZBM, ZBH, ZCM, ZCH, ZDUDZ
                    
LOGICAL ::         LLRICU
                    
REAL(KIND=JPRB) :: ZDRADFLX(KLON), ZRADKBASE(KLON), ZRADKDEPTH(KLON), &
                 & ZRADKFAC(KLON)

REAL(KIND=JPRB) :: ZWECUTOP(KLON)  , ZDTHVCUTOP, &
                 & ZTHVEN(KLON,KLEV), ZTHEN(KLON,KLEV), ZFAC

REAL(KIND=JPRB) ::    ZTMP1(KFDIA-KIDIA+1+N_VMASS)
REAL(KIND=JPRB) ::    ZTMP2(KFDIA-KIDIA+1+N_VMASS)
REAL(KIND=JPRB) ::    ZTMP3(KFDIA-KIDIA+1+N_VMASS)
REAL(KIND=JPRB) ::    ZTMP4(KFDIA-KIDIA+1+N_VMASS)
REAL(KIND=JPRB) ::    ZTMP5(KFDIA-KIDIA+1+N_VMASS)
REAL(KIND=JPRB) ::    ZHOOK_HANDLE
! 
! INTERFACE
! #include "surf_inq.h"
! END INTERFACE
! 
! #include "fcvdfs.h"



!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 --------------------

! IF (LHOOK) CALL DR_HOOK('VDFEXCU',0,ZHOOK_HANDLE)

ZENTRSFC  = 0.2_JPRB       ! factor for surface based top entrainment 
ZENTRRAD  = 0.2_JPRB       ! factor for radiative based top entrainment 
ZENTRTOP  = 0.4_JPRB       ! entrainment efficiency factor at cumulus PBL top, as proposed by Wyant et al (JAS, 1997)

ZCD       = 1.0_JPRB
ZCB       = 5.0_JPRB
ZEPS      = 1.E-10_JPRB

LLRICU = .TRUE.   ! switch for top-entrainment efficiency closure using Ri^cu at cumulus PBL top
!LLRICU = .FALSE.  

!ZKFACEDMF = 1.0_JPRB
ZKFACEDMF = 0.8_JPRB     !aup = 5%   !cy32r3
!ZKFACEDMF = 0.692_JPRB   !aup = 10%

! optimization
ZRG       = 1.0_JPRB/RG
ZCONS13   = 1.0_JPRB/3._JPRB
ZCONS1    = 0.5_JPRB*RKAP*ZRG/RLAM
ZHLM2     = 1.0_JPRB / ((2.0_JPRB*RG*RLPMIXL)**2)

IF(N_VMASS > 0) THEN
  JLEN=KFDIA-KIDIA+N_VMASS-MOD(KFDIA-KIDIA,N_VMASS)
ENDIF


!     ------------------------------------------------------------------

!*         2.     PREPARE SCALING COEFFICIENTS
!                 ----------------------------

  DO JL=KIDIA,KFDIA
    ZUST  (JL)=SQRT(MAX(PKMFL(JL),REPUST**2))
    ZKHVFL(JL)=PKHFL(JL)+RETV*PTM1(JL,KLEV)*PKQFL(JL)

    PTAUXCG(JL,KLEV)=0.0_JPRB
    PTAUYCG(JL,KLEV)=0.0_JPRB

    IF (ZKHVFL(JL)  <  0.0_JPRB) THEN
      ZWST3=-ZKHVFL(JL)*PZINV(JL)*RG/PTM1(JL,KLEV)
    ELSE
      ZWST3=0.0_JPRB
    ENDIF
    ZRHO =PAPHM1(JL,KLEV)/( RD*PTM1(JL,KLEV)*(1.0_JPRB+RETV*PQM1(JL,KLEV)) )
    ZUABS=MAX(0.1_JPRB,SQRT(PUM1(JL,KLEV)**2+PVM1(JL,KLEV)**2))
    ZHU1 =ZRHO*PKMFL(JL)/ZUABS
    ZTAUX=ZHU1*PUM1(JL,KLEV)
    ZTAUY=ZHU1*PVM1(JL,KLEV)

    ZHU2=2.7_JPRB*ZWST3/(ZUST(JL)**3+0.6_JPRB*ZWST3)
    ZTAUXCG(JL)=ZTAUX*ZHU2
    ZTAUYCG(JL)=ZTAUY*ZHU2
  ENDDO
  
  !calculate full level mean theta_v profile for later use
  DO JK=KLEV,1,-1
    DO JL=KIDIA,KFDIA
      ZTHEN(JL,JK)   = ( PAPM1(JL,JK)/RATM )**(-RD/RCPD) * PTM1(JL,JK)
      ZTHVEN(JL,JK)  = ZTHEN(JL,JK) * &
                     & ( 1.0_JPRB + RETV * PQM1(JL,JK)   - PLM1(JL,JK)   - PIM1(JL,JK)   )
    ENDDO
  ENDDO

!  IF (LLDIAG) THEN
!    DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
!      IF (JK>=KVARTOP(JL)-1) THEN
!        PEXTRA(JL,JK,4) = ZTHVEN(JL,JK) - ZTHVEN(JL,KLEV)
!        PEXTRA(JL,JK,5) = ZTHEN(JL,JK)  - ZTHEN(JL,KLEV)
!      ENDIF  
!    ENDDO
!    ENDDO
!  ENDIF  
  

!          Calculate PBL cloud top radiative flux jump [Km/s] (cooling)
!          for top-driven K and entrainment velocity formulations.

  DO JL=KIDIA,KFDIA
    SELECT CASE (KPBLTYPE(JL))
      CASE(2)
        ZRADKDEPTH(JL) = PZINV(JL)                               
        ZRADKBASE(JL)  = 0._JPRB    
        ZRADKFAC(JL)   = 1._JPRB                        
      CASE(3)
        ZRADKDEPTH(JL) = MAX( PZCLDTOP(JL) - PZCLDBASE(JL),  0._JPRB )
        ZRADKBASE(JL)  = PZCLDTOP(JL) - ZRADKDEPTH(JL)
        ZRADKFAC(JL)   = 0.25_JPRB                        
      CASE DEFAULT
        ZRADKDEPTH(JL) = -100._JPRB
        ZRADKBASE(JL)  = -100._JPRB
        ZRADKFAC(JL)   = 0._JPRB                        
    END SELECT 
  ENDDO


  ZDRADFLX(:) = 0.0_JPRB
  DO JK=KLEV-1,1,-1
    DO JL=KIDIA,KFDIA

!      IF ( KPBLTYPE(JL) == 2 ) THEN  !RN for now, only stratocumulus: extension to cumulus planned (for intermediate scenarios like ATEX)

        IF ( PGEOH(JL,JK)*ZRG <= ZRADKBASE(JL)+ZRADKDEPTH(JL) .AND. ZRADKBASE(JL)+ZRADKDEPTH(JL) < PGEOH(JL,JK-1)*ZRG ) THEN
!          ZDRADFLX(JL) = -  PHRLW(JL,JK+1)                 * (PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZRG   !use LW divergence only
          ZDRADFLX(JL) = - (PHRLW(JL,JK+1)+PHRSW(JL,JK+1)) * (PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZRG
!         ... add solar heating at 2nd cld level
          IF ( PZCLDBASE(JL) < PGEOH(JL,JK+1)*ZRG .AND. JK < KLEV-1 ) THEN 
            ZDRADFLX(JL) = ZDRADFLX(JL) - PHRSW(JL,JK+2)   * (PGEOH(JL,JK+1)-PGEOH(JL,JK+2))*ZRG
          ENDIF
          ZDRADFLX(JL) = MAX( ZDRADFLX(JL), 0.0_JPRB )    !safety against rad. heating cases
        ENDIF

!      ENDIF

    ENDDO
  ENDDO



!     ------------------------------------------------------------------

!*         3.     VERTICAL LOOP - non-linear physics
!                 ----------------------------------

  IF(N_VMASS > 0) THEN
    IF(KFDIA-KIDIA+1 /= JLEN) THEN
      ZTMP1(KFDIA-KIDIA+2:JLEN)=1.0_JPRB 
      ZTMP2(KFDIA-KIDIA+2:JLEN)=1.0_JPRB 
      ZTMP3(KFDIA-KIDIA+2:JLEN)=1.0_JPRB 
      ZDU2(KFDIA+1:KIDIA+JLEN-1)=1.0_JPRB 
    ENDIF
  ENDIF


!***
  DO JK=KLEV-1,1,-1
!***

    DO JL=KIDIA,KFDIA
      PCFM(JL,JK)=0.0_JPRB
      PCFH(JL,JK)=0.0_JPRB
      PKH(JL,JK) =0.0_JPRB
      PTAUXCG(JL,JK)=0.0_JPRB
      PTAUYCG(JL,JK)=0.0_JPRB
    ENDDO

    IF(N_VMASS <= 0) THEN   ! efficiency of exponentials

!          COMPUTE RI-NUMBER

      DO JL=KIDIA,KFDIA
        ZDU2(JL)=MAX(REPDU2,(PUM1(JL,JK)-PUM1(JL,JK+1))**2&
                         & +(PVM1(JL,JK)-PVM1(JL,JK+1))**2)  
        ZDRORO= 2.0_JPRB * (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
         & / ( PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)&
         &   - PGEOM1(JL,JK)-PGEOM1(JL,JK+1))&
         & - (RVTMP2-RETV)*(PQM1(JL,JK)-PQM1(JL,JK+1))
        ZDTV(JL)=( (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
         & - (RVTMP2-RETV)*0.5_JPRB * (PQM1(JL,JK)-PQM1(JL,JK+1))&
         & * (PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)) ) * (1.0_JPRB/RCPD)  
        ZMGEOM(JL)=PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
        ZRI(JL)=ZMGEOM(JL)*ZDRORO/ZDU2(JL)
      ENDDO

    ELSE

      DO JL=KIDIA,KFDIA
        ZDU2(JL)=MAX(REPDU2,(PUM1(JL,JK)-PUM1(JL,JK+1))**2&
                         & +(PVM1(JL,JK)-PVM1(JL,JK+1))**2)
        ZTMP2(JL-KIDIA+1)= PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)&
                        & -PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
        ZDTV(JL)=( (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
         & -(RVTMP2-RETV)*0.5_JPRB* (PQM1(JL,JK)-PQM1(JL,JK+1))&
         & *(PCPTGZ(JL,JK)+PCPTGZ(JL,JK+1)) ) * (1.0_JPRB/RCPD)
        ZMGEOM(JL)=PGEOM1(JL,JK)-PGEOM1(JL,JK+1)
      ENDDO
! 
!       CALL VREC(ZTMP4,ZDU2(KIDIA),JLEN)
!       CALL VREC(ZTMP5,ZTMP2,JLEN)

      DO JL=KIDIA,KFDIA
        ZDRORO= 2.0_JPRB * (PCPTGZ(JL,JK)-PCPTGZ(JL,JK+1))&
         & *ZTMP5(JL-KIDIA+1)&
         & - (RVTMP2-RETV)*(PQM1(JL,JK)-PQM1(JL,JK+1))  
        ZRI(JL)=ZMGEOM(JL)*ZDRORO*ZTMP4(JL-KIDIA+1)
      ENDDO

    ENDIF

    DO JL=KIDIA,KFDIA

!          COMPUTE STABILITY FUNCTIONS

      IF (ZRI(JL)  >  0.0_JPRB) THEN

!        INTERPLOLATE ETA WITH SPLINES FOR POSITIVE RICHARDSON
!        NUMBERS

        IRIB=INT(ZRI(JL)*(1.0_JPRB/DRITBL))+1
        
        IF (IRIB  >=  JPRITBL) THEN
!           LINEAR EXTENSION OF LOOK-UP TABLE
          ZETA = RITBL(JPRITBL)*(ZRI(JL)*(1.0_JPRB/RIMAX))
        ELSE
          ZX2  = IRIB*DRITBL
          ZA   = (ZX2-ZRI(JL))*(1.0_JPRB/DRITBL)
          ZB   = 1.0_JPRB-ZA
          ZETA = ZA*RITBL(IRIB) + ZB*RITBL(IRIB+1)&
           & +( (ZA**3-ZA)*ARITBL(IRIB)&
           & +(  ZB**3-ZB)*ARITBL(IRIB+1) )*DRI26  
        ENDIF

!        STABLE PHI-FUNCTIONS

        ZPHIM(JL) = PHIMS(ZETA)
        ZPHIH(JL) = PHIHS(ZETA)
      ELSE

!        UNSTABLE SITUATIONS

        ZETA  = ZRI(JL)
        ZPHIM(JL) = PHIMU(ZETA)
        ZPHIH(JL) = PHIHU(ZETA)
      ENDIF
    ENDDO

!   IF(N_VMASS <= 0) THEN ! Vector MASS taken out because completely changed code

    DO JL=KIDIA,KFDIA

!-------- up to CY32R3 --------
!
!!          COMMON FACTORS FOR STABLE AND UNSTABLE
!
!     Z2GEOMF=PGEOM1(JL,JK)+PGEOM1(JL,JK+1)+2.0_JPRB*RG*PZ0MM(JL)
!     ZLIM=RLPBETA + (1.0_JPRB-RLPBETA)/(1.0_JPRB+Z2GEOMF*Z2GEOMF*ZHLM2)
!     ZLIM2=ZLIM*ZLIM
!     ZALM2=ZLIM2*(0.5_JPRB*RKAP*ZRG*Z2GEOMF/(1.0_JPRB+ZCONS1*Z2GEOMF))**2
!     ZALH2=ZALM2
!     ZCFNC1=RVDIFTS*PTMST*RG**2 * PAPHM1(JL,JK)&
!      & /( 0.5_JPRB*RD * ZMGEOM(JL)&
!      & *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  ))&
!      & +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1))))  
!     ZCFNC=RG*ZCFNC1*SQRT(ZDU2(JL))/ZMGEOM(JL)
!
!!          DIMENSIONLESS COEFFICIENTS MULTIPLIED BY PRESSURE
!!          THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE.
!
!     IF (ZRI(JL)  >  0.0_JPRB) THEN  ! statically stable
!       ZSCF=SQRT(1.0_JPRB+ZCD*ZRI(JL))
!       PCFM(JL,JK)=ZCFNC*ZALM2/(1.0_JPRB+2.0_JPRB*ZCB*ZRI(JL)/ZSCF)
!       PCFH(JL,JK)=ZCFNC*ZALH2/(1.0_JPRB+2.0_JPRB*ZCB*ZRI(JL)*ZSCF)
!     ELSE                            ! statically unstable
!       PCFM(JL,JK)=ZCFNC*ZALM2/(ZPHIM(JL)**2)
!       PCFH(JL,JK)=ZCFNC*ZALM2/(ZPHIM(JL)*ZPHIH(JL))
!     ENDIF
!-----------------------------

!-------- CY32R3 -------------
! new K: K,LTG scales with l=kappa*z 
!        K,MO  scales with l=150m above surface layer

!          COMMON FACTORS FOR STABLE AND UNSTABLE

      ZIFMOM  = 1.0_JPRB / (ZPHIM(JL)**2)                              !F(MO),M
      ZIFMOH  = 1.0_JPRB / (ZPHIM(JL)*ZPHIH(JL))                       !F(MO),H
      ZDUDZ   = SQRT(ZDU2(JL))/ZMGEOM(JL)*RG                           !shear
      ZCFNC1=RVDIFTS*PTMST*RG**2 * PAPHM1(JL,JK)&                      !factor for vdfdifh/m
       & /( 0.5_JPRB*RD * ZMGEOM(JL)&
       & *( PTM1(JL,JK  )*(1.0_JPRB+RETV*PQM1(JL,JK  ))&
       & +  PTM1(JL,JK+1)*(1.0_JPRB+RETV*PQM1(JL,JK+1))))  

!          DIMENSIONLESS COEFFICIENTS MULTIPLIED BY PRESSURE
!          THICKNESSES FOR MOMENTUM AND HEAT EXCHANGE.

      IF ( ZRI(JL) > 0.0_JPRB ) THEN  ! statically stable
        ZZH     = 0.5_JPRB * ZRG * (PGEOM1(JL,JK)+PGEOM1(JL,JK+1)) + PZ0MM(JL)
        ZSCF    = SQRT(1.0_JPRB+ZCD*ZRI(JL))
        ZIFLTGM = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)/ZSCF) !F(LTG),M
        ZIFLTGH = 1.0_JPRB / (1.0_JPRB + 2.0_JPRB *ZCB * ZRI(JL)*ZSCF) !F(LTG),H
        ZBM     = RKAP * ZZH * SQRT(ZIFLTGM) 
        ZBH     = RKAP * ZZH * SQRT(ZIFLTGH) 
        ZCM     = 150.0_JPRB * SQRT(ZIFMOM )
        ZCH     = 150.0_JPRB * SQRT(ZIFMOH )

        PCFM(JL,JK) = ZCFNC1 * ZDUDZ * (ZBM*ZCM/(ZBM+ZCM))**2
        PCFH(JL,JK) = ZCFNC1 * ZDUDZ * (ZBH*ZCH/(ZBH+ZCH))**2
      ELSE                            ! statically unstable
        PCFM(JL,JK) = ZCFNC1 * ZDUDZ * 150.0_JPRB**2 * ZIFMOM
        PCFH(JL,JK) = ZCFNC1 * ZDUDZ * 150.0_JPRB**2 * ZIFMOH
      ENDIF
!-----------------------------


      !  overwrite Ri-based K values within boundary layer
      !IF ( JK >= KVARTOP(JL) .AND. KVARTOP(JL)>0 ) THEN 
      !  PCFH(JL,JK) = 0._JPRB
      !  PCFM(JL,JK) = 0._JPRB
      !ENDIF

      
      
      !------------------------------------------------------------------
      !
      !   3.1   INTERNAL K MODE OF MIXED LAYER
      !         ------------------------------
      !
      !         Using a prescribed vertical structure.
      !
      
      ZZ     = PGEOH(JL,JK)*ZRG

      !IF (PGEOH(JL,JK-1)*ZRG <= PZINV(JL) .AND. ZKHVFL(JL)<0._JPRB ) THEN  ! up to level below entr. level
      IF ( JK > KHPBL(JL) .AND. ZKHVFL(JL)<0._JPRB ) THEN 

        ZL(JL) = ZUST (JL)**3*PTM1(JL,KLEV)/(RKAP*RG*(ZKHVFL(JL)-ZEPS))

        ZDH    = ZZ/PZINV(JL)
        ZDL    = ZZ/ZL(JL)
        ZETA   = ZZ/ZL(JL)
        ZPHIKH = (1-39._JPRB*ZDL)**(-ZCONS13)
        ZPHIKM = (1-15._JPRB*ZDL)**(-ZCONS13)

        !   K,surface
        
        PCFH(JL,JK)  = ZKFACEDMF * ZCFNC1 * RKAP / ZPHIKH * ZUST(JL) * ZZ * (1.0_JPRB-ZDH)**2
        PCFM(JL,JK)  = ZKFACEDMF * ZCFNC1 * RKAP / ZPHIKM * ZUST(JL) * ZZ * (1.0_JPRB-ZDH)**2

        PTAUXCG(JL,JK)=ZDH*(1._JPRB-ZDH)**2*ZTAUXCG(JL)
        PTAUYCG(JL,JK)=ZDH*(1._JPRB-ZDH)**2*ZTAUYCG(JL)

      ENDIF
      
      
      
      !------------------------------------------------------------------
      !
      !   3.1   INTERNAL K MODE IN CUMULUS CLOUD LAYER
      !         ------------------------------
      !
      
      IF ( KPBLTYPE(JL) == 3 .OR. KPBLTYPE(JL) == 4) THEN   

        IF ( JK < KHPBL(JL) .AND. JK >= KVARTOP(JL) ) THEN
          IF ( JK > KVARTOP(JL) ) THEN 
            !  K diffusion within cumulus layer
            !PCFH(JL,JK)  = PCFH(JL,JK)          !testing: Ri diffusion 
            !PCFM(JL,JK)  = PCFM(JL,JK)
            PCFH(JL,JK)  = PCFH(JL,KHPBL(JL))   !cy32r3
            PCFM(JL,JK)  = PCFM(JL,KHPBL(JL))
          ELSE
            !  Reset K 
            PCFH(JL,JK)  = 0._JPRB 
            PCFM(JL,JK)  = 0._JPRB
          ENDIF
        ENDIF
        
        !  Protect top-entrainment of shallow cu topped mixed layers against a zero
        !    buoyancy jump (dthv) through the cloud base transition layer.
        !    This can easily occur just before h grows one layer,due to dq
        !    cancelling ds in dthv.    -RN
        !
        ZDTV(JL) = MAX( ZDTV(JL), 0.2_JPRB ) 
        
      ENDIF



      !------------------------------------------------------------------
      !
      !   3.2    INTERNAL K MODE DRIVEN BY CLOUD-TOP COOLING
      !          -------------------------------------------
      !
      !          As in Lock et al. (2000, MWR p3187f), equ. 5:
      !          Using simplified radiative velocity scale as in Lock, 1998, equ. 12
      !                 and ignore buoyancy reversal velocity scale
      !

      IF ( ZKHVFL(JL)<0._JPRB .AND. ZRADKDEPTH(JL)>0._JPRB) THEN
        IF ( ZZ >= ZRADKBASE(JL)  .AND.  ZZ <= ZRADKBASE(JL)+ZRADKDEPTH(JL) ) THEN  
          ZVSC  = ( RG / PTM1(JL,JK) * ZRADKDEPTH(JL) * ZDRADFLX(JL) ) ** ZCONS13 
          !ZKCLD = 0.85_JPRB * RKAP * ZVSC &
          ZKCLD = ZRADKFAC(JL) * 0.85_JPRB * RKAP * ZVSC &
              & * (ZZ-ZRADKBASE(JL)) ** 2 / ZRADKDEPTH(JL) &
              & * ( 1 - (ZZ-ZRADKBASE(JL)) / ZRADKDEPTH(JL) ) ** 0.5_JPRB 
          IF (KPBLTYPE(JL)==2) THEN     
            PCFH(JL,JK)  = PCFH(JL,JK) + ZCFNC1 * ZKCLD
            PCFM(JL,JK)  = PCFM(JL,JK) + ZCFNC1 * ZKCLD * 0.75_JPRB
          ENDIF  
        ENDIF
      ENDIF



      !------------------------------------------------------------------
      !
      !   3.3   INTERFACE K AT CUMULUS PBL TOP
      !         ------------------------------
      !
      !         At top level, mass flux is replaced by diffusion, using an entrainment efficiency formulation.
      !         This is important for representing the intermediate regime (StCu->Cu transitions)    -RN
      !

      IF ( LLRICU .AND. JK == KVARTOP(JL) .AND. KPBLTYPE(JL) == 3) THEN 

        !entrainment efficiency - after Wyant et al. (JAS, 1997)
        !ZWECUTOP(JL) = PMCU(JL) * ZENTRTOP * PRICUI(JL)
        ZWECUTOP(JL) = 2._JPRB * PMCU(JL) * ZENTRTOP * PRICUI(JL)
        ZWECUTOP(JL) = MAX(0.0_JPRB,ZWECUTOP(JL))
          
        !  translation into K [m2/s] at this level: K = entrainment velocity * mixing-length (dz)
        ZKH = ZWECUTOP(JL)                                !top-entrainment by overshooting surface-driven thermals
        ZKH = ZKH + ZENTRRAD * ZDRADFLX(JL) / ZDTV(JL )   !add cloud top cooling driven entrainment
        ZKH = ZKH * ZMGEOM(JL) * ZRG
        
        ZKH     = MAX(0.0_JPRB,ZKH)
        ZCFHNEW = ZCFNC1 * ZKH
            
        PCFH(JL,JK) = ZCFHNEW
        PCFM(JL,JK) = ZCFHNEW * 0.75_JPRB
            
        !  reset any updraft M
        PMFLX(JL,JK,2) = 0._JPRB
        PMFLX(JL,JK,3) = 0._JPRB

      ENDIF



      !---------------------------------------------------------------------
      !
      !   3.4   INTERFACE K AT MIXED LAYER TOP 
      !         ------------------------------
      ! 
      !    This is cloud top for dry & stratocu PBL, and cloud base in 
      !    shallow cu PBL (PZINV).
      !
      !    w_e is the mixed layer top-entrainment rate, defined as
      !
      !          w_e = w'thv'_h / dthv_h = -0.2 w'thv'_s / dthv_h = A/Ri w_*.
      !
      !    This w_e is here translated into K at this level.
      !

      !IF ( PGEOH(JL,JK)*ZRG <= PZINV(JL)  .AND.  PZINV(JL) < PGEOH(JL,JK-1)*ZRG ) THEN
      IF ( JK == KHPBL(JL) ) THEN 
        
        !wthv_h = -0.2 wthv_s
        ZWTVENTR = -ZENTRSFC * ZKHVFL(JL) * 0.1_JPRB
        IF (KPBLTYPE(JL) == 2) THEN
          ZWTVENTR = -ZENTRSFC * ZKHVFL(JL)
        ENDIF  

        !--- Special stratocumulus treatment: radiation impacts on entrainment
        !---
        !--- ENTRAINMENT VELOCITY * T,v JUMP DUE TO LW-RADIATIVE COOLING & SFC FLUX
        !---    (Lock & MacVean, 1999, equ. 11)
        IF ( KPBLTYPE(JL) == 2 ) THEN                     
          ZWTVENTR = ZWTVENTR + ZENTRRAD * ZDRADFLX(JL) !radiation flux jump
        ENDIF

        
        ZWTVENTR = MAX(0.0_JPRB,ZWTVENTR)

        !ZDTV(JL) = 1.0_JPRB
        !ZDTV(JL) = 10.*ZDTV(JL)
        !ZDTV(JL) = 2._JPRB * PDTHV(JL)

        ZKH     = ZWTVENTR * ZMGEOM(JL) / ( RG * ZDTV(JL) )

        ZKH     = MAX(0.0_JPRB,ZKH)
        ZCFHNEW = ZCFNC1 * ZKH

        !RN PCFH(JL,JK)  = MAX(PCFH(JL,JK),ZCFHNEW)           !protection against K=0
        !RN PCFM(JL,JK)  = MAX(PCFM(JL,JK),ZCFHNEW * 0.75_JPRB)
        PCFH(JL,JK)  = ZCFHNEW           
        PCFM(JL,JK)  = ZCFHNEW * 0.75_JPRB

!        IF (LLDIAG) THEN
!          PEXTR2(JL,4) = - ZDRADFLX(JL)         ! radiative flux jump         [K m/s]
!          PEXTR2(JL,5) = ZWTVENTR               ! entrainment flux = we * dTv [K m/s]
!          PEXTR2(JL,7) = ZKH / ZMGEOM(JL)*RG    ! we                          [m/s]
!          PEXTR2(JL,6) = ZDTV(JL)               ! dTv                         [K]
!          PEXTR2(JL,7) = ZDTV(JL) * RG / ZMGEOM(JL)
!        ENDIF
          
      ENDIF



!          DIFFUSION COEFFICIENT FOR HEAT FOR POSTPROCESSING ONLY IN (M2/S)

      PKH(JL,JK) = PCFH(JL,JK) / ZCFNC1


      !RN output
      IF (LLDIAG) THEN
!        PEXTR2(JL,19) = ZWECUTOP(JL)           !top entrainment rate
!        PEXTR2(JL,20) = ZENTRTOP * PRICUI(JL)  !entrainment efficiency
        PEXTRA(JL,JK,3) = PCFH(JL,JK)   / ZCFNC1   !K [m2/s]
      ENDIF  
      
    ENDDO

!***
  ENDDO !JK
!***


! IF (LHOOK) CALL DR_HOOK('VDFEXCU',1,ZHOOK_HANDLE)
END SUBROUTINE VDFEXCU
