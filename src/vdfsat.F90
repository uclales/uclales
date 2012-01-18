!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFSAT (s1,s0,x1,x0,smin,xdir)
                   
                     
!     ------------------------------------------------------------------

!**   *VDFSAT* - DETERMINES A SATURATION SPECIFIC HUMIDITY

!     Roel Neggers     02/11/2005  


!     PURPOSE     
!     -------     

!     DETERMINE THE INTERSECTION POINT OF THE LINEARIZED DRY SATURATION CURVE
!     AND A JOINT-PDF AXIS.
!                

!     INTERFACE
!     ---------

!     *VDFSAT* IS CALLED BY *VDFCLOUD*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):


!     INPUT PARAMETERS (REAL):


!     INPUT PARAMETERS (LOGICAL):


!     OUTPUT PARAMETERS (REAL):


!     METHOD
!     ------

      !--------- Vector routine to retrieve qtsatmin ------------------------!
      !                                                                      !
      !  This routine returns the intersection point smin of the two         !
      !  vectors s1-s0 (the linearized dry saturation curve) and x1-x0       !
      !  (the linearized mixing line).                                       !
      !                                                                      !
      !  The vectors are defined in {thl,qt} space. Use is made of a vector  !
      !  calculus method: the mixing line vector is projected onto a new     !
      !  system defined by the axes parallel and perpendicular to the        !
      !  linearized saturation curve.                                        !
      !                                                                      !
      !  The benefit of this vector method over a gradient method is that it !
      !  avoids awkward gradients (e.g. dqt/dthl with dthl=0).               !
      !                                                                      !
      !                                                                      !
      !  Roel Neggers, ECMWF    22-02-2006                                   !
      !                                                                      !
      !----------------------------------------------------------------------!


!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

REAL(KIND=JPRB)   ,INTENT(IN)    :: S1(3),S0(3),X1(3),X0(3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: SMIN(3),XDIR


!*         0.2    LOCAL VARIABLES
      
REAL(KIND=JPRB) ::  n1(3),n2(3),n3(3),nsat(3),nsatperp(3),nmeantosat(3),nmix(3)

REAL(KIND=JPRB) ::  norm,pA,pB,pC,pD,dist

! REAL(KIND=JPRB) ::  ZHOOK_HANDLE


 

!     ------------------------------------------------------------------

!*         1.     
!                 -------

! IF (LHOOK) CALL DR_HOOK('VDFSAT',0,ZHOOK_HANDLE)


smin(:) = (/0._JPRB,0._JPRB,0._JPRB/)
  
  
!-- vector perpendicular to dry saturation curve --
nsat = s1 - s0
norm = SQRT(DOT_PRODUCT(nsat,nsat))
nsat = nsat / norm     !normalize
        
        
!-- vector coming out of {thl,qt} plane --
n3 = (/0._JPRB, 0._JPRB, 1._JPRB/)
        
        
!-- vector perpendicular on dry saturation vector (cross product with n3) --
nsatperp =  (/ nsat(2) * n3(3) - nsat(3) * n3(2), &
         & -nsat(1) * n3(3) + nsat(3) * n3(1), &
         &  nsat(1) * n3(2) - nsat(2) * n3(1) /)


!-- mixing line vector --
nmix = x1 - x0
        
        
!-- vector connecting mean to dry saturation curve --
nmeantosat = s0 - x0
        
        
!-- do some projections --
pA =  DOT_PRODUCT(nmix,nsatperp)
pB =  DOT_PRODUCT(nmeantosat,nsatperp)
pC = -DOT_PRODUCT(nmeantosat,nsat)
pD =  DOT_PRODUCT(nmix,nsat)

xdir = pA        


IF ( ABS(pD)>100._JPRB*ABS(pA) ) THEN

  !-- lines are almost parallel: their angle is < 0.6 degrees (=invtan(0.01)) --
  smin = s0
  
  !smin(:) = (/-999._JPRB,-999._JPRB,-999._JPRB/)
  smin(3) = -999._JPRB   !send warning back to VDFCLOUD that no intersection point is found
  
  !print '(a,6f)','vdfsat warning: ',ABS(pD),ABS(pA),s1(1)-s0(1),s1(2)-s0(2),x1(1)-x0(1),x1(2)-x0(2)
  
ELSE

  !-- there is a reasonable intersection point --
  dist = (pC + pB * pD/pA)
  smin = s0  + dist * nsat

ENDIF  



! IF (LHOOK) CALL DR_HOOK('VDFSAT',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE VDFSAT
