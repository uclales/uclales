module parkind1
 implicit none
 integer, parameter :: jprb = selected_real_kind(p=13,r=307) !< the kind of a 64bit real
  integer, parameter :: jpim = selected_int_kind( 9)          !< the kind of a 32bit integer


end module parkind1

module yos_cst

use parkind1  ,only : jpim     ,jprb
  
implicit none

save
real(kind=jprb),parameter :: zrclum=299792458._jprb
real(kind=jprb),parameter :: zrhpla=6.6260755e-34_jprb
real(kind=jprb),parameter :: zrkbol=1.380658e-23_jprb
real(kind=jprb),parameter :: zrnavo=6.0221367e+23_jprb
real(kind=jprb),parameter :: zrmd=28.9644_jprb
real(kind=jprb),parameter :: zrmv=18.0153_jprb
real(kind=jprb),parameter :: rtt =273.16_jprb    ! rt0=273.16= triple point temperature
real(kind=jprb),parameter :: rpi = 2.0_jprb*asin(1.0_jprb)     ! pi=3.14...
real(kind=jprb),parameter :: rday=86400._jprb   ! astronomical constants
real(kind=jprb),parameter :: r =zrnavo*zrkbol      ! thermodynamic gas phase
real(kind=jprb),parameter :: rd =1000._jprb*r/zrmd     ! thermodynamic gas phase
real(kind=jprb),parameter :: rv =1000._jprb*r/zrmv     ! thermodynamic gas phase
real(kind=jprb),parameter :: retv  =rv/rd-1.0_jprb  ! thermodynamic gas phase
real(kind=jprb),parameter :: rlstt =2.8345e+6_jprb  ! thermodynamic transition of phase
real(kind=jprb),parameter :: rlvtt =2.5008e+6_jprb  ! thermodynamic transition of phase
real(kind=jprb),parameter :: rlmlt =rlstt-rlvtt  ! thermodynamic transition of phase
real(kind=jprb),parameter :: rcpd  =3.5_jprb*rd  ! thermodynamic gas phase
real(kind=jprb),parameter :: rcpv  =4._jprb *rv  ! thermodynamic gas phase
real(kind=jprb),parameter :: rsigma =2.0_jprb * rpi**5 * zrkbol**4 /(15._jprb* zrclum**2 * zrhpla**3) ! radiation
real(kind=jprb),parameter :: rg    =9.80665_jprb  ! geoide
real(kind=jprb),parameter :: ratm =100000._jprb   ! 

end module yos_cst

module yomct0
  use parkind1, only : jprb
implicit none
  logical, parameter :: lscmec = .true., lsfcflx = .true.
  real(kind=jprb)    :: rextshf = 1, rextlhf = 1
end module yomct0

module yomjfh
  use parkind1, only : jpim
implicit none
  integer(kind=jpim), parameter :: n_vmass = 0
end module yomjfh

module yoethf
  use yos_cst, only : rd, rv, rtt, rlvtt, rcpd, rlstt
  use parkind1, only : jprb
implicit none
  real(kind=jprb), parameter :: r2es = 611.21_jprb*rd/rv, r3les = 17.502_jprb, r4les = 32.19_jprb, r3ies = 22.587_jprb, r4ies = -0.7_jprb, r5ies = r3ies*(rtt-r4ies), r5les = r3les*(rtt-r4les), r5alvcp = 1, r5alscp = 1
  real(kind=jprb), parameter :: ralvdcp = rlvtt/rcpd, ralsdcp = rlstt/rcpd, rtwat = rtt, rtice = rtt-23._jprb, rticecu = rtt-23._jprb, rtwat_rtice_r = 1.0_jprb/(rtwat-rtice), rtwat_rticecu_r =1.0_jprb/(rtwat-rticecu)
  real(kind=jprb), parameter :: rvtmp2 = 0._jprb
end module yoethf
module yoecumf
  use parkind1, only : jprb
implicit none
  real(kind=jprb), parameter :: rtaumel=5._jprb*3.6e3_jprb
end module yoecumf

module garbage
  use parkind1, only : jprb, jpim
  use yoethf, only : r2es, r3ies, r3les, r4ies, r4les, rtt, ralsdcp, ralvdcp, r5alscp, r5alvcp, rtice, rtwat, rtwat_rtice_r
  implicit none
contains
  real(kind=jprb) function foedelta (ptare)
  real(kind=jprb) :: ptare
      foedelta = max (0.0_jprb,sign(1.0_jprb,ptare-rtt))
 end function foedelta
  real(kind=jprb) function foealfa (ptare)
  real(kind=jprb) :: ptare
foealfa = min(1.0_jprb,((max(rtice,min(rtwat,ptare))-rtice)&
 &*rtwat_rtice_r)**2) 
 end function foealfa
  real(kind=jprb) function foedem (ptare)
  real(kind=jprb) :: ptare
      foedem = foealfa(ptare)*r5alvcp*(1.0_jprb/(ptare-r4les)**2)+&
             &(1.0_jprb-foealfa(ptare))*r5alscp*(1.0_jprb/(ptare-r4ies)**2)
 end function foedem
  real(kind=jprb) function foeldcpm (ptare)
  real(kind=jprb) :: ptare
foeldcpm = foealfa(ptare)*ralvdcp+&
            &(1.0_jprb-foealfa(ptare))*ralsdcp
 end function foeldcpm
 
 real(kind=jprb) function foeewm (ptare)
  real(kind=jprb) :: ptare
    foeewm = r2es *&
     &(foealfa(ptare)*exp(r3les*(ptare-rtt)/(ptare-r4les))+&
  &(1.0_jprb-foealfa(ptare))*exp(r3ies*(ptare-rtt)/(ptare-r4ies)))
 end function foeewm

 subroutine cuadjtq &
 & (kidia,    kfdia,    klon,     ktdia,    klev,&
 & kk,&
 & psp,      pt,       pq,       ldflag,   kcall)  

!          m.tiedtke         e.c.m.w.f.     12/89

!          modifications
!          -------------
!          d.salmond         cray(uk))      12/8/91
!          j.j. morcrette    ecmwf          92-09-18   update to cy44
!          j.f. mahfouf      ecmwf          96-06-11   smoothing option
!          d.salmond & m.hamrud ecmwf       99-06-04   optimisation
!          j.hague                          03-01-13   mass vector functions
!          j.hague                          03-07-07   more mass v.f.
!        m.hamrud              01-oct-2003 cy28 cleaning
!        j.hague & d.salmond   22-nov-2005 optimisations 

!          purpose.
!          --------
!          to produce t,q and l values for cloud ascent

!          interface
!          ---------
!          this routine is called from subroutines:
!              *cond*     (t and q at condensation level)
!              *cubase*   (t and q at condensation level)
!              *cuasc*    (t and q at cloud levels)
!              *cuini*    (environmental t and qs values at half levels)
!              *custrat*  (t and q at condensation level)
!          input are unadjusted t and q values,
!          it returns adjusted values of t and q

!     parameter     description                                   units
!     ---------     -----------                                   -----
!     input parameters (integer):

!    *kidia*        start point
!    *kfdia*        end point
!    *klon*         number of grid points per packet
!    *ktdia*        start of the vertical loop
!    *klev*         number of levels
!    *kk*           level
!    *kcall*        defines calculation as
!                      kcall=0  env. t and qs in*cuini*
!                      kcall=1  condensation in updrafts  (e.g. cubase, cuasc)
!                      kcall=2  evaporation in downdrafts (e.g. cudlfs,cuddraf)

!     input parameters (logical):

!    *ldland*       land-sea mask (.true. for land points)

!     input parameters (real):

!    *psp*          pressure                                        pa

!     updated parameters (real):

!    *pt*           temperature                                     k
!    *pq*           specific humidity                             kg/kg

!          externals   
!          ---------
!          3 lookup tables ( tlucua, tlucub, tlucuc )
!          for condensation calculations.
!          the tables are initialised in *suphec*.

!----------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
! ! use yomhook   ,only : lhook,   dr_hook
! 
! use yos_cst   , only : retv     ,rlvtt    ,rlstt    ,rtt
! use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les    ,&
!  & r4ies    ,r5les    ,r5ies    ,r5alvcp  ,r5alscp  ,&
!  & ralvdcp  ,ralsdcp  ,rtwat    ,rtice    ,rticecu  ,&
!  & rtwat_rtice_r      ,rtwat_rticecu_r  
! use yoephli  , only : lphylin  ,rlptrc   ,rlpal1   ,rlpal2
use yomjfh   , only : n_vmass
use yos_cst, only : retv

implicit none

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim)               :: ktdia ! argument not used
integer(kind=jpim),intent(in)    :: kk 
real(kind=jprb)   ,intent(in)    :: psp(klon) 
real(kind=jprb)   ,intent(inout) :: pt(klon,klev) 
real(kind=jprb)   ,intent(inout) :: pq(klon,klev) 
logical           ,intent(in)    :: ldflag(klon) 
integer(kind=jpim),intent(in)    :: kcall 
integer(kind=jpim) :: jl, jlen

real(kind=jprb) :: ztmp0(kfdia-kidia+1+n_vmass)
real(kind=jprb) :: ztmp1(kfdia-kidia+1+n_vmass)
real(kind=jprb) :: ztmp2(kfdia-kidia+1+n_vmass)
real(kind=jprb) :: ztmp3(kfdia-kidia+1+n_vmass)
real(kind=jprb) :: ztmp4(kfdia-kidia+1+n_vmass)
real(kind=jprb) :: ztmp5(kfdia-kidia+1+n_vmass)
real(kind=jprb) :: ztmp6(kfdia-kidia+1+n_vmass)

real(kind=jprb) :: z1s, z2s, zcond,zcond1, zcor, zfoeewi, zfoeewl,&
 & zoealfa, zqmax, zqsat, ztarg, zqp

!dir$ vfunction exphf
! #include "fcttre.h"

!     statement functions
!real_b :: foealfaj,foedemj,foeldcpmj,foeewmj

real(kind=jprb) :: minj, maxj, x, y
! real(kind=jprb) :: zhook_handle

minj(x,y) = y - 0.5_jprb*(abs(x-y)-(x-y))
maxj(x,y) = y + 0.5_jprb*(abs(x-y)+(x-y))

!----------------------------------------------------------------------

!     1.           define constants
!                  ----------------

zqmax=0.5_jprb


!dir$    ivdep
!ocl novrec
    do jl=kidia,kfdia
      if(ldflag(jl)) then
        zqp    =1.0_jprb/psp(jl)
        zqsat=foeewm(pt(jl,kk))*zqp    
        zqsat=min(0.5_jprb,zqsat)
        zcor=1.0_jprb/(1.0_jprb-retv  *zqsat)
        zqsat=zqsat*zcor
        zcond=(pq(jl,kk)-zqsat)/(1.0_jprb+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond
        pq(jl,kk)=pq(jl,kk)-zcond
        zqsat=foeewm(pt(jl,kk))*zqp    
        zqsat=min(0.5_jprb,zqsat)
        zcor=1.0_jprb/(1.0_jprb-retv  *zqsat)
        zqsat=zqsat*zcor
        zcond1=(pq(jl,kk)-zqsat)/(1.0_jprb+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk)=pt(jl,kk)+foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk)=pq(jl,kk)-zcond1
      endif
    enddo

end subroutine cuadjtq

subroutine satadj(pt, psp, pqt, ql)
  use thrm, only : rslf
  real, intent(in)  :: psp, pqt
  real, intent(inout) :: pt
  real, intent(out)   :: ql
  real :: qs, T, Told
  integer :: n
    T = pt
    Told = 0
    n = 0
    do while (abs(T - Told)/Told > 1e-5)
      n = n + 1
      Told = T
      qs = rslf(psp, T)
      ql=max(0._jprb,pqt-qs)
      T = 0.5*(pt +ralvdcp*ql+Told)
! print *, n, t  , told ,  abs(T - Told)/Told , ql
    end do
! print *, pt, T, pqt, psp, ql    
    pt = T
!     stop
end subroutine satadj

end module garbage

module yoevdfs

use parkind1  ,only : jpim     ,jprb

implicit none

save

!     ------------------------------------------------------------------
!*    ** *yoevdfs* contains stability function tables for *vdf...*
!     ------------------------------------------------------------------
logical :: firsttime = .true.
integer(kind=jpim), parameter :: jpritbl=101
real(kind=jprb) :: ritbl(jpritbl)
real(kind=jprb) :: aritbl(jpritbl)
real(kind=jprb) :: rchba
real(kind=jprb) :: rchbb
real(kind=jprb) :: rchbc
real(kind=jprb) :: rchbd
real(kind=jprb) :: rchb23a
real(kind=jprb) :: rchbbcd
real(kind=jprb) :: rchbcd
real(kind=jprb) :: rcheta
real(kind=jprb) :: rchetb
real(kind=jprb) :: rchetc
real(kind=jprb) :: rchbhdl
real(kind=jprb) :: rcdhalf
real(kind=jprb) :: rcdhpi2
real(kind=jprb) :: rimax
real(kind=jprb) :: dritbl
real(kind=jprb) :: dri26

!**   ** *yoedfs* contains stability function tables for *vdf...*

!     a.c.m. beljaars   e.c.m.w.f.       26/03/90.

!      name      type        purpose
!      ----      ----        -------

!     *rchba*     real       *constant a in *holtslag and *debruin
!                            functions for stable situations
!     *rchbb*     real       *constant b in *hb* functions
!     *rchbc*     real       *constant c in *hb* functions
!     *rchbd*     real       *constant d in *hb* functions
!     *rchb23a    real       2./3.*a in *hb* functions
!     *rchbbcd    real       b*c/d in *hb* functions
!     *rchbcd     real       c/d in *hb* functions
!     *rcheta*    real       constant in the *hogstrom *ellison *turner
!                            functions for stably stratified turbulence 
!     *rchetb*    real       constant in the *het* functions     
!     *rchetc*    real       constant in *het* functions  
!     *rchbhdl*   real       maxim znlev/l for stable boundary layer  
!     *rcdhalf    real       constant in *dyer and *hicks formulae
!                            for unstable situations
!     *rcdhpi2    real       pi/2.
!     *rimax*     real       *maximim richardson number tabulated
!     *dritbl*    real       *increment of the richardson number
!                            between tabulated values.
!     *dri26*     real       dritbl**2/6.
!     *ritbl*     real array *tabulated eta-values (z/l) as a function
!                            of the richardson number for stable cases.
!     *aritbl*    real array *second derivatives of tabulated function
!                            for spline interpolation.
!     ------------------------------------------------------------------
contains
  subroutine suvdfs

  !**   *subroutine* *suvdfs* initializes common block *yoevdfs*

  !      a. beljaars   e.c.m.w.f.   26/03/90

  !      purpose
  !      -------

  !           *subroutine *suvdfs* initializes the constants needed by
  !      the empirical stability functions and sets up the table that
  !      gives the stability parameter eta (height devided by
  !      *obukhov length) as a function of the gradient *richardson number
  !      for stable situations. also the second derivatives are
  !      tabulated for later use by spline interpolation.
  !           *inivdfs* checks wether the *phi* and *psi* expressions are
  !      consistent (by means of numerical differentiation). if the
  !      functions are not consistent, the routine aborts.

  !      interface
  !      ---------

  !           *call* *suvdfs* from *suphec*

  !      method
  !      ------

  !           *the algebraic equation to be solved is
  !      ri=(phih/phim**2)*eta, where *phih* and *phim* are the gradient
  !      stability functions for heat and momentum respectively.
  !           *to solve the implicit algebraic equation for *eta* as a
  !      function of *ri*, *newton's method is employed with a fixed
  !      number of iterations. the necessary derivatives are evaluated
  !      numerically.
  !           *after completion of the function table, the second
  !      derivatives are derived for later use by the spline
  !      interpolation.

  !      externals
  !      ---------

  !           *statement functions *phims*, *phihs*, *psims* and *psihs*

  !      reference
  !      ---------

  !           *see *press et al. (1986; numerical recipes - the art of
  !      scientific computing) for details on the spline interpolation.

  !      -----------------------------------------------------------------

  use parkind1  ,only : jpim     ,jprb
  real(kind=jprb) :: zu(jpritbl)

  integer(kind=jpim) :: itmax, jit, jjp

  real(kind=jprb) :: zdrv1, zdrvn, zeps, zeta, zetam, zetap, zfxm,&
  & zfxp, zp, zqn, zrib, zun  
  real(kind=jprb) :: zhook_handle



  !*     1. initialize constants in common block
  !         ---------- --------- -- ------ -----

  !     1.1 constants related to eta(ri)-table

  rimax   = 10._jprb
  dritbl  = rimax/(jpritbl-1)
  dri26   = dritbl**2/6._jprb
  ritbl(1)= 0._jprb

  !     1.2 constants for the ellison turner functions (stable)

  rcheta=5._jprb
  rchetb=4._jprb
  rchetc=8._jprb
  rchbhdl = 5._jprb

  !     1.21 constants for holtslag and debruin functions (stable psi)

  rchba   = 1._jprb
  rchbb   = 2.0_jprb/3._jprb
  rchbc   = 5._jprb
  rchbd   = 0.35_jprb
  rchb23a = (2.0_jprb/3._jprb)*rchba
  rchbbcd = rchbb*rchbc/rchbd
  rchbcd  = rchbc/rchbd

  !     1.3 constants for dyer and hicks expressions (unstable)

  rcdhalf = 16._jprb
  rcdhpi2 = 2.0_jprb*atan(1.0_jprb)

  !*    3. loop over table index
  !        ---- ---- ----- -----

  do jjp=2, jpritbl
    zrib=(jjp-1)*dritbl
  !     3.1 initial guess

    zeta=ritbl(jjp-1)
    if (zrib  <  0.5_jprb) then
      itmax=5
    else
      itmax=3
    endif

  !     3.2 newton's iteration loop with derivative from finite difference

    do jit=1, itmax
      zeps=(zeta+1.0_jprb)*0.001_jprb
      zetap=zeta+zeps
      zetam=zeta-zeps
      zfxp=zrib-zetap*phihs(zetap)/phims(zetap)**2
      zfxm=zrib-zetam*phihs(zetam)/phims(zetam)**2
      zeta=zeta-zeps*(zfxp+zfxm)/(zfxp-zfxm)
    enddo
  !     3.3 store result

    ritbl(jjp)=zeta
  enddo

  !     4. compute second derivatives from tabulated results
  !        ------- ------ ----------- ---- --------- -------

  !     4.1 derivative at ri=0. and ri=rimax

  zdrv1=1._jprb
  zdrvn=(ritbl(jpritbl)-ritbl(jpritbl-1))/dritbl

  aritbl(1)=-0.5_jprb
  zu(1)=(3._jprb/dritbl)*((ritbl(2)-ritbl(1))/dritbl-zdrv1)

  !     4.2. decomposition loop of tridiagonal matrix

  do jjp=2,jpritbl-1
    zp=0.5_jprb*aritbl(jjp-1)+2.0_jprb
    aritbl(jjp)=-0.5_jprb/zp
    zu(jjp)=(6._jprb*((ritbl(jjp+1)-ritbl(jjp))/dritbl &
    & -(ritbl(jjp)-ritbl(jjp-1))&
    & /dritbl)/(2*dritbl)-0.5_jprb*zu(jjp-1))/zp  
  enddo

  zqn=0.5_jprb
  zun=(3._jprb/dritbl)*(zdrvn-(ritbl(jpritbl)-ritbl(jpritbl-1))/dritbl)
  aritbl(jpritbl)=(zun-zqn*zu(jpritbl-1))/(zqn*aritbl(jpritbl-1)+1.0_jprb)

  !     4.3 back substitution

  do jjp=jpritbl-1,1,-1
    aritbl(jjp)=aritbl(jjp)*aritbl(jjp+1)+zu(jjp)
  enddo

  end subroutine suvdfs

  real function phihs(peta)
    real(kind=jprb) :: peta
    phihs=(1.0_jprb+rchetb*peta)**2
  end function phihs
  
  real function phims(peta)
    real(kind=jprb) :: peta
    phims=1.0_jprb+rcheta*peta
  end function phims
  
  real function phihu(peta)
    real(kind=jprb) :: peta
    phihu= 1.0_jprb/     sqrt(1.0_jprb-rcdhalf*peta)
  end function phihu
  
  real function phimu(peta)
    real(kind=jprb) :: peta
  phimu= 1.0_jprb/sqrt(sqrt(1.0_jprb-rcdhalf*peta))
  end function phimu
  
end module yoevdfs
module yoevdf

use parkind1  ,only : jpim     ,jprb

implicit none

save

!     ------------------------------------------------------------------
!*    ** *yoevdf* contains constants needed by *vdf....*
!     ------------------------------------------------------------------

! integer(kind=jpim) ,parameter :: nvtypes
real(kind=jprb) ,parameter :: rlam =150._jprb
real(kind=jprb) ,parameter :: rkap=0.4_jprb
real(kind=jprb) ,parameter :: rvdifts=1.5_jprb
real(kind=jprb) ,parameter :: repdu2=(0.1_jprb)**2
! real(kind=jprb) ,parameter :: rentr=0.20_jprb
! real(kind=jprb) ,parameter :: rpar=2._jprb
! real(kind=jprb) ,parameter :: rpar1=0.6_jprb
! real(kind=jprb) ,parameter :: rparsrf=0.1_jprb
! 
logical :: lldiag = .true.

!*     *yoevdf* contains constants needed by *vdf....*
!     for the computation of vertical diffusion

!     a.c.m. beljaars      e.c.m.w.f.    14/12/89

!     obukhov-l update     acmb          26/03/90.   

!     name        type     description
!     ----        ----     -----------

!     *nvtypes*   integer  number of vegetation (surface cover) types
!     *rlam*      real     *asymptotic mixing length for momentum
!     *rkap*      real     *vonkarman constant
!     *rvdifts*   real     *factor for time step weighting in *vdf....*
!     *repdu2*    real     *minimum velocity difference in ri-number
!     *rentr*     real     *entrainment constant          
!     *rpar*      real     *parameter for temperature excess in thermal 
!                           at boundary layer top      
!     *rpar1*     real     *coefficient of (w*)**3 in ws         
!     *rparsrf*   real     *depth of surface layer as fraction of pbl-h 
!     *lldiag*    logical  *switch for vdf extra variables output
!     ------------------------------------------------------------------
end module yoevdf

module yoephli

use parkind1  ,only : jpim     ,jprb
use yoethf, only : rtice, rtwat
implicit none

save

!     ------------------------------------------------------------------
!*    ** *yoephli* contains constants for the linearized physics
!     ------------------------------------------------------------------

logical :: lphylin = .false.
logical :: lenopert = .false.
logical :: leppcfls = .false.
logical :: lraisanen = .false.

real(kind=jprb) ,parameter :: rlptrc=rtice+(rtwat-rtice)/sqrt(2.0_jprb)
real(kind=jprb) ,parameter :: rlpal1=0.15_jprb
real(kind=jprb) ,parameter :: rlpal2=20._jprb
real(kind=jprb) ,parameter :: rlpbb  =5._jprb
real(kind=jprb) ,parameter :: rlpcc  =5._jprb
real(kind=jprb) ,parameter :: rlpdd  =5._jprb
real(kind=jprb) ,parameter :: rlpmixl=4000._jprb
real(kind=jprb) ,parameter :: rlpbeta=0.2_jprb
real(kind=jprb) ,parameter :: rlpdrag=0._jprb
real(kind=jprb) ,parameter :: rlpevap=0.0_jprb
real(kind=jprb) ,parameter :: rlpp00=30000._jprb

!*     *yoephli* contains constants needed by 
!     the linearized physics

!     j.f. mahfouf        e.c.m.w.f.    23/06/96

!     name        type     description
!     ----        ----     -----------

!     *rlptrc*    real     critical temperature for mixed phase properties
!                          of water 
!     *rlpal1*    real     smoothing coefficient
!     *rlpal2*    real     smoothing coefficient
!     *rlpbb*     real     constant from the louis et al. formulation
!     *rlpcc*     real     constant from the louis et al. formulation
!     *rlpdd*     real     constant from the louis et al. formulation
!     *rlpmixl*   real     pseudo depth of the planetary boundary layer
!     *rlpbeta*   real     reduction factor of the asymptotic mixing length
!     *rlpdrag*   real     coefficient for the estimation of surface drag
!     *rlpevap*   real     fraction of possible rainfall evaporation
!     *rlpp00*    real     pressure above which radiation is not applied
!     *lphylin*   logical  true when linearized physics is activated 
!     *lenopert   logical  true when no perturbation is required
!                          for surface arrays
!     *leppcfls   logical  true when post-processing of surface fields 
!                          required
!     *lraisanen  logical  true when raisanen overlap scheme is
!                          activated
!     ------------------------------------------------------------------
end module yoephli

module yoephy
  use parkind1, only : jprb
implicit none
  logical, parameter :: leocwa = .false., leocco = .false., lvdftrac = .false.
end module yoephy
module yos_exc
 
use parkind1  ,only : jpim     ,jprb

implicit none

save

real(kind=jprb), parameter :: repust =0.0001_jprb   ! minimum friction velocity (security parameter)

end module yos_exc

