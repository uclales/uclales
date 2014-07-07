subroutine vdfincr(kidia  , kfdia  , klon   , klev  , ktop  ,ptmst  , &
 & pum1   , pvm1   , pslgm1 , ptm1   , pqtm1  , paphm1 , pgeom1 , &
 & pcfm   , ptodc  , psoteu , psotev , psoc   ,&
 & pudif  , pvdif  ,pucurr ,pvcurr ,pslgdif, pqtdif , &
 & pvom   , pvol   , pslge  , pqte   , pslgewodis, &
 & pvdis  , pvdisg, pstrtu , pstrtv , pstrsou , pstrsov , ptofdu , ptofdv)  
!     ------------------------------------------------------------------

!**   *vdfincr* - increments u,v,t and q-tendencies; compute multilevel
!                 fluxes and dissipation.

!     a.c.m. beljaars  18/01/90   derived from vdiff (cy34)
!     a.c.m. beljaars  26/03/90   obukhov-l update
!     m. ko"hler        3/12/2004 conserved variables (qt and slg)
!     a. beljaars       4/04/2005 turbulent orogr. drag acmb
!     a  beljaars      30/09/2005 include subgr. oro. in solver  
!     ocean current b.c.    acmb          12/11/02.
!     p. lopez         02/06/2005 removed option for linearized
!                                 physics (now called separately)
!                                  
!     purpose
!     -------

!     increment u,v,t and q; compute multilevel fluxes and dissipation

!     interface
!     ---------

!     *vdfincr* is called by *vdfmain*

!     input parameters (integer):

!     *kidia*        start point
!     *kfdia*        end point
!     *klev*         number of levels
!     *klon*         number of grid points per packet
!     *ktop*         first level index without zero-diffusion


!     input parameters (real):

!     *ptmst*        double time step (single at 1th step)
!     *pum1*         x-velocity component at t-1
!     *pvm1*         y-velocity component at t-1
!     *pslgm1*       generalized liquid water static energy (slg) at t-1
!     *ptm1*         temperature at t-1
!     *pqtm1*        total water at t-1
!     *paphm1*       pressure at t-1
!     *pgeom1*       geopotential at t-1
!     *pcfm*         prop. to exch. coeff. for momentum (c-star in doc.)
!     *ptodc*        turbulent orographic drag coefficient
!     *psoteu*       explicit part of u-tendency from subgrid orography scheme    
!     *psotev*       explicit part of v-tendency from subgrid orography scheme     
!     *psoc*         implicit part of subgrid orography (df/dt=psote-psoc*f/alpha)
!     *pudif*        u-double tilde devided by alfa
!     *pvdif*        v-double tilde devided by alfa
!     *pucurr*       u-ocean current
!     *pvcurr*       v-ocean current

!     updated parameters (real):

!     *pslgdif*      slg-double tilde devided by alfa (on entry)
!                    slg-single tilde                 (on exit)
!     *pqtdif*       qt-double tilde devided by alfa  (on entry)
!                    qt-single tilde                  (on exit)
!     *pvom*         u-tendency
!     *pvol*         v-tendency
!     *pslge*        slg-tendency
!     *pqte*         qt-tendency

!     output parameters (real):

!     *pvdis*        turbulent dissipation
!     *pvdisg*       subgrid orography dissipation
!     *pstrtu*       turbulent flux of u-momemtum         kg*(m/s)/(m2*s)
!     *pstrtv*       turbulent flux of v-momemtum         kg*(m/s)/(m2*s)
!     *pstrsou*      subgrid orography flux of u-momemtum kg*(m/s)/(m2*s)
!     *pstrsov*      subgrid orography flux of v-momemtum kg*(m/s)/(m2*s)
!     *pslgewodis*   slg-tendency minus (total) dissipation
!     *ptofdu*       tofd comp. of turbulent flux of u-momemtum    kg*(m/s)/(m2*s)
!     *ptofdv*       tofd comp. of turbulent flux of v-momemtum    kg*(m/s)/(m2*s)

!     method
!     ------

!     see documentation

!     ------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
! use yomhook   ,only : lhook    ,dr_hook

use yos_cst   , only : rg
use yoevdf   , only : rvdifts

implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim),intent(in)    :: ktop 
real(kind=jprb)   ,intent(in)    :: ptmst 
real(kind=jprb)   ,intent(in)    :: pum1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pvm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pslgm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: ptm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqtm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: paphm1(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: pgeom1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pcfm(klon,klev) 
real(kind=jprb)   ,intent(in)    :: ptodc(klon,klev) 
real(kind=jprb)   ,intent(in)    :: psoteu(klon,klev) 
real(kind=jprb)   ,intent(in)    :: psotev(klon,klev) 
real(kind=jprb)   ,intent(in)    :: psoc(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pudif(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pvdif(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pucurr(klon) 
real(kind=jprb)   ,intent(in)    :: pvcurr(klon) 
real(kind=jprb)   ,intent(inout) :: pslgdif(klon,klev) 
real(kind=jprb)   ,intent(inout) :: pqtdif(klon,klev) 
real(kind=jprb)   ,intent(inout) :: pvom(klon,klev) 
real(kind=jprb)   ,intent(inout) :: pvol(klon,klev) 
real(kind=jprb)   ,intent(out)   :: pslge(klon,klev) 
real(kind=jprb)   ,intent(out)   :: pqte(klon,klev) 
real(kind=jprb)   ,intent(out)   :: pslgewodis(klon,klev) 
real(kind=jprb)   ,intent(out)   :: pvdis(klon)
real(kind=jprb)   ,intent(out)   :: pvdisg(klon) 
real(kind=jprb)   ,intent(out)   :: pstrtu(klon,0:klev) 
real(kind=jprb)   ,intent(out)   :: pstrtv(klon,0:klev) 
real(kind=jprb)   ,intent(out)   :: pstrsou(klon,0:klev) 
real(kind=jprb)   ,intent(out)   :: pstrsov(klon,0:klev) 
real(kind=jprb)   ,intent(out)   :: ptofdu(klon) 
real(kind=jprb)   ,intent(out)   :: ptofdv(klon) 

!*         0.2    local variables

real(kind=jprb) ::    ztofdu(klon),ztofdv(klon),&
                     &zsou(klon),zsov(klon)

integer(kind=jpim) :: jk, jl

real(kind=jprb) ::    zcons1, zcons2, &
                    & zdudt, zdvdt, zvdfdis, zsodis, ztpfac2, ztpfac3, ztpfac4,&
                    & zdp, zrg, zgdph, ztetofdu,ztetofdv,&
                    & ztesou,ztesov,zhu2,zhu3,zu1,zu2,zu3,zv1,zv2,zv3
! real(kind=jprb) ::    zhook_handle

logical         ::    llocons


!     ------------------------------------------------------------------

!*         1.     initialize constants
!                 --------------------

! if (lhook) call dr_hook('vdfincr',0,zhook_handle)
ztpfac2 = 1.0_jprb/rvdifts
ztpfac3 = 1.0_jprb-ztpfac2
ztpfac4 = 1.0_jprb+ztpfac3

zcons1  = 1.0_jprb/ptmst
zcons2  = 1.0_jprb/(rg*ptmst)

zrg     = 1.0_jprb/rg
llocons  = .true.

!     ------------------------------------------------------------------

!*         2.    compute tendencies and budgets
!                ------------------------------

do jl=kidia,kfdia
  pvdis(jl)    =0.0_jprb
  pvdisg(jl)   =0.0_jprb
  pstrtu(jl,0) =0.0_jprb
  pstrtv(jl,0) =0.0_jprb
  pstrsou(jl,0)=0.0_jprb
  pstrsov(jl,0)=0.0_jprb
  ztofdu(jl)   =0.0_jprb
  ztofdv(jl)   =0.0_jprb
  zsou(jl)   =0.0_jprb
  zsov(jl)   =0.0_jprb
enddo

!*         2.1  vertical loop

do jk=1,klev
  do jl=kidia,kfdia

!   compute total tendencies (dynamics + vertical diffusion + so) 
    zdudt          = ( pudif(jl,jk) - ztpfac2 * pum1(jl,jk) ) * zcons1
    zdvdt          = ( pvdif(jl,jk) - ztpfac2 * pvm1(jl,jk) ) * zcons1
   
    zdp            = paphm1(jl,jk)-paphm1(jl,jk-1)
    zgdph          = -zdp*zrg

!   tofd-tendencies 
    zhu2=-ptodc(jl,jk)*zcons1
    ztetofdu=pudif(jl,jk)*zhu2
    ztetofdv=pvdif(jl,jk)*zhu2
    ztofdu(jl)=ztofdu(jl)+zgdph*ztetofdu
    ztofdv(jl)=ztofdv(jl)+zgdph*ztetofdv
    ptofdu(jl)=ztofdu(jl)
    ptofdv(jl)=ztofdv(jl)

!   implicit part of so-tendencies 
    zhu3=-psoc(jl,jk)*zcons1
    ztesou=pudif(jl,jk)*zhu3
    ztesov=pvdif(jl,jk)*zhu3
    zsou(jl)=zsou(jl)+zgdph*ztesou
    zsov(jl)=zsov(jl)+zgdph*ztesov
   
!   velocity before vdf 
    zu1=pum1(jl,jk)+pvom(jl,jk)*ptmst
    zv1=pvm1(jl,jk)+pvol(jl,jk)*ptmst
!   velocity after vdf 
    zu2=pum1(jl,jk)+(zdudt-psoteu(jl,jk)-ztesou)*ptmst
    zv2=pvm1(jl,jk)+(zdvdt-psotev(jl,jk)-ztesov)*ptmst

!   velocity after so 
    zu3=pum1(jl,jk)+zdudt*ptmst   
    zv3=pvm1(jl,jk)+zdvdt*ptmst   
    zvdfdis=0.5_jprb*(zu1-zu2)*(zu1+zu2) + 0.5_jprb*(zv1-zv2)*(zv1+zv2)
    zsodis =0.5_jprb*(zu2-zu3)*(zu2+zu3) + 0.5_jprb*(zv1-zv2)*(zv1+zv2)

!   integrate vdf-tendencies (including tofd) to find vdf-stress profile
    pstrtu(jl,jk)  = (zdudt-pvom(jl,jk)-psoteu(jl,jk)-ztesou)*zgdph+pstrtu(jl,jk-1)
    pstrtv(jl,jk)  = (zdvdt-pvol(jl,jk)-psotev(jl,jk)-ztesov)*zgdph+pstrtv(jl,jk-1)

!   integrate so-tendencies to find so-stress profile
    pstrsou(jl,jk)  = (psoteu(jl,jk)+ztesou)*zgdph+pstrsou(jl,jk-1)
    pstrsov(jl,jk)  = (psotev(jl,jk)+ztesov)*zgdph+pstrsov(jl,jk-1)

    
    pvom(jl,jk)    = zdudt
    pvol(jl,jk)    = zdvdt
    pvdis(jl)      = pvdis(jl) +zvdfdis*zdp
    pvdisg(jl)     = pvdisg(jl)+zsodis*zdp
   
    pqtdif(jl,jk)  =   pqtdif(jl,jk) + ztpfac3 * pqtm1(jl,jk)
    pqte(jl,jk)    = ( pqtdif(jl,jk) - pqtm1(jl,jk) ) * zcons1

!----------------------------------------------------
!  mpbl - liquid water static energy
!    input:  pslgdif liquid static energy first guess
!            pqtdif  total water first guess
!    output: pslge   liquid static energy tendency
!----------------------------------------------------

      pslgdif(jl,jk)   =   pslgdif(jl,jk) + ztpfac3 * pslgm1(jl,jk)
      pslge(jl,jk)     = ( pslgdif(jl,jk) + zvdfdis+zsodis-pslgm1(jl,jk) ) * zcons1
      pslgewodis(jl,jk)= ( pslgdif(jl,jk)          - pslgm1(jl,jk) ) * zcons1

  enddo
enddo


do jl=kidia,kfdia
  pvdis(jl) =pvdis(jl) *zcons2
  pvdisg(jl)=pvdisg(jl)*zcons2
enddo


! if (lhook) call dr_hook('vdfincr',1,zhook_handle)
end subroutine vdfincr
