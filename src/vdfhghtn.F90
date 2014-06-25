!options xopt(hsfun)
subroutine vdfhghtn (pcog    , pnog    , kidia   , kfdia   , klon    , klev   , kdraft  , ptmst  , kstep    , &
                   & pum1    , pvm1    , ptm1    , pqm1   , plm1    , pim1   , pam1     ,&
                   & paphm1  , papm1   , pgeom1  , pgeoh  , pvervel , pqe    , pte      , &
                   & pkmfl   , pkhfl   , pkqfl   , pmflx  , &
                   & pextr2  , kfldx2  , pextra  , klevx  , kfldx , &
                   & puuh    , pvuh    , pslguh  , pqtuh  , pfracb  , pwuh  , &
                   & pzptop  , kptop   , pzplcl  , kplcl  , kplzb   , &
                   & pwuavg  , pricui  , pmcu    , pdthv  , &
                   & pfplvl  , pfplvn  , pdetr   , &
                   & pbir    , ldnodecp, ldrundry, kpbltype, pwqt2)
                   
                     
!     ------------------------------------------------------------------

!**   *vdfhghtn* - determines the pbl-height and strong updraft fields
!                  using a entraining parcel ascent method.

!     a.p. siebesma    30/06/99    original (dry)
!     m. ko"hler        3/12/2004  moist version
!     roel neggers     12/04/2005  multiple updraft extension


!     purpose
!     -------

!     determine pbl height and updraft fields

!     interface
!     ---------

!     *vdfhghtn* is called by *vdfmain*

!     parameter     description                                   units
!     ---------     -----------                                   -----
!     input parameters (integer):

!     *kidia*        start point
!     *kfdia*        end point
!     *klev*         number of levels
!     *klon*         number of grid points per packet
!     *kdraft*       number of explicitly modeled drafts - currently 3:
!                    1: test parcel
!                    2: rising dry thermals which stop at cloud base or inversion
!                    3: rising dry thermals which become cloudy
!                    (4: downdrafts .. to be done?)

!     input parameters (real):

!     *ptmst*        double time step (single at 1th step)        s
!     *pum1*         x-velocity component at t-1                  m/s
!     *pvm1*         y-velocity component at t-1                  m/s
!     *ptm1*         temperature at t-1                           k
!     *pqm1*         specific humudity at t-1                     kg/kg
!     *plm1*         specific cloud liquid water at t-1           kg/kg
!     *pim1*         specific cloud ice at t-1                    kg/kg
!     *pam1*         cloud fraction at t-1                        kg/kg
!     *paphm1*       pressure at half level at t-1                pa
!     *papm1*        pressure at full level at t-1                pa
!     *pgeom1*       geopotential at t-1                          m2/s2
!     *pgeoh*        geopotential at half level                   m2/s2
!     *pkmfl*        surface kinematic momentum flux              m2/s2  
!     *pkhfl*        surface kinematic heat flux                  k*m/s
!     *pkqfl*        surface kinematic moisture flux              m/s
!     *pbir*         buoyancy-flux integral ratio (-n/p)
!                    used for decoupling criteria

!     *pvervel*      vertical velocity

!     input parameters (logical):

!     *ldnodecp*     true:  never decouple
!                    false: maybe decouple
!     *ldrundry*     true:  run parcel without condensation
!                    false: run parcel with condensation

!     output parameters (real):

!     *pfplvl*       pbl precipitation flux as rain                kg/(m**2*s)
!     *pfplvn*       pbl precipitation flux as snow                kg/(m**2*s)

!     *puuh*         updraft x-momentum
!     *pvuh*         updraft y-momentum
!     *pslguh*       updraft generalized liquid static energy (slg)
!                    at half level                                m2/s2
!     *pqtuh*        updraft specific total water at half level   kg/kg
!     *pmflx*        pbl mass flux                                m/s
!     *pzplcl*       height of lifting condensation level of updraft          m
!     *pzptop*       height of level of zero kinetic energy (w=0) of updraft  m

!     output parameters (integer):

!     *kplcl*         first half level above real height of upraft lcl
!     *kptop*         highest half level below pztop, and
!                       updraft top full level (pztop is within that layer)
!     *kplzb*         level of upraft zero buoyancy (last full level that is pos. buoyant)
!     *kpbltype*    -1: not defined yet
!                    0: stable pbl
!                    1: dry convective pbl (no cloud below parcel top)
!                    2: stratocumulus
!                    3: shallow cumulus
!                    4: deep cumulus

!     method
!     ------

!     see documentation

!     ------------------------------------------------------------------

! use garbage, only : foeewm
use thrm, only : rslf
use parkind1  ,only : jpim     ,jprb

! use yomhook   ,only : lhook,   dr_hook

use yos_cst   , only : rg       ,rd       ,rcpd     ,retv     ,rlvtt    ,&
                    & rlstt    ,ratm     ,rtt      ,rlmlt

use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les    ,&
                    & r4ies    ,r5les    ,r5ies    ,rvtmp2   ,r5alvcp  ,&
                    & r5alscp  ,ralvdcp  ,ralsdcp  ,rtwat    ,rtice    ,&
                    & rticecu  ,rtwat_rtice_r      ,rtwat_rticecu_r  
	      
use yoevdf   , only : rkap     ,rvdifts  ,lldiag
use yoecumf  , only : rtaumel
use yos_exc, only : repust
!use yomgf1c  , only : nc
!use yomlog1c , only : lccn

implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon, pcog, pnog
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: kdraft
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim),intent(in)    :: kstep 
integer(kind=jpim),intent(inout)   :: kplcl(klon,kdraft) 
integer(kind=jpim),intent(inout)   :: kptop(klon,kdraft) 
integer(kind=jpim),intent(inout)   :: kplzb(klon,kdraft) 
real(kind=jprb)   ,intent(in)    :: ptmst 
real(kind=jprb)   ,intent(in)    :: pum1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pvm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: ptm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: plm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pim1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pam1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: paphm1(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: papm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pgeom1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pgeoh(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: pkmfl(klon) 
real(kind=jprb)   ,intent(in)    :: pkhfl(klon) 
real(kind=jprb)   ,intent(in)    :: pkqfl(klon) 
real(kind=jprb)   ,intent(inout)   :: pmflx(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(out)     :: puuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(out)     :: pvuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pslguh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pqtuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pwuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(out)     :: pfracb(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pzplcl(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pzptop(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pwuavg(klon) 
real(kind=jprb)   ,intent(inout) :: pfplvl(klon,0:klev)
real(kind=jprb)   ,intent(inout) :: pfplvn(klon,0:klev)
real(kind=jprb)   ,intent(out)   :: pdetr(klon,klev)
real(kind=jprb)   ,intent(in)    :: pbir(klon) 
real(kind=jprb)   ,intent(out)   :: pricui(klon)
real(kind=jprb)   ,intent(out)   :: pdthv(klon)
real(kind=jprb)   ,intent(out)   :: pmcu(klon) 
logical           ,intent(in)    :: ldnodecp(klon) 
!ldrundry not used now
logical           ,intent(in)    :: ldrundry(klon) 
integer(kind=jpim),intent(inout)   :: kpbltype(klon) 
real(kind=jprb)   ,intent(out)   :: pwqt2(klon,0:klev)
real(kind=jprb)   ,intent(in)    :: pvervel(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pte(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqe(klon,klev) 

!diagnostic output
integer(kind=jpim),intent(in) :: kfldx2, klevx, kfldx
real(kind=jprb)   ,intent(inout) :: pextr2(pnog,kfldx2), pextra(pnog,klevx,kfldx)



!*         0.2    local variables

!--- mean & environmental properties ---
real(kind=jprb) ::    zustar  (klon)     , zwstar(klon)       , zkhvfl(klon)       , &
                    & zusigma(klon)      , zwsigma(klon)      , &
                    & zslgenh(klon,0:klev),zqlenh(klon,0:klev), zqienh(klon,0:klev), &
                    & zqtenh(klon,0:klev), zuenh(klon,0:klev) , zvenh(klon,0:klev) , &
                    & ztven(klon,klev)   , zqtm1 (klon,klev)  , &
                    & zslgm1(klon,klev)  , zmgeom(klon,0:klev), &
                    & ztenh(klon,0:klev) , zrhoh (klon,0:klev), zthven(klon,klev)

!--- updraft parameters ---
real(kind=jprb) ::    zwu2h (klon,0:klev,kdraft), zwuh, &
                    & zqcuh (klon,0:klev,kdraft), zquh  (klon,0:klev,kdraft), &
        & ztuh  (klon,0:klev,kdraft), zeps  (klon,0:klev,kdraft), &
        & zfrac (klon,0:klev,kdraft), &
        & zbuof (klon,klev,kdraft)  , zmcld (klon)              , &
                    & zabulk(klon,0:klev) , zwbulk(klon,0:klev)  , &
                    & zqtbulk(klon,0:klev), zslgbulk(klon,0:klev), &
                    & zubulk(klon,0:klev) , zvbulk(klon,0:klev)  , &
                    & zzplzb(klon,kdraft) , zcape1(klon)
                    
real(kind=jprb) ::    zqsatm, zsatdef, &
                    & zupflxl(klon,0:klev,kdraft), zupflxn(klon,0:klev,kdraft), &
                    & zupgenl(klon,klev,kdraft), zupgenn(klon,klev,kdraft), &
                    & zalfaw, zdzrho, zpflxtot, zpevapup, zfac, zupmelt, &
                    & zaprecevap, zqtep

real(kind=jprb) ::    zfracb(klon,kdraft), zmflxb(klon,kdraft)
      
real(kind=jprb) ::    zfracmax , zfacmaxexc , zfractest , zfactestexc , &
                    & zfacexc(klon,kdraft), zdumfrac, zdumr, zmasscapdepth, &
                    & zpdffacphi(klon), zpdffacw(klon), zlobukhov, &
                    & zstability(klon), zcoupling(klon)

logical ::            lldone(klon,kdraft), llmasscap, llmcind(klon), & 
                    & llwipe, llstcu


integer(kind=jpim) :: is, jk, jl, jd, ibase, jkh, itop

real(kind=jprb) ::    zqexc   , ztexc   , zdz     ,  zdb   , &
                    & zspeedenv         , zspeedup, zwindir , &
                    & zz      , ztotw2(klon)      , ztotp(klon) , &
                    & zcons10 , zcfnc1(klon,0:klev)      , ztvmean     , &
                    & zrg     , zmfmax  , zmfs(klon,kdraft)

!          remaining model parameters

real(kind=jprb) ::    ztaueps , zclddepth     , &
                    & zw2thresh               , zstabthresh    , zbirthresh , &
                    & ztvlim  , zclddepthdp   , zdzcloud(klon) , zw2h       , &
                    & zzfunc  , zzfunc3(klon) , zzfunc4(klon)  , zdhri(klon)    , zdhcl      , &
                    & zdthvdz        , zghm1
		    
real(kind=jprb) ::    zzi(klon)

integer(kind=jpim) :: izi(klon,kdraft)

real(kind=jprb) ::    zbuoycu, zdthvcutop, zricuinorm, zdmdz

real(kind=jprb) ::    zfracbplus(klon,0:klev), zwuhbplus(klon,0:klev), &
                    & zqtuhbplus(klon,0:klev), zslguhbplus(klon,0:klev)

logical ::            llcape, llcapetest

real(kind=jprb) ::    zrhs(klon), zmasscape(klon), ztaubm, &
                    & zfracbcong(klon), zmscale(klon)

! real(kind=jprb) :: zhook_handle
! 
! 
! interface
! #include "surf_inq.h"
! end interface
! 
! 
! #include "vdfparcel.intfb.h"
! #include "vdfstcucrit.intfb.h"
! #include "vdfpdftable.intfb.h"
! #include "vdfbuoysort.intfb.h"
! 
! #include "fcttre.h"
! 


!     ------------------------------------------------------------------

!*         1.     initialization
!                 --------------

! if (lhook) call dr_hook('vdfhghtn',0,zhook_handle)

!if (lccn) then
!  write(0,*) 'vdfhghtn: nc=',nc
!endif


!-- top % of the pdf associated with the test parcel --
!zfractest   = 0.0002_jprb    
!zfractest   = 0.001_jprb    
zfractest   = 0.002_jprb      !cy32r3 
!zfractest   = 0.005_jprb    
call vdfpdftable (zfractest, zfactestexc, zdumr, zdumr, 0) ! associated pdf scaling factor

!-- total convective area fraction that is done with mass flux --
!zfracmax    = 0.05_jprb     
zfracmax    = 0.075_jprb     !cy32r3
!zfracmax    = 0.1_jprb     
!zfracmax    = 0.15_jprb     
call vdfpdftable (zfracmax, zfacmaxexc, zdumr, zdumr, 0) ! associated pdf scaling factor

!-- eddy turnover time scale used in parcel entrainment [s]  (neggers, siebesma & jonker, jas 2002) --
!ztaueps     = 300._jprb    
ztaueps     = 400._jprb    ! cy32 
!ztaueps     = 500._jprb    

!-- threshold parcel vertical velocity squared [m2/s2] --
!zw2thresh  = -1._jprb     
zw2thresh   = 0._jprb      

!-- threshold cloud thickness for stcu/cu transition [m] --
zclddepth   = 2000._jprb   

!-- threshold cloud thickness used in shallow/deep decision [m] --
!zclddepthdp = 2000._jprb   
!zclddepthdp = 3000._jprb   !cy32
zclddepthdp = 100000._jprb   

zstabthresh = 20._jprb     ! threshold stability (klein & hartmann criteria) [k]
zbirthresh  = 0.1_jprb     ! threshold bir (tke decoupling criteria) [1]
ztvlim      = 0.1_jprb     ! cloud fraction limit in tv,env calculation

!-- switch for moist mass flux depth limiter --      
!llmasscap     = .true.
llmasscap     = .false.    
zmasscapdepth = 3000._jprb 
!zmasscapdepth = 50000._jprb 

!-- switch for applying klein-hartmann criterion for stratocumulus --
!llstcu = .true.
llstcu = .false.

!-- factor used in updraft initialization --
do jl=kidia,kfdia
  !zpdffacw(jl)   = 0.2_jprb
  !zpdffacphi(jl) = 0.2_jprb
  zpdffacw(jl)   = 1.0_jprb
  zpdffacphi(jl) = 1.0_jprb
enddo  

!-- updraft precip evaporation constant --
zaprecevap = 0.001_jprb
!zaprecevap = 0.000544_jprb
!zaprecevap = 0.0001_jprb

!-- settings for fritsch-chappell cape-removal --
!llcape     = .true.                  !activation switch
llcape     = .false.
ztaubm     = 3600._jprb * 1._jprb    !the associated cape adjustment timescale
!ztaubm     = 3600._jprb * 2._jprb
llcapetest = .true.                  !option i: test updraft carries the required additional transport
!llcapetest = .false.                 !option ii: moist updraft carries the required additional transport
  
!-- optimization --
zrg    = 1.0_jprb/rg


! set some stuff to zero
do jl=kidia,kfdia
  
  pwuavg(jl)     = 0.0_jprb
  kpbltype(jl)   = -1          ! -1 means: yet unknown
  
  zzi(jl)        = 0._jprb      !mixed layer scalings
  zwstar(jl)     = 0._jprb        
  
  pricui(jl)  = 1._jprb       ! 1 / cumulus inversion richardson number
  pdthv(jl)   = 0._jprb
   
  zcape1(jl)     = 0._jprb

  zmcld(jl)       = 0._jprb
  pmcu(jl)        = 0._jprb       !cloud-depth average moist updraft mass flux
  
  zcoupling(jl) = 0._jprb
  
  zrhs(jl)      = 0._jprb 
  
enddo


do jd=1,kdraft
  do jl=kidia,kfdia
    pzplcl(jl,jd)  = -100._jprb  ! default value: -100 (no lcl)
    pzptop(jl,jd)  = 0._jprb     
    kplcl(jl,jd)   = 0           ! default value: 0 (no pbl cloud)
    kptop(jl,jd)   = 0          
    kplzb(jl,jd)   = 0          
    lldone(jl,jd)  = .true.       ! default: true (don't launch the parcel)
    zfracb(jl,jd)  = 0._jprb 
    pfracb(jl,jd)  = 0._jprb 
    zfacexc(jl,jd) = 0._jprb 
    zfacexc(jl,jd) = 0._jprb 
    zmflxb(jl,jd)  = 0._jprb 
    izi(jl,jd)     = 0._jprb        
  enddo
enddo


do jk=0,klev
  do jl=kidia,kfdia
    pwqt2(jl,jk) = 0._jprb  
  enddo
enddo


do jk=1,klev
  do jl=kidia,kfdia
    pdetr(jl,jk) = 0._jprb
  enddo
enddo    


!--- parcel half level parameters ---
do jd=1,kdraft
  do jk=0,klev
    do jl=kidia,kfdia
    puuh(jl,jk,jd)  = 0.0_jprb
    pvuh(jl,jk,jd)  = 0.0_jprb
    pslguh(jl,jk,jd)= 0.0_jprb
    pqtuh(jl,jk,jd) = 0.0_jprb
    pmflx(jl,jk,jd) = 0.0_jprb
    ztuh(jl,jk,jd)  = 0.0_jprb
    zquh(jl,jk,jd)  = 0.0_jprb
    zqcuh(jl,jk,jd) = 0.0_jprb
    zeps(jl,jk,jd)  = 0.0_jprb
    zwu2h(jl,jk,jd) = 0.0_jprb
    zfrac(jl,jk,jd) = 0.0_jprb
    zupflxl(jl,jk,jd)  = 0.0_jprb
    zupflxn(jl,jk,jd)  = 0.0_jprb
    enddo
  enddo
enddo


!--- parcel full level parameters ---
do jd=1,kdraft
  do jk=1,klev
    do jl=kidia,kfdia
      zbuof(jl,jk,jd)    = 0.0_jprb
      zupgenl(jl,jk,jd)  = 0.0_jprb
      zupgenn(jl,jk,jd)  = 0.0_jprb
    enddo
  enddo
enddo


!--- reset output stuff ---
if (lldiag) then

!  do jl=kidia,kfdia
!    pextr2(jl,1:49) = 0._jprb
!  enddo

  do jk=1,klevx
  do jl=kidia,kfdia
    pextra(jl,jk,49) = 0._jprb
    pextra(jl,jk,50) = 0._jprb
  enddo
  enddo

endif


!     -----------------------------------------------------------------
!
!*         2.     prepare fields on half levels by linear interpolation
!*                of conserved variables
!                 -----------------------------------------------------

!*         2.1  full level cpm, slg, qt and tv
!*

  do jk=1,klev
    do jl=kidia,kfdia
      zslgm1(jl,jk) = rcpd * ptm1(jl,jk) + pgeom1(jl,jk) &
                  & - rlvtt * plm1(jl,jk) - rlstt * pim1(jl,jk)  
      zqtm1 (jl,jk) = pqm1(jl,jk) + plm1(jl,jk) + pim1(jl,jk)

!          parcel goes through cloud portion of environment
!          (added ql loading; ql,cld=ql,mean/fc; qv = qsat) 
!          safety: fc>0.1; linear interpolation between overcast 
!                  and cloudy portion for 0<fc<0.1
!                  guaranteed to be < tv from mean conditions

!          grid box mean virtual effect
      ztvmean       = ptm1(jl,jk) * ( 1.0_jprb + retv * pqm1(jl,jk) &
                  & - plm1(jl,jk) - pim1(jl,jk) )       !qli loading  
      ztven(jl,jk) = ztvmean
      zthven(jl,jk) = ( papm1(jl,jk)/ratm )**(-rd/rcpd) * ztven(jl,jk)
    enddo
  enddo


!*         2.2  half-level environment interpolation (qt, ql, qi, slg)
!*

  do jk=1,klev-1
  do jl=kidia,kfdia

    if (jk==1) then
      zghm1 = pgeoh(jl,jk) + 50000._jprb*rg   !avoid using top half level (=inf)
    else
      zghm1 = pgeoh(jl,jk-1)
    endif  
    
    zqtenh(jl,jk) = ( zqtm1(jl,jk+1) *(zghm1-pgeoh(jl,jk  )) &
                & +   zqtm1(jl,jk)   *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                &   )                /(zghm1-pgeoh(jl,jk+1))
    zqlenh(jl,jk) = ( plm1(jl,jk+1)  *(zghm1-pgeoh(jl,jk  )) &
                & +   plm1(jl,jk)    *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                &   )                /(zghm1-pgeoh(jl,jk+1))
    zqienh(jl,jk) = ( pim1(jl,jk+1)  *(zghm1-pgeoh(jl,jk  )) &
                & +   pim1(jl,jk)    *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                &   )                /(zghm1-pgeoh(jl,jk+1))
    zslgenh(jl,jk)= ( zslgm1(jl,jk+1)*(zghm1-pgeoh(jl,jk  )) &
                & +   zslgm1(jl,jk)  *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                &   )                /(zghm1-pgeoh(jl,jk+1))
    zuenh(jl,jk)  = ( pum1(jl,jk+1)  *(zghm1-pgeoh(jl,jk  )) &
                & +   pum1(jl,jk)    *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                &   )                /(zghm1-pgeoh(jl,jk+1))
    zvenh(jl,jk)  = ( pvm1(jl,jk+1)  *(zghm1-pgeoh(jl,jk  )) &
                & +   pvm1(jl,jk)    *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                &   )                /(zghm1-pgeoh(jl,jk+1))

    !calculate t at half levels from sl, for later use in density calculations    
    !ztenh        = ( zslgenh (jl,jk) - pgeoh(jl,jk) &
    !                   & + rlvtt*zqlenh(jl,jk) + rlstt*zqienh(jl,jk) &
    !                   & ) / rcpd
    
    ztenh(jl,jk)  =  ( ptm1(jl,jk+1)  *(zghm1-pgeoh(jl,jk  )) &
                 & +   ptm1(jl,jk)    *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                 &   )                /(zghm1-pgeoh(jl,jk+1))
                
    !air density at half levels
    zrhoh(jl,jk) = paphm1(jl,jk)/(rd*ztenh(jl,jk))
      
  enddo
  enddo



!     -----------------------------------------------------------------

!*         3.     release the first (test) updraft to get pbl heights
  
  
  !* set updraft index to 1
  jd = 1   

  do jl=kidia,kfdia
    
 
    pfracb(jl,jd) = zfractest
 
 
    !* 3.1    determine stability of bl using the surface buoyancy flux
    !*
    zkhvfl(jl)  = ( 1.0_jprb + retv *  zqtm1(jl,klev) ) * pkhfl(jl) + &
                & ( retv * zslgm1(jl,klev) / rcpd )     * pkqfl(jl) 

    if ( zkhvfl(jl) >= 0.0_jprb ) then
      
      ! stable bl (no updrafts expected/needed)
      kpbltype(jl)  = 0

    else

      lldone(jl,jd) = .false.  !confirm launch
     
     
      !* 3.2    surface layer scaling
      !*
      zustar(jl)  = max( sqrt(pkmfl(jl)), repust )               !u* (repust=10e-4)
      zwstar(jl)  = (- zkhvfl(jl) * rg / ptm1(jl,klev) * 1000._jprb ) &   !zi=1000m
                       & ** ( 1._jprb/3._jprb) 
      zwsigma(jl)      = 1.2_jprb &
       & * ( zustar(jl)**3 &
       & - 1.5_jprb * rkap * zkhvfl(jl) * pgeoh(jl,klev-1) / ptm1(jl,klev-1) &
       & ) ** ( 1.0_jprb/3._jprb )                         ! kolmogorov 1/3-power
      zusigma(jl)      = 2.29_jprb  &
       & * ( zustar(jl)**3 &
       & + 0.5_jprb / 12.0_jprb * rkap * zwstar(jl)**3 &
       &   ) ** ( 1._jprb/3._jprb)
      
      !  scaling factors between phi* and initial updraft phi excesses
      !zlobukhov = -zustar(jl)**3 * ptm1(jl,klev) / (rg * rkap * zkhvfl(jl))
      !zpdffacphi(jl) = zrg*pgeoh(jl,klev-1)/zlobukhov
      
      
      !* 3.3    initialize test updraft
      !*
      
      !get the constant associated with the top zfractest % of the pdf
      zfacexc(jl,1) = zfactestexc
      
      !calculate the initial excess values
      zwu2h(jl,klev-1,jd) = ( zpdffacw(jl) * zfacexc(jl,1) * zwsigma(jl) )**2         
      ztexc            = - zpdffacphi(jl) * zfacexc(jl,1) * pkhfl(jl) / zwsigma(jl) 
      zqexc            = - zpdffacphi(jl) * zfacexc(jl,1) * pkqfl(jl) / zwsigma(jl) 
      ztexc            = max(ztexc, 0.0_jprb)
      zqexc            = max(zqexc, 0.0_jprb)
      pqtuh(jl,klev-1,jd) = zqtenh(jl,klev-1) + zqexc
      zqcuh(jl,klev-1,jd) = zqlenh(jl,klev-1) + zqienh(jl,klev-1)
      zquh (jl,klev-1,jd) = pqtuh(jl,klev-1,jd)  - zqcuh(jl,klev-1,jd)
      pslguh(jl,klev-1,jd)= zslgenh(jl,klev-1) + rcpd * ztexc
      ztuh (jl,klev-1,jd) = ( pslguh (jl,klev-1,jd) - pgeoh(jl,klev-1) &
                       & + rlvtt*zqlenh(jl,klev-1) + rlstt*zqienh(jl,klev-1) &
                       & ) / rcpd

! ... u & v: (wind speed assumed to be negatively correlated with t and q excesses)
      zspeedenv        = sqrt( zuenh(jl,klev-1)**2 + zvenh(jl,klev-1)**2 )
      zspeedup         = max( zspeedenv - zfacexc(jl,1) * zusigma(jl), 0._jprb )
!     zspeedup         = max( zspeedenv - zfacexc(jl,1) * zustar(jl)**2/ zwsigma(jl) , 0._jprb )

!      puuh(jl,klev-1,jd)= zuenh(jl,klev-1) * zspeedup/zspeedenv
!      pvuh(jl,klev-1,jd)= zvenh(jl,klev-1) * zspeedup/zspeedenv
      puuh(jl,klev-1,jd)= zuenh(jl,klev-1)  !mean wind at this half level
      pvuh(jl,klev-1,jd)= zvenh(jl,klev-1)
!      puuh(jl,klev-1,jd)= pum1(jl,klev)     !10m wind instead of 20m
!      pvuh(jl,klev-1,jd)= pvm1(jl,klev)
!      puuh(jl,klev-1,jd)= 0._jprb           !using this is really bad for the scores!
!      pvuh(jl,klev-1,jd)= 0._jprb
!      puuh(jl,klev-1,jd)= pum1(jl,klev-1)   !full level wind
!      pvuh(jl,klev-1,jd)= pvm1(jl,klev-1)
      
    endif

  enddo !jl


  !* 3.4   release the test updraft #1
  !*          - used to make a first guess of the heights of cloud base & inversion,
  !*            and to determine pbl type.
  !*
  call vdfparcel (kidia   , kfdia   , klon    , klev    , kdraft  , &
                & pgeoh   , pgeom1  , paphm1  , &
    & pum1    , pvm1    , zqtm1   , zslgm1  , ztven   , &
    & puuh    , pvuh    , pslguh  , pqtuh   , zwu2h   , zqcuh  , zbuof , & 
    & zquh    , ztuh    , zeps    , zfacexc , &
    & pzplcl  , kplcl   , pzptop  , kptop   , kplzb   , &
    & jd      , zupgenl , zupgenn , &
    & ztaueps , zw2thresh, lldone , kpbltype)  



!     -----------------------------------------------------------------

!*         4.     classification of the convective pbl
!                 ------------------------------------


  !* 4.1    classify the convective pbl
  !*
  do jl=kidia,kfdia
    if ( kpbltype(jl)/=0 ) then
 
      if ( pzplcl(jl,1) > pzptop(jl,1) .or. kplcl(jl,1) == 0 ) then
      
        !dry convective pbl
        kpbltype(jl)  = 1                   !dry convective pbl
        zdzcloud(jl)  = 0.0_jprb            !cloud thickness

      else

        !moist convective pbl
        zdzcloud(jl)  = pzptop(jl,1) - pzplcl(jl,1) !cloud thickness
            
        if (zdzcloud(jl)>zclddepthdp .and. .not.llmasscap ) then
        
	  !deep convection
    kpbltype(jl) = 4
          
  else
        
          if (llstcu) then
            kpbltype(jl) = 2   !set the type to stratocumulus for the moment
          else
            kpbltype(jl) = 3   !rn run without klein-hartmann criterion!
          endif  
          
        endif
	
      endif

    endif !kpbltype /=0
  enddo !jl
  
  
  
  !* 4.2    check the stratocumulus/shallow cumulus criterion (trigger function)
  !*        if shallow cumulus is diagnosed, kpbltype will be set to 3
  !*
  call vdfstcucrit ( kidia , kfdia  , klon  , klev , kdraft , &
      &    ptm1  , zslgm1 , zqtm1 , papm1 , &
      &    zstabthresh, zclddepth, zbirthresh, zdzcloud, &
      &    kptop , kpbltype, ldnodecp, &
                  &    zstability )
  
  !-- formulate decoupling constraint on moist updraft area fraction --                    
  do jl=kidia,kfdia
    
    zcoupling(jl) = max( 0._jprb, (zstability(jl) - 20._jprb) / 2._jprb )
    zcoupling(jl) = min( zcoupling(jl), 1._jprb )
    
    zzfunc4(jl) = zfracmax * zcoupling(jl)
    
    !zzfunc4(jl) = zfracmax
    !zzfunc4(jl) = 0._jprb
  
!    if (lldiag) then
!      pextr2(jl,10) = zstability(jl)
!    endif 

  enddo

  
  
!     -----------------------------------------------------------------

!*         5.     closure for organized updrafts (jd=2,3)
!                 ---------------------------------------


  !* 5.1    mixed layer scalings
  !*

  do jl=kidia,kfdia
    
    if ( kpbltype(jl)/=0 ) then    !don't do this for stable pbl
      
      !--- mixed layer scaling depth ---
      select case (kpbltype(jl))
    
        case(1)
          !dry convective pbl - inversion height
          zzi(jl)   = pzptop(jl,1)
          izi(jl,1) = kptop(jl,1)
		
        case(2)
          !stratocumulus - inversion height
          !caution: during decoupling in the intermediate regime (e.g. astex/atex) the
    !   relevant ml scaling height changes from pbl inversion to level of minimum
          !   buoyancy flux. in the current setup this is not modelled yet!
          zzi(jl)   = pzptop(jl,1)
          izi(jl,1) = kptop(jl,1)
	  
        case(3)
          !shallow cumulus - level of minimum buoyancy flux
    !assume that the moist updraft lcl is very close to this level
          zzi(jl)   = pzplcl(jl,1)
          izi(jl,1) = kplcl(jl,1)
	  
        case(4)
          !deep cumulus - only do a dry parcel up to cloud base
          zzi(jl)   = pzplcl(jl,1)
          izi(jl,1) = kplcl(jl,1)
		
      end select 

      !--- mixed layer convective velocity scale ---
      zwstar(jl) = ( -zkhvfl(jl) * rg * zzi(jl) / zthven(jl,klev)  ) ** (1._jprb/3._jprb)
      
    endif
    
  enddo  

    
    
  !*  5.2    ri number of cumulus inversion
  !* 
  
  !-- test-updraft cloudy cape --
  do jk=klev-1,1,-1
    do jl=kidia,kfdia
      if ( zqcuh(jl,jk,1)>0._jprb .and. zbuof(jl,jk,1)>0._jprb .and. jk<=kplcl(jl,1) ) then  
        zdz = zrg*( pgeoh(jl,jk-1) - pgeoh(jl,jk) )
        zcape1(jl) = zcape1(jl) + zdz * zbuof(jl,jk,1) 
      endif
    enddo
  enddo
  
  do jl=kidia,kfdia
    if ( kpbltype(jl)==2 .or. kpbltype(jl)==3 ) then
    
      !-- interpolate lzb height --
      zzplzb(jl,1) = pzptop(jl,1)
      if (kplzb(jl,1)>2) then
        zdz = (pgeoh(jl,kplzb(jl,1)-1) - pgeoh(jl,kplzb(jl,1)))*zrg
        zdb = zbuof(jl,kplzb(jl,1),1)  - zbuof(jl,kplzb(jl,1)-1,1)
        if (zdb>0._jprb) then
          zzplzb(jl,1) = pgeoh(jl,kplzb(jl,1)) * zrg + &
                       & zdz * zbuof(jl,kplzb(jl,1),1) / zdb 
          zzplzb(jl,1) = min( zzplzb(jl,1), pgeoh(jl,kplzb(jl,1)-1)*zrg )
        endif               
      endif
      
      !-- cloud layer average test-updraft positive buoyancy --
      zbuoycu = 0._jprb
      if ( zzplzb(jl,1)-pzplcl(jl,1)>0._jprb ) then  
  zbuoycu = zcape1(jl) / ( zzplzb(jl,1) - pzplcl(jl,1) )
      endif  
  
      !-- inversion theta_v jump --
      !jk = kplzb(jl,1)    !use level of zero buoyancy (lzb) of test-updraft
      jk = kptop(jl,1)    !use top level of test-updraft
      zdthvcutop = 0
      if (jk>2) then
        zdthvcutop = max( zthven(jl,jk-1)-zthven(jl,jk), zthven(jl,jk)-zthven(jl,jk+1) )
      
        !-- cumulus ri number - used again in vdfexcu --
        if ( zdthvcutop > 0._jprb ) then 
          pricui(jl) = zbuoycu * zrg * zthven(jl,klev) / zdthvcutop   
        endif  
      endif  
      
      !-- rn testing: no top entrainment for stcu (yikes) --
      !if ( zstability(jl) > zstabthresh ) then
      !  pricui(jl) = 0._jprb
      !endif
      !pricui(jl) = ( 1._jprb - zcoupling(jl) ) * pricui(jl)
      
!      if (lldiag) then
!        pextr2(jl,14) = zbuoycu
!        pextr2(jl,15) = zdthvcutop
!        pextr2(jl,16) = rg * zdthvcutop / zthven(jl,klev)
!        pextr2(jl,17) = min(pricui(jl),1._jprb)
!      endif  

    endif
  enddo
  
  

  !* 5.3    closure of updraft area fractions (jd=2,3)
  !*
    
  do jl=kidia,kfdia
    
    if ( kpbltype(jl)/=0 ) then    !don't do this for stable pbl


      if (kpbltype(jl)>1) then
      
        !--- transition layer depth, scale i: dry entrainment layer depth ---
        !       using thv gradient averaged over a number of layers above h
        !
        zdthvdz  = max( 0.01_jprb,zthven(jl,izi(jl,1)-2)-zthven(jl,izi(jl,1)) ) * rg / &
               & ( pgeom1(jl,izi(jl,1)-2) - pgeom1(jl,izi(jl,1)) )

        zw2h      = 0.5_jprb * zwstar(jl)**2.
        zdhri(jl) = (   zw2h * zrg * zthven(jl,klev) / (0.5_jprb * zdthvdz)  )**0.5
        pdthv(jl) = zdhri(jl) * zdthvdz   !used again in vdfexcu
      
        zzfunc3(jl) = max( 0._jprb, 0.2_jprb * zdhri(jl) / pzplcl(jl,1) )   !cy32r3
      
!        if (lldiag) then
!          pextr2(jl,50) = pzplcl(jl,1)
!          pextr2(jl,51) = 100._jprb * zzfunc3(jl)
!          pextr2(jl,52) = zdhri(jl)
!          pextr2(jl,48) = zdthvdz
!          pextr2(jl,49) = pdthv(jl)
!          pextr2(jl,55) = kplcl(jl,1)-2 - kptop(jl,1)
!        endif      

      endif
      
      
      !--- calculation of moist updraft area fraction ---
      select case (kpbltype(jl))
    
        case(1)
	
          !dry convective pbl
    !set area fraction of moist group to zero
    zfracb(jl,3) = 0._jprb

        case(2)
	
          !stratocumulus
    !set area fraction of moist group to zfracmax
    zfracb(jl,3) = zfracmax

        case(3)
	  
          !shallow cumulus
          !flexible updraft area fractions

    !-- transition layer depth, scale ii: cumulus condensation depth-scale --
          !
          zdhcl = min(200._jprb, 0.1_jprb * zdzcloud(jl))
          zdhcl = max(zdhcl,0._jprb)
          
          !use m/w* ~ dh / h  (neggers et al., qj, 2007)
          if (pzplcl(jl,1).gt.0._jprb) then
            zzfunc =  0.15_jprb * ( zdhcl/pzplcl(jl,1) )     !cy32r3
            !zzfunc =  0.2_jprb * ( zdhcl/pzplcl(jl,1) )   
          else
            zzfunc = 0._jprb
          endif
          
!          if (lldiag) then
!            pextr2(jl,53) = zdhcl
!            pextr2(jl,44) = 100._jprb*zzfunc
!          endif  
          
          !-- choose the minimum of scales i and ii --
          zfracb(jl,3)  = min( zfracmax, zzfunc, zzfunc3(jl) )
          zfracb(jl,3)  = max( zfracb(jl,3), zzfunc4(jl) )  !superimpose decoupling criterion
              
          !-- 1st guess for the lcl mass flux --
          zmflxb(jl,3)  = zwstar(jl) * zfracb(jl,3) * zrhoh(jl,kplcl(jl,1))
          
          !-- switch kpbltype to dry convective if moist updraft is not launched --
          if (zfracb(jl,3).eq.0._jprb) then
            kpbltype(jl)=1
          endif  
	  
	  
        case(4)
	
          !deep cumulus
    !set area fraction of moist group to zero (only allow updraft transport in dry mixed layer)
    zfracb(jl,3) = 0._jprb


      end select !kpbltype
      
      
      !--- dry updraft area fraction (jd=2) ---
      zfracb(jl,2) = max( 0._jprb, zfracmax - zfracb(jl,3) )
      
      pfracb(jl,2) = zfracb(jl,2)
      pfracb(jl,3) = zfracb(jl,3)


    endif !kpbltype /=0
    
  enddo !jl


      
!     -----------------------------------------------------------------

!*         6.     calculate vertical profiles of all updrafts (jd=2,3)
!                 ----------------------------------------------------


  !*       6.1    calculate the scaling factors of the updraft excess with the surface joint pdfs
  !*
  do jd = 2,kdraft
    do jl=kidia,kfdia
      
      if ( kpbltype(jl)/=0 .and. zfracb(jl,jd)>0._jprb ) then
        
        !-- get the pdf scaling factor --
  select case (jd)
  	  
    case(2)
      !lower part of top zfracmax %
      zdumfrac = zfracmax - zfracb(jl,2)
            call vdfpdftable(zdumfrac , zfacexc(jl,2), zdumr, zdumr, 0)
      zfacexc(jl,2) = ( zfracmax * zfacmaxexc - zdumfrac * zfacexc(jl,2) ) / zfracb(jl,2)
		    
    case(3)
      !upper part of top zfracmax %
      zdumfrac = zfracb(jl,jd)
            call vdfpdftable(zdumfrac , zfacexc(jl,3), zdumr, zdumr, 0)
	    
  end select
	
      endif !kpbltype & zfracb

    enddo !jl
  enddo !jd
    
          
    
  !*       6.2    vertical integration of dry & moist updraft budgets (jd=2,3)
  !*
  do jd = 2,kdraft
  
    !-- initialize updraft --
    do jl=kidia,kfdia
      
      if ( kpbltype(jl)/=0 .and. zfracb(jl,jd)>0._jprb ) then
        
        lldone(jl,jd) = .false. !confirm launch
          
        zwu2h(jl,klev-1,jd) = ( zpdffacw(jl) * zfacexc(jl,jd) * zwsigma(jl) )**2 
        ztexc            = - zpdffacphi(jl) * zfacexc(jl,jd) * pkhfl(jl) / zwsigma(jl) 
        zqexc            = - zpdffacphi(jl) * zfacexc(jl,jd) * pkqfl(jl) / zwsigma(jl) 
        ztexc            = max(ztexc, 0.0_jprb)
        zqexc            = max(zqexc, 0.0_jprb)
        pqtuh(jl,klev-1,jd) = zqtenh(jl,klev-1) + zqexc
        zqcuh(jl,klev-1,jd) = zqlenh(jl,klev-1) + zqienh(jl,klev-1)
        zquh (jl,klev-1,jd) = pqtuh(jl,klev-1,jd)  - zqcuh(jl,klev-1,jd)
        pslguh(jl,klev-1,jd)= zslgenh(jl,klev-1) + rcpd * ztexc
        ztuh (jl,klev-1,jd) = ( pslguh (jl,klev-1,jd) - pgeoh(jl,klev-1) &
                       & + rlvtt*zqlenh(jl,klev-1) + rlstt*zqienh(jl,klev-1) &
                       & ) / rcpd

!   ... u & v: (wind speed assumed to be negatively correlated with t and q excesses)
        zspeedenv        = sqrt( zuenh(jl,klev-1)**2 + zvenh(jl,klev-1)**2 )
        zspeedup         = max( zspeedenv - zfacexc(jl,jd) * zusigma(jl), 0._jprb )
!       zspeedup         = max( zspeedenv - zfacexc(jl,jd) * zustar(jl)**2/ zwsigma(jl) , 0._jprb )

!        puuh(jl,klev-1,jd)= zuenh(jl,klev-1) * zspeedup/zspeedenv
!        pvuh(jl,klev-1,jd)= zvenh(jl,klev-1) * zspeedup/zspeedenv
        puuh(jl,klev-1,jd)= zuenh(jl,klev-1)  !mean wind at this half level
        pvuh(jl,klev-1,jd)= zvenh(jl,klev-1)
!        puuh(jl,klev-1,jd)= pum1(jl,klev)     !10m wind instead of 20m
!        pvuh(jl,klev-1,jd)= pvm1(jl,klev)
!        puuh(jl,klev-1,jd)= 0._jprb           !using this is really bad for the scores!
!        pvuh(jl,klev-1,jd)= 0._jprb
!        puuh(jl,klev-1,jd)= pum1(jl,klev-1)   !full level wind
!        pvuh(jl,klev-1,jd)= pvm1(jl,klev-1)
 
      endif !kpbltype & zfracb
      
    enddo !jl
    
    
    
    
    !-- release the updraft --
    call vdfparcel (kidia   , kfdia   , klon    , klev    , kdraft  , &
                  & pgeoh   , pgeom1  , paphm1  , &
      & pum1    , pvm1    , zqtm1   , zslgm1  , ztven   , &
      & puuh    , pvuh    , pslguh  , pqtuh   , zwu2h   , zqcuh  , zbuof , & 
      & zquh    , ztuh    , zeps    , zfacexc , &
      & pzplcl  , kplcl   , pzptop  , kptop   , kplzb   , &
      & jd      , zupgenl , zupgenn , &
      & ztaueps , zw2thresh, lldone , kpbltype)  
    
  enddo !jd 



  !*        6.3 some protectional measures against updraft #3 failing (only just) to reach condensation.
  !*
  do jl=kidia,kfdia

    if ( kpbltype(jl)==2 .or. kpbltype(jl)==3) then
      if ( zfracb(jl,3)>0._jprb .and. pzplcl(jl,3)<0._jprb) then
      
        !set cloud base height to updraft #3 top (cloud layer exists but has zero depth)
        kplcl(jl,3)  = kptop(jl,3)
        pzplcl(jl,3) = pzptop(jl,3)
        
      endif
      
    endif
     
  enddo !jl 


    
  !*        6.4  limiter for cloudy depth of moist updraft  *testing*
  !* 
  if (llmasscap) then
  
    do jl=kidia,kfdia
      llmcind(jl) = .false.
      if (kpbltype(jl)==3) then
        if ( pzplcl(jl,3)>0._jprb .and. (pzptop(jl,3) - pzplcl(jl,3))>zmasscapdepth ) then
          pzptop(jl,3) = pzplcl(jl,3) + zmasscapdepth
          llmcind(jl) = .true.
        endif
      endif
    enddo
    
    do jk=klev-1,1,-1
      do jl=kidia,kfdia
      
        if ( llmcind(jl) ) then
          if ( pgeoh(jl,jk+1)*zrg<=pzptop(jl,3) .and. pgeoh(jl,jk)*zrg>pzptop(jl,3) ) then
            kptop(jl,3) = jk+1
          endif
        endif
        
      enddo
    enddo
    
  endif


  !*        6.5  updraft precipitation fluxes (rain and snow)
  !* 
  do jd = 3,kdraft  !moist updrafts only
    
    do jk=2,klev
      do jl=kidia,kfdia
      
        zdzrho = zrg * ( paphm1(jl,jk)-paphm1(jl,jk-1) ) 
        
        !-- add precip generation to flux [kg /m2 /s: tendency * layer depth * air density] --
        zupflxl(jl,jk,jd) = zupflxl(jl,jk-1,jd) + zupgenl(jl,jk,jd) * zdzrho
        zupflxn(jl,jk,jd) = zupflxn(jl,jk-1,jd) + zupgenn(jl,jk,jd) * zdzrho
        
        !-- do some melting at freezing level (snow->rain) --
        if (zupflxn(jl,jk,jd)>0._jprb .and. ptm1(jl,jk) > rtt) then
          zupmelt = (1.0_jprb+0.5_jprb*(ptm1(jl,jk)-rtt)) * &
                  & (ptm1(jl,jk)-rtt) * rcpd/(rlmlt*rtaumel) * zdzrho
          zupmelt = min(zupflxn(jl,jk,jd),zupmelt)
          zupflxl(jl,jk,jd) = zupflxl(jl,jk,jd) + zupmelt
          zupflxn(jl,jk,jd) = zupflxn(jl,jk,jd) - zupmelt
        endif  
        
        !-- saturation deficit of mean state t --
!         zqsatm = foeewm(ptm1(jl,jk))/papm1(jl,jk)
!         zqsatm = min(0.5_jprb,zqsatm)
!         zqsatm = zqsatm/(1.0_jprb-retv*zqsatm)
        zqsatm = rslf(papm1(jl,jk),ptm1(jl,jk))
        zsatdef = max( 0._jprb, zqsatm-pqm1(jl,jk) )
        
        zpflxtot = zupflxl(jl,jk,jd) + zupflxn(jl,jk,jd)
        if (zpflxtot>0._jprb) then
         
          !-- precip evaporation tendency [kg/kg /s] (kessler 1969, tiedtke 1993) --
          zpevapup = zaprecevap * zsatdef * ( &
              & ( zpflxtot / 0.00509_jprb ) * &
              & ( papm1(jl,jk)/paphm1(jl,klev) )**0.5_jprb &
              & )**0.5777_jprb

          !rn testing: instantaneous evaporation
          !zpevapup = ( zupflxl(jl,jk,jd) + zupflxn(jl,jk,jd) ) / zdzrho

          !rn testing: no evap in subcloud layer
          !if (jk>=kplcl(jl,3)) then
          !  zpevapup = 0._jprb
          !endif
        
          !-- back-partition evaporation and substract from fluxes --
          zfac = zupflxl(jl,jk,jd) / zpflxtot
          zupflxl(jl,jk,jd) = zupflxl(jl,jk,jd) - zpevapup * zdzrho * zfac
          zupflxn(jl,jk,jd) = zupflxn(jl,jk,jd) - zpevapup * zdzrho * (1._jprb - zfac)
          zupflxl(jl,jk,jd) = max(0._jprb,zupflxl(jl,jk,jd))
          zupflxn(jl,jk,jd) = max(0._jprb,zupflxn(jl,jk,jd))
        endif
        
      enddo
    enddo
    
    !add contribution to total flux - weight by updraft area fraction
    do jk=0,klev
      do jl=kidia,kfdia
        pfplvl(jl,jk) = pfplvl(jl,jk) + zfracb(jl,jd) * zupflxl(jl,jk,jd)
        pfplvn(jl,jk) = pfplvn(jl,jk) + zfracb(jl,jd) * zupflxn(jl,jk,jd)
      enddo
    enddo
    
  enddo !jd
        
        

!     -----------------------------------------------------------------

!*         7.     construct mass flux profiles (jd=2,3)
!                 -------------------------------------


  !*         7.1  determine the mixed layer scaling height for jd=2,3
  !*
  do jl=kidia,kfdia
    
    select case (kpbltype(jl))
    
      case(1)
        !dry convective pbl - no moist parcel
        izi(jl,2) = kptop(jl,2)      !half level below level of zero kinetic energy
                
      case(2)
        !stratocumulus - no dry parcel
        izi(jl,3) = kplcl(jl,3)+1      !half level below lcl
	  
      case(3)
        !shallow cumulus - both dry and moist
        izi(jl,3) = kplcl(jl,3)+1    !half level below lcl
        izi(jl,2) = kptop(jl,2)      !half level below level of zero kinetic energy

      case(4)
        !deep cumulus - no moist parcel
        izi(jl,2) = kptop(jl,2)  

    end select 

  enddo !jl 



  !*         7.2  construct mixed layer mass fluxes
  !*               - use constant area fraction, and multiply by parcel w
  !*
  do jd = 2,kdraft

    do jk=klev-1,1,-1

      do jl=kidia,kfdia

        if ( kpbltype(jl)/=0 .and. zfracb(jl,jd)>0._jprb) then
   if (jk>=izi(jl,jd) ) then
	  
          zwuh = max( zwu2h(jl,jk,jd),0._jprb )
          zwuh = zwuh**0.5_jprb
          pmflx(jl,jk,jd)  = zfracb(jl,jd) * zwuh * zrhoh(jl,jk)
          
    if (zwu2h(jl,jk,jd)>0._jprb) then
      zfrac(jl,jk,jd) = zfracb(jl,jd)      
    else
      zfrac(jl,jk,jd) = 0._jprb
    endif
          
   endif
  endif
	
      enddo !jl
    
    enddo !jk
    
  enddo !jd
  
    
  
  !*         7.3    buoyancy sorting on moist updraft in cloud layer
  !*
  
  call vdfbuoysort( kidia     , kfdia   , klon    , klev   , kdraft , &
                  & papm1     , pgeom1  , pgeoh   , &
                  & zqtm1     , zslgm1  , &
                  & pfracb    , kplcl   , kptop   , kplzb  , &
                  & pqtuh     , pslguh  , zwu2h   , puuh   , pvuh   , &
                ! diagnostic output
                  & pextr2    , kfldx2  , pextra  , klevx  , kfldx  , &
                !              
                  & zabulk    , zwbulk  , zqtbulk , zslgbulk , zubulk , zvbulk )
  


  !*         7.4    construct cloudy mass flux profile (jd=3 only)
  !*
  do jk=klev-2,1,-1

    do jl=kidia,kfdia
   
      if ( kpbltype(jl)/=0 .and. kpbltype(jl)/=4 .and. zfracb(jl,3)>0._jprb ) then
        
        if (jk>=kptop(jl,1) .and. jk<=kplcl(jl,3) ) then

        zwuh = max(0._jprb,zwu2h(jl,jk,3))**0.5_jprb
        
        if( jk==kplcl(jl,3) ) then 

          !  special treatment for cloud base level

          ! -- moist updraft --
!          pmflx(jl,jk,3) = zfracb(jl,3)*zwuh*zrhoh(jl,jk) 
          pmflx(jl,jk,3) = min( zfracb(jl,3)*zwuh*zrhoh(jl,jk), 2._jprb*pmflx(jl,jk+1,3))   !cap acceleration through cloud base (can cause instability)
          zfrac(jl,jk,3) = zfracb(jl,3)
          
          ! -- dry updraft --
          !    tie dry updraft flux to its value at layer below, to prevent too sharp gradients
          !pmflx(jl,jk,2) = 0.3_jprb*pmflx(jl,jk+1,2)
          !zfrac(jl,jk,2) = 0.3_jprb*zfrac(jl,jk+1,2)
          pmflx(jl,jk,2) = 1.0_jprb*pmflx(jl,jk+1,2)
          zfrac(jl,jk,2) = 1.0_jprb*zfrac(jl,jk+1,2)
          !pmflx(jl,jk,2) = 0.0_jprb
          !zfrac(jl,jk,2) = 0.0_jprb
          
          pqtuh(jl,jk,2)  = 0.3_jprb*( pqtuh(jl,jk+1,2)  - zqtenh(jl,jk+1)  ) + zqtenh(jl,jk)
          pslguh(jl,jk,2) = 0.3_jprb*( pslguh(jl,jk+1,2) - zslgenh(jl,jk+1) ) + zslgenh(jl,jk)
          puuh(jl,jk,2)   = 0.3_jprb*( puuh(jl,jk+1,2)   - zuenh(jl,jk+1)   ) + zuenh(jl,jk)
          pvuh(jl,jk,2)   = 0.3_jprb*( pvuh(jl,jk+1,2)   - zvenh(jl,jk+1)   ) + zvenh(jl,jk)
          
  elseif ( pmflx(jl,jk+1,3) > 0._jprb ) then
          
          if ( jk>=kptop(jl,3) ) then
            
            !-- convective cloud layer --
            !   buoyancy sorting
            
            zfrac(jl,jk,3)  = zabulk(jl,jk)
            zwu2h(jl,jk,3)  = zwbulk(jl,jk) ** 2._jprb
            pmflx(jl,jk,3)  = zabulk(jl,jk) * zwbulk(jl,jk) * zrhoh(jl,jk)
            pqtuh(jl,jk,3)  = zqtbulk(jl,jk)
            pslguh(jl,jk,3) = zslgbulk(jl,jk)
            puuh(jl,jk,3)   = zubulk(jl,jk)
            pvuh(jl,jk,3)   = zvbulk(jl,jk)
            
          else 
            
            !-- inversion (layer between tops of test parcel and moist parcel) --
            !   prescribed linear decay of flux
            
            !zfrac(jl,jk,3)  = zfrac(jl,jk+1,3) * 0.25_jprb
            !pmflx(jl,jk,3)  = pmflx(jl,jk+1,3) * 0.25_jprb
            zfac = ( pzptop(jl,1) - pgeoh(jl,jk)*zrg  ) / ( pzptop(jl,1) - pzptop(jl,3) )
            zfac = max( 0._jprb, min(1._jprb,zfac) )
            zfrac(jl,jk,3)  = zfrac(jl,kptop(jl,3),3) * zfac
            pmflx(jl,jk,3)  = pmflx(jl,kptop(jl,3),3) * zfac
            
            pqtuh(jl,jk,3)  = ( pqtuh(jl,jk+1,3)  - zqtenh(jl,jk+1)  ) + zqtenh(jl,jk)
            pslguh(jl,jk,3) = ( pslguh(jl,jk+1,3) - zslgenh(jl,jk+1) ) + zslgenh(jl,jk)
            puuh(jl,jk,3)   = ( puuh(jl,jk+1,3)   - zuenh(jl,jk+1)   ) + zuenh(jl,jk)
            pvuh(jl,jk,3)   = ( pvuh(jl,jk+1,3)   - zvenh(jl,jk+1)   ) + zvenh(jl,jk)
            
            !zwu2h(jl,jk,3)  = zwu2h(jl,kptop(jl,3),3)
            zwu2h(jl,jk,3)  = zwu2h(jl,kptop(jl,3),3) * (zfac**2._jprb)
            !zwu2h(jl,jk,3)  = zwu2h(jl,jk,1)
            
          endif
	  
          !make sure that updraft #2 does not do any flux here
          pmflx(jl,jk,2)  = 0._jprb
          zfrac(jl,jk,2)  = 0._jprb
          pqtuh(jl,jk,2)  = 0._jprb
          pslguh(jl,jk,2) = 0._jprb
          puuh(jl,jk,2)   = 0._jprb
          pvuh(jl,jk,2)   = 0._jprb
          
  endif 
	
	
        !--- limit mass flux covering 50% area (m<rho*w,up*0.5) ---
        !    (detrainment is initiated if strong w,up slowdown)
        !inv: comment out
!inv        pmflx(jl,jk,3) =   min( pmflx(jl,jk,3) , 0.5_jprb * zwuh * zrhoh(jl,jk) )   

  endif 
        
      endif
     
    enddo !jl
     
  enddo !jk

  

  !*        7.5  compute i) variance transport flux and ii) bulk updraft detrainment, as used in vdfmain
  !*
  do jk=klev-1,1,-1
    do jl=kidia,kfdia
!      if (jk >= kptop(jl,1) ) then  
      if (jk > kptop(jl,1) ) then  
      
        !-- calculate w'qt'qt' at half levels --
        pwqt2(jl,jk) =  pmflx(jl,jk,3) * (pqtuh(jl,jk,3) - zqtenh(jl,jk))**2._jprb + &
                     &  pmflx(jl,jk,2) * (pqtuh(jl,jk,2) - zqtenh(jl,jk))**2._jprb

      endif 
    enddo !jl
  enddo !jk
  
  do jk=1,klev
    do jl=kidia,kfdia
      if ( kpbltype(jl)/=0 .and. kpbltype(jl)/=4 .and. zfracb(jl,3)>0._jprb ) then
        
!      if (jk >= kptop(jl,1) .and. jk<=kplcl(jl,3) ) then  
      if (jk > kptop(jl,1) .and. jk<=kplcl(jl,3) ) then  
	
        zdz   =  (   pgeoh(jl,jk-1)   - pgeoh(jl,jk)   ) * zrg 
        
  zdmdz = -(   pmflx(jl,jk-1,3) - pmflx(jl,jk,3) &
!           & + pmflx(jl,jk-1,2) - pmflx(jl,jk,2) &  
     & ) / zdz
                 
        if (jk==kptop(jl,1)+1) then
    zdmdz = pmflx(jl,jk,3) / zdz
        endif
        
        zdmdz = max( 0._jprb, zdmdz )   
           
        pdetr(jl,jk) = zdmdz + pmflx(jl,jk,3) / ( ztaueps * max(0.0001_jprb,zwu2h(jl,jk,3))**0.5_jprb )  
        
!        if (lldiag) then
!          pextra(jl,jk,98) = zdmdz
!          pextra(jl,jk,99) = pdetr(jl,jk) - zdmdz
!        endif  
	
      endif
      	
      endif 
    enddo !jl
  enddo !jk
          
  

  !*        7.6  estimate mass flux needed for cape removal (fritsch & chappell)
  !*
  
  if ( llcape ) then
  
  !-- calculate right-hand-side --
  do jk=klev-1,2,-1
    do jl=kidia,kfdia
      if (kpbltype(jl) ==3 .and. jk >= kptop(jl,3) .and. jk <= kplcl(jl,3) .and. pmflx(jl,kplcl(jl,3),3)>0._jprb ) then 
        zrhs(jl) = zrhs(jl) + &
      & ( pmflx(jl,jk,3) / pmflx(jl,kplcl(jl,3),3) ) * &
      & ( rg/zthven(jl,klev) ) * &
      & ( zthven(jl,jk-1)-zthven(jl,jk) ) 
      endif
    enddo !jl
  enddo !jk
  
  !-- calculate m_cape --
  do jl=kidia,kfdia
    zfracbcong(jl) = 0._jprb
    zmscale(jl)    = 1._jprb
    zmasscape(jl)  = 0._jprb
    if (zrhs(jl)>0._jprb) then
      zmasscape(jl)  = zcape1(jl) / (ztaubm * zrhs(jl))
      if (kplcl(jl,3).gt.0) then
        jk = kplcl(jl,3)
        zwuh = sqrt( max( 0._jprb, zwu2h(jl,jk,1) ) )
        if (zwuh > 0._jprb) then
          zfracbcong(jl) = max( 0._jprb, zmasscape(jl)-pmflx(jl,jk,3) ) / zwuh
        endif  
        if (pmflx(jl,jk,3) > 0._jprb) then
          !zmscale(jl) = zmasscape(jl)/pmflx(jl,jk,3)                  !pure cape closure
          zmscale(jl) = max(1._jprb, zmasscape(jl)/pmflx(jl,jk,3) )   !take the maximum of the two
        endif  
      endif  
    endif  
  enddo !jl

  do jk=klev-1,2,-1
    do jl=kidia,kfdia
      if (llcapetest) then
        ! option i: assign mass to the test updraft
        pmflx(jl,jk,1) = zfracbcong(jl) * sqrt( max( 0._jprb, zwu2h(jl,jk,1) ) )  
      else
        ! option ii: rescale moist mass flux
        pmflx(jl,jk,3) = zmscale(jl) * pmflx(jl,jk,3)
      endif  
    enddo !jl
  enddo !jk
  
  !-- diagnostics --
!  if ( lldiag ) then
!    do jl=kidia,kfdia
!      pextr2(jl,5) = zcape1(jl)
!      pextr2(jl,6) = zfracbcong(jl)
!      pextr2(jl,7) = zmasscape(jl)
!      pextr2(jl,8) = pmflx(jl,kplcl(jl,3),3)
!    enddo
!    do jk=1,klev
!      do jl=kidia,kfdia
!        pextra(jl,jk,98) = sqrt( max( 0._jprb, zwu2h(jl,jk,1) ) )
!        pextra(jl,jk,99) = pmflx(jl,jk,1)
!      enddo
!    enddo
!  endif

  endif



  !*        7.6  mass flux limit according to cfl criterion
  !*
 
  zcons10 = 1.0_jprb/(rg*ptmst)

!  do jd = 1,kdraft
!    do jl=kidia,kfdia
!      zmfs(jl,jd) = 1.0_jprb    ! default reduction factor
!    enddo
!  enddo
  
  do jd = 1,kdraft
    do jk=1,klev-1
      do jl=kidia,kfdia
        if ( jk >= kptop(jl,jd) .and. kptop(jl,jd)>0) then
          zmfmax = (papm1(jl,jk+1)-papm1(jl,jk)) * zcons10

!          pextra(jl,jk,6) = zmfmax

          if (jd==1) then
            zmfmax = 1.0_jprb * zmfmax
          else
            zmfmax = 3.0_jprb * zmfmax
          endif  
          
!          zmfmax = 4.0_jprb * zmfmax
!          zmfmax = 2.0_jprb * zmfmax
!          zmfmax = 1.0_jprb * zmfmax
!          zmfmax = 0.3_jprb * zmfmax

!          pextra(jl,jk,7) = zmfmax

          !option i: preserve vertical structure
!          if ( pmflx(jl,jk,jd) > zmfmax ) then
!            zmfs(jl,jd) = min(zmfs(jl,jd),zmfmax/pmflx(jl,jk,jd))
!          endif
          
          !option ii: correct level by level
    pmflx(jl,jk,jd) = min( zmfmax, pmflx(jl,jk,jd) )
	  
        endif
      enddo
    enddo
  enddo
  
!  do jd = 1,kdraft
!    do jk=1,klev
!      do jl=kidia,kfdia
!        pmflx(jl,jk,jd) = pmflx(jl,jk,jd)*zmfs(jl,jd)
!      enddo
!    enddo
!  enddo
          

  !rn sensitivity test: no dry updraft
!  do jk=1,klev
!    do jl=kidia,kfdia
!        pmflx(jl,jk,2) = 0.
!    enddo
!  enddo
          


  !*        7.7  cloud-layer average moist updraft mass flux.
  !*               (used in vdfexcu for estimating entrainment-k at cumulus pbl top)  
  !*
  do jk=klev-1,1,-1
    do jl=kidia,kfdia
      !-- cloudy mass flux: use moist updraft m (depth average) --
      if ( zqcuh(jl,jk,3)>0._jprb .and. pmflx(jl,jk,3)>0._jprb ) then  
        zmcld(jl) = zmcld(jl) + pmflx(jl,jk,3) * &
                &     zrg*( pgeoh(jl,jk-1) - pgeoh(jl,jk) )
      endif
    enddo
  enddo
  
  do jl=kidia,kfdia
    if (pzptop(jl,3)-pzplcl(jl,3)>0._jprb ) then  
!      pmcu(jl)    = zmcld(jl) / (pzptop(jl,3)-pzplcl(jl,3))
      pmcu(jl) = zmflxb(jl,3)
    endif
  enddo


    
  !*        7.8  scaling (see vdfexcu)
  !* 
  do jd = 1,kdraft
    do jk=1,klev-1
      do jl=kidia,kfdia
      
        zmgeom(jl,jk)  = pgeom1(jl,jk)-pgeom1(jl,jk+1)      
        zcfnc1(jl,jk)  = rvdifts * ptmst * rg**2 * paphm1(jl,jk) &
                 & /( zmgeom(jl,jk) * rd * 0.5_jprb &
                 & *( ptm1(jl,jk  )*(1.0_jprb+retv*pqm1(jl,jk  )) &
                 & +  ptm1(jl,jk+1)*(1.0_jprb+retv*pqm1(jl,jk+1))))  
        pmflx(jl,jk,jd) = zcfnc1(jl,jk) * zmgeom(jl,jk) * zrg * pmflx(jl,jk,jd)
        
      enddo
    enddo
  enddo



!     -----------------------------------------------------------------

!*         8.     set some w-scales for use in vdfmain
!                 ------------------------------------

  do jl=kidia,kfdia
    pwuavg(jl) = zwstar(jl)
    !pwuavg(jl) = 2.0_jprb * zwstar(jl)
  enddo

  !for use in buoyancy sorting scheme
  do jd = 1,kdraft
    do jk=0,klev
      do jl=kidia,kfdia
        pwuh(jl,jk,jd) = sqrt( max( 0._jprb, zwu2h(jl,jk,jd) ) )
      enddo
    enddo
  enddo
  
  

!     -----------------------------------------------------------------

!*         9.     advective flux adjustments at updraft top-level
!                 -----------------------------------------------
!
!                 vertical mixing at the top-level is prescribed
!                 and controlled in vdfexcu.
!

  do jd = 2,kdraft
    do jk=0,klev
      do jl=kidia,kfdia
        
        !  remove all mass flux at and above top layer of updraft
  itop = max( kptop(jl,jd), kptop(jl,3) )
	
        if ( jk <= itop ) then
          
          llwipe=.true.    
          
          ! protect cumulus cloud top: this is treated later in vdfexcu, dept. on llricu=t
          !   note: do wipe in case of single layer moist convection
          if ( kpbltype(jl)==3 .and. kptop(jl,3)<kplcl(jl,3)+1 ) then   
            llwipe=.false. 
          endif
          
          if (llwipe) then
            zfrac(jl,jk,jd)  = 0.0_jprb
            pmflx(jl,jk,jd)  = 0.0_jprb
            pwqt2(jl,jk)     = 0.0_jprb
          endif
          
        endif
        
      enddo
    enddo
  enddo
  
  
  
  !---------------- some output ------------------------
  if (lldiag) then
  
  do jl=kidia,kfdia
  
    !  boundary layer classification
    pextr2(pcog,30) = kpbltype(jl)

    !  updraft heights
    pextr2(pcog,31) = pzptop(jl,1)
    pextr2(pcog,32) = pzptop(jl,2)
    pextr2(pcog,33) = pzptop(jl,3)
    pextr2(pcog,34) = pzplcl(jl,1)
    pextr2(pcog,35) = pzplcl(jl,2)
    pextr2(pcog,36) = pzplcl(jl,3)
    
    !pextr2(jl,37) = pwuavg(jl)
    !pextr2(jl,38) = sqrt ( zwu2h(jl,klev-1,1) )
   
    !  updraft levels
    pextr2(pcog,24) = kplcl(jl,1)
    pextr2(pcog,25) = kplcl(jl,2)
    pextr2(pcog,26) = kplcl(jl,3)
    pextr2(pcog,27) = kptop(jl,1)
    pextr2(pcog,28) = kptop(jl,2)
    pextr2(pcog,29) = kptop(jl,3)
    
    !  various scalings
    !pextr2(jl,45) = 100._jprb*zfracb(jl,3)
        
    !  test updraft properties   
    zfrac(jl,0:kptop(jl,1)-1,1) = 0._jprb
    zfrac(jl,kptop(jl,1):klev-1,1) = 0.0001_jprb
    pextra(pcog,:,15) = 1000._jprb * ceiling(zfrac(jl,:,1)) * ( pqtuh(jl,:,1)  - zqtenh(jl,:) )
    pextra(pcog,:,16) = ceiling(zfrac(jl,:,1)) * ( pslguh(jl,:,1) - zslgenh(jl,:) )/rcpd 
    pextra(pcog,:,17) = max(0._jprb,zwu2h(jl,:,1))**0.5
    pextra(pcog,:,18) = zeps(jl,:,1)
    pextra(pcog,:,19) = 1000._jprb * zqcuh(jl,:,1)
    pextra(pcog,1:klev,20) = zbuof(jl,:,1)
    
    !pextr2(jl,13) = zcape1(jl)
    !pextr2(jl,18) = pmcu(jl)
    
    !  updraft precip generation
    pextra(pcog,:,21) = pfplvl(jl,:)
    pextra(pcog,:,22) = pfplvn(jl,:)

    !  updraft buoyancy
    pextra(pcog,1:klev,28) = zbuof(jl,:,2)
    pextra(pcog,1:klev,29) = zbuof(jl,:,3)

    !  updraft mass flux
    pextra(pcog,1:klev,30) = pmflx(jl,0:klev-1,2) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )  
    pextra(pcog,1:klev,31) = pmflx(jl,0:klev-1,3) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )
    
    !  updraft excesses
    pextra(pcog,:,32) = 1000._jprb * ceiling(zfrac(jl,:,2)) * ( pqtuh(jl,:,2)  - zqtenh(jl,:)  )
    pextra(pcog,:,33) = 1000._jprb * ceiling(zfrac(jl,:,3)) * ( pqtuh(jl,:,3)  - zqtenh(jl,:)  )
    pextra(pcog,:,34) = ceiling(zfrac(jl,:,2)) * ( pslguh(jl,:,2) - zslgenh(jl,:) )/rcpd 
    pextra(pcog,:,35) = ceiling(zfrac(jl,:,3)) * ( pslguh(jl,:,3) - zslgenh(jl,:) )/rcpd 
    pextra(pcog,:,36) = zwu2h(jl,:,2)
    pextra(pcog,:,37) = zwu2h(jl,:,3)
    pextra(pcog,:,26) = max(0._jprb,zwu2h(jl,:,2))**0.5
    pextra(pcog,:,27) = max(0._jprb,zwu2h(jl,:,3))**0.5

    !  updraft fractions
    pextra(pcog,:,38) = 100._jprb * zfrac(jl,:,2)
    pextra(pcog,:,39) = 100._jprb * zfrac(jl,:,3)
    
    !  updraft entrainment / detrainment
    pextra(pcog,:,40) = zeps(jl,:,2)
    pextra(pcog,:,41) = zeps(jl,:,3)
    
    !  updraft qt flux
    pextra(pcog,1:klev,44) =  rlvtt * ( pmflx(jl,0:klev-1,2) * (pqtuh(jl,0:klev-1,2) - zqtenh(jl,0:klev-1)) ) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )
    pextra(pcog,1:klev,45) =  rlvtt * ( pmflx(jl,0:klev-1,3) * (pqtuh(jl,0:klev-1,3) - zqtenh(jl,0:klev-1)) ) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )
    pextra(pcog,1:klev,46) =  rlvtt * ( pmflx(jl,0:klev-1,2) * (pqtuh(jl,0:klev-1,2) - zqtenh(jl,0:klev-1)) + &
                    &  pmflx(jl,0:klev-1,3) * (pqtuh(jl,0:klev-1,3) - zqtenh(jl,0:klev-1)) ) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg ) 

    !  updraft condensate
    pextra(pcog,:,47) = 1000._jprb * zqcuh(jl,:,2)
    pextra(pcog,:,48) = 1000._jprb * zqcuh(jl,:,3)
            
    !  updraft thl flux
    pextra(pcog,1:klev,51) =  ( pmflx(jl,0:klev-1,2) * (pslguh(jl,0:klev-1,2) - zslgenh(jl,0:klev-1)) ) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )
    pextra(pcog,1:klev,52) =  ( pmflx(jl,0:klev-1,3) * (pslguh(jl,0:klev-1,3) - zslgenh(jl,0:klev-1)) ) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )
    pextra(pcog,1:klev,53) =  ( pmflx(jl,0:klev-1,2) * (pslguh(jl,0:klev-1,2) - zslgenh(jl,0:klev-1)) + &
                    &    pmflx(jl,0:klev-1,3) * (pslguh(jl,0:klev-1,3) - zslgenh(jl,0:klev-1)) ) / (zcfnc1(jl,0:klev-1) * zmgeom(jl,0:klev-1) * zrg )

    !  qt and th(environment)
    pextra(pcog,:,42) = zqtenh(jl,:)
    pextra(pcog,:,43) = zslgenh(jl,:)

  enddo !jl
  
  endif !lldiag
!print*, 'pextra (2,:,54)',pextra (1,:,43)
!print*, 'zqtenh(1,:)',maxval(zslgenh(:,:))

! if (lhook) call dr_hook('vdfhghtn',1,zhook_handle)

end subroutine vdfhghtn
