!options xopt(hsfun)
subroutine vdfhghtn (kidia   , kfdia   , klon    , klev   , kdraft  , ptmst  , kstep    , &
                   & pum1    , pvm1    , pwm1    , ptm1    , pqm1   , plm1    , pim1   , pam1     ,&
                   & paphm1  , papm1   , pgeom1  , pgeoh  , pvervel , pqe    , pte      , &
                   & pkmfl   , pkhfl   , pkqfl   , pmflx  , pfrac   ,&
                   & pr_max  , ptkeint , &
                   & pextr2  , kfldx2  , pextra  , klevx  , kfldx , &
                   & puuh    , pvuh    , pslguh  , pqtuh  , pwuh  , &
                   & pzptop  , kptop   , pzplcl  , kplcl  , kplzb   , &
                   & pwuavg  , &
                   & pfracc  , pqcc    , pfplvl  , pfplvn  , &
                   & ldnodecp, ldrundry, kpbltype, pwqt2 , pqt2uh)
                   
                     
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
!     *pwm1*         w-velocity component at t-1                  m/s
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
!     *pr_max*       maximum cloud size

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
!                    1: dry cbl
!                    2: stratocumulus
!                    3: shallow cumulus
!                   (4: deep cumulus)

!     method
!     ------

!     see documentation

!     ------------------------------------------------------------------

! use garbage, only : foeewm
use thrm, only : rslf
use grid, only: edmfnfilter
use parkind1  ,only : jpim     ,jprb

! use yomhook   ,only : lhook,   dr_hook

use yos_cst   , only : rg       ,rd       ,rcpd     ,retv     ,rlvtt    ,&
                    & rlstt    ,ratm     ,rtt      ,rlmlt

use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les    ,&
                    & r4ies    ,r5les    ,r5ies    ,rvtmp2   ,r5alvcp  ,&
                    & r5alscp  ,ralvdcp  ,ralsdcp  ,rtwat    ,rtice    ,&
                    & rticecu  ,rtwat_rtice_r      ,rtwat_rticecu_r  
	      
use yoevdf   , only : rvdifts  ,lldiag
use yoecumf  , only : rtaumel
use yos_exc  , only : rkap     ,repust
!use yomgf1c  , only : nc
!use yomlog1c , only : lccn

implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
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
real(kind=jprb)   ,intent(in)    :: pwm1(klon,klev) 
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
!rn
real(kind=jprb)   ,intent(inout) :: pr_max(klon)
real(kind=jprb)   ,intent(inout) :: ptkeint(klon,kdraft)
!rn
real(kind=jprb)   ,intent(inout)   :: pmflx(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pfrac(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(out)     :: puuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(out)     :: pvuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pslguh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pqtuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pwuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pzplcl(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pzptop(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pwuavg(klon) 
real(kind=jprb)   ,intent(inout) :: pfplvl(klon,0:klev)
real(kind=jprb)   ,intent(inout) :: pfplvn(klon,0:klev)
real(kind=jprb)   ,intent(inout) :: pfracc(klon,0:klev)
real(kind=jprb)   ,intent(inout) :: pqcc  (klon,0:klev)
logical           ,intent(in)    :: ldnodecp(klon) 
!ldrundry not used now
logical           ,intent(in)    :: ldrundry(klon) 
integer(kind=jpim),intent(inout)   :: kpbltype(klon) 
real(kind=jprb)   ,intent(out)   :: pwqt2(klon,0:klev)
real(kind=jprb)   ,intent(out)   :: pqt2uh(klon,0:klev)
real(kind=jprb)   ,intent(in)    :: pvervel(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pte(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqe(klon,klev) 

!diagnostic output
integer(kind=jpim),intent(in) :: kfldx2, klevx, kfldx
real(kind=jprb)   ,intent(inout) :: pextr2(klon,kfldx2), pextra(klon,klevx,kfldx)



!*         0.2    local variables

!--- mean & environmental properties ---
real(kind=jprb) ::    zustar  (klon)     , zwstar(klon)       , zkhvfl(klon)       , &
                    & zusigma(klon)      , zwsigma(klon)      , &
                    & zslgenh(klon,0:klev),zqlenh(klon,0:klev), zqienh(klon,0:klev), &
                    & zqtenh(klon,0:klev), zuenh(klon,0:klev) , zvenh(klon,0:klev) , &
                    & ztven(klon,klev)   , zqtm1 (klon,klev)  , &
                    & zslgm1(klon,klev)  , zmgeom(klon,0:klev), &
                    & ztenh(klon,0:klev) , zrhoh (klon,0:klev), &
                    & zthven(klon,klev)  , zthvenh(klon,0:klev)

!--- updraft parameters ---
real(kind=jprb) ::    zwu2h (klon,0:klev,kdraft), zwuh, &
                    & zqcuh (klon,0:klev,kdraft), zquh  (klon,0:klev,kdraft), &
		    & ztuh  (klon,0:klev,kdraft), zeps  (klon,0:klev,kdraft), &
		    & zbuof (klon,klev,kdraft)  , &
                    & zzplzb(klon,kdraft)
                    
real(kind=jprb) ::    zabulk(klon,0:klev)   , zwbulk(klon,0:klev)    , zmbulk   (klon,0:klev), &
                    & zabulk_c(klon,0:klev) , zmbulk_c(klon,0:klev)  , zqcbulk_c(klon,0:klev), &
                    & zqtexbulk(klon,0:klev), zslgexbulk(klon,0:klev), &
                    & zuexbulk(klon,0:klev) , zvexbulk(klon,0:klev)
                    
real(kind=jprb) ::    zqsatm, zsatdef, &
                    & zupflxl(klon,0:klev,kdraft), zupflxn(klon,0:klev,kdraft), &
                    & zupgenl(klon,klev,kdraft), zupgenn(klon,klev,kdraft), &
                    & zalfaw, zdzrho, zpflxtot, zpevapup, zfac, zupmelt, &
                    & zaprecevap, zqtep

real(kind=jprb) ::    zmflxb(klon,kdraft)
      
real(kind=jprb) ::    zfracmax, zfacmaxexc , zfractest , zfactestexc , &
                    & zfacexc(klon,kdraft), zdumfrac, zdumr, &
                    & zpdffacphi(klon), zpdffacw(klon), zlobukhov

logical ::            lldone(klon,kdraft), llmcind(klon), & 
                    & llwipe, llstcu


integer(kind=jpim) :: is, jk, jl, jd, ibase, jkh, itop

real(kind=jprb) ::    zqexc   , ztexc   , zdz     ,  zdb   , &
                    & zspeedenv         , zspeedup, zwindir , &
                    & zz      , ztotw2(klon)      , ztotp(klon) , &
                    & zcons10 , zcfnc1(klon,0:klev)      , ztvmean     , &
                    & zrg     , zmfmax  , zmfs(klon,kdraft)

!          remaining model parameters

real(kind=jprb) ::    ztaueps , ztauepstemp(klon), zclddepth     , &
                    & zw2thresh               , zstabthresh    , zbirthresh , &
                    & ztvlim  , zclddepthdp   , zdzcloud(klon) , zw2h       , &
                    & zzfunc  , zzfunc3(klon) , zzfunc4(klon)  , zdhri(klon)    , zdhcl      , &
                    & zdthvdz        , zghm1
		    
real(kind=jprb) ::    zzi(klon)

integer(kind=jpim) :: izi(klon,kdraft)

real(kind=jprb) ::    zbuoycu, zdthvcutop

real(kind=jprb) ::    zmscale(klon)


!------- size densities --------

real(kind=jprb) ::   sd_r(klon,kdraft), sd_n(klon,kdraft), &
                   & zwthvint(klon,kdraft), zw2int(klon,kdraft), ztrans(klon,kdraft)
                   
real(kind=jprb) ::   sd_w    (klon,0:klev,kdraft), &
                   & sd_wmask(klon,0:klev,kdraft), &
                   & sd_cmask(klon,0:klev,kdraft)

integer(kind=jpim) :: jdmax(klon),i0,i1,ik0,ikb,ikt,ientrain,ihisto

logical :: llndensprescr

real(kind=jprb) ::   za1,za2,zb1,zb2, &
                   & zdr(klon),zntot(klon),zatot(klon), &
                   & zr_filter(klon), zr_break(klon)

!-----------------------------


		    
! real(kind=jprb) :: zhook_handle
! 
! 
! interface
! #include "surf_inq.h"
! end interface
! 
! 
! #include "vdfparcel.intfb.h"
! #include "vdfpdftable.intfb.h"
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

!-- total convective area fraction that is done with mass flux --
!zfracmax    = 0.2_jprb     
!zfracmax    = 0.075_jprb
zfracmax    = 0.1_jprb     
!zfracmax    = 0.15_jprb     
call vdfpdftable (zfracmax, zfacmaxexc, zdumr, zdumr, 0) ! associated pdf scaling factor

!-- top % of the pdf associated with the test parcel --
!zfractest   = 0.0002_jprb    
!zfractest   = 0.001_jprb    
zfractest   = 0.002_jprb      !cy32r3 
!zfractest   = 0.005_jprb    
call vdfpdftable (zfractest, zfactestexc, zdumr, zdumr, 0) ! associated pdf scaling factor

!-- switch for using prescribed cloud size density (powerlaw) --
llndensprescr = .true.
!llndensprescr = .false.

!-- power-law constants for prescribed cloud size density --
za1 =  ( 10._jprb**1.121_jprb ) / alog(10.)  !neggers et al.  , jas 2003
!zb1 = -1.70                                       !neggers et al.  , jas 2003
zb1 = -1.98_jprb                                  !benner and curry, jgr 1998

zb2 = -3.0_jprb                           !benner and curry, jgr 1998

!-- switch for discretization of histograms in size-space --
ihisto = 1   !option i: flexible bin-width, range up to maximum size 
!ihisto = 2   !option ii: fixed bin-width, range up to filter size

!-- parcel entrainment options --
!ientrain = 0   !default: epsilon = 1 / ( w_up * ztaueps )
ientrain = 1   !epsilon = 1 / cloud size

!-- eddy turnover time scale [s] used in parcel entrainment option 0  (neggers, siebesma & jonker, jas 2002) --
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

!call surf_inq(prepust=zrepust)
                 
!-- switch for applying klein-hartmann criterion for stratocumulus --
!llstcu = .true.
llstcu = .false.

!-- factor used in updraft initialization --
do jl=kidia,kfdia
  zpdffacw(jl)   = 1.0_jprb
  zpdffacphi(jl) = 1.0_jprb
enddo  

!-- updraft precip evaporation constant --
zaprecevap = 0.001_jprb
!zaprecevap = 0.000544_jprb
!zaprecevap = 0.0001_jprb

!-- optimization --
zrg    = 1.0_jprb/rg


!-- reset extra variables --
if (lldiag) then

  do jl=kidia,kfdia
    pextr2(jl,1:49) = 0._jprb
  enddo

  do jk=1,klevx
  do jl=kidia,kfdia
    pextra(jl,jk,49) = 0._jprb
    pextra(jl,jk,50) = 0._jprb
  enddo
  enddo

endif


!-- reset variables --
do jl=kidia,kfdia
  
  pwuavg(jl)     = 0.0_jprb
  kpbltype(jl)   = -1          ! -1 means: yet unknown
  
  zzi(jl)        = 0._jprb      !mixed layer scalings
  zwstar(jl)     = 0._jprb        
  
enddo


do jd=1,kdraft
  do jl=kidia,kfdia
    pzplcl(jl,jd)  = -100._jprb  ! default value: -100 (no lcl)
    pzptop(jl,jd)  = 0._jprb     
    kplcl(jl,jd)   = 0           ! default value: 0 (no pbl cloud)
    kptop(jl,jd)   = 0          
    kplzb(jl,jd)   = 0          
    lldone(jl,jd)  = .true.       ! default: true (don't launch the parcel)
    zfacexc(jl,jd) = 0._jprb 
    zfacexc(jl,jd) = 0._jprb 
    zmflxb(jl,jd)  = 0._jprb 
    izi(jl,jd)     = 0._jprb     
    
    sd_r(jl,jd)   = 0._jprb
    sd_n(jl,jd)   = 0._jprb
    
    zwthvint(jl,jd)   = 0._jprb   !integrated buoyancy flux
    zw2int(jl,jd)     = 0._jprb   !integrated plume-tke
    ztrans(jl,jd)     = 0._jprb   !transport
    
  enddo
enddo


do jk=0,klev
  do jl=kidia,kfdia
  
    pwqt2(jl,jk)  = 0._jprb  
    pqt2uh(jl,jk)  = 0._jprb  
    
    zabulk    (jl,jk) = 0._jprb  
    zwbulk    (jl,jk) = 0._jprb  
    zmbulk    (jl,jk) = 0._jprb  
    zqtexbulk (jl,jk) = 0._jprb  
    zslgexbulk(jl,jk) = 0._jprb  
    zuexbulk  (jl,jk) = 0._jprb  
    zvexbulk  (jl,jk) = 0._jprb  
    zabulk_c  (jl,jk) = 0._jprb  
    zmbulk_c  (jl,jk) = 0._jprb  
    zqcbulk_c (jl,jk) = 0._jprb  
    
    pfracc     (jl,jk) = 0._jprb
    pqcc       (jl,jk) = 0._jprb
    
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
    pfrac(jl,jk,jd) = 0.0_jprb
    ztuh(jl,jk,jd)  = 0.0_jprb
    zquh(jl,jk,jd)  = 0.0_jprb
    zqcuh(jl,jk,jd) = 0.0_jprb
    zeps(jl,jk,jd)  = 0.0_jprb
    zwu2h(jl,jk,jd) = 0.0_jprb
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



!--- set up size densities ---

do jl=kidia,kfdia

  !-- filter size --
  zr_filter(jl) = edmfnfilter
  !zr_filter(jl) = 100._jprb

  !-- scale break size --
  !write(0,*) 'vdfhghtn 1.1: prmax=',pr_max(jl),za1,zb1
  zr_break(jl) =             pr_max(jl)  !no scale break
  !zr_break(jl) = 0.75_jprb * pr_max(jl)  !scale break at 75 % of max size
  !zr_break(jl) = 200._jprb
  !zr_break(jl) = 600._jprb
  !zr_break(jl) = 1000._jprb

  !-- power-law constants for cloud size density --
  za2 = za1 * ( zr_break(jl)**zb1 ) / ( zr_break(jl)**zb2 )

  !-- discretization of the size densities: bin-width --
  select case (ihisto)
    case(1)  ! histogram option i: flexible bin-width, range up to maximum size 
      zdr(jl) = pr_max(jl) / (kdraft-1)
    case(2)  ! histogram option ii: fixed bin-width, range up to filter size
      zdr(jl) = zr_filter(jl) / (kdraft-1)
  end select

  if (lldiag) then
    pextr2(jl,2) = zdr(jl)
  endif  
  
  !-- construct size-array and determine jdmax (= index of largest bin to be calculated) --    
  jdmax(jl) = 1   !Note: plume 1 is a bulk plume - plume 2 is the first 
  do jd=2,kdraft
    sd_r(jl,jd)  = zdr(jl) /2. + (jd-2) * zdr(jl)
    !if (sd_r(jl,jd) .le. pr_max(jl)) then                       !only calculate bins of sizes < max size
    if ( sd_r(jl,jd) .le. min(pr_max(jl),zr_filter(jl)) ) then  !only calculate bins of sizes < max size AND filter size
      jdmax(jl) = jd 
    endif
  enddo
  !write(0,*) 'vdfhghtn 2.0: prmax=',pr_max(jl), ' zrfilter=',zr_filter(jl), ' jdmax=',jdmax(jl)

  !-- reset some size-densities --
  do jd=2,kdraft
    do jk=0,klev
      sd_w    (jl,jk,jd) = 0._jprb
      sd_wmask(jl,jk,jd) = 0._jprb
      sd_cmask(jl,jk,jd) = 0._jprb
    enddo
  enddo

enddo



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
      
    zthvenh(jl,jk)  =  ( zthven(jl,jk+1)  *(zghm1-pgeoh(jl,jk  )) &
                   & +   zthven(jl,jk)    *(pgeoh(jl,jk  )-pgeoh(jl,jk+1)) &
                   &   )                  /(zghm1-pgeoh(jl,jk+1))
                
  enddo
  enddo



!     -----------------------------------------------------------------

!*         3.     launch <kdraft> plumes to construct (z,l) fields
  
  
  do jl=kidia,kfdia
    
 
    !* 3.1    determine stability of bl using the surface buoyancy flux
    !*
    zkhvfl(jl)  = ( 1.0_jprb + retv *  zqtm1(jl,klev) ) * pkhfl(jl) + &
                & ( retv * zslgm1(jl,klev) / rcpd )     * pkqfl(jl) 


    !* 3.2    surface layer scaling
    !*
    if ( zkhvfl(jl) >= 0.0_jprb ) then
      
      ! stable bl (no updrafts expected/needed)
      kpbltype(jl) = 0

    else

      !convective bl
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
      
    endif

  enddo !jl
  
  
  
  do jd = 2,kdraft     !index 1 is reserved for ensemble-mean (bulk)

  !write(0,'(a,i3)' ) 'vdfhghtn: plume jd = ',jd

  !* 3.3  updraft initialization

  do jl=kidia,kfdia
    
    if ( zkhvfl(jl) < 0.0_jprb .and. jd.le.jdmax(jl) ) then
      
      !write(0,'(a,i,a,i)' ) 'vdfhghtn 3.0: ',jl,'  launching parcel ', jd
  
      lldone(jl,jd) = .false.  !confirm launch
     
      !-- set factor --
      !zfacexc(jl,jd) =              ( zfactestexc              /(kdraft-1)) * (jd - 1)   !between 0 and that of top-0.2%
      zfacexc(jl,jd) = zfacmaxexc + ((zfactestexc - zfacmaxexc)/(kdraft-1)) * (jd - 1)   !between that of top-10% and top-0.2%
      !zfacexc(jl,jd) = zfacmaxexc + (0.5/(kdraft-1)) * (jd - 1)
      
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

      puuh(jl,klev-1,jd)= zuenh(jl,klev-1)  !mean wind at this half level
      pvuh(jl,klev-1,jd)= zvenh(jl,klev-1)
      
      !write(0,'(a,i,f,a,f)') 'vdfhghtn 3.1: zfacexc:',jd,zfacexc(jl,jd),' zqexc:',zqexc
      
    endif

  enddo !jl


  !* 3.4   launch updraft
  
  select case (ientrain)
    case(0)
      do jl=kidia,kfdia
        ztauepstemp(jl) = ztaueps
      enddo !jl
    case(1)
      do jl=kidia,kfdia
        ztauepstemp(jl) = 1. / sd_r(jl,jd)
      enddo !jl
  end select 
  
  call vdfparcel (kidia   , kfdia   , klon    , klev    , kdraft  , ientrain, &
                & pgeoh   , pgeom1  , paphm1  , &
		& pum1    , pvm1    , pwm1    , zqtm1   , zslgm1  , ztven   , &
		& puuh    , pvuh    , pslguh  , pqtuh   , zwu2h   , zqcuh  , zbuof , & 
		& zquh    , ztuh    , zeps    , zfacexc , &
		& pzplcl  , kplcl   , pzptop  , kptop   , kplzb   , &
		& jd      , zupgenl , zupgenn , &
		& ztauepstemp , zw2thresh, lldone)  

  !do jl=kidia,kfdia
  !  write(0,'(a,i3,a,f15.7)' ) 'vdfhghtn:       jl=',jl,' pzptop=',pzptop(jl,jd)
    !do jk=kptop(jl,jd),klev-1
    !  write(0,'(a,i3,a,f10.3,a,f8.3,a,f8.3)' ) 'vdfhghtn:           jk=',jk,&
    !    & '  zh=',pgeoh(jl,jk)*zrg,'  wh=',max(0._jprb,zwu2h(jl,jk,jd))**0.5,&
    !    & '  qlh=',zqcuh(jl,jk,jd)*1000._jprb
    !enddo  !jk
  !enddo !jl

  
  enddo !jd


  !write(0,'(a)' ) 'vdfhghtn: plumes done'

      



!     -----------------------------------------------------------------

!*         4.     updraft precipitation
!                 ---------------------

  
  !*        4.1  updraft precipitation fluxes (rain and snow)
  !* 
  do jl=kidia,kfdia

    if (jdmax(jl).gt.1) then
      
    do jd = 2,jdmax(jl)
    
      do jk=2,klev

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

      enddo !jk

    enddo !jd
    
    endif  !jdmax>1
       
  enddo !jl
        
        

!     -----------------------------------------------------------------

!*         5.     integrate densities
!                 -------------------

  
  !* 5.1   prepare variables
  
  do jd = 2,kdraft
  do jl = kidia,kfdia
  do jk = 0,klev
    sd_w(jl,jk,jd) = max(0._jprb, zwu2h(jl,jk,jd))**0.5
    if (sd_w(jl,jk,jd).gt.0.) then
      sd_wmask(jl,jk,jd) = 1.
      if (zqcuh(jl,jk,jd).gt.0.) then
        sd_cmask(jl,jk,jd) = 1.
      endif  
    endif  
  enddo
  enddo
  enddo !jd
  
 
  !* 5.2  construct the number density N(l,z)
  
  do jl = kidia,kfdia
  
    !if ( zkhvfl(jl) < 0.0_jprb .and. jdmax(jl).gt.1) then
    if ( zkhvfl(jl) < 0.0_jprb ) then
    
    !do jd = 2,jdmax(jl)
    do jd = 2,kdraft   

      !calculate N(l,z) up to maximum size pr_max, not filter size zr_filter:
      !    this is important to get fractions right
      if ( sd_r(jl,jd).le.pr_max(jl) ) then

      !-- column-integrated properties of individual size bins --
      do jk = 0,klev
        if (jk.ge.2) then
          zdz  = ( pgeoh(jl,jk-1) - pgeoh(jl,jk) ) * zrg     
          zwthvint(jl,jd)   = zwthvint(jl,jd) + sd_w(jl,jk,jd) * max( 0._jprb, zbuof(jl,jk,jd)) * zdz
          zw2int  (jl,jd)   = zw2int  (jl,jd) + (sd_w(jl,jk,jd)**2)                             * zdz
        endif  
      enddo  !jk
      zw2int(jl,jd) = zw2int(jl,jd) * zdr(jl)
    
      if (llndensprescr) then
        
        !-- prescribe number using a powerlaw --
        if (sd_r(jl,jd).lt.zr_break(jl)) then
          sd_n(jl,jd)  = za1 * ( sd_r (jl,jd)**zb1 )
        else
          sd_n(jl,jd)  = za2 * ( sd_r (jl,jd)**zb2 )
        endif  

      else
        
        !-- number = total tkeint in bin / tkeint per plume --
        
        if ( jd.eq.jdmax(jl) .and. jd.gt.2 .and. ptkeint(jl,jd).eq.0._jprb) then
          ptkeint(jl,jd) = 0.1_jprb * ptkeint(jl,jd-1)   !in case of a new (i.e.larger) bin: supply a wee bit of initial tke
        endif
        
        sd_n(jl,jd)  = ptkeint(jl,jd) / zw2int(jl,jd)
        
      endif  !prescribed number density
      
      !write(0,'(2a,f,a,f,a,f,a,f,a,f)') 'vdfhghtn 4.0: ', &
      !                       & ' sd_r=',sd_r(jl,jd), &
      !                       & ' sd_n=',sd_n(jl,jd), &
      !                       & ' ptkeint=', ptkeint(jl,jd), &
      !                       & ' zw2int=',zw2int(jl,jd), &
      !                       & ' zwthvint=',zwthvint(jl,jd)

      endif  !sd_r<pr_max

    enddo  !jd
    
    endif  !wthv_s<0
    
  
    !-- integrate number density to get total number and coverage up to max size (pr_max) --
    zatot(jl) = 0._jprb
    zntot(jl) = 0._jprb
    !if (jdmax(jl).gt.1) then
    !  do jd=2,jdmax(jl)
      do jd=2,kdraft
        zntot(jl) = zntot(jl) + sd_n(jl,jd)                    * zdr(jl)
        zatot(jl) = zatot(jl) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * zdr(jl)
      enddo
    !endif
    !write(0,'(a,f15.7,a,f15.7)') 'vdfhghtn 4.1:  zntot=',zntot(jl),'  zatot=',zatot(jl)

  enddo  !jl
 
 
  !* 5.3  integrate densities with size
  
  !-- integrate --
  do jl = kidia,kfdia
  
    if ( zntot(jl) > 0.0_jprb .and. jdmax(jl).gt.1 ) then
    
    do jd = 2,jdmax(jl)

      do jk = 0,klev
        
        !- fractional quantities of individual bins -
        pfrac  (jl,jk,jd) =                     sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask (jl,jk,jd)
        pmflx  (jl,jk,jd) =                     sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_w     (jl,jk,jd)
  
        !- fractional quantities of bulk plume -
        zabulk    (jl,jk) = zabulk    (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask (jl,jk,jd)
        zmbulk    (jl,jk) = zmbulk    (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_w     (jl,jk,jd)
        zqtexbulk (jl,jk) = zqtexbulk (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask (jl,jk,jd) * ( pqtuh (jl,jk,jd) - zqtenh (jl,jk) )
        zslgexbulk(jl,jk) = zslgexbulk(jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask (jl,jk,jd) * ( pslguh(jl,jk,jd) - zslgenh(jl,jk) )
        zuexbulk  (jl,jk) = zuexbulk  (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask (jl,jk,jd) * ( puuh  (jl,jk,jd) - zuenh  (jl,jk) )
        zvexbulk  (jl,jk) = zvexbulk  (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask (jl,jk,jd) * ( pvuh  (jl,jk,jd) - zvenh  (jl,jk) )
      
        zabulk_c  (jl,jk) = zabulk_c  (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_cmask(jl,jk,jd)
        zmbulk_c  (jl,jk) = zmbulk_c  (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_cmask(jl,jk,jd) * sd_w (jl,jk,jd)
        zqcbulk_c (jl,jk) = zqcbulk_c (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_cmask(jl,jk,jd) * zqcuh(jl,jk,jd)
      
        pwqt2     (jl,jk) = pwqt2     (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_w    (jl,jk,jd) * ( pqtuh (jl,jk,jd) - zqtenh (jl,jk) )**2.
        pqt2uh    (jl,jk) = pqt2uh    (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * sd_wmask(jl,jk,jd) * pqtuh (jl,jk,jd)**2.
    
        pfplvl    (jl,jk) = pfplvl    (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * zupflxl(jl,jk,jd)
        pfplvn    (jl,jk) = pfplvn    (jl,jk) + sd_n(jl,jd) * (sd_r(jl,jd)**2) * zupflxn(jl,jk,jd)
      
        !if (lldiag) then
          !pextra(jl,1+jk,28+jd) = zabulk(jl,jk) * zdr(jl) * zfracmax / zatot(jl)
          !pextra(jl,1+jk,28+jd) = zmbulk(jl,jk) * zdr(jl) * zfracmax / zatot(jl)
          !pextra(jl,1+jk,28+jd) = pwqt2 (jl,jk) * zdr(jl) * zfracmax / zatot(jl)
          !pextra(jl,1+jk,28+jd) = sd_w  (jl,jk,jd)
          !pextra(jl,1+jk,28+jd) = zqcuh(jl,jk,jd)
          !pextra(jl,1+jk,28+jd) = zupflxl(jl,jk,jd)
        !endif  
  
      enddo  !jk
  
    enddo  !jd
  
    endif  !zntot>0, jdmax>1
    
  enddo  !jl

  
  !-- normalize --
  do jl = kidia,kfdia
  
    if ( zntot(jl) > 0.0_jprb .and. jdmax(jl).gt.1) then
    
    do jk = 0,klev

    !- fractional quantities of individual bins -
    do jd = 2,jdmax(jl)
      pmflx(jl,jk,jd) = pmflx(jl,jk,jd) * zdr(jl) * zfracmax / zatot(jl)
      pfrac(jl,jk,jd) = pfrac(jl,jk,jd) * zdr(jl) * zfracmax / zatot(jl)
    enddo
    
    if (zabulk(jl,jk).gt.0._jprb) then
    
      !- bulk averages -
      zwbulk    (jl,jk)  = zmbulk    (jl,jk) / zabulk(jl,jk)
      zqtexbulk (jl,jk)  = zqtexbulk (jl,jk) / zabulk(jl,jk)
      zslgexbulk(jl,jk)  = zslgexbulk(jl,jk) / zabulk(jl,jk)
      zuexbulk  (jl,jk)  = zuexbulk  (jl,jk) / zabulk(jl,jk)
      zvexbulk  (jl,jk)  = zvexbulk  (jl,jk) / zabulk(jl,jk)
      
      !- bulk variances among ensemble -
      pqt2uh(jl,jk)      = pqt2uh(jl,jk) / zabulk(jl,jk)  - (zqtexbulk (jl,jk) + zqtenh (jl,jk))**2.
      
    endif  
    
    !- cloud bulk averages -
    if (zabulk_c(jl,jk).gt.0._jprb) then
      zqcbulk_c (jl,jk)  = zqcbulk_c (jl,jk) / zabulk_c(jl,jk)
    endif  
    
    !- fractional quantities of bulk plume -
    zmbulk  (jl,jk) = zmbulk  (jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    zabulk  (jl,jk) = zabulk  (jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    pwqt2   (jl,jk) = pwqt2   (jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    pfplvl  (jl,jk) = pfplvl  (jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    pfplvn  (jl,jk) = pfplvn  (jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    
    zmbulk_c(jl,jk) = zmbulk_c(jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    zabulk_c(jl,jk) = zabulk_c(jl,jk) * zdr(jl) * zfracmax / zatot(jl)
    
    enddo !jk

    endif   !zntot>0, jdmax>1
    
  enddo !jl


  !* 5.4  time-integration of tkeint
  
  do jl = kidia,kfdia
  
    if (jdmax(jl).gt.1) then

    do jd = 2,jdmax(jl)
      if ( zntot(jl) > 0.0_jprb ) then
        zwthvint(jl,jd)   = zwthvint(jl,jd) * sd_n(jl,jd) * (sd_r(jl,jd)**2) * zdr(jl) * zfracmax / zatot(jl)
      endif  
      if (jd.lt.jdmax(jl)) then
        ztrans(jl,jd)     = ( ptkeint(jl,jd+1) - ptkeint(jl,jd) ) / 900._jprb
      else
        ztrans(jl,jd)     =                    - ptkeint(jl,jd)   / 900._jprb
      endif
    enddo
        
    do jd = 2,jdmax(jl)
      ptkeint(jl,jd) = ptkeint(jl,jd) + ptmst * ( zwthvint(jl,jd) + ztrans(jl,jd) )
    enddo
      
    endif  !jdmax>1
    
  enddo
  

  !- fill first index of histogram with bulk values -
  
  jd=1
  
  do jk = 0,klev
  do jl=kidia,kfdia

    pfrac (jl,jk,jd) = zabulk    (jl,jk)
    pmflx (jl,jk,jd) = zmbulk    (jl,jk)
    pwuh  (jl,jk,jd) = zwbulk    (jl,jk)
    pqtuh (jl,jk,jd) = zqtexbulk (jl,jk) + zqtenh  (jl,jk)
    pslguh(jl,jk,jd) = zslgexbulk(jl,jk) + zslgenh (jl,jk)
    puuh  (jl,jk,jd) = zuexbulk  (jl,jk) + zuenh   (jl,jk)
    pvuh  (jl,jk,jd) = zvexbulk  (jl,jk) + zvenh   (jl,jk)

    pfracc(jl,jk) = zabulk_c (jl,jk)     !area fraction covered by cloudy updrafts
    pqcc(jl,jk)   = zqcbulk_c(jl,jk)     !condensate averaged over cloudy updrafts

  enddo
  enddo


    
!     -----------------------------------------------------------------

!*         6.     postprocessing
!                 --------------

  
  !* 6.1   determine lcl and top heights of bulk updraft (for k-profile diffusion calculations in vdfexcu)
  !*              &
  !*       classify pbl type (kpbltype)
  
  jd=1
  
  do jl=kidia,kfdia
  
    if ( zntot(jl) == 0.0_jprb ) then
    
      kpbltype(jl)  = 0   !revert to non-convective state when size-density is zero
      
    else  
    
    ikt = klev+1
    ikb = -1
   
    do jk = 0,klev
      if (zabulk  (jl,jk).gt.0. .and. jk.lt.ikt) then
        ikt = jk
      endif  
      if (zabulk_c(jl,jk).gt.0. .and. jk.gt.ikb) then
        ikb = jk
      endif  
    enddo
    
    if (ikt.le.klev) then
    
      ! dry cbl
      kpbltype(jl)  = 1
      
      kptop (jl,jd) = ikt
      pzptop(jl,jd) = pgeoh(jl,ikt) * zrg
      
      if (ikb.eq.-1) then
        !ikb = max(ikb,ikt-1)   !set lcl at level above top
      else
        ! lcl found: shallow cumulus
        kpbltype(jl)  = 3
        kplcl (jl,jd) = ikb
        pzplcl(jl,jd) = pgeoh(jl,ikb) * zrg
      endif
      
    endif  
  
    !write(0,'(a,i,f,i,f,i,f)') 'vdfhghtn:  5.0:',jk, zkhvfl(jl), kptop(jl,jd), pzptop(jl,jd), kplcl(jl,jd), pzplcl(jl,jd)
    
    endif    !zntot>0
     
  enddo
  
  
  !* 6.2   w-scale for use in qt-variance budget (used in vdfmain)

  do jl=kidia,kfdia
    pwuavg(jl) = zwstar(jl)
    !pwuavg(jl) = 2.0_jprb * zwstar(jl)
  enddo


  !* 6.5   scaling (for adv-diff solver, see vdfdifh & vdfdifm)

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
  
  
  !* 6.6   update maximum size
  
  do jl=kidia,kfdia
  
    if (jdmax(jl).gt.1) then

    select case (ihisto)

      case(1)  ! histogram option i: flexible bin-width, range up to maximum size 
        !write(0,'(a,f,a,2f,a,f)') 'vdfhghtn 6.0: zdr:',zdr(jl),' ztop:',pzptop(jl,1), pzptop(jl,jdmax(jl)), ' prmax:', pr_max(jl)
        !pr_max(jl) = max( 1000._jprb, pzptop(jl,1) )
        pr_max(jl) = max( 1000._jprb,            pzptop(jl,jdmax(jl)) )
        !pr_max(jl) = 1000._jprb

      case(2)  ! histogram option ii: fixed bin-width, range up to filter size
        pr_max(jl) = min( zr_filter(jl), max( zdr(jl), 0.5_jprb*pzptop(jl,jdmax(jl))) )

    end select
    
    endif !jdmax
    
  enddo !jl



  !---------------- some output ------------------------
  if (lldiag) then
  
  do jl=kidia,kfdia

    if (jdmax(jl).gt.1) then

    do jd=1,jdmax(jl)
      if (88+jd<=kfldx) then
        pextra(jl,:,88+jd) = sd_w(jl,:,jd)
        !pextra(jl,:,88+jd) = zqcuh(jl,:,jd)
        !pextra(jl,:,88+jd) = pqtuh(jl,:,jd)
        !pextra(jl,:,88+jd) = pmflx(jl,:,jd)
        !pextra(jl,:,88+jd) = zupflxl(jl,:,jd)
        pextra(jl,:,99+jd) = pfrac(jl,:,jd)
      endif 
    enddo
    
    do jd=2,jdmax(jl)
      
      if (klev-kdraft+jd>0) then
        pextra(jl,klev-kdraft+jd, 8)   = sd_r    (jl,jd)
        pextra(jl,klev-kdraft+jd, 9)   = sd_n    (jl,jd)
        pextra(jl,klev-kdraft+jd,11)   = zwthvint(jl,jd)
        pextra(jl,klev-kdraft+jd,12)   = ztrans  (jl,jd)
        pextra(jl,klev-kdraft+jd,13)   = zw2int(jl,jd)
        if (zw2int(jl,jd).gt.0._jprb) then
          pextra(jl,klev-kdraft+jd,14)   = ptkeint(jl,jd) / zw2int(jl,jd)
        endif  
        pextra(jl,klev-kdraft+jd,15)   = pzplcl(jl,jd)
        pextra(jl,klev-kdraft+jd,16)   = pzptop(jl,jd)
      endif 
       
    enddo
    
    !pextra(jl,1:klev,88) = zbuof(jl,:,kdraft)
      
    pextra(jl,:,50) = zabulk_c  (jl,:)
    pextra(jl,:,51) = zmbulk_c  (jl,:)
    pextra(jl,:,52) = zqcbulk_c (jl,:) * 1000._jprb   !in g/kg
    
    !pextra(jl,:,79) = 1000000._jprb * pqt2uh(jl,:)
    pextra(jl,:,80) = zabulk    (jl,:)
    pextra(jl,:,81) = zwbulk    (jl,:)
    pextra(jl,:,82) = zmbulk    (jl,:)
    pextra(jl,:,83) = zqtexbulk (jl,:)
    pextra(jl,:,84) = zslgexbulk(jl,:)
    pextra(jl,:,85) = zuexbulk  (jl,:)
    pextra(jl,:,86) = zvexbulk  (jl,:)
    pextra(jl,:,87) = pfplvl  (jl,:)
    
    endif  !jdmax>1

  enddo !jl
  
  endif !lldiag
  


  !write(0,'(a)' ) 'vdfhghtn: end'


! if (lhook) call dr_hook('vdfhghtn',1,zhook_handle)

end subroutine vdfhghtn
