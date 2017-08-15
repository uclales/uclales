!options xopt(hsfun)
subroutine vdfparcel (kidia   , kfdia   , klon    , klev    , kdraft  , kentrain, &
                    & pgeoh   , pgeom1  , paphm1  , &
                    & pum1    , pvm1    , pwm1    , pqtm1   , pslgm1  , ptven   , &
                    & puuh    , pvuh    , pslguh  , pqtuh   , pwu2h   , pqcuh  , pbuof , & 
                    & pquh    , ptuh    , peps    , pfacexc , &
                    & pzplcl  , kplcl   , pzptop  , kptop   , kplzb   , &
                    & kd      , pupgenl , pupgenn , &
                    & ptaueps , pw2thresh, lddone )  
!     ------------------------------------------------------------------

!**   *vdfparcel* - vertical integration for parcel ascent
!
!     based on original vdfhghtn.f90 (cy29r1 and earlier) by
!             a.p. siebesma    30/06/99  
!             m. ko"hler       3/12/2004 
!             roel neggers     12/04/2005     separated from vdfhghtn
!                              15/10/2005     pressure term added
!                              12/04/2006     debugged interpolation of lcl height
!                              30/11/2006     updraft precipitation generation added
!                                             level of zero buoyancy determination


!     purpose
!     -------

!     determine pbl height and updraft fields

!     interface
!     ---------

!     *vdfparcel* is called by *vdfhghtn*

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
!                    (4: downdrafts .. to be done)
!     *kd*           draft index



!     input parameters (real):

!     *pgeoh*        geopotential at half level                   m2/s2
!     *pgeom1*       geopotential at t-1                          m2/s2
!     *paphm1*       pressure at half level at t-1                pa
!     *pum1*         x-velocity component at t-1                  m/s
!     *pvm1*         y-velocity component at t-1                  m/s
!     *pwm1*         w-velocity component at t-1                  m/s
!     *pqtm1*        mean specific total water at full level      kg/kg
!     *pslgm1*       mean liquid static energy at full level      m2/s2
!     *ptven*        environmental virtual temperature            k
!     *ptaueps*      updraft entrainment timescale                s
!     *pw2thresh*    threshold updraft velocity squared           m2/s2


!     input parameters (logical):
!     *lddone*       parcel liftoff confirmation (.true. means 'don't launch parcel')


!     output parameters (real):

!     *puuh*         updraft x-momentum
!     *pvuh*         updraft y-momentum
!     *pslguh*       updraft generalized liquid static energy (slg)
!                    at half level                                   m2/s2
!     *pqtuh*        updraft total specific humidity at half level   kg/kg
!     *pwu2h*        updraft vertical velocity square at half level  m2/s2
!     *pqcuh*        updraft liquid water at half level              kg/kg
!     *pquh*         updraft specific humidity at half level         kg/kg
!     *ptuh*         updraft temperature at half level               k
!     *pbuof*        updraft buoyancy at full level                  m/s2
!     *peps*         updraft entrainment rate                        1/m
!     *pzplcl*       height of lifting condensation level of updraft          m
!     *pzptop*       height of level of zero kinetic energy (w=0) of updraft  m
!     *pupgenl*      updraft rain generation                         kg/kg /s
!     *pupgenn*      updraft snow generation                         kg/kg /s
!
!     output parameters (integer):

!     *kplcl*        first half level above real height of upraft lcl
!     *kptop*        highest half level below pztop, and
!                    updraft top full level (pztop is within that layer)
!     *kplzb*        level of upraft zero buoyancy (highest full level that is pos. buoyant)

!     method
!     ------

!     see documentation

!     ------------------------------------------------------------------

!#include "tsmbkind.h"
use garbage, only : foealfa, satadj
use parkind1  ,only : jpim     , jprb

! use yomhook   ,only : lhook,   dr_hook

use yoethf   , only : r2es     , r3les    , r3ies    , &
                     &r4les    , r4ies    , r5les    , r5ies     , rvtmp2  , & 
                     &ralvdcp  , ralsdcp  , &
		     &rtwat    , rtice    , rticecu  , r5alvcp   , r5alscp , & 
                     &rtwat_rtice_r       , rtwat_rticecu_r
                     
use yos_cst   , only : rg  , rlstt , rcpd , rlvtt , retv ,rtt         



implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: kdraft
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim),intent(in)    :: kd
integer(kind=jpim),intent(inout) :: kplcl(klon,kdraft)
integer(kind=jpim),intent(inout) :: kptop(klon,kdraft)
integer(kind=jpim),intent(inout) :: kplzb(klon,kdraft)
integer(kind=jpim),intent(in)    :: kentrain
real(kind=jprb)   ,intent(in)    :: pgeoh(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: pgeom1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: paphm1(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: pum1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pvm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pwm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqtm1 (klon,klev)
real(kind=jprb)   ,intent(in)    :: pslgm1(klon,klev)
real(kind=jprb)   ,intent(in)    :: ptven(klon,klev)
real(kind=jprb)   ,intent(inout)   :: puuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pvuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pslguh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pqtuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pwu2h(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pqcuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pquh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: ptuh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: peps(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(in)      :: pfacexc(klon,kdraft)
real(kind=jprb)   ,intent(inout)   :: pzptop(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pzplcl(klon,kdraft) 
real(kind=jprb)   ,intent(in)      :: ptaueps(klon) 
real(kind=jprb)   ,intent(in)      :: pw2thresh
logical           ,intent(inout)   :: lddone(klon,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pbuof(klon,klev,kdraft) 
real(kind=jprb)   ,intent(inout)   :: pupgenl(klon,klev,kdraft)
real(kind=jprb)   ,intent(inout)   :: pupgenn(klon,klev,kdraft)

!*         0.2    local variables

logical ::            lldoit(klon), llcloud(klon), llprec

real(kind=jprb) ::    zrg   , zdz   , ztvuf , zquf    , zqcuf   , zmu, zb, &
                    & zqluf , zqiuf , zqtuf , ztuf  , zslguf  , &
		    & zmix  (klon,0:klev), zmixw (klon,0:klev) !,zbuof(klon,klev)

real(kind=jprb) ::    zwuh  (klon,0:klev)   , zquh  (klon,0:klev)    , zph  (klon), &
                    & zttemp(klon,klev)     , zqtemp(klon,klev)      
		    
real(kind=jprb) ::    zalfaw  , zfacw   , zfaci   , zfac    , ztemp ,&
                    & zesdp   , zcor    , zdqsdtemp(klon)   , zqs(klon,0:klev), &
                    & zpgenup , zwuhtemp(klon), zcepsz(klon)

integer(kind=jpim) :: is, jk, jl, jkmax, jkmin

real(kind=jprb) ::    zdqsudz, zdqtudz, zlclfac(klon), zepscflfac, zqlwork, zlcrit

! real(kind=jprb) ::    zhook_handle

!dir$ vfunction exphf
! #include "fcttre.h"
! #include "cuadjtq.intfb.h"


!     -----------------------------------------------------------------

!*         1.     set some constants
!                 --------------------

! if (lhook) call dr_hook('vdfparcel',0,zhook_handle)

!  constant of proportionality between updraft induced pressure term and
!     updraft vertical acceleration (siebesma, soares and teixeira, jas 2007)
!zmu = 0._jprb  
zmu = 0.15_jprb   !cy32r3

!  constant of proportionality between kinematic and thermodynamic entrainment
!zb = 1.0_jprb
zb = 0.5_jprb   !cy32r3

!  cfl criterion factor for updraft lateral entrainment
zepscflfac = 0.6_jprb   
!zepscflfac = 1.0_jprb

!  optimization
zrg         = 1.0_jprb/rg 

!  updraft precipitation switch
llprec=.true.
!llprec=.false.

!  critical updraft condensate [kg/kg] in sundqvist precipitation generation
!zlcrit = 0.0005_jprb  
zlcrit = 0.001_jprb     !cy32r3
!zlcrit = 0.0015_jprb   



!     -----------------------------------------------------------------

!*         2.     some final initialization
!                 ------------------------

  do jl=kidia,kfdia
    llcloud(jl)        = .false.
    pbuof (jl,klev,kd) = 0.0_jprb
    zlclfac(jl)        = 0.0_jprb
    kplzb(jl,kd)       = klev-1
    
    !  1/z scaling factor in eps
    zcepsz(jl) = 0._jprb
    !if (kd==2 .and. kpbltype(jl)==1) then
    !  zcepsz(jl) = 0.4_jprb
    !endif
    !if (.not.lddone(jl,kd)) then
    !  zcepsz(jl) = ( 2._jprb * pfacexc(jl,kd) )**(-1._jprb)
    !endif  
    !write(0,'(a,i,2f)') '    ',kd,pfacexc(jl,kd),zcepsz(jl)
    
  enddo
    
  !integration over total depth for now: 
  !    note that limiting jkmin to kplcl/kptop for kd=2 could speed things up a bit..
  jkmax = klev-2
  jkmin = 1


!     -----------------------------------------------------------------

!*         3.     vertical ascent until velocity becomes negative
!                 -----------------------------------------------
  
  do jk=jkmax,jkmin,-1
  
    is=0
    do jl=kidia,kfdia
      if (.not.lddone(jl,kd)) then
        is            = is+1


!*         3.1  updraft entrainment
        
        select case (kentrain)
        
          case(0)
            zwuh(jl,jk+1) = sqrt( max( pwu2h(jl,jk+1,kd), 0.01_jprb) ) ! w,up > 0.1 m/s (safety)
            peps(jl,jk+1,kd) = 1.0_jprb / ( zwuh(jl,jk+1) * ptaueps(jl) ) !& ! eps=1/(w,up*tau)
        
          case(1)
            peps(jl,jk+1,kd) = ptaueps(jl)
          
        end select 
        
        !rn numerical entrainment limiter: maximally c_e/dz
        zdz           = (pgeoh(jl,jk) - pgeoh(jl,jk+1))*zrg
        peps(jl,jk+1,kd) = min( zepscflfac/zdz, peps(jl,jk+1,kd) )
        peps(jl,jk+1,kd) = max( peps(jl,jk+1,kd), zcepsz(jl) * rg/pgeoh(jl,jk+1) )
        

!*         3.2  ascent of slg and qt (exact)
 
        zmix(jl,jk+1) = exp( - zdz * peps(jl,jk+1,kd) )
        pqtuh(jl,jk,kd)  = ( pqtuh (jl,jk+1,kd) - pqtm1 (jl,jk+1) ) * zmix(jl,jk+1) &
                    & + pqtm1 (jl,jk+1)
        pslguh(jl,jk,kd) = ( pslguh(jl,jk+1,kd) - pslgm1(jl,jk+1) ) * zmix(jl,jk+1) &
                    & + pslgm1(jl,jk+1)
        puuh(jl,jk,kd)   = ( puuh (jl,jk+1,kd)  - pum1  (jl,jk+1) ) * zmix(jl,jk+1) &
                    & + pum1 (jl,jk+1)
        pvuh(jl,jk,kd)   = ( pvuh (jl,jk+1,kd)  - pvm1  (jl,jk+1) ) * zmix(jl,jk+1) &
                    & + pvm1 (jl,jk+1)
                    

!*         3.3  condensation - diagnose t, qv, ql
!*         (taylor, becomes more inaccurate for states far from q*
!*          -> double condensation step)

!          cuadjtq initialization (assume qv=qt)

        pquh(jl,jk,kd)   = pqtuh(jl,jk,kd)
        pqcuh(jl,jk,kd)  = 0.0_jprb

!          cuadjtq initialization (assume qv=qv(jk+1) for speed)

!       if ( zquh(jl,jk+1) < pqtuh(jl,jk) ) then
!         zquh(jl,jk) = zquh(jl,jk+1)
!         zqcuh(jl,jk)= pqtuh(jl,jk) - zquh(jl,jk+1)
!       endif

        ptuh(jl,jk,kd)   = ( pslguh(jl,jk,kd) - pgeoh(jl,jk) + rlvtt*pqcuh(jl,jk,kd) ) &
                    & / rcpd      ! assume liquid phase!
!         zph(jl)       = paphm1(jl,jk)
!         zqtemp(jl,jk) = pquh(jl,jk,kd)
!         zttemp(jl,jk) = ptuh(jl,jk,kd)
! 	
      endif
      
      lldoit(jl)      = .not. lddone(jl,kd)
      
      if ( kd==2 ) then    ! condensation not done for dry subcloud thermal
        lldoit(jl)    = .false.
      endif
      if (lldoit(jl)) then
! print *, jl, jk, kd , paphm1(jl,jk)     
        call satadj(ptuh(jl,jk,kd), paphm1(jl,jk), pqtuh(jl,jk,kd), pqcuh(jl,jk,kd))
        pquh(jl,jk,kd) = pqtuh(jl,jk,kd) - pqcuh(jl,jk,kd)
      end if
      
    enddo

! 
!     call cuadjtq &
!      & ( kidia,    kfdia,    klon,     0,       klev,&
!      &   jk,&
!      &   zph,      zttemp,   zqtemp,   lldoit,  4)  
! 
! 
!     do jl=kidia,kfdia
!       if ( lldoit(jl) ) then
!         if ( zqtemp(jl,jk) < pqtuh(jl,jk,kd) ) then !allow evaporation up to qt
!           pquh(jl,jk,kd) = zqtemp(jl,jk)
!           pqcuh(jl,jk,kd)= pqtuh(jl,jk,kd) - pquh(jl,jk,kd)
!           ptuh(jl,jk,kd) = zttemp(jl,jk)
!         else                          !case where qv(initial)<qt but qv(final)>qt
!           pquh(jl,jk,kd) = pqtuh(jl,jk,kd)  !(unusual!)
!           pqcuh(jl,jk,kd)= 0.0_jprb
!           ptuh(jl,jk,kd) = ( pslguh(jl,jk,kd) - pgeoh(jl,jk) + rlvtt*pqcuh(jl,jk,kd) ) &
!                     & / rcpd
!         endif
!       endif
!     enddo


    do jl=kidia,kfdia


!*         3.4  updraft microphysics                *experimental* rn 
!*
      !precip generation tendency (sundqvist,1978) 
      !    [kg/kg /s]    in full level below current half level
      if (llprec .and. llcloud(jl)) then

        zdz = (pgeoh(jl,jk) - pgeoh(jl,jk+1)) * zrg
        zqlwork = pqcuh(jl,jk,kd)
        !zqlwork = pqcuh(jl,jk+1,kd)
        !zqlwork = ( pqcuh(jl,jk,kd) + pqcuh(jl,jk+1,kd) ) / 2._jprb
        zpgenup = 0.0015_jprb * zqlwork * (1._jprb - exp(-(zqlwork/zlcrit)**2) )
        zpgenup = min( zpgenup, pqcuh(jl,jk,kd)/zdz )              

        zalfaw = foealfa(ptuh(jl,jk,kd))
        pupgenl(jl,jk+1,kd) = zpgenup * zalfaw
        pupgenn(jl,jk+1,kd) = zpgenup * (1._jprb - zalfaw)
        
        !adjust the associated updraft state variables (integrate tendency over layer)      
        pqcuh(jl,jk,kd)  = pqcuh(jl,jk,kd)  - zpgenup * zdz
        pqtuh(jl,jk,kd)  = pqtuh(jl,jk,kd)  - zpgenup * zdz
        pslguh(jl,jk,kd) = pslguh(jl,jk,kd) + zdz * &
           & (rlvtt * pupgenl(jl,jk+1,kd) + rlstt * pupgenn(jl,jk+1,kd) )

      endif


!*         3.5  interpolation of updraft lcl

      if ( pqcuh(jl,jk,kd) > 0.0_jprb  .and.  .not. llcloud(jl) ) then

        llcloud(jl)   = .true.
	
	!cloud base level is first level with ql>0
        kplcl(jl,kd)  = jk

        !rn --- new interpolation method incorporating dqt/dz *and* dqsat/dz ---
        zdqsudz = ( zqtemp(jl,jk)   - zqtemp(jl,jk+1)   ) * rg / ( pgeoh(jl,jk) - pgeoh(jl,jk+1) )
        zdqtudz = ( pqtuh(jl,jk,kd) - pqtuh(jl,jk+1,kd) ) * rg / ( pgeoh(jl,jk) - pgeoh(jl,jk+1) )
        
        pzplcl(jl,kd) = pgeoh(jl,jk)*zrg 
        if (zdqsudz-zdqtudz.lt.0._jprb) then
          pzplcl(jl,kd) = pgeoh(jl,jk)*zrg + pqcuh(jl,jk,kd)/(zdqsudz-zdqtudz)
        endif
        
        if ( kplcl(jl,kd) < klev ) then
          pzplcl(jl,kd) = max( pzplcl(jl,kd), pgeoh(jl,kplcl(jl,kd)+1)*zrg )
        else
          pzplcl(jl,kd) = max( pzplcl(jl,kd), 0.0_jprb )
        endif
        
        zlclfac(jl) = ( pgeoh(jl,jk) - pzplcl(jl,kd)*rg  ) / ( pgeoh(jl,jk) - pgeoh(jl,jk+1) )
        zlclfac(jl) = max(0._jprb, min(1._jprb, zlclfac(jl)) )

      endif


!*         3.6  updraft buoyancy
!*         (at full level k+1 from interpolation of slg, q, qc)

      if ( .not. lddone(jl,kd) ) then

        zslguf        = 0.5_jprb * ( pslguh(jl,jk,kd) + pslguh(jl,jk+1,kd) )
        zqtuf         = 0.5_jprb * ( pqtuh (jl,jk,kd) + pqtuh (jl,jk+1,kd) )
        zqcuf         = 0.5_jprb * ( pqcuh (jl,jk,kd) + pqcuh (jl,jk+1,kd) )
	
	! at first full level above real cloud base, interpolate ql between 
	! height of real cloud base and the half level of kplcl
        if ( jk == kplcl(jl,kd) ) then
          zqcuf = zqcuf * zlclfac(jl)
        endif
        
	zquf          = zqtuf - zqcuf
        ztuf          = ( zslguf - pgeom1(jl,jk+1) & ! preliminary estimate:
                    & + rlvtt * zqcuf ) / rcpd       ! all liquid 
        zalfaw        = foealfa( ztuf )
        zqluf         = zalfaw            * zqcuf
        zqiuf         = (1.0_jprb-zalfaw) * zqcuf
        ztuf          = ( zslguf - pgeom1(jl,jk+1) &
                    & + rlvtt * zqluf + rlstt * zqiuf ) / rcpd  
        ztvuf         = ztuf * ( 1.0_jprb + retv * zquf - zqcuf )  
        pbuof(jl,jk+1,kd) = rg * ( ztvuf - ptven(jl,jk+1) ) / ptven(jl,jk+1)
        
        !update level of zero buoyancy
        if ( pbuof(jl,jk+1,kd)>0._jprb ) then
          kplzb(jl,kd) = jk+1
        endif


!*         3.7  kinetic energy equation (exact)
        
        zdz           = (pgeoh(jl,jk) - pgeoh(jl,jk+1))*zrg
        ztemp         = pbuof(jl,jk+1,kd) / (zb*peps(jl,jk+1,kd))
!        pwu2h(jl,jk,kd)  = ( pwu2h(jl,jk+1,kd) - ztemp ) * zmix(jl,jk+1)**2 + ztemp
        zmixw(jl,jk+1) = exp( - zdz * (zb*peps(jl,jk+1,kd)) / (1._jprb - 2._jprb*zmu) )
!        pwu2h(jl,jk,kd)  = ( pwu2h(jl,jk+1,kd) - ztemp ) * zmixw(jl,jk+1)**2 + ztemp
        !RN maintaining the gridbox-average vertical velocity:
        pwu2h(jl,jk,kd)  = ( pwu2h(jl,jk+1,kd) - pwm1(jl,jk+1)**2 - ztemp ) * zmixw(jl,jk+1)**2 + ztemp + pwm1(jl,jk+1)**2
        

!*         3.8  inversion height at w=0  (lin. interpolation in w^2)
        
        if ( pwu2h(jl,jk,kd) < 0.0_jprb  .and.  pwu2h(jl,jk+1,kd) > 0.0_jprb ) then 
	
	  !set top level to last level with positive w
          kptop(jl,kd)   = jk+1  
	  
          pzptop(jl,kd)   = pgeoh(jl,jk+1) * zrg &
                    & + zdz * pwu2h(jl,jk+1,kd) / ( pwu2h(jl,jk+1,kd) - pwu2h(jl,jk,kd) ) 
                     
        endif

        if ( pwu2h(jl,jk,kd) < pw2thresh ) then   !allow parcel to overcome layers 
          lddone(jl,kd)  = .true.                 !with small negative kin. energy
        endif                                  !but remember last w(z)=0 (z=pzptop)


      endif
    enddo
    if (is == 0) exit
    
        
  enddo !jk

  
  !protect for updrafts that have reached top level (let's hope this is unnecessary!)
  do jl=kidia,kfdia
    if ( .not. lddone(jl,kd) ) then
      lddone(jl,kd) = .true.
      kptop(jl,kd)  = jkmin+1
      kplzb(jl,kd)  = jkmin+1
      pzptop(jl,kd) = pgeoh(jl,jkmin+1) * zrg
    endif
    
    !rn testing rico: 
    !kplcl(jl,kd) = 81
    
  enddo


! if (lhook) call dr_hook('vdfparcel',1,zhook_handle)
end subroutine vdfparcel
