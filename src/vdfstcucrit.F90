!options xopt(hsfun)
subroutine vdfstcucrit (kidia   , kfdia   , klon    , klev    , kdraft  , &
                   &    ptm1    , pslgm1  , pqtm1   , papm1   , &
                   &    pstabthresh, pclddepth, pbirthresh, pdzcloud, &
                   &    kptop   , kpbltype , ldnodecp, &
                   &    pstability  )  
!     ------------------------------------------------------------------

!**   *vdfstcucrit* - criteria for stratocumulus occurrence
!
!     based on original vdfhghtn.f90 (cy29r1 and earlier) by
!             a.p. siebesma    30/06/99  
!             m. ko"hler       3/12/2004 
!     put into separate file by   
!             roel neggers     12/04/2005 


!     purpose
!     -------

!     determine stratocumulus occurrence

!     interface
!     ---------

!     *vdfstcucrit* is called by *vdfhghtn*

!     parameter     description                                   units
!     ---------     -----------                                   -----
!     input parameters (integer):

!     *kidia*        start point
!     *kfdia*        end point
!     *klev*         number of levels
!     *klon*         number of grid points per packet
!     *kptop*        highest half level below pbl height, and
!                    pbl top full level (pzptop is within that layer)
!     *kdraft*       number of explicitly modeled drafts - currently 3:
!                    1: test parcel
!                    2: rising dry thermals which stop at cloud base or inversion
!                    3: rising dry thermals which become cloudy
!                    (4: downdrafts .. to be done)


!     input parameters (real):

!     *ptm1*         temperature at t-1                               k
!     *pslgm1*       liquid static energy (slg) at t-1                k
!     *pqtm1*        total specific humidity at t-1                   kg/kg
!     *papm1*        pressure at full level at t-1                    pa
!     *pstabthresh*  stability criterion (klein & hartmann criteria)  k                      
!     *pclddepth*    threshold cloud thickness for stcu/cu transition m
!     *pbirththrash* threshold bir (tke decoupling criteria)          1
!     *pdzcloud*     cloud thickness                                  m


!     input parameters (logical):

!     *ldnodecp*     true:  never decouple
!                    false: maybe decouple



!     output parameters (integer):

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

!#include "tsmbkind.h"
use thrm, only : rslf, esl
use garbage, only : foealfa
use parkind1  ,only : jpim     , jprb

! use yomhook   ,only : lhook    , dr_hook

use yos_cst   , only : rd       , rg      , rcpd     , retv     , rlvtt, &
                     &rlstt    ,rtt

use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les    ,r4ies   ,r5les    ,r5ies, &
                     &r5alvcp  ,r5alscp  ,ralvdcp  ,ralsdcp  ,rtwat   ,rtice    ,rticecu, &
                     &rtwat_rtice_r      ,rtwat_rticecu_r

implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim),intent(in)    :: kdraft 
integer(kind=jpim),intent(in)    :: kptop(klon,kdraft) 
integer(kind=jpim),intent(inout) :: kpbltype(klon) 
real(kind=jprb)   ,intent(in)    :: ptm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pslgm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqtm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: papm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pstabthresh
real(kind=jprb)   ,intent(in)    :: pclddepth
real(kind=jprb)   ,intent(in)    :: pbirthresh
real(kind=jprb)   ,intent(in)    :: pdzcloud(klon) 
logical           ,intent(in)    :: ldnodecp(klon) 
real(kind=jprb)   ,intent(out)   :: pstability(klon)


!*         0.2    local variables

real(kind=jprb) ::   zrg,  zdctei(klon)

!          variables for cloud base estimation

real(kind=jprb) ::    zqs(klon,0:klev)  , zalfaw  , zfacw   , zfaci   , zfac        , &
                    & zesdp   , zcor    , zdqsdtemp(klon)   , zbeta   , zdsl(klev)  , &
        & zdqt(klev)
		      
integer(kind=jpim) :: is, jk, jl, jd

integer(kind=jpim) :: i700(klon)

! real(kind=jprb) ::    zhook_handle


!#include "cuadjtq.intfb.h"

!dir$ vfunction exphf
! #include "fcttre.h"



!     -----------------------------------------------------------------

!*         1.     set some constants
!                 --------------------

! if (lhook) call dr_hook('vdfstcucrit',0,zhook_handle)

! optimization
zrg         = 1.0_jprb/rg




!     -----------------------------------------------------------------

!*         2.     prepare variables appearing in some criteria 
!                 --------------------------------------------

!*         2.1  stability criteria == theta(700hpa) - theta(sfc)

!          find index i700 of pressure closest to 700hpa

  do jl=kidia,kfdia
    i700(jl) = 0
  enddo
  do jk=1,klev
    do jl=kidia,kfdia
      if ( i700(jl) == 0  .and.  papm1(jl,jk) > 70000.0_jprb ) then
        i700(jl) = jk
      endif
    enddo
  enddo
  do jl=kidia,kfdia
    if ( i700(jl) > 1 ) then
      if ( abs(papm1(jl,i700(jl)-1)-70000.0_jprb)  <  &
         & abs(papm1(jl,i700(jl))  -70000.0_jprb)    ) then  
        i700(jl) = i700(jl)-1
      endif

      pstability(jl) = ptm1(jl,i700(jl)) * ( 1.0e5_jprb/papm1(jl,i700(jl)) ) ** (rd/rcpd) &
                   & - ptm1(jl,klev)     * ( 1.0e5_jprb/papm1(jl,klev) )     ** (rd/rcpd)  
    else
      pstability(jl) = 0.0_jprb
    endif

  enddo


!*         2.2  cloud top entrainment instability (ctei) criteria

  do jl=kidia,kfdia

   if ( .false. ) then
    if ( kpbltype(jl) == 2) then

      jk            = kptop(jl,1)  ! pbl top full level taken as k+1
                                   ! full level above pbl taken as k-1

!          qsat (full level)
!       zqs(jl,jk+1)  = foeewm(ptm1(jl,jk+1))/papm1(jl,jk+1)
!       zqs(jl,jk+1)  = min(0.5_jprb,zqs(jl,jk+1))
!       zqs(jl,jk+1)  = zqs(jl,jk+1)/(1.0_jprb-retv*zqs(jl,jk+1))
      zqs(jl,jk+1) = rslf(ptm1(jl,jk+1),papm1(jl,jk+1))
!          calculate dqs/dt correction factor (full level)
      zalfaw        = foealfa(ptm1(jl,jk+1))
      zfacw         = r5les/((ptm1(jl,jk+1)-r4les)**2)
      zfaci         = r5ies/((ptm1(jl,jk+1)-r4ies)**2)
      zfac          = zalfaw*zfacw+(1.0_jprb-zalfaw)*zfaci
      zesdp         = esl(ptm1(jl,jk+1))/papm1(jl,jk+1)
      zcor          = 1.0_jprb/(1.0_jprb-retv*zesdp)
      zdqsdtemp(jl) = zfac*zcor*zqs(jl,jk+1)

!          ctei
      zbeta = ( 1.0_jprb + (1.0_jprb+retv) * ptm1(jl,jk+1) * zdqsdtemp(jl) ) &
          & / ( 1.0_jprb + rlvtt/rcpd * zdqsdtemp(jl) )  
      zdsl(jl)  = pslgm1(jl,jk-1) - pslgm1(jl,jk+1)
      zdqt(jl)  = pqtm1 (jl,jk-1) - pqtm1 (jl,jk+1)

      zdctei(jl) = zbeta * zdsl(jl) &
               & + ( zbeta - rcpd/rlvtt * ptm1(jl,jk+1) ) * rlvtt * zdqt(jl)  

    else

      zdctei(jl) = 1.0

    endif
   endif




!     -----------------------------------------------------------------

!*         3      stratocumulus - shallow cumulus criteria: 
!*                * cloud thickness = 1000m
!*                * stability = 15k
!*                * tke decoupling, bir = 0.1
!                 -----------------------------------------


    if ( .not. ldnodecp(jl) .and.  kpbltype(jl) == 2 ) then

!..........stability criteria (klein & hartmann 1993)
      if ( pstability(jl) < pstabthresh ) then 

!..........cloud thickness criteria
!     if ( pdzcloud(jl) > pclddepth ) then

!..........ctei...
!     if ( zdctei(jl) < 0 ) then

!..........tke decoupling (or cloud thickness criteria for safety)
!     if ( pbir(jl) > pbirthresh .or. pdzcloud(jl) > pclddepth ) then

!..........always decouple
!     if ( .true. ) then

        kpbltype(jl) = 3   !decouple: pbl type 3 (shallow cumulus)

      endif
    endif


  enddo !jl




! if (lhook) call dr_hook('vdfstcucrit',1,zhook_handle)
end subroutine vdfstcucrit
