subroutine vdfexcu(kidia  , kfdia  , klon   , klev   , kdraft , ptmst  , pz0mm  , &
                  &phrlw  , phrsw  , &
                  &pum1   , pvm1   , ptm1   , pqm1   , plm1   , pim1   , &
                  &paphm1 , papm1  , pgeom1 , pgeoh  , pcptgz , &
 ! diagnostic output
                  &pextr2 , kfldx2 , pextra , klevx  , kfldx  , &
 !
                  &pkmfl  , pkhfl  , pkqfl  , &
                  &pcfm   , pcfh   , ptauxcg, ptauycg, &
                  &pricui , pmcu   , pdthv  , pmflx  , kvartop, &
                  &pzinv  , khpbl  , pkh    , pzcldbase , pzcldtop     , kpbltype)
!     ------------------------------------------------------------------

!**   *vdfexcu* - determines the exchange coefficients between the
!                 upper model levels with stability as a function of
!                 obukhov-l

!     a.c.m. beljaars  26/03/90.  original
!     a.c.m. beljaars  26/03/99   tiling of the land surface.
!     j.hague          13/01/2003 mass vector functions
!     m. ko"hler        3/12/2004 moist advection-diffusion incl.
!                                 k,cloud and cloud top entrainment
!     p. lopez         02/06/2005 removed option for linearized
!                                 physics (now called separately)
!     r. neggers       01/06/2006 reorganization (into internal & interface k-modes)
!                                 entrainment efficiency closure at cumulus pbl top
!     a. beljaars      29/03/2006 counter gradient stresses (brown and grant, 1997)
!
!     purpose
!     -------

!     determine exchange coefficients between the upper model levels

!     interface
!     ---------

!     *vdfexcu* is called by *vdfmain*

!     input parameters (integer):

!     *kidia*        start point
!     *kfdia*        end point
!     *klev*         number of levels
!     *klon*         number of grid points per packet

!     input parameters (real):

!     *ptmst*        double time step (single at 1th step)
!     *pum1*         x-velocity component at t-1
!     *pvm1*         y-velocity component at t-1
!     *ptm1*         temperature at t-1
!     *pqm1*         specific humudity at t-1
!     *paphm1*       pressure at half levels at t-1
!     *papm1*        pressure at full levels at t-1
!     *pgeom1*       geopotential at t-1
!     *pcptgz*       dry static energy
!     *pkmfl*        kinematic momentum flux                [#]
!     *pkhfl*        kinematic heat flux                    [#]
!     *pkqfl*        kinematic moisture flux                [#]
!     *pzinv*        inversion height                   [m]
!     *pkh*          turb. diff. coeff. for heat above surf. lay.  (m2/s)

!     output parameters (real):

!     *pcfm*         prop. to exch. coeff. for momentum (c-star in doc.)
!     *pcfh*         prop. to exch. coeff. for heat     (c-star in doc.)
!                    (only pcfm(*,1:klev-1) and
!                          pcfh(*,1:klev-1) are computed)
!     *ptauxcg*      counter gradient stress x-component    (n/m2)
!     *ptauycg*      counter gradient stress y-component    (n/m2)

!     remark: [#] unused parameters in tangent linear and adjoint versions
!     ------

!     method
!     ------

!     see documentation

!     ------------------------------------------------------------------
use parkind1  ,only : jpim     ,jprb
! ! use yomhook   ,only : lhook,   dr_hook

use yos_cst   , only : rg       ,rd       ,rcpd     ,retv     ,ratm
use yoethf   , only : rvtmp2
use yoevdf   , only : rlam     ,rkap     ,rvdifts  ,repdu2   ,lldiag   
use yoevdfs  , only : jpritbl  ,ritbl    ,aritbl   ,rchba    ,&
                    & rchbb    ,rchbd    ,rchb23a  ,rchbbcd  ,rchbcd   ,&
                    & rcheta   ,rchetb   ,rcdhalf  ,rcdhpi2  ,rimax    ,&
                    & dritbl   ,dri26  ,phims, phihs, phimu, phihu
use yoephli  , only : rlpmixl  ,rlpbeta
use yomjfh   , only : n_vmass
use yos_exc, only : repust

implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim),intent(in)    :: kdraft
real(kind=jprb)   ,intent(in)    :: ptmst 
real(kind=jprb)   ,intent(in)    :: pz0mm(klon) 
real(kind=jprb)   ,intent(in)    :: phrlw(klon,klev) 
real(kind=jprb)   ,intent(in)    :: phrsw(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pum1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pvm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: ptm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: plm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pim1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: paphm1(klon,0:klev) 
real(kind=jprb)                  :: papm1(klon,klev) ! argument not used
real(kind=jprb)   ,intent(in)    :: pgeom1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pgeoh(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: pcptgz(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pkmfl(klon) 
real(kind=jprb)   ,intent(in)    :: pkhfl(klon) 
real(kind=jprb)   ,intent(in)    :: pkqfl(klon) 
real(kind=jprb)   ,intent(inout) :: pcfm(klon,klev) 
real(kind=jprb)   ,intent(inout) :: pcfh(klon,klev) 
real(kind=jprb)   ,intent(inout) :: ptauxcg(klon,klev) 
real(kind=jprb)   ,intent(inout) :: ptauycg(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pzinv(klon) 
integer(kind=jpim),intent(in)    :: khpbl(klon)
real(kind=jprb)   ,intent(out)   :: pkh(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pzcldbase(klon)
real(kind=jprb)   ,intent(in)    :: pzcldtop(klon)
integer(kind=jpim),intent(in)    :: kpbltype(klon) 
real(kind=jprb)   ,intent(in)    :: pricui(klon) 
real(kind=jprb)   ,intent(in)    :: pdthv(klon) 
real(kind=jprb)   ,intent(in)    :: pmcu(klon) 
real(kind=jprb)   ,intent(inout) :: pmflx(klon,0:klev,kdraft) 
integer(kind=jpim),intent(in)    :: kvartop(klon)
!          diagnostic output
integer(kind=jpim),intent(in)     :: kfldx2, klevx, kfldx
real(kind=jprb)   ,intent(inout)  :: pextr2(klon,kfldx2), pextra(klon,klevx,kfldx)

!*         0.2    local variables

real(kind=jprb) ::    zri(klon),zmgeom(klon),zust(klon),&
                    & zdtv(klon),zkhvfl(klon),zl(klon),zphim(klon),&
                    & zphih(klon),ztauxcg(klon),ztauycg(klon)    
real(kind=jprb) ::    zdu2(klon+n_vmass)

integer(kind=jpim) :: irib, jk, jl, jlen
real(kind=jprb) ::    zentrsfc, zentrrad, zentrtop, &
                    & z2geomf, za, zalh2, zalm2, zb, zcb, &
                    & zcd, zcfnc1, zrho, zuabs, &
                    & zcons13, zcons1, zhu1, zhu2, zwst3, zrg, &
                    & zdh, zdl, zdroro, zeps, zeta, zhlm2, &
                    & zlim, zlim2, zphikh, zphikm, zscf, &
                    & zx2, zz, zwtventr, zkh, zcfhnew, &
                    & zml, zbase, zvsc, zkcld, &
                    & zkfacedmf, ztaux, ztauy

real(kind=jprb) ::    zzh, zifltgm, zifltgh, zifmom, zifmoh, zbm, zbh, zcm, zch, zdudz
                    
logical ::         llricu
                    
real(kind=jprb) :: zdradflx(klon), zradkbase(klon), zradkdepth(klon), &
                 & zradkfac(klon)

real(kind=jprb) :: zwecutop(klon)  , zdthvcutop, &
                 & zthven(klon,klev), zthen(klon,klev), zfac

real(kind=jprb) ::    ztmp1(kfdia-kidia+1+n_vmass)
real(kind=jprb) ::    ztmp2(kfdia-kidia+1+n_vmass)
real(kind=jprb) ::    ztmp3(kfdia-kidia+1+n_vmass)
real(kind=jprb) ::    ztmp4(kfdia-kidia+1+n_vmass)
real(kind=jprb) ::    ztmp5(kfdia-kidia+1+n_vmass)
real(kind=jprb) ::    zhook_handle
! 
! interface
! #include "surf_inq.h"
! end interface
! 
! #include "fcvdfs.h"



!     ------------------------------------------------------------------

!*         1.     initialize constants
!                 --------------------

! if (lhook) call dr_hook('vdfexcu',0,zhook_handle)

zentrsfc  = 0.2_jprb       ! factor for surface based top entrainment 
zentrrad  = 0.2_jprb       ! factor for radiative based top entrainment 
zentrtop  = 0.4_jprb       ! entrainment efficiency factor at cumulus pbl top, as proposed by wyant et al (jas, 1997)

zcd       = 1.0_jprb
zcb       = 5.0_jprb
zeps      = 1.e-10_jprb

llricu = .true.   ! switch for top-entrainment efficiency closure using ri^cu at cumulus pbl top
!llricu = .false.  

!zkfacedmf = 1.0_jprb
zkfacedmf = 0.8_jprb     !aup = 5%   !cy32r3
!zkfacedmf = 0.692_jprb   !aup = 10%

! optimization
zrg       = 1.0_jprb/rg
zcons13   = 1.0_jprb/3._jprb
zcons1    = 0.5_jprb*rkap*zrg/rlam
zhlm2     = 1.0_jprb / ((2.0_jprb*rg*rlpmixl)**2)

if(n_vmass > 0) then
  jlen=kfdia-kidia+n_vmass-mod(kfdia-kidia,n_vmass)
endif


!     ------------------------------------------------------------------

!*         2.     prepare scaling coefficients
!                 ----------------------------

  do jl=kidia,kfdia
    zust  (jl)=sqrt(max(pkmfl(jl),repust**2))
    zkhvfl(jl)=pkhfl(jl)+retv*ptm1(jl,klev)*pkqfl(jl)

    ptauxcg(jl,klev)=0.0_jprb
    ptauycg(jl,klev)=0.0_jprb

    if (zkhvfl(jl)  <  0.0_jprb) then
      zwst3=-zkhvfl(jl)*pzinv(jl)*rg/ptm1(jl,klev)
    else
      zwst3=0.0_jprb
    endif
    zrho =paphm1(jl,klev)/( rd*ptm1(jl,klev)*(1.0_jprb+retv*pqm1(jl,klev)) )
    zuabs=max(0.1_jprb,sqrt(pum1(jl,klev)**2+pvm1(jl,klev)**2))
    zhu1 =zrho*pkmfl(jl)/zuabs
    ztaux=zhu1*pum1(jl,klev)
    ztauy=zhu1*pvm1(jl,klev)

    zhu2=2.7_jprb*zwst3/(zust(jl)**3+0.6_jprb*zwst3)
    ztauxcg(jl)=ztaux*zhu2
    ztauycg(jl)=ztauy*zhu2
  enddo
  
  !calculate full level mean theta_v profile for later use
  do jk=klev,1,-1
    do jl=kidia,kfdia
      zthen(jl,jk)   = ( papm1(jl,jk)/ratm )**(-rd/rcpd) * ptm1(jl,jk)
      zthven(jl,jk)  = zthen(jl,jk) * &
                     & ( 1.0_jprb + retv * pqm1(jl,jk)   - plm1(jl,jk)   - pim1(jl,jk)   )
    enddo
  enddo

!  if (lldiag) then
!    do jk=1,klev
!    do jl=kidia,kfdia
!      if (jk>=kvartop(jl)-1) then
!        pextra(jl,jk,4) = zthven(jl,jk) - zthven(jl,klev)
!        pextra(jl,jk,5) = zthen(jl,jk)  - zthen(jl,klev)
!      endif  
!    enddo
!    enddo
!  endif  
  

!          calculate pbl cloud top radiative flux jump [km/s] (cooling)
!          for top-driven k and entrainment velocity formulations.

  do jl=kidia,kfdia
    select case (kpbltype(jl))
      case(2)
        zradkdepth(jl) = pzinv(jl)                               
        zradkbase(jl)  = 0._jprb    
        zradkfac(jl)   = 1._jprb                        
      case(3)
        zradkdepth(jl) = max( pzcldtop(jl) - pzcldbase(jl),  0._jprb )
        zradkbase(jl)  = pzcldtop(jl) - zradkdepth(jl)
        zradkfac(jl)   = 0.25_jprb                        
      case default
        zradkdepth(jl) = -100._jprb
        zradkbase(jl)  = -100._jprb
        zradkfac(jl)   = 0._jprb                        
    end select 
  enddo


  zdradflx(:) = 0.0_jprb
  do jk=klev-1,1,-1
    do jl=kidia,kfdia

!      if ( kpbltype(jl) == 2 ) then  !rn for now, only stratocumulus: extension to cumulus planned (for intermediate scenarios like atex)

        if ( pgeoh(jl,jk)*zrg <= zradkbase(jl)+zradkdepth(jl) .and. zradkbase(jl)+zradkdepth(jl) < pgeoh(jl,jk-1)*zrg ) then
!          zdradflx(jl) = -  phrlw(jl,jk+1)                 * (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg   !use lw divergence only
          zdradflx(jl) = - (phrlw(jl,jk+1)+phrsw(jl,jk+1)) * (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
!         ... add solar heating at 2nd cld level
          if ( pzcldbase(jl) < pgeoh(jl,jk+1)*zrg .and. jk < klev-1 ) then 
            zdradflx(jl) = zdradflx(jl) - phrsw(jl,jk+2)   * (pgeoh(jl,jk+1)-pgeoh(jl,jk+2))*zrg
          endif
          zdradflx(jl) = max( zdradflx(jl), 0.0_jprb )    !safety against rad. heating cases
        endif

!      endif

    enddo
  enddo



!     ------------------------------------------------------------------

!*         3.     vertical loop - non-linear physics
!                 ----------------------------------

  if(n_vmass > 0) then
    if(kfdia-kidia+1 /= jlen) then
      ztmp1(kfdia-kidia+2:jlen)=1.0_jprb 
      ztmp2(kfdia-kidia+2:jlen)=1.0_jprb 
      ztmp3(kfdia-kidia+2:jlen)=1.0_jprb 
      zdu2(kfdia+1:kidia+jlen-1)=1.0_jprb 
    endif
  endif


!***
  do jk=klev-1,1,-1
!***

    do jl=kidia,kfdia
      pcfm(jl,jk)=0.0_jprb
      pcfh(jl,jk)=0.0_jprb
      pkh(jl,jk) =0.0_jprb
      ptauxcg(jl,jk)=0.0_jprb
      ptauycg(jl,jk)=0.0_jprb
    enddo

    if(n_vmass <= 0) then   ! efficiency of exponentials

!          compute ri-number

      do jl=kidia,kfdia
        zdu2(jl)=max(repdu2,(pum1(jl,jk)-pum1(jl,jk+1))**2&
                         & +(pvm1(jl,jk)-pvm1(jl,jk+1))**2)  
        zdroro= 2.0_jprb * (pcptgz(jl,jk)-pcptgz(jl,jk+1))&
         & / ( pcptgz(jl,jk)+pcptgz(jl,jk+1)&
         &   - pgeom1(jl,jk)-pgeom1(jl,jk+1))&
         & - (rvtmp2-retv)*(pqm1(jl,jk)-pqm1(jl,jk+1))
        zdtv(jl)=( (pcptgz(jl,jk)-pcptgz(jl,jk+1))&
         & - (rvtmp2-retv)*0.5_jprb * (pqm1(jl,jk)-pqm1(jl,jk+1))&
         & * (pcptgz(jl,jk)+pcptgz(jl,jk+1)) ) * (1.0_jprb/rcpd)  
        zmgeom(jl)=pgeom1(jl,jk)-pgeom1(jl,jk+1)
        zri(jl)=zmgeom(jl)*zdroro/zdu2(jl)
      enddo

    else

      do jl=kidia,kfdia
        zdu2(jl)=max(repdu2,(pum1(jl,jk)-pum1(jl,jk+1))**2&
                         & +(pvm1(jl,jk)-pvm1(jl,jk+1))**2)
        ztmp2(jl-kidia+1)= pcptgz(jl,jk)+pcptgz(jl,jk+1)&
                        & -pgeom1(jl,jk)-pgeom1(jl,jk+1)
        zdtv(jl)=( (pcptgz(jl,jk)-pcptgz(jl,jk+1))&
         & -(rvtmp2-retv)*0.5_jprb* (pqm1(jl,jk)-pqm1(jl,jk+1))&
         & *(pcptgz(jl,jk)+pcptgz(jl,jk+1)) ) * (1.0_jprb/rcpd)
        zmgeom(jl)=pgeom1(jl,jk)-pgeom1(jl,jk+1)
      enddo
! 
!       call vrec(ztmp4,zdu2(kidia),jlen)
!       call vrec(ztmp5,ztmp2,jlen)

      do jl=kidia,kfdia
        zdroro= 2.0_jprb * (pcptgz(jl,jk)-pcptgz(jl,jk+1))&
         & *ztmp5(jl-kidia+1)&
         & - (rvtmp2-retv)*(pqm1(jl,jk)-pqm1(jl,jk+1))  
        zri(jl)=zmgeom(jl)*zdroro*ztmp4(jl-kidia+1)
      enddo

    endif

    do jl=kidia,kfdia

!          compute stability functions

      if (zri(jl)  >  0.0_jprb) then

!        interplolate eta with splines for positive richardson
!        numbers

        irib=int(zri(jl)*(1.0_jprb/dritbl))+1
        
        if (irib  >=  jpritbl) then
!           linear extension of look-up table
          zeta = ritbl(jpritbl)*(zri(jl)*(1.0_jprb/rimax))
        else
          zx2  = irib*dritbl
          za   = (zx2-zri(jl))*(1.0_jprb/dritbl)
          zb   = 1.0_jprb-za
          zeta = za*ritbl(irib) + zb*ritbl(irib+1)&
           & +( (za**3-za)*aritbl(irib)&
           & +(  zb**3-zb)*aritbl(irib+1) )*dri26  
        endif

!        stable phi-functions

        zphim(jl) = phims(zeta)
        zphih(jl) = phihs(zeta)
      else

!        unstable situations

        zeta  = zri(jl)
        zphim(jl) = phimu(zeta)
        zphih(jl) = phihu(zeta)
      endif
    enddo
!   if(n_vmass <= 0) then ! vector mass taken out because completely changed code

    do jl=kidia,kfdia

!-------- up to cy32r3 --------
!
!!          common factors for stable and unstable
!
!     z2geomf=pgeom1(jl,jk)+pgeom1(jl,jk+1)+2.0_jprb*rg*pz0mm(jl)
!     zlim=rlpbeta + (1.0_jprb-rlpbeta)/(1.0_jprb+z2geomf*z2geomf*zhlm2)
!     zlim2=zlim*zlim
!     zalm2=zlim2*(0.5_jprb*rkap*zrg*z2geomf/(1.0_jprb+zcons1*z2geomf))**2
!     zalh2=zalm2
!     zcfnc1=rvdifts*ptmst*rg**2 * paphm1(jl,jk)&
!      & /( 0.5_jprb*rd * zmgeom(jl)&
!      & *( ptm1(jl,jk  )*(1.0_jprb+retv*pqm1(jl,jk  ))&
!      & +  ptm1(jl,jk+1)*(1.0_jprb+retv*pqm1(jl,jk+1))))  
!     zcfnc=rg*zcfnc1*sqrt(zdu2(jl))/zmgeom(jl)
!
!!          dimensionless coefficients multiplied by pressure
!!          thicknesses for momentum and heat exchange.
!
!     if (zri(jl)  >  0.0_jprb) then  ! statically stable
!       zscf=sqrt(1.0_jprb+zcd*zri(jl))
!       pcfm(jl,jk)=zcfnc*zalm2/(1.0_jprb+2.0_jprb*zcb*zri(jl)/zscf)
!       pcfh(jl,jk)=zcfnc*zalh2/(1.0_jprb+2.0_jprb*zcb*zri(jl)*zscf)
!     else                            ! statically unstable
!       pcfm(jl,jk)=zcfnc*zalm2/(zphim(jl)**2)
!       pcfh(jl,jk)=zcfnc*zalm2/(zphim(jl)*zphih(jl))
!     endif
!-----------------------------

!-------- cy32r3 -------------
! new k: k,ltg scales with l=kappa*z 
!        k,mo  scales with l=150m above surface layer

!          common factors for stable and unstable

      zifmom  = 1.0_jprb / (zphim(jl)**2)                              !f(mo),m
      zifmoh  = 1.0_jprb / (zphim(jl)*zphih(jl))                       !f(mo),h
      zdudz   = sqrt(zdu2(jl))/zmgeom(jl)*rg                           !shear
      zcfnc1=rvdifts*ptmst*rg**2 * paphm1(jl,jk)&                      !factor for vdfdifh/m
       & /( 0.5_jprb*rd * zmgeom(jl)&
       & *( ptm1(jl,jk  )*(1.0_jprb+retv*pqm1(jl,jk  ))&
       & +  ptm1(jl,jk+1)*(1.0_jprb+retv*pqm1(jl,jk+1))))  

!          dimensionless coefficients multiplied by pressure
!          thicknesses for momentum and heat exchange.

      if ( zri(jl) > 0.0_jprb ) then  ! statically stable
        zzh     = 0.5_jprb * zrg * (pgeom1(jl,jk)+pgeom1(jl,jk+1)) + pz0mm(jl)
        zscf    = sqrt(1.0_jprb+zcd*zri(jl))
        zifltgm = 1.0_jprb / (1.0_jprb + 2.0_jprb *zcb * zri(jl)/zscf) !f(ltg),m
        zifltgh = 1.0_jprb / (1.0_jprb + 2.0_jprb *zcb * zri(jl)*zscf) !f(ltg),h
        zbm     = rkap * zzh * sqrt(zifltgm) 
        zbh     = rkap * zzh * sqrt(zifltgh) 
        zcm     = 150.0_jprb * sqrt(zifmom )
        zch     = 150.0_jprb * sqrt(zifmoh )

        pcfm(jl,jk) = zcfnc1 * zdudz * (zbm*zcm/(zbm+zcm))**2
        pcfh(jl,jk) = zcfnc1 * zdudz * (zbh*zch/(zbh+zch))**2
      else                            ! statically unstable
        pcfm(jl,jk) = zcfnc1 * zdudz * 150.0_jprb**2 * zifmom
        pcfh(jl,jk) = zcfnc1 * zdudz * 150.0_jprb**2 * zifmoh
       
      endif
!-----------------------------


      !  overwrite ri-based k values within boundary layer
      !if ( jk >= kvartop(jl) .and. kvartop(jl)>0 ) then 
      !  pcfh(jl,jk) = 0._jprb
      !  pcfm(jl,jk) = 0._jprb
      !endif

      
      
      !------------------------------------------------------------------
      !
      !   3.1   internal k mode of mixed layer
      !         ------------------------------
      !
      !         using a prescribed vertical structure.
      !
      
      zz     = pgeoh(jl,jk)*zrg

      !if (pgeoh(jl,jk-1)*zrg <= pzinv(jl) .and. zkhvfl(jl)<0._jprb ) then  ! up to level below entr. level
      if ( jk > khpbl(jl) .and. zkhvfl(jl)<0._jprb ) then 

        zl(jl) = zust (jl)**3*ptm1(jl,klev)/(rkap*rg*(zkhvfl(jl)-zeps))

        zdh    = zz/pzinv(jl)
        zdl    = zz/zl(jl)
        zeta   = zz/zl(jl)
        zphikh = (1-39._jprb*zdl)**(-zcons13)
        zphikm = (1-15._jprb*zdl)**(-zcons13)

        !   k,surface
        
        pcfh(jl,jk)  = zkfacedmf * zcfnc1 * rkap / zphikh * zust(jl) * zz * (1.0_jprb-zdh)**2
        pcfm(jl,jk)  = zkfacedmf * zcfnc1 * rkap / zphikm * zust(jl) * zz * (1.0_jprb-zdh)**2

        ptauxcg(jl,jk)=zdh*(1._jprb-zdh)**2*ztauxcg(jl)
        ptauycg(jl,jk)=zdh*(1._jprb-zdh)**2*ztauycg(jl)

      endif
      
      
      
      !------------------------------------------------------------------
      !
      !   3.1   internal k mode in cumulus cloud layer
      !         ------------------------------
      !
      
      if ( kpbltype(jl) == 3 .or. kpbltype(jl) == 4) then   

        if ( jk < khpbl(jl) .and. jk >= kvartop(jl) ) then
          if ( jk > kvartop(jl) ) then 
            !  k diffusion within cumulus layer
            !pcfh(jl,jk)  = pcfh(jl,jk)          !testing: ri diffusion 
            !pcfm(jl,jk)  = pcfm(jl,jk)
            pcfh(jl,jk)  = pcfh(jl,khpbl(jl))   !cy32r3
            pcfm(jl,jk)  = pcfm(jl,khpbl(jl))
          else
            !  reset k 
            pcfh(jl,jk)  = 0._jprb 
            pcfm(jl,jk)  = 0._jprb
          endif
        endif
        
        !  protect top-entrainment of shallow cu topped mixed layers against a zero
        !    buoyancy jump (dthv) through the cloud base transition layer.
        !    this can easily occur just before h grows one layer,due to dq
        !    cancelling ds in dthv.    -rn
        !
        zdtv(jl) = max( zdtv(jl), 0.2_jprb ) 
        
      endif



      !------------------------------------------------------------------
      !
      !   3.2    internal k mode driven by cloud-top cooling
      !          -------------------------------------------
      !
      !          as in lock et al. (2000, mwr p3187f), equ. 5:
      !          using simplified radiative velocity scale as in lock, 1998, equ. 12
      !                 and ignore buoyancy reversal velocity scale
      !

      if ( zkhvfl(jl)<0._jprb .and. zradkdepth(jl)>0._jprb) then
        if ( zz >= zradkbase(jl)  .and.  zz <= zradkbase(jl)+zradkdepth(jl) ) then  
          zvsc  = ( rg / ptm1(jl,jk) * zradkdepth(jl) * zdradflx(jl) ) ** zcons13 
          !zkcld = 0.85_jprb * rkap * zvsc &
          zkcld = zradkfac(jl) * 0.85_jprb * rkap * zvsc &
              & * (zz-zradkbase(jl)) ** 2 / zradkdepth(jl) &
              & * ( 1 - (zz-zradkbase(jl)) / zradkdepth(jl) ) ** 0.5_jprb 
          if (kpbltype(jl)==2) then     
            pcfh(jl,jk)  = pcfh(jl,jk) + zcfnc1 * zkcld
            pcfm(jl,jk)  = pcfm(jl,jk) + zcfnc1 * zkcld * 0.75_jprb
          endif  
        endif
      endif



      !------------------------------------------------------------------
      !
      !   3.3   interface k at cumulus pbl top
      !         ------------------------------
      !
      !         at top level, mass flux is replaced by diffusion, using an entrainment efficiency formulation.
      !         this is important for representing the intermediate regime (stcu->cu transitions)    -rn
      !

      if ( llricu .and. jk == kvartop(jl) .and. kpbltype(jl) == 3) then 

        !entrainment efficiency - after wyant et al. (jas, 1997)
        !zwecutop(jl) = pmcu(jl) * zentrtop * pricui(jl)
        zwecutop(jl) = 2._jprb * pmcu(jl) * zentrtop * pricui(jl)
        zwecutop(jl) = max(0.0_jprb,zwecutop(jl))
          
        !  translation into k [m2/s] at this level: k = entrainment velocity * mixing-length (dz)
        zkh = zwecutop(jl)                                !top-entrainment by overshooting surface-driven thermals
        zkh = zkh + zentrrad * zdradflx(jl) / zdtv(jl )   !add cloud top cooling driven entrainment
        zkh = zkh * zmgeom(jl) * zrg
        
        zkh     = max(0.0_jprb,zkh)
        zcfhnew = zcfnc1 * zkh
            
        pcfh(jl,jk) = zcfhnew
        pcfm(jl,jk) = zcfhnew * 0.75_jprb
            
        !  reset any updraft m
        pmflx(jl,jk,2) = 0._jprb
        pmflx(jl,jk,3) = 0._jprb

      endif



      !---------------------------------------------------------------------
      !
      !   3.4   interface k at mixed layer top 
      !         ------------------------------
      ! 
      !    this is cloud top for dry & stratocu pbl, and cloud base in 
      !    shallow cu pbl (pzinv).
      !
      !    w_e is the mixed layer top-entrainment rate, defined as
      !
      !          w_e = w'thv'_h / dthv_h = -0.2 w'thv'_s / dthv_h = a/ri w_*.
      !
      !    this w_e is here translated into k at this level.
      !

      !if ( pgeoh(jl,jk)*zrg <= pzinv(jl)  .and.  pzinv(jl) < pgeoh(jl,jk-1)*zrg ) then
      if ( jk == khpbl(jl) ) then 
        
        !wthv_h = -0.2 wthv_s
        zwtventr = -zentrsfc * zkhvfl(jl) * 0.1_jprb
        if (kpbltype(jl) == 2) then
          zwtventr = -zentrsfc * zkhvfl(jl)
        endif  

        !--- special stratocumulus treatment: radiation impacts on entrainment
        !---
        !--- entrainment velocity * t,v jump due to lw-radiative cooling & sfc flux
        !---    (lock & macvean, 1999, equ. 11)
        if ( kpbltype(jl) == 2 ) then                     
          zwtventr = zwtventr + zentrrad * zdradflx(jl) !radiation flux jump
        endif

        
        zwtventr = max(0.0_jprb,zwtventr)

        !zdtv(jl) = 1.0_jprb
        !zdtv(jl) = 10.*zdtv(jl)
        !zdtv(jl) = 2._jprb * pdthv(jl)

        zkh     = zwtventr * zmgeom(jl) / ( rg * zdtv(jl) )

        zkh     = max(0.0_jprb,zkh)
        zcfhnew = zcfnc1 * zkh

        !rn pcfh(jl,jk)  = max(pcfh(jl,jk),zcfhnew)           !protection against k=0
        !rn pcfm(jl,jk)  = max(pcfm(jl,jk),zcfhnew * 0.75_jprb)
        pcfh(jl,jk)  = zcfhnew           
        pcfm(jl,jk)  = zcfhnew * 0.75_jprb

!        if (lldiag) then
!          pextr2(jl,4) = - zdradflx(jl)         ! radiative flux jump         [k m/s]
!          pextr2(jl,5) = zwtventr               ! entrainment flux = we * dtv [k m/s]
!          pextr2(jl,7) = zkh / zmgeom(jl)*rg    ! we                          [m/s]
!          pextr2(jl,6) = zdtv(jl)               ! dtv                         [k]
!          pextr2(jl,7) = zdtv(jl) * rg / zmgeom(jl)
!        endif
          
      endif



!          diffusion coefficient for heat for postprocessing only in (m2/s)

      pkh(jl,jk) = pcfh(jl,jk) / zcfnc1


      !rn output
      if (lldiag) then
!        pextr2(jl,19) = zwecutop(jl)           !top entrainment rate
!        pextr2(jl,20) = zentrtop * pricui(jl)  !entrainment efficiency
        pextra(jl,jk,3) = pcfh(jl,jk)   / zcfnc1   !k [m2/s]
      endif  
      
    enddo

!***
  enddo !jk
!***

! if (lhook) call dr_hook('vdfexcu',1,zhook_handle)
end subroutine vdfexcu
