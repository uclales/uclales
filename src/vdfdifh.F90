subroutine vdfdifh(&
 & kidia  , kfdia  , klon   , klev   , kdraft , ktop   , ktiles, &
 & ptmst  , pextshf, pextlhf, ldsfcflx, &
 & pfrti  , pssrflti,pslrfl , pemis  , pevapsnw, &
 & pslm1  , ptm1   , pqm1   , pqtm1  , paphm1 , &
 & pcfh   , pcfhti , pcfqti , pmflx  , psluh  , pqtuh  , &
 & ptdif  , pqdif  , pcptsti, pqsti  , pcairti, pcsatti, &
 & pdqsti , ptskti , ptskrad, ptsm1m , ptsnow , ptice  , psst, &
 & ptsktip1, pslge , pte    , pqte, &
 & pjq    ,pssh    ,pslh    , pstr   ,pg0, &
 & ldmultisolv)  
!     ------------------------------------------------------------------

!**   *vdfdifh* - does the implicit calculation for diffusion of s. l.

!     derived from vdiff (cy34) by
!     a.c.m. beljaars       e.c.m.w.f.    10-11-89

!     obukhov-l update      acmb          26/03/90.
!     skin t cleaning       p. viterbo    15-11-96.
!     tile boundary cond.   acmb          20-11-98.
!     surface ddh for tiles p. viterbo    17-05-2000.
!     new tile coupling, 
!     ddh moved to vdfmain  a. beljaars   2-05-2003.
!     mass flux terms,
!     flux b.c. for scm,
!     moist generalization  a. beljaars/m. ko"hler 3-12-2004. 
!     removed option for linearized    p. lopez  02/06/2005
!     physics (now called separately)   

!     purpose
!     -------

!     solve tridiagonal matrices for diffusion of dry static energy
!     and moisture; in so doing, it also solves the skin temperature
!     equation.

!     interface
!     ---------

!     *vdfdifh* is called by *vdfmain*

!     input parameters (integer):

!     *kidia*        start point
!     *kfdia*        end point
!     *klev*         number of levels
!     *klon*         number of grid points per packet
!     *ktop*         index for boundary layer top

!     input parameters (real):

!     *ptmst*        double time step (single at 1th step)
!     *pfrti*        fraction of surface area covered by tiles
!     *pslm1*        generalized liquid water static energy    at t-1
!                    (note: in lin/adj physics = dry static energy)
!     *ptm1*         temperature at t-1
!     *pqm1*         specific humidity      at t-1
!     *pqtm1*        specific total water   at t-1
!     *paphm1*       pressure at t-1
!     *pcfh*         prop. to exch. coeff. for heat (c,k-star in doc.)
!     *pcfhti*       idem for heat (surface layer only)
!     *pcfqti*       idem for moisture (surface layer only)
!     *pmflx*        massflux
!     *psluh*        updraft generalized liquid water static energy at half level
!     *pqtuh*        updraft specific total water at half level
!     *pcptsti*      dry static enrgy at surface
!     *pqsti*        saturation q at surface
!     *pcairti*      multiplication factor for q at lowest model level
!                    for surface flux computation
!     *pcsatti*      multiplication factor for qs at surface
!                    for surface flux computation
!     *pdqsti*       d/dt (pqs)
!     *pssrflti*     net solar radiation at the surface, for each tile
!     *pslrfl*       net thermal radiation at the surface
!     *pemis*        model surface longwave emissivity
!     *pevapsnw*     evaporation from snow under forest
!     *ptskti*       skin temperature at t-1
!     *ptskrad*      skin temperature of latest full radiation timestep
!     *ptsm1m*       top soil layer temperature
!     *ptsnow*       snow temperature 
!     *ptice*        ice temperature (top slab)
!     *psst*         (open) sea surface temperature
!     *pslge*        generalized dry static energy tendency
!                    (note: in lin/adj physics = temperature tendency)
!     *pqte*         total water tendency
!                    (note: in lin/adj physics = humidity tendency)

!     output parameters (real):

!     *ptdif*        slg-double-tilde devided by alfa
!     *pqdif*        qt-double-tilde devided by alfa
!     *ptsktip1*     skin temperature at t+1
!     *pjq*          surface moisture flux                      (kg/m2s)
!     *pssh*         surface sensible heat flux                 (w/m2)
!     *pslh*         surface latent heat flux                   (w/m2)
!     *pstr*         surface net thermal radiation              (w/m2)
!     *pg0*          surface ground heat flux (solar radiation  (w/m2)
!                    leakage is not included in this term)

!     additional parameters for flux boundary condtion (in scm model):

!     *ldsfcflx*     if .true. flux boundary condtion is used 
!     *pextshf*      specified sensible heat flux (w/m2)
!     *pextlhf*      specified latent heat flux (w/m2)

!     method
!     ------

!     *lu*-decomposition (downward scan), followed by skin-temperature
!     solver, and back substitution (upward scan).

!     externals.
!     ----------

!     *vdfdifh* calls:
!         *surfseb*

!     ------------------------------------------------------------------
! use garbage, only : surfseb
use parkind1  ,only : jpim     ,jprb
! use yomhook   ,only : lhook    ,dr_hook
use yoevdf   , only : rvdifts
use yos_cst  , only : rcpd     ,rg     ,rlstt  ,rlvtt
use yoethf   , only : rvtmp2
use yos_exc  , only : leocwa   ,leocco

implicit none
! 
! interface
! #include "surfseb.h"
! end interface


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: ktiles 
integer(kind=jpim),intent(in)    :: kidia 
integer(kind=jpim),intent(in)    :: kfdia 
integer(kind=jpim),intent(in)    :: kdraft
integer(kind=jpim),intent(in)    :: ktop 
real(kind=jprb)   ,intent(in)    :: ptmst 
real(kind=jprb)   ,intent(in)    :: pextshf(klon) 
real(kind=jprb)   ,intent(in)    :: pextlhf(klon) 
logical           ,intent(in)    :: ldsfcflx

!-- tiles, input --
real(kind=jprb)   ,intent(in)    :: pfrti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pssrflti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pcfhti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pcfqti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pcptsti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pqsti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pcairti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pcsatti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: pdqsti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    :: ptskti(klon,ktiles) 

real(kind=jprb)   ,intent(in)    :: pslrfl(klon) 
real(kind=jprb)   ,intent(in)    :: pemis(klon) 
real(kind=jprb)   ,intent(in)    :: pevapsnw(klon) 
real(kind=jprb)   ,intent(in)    :: pslm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: ptm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqtm1(klon,klev) 
real(kind=jprb)   ,intent(in)    :: paphm1(klon,0:klev) 
real(kind=jprb)   ,intent(in)    :: pcfh(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pmflx(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(in)    :: psluh(klon,0:klev,kdraft) 
real(kind=jprb)   ,intent(in)    :: pqtuh(klon,0:klev,kdraft) 

real(kind=jprb)   ,intent(in)    :: ptskrad(klon) 
real(kind=jprb)   ,intent(in)    :: ptsm1m(klon) 
real(kind=jprb)   ,intent(in)    :: ptsnow(klon) 
real(kind=jprb)   ,intent(in)    :: ptice(klon) 
real(kind=jprb)   ,intent(in)    :: psst(klon) 
  
real(kind=jprb)   ,intent(in)    :: pslge(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pte(klon,klev) 
real(kind=jprb)   ,intent(in)    :: pqte(klon,klev) 

!-- non-tiles, output --
real(kind=jprb)   ,intent(out)   :: ptdif(klon,klev) 
real(kind=jprb)   ,intent(out)   :: pqdif(klon,klev) 

!-- tiles, output --
real(kind=jprb)   ,intent(out)   :: ptsktip1(klon,ktiles)
real(kind=jprb)   ,intent(out)   :: pjq(klon,ktiles) 
real(kind=jprb)   ,intent(out)   :: pssh(klon,ktiles) 
real(kind=jprb)   ,intent(out)   :: pslh(klon,ktiles) 
real(kind=jprb)   ,intent(out)   :: pstr(klon,ktiles) 
real(kind=jprb)   ,intent(out)   :: pg0(klon,ktiles) 

logical           ,intent(in)    :: ldmultisolv

!*         0.2    local variables

real(kind=jprb) ::    zaa(klon,klev) ,zbb(klon,klev) ,zcc(klon,klev) ,&
                    & ztyy(klon,klev),zqyy(klon,klev),zgam(klon,klev),&
                    & z1dp(klon,klev)
real(kind=jprb) ::    z1bet(klon)    ,&
                    & zaql(klon)     ,zbql(klon)     ,zasl(klon)     ,&
                    & zbsl(klon)     ,zsl(klon)      ,zql(klon)
real(kind=jprb) ::    ztsrf(klon,ktiles)  ,zrhochu(klon,ktiles)  ,& 
                    & zrhocqu(klon,ktiles),zjs(klon,ktiles)      ,&
                    & zssk(klon,ktiles)   ,ztsk(klon,ktiles)  

integer(kind=jpim) :: jk, jl, jt, jd, jd0, jd1

real(kind=jprb) ::    zqdp, zqsp1, ztpfac2, ztpfac3,&
                    & zcsnq, zcsns, zcoef1,zcoef2
! real(kind=jprb) ::    zhook_handle

real(kind=jprb) ::    zmflx(klon,0:klev), zmsl(klon,0:klev), zmqt(klon,0:klev)

!     ------------------------------------------------------------------

! if (lhook) call dr_hook('vdfdifh',0,zhook_handle)

ztpfac2=1.0_jprb/rvdifts
ztpfac3=1-ztpfac2


!     ------------------------------------------------------------------

!*         1.     full model physics with moist mass flux pbl
!                 -------------------------------------------

!*         0.9    precalculation of multiple mass-flux terms

  do jk=ktop,klev
    do jl=kidia,kfdia
      zmflx(jl,jk) = 0.0_jprb
      zmsl (jl,jk) = 0.0_jprb
      zmqt (jl,jk) = 0.0_jprb
    enddo
  enddo
  
  if (ldmultisolv) then
    jd0 = 2      !multiple resolved updrafts
    jd1 = kdraft
  else
    jd0 = 1      !only updraft #1 (bulk)
    jd1 = 1
  endif
  
  do jd=jd0,jd1
    do jk=ktop,klev-1
      do jl=kidia,kfdia
        zmflx(jl,jk) = zmflx(jl,jk) + pmflx(jl,jk,jd)
        zmsl (jl,jk) = zmsl (jl,jk) + pmflx(jl,jk,jd)*psluh(jl,jk,jd)
        zmqt (jl,jk) = zmqt (jl,jk) + pmflx(jl,jk,jd)*pqtuh(jl,jk,jd)
      enddo
    enddo
  enddo

!*         1.0    setting of the matrix a, b and c.

do jk=ktop+1,klev-1
  do jl=kidia,kfdia
    z1dp(jl,jk)=1.0_jprb/(paphm1(jl,jk)-paphm1(jl,jk-1))
    zaa(jl,jk) =(-pcfh(jl,jk-1)-0.5_jprb*zmflx(jl,jk-1))*z1dp(jl,jk)
    zcc(jl,jk) =(-pcfh(jl,jk)  +0.5_jprb*zmflx(jl,jk)  )*z1dp(jl,jk)
    zbb(jl,jk) =1.0_jprb+(pcfh(jl,jk-1)+pcfh(jl,jk)&
     & -0.5_jprb*(zmflx(jl,jk-1)-zmflx(jl,jk)))*z1dp(jl,jk)  
  enddo
enddo

!          1.0a   the surface boundary condition

do jl=kidia,kfdia
  z1dp(jl,klev)=1.0_jprb/(paphm1(jl,klev)-paphm1(jl,klev-1))
  zcc(jl,klev) =0.0_jprb
  zaa(jl,klev) =        (-pcfh(jl,klev-1)-0.5_jprb*zmflx(jl,klev-1))*z1dp(jl,klev)
  zbb(jl,klev) =1.0_jprb+(pcfh(jl,klev-1)-0.5_jprb*zmflx(jl,klev-1))*z1dp(jl,klev)  
enddo

!          1.0b   the top boundary condition    

  do jl=kidia,kfdia
    z1dp(jl,ktop)=1.0_jprb/(paphm1(jl,ktop)-paphm1(jl,ktop-1))
    zaa(jl,ktop) =0.0_jprb
    zcc(jl,ktop) =         (-pcfh(jl,ktop)+0.5_jprb*zmflx(jl,ktop))*z1dp(jl,ktop)
    zbb(jl,ktop) =1.0_jprb+( pcfh(jl,ktop)+0.5_jprb*zmflx(jl,ktop))*z1dp(jl,ktop)
  enddo


!*         1.1    setting of right hand sides.

do jk=ktop+1,klev-1
  do jl=kidia,kfdia
    ztyy(jl,jk) = ztpfac2 * pslm1(jl,jk) &
     & + ptmst * pslge(jl,jk) &
     & + ztpfac2 *(zmsl(jl,jk)-zmsl(jl,jk-1)) *z1dp(jl,jk)  
    zqyy(jl,jk) = ztpfac2 * pqtm1(jl,jk) &
     & + ptmst * pqte(jl,jk) &
     & + ztpfac2 *(zmqt(jl,jk)-zmqt(jl,jk-1)) *z1dp(jl,jk)  
  enddo
enddo

!          1.1a   surface

jk=klev
do jl=kidia,kfdia
  ztyy(jl,jk) = ztpfac2 * pslm1(jl,jk) &
   & + ptmst * pslge(jl,jk) &
   & + ztpfac2 *(-zmsl(jl,jk-1)) *z1dp(jl,jk)  
  zqyy(jl,jk) = ztpfac2 * pqtm1(jl,jk) &
   & + ptmst * pqte(jl,jk) &
   & + ztpfac2 *(-zmqt(jl,jk-1)) *z1dp(jl,jk)  
enddo

!          1.1b   top

jk=ktop
do jl=kidia,kfdia
  ztyy(jl,jk) = ztpfac2 * pslm1(jl,jk) &
   & + ptmst * pslge(jl,jk) &
   & + ztpfac2 *(zmsl(jl,jk)) *z1dp(jl,jk)  
  zqyy(jl,jk) = ztpfac2 * pqtm1(jl,jk) &
   & + ptmst * pqte(jl,jk) &
   & + ztpfac2 *(zmqt(jl,jk)) *z1dp(jl,jk)  
enddo

! 
! !*         1.2    add moisture flux from snow from tile 7 as explicit term
! 
! jk=klev
! do jl=kidia,kfdia
!   zcsnq=rg*ptmst*pfrti(jl,7)*pevapsnw(jl)*z1dp(jl,jk)
!   zcsns=rcpd*rvtmp2*ptskti(jl,7)*zcsnq
! 
!   ztyy(jl,jk)=ztyy(jl,jk)-zcsns
!   zqyy(jl,jk)=zqyy(jl,jk)-zcsnq
! enddo

!*         1.3    top layer elimination.

do jl=kidia,kfdia
  z1bet(jl)=1.0_jprb/zbb(jl,ktop)
  ptdif(jl,ktop)=ztyy(jl,ktop)*z1bet(jl)
  pqdif(jl,ktop)=zqyy(jl,ktop)*z1bet(jl)
enddo

!*         1.4    elimination for middle layers.

do jk=ktop+1,klev-1
  do jl=kidia,kfdia
    zgam(jl,jk)=zcc(jl,jk-1)*z1bet(jl)
    z1bet(jl)=1.0_jprb/(zbb(jl,jk)-zaa(jl,jk)*zgam(jl,jk))
    ptdif(jl,jk)=(ztyy(jl,jk)-zaa(jl,jk)*ptdif(jl,jk-1))*z1bet(jl)
    pqdif(jl,jk)=(zqyy(jl,jk)-zaa(jl,jk)*pqdif(jl,jk-1))*z1bet(jl)
  enddo
enddo

!*         1.5    bottom layer, linear relation between lowest
!                 model level s and q and fluxes.

do jl=kidia,kfdia
  zgam(jl,klev)=zcc(jl,klev-1)*z1bet(jl)
  z1bet(jl)=1.0_jprb/(zbb(jl,klev)-zaa(jl,klev)*zgam(jl,klev))
  zbsl(jl)=(ztyy(jl,klev)-zaa(jl,klev)*ptdif(jl,klev-1))*rvdifts*z1bet(jl)
  zbql(jl)=(zqyy(jl,klev)-zaa(jl,klev)*pqdif(jl,klev-1))*rvdifts*z1bet(jl)
  zqdp=1.0_jprb/(paphm1(jl,klev)-paphm1(jl,klev-1))
  zasl(jl)=-rg*ptmst*zqdp*rvdifts*z1bet(jl)
  zaql(jl)=zasl(jl)
enddo
 
!*         1.6    prepare array's for call to surface energy
!                 balance routine
 
if (leocwa .or. leocco) then
  ztsrf(kidia:kfdia,1)=ptskti(kidia:kfdia,1)
else
  ztsrf(kidia:kfdia,1)=psst(kidia:kfdia)
endif
ztsrf(kidia:kfdia,2)=ptice(kidia:kfdia)
ztsrf(kidia:kfdia,3)=ptsm1m(kidia:kfdia)
ztsrf(kidia:kfdia,4)=ptsm1m(kidia:kfdia)
ztsrf(kidia:kfdia,5)=ptsnow(kidia:kfdia)
ztsrf(kidia:kfdia,6)=ptsm1m(kidia:kfdia)
ztsrf(kidia:kfdia,7)=ptsnow(kidia:kfdia)
ztsrf(kidia:kfdia,8)=ptsm1m(kidia:kfdia)

zcoef1=1.0_jprb/(rg*rvdifts*ptmst)
do jt=1,ktiles
  do jl=kidia,kfdia
    zrhochu(jl,jt)=pcfhti(jl,jt)*zcoef1
    zrhocqu(jl,jt)=pcfqti(jl,jt)*zcoef1

    !if (pfrti(jl,jt).gt.0._jprb) then
    !if (jl.eq.kidia) then
    !write(0,'(a,i3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3)') &
    !  & 'vdfdifh: tile props before call to vdfsurfseb: jt=',jt, &
    !  & ' pcptsti=',pcptsti(jl,jt), &
    !  & ' ptskti=',ptskti(jl,jt), &
    !  & ' pqsti=',pqsti(jl,jt), &
    !  & ' pdqsti=',pdqsti(jl,jt), &
    !  & ' zrhochu=',zrhochu(jl,jt), &
    !  & ' zrhocqu=',zrhocqu(jl,jt), &
    !  & ' pcairti=',pcairti(jl,jt), &
    !  & ' pcsatti=',pcsatti(jl,jt), &
    !  & ' pssrflti=',pssrflti(jl,jt), &
    !  & ' pfrti=',pfrti(jl,jt), &
    !  & ' ztsrf=',ztsrf(jl,jt)
    !endif
    !endif

  enddo
enddo
! 
! !*         1.7    call to surface energy balance routine
! !                 remember: output is extrapolated in time

! call vdfsurfseb   (kidia=kidia,kfdia=kfdia,klon=klon,ktiles=ktiles,&
!  & psskm1m=pcptsti,ptskm1m=ptskti,pqskm1m=pqsti,&
!  & pdqsdt=pdqsti,prhochu=zrhochu,prhocqu=zrhocqu,&
!  & palphal=pcairti,palphas=pcsatti,&
!  & pssrfl=pssrflti,pfrti=pfrti,ptsrf=ztsrf,&
!  & pslrfl=pslrfl,ptskrad=ptskrad,pemis=pemis,&
!  & pasl=zasl,pbsl=zbsl,paql=zaql,pbql=zbql,&
!  !out
!  & pjs=zjs,pjq=pjq,pssk=zssk,ptsk=ztsk,&
!  & pssh=pssh,pslh=pslh,pstr=pstr,pg0=pg0,&
!  & psl=zsl,pql=zql)  

! !*         1.7a   add snow evaporation to fluxes
! 
! pjq (kidia:kfdia,7)=pjq (kidia:kfdia,7)+pevapsnw(kidia:kfdia)
! pslh(kidia:kfdia,7)=pslh(kidia:kfdia,7)+pevapsnw(kidia:kfdia)*rlstt
! 
! !*         1.7b   flux boundary condition for 1d model (fluxes in w/m2)
! !                 (over-write output of surfseb)

!write(0,*) 'vdfdifh:',ldsfcflx

if (ldsfcflx) then
  do jt=1,ktiles
    do jl=kidia,kfdia
      zjs(jl,jt)=pextshf(jl)+rcpd*ptskti(jl,jt)*rvtmp2*pextlhf(jl)/rlvtt
      pjq(jl,jt)=pextlhf(jl)/rlvtt

!rn bug      zssk(jl,jt)=zbsl(jl)-zjs(jl,jt)*(zasl(jl)-1.0_jprb/zrhochu(jl,jt)) 
      zssk(jl,jt)=zbsl(jl)+zjs(jl,jt)*(zasl(jl)-1.0_jprb/zrhochu(jl,jt)) 
      ztsk(jl,jt)=zssk(jl,jt)/(rcpd*(1.+rvtmp2*pqsti(jl,jt)))
      pssh(jl,jt)=pextshf(jl)
      pslh(jl,jt)=pextlhf(jl)
      pstr(jl,jt)=pslrfl(jl)
      pg0 (jl,jt)=pextshf(jl)+pextlhf(jl)+pslrfl(jl)+pssrflti(jl,jt)
    enddo
  enddo
  do jl=kidia,kfdia
    zsl(jl)=zjs(jl,1)*zasl(jl)+zbsl(jl)   !apparently tile 1 is used for prescribed fluxes!
    zql(jl)=pjq(jl,1)*zaql(jl)+zbql(jl)
  enddo
endif

!*         1.7c   compute parameters at new time level 

do jt=1,ktiles

  do jl=kidia,kfdia
    ptsktip1(jl,jt)=ztpfac2*ztsk(jl,jt)+ztpfac3*ptskti(jl,jt)
    zqsp1=pqsti(jl,jt)+pdqsti(jl,jt)*(ztsk(jl,jt)-ptskti(jl,jt))

    !if (pfrti(jl,jt).gt.0._jprb) then
    !if (jl.eq.kidia) then
    !write(0,'(a,i3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3)') &
    !  & 'vdfdifh: sfc en budg: jt=',jt, &
    !  & ' pfrti=',pfrti(jl,jt),' pssh=',pssh(jl,jt),' pslh=',pslh(jl,jt), &
    !  & ' pslrfl=',pslrfl(jl)+pssrflti(jl,jt),' pg0=',pg0(jl,jt), &
    !  & ' izrhochu=',1.0_jprb/zrhochu(jl,jt), ' zjs=',zjs(jl,jt)
    !endif     
    !endif

  enddo
enddo

!*         1.7d   copy lowest model solution from surfseb

zcoef2=1.0_jprb/rvdifts
do jl=kidia,kfdia
  ptdif(jl,klev)=zsl(jl)*zcoef2
  pqdif(jl,klev)=zql(jl)*zcoef2
enddo

!*         1.8    back-substitution.

do jk=klev-1,ktop,-1
  do jl=kidia,kfdia
    ptdif(jl,jk)=ptdif(jl,jk)-zgam(jl,jk+1)*ptdif(jl,jk+1)
    pqdif(jl,jk)=pqdif(jl,jk)-zgam(jl,jk+1)*pqdif(jl,jk+1)
  enddo
enddo
  
! if (lhook) call dr_hook('vdfdifh',1,zhook_handle)
end subroutine vdfdifh
