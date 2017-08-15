!options xopt(hsfun)
! 
module vdf
implicit none
!logical :: ledmf = .false.
logical :: ledmf = .true.


contains
function flip(a)
  real, dimension(:), intent(in) :: a
  real, dimension(size(a)) :: flip
  integer :: k
  do k=1,size(a)
    flip(k) = a(1+size(a)-k)
  end do
  
end function flip

function interpolate(a)
  real, dimension(:), intent(in) :: a
  real, dimension(size(a)-1) :: interpolate
  integer :: k
  do k=1,size(a) - 1
    interpolate(k) = 0.5*(a(k) + a(k+ 1) )
  end do
  
end function interpolate
  
subroutine vdfouter(sst,kstp)
use parkind1, only : jprb, jpim
use grid,  only : nzp, nxp, nyp, &
                & liquid, a_up, a_vp, a_wp, a_rp, a_tp, &    !watch out! a_rp = q_t, a_tp = theta_l
                & a_theta, a_pexnr, pi0, pi1, dt, press,zm, zt, &
                & a_tt, a_rt, &
                & wt_sfc, wq_sfc, uw_sfc, vw_sfc, &
                & dn0, th00, &
                & iradtyp, sfc_albedo, a_lflxu, a_lflxd, a_sflxu, a_sflxd, albedo,&
                & a_tskin, kh_les!, a_tsoil, a_phiw
use defs, only : g,alvl, cp, cpr, p00, tmelt
use srfc, only : zrough
use util, only : get_avg3, get_avg
use stat, only : sflg, updtst,updtst_ts, stat_edmf
use yomct0, only : rextlhf, rextshf
use yoevdfs, only : firsttime, suvdfs
use yos_veg, only : ntilesifs, susveg
use yos_cst, only : rlvtt
use thrm, only : fll_tkrs
use lsmdata, only: ksoilmax !,z0m,z0h,Qnetm,tsoilm,phiwm,tskinm
!use step, only: istp
! implicit none

!-- declare local processing variables --
real, dimension(nzp, nxp, nyp) :: a_scr1
real, dimension(nzp) :: a1
real :: a2
!real, dimension((nxp-4)*(nyp-4),nzp) :: a1
!real, dimension((nxp-4)*(nyp-4)) :: a2
real, dimension(nzp)   :: presf
real, dimension(nzp-1) :: presh
integer :: i, j, k, pcog, pnog, klevx,kcog

real, intent(in) :: sst
integer, intent(in) :: kstp

!-- declare input variables for vdfmain --
real(kind=jprb), dimension((nxp-4)*(nyp-4),nzp-2) :: &  !x,y,full levels
  & pum1,pvm1,pwm1,ptm1,pqm1,plm1,pim1,pam1,papm1,psobeta,psoteu,psotev, &
  & pvdis,pvdisg,ptofdu,ptofdv, &
  & pgeom1,phrlw,phrsw,pie,pkh,ple,pqe,pae,pte,pvar,pvol, pthle, &
  & pvom, pzinv,pldiff,pwuavg
real(kind=jprb), dimension((nxp-4)*(nyp-4),ntilesifs) :: &    !x,y,ntilesifs -- this is the IFS surface tile system
  & pfrti,palbti,pssrflti,pustrti, pvstrti, pahfsti, pevapti, ptskti
real(kind=jprb), dimension((nxp-4)*(nyp-4),0:nzp-2) :: &   !x,y,half levels
  & paphm1,pgeoh, &                                        !half level pressure and geopotential height
  & pfplvl,pfplvn,pfhpvl,pfhpvn, &                         !EDMF precipitation fluxes
  & pdiftq,pdiftl,pdifti,pdifts,pstrtu,pstrtv, &           !EDMF turbulent fluxes
  & pstrsou,pstrsov                                        !subgrid orography momentum flux
real(kind=jprb), dimension((nxp-4)*(nyp-4),nzp-1,0):: &   !tracers. last dimension: ktrac (# of tracers, currently 0) 
  & pcm1   
real(kind=jprb), dimension((nxp-4)*(nyp-4)):: &   !x,y (real scalars)  surface properties
  & psst , ptsnow, ptice, pz0m, pz0h, ptskm1m, ptskrad, pevapsnw, &          
  & pslrfl,pssrfl,pemis, &
  & pkhfl, pkqfl, pkmfl, &                        !local surface fluxes
  & pcvl, pcvh, &                                 !vegetation fractions
  & psigflt, pchar, pucurr, pvcurr
real(kind=jprb), dimension((nxp-4)*(nyp-4),ksoilmax) :: &   !x,y,ksoilmax -- properties at soil levels
  & ptsam1m, pwsam1m
integer(kind=jpim), dimension((nxp-4)*(nyp-4)) :: &   !x,y (integers scalars)
  & kpbltype, khpbln, kvartop,&
  & ktvl,ktvh                                         !vegetation types
integer(kind=jpim),parameter :: &   !diagnostics - # of fields
  & kfldx=152, kfldx2=100
real(kind=jprb) :: &    !diagnostics - x,y and x,y,z arrays
  & pextr2((nxp-4)*(nyp-4),kfldx2), pextra((nxp-4)*(nyp-4),nzp-1,kfldx), pk_les((nxp-4)*(nyp-4),nzp-1)

!write(0,'(a)' ) 'vdfouter: start'

pextr2 = 0.
pextra = 0.
pk_les = 0.
pnog = (nxp-4)*(nyp-4)
klevx = nzp-1    !output array should be large enough to accomodate half levels

if(ledmf) then

  !if(.not. sflg) return

  !-- fill pressure arrays --

  presf=p00*((pi0+pi1)/cp)**cpr    !LES height and pressure arrays are 1D (see grid.f90)
  presh = 0.
  presh = interpolate(presf)
  
  !do k=1,nzp
  !  write(0,'(a,i3,a,f10.3,a,f8.3,a,f8.3)' ) 'vdfouter:   k=',k,&
  !    & '  presf=',presf(k), ' ', zt(k), ' ',zm(k)
  !enddo
  !do k=1,nzp-1
  !  write(0,'(a,i3,a,f10.3,a,f8.3,a,f8.3)' ) 'vdfouter:   k=',k,&
  !    & '  presh=',presh(k)
  !enddo
  

  !-- fill x,y dependent fields --
  do j=3,nyp-2
  do i=3,nxp-2

     pcog = (i-2) + (j-3) * (nxp-4)

     !-- 3D state variables --
     pum1(pcog,:) = flip(a_up(2:nzp-1,i,j))
     pvm1(pcog,:) = flip(a_vp(2:nzp-1,i,j))
     pwm1(pcog,:) = flip(a_wp(2:nzp-1,i,j))
     
     !from theta_l to T
     call fll_tkrs(nzp,nxp,nyp,a_tp+th00,a_pexnr,pi0,pi1,a_scr1)                   !T_l = theta_l * Exner
     a_scr1(2:nzp-1,i,j) = a_scr1(2:nzp-1,i,j) + liquid(2:nzp-1,i,j) * alvl / cp   !T   = T_l + L * q_l / c_p
     ptm1(pcog,:) = flip(a_scr1(2:nzp-1,i,j))

     !from qt to qv
     a_scr1 = 0
     a_scr1(2:nzp-1,i,j) = a_rp(2:nzp-1,i,j) - liquid(2:nzp-1,i,j)
     pqm1(pcog,:) = flip(a_scr1(2:nzp-1,i,j))    

     !condensate
     plm1(pcog,:) = flip(liquid(2:nzp-1,i,j))    
     pim1(pcog,:) = 0.

     !cloud fraction
     a_scr1 = 0
     where(liquid>0) a_scr1 = 1
     pam1(pcog,:) = 0.!flip(a_scr1(2:nzp-1,i,j))    !no cloud in column

     !-- full level pressure and height --
     papm1(pcog,:) = flip(presf(2:nzp-1))    !presf is at zt-levels
     pgeom1(pcog,:) = g*flip(zt(2:nzp-1))    ! zt(1)=-20, zt(2)=20, zt(3)=60., etc

     !-- half level pressure and height --
     paphm1(pcog,0:nzp-2) = flip(presh(1:nzp-1))   ! presh is at zm-levels: 1st is in between presf(1) and presf(2)
     pgeoh (pcog,0:nzp-2) = g*flip(zm(1:nzp-1))    ! zm(1)=0., zm(2)=40., zm(3)=80., etc.

     !-- local surface fluxes, for EDMF plume initialization --
     !   Note: upwards = negative sign, ECMWF IFS convention
     pkhfl(pcog) =  - wt_sfc(i,j)*cp*(dn0(1)+dn0(2))*0.5       !sensible heat   
     pkqfl(pcog) =  - wq_sfc(i,j)*alvl*(dn0(1)+dn0(2))*0.5     !latent heat
     pkmfl(pcog) =  ( uw_sfc(i,j)**2. + vw_sfc(i,j)**2. )**0.5 !momentum

     pk_les(pcog,:) = flip(kh_les(2:nzp-1,i,j))

  end do  !i
  end do  !j

  !-- prescribed domain-averaged surface fluxes (option in EDMF)
  !   Note: upwards = negative sign, ECMWF IFS convention
  rextshf = - get_avg(1,nxp,nyp,1,wt_sfc)*cp*(dn0(1)+dn0(2))*0.5     !sensible heat
  rextlhf = - get_avg(1,nxp,nyp,1,wq_sfc)*alvl*(dn0(1)+dn0(2))*0.5   !latent heat

  !-- initialize/reset --
  ptsnow   = 0._jprb    !d  snow temperature
  ptice    = 0._jprb    !d  ice temperature (top slab)
  psst     = 0._jprb    !d  (open) sea surface temperature
  pevapsnw = 0._jprb    !d  evaporation of snow under trees

  pz0m    = 0._jprb
  pz0h    = 0._jprb
  ptskm1m = 0._jprb
  ptskrad = 0._jprb  

  pcm1    = 0._jprb

  psobeta = 0._jprb
  psoteu  = 0._jprb
  psotev  = 0._jprb

  pte     = 0._jprb
  pqe     = 0._jprb
  ple     = 0._jprb
  pie     = 0._jprb
  pae     = 0._jprb
  pvom    = 0._jprb
  pvol    = 0._jprb  

  pstrtu  = 0._jprb
  pstrtv  = 0._jprb
  ptofdu  = 0._jprb
  ptofdv  = 0._jprb

  !radiation
  pssrfl  = 0._jprb
  pslrfl  = 0._jprb    !d  net thermal radiation at the surface
  pemis   = 0._jprb    !d  model surface longwave emissivity
  phrlw   = 0._jprb
  phrsw   = 0._jprb

  !tiles
  pfrti    = 0._jprb   !d  fraction
  palbti   = 0._jprb   !   albedo
  ptskti   = 0._jprb   !d  skin temperature at t-1  
  pssrflti = 0._jprb   !d  net sw rad at sfc
  pustrti  = 0._jprb
  pvstrti  = 0._jprb
  pahfsti  = 0._jprb
  pevapti  = 0._jprb

  !soil levels
  ptsam1m  = 0._jprb   !d  temperature of soil layers   
  pwsam1m  = 0._jprb                                    

  !vegetation
  ktvl = 1          !low vegetation type
  ktvh = 1          !high veg type
  pcvl = 0.5_jprb   !low veg cover
  pcvh = 0.5_jprb   !high veg cover
  
  !orography
  psigflt = 0._jprb

  !ocean properties
  pchar  = 0._jprb
  pucurr = 0._jprb
  pvcurr = 0._jprb

  !##### temporary for running over ocean (should be generalized later) ##################
  !      variables used in vdfdifh.F90 and vdfsurfexcdriver.F90
  !      to be made dependent on LES sfc scheme!

  pz0m        = zrough    !ocean values, see lsm in src.f90  !z0m
  pz0h        = zrough                                       !z0h
  psst        = sst
  ptsnow      = tmelt     !ze smelten de kazen!
  ptice       = tmelt
  ptskm1m     = sst       ! skin T at t-1  lsmdata: tskinm  !grid: a_tskin   (tskinm=a_tskin at t-1, see line 951 in srfc.f90, lsm)
  ptskrad     = sst       ! skin T used in calculation of radiative fluxes (pslrfl)

  !-- tile properties  --
  pfrti (:,1) = 1._jprb   !tile fraction (sea) = 1
  palbti(:,1) = sfc_albedo          !tile broadband albedo (needed in vdfsurfexcdriver)   !grid: albedo, sfc_albedo
  ptskti(:,1) = sst   
  !do j=3,nyp-2
  !do i=3,nxp-2
  !  pcog = (i-2) + (j-3) * (nxp-4)
  !  ptskti(pcog,1) = a_tskin(i,j)     !RN a_tskin etc only allocated if isfctyp==5 (see define_vars in grid.f90) 
  !enddo
  !enddo
  !write(0,'(a,f10.3,a,f10.3)' ) &
  !  & 'vdfouter:   sfc tiles: sfc_albedo=', sfc_albedo!, ' a_tskin=',a_tskin(3,3)   

  !-- soil at t-1 --
  ptsam1m     = sst       !lsmdata: tsoilm  !grid: a_tsoil     (tsoilm = a_tsoil at previous timestep, see line 1036 in srfc.f90, lsm)
  !pwsam1m     =           !lsmdata: phiwm   !grid: a_phiw     only used in vdfsurfexcdriver

  !-- radiative fluxes --
  !note: (a_lflxu, a_lflxd) only defined for iradtyp > 0, and 
  !      (albedo, a_sflxu, a_sflxd) only defined for iradtyp > 2     
  !           see grid.f90 (subr define_vars) and forc.f90 (subr forcings)
  if (iradtyp>2) then
    do j=3,nyp-2
    do i=3,nxp-2
      pcog = (i-2) + (j-3) * (nxp-4)
      pssrfl(pcog) = a_sflxd(2,i,j) - a_sflxu(2,i,j)   !use level 2, see e.g. subroutine surfacerad (in rad_driver.f90)
      pslrfl(pcog) = a_lflxd(2,i,j) - a_lflxu(2,i,j)
    enddo
    enddo
  endif

  !k=2  !zt(1) = -20, zm(1) = 0, zt(2) = 20
  !i=3
  !j=3
  !write(0,'(a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3)' ) &
  !  & 'vdfouter:   sfc: zt=',zt(k), ' zm=',zm(k) &
  !  & , ' a_lflxd', a_lflxd(k,i,j), ' a_lflxu=', a_lflxu(k,i,j)
  !  & !, ' albedo=',sfc_albedo, ' a_sflxd', a_sflxd(k,i,j), ' a_sflxu=', a_sflxu(k,i,j)
 
  !-- vegetation props --
  ktvl(:) = 15         !ocean
  ktvh(:) = 15         !ocean
  pcvl(:) = 0.5_jprb   !low veg cover
  pcvh(:) = 0.5_jprb   !high veg cover

  !-- tile surface heat fluxes and stresses --
  do j=3,nyp-2
  do i=3,nxp-2
    pcog = (i-2) + (j-3) * (nxp-4)
    pustrti(:,1)  = uw_sfc(i,j)
    pvstrti(:,1)  = vw_sfc(i,j)
  enddo
  enddo
  pahfsti(:,1)  = pkhfl(:)           !should be in W/m2
  pevapti(:,1)  = pkqfl(:) / rlvtt   !should be in kg m-2 s-1

  !##### temporary #################################################################
  
else
  stop
end if   !ledmf

if (firsttime) then
  call suvdfs
  call susveg
endif
firsttime = .false.

call vdfmain(   &
  cdconf = 'a', &
  kidia  = 1,     &
  kfdia  = pnog,     &
  klon   = pnog,     &
  klev   = nzp-2,    &
  klevs  = 1,     & !or 0?
!  kstep  = istp,    & 
  kstep  = kstp,    & 
  ktiles = ntilesifs,     &
  ktrac  = 0, & !number of tracers does this include qt, thl? - no.
  klevsn = 0,     & 
  klevi  = 0,     &
  kdhvtls = 0,    &
  kdhftls = 0,    &
  kdhvtss = 0, &
  kdhftss = 0, &
  kdhvtts = 0, &
  kdhftts = 0, &
  kdhvtis = 0, & 
  kdhftis = 0, &
  ptsphy  = dt, &
! 
  ktvl = ktvl, &   !low vegetation type
  ktvh = ktvh, &   !high veg type
  pcvl = pcvl, &   !low veg cover
  pcvh = pcvh, &   !high veg cover
!
  kcnt = 1, &
!
  psigflt = psigflt, &
!
  pum1    = pum1, & 
  pvm1    = pvm1, & 
  pwm1    = pwm1, &
  ptm1    = ptm1, & 
  pqm1    = pqm1, & 
  plm1    = plm1, & 
  pim1    = pim1, & 
  pam1    = pam1, &
  pcm1    = pcm1, & 
!
  paphm1  = paphm1, & 
  papm1   = papm1, & 
  pgeom1  = pgeom1, & 
  pgeoh   = pgeoh, & 
!
  ptskm1m = ptskm1m, &
  ptsam1m = ptsam1m, & 
  pwsam1m = pwsam1m, &
!
  pssrfl  = pssrfl, &
  pslrfl  = pslrfl, & 
  pemis   = pemis, &
  phrlw   = phrlw, & 
  phrsw   = phrsw, & 
!
  ptsnow  = ptsnow, &
  ptice   = ptice, &
  pevapsnw = pevapsnw, &
!
  pkhfl   = pkhfl, &
  pkqfl   = pkqfl, &
  pkmfl   = pkmfl, &
!
  psst    = psst, & 
!
  pchar  = pchar,  & 
  pucurr = pucurr,  & 
  pvcurr = pvcurr,  &  
!
  ptskrad = ptskrad , & 
  psobeta = psobeta, &
  psoteu  = psoteu, &
  psotev = psotev, &
  ldnodecp = (/.false./), &
! roughness lengths
  pz0m    = pz0m, & !inout
  pz0h    = pz0h, & 
! tendencies
  pte     = pte, & 
  pqe     = pqe, & 
  ple     = ple, & 
  pie     = pie, & 
  pae     = pae, & 
  pvom    = pvom, & 
  pvol    = pvol, &  
  pthle   = pthle, &
! tiles
  pfrti    = pfrti, &   
  palbti   = palbti, &
  pssrflti = pssrflti, &   
  pustrti  = pustrti, &
  pvstrti  = pvstrti, &
  ptskti   = ptskti, &
  pevapti  = pevapti, &
  pahfsti  = pahfsti, &
!
  pvar    = pvar, & !gridbox average qt variance
  pzinv   = pzinv, & 
  kpbltype = kpbltype, & 
! precip fluxes
  pfplvl  = pfplvl, & 
  pfplvn  = pfplvn, & 
  pfhpvl  = pfhpvl, & 
  pfhpvn  = pfhpvn, & 
!
  pwuavg  = pwuavg, & 
  khpbln   = khpbln, & 
  kvartop   = kvartop, & 
  pldiff   = pldiff, & 
  kfldx2 = kfldx2, klevx = klevx, kfldx = kfldx, &
  pextr2 = pextr2, &
  pextra = pextra, &
  pk_les = pk_les, &
  pvdis = pvdis , pvdisg = pvdisg, &
  pdiftq = pdiftq, pdiftl = pdiftl, pdifti = pdifti, pdifts = pdifts, &
  pstrtu = pstrtu, pstrtv = pstrtv, ptofdu = ptofdu , ptofdv = ptofdv, &
  pstrsou = pstrsou, pstrsov =pstrsov, &
  pkh   = pkh & 
  )


!--- add vdf tendencies to uclales tendencies ---
do j=3,nyp-2
do i=3,nxp-2
!j=3
!i=3

    pcog = (i-2) + (j-3) * (nxp-4)

    do k=2,nzp-1

      kcog = nzp-k

      !write(0,'(a,f10.3,a,f10.3,a,f15.7,a,f15.7,a,f15.7,a,f15.7,a,f15.7)' ) &
      !      & 'vdfouter end: zt=',zt(k), ' zm=',zm(k), &
      !      & ' pte=' , 3600.*pte(pcog,kcog), &
      !      & ' pthle=' , 3600.*pthle(pcog,kcog), &
      !      & ' pqe=' , 3600.*1000.*pqe(pcog,kcog), &
      !      & ' ple=' , 3600.*1000.*ple(pcog,kcog), &
      !      & ' pqte=', 3600.*1000.*(pqe(pcog,kcog)+ple(pcog,kcog))
      
      a_tt(k,i,j) = a_tt(k,i,j) + pthle(pcog,kcog)
      a_rt(k,i,j) = a_rt(k,i,j) + pqe  (pcog,kcog) + ple(pcog,kcog)

    end do

enddo
enddo


!output of edmf variables
  if (sflg) then
    call stat_edmf(pnog, pextr2, pextra, kfldx2, kfldx, klevx)
  end if

!write(0,'(a)' ) 'vdfouter: end'

end subroutine vdfouter


subroutine vdfmain    ( cdconf, &
 & kidia  , kfdia  , klon   , klev   , klevs  , kstep  , ktiles , &
 & ktrac  , klevsn , klevi  , kdhvtls, kdhftls, kdhvtss, kdhftss, &
 & kdhvtts, kdhftts, kdhvtis, kdhftis, &
 & ptsphy , ktvl   , ktvh   , kcnt   , pcvl   , pcvh   , psigflt, &
 & pum1   , pvm1   , pwm1   , ptm1   , pqm1   , plm1   , pim1   , pam1   , pcm1   , &
 & paphm1 , papm1  , pgeom1 , pgeoh  , ptskm1m, ptsam1m, pwsam1m, &
 & pssrfl , pslrfl , pemis  , phrlw  , phrsw  , &
 & ptsnow , ptice  , &
 & pkhfl  , pkqfl  , pkmfl  , &
 & psst   , pfrti  , palbti , pwlmx  , &
 & pchar  , pucurr , pvcurr , ptskrad, pcflx  , &
 & psoteu , psotev , psobeta, pvervel, &
 & pz0m   , pz0h   , &
 & pvdis  , pvdisg , pvar   , &
 & pzinv  , khpbln , kvartop, &
 & pssrflti,pevapsnw, pwuavg , ldnodecp,kpbltype, pldiff, &
 & pfplvl , pfplvn , pfhpvl , pfhpvn , &
 & pextr2 , kfldx2 , pextra , pk_les , klevx  , kfldx  , &
 & pte    , pqe    , ple    , pie    , pae    , pvom   , pvol   , pthle, &
 & ptenc  , ptske1 , &
 & pustrti, pvstrti, pahfsti, pevapti, ptskti , &
 & pdifts , pdiftq , pdiftl , pdifti , pstrtu , pstrtv , ptofdu , ptofdv, &
 & pstrsou, pstrsov,   pkh   &
 & )

!***

!**   *vdfmain* - does the vertical exchange of u,v,slg,qt by turbulence.

!     j.f.geleyn       20/04/82   original  
!     c.a.blondin      18/12/86
!     a.c.m. beljaars  20/10/89   ifs-version (technical revision of cy34)
!     a.c.m. beljaars  26/03/90   obukhov-l update 
!     a.c.m. beljaars  30/09/98   surface tiles 
!     p. viterbo       17/05/2000 surface ddh for tiles
!     d. salmond       15/10/2001 fullimp mods
!     s. abdalla       27/11/2001 passing zi/l to waves
!     a. beljaars       2/05/2003 new tile coupling     
!     p.viterbo        24/05/2004 change surface units
!     m. ko"hler        3/12/2004 moist advection-diffusion
!     a. beljaars       4/04/2005 turb. orogr. drag
!     a. beljaars      30/09/2005 include subgr. oro. in solver
!     a.beljaars       31/03/2005 introduction of ocean current b.c.
!     p. viterbo       17/06/2005 surf external library
!     m. ko"hler        6/06/2006 single column model option (lscmec) 
!     r. neggers       24/05/2006 pbl dualm scheme + bimodal cloud scheme
 
!     purpose.
!     --------

!          this routine computes the physical tendencies of the four
!     prognostic variables u,v,t and q due to the vertical exchange by
!     turbulent (= non-moist convective) processes. these tendencies are
!     obtained as the difference between the result of an implicit
!     time-step starting from values at t-1 and these t-1 values. all
!     the diagnostic computations (exchange coefficients, ...) are done
!      from the t-1 values. as a by-product the roughness length over sea
!     is updated accordingly to the *charnock formula. heat and moisture
!     surface fluxes and their derivatives against ts, ws and wl
!     (the latter will be later weighted with the snow factor in
!     *vdiff*), later to be used for soil processes treatment, are also
!     computed as well as a stability value to be used as a diagnostic
!     of the depth of the well mixed layer in convective computations.

!     interface.
!     ----------
!          *vdiff* takes the model variables at t-1 and returns the values
!     for the prognostic time t+1 due to vertical diffusion.
!     the model variables, the model dimensions and the diagnostics data
!     are passed as subroutine arguments. constants that do not change
!     during a model run (e.g. physical constants, switches etc.) are
!     stored in a single common block *yomvdf*, which is initialized
!     by set-up routine *suvdf*.

!     parameter     description                                   units
!     ---------     -----------                                   -----
!     input parameters (integer):

!    *kidia*        start point
!    *kfdia*        end point
!    *klev*         number of levels
!    *klon*         number of grid points per packet
!    *klevs*        number of soil layers
!    *kstep*        current time step index
!    *ktiles*       number of tiles (i.e. subgrid areas with different 
!                   of surface boundary condition)
!    *ktrac*        number of tracers
!    *klevsn*       number of snow layers (diagnostics) 
!    *klevi*        number of sea ice layers (diagnostics)
!    *kdhvtls*      number of variables for individual tiles
!    *kdhftls*      number of fluxes for individual tiles
!    *kdhvtss*      number of variables for snow energy budget
!    *kdhftss*      number of fluxes for snow energy budget
!    *kdhvtts*      number of variables for soil energy budget
!    *kdhftts*      number of fluxes for soil energy budget
!    *kdhvtis*      number of variables for sea ice energy budget
!    *kdhftis*      number of fluxes for sea ice energy budget

!    *ktvl*         vegetation type for low vegetation fraction
!    *ktvh*         vegetation type for high vegetation fraction

!    *kcnt*         index of vdf sub steps.

!     input parameters (logical)

!     input parameters (real)

!    *ptsphy*       time step for the physics

!     input parameters at t-1 or constant in time (real):

!    *pcvl*         low vegetation cover                          -  
!    *pcvh*         high vegetation cover                         -  
!    *psigflt*      standard deviation of filtered orography      m
!    *pum1*         x-velocity component                          m/s
!    *pvm1*         y-velocity component                          m/s
!    *pwm1*         w-velocity component                          m/s
!    *ptm1*         temperature                                   k
!    *pqm1*         specific humidity                             kg/kg
!    *plm1*         specific cloud liquid water                   kg/kg
!    *pim1*         specific cloud ice                            kg/kg
!    *pam1*         cloud fraction                                1
!    *pcm1*         tracer concentration                          kg/kg
!    *papm1*        pressure on full levels                       pa
!    *paphm1*       pressure on half levels                       pa
!    *pgeom1*       geopotential                                  m2/s2
!    *pgeoh*        geopotential at half levels                   m2/s2
!    *ptskm1m*      skin temperature at t-1                       k
!    *ptsam1m*      surface temperature                           k
!    *pwsam1m*      soil moisture all layers                      m**3/m**3
!    *pssrfl*       net shortwave radiation flux at surface       w/m2
!    *pslrfl*       net longwave radiation flux at surface        w/m2
!    *pemis*        model surface longwave emissivity
!    *phrlw*        longwave heating rate                         k/s
!    *phrsw*        shortwave heating rate                        k/s
!    *ptsnow*       snow temperature                              k
!    *ptice*        ice temperature (top slab)                    k
!    *pkhfl*        surface sensible heat flux                    w/m2
!    *pkqfl*        surface latent heat flux                      w/m2
!    *pkmfl*        surface momentum flux                         m2/s2
!    *psst*         (open) sea surface temperature                k
!    *pfrti*        tile fractions                                (0-1)
!            1 : water                  5 : snow on low-veg+bare-soil
!            2 : ice                    6 : dry snow-free high-veg
!            3 : wet skin               7 : snow under high-veg
!            4 : dry snow-free low-veg  8 : bare soil
!    *palbti*       broadband albedo for tile fractions
!    *pwlmx*        maximum skin reservoir capacity               kg/m**2
!    *pchar*        "equivalent" charnock parameter
!    *pucurr*       ocean current x_component      
!    *pvcurr*       ocean current y_component      
!    *ptskrad*      skin temperature of latest full radiation
!                      timestep                                   k
!    *pcflx*        tracer surface flux                           kg/(m2 s)
!    *psoteu*       explicit part of u-tendency from subgrid orography scheme    
!    *psotev*       explicit part of v-tendency from subgrid orography scheme     
!    *psobeta*      implicit part of subgrid orography 

!    *pvervel*      vertical velocity

!     input parameters (logical):

!     contributions to budgets (output,real):

!    *pvdis*        turbulent dissipation                         w/m2
!    *pvdisg*        subgrid orography dissipation                 w/m2
!    *pahflev*      latent heat flux  (snow/ice free part)        w/m2
!    *pahflsb*      latent heat flux  (snow/ice covered part)     w/m2

!     updated parameters (real):

!    *pte*          temperature tendency                          k/s
!    *pqe*          moisture tendency                             kg/(kg s)
!    *ple*          liquid water tendency                         kg/(kg s)
!    *pie*          ice water tendency                            kg/(kg s)
!    *pae*          cloud fraction tendency                       1/s)
!    *pvom*         meriodinal velocity tendency (du/dt)          m/s2
!    *pvol*         latitude tendency            (dv/dt)          m/s2
!    *ptenc*        tracer tendency                               kg/(kg s)
!    *ptske1*       skin temperature tendency                     k/s
!    *pz0m*         aerodynamic roughness length                  m
!    *pz0h*         roughness length for heat                     m

!     updated parameters for tiles (real): 

!    *pustrti*      surface u-stress                              n/m2 
!    *pvstrti*      surface v-stress                              n/m2
!    *pahfsti*      surface sensible heat flux                    w/m2
!    *pevapti*      surface moisture flux                         kg/m2/s
!    *ptskti*       skin temperature                              k

!     output parameters (real):

!    *pfplvl*       pbl precipitation flux as rain                kg/(m**2*s)
!    *pfplvn*       pbl precipitation flux as snow                kg/(m**2*s)
!    *pfhpvl*       enthalpy flux of pbl precipitation as rain    j/(m**2*s)
!    *pfhpvn*       enthalpy flux of pbl precipitation as snow    j/(m**2*s)

!    *pldiff*       contrib to pbl condensate by passive clouds   kg/kg

!    *pfwsb*        evaporation of snow                           kg/(m**2*s)
!    *pu10m*        u-component wind at 10 m                      m/s
!    *pv10m*        v-component wind at 10 m                      m/s
!    *pt2m*         temperature at 2m                             k
!    *pd2m*         dew point temperature at 2m                   k
!    *pq2m*         specific humidity at 2m                       kg/kg
!    *pgust*        gust at 10 m                                  m/s
!    *pzidlwv*      zi/l used for gustiness in wave model         m/m
!                   (note: positive values of zi/l are set to zero)
!    *pblh*         pbl height (dry diagnostic based on ri#)      m
!    *pzinv*        pbl height (moist parcel, not for stable pbl) m
!    *pssrflti*     net shortwave radiation flux at surface, for
!                      each tile                                  w/m2
!    *pevapsnw*     evaporation from snow under forest            kg/(m2*s)
!    *pstrtu*       turbulent flux of u-momemtum            kg*(m/s)/(m2*s)
!    *pstrtv*       turbulent flux of v-momemtum            kg*(m/s)/(m2*s)
!    *ptofdu*       tofd comp. of turbulent flux of u-momemtum   kg*(m/s)/(m2*s)
!    *ptofdv*       tofd comp. of turbulent flux of v-momemtum   kg*(m/s)/(m2*s)
!    *pdifts*       turbulent flux of heat                         j/(m2*s)
!    *pdiftq*       turbulent flux of specific humidity           kg/(m2*s)
!    *pdiftl*       turbulent flux of liquid water                kg/(m2*s)
!    *pdifti*       turbulent flux of ice water                   kg/(m2*s)
!    *pstrsou*      subgrid orography flux of u-momemtum    kg*(m/s)/(m2*s)
!    *pstrsov*      subgrid orography flux of v-momemtum    kg*(m/s)/(m2*s)

!    *pkh*          turb. diff. coeff. for heat above surf. lay.  (m2/s)
!                   in surface layer: ch*u                        (m/s)

!     additional parameters for flux boundary condition (in scm model):

!    *lsfcflx*      if .true. flux boundary condtion is used 
!    *rextshf*      specified sensible heat flux [w/m2]
!    *rextlhf*      specified latent heat flux [w/m2]

!     method.
!     -------

!          first an auxialiary variable cp(q)t+gz is created on which
!     the vertical diffusion process will work like on u,v and q. then
!     along the vertical and at the surface, exchange coefficients (with
!     the dimension of a pressure thickness) are computed for momentum
!     and for heat (sensible plus latent). the letters m and h are used
!     to distinguish them and the computation is the result of a
!     conditional merge between the stable and the unstable case
!     (depending on the sign of the *richardson bulk number).
!          in the second part of the routine the implicit linear
!     systems for u,v first and t,q second are solved by a *gaussian
!     elimination back-substitution method. for t and q the lower
!     boundary condition depends on the surface state.
!     over land, two different regimes of evaporation prevail:
!     a stomatal resistance dependent one over the vegetated part
!     and a soil relative humidity dependent one over the
!     bare soil part of the grid mesh.
!     potential evaporation takes place over the sea, the snow
!     covered part and the liquid water covered part of the
!     grid mesh as well as in case of dew deposition.
!          finally one returns to the variable temperature to compute
!     its tendency and the later is modified by the dissipation's effect
!     (one assumes no storage in the turbulent kinetic energy range) and
!     the effect of moisture diffusion on cp. z0 is updated and the
!     surface fluxes of t and q and their derivatives are prepared and
!     stored like the difference between the implicitely obtained
!     cp(q)t+gz and cp(q)t at the surface.

!     externals.
!     ----------

!     *vdfmain* calls sucessively:
!         *surfexcdriver*
!         *vdfexcu*
!         *vdftofdc*
!         *vdfdifm*
!         *vdfdifh*
!         *vdfdifc*
!         *vdfincr*
!         *vdfsdrv*
!         *vdfppcfl*
!         *vdfupdz0*

!     reference.
!     ----------

!          see vertical diffusion's part of the model's documentation
!     for details about the mathematics of this routine.

!     ------------------------------------------------------------------
use garbage, only : foealfa!, foeewm, foeldcpm
use parkind1  ,only : jpim     ,jprb
! use yomhook   ,only : lhook    ,dr_hook

use yomct0   , only : lsfcflx  ,lsfcflx_les, rextshf  ,rextlhf,  rextmomf 
use yoevdf   , only : rvdifts  ,lldiag
use yos_cst   , only : rg       ,rd       ,&
                    & rcpd     ,retv     ,rlvtt    ,rlstt    ,rtt 
use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les    ,&
                    & r4ies    ,r5les    ,r5ies    ,rvtmp2   ,r5alvcp  ,&
                    & r5alscp  ,ralvdcp  ,ralsdcp  ,rtwat    ,rtice    ,&
                    & rticecu  ,rtwat_rtice_r      ,rtwat_rticecu_r  
use yomjfh   , only : n_vmass
use yoephy   , only : lvdftrac

use defs, only : cpr, p00

implicit none


!*         0.1    global variables

integer(kind=jpim),intent(in)    :: klon 
integer(kind=jpim),intent(in)    :: klev 
integer(kind=jpim),intent(in)    :: klevs 
integer(kind=jpim),intent(in)    :: ktiles 
integer(kind=jpim),intent(in)    :: ktrac
integer(kind=jpim),intent(in)    :: klevsn 
integer(kind=jpim),intent(in)    :: klevi 
integer(kind=jpim),intent(in)    :: kdhvtls 
integer(kind=jpim),intent(in)    :: kdhftls 
integer(kind=jpim),intent(in)    :: kdhvtss 
integer(kind=jpim),intent(in)    :: kdhftss 
integer(kind=jpim),intent(in)    :: kdhvtts 
integer(kind=jpim),intent(in)    :: kdhftts 
integer(kind=jpim),intent(in)    :: kdhvtis 
integer(kind=jpim),intent(in)    :: kdhftis 
character(len=1)  ,intent(in)    ,optional:: cdconf 
integer(kind=jpim),intent(in)    ,optional:: kidia
integer(kind=jpim),intent(in)    ,optional:: kfdia
integer(kind=jpim),intent(in)    ,optional:: kstep 
real(kind=jprb)   ,intent(in)    ,optional:: ptsphy 
integer(kind=jpim),intent(in)    ,optional:: ktvl(klon) 
integer(kind=jpim),intent(in)    ,optional:: ktvh(klon) 
integer(kind=jpim),intent(in)    ,optional:: kcnt
real(kind=jprb)   ,intent(in)    ,optional:: pcvl(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pcvh(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: psigflt(klon)
real(kind=jprb)   ,intent(in)    ,optional:: pum1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pvm1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pwm1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: ptm1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pqm1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: plm1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pim1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pam1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pcm1(klon,klev,ktrac) 
real(kind=jprb)   ,intent(in)    ,optional:: paphm1(klon,0:klev)
real(kind=jprb)   ,intent(in)    ,optional:: papm1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pgeom1(klon,klev)
real(kind=jprb)   ,intent(in)    ,optional:: pgeoh(klon,0:klev)
real(kind=jprb)   ,intent(in)    ,optional:: ptskm1m(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: ptsam1m(klon,klevs) 
real(kind=jprb)   ,intent(in)    ,optional:: pwsam1m(klon,klevs) 
real(kind=jprb)   ,intent(in)    ,optional:: pssrfl(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pslrfl(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pemis(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: phrlw(klon,klev) 
real(kind=jprb)   ,intent(in)    ,optional:: phrsw(klon,klev) 
real(kind=jprb)   ,intent(in)    ,optional:: ptsnow(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: ptice(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pkhfl(klon)
real(kind=jprb)   ,intent(in)    ,optional:: pkqfl(klon)
real(kind=jprb)   ,intent(in)    ,optional:: pkmfl(klon)
real(kind=jprb)   ,intent(in)    ,optional:: psst(klon) 

real(kind=jprb)   ,intent(in)    ,optional:: pfrti(klon,ktiles) 
real(kind=jprb)   ,intent(in)    ,optional:: palbti(klon,ktiles) 
real(kind=jprb)   ,intent(out)   ,optional:: pssrflti(klon,ktiles) 
real(kind=jprb)   ,intent(inout) ,optional:: pustrti(klon,ktiles) 
real(kind=jprb)   ,intent(inout) ,optional:: pvstrti(klon,ktiles) 
real(kind=jprb)   ,intent(inout) ,optional:: pahfsti(klon,ktiles) 
real(kind=jprb)   ,intent(inout) ,optional:: pevapti(klon,ktiles) 
real(kind=jprb)   ,intent(inout) ,optional:: ptskti(klon,ktiles) 

real(kind=jprb)   ,intent(in)    ,optional:: pwlmx(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pchar(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pucurr(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pvcurr(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: ptskrad(klon) 
real(kind=jprb)   ,intent(in)    ,optional:: pcflx(klon,ktrac)
real(kind=jprb)   ,intent(in)    ,optional:: psoteu(klon,klev) 
real(kind=jprb)   ,intent(in)    ,optional:: psotev(klon,klev) 
real(kind=jprb)   ,intent(in)    ,optional:: psobeta(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pz0m(klon) 
real(kind=jprb)   ,intent(inout) ,optional:: pz0h(klon) 
real(kind=jprb)   ,intent(out)   ,optional:: pvdis(klon) 
real(kind=jprb)   ,intent(out)   ,optional:: pvdisg(klon) 
real(kind=jprb)   ,intent(out)   ,optional:: pvar(klon,klev)
real(kind=jprb)   ,intent(out)   ,optional:: pzinv(klon)
integer(kind=jpim),intent(out)   ,optional:: khpbln(klon)
real(kind=jprb)   ,intent(out)   ,optional:: pevapsnw(klon) 
real(kind=jprb)   ,intent(out)   ,optional:: pfplvl(klon,0:klev)
real(kind=jprb)   ,intent(out)   ,optional:: pfplvn(klon,0:klev)
real(kind=jprb)   ,intent(out)   ,optional:: pfhpvl(klon,0:klev)
real(kind=jprb)   ,intent(out)   ,optional:: pfhpvn(klon,0:klev)
real(kind=jprb)   ,intent(out)   ,optional:: pwuavg(klon)
logical           ,intent(in)    ,optional:: ldnodecp(klon)
integer(kind=jpim),intent(out)   ,optional:: kpbltype(klon)
real(kind=jprb)   ,intent(out)   ,optional:: pldiff(klon,klev) 
integer(kind=jpim),intent(out)   ,optional:: kvartop(klon)
real(kind=jprb)   ,intent(inout) ,optional:: pte(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pqe(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: ple(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pie(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pae(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pvom(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pvol(klon,klev) 
real(kind=jprb)   ,intent(out)   ,optional:: pthle(klon,klev) 
real(kind=jprb)   ,intent(inout) ,optional:: ptenc(klon,klev,ktrac)
real(kind=jprb)   ,intent(inout) ,optional:: ptske1(klon) 
real(kind=jprb)   ,intent(out)   ,optional:: pdifts(klon,0:klev) 
real(kind=jprb)   ,intent(out)   ,optional:: pdiftq(klon,0:klev) 
real(kind=jprb)   ,intent(out)   ,optional:: pdiftl(klon,0:klev) 
real(kind=jprb)   ,intent(out)   ,optional:: pdifti(klon,0:klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pstrtu(klon,0:klev) 
real(kind=jprb)   ,intent(inout) ,optional:: pstrtv(klon,0:klev) 
real(kind=jprb)   ,intent(inout) ,optional:: ptofdu(klon) 
real(kind=jprb)   ,intent(inout) ,optional:: ptofdv(klon)
real(kind=jprb)   ,intent(out)   ,optional:: pstrsou(klon,0:klev) 
real(kind=jprb)   ,intent(out)   ,optional:: pstrsov(klon,0:klev) 
real(kind=jprb)   ,intent(out)   ,optional:: pkh(klon,klev) 
real(kind=jprb)   ,intent(in)    ,optional:: pvervel(klon,klev) 
!          diagnostic output
integer(kind=jpim),intent(in)     :: kfldx2, klevx, kfldx
real(kind=jprb)   ,intent(inout)  ,optional:: pextr2(klon,kfldx2), pextra(klon,klevx,kfldx), pk_les(klon,klevx)

!*         0.2    local variables

real(kind=jprb) ::    zdiftqt(klon,0:klev), zdiftslg(klon,0:klev) 

real(kind=jprb) ::    zcptgz(klon,klev) , zcfm(klon,klev)   , zcfh(klon,klev)   ,&
                    & zudif(klon,klev)  , zvdif(klon,klev)  ,&
                    & zqtdif(klon,klev) , zslgdif(klon,klev),&
                    & zslgm1(klon,klev) , zqtm1(klon,klev)  , zqte(klon,klev)   ,&
                    & zslge(klon,klev)  , ztofdc(klon,klev) , zsoc(klon,klev)   , &
                    & zthlm1(klon,klev) , zexni(klon,klev)
real(kind=jprb) ::    zkhfl(klon)       , zkqfl(klon)       , zkmfl(klon)  
real(kind=jprb) ::    zqea(klon,klev)   , zlea(klon,klev)   , ziea(klon,klev)   ,&
                    & zqtea(klon,klev)  , zslgea(klon,klev) , zaea(klon,klev)   ,&
                    & ztea(klon,klev)   , zuea(klon,klev)   , zvea(klon,klev)   ,&
                    & zslgewodis(klon,klev)  
real(kind=jprb) ::    zqtep(klon,klev)  , zslgep(klon,klev), zdzrhoi

!rn --- variables associated with ED(MF)^n scheme -------------------------

integer(kind=jpim), parameter   ::   idraft = 11    ! nr of resolved bins in size density
!integer(kind=jpim), parameter   ::   idraft = 21    ! nr of resolved bins in size density

real(kind=jprb) :: zr_max(klon)          !prognostic: maximum size in density
real(kind=jprb) :: ztkeint(klon,idraft)  !prognostic: size density of vertically integrated tke

logical :: llmultisolv    !switch for multi-plume mode in solver (recommended: T)

real(kind=jprb) ::  zfrac (klon,0:klev,idraft) , &
                  & zmflx (klon,0:klev,idraft) , &
                  & zqtuh (klon,0:klev,idraft) , &   
                  & zslguh(klon,0:klev,idraft) , &
                  & zwuh  (klon,0:klev,idraft) , &
                  & zmflxm(klon,0:klev,idraft) , &
                  & zuuh  (klon,0:klev,idraft) , &
                  & zvuh  (klon,0:klev,idraft)
                  
real(kind=jprb) ::  zfracc(klon,0:klev) , &
                  & zqcc  (klon,0:klev)

real(kind=jprb) ::  zucurr(klon) , &
                  & zvcurr(klon) , &
                  & ztaux (klon) , &
                  & ztauy (klon)

real(kind=jprb) ::   zzptop(klon,idraft) , &
                   & zzplcl(klon,idraft)

integer(kind=jpim) ::   iptop(klon,idraft) , &
                      & iplcl(klon,idraft) , &
                      & iplzb(klon,idraft)  

!rn --------------------------------------------------------------------

real(kind=jprb) ::    zz0mw(klon)       , zz0hw(klon)       , zz0qw(klon)       ,&
                    & zblend(klon)      , zfblend(klon)
real(kind=jprb) ::    zzcpts(klon)      , zzqsa(klon)       , zzbuom(klon)      ,&
                    & zzzdl(klon)
real(kind=jprb) ::    ztupd(klon,klev)  , zqupd(klon,klev)  , zlupd(klon,klev)  ,&
                    & ziupd(klon,klev)  , zqtupd(klon,klev) , zliupd(klon,klev) ,&
                    & zslgupd(klon,klev), zaupd(klon,klev)  , zthlupd(klon,klev)
real(kind=jprb) ::    ztauxcg(klon,klev), ztauycg(klon,klev)
real(kind=jprb) ::    zqsvar(klon,klev) , zanew(klon,klev)  , zlnew(klon,klev)  
real(kind=jprb) ::    zsvfluxcld(klon,0:klev)               , zsvfluxsub(klon,0:klev),&
                    & zsvflux(klon,0:klev),zbuoypos(klon)   , zbuoyneg(klon)    ,&
                    & zdzh(klon,0:klev)
real(kind=jprb) ::    zalfa1            , zalfa2            , zdelq             ,&
                    & zcorqs(klon,klev) , zdqsdtemp(klon,klev)

real(kind=jprb) ::    zcldbase(klon)    , zcldtop(klon)

real(kind=jprb) ::    zpflxusum(klon,0:klev)

!rn --- vdf qt variance budget & bimodal cloud scheme variables ---------

real(kind=jprb) ::    zvargen           , zvartrans         , zvardiss          ,&
                    & ztau(klon,klev)   , ztaunew           , &
                    & zdqtdz            , zwqtf             , &
                    & zwqt2(klon,0:klev), zqt2uh(klon,0:klev), & 
                    & zqtup(klon,klev)  , zslgup(klon,klev) , &
                    & zaup(klon,klev)   , zsigqt2up(klon,klev)
                 
real(kind=jprb) ::    zcldfrac(klon,klev),zqlav(klon,klev)

!rn ---------------------------------------------------------------------

real(kind=jprb) ::    ztskintiold(klon,ktiles)

real(kind=jprb) ::    zextshf(klon)     , zextlhf(klon)    , zextmomf(klon)

real(kind=jprb) ::    zcptsti(klon,ktiles), zqsti(klon,ktiles)  ,&
                    & zdqsti(klon,ktiles) , zcsatti(klon,ktiles),&
                    & zcairti(klon,ktiles), zcfhti(klon,ktiles) ,&
                    & zcfqti(klon,ktiles) , zahflti(klon,ktiles),&
                    & ztsktip1(klon,ktiles)
real(kind=jprb) ::    zstr(klon,ktiles)   , zg0(klon,ktiles)

logical ::            llrundry(klon), llpbl(klon,klev), &
                      lltropdamp

integer(kind=jpim) :: itop, jk, jl, jd, kcap,jkk, jt

real(kind=jprb) ::    zgdph, zrhoh(klon,0:klev),ztmst, zrg, zrtmst

logical ::            llsfcflx

real(kind=jprb) ::    zhu1

real(kind=jprb) ::    zalfaw(klon,klev), zfacw, zfaci, zfac, zesdp, zcor

real(kind=jprb) ::    zdummy2, zdummy3, zcfnc1, zmgeom

! real(kind=jprb) ::    zhook_handle

! interface
! #include "surfexcdriver.h"
! #include "surfpp.h"
! end interface
! 
! #include "vdfdifh.intfb.h"
! #include "vdfdifh5.intfb.h"
! #include "vdfdifm.intfb.h"
! #include "vdfdifm2.intfb.h"
! #include "vdfdifc.intfb.h"
! #include "vdfdpbl.intfb.h"
! #include "vdfexcu.intfb.h"
! #include "vdfhghtn.intfb.h"
! #include "vdfincr.intfb.h"
! #include "vdffblend.intfb.h"
! #include "vdfcloud.intfb.h"
! #include "vdftofdc.intfb.h"
! 
! #include "fcttre.h"



!     ------------------------------------------------------------------

!*         1.     initialize constants
!                 --------------------

! if (lhook) call dr_hook('vdfmain',0,zhook_handle)

!write(0,'(a)' ) 'vdfmain: start'

ztmst       = ptsphy
zrtmst      = 1.0_jprb/ptsphy    ! optimization
zrg         = 1.0_jprb/rg        !     -"-
llrundry(:) = .false.  ! option to run dry updrafts with no condensation
!do jl=kidia,kfdia
! if ( .not. ldnodecp(jl) )  llrundry(jl) = .true. ! run dry for real vdfmain's
!enddo

!-- switch for eliminating ls tendencies in troposphere (for research purposes only!) --
!lltropdamp = .true.
lltropdamp = .false.

!-- switch for multiple updraft solver mode --
!llmultisolv = .false.   !solve bulk updraft, integrated over size density
llmultisolv = .true.    !solve multiple resolved updrafts (i.e. the ones coming from vdfparcel) -- recommended

!-- temporary: use fixed maximum size in density --
zr_max(:)    = 1000._jprb
ztkeint(:,:) = 0._jprb

!*         1.0  set the external fluxes

!if (lscmec) then
  llsfcflx   = lsfcflx   
  if (lsfcflx_les) then
    ! use t-1 fluxes from LES (vary with time- and location)
    zextshf(:)  = pkhfl(:)
    zextlhf(:)  = pkqfl(:)
    zextmomf(:) = pkmfl(:)
  else
    ! use prescribed fluxes from scm namelist (constant & homogeneous)
    zextshf(:)  = rextshf
    zextlhf(:)  = rextlhf
    zextmomf(:) = rextmomf
  endif
!else
!  llsfcflx   = .false.
!  zextshf(:) = 0.0_jprb
!  zextlhf(:) = 0.0_jprb
!endif


!*         1.1  store initial tendencies for flux calculation
!*              and initialize variable.

do jk=1,klev
  do jl=kidia,kfdia
    zqea(jl,jk)=pqe(jl,jk)
    zlea(jl,jk)=ple(jl,jk)
    ziea(jl,jk)=pie(jl,jk)
    zaea(jl,jk)=pae(jl,jk)
    ztea(jl,jk)=pte(jl,jk)
    zuea(jl,jk)=pvom(jl,jk)
    zvea(jl,jk)=pvol(jl,jk)
  enddo
enddo
zqtep  = 0._jprb
zslgep = 0._jprb
pldiff = 0._jprb

zmflx  = 0._jprb
zslguh = 0._jprb
zqtuh  = 0._jprb
zwuh   = 0._jprb

pfplvl = 0._jprb
pfplvn = 0._jprb
pfhpvl = 0._jprb
pfhpvn = 0._jprb

zcfm = 0._jprb 
zcfh = 0._jprb
pkh  = 0._jprb


!     ------------------------------------------------------------------

!*         2.     new variables s, slg, qt
!*                (at initial time level)
!                 ------------------------

do jk=1,klev
  do jl=kidia,kfdia

!*         2.1  dry static energy cp(q)*t + gz

    zcptgz(jl,jk)  =pgeom1(jl,jk) + ptm1(jl,jk)*rcpd * (1.0_jprb + rvtmp2*pqm1(jl,jk))

!*         2.2  total water and generalized liquid water static energy 
!*              slg = cp*t + gz - lcond*ql - ldep*qi

    zslgm1(jl,jk) = zcptgz(jl,jk) - rlvtt * plm1(jl,jk) - rlstt * pim1(jl,jk)
    zslge(jl,jk)  = rcpd * ( ( 1.0_jprb + rvtmp2 * pqm1(jl,jk) ) * pte(jl,jk)   &  !dcpt/dt
                &                       + rvtmp2 * ptm1(jl,jk)   * pqe(jl,jk) ) &  
                & - rlvtt * ple(jl,jk) - rlstt * pie(jl,jk)                        !dlqli/dt  

    zexni (jl,jk) = ( p00 / papm1(jl,jk) ) ** (rd/rcpd)
    zthlm1(jl,jk) = ( ptm1(jl,jk) - rlvtt * plm1(jl,jk) / rcpd ) * zexni(jl,jk)    !old thl

    zqtm1(jl,jk ) = pqm1(jl,jk) + plm1(jl,jk) + pim1(jl,jk)
    zqte(jl,jk)   = pqe(jl,jk)  + ple(jl,jk)  + pie(jl,jk)             !rad+dyn. qt tendency
    
    zslgea(jl,jk) = zslge(jl,jk)
    zqtea(jl,jk)  = zqte(jl,jk)
    
    if (jl==kfdia ) then
      !write(0,'(a,i3,a,f15.7,a,f15.7)') 'vdfmain 1.0:  jk=',jk,&
      !         & ' zslge:', zslge(jl,jk) * ztmst / rcpd, &
      !         & ' zqte:', 1000._jprb * zqte(jl,jk) * ztmst
      !write(0,'(a,i3,a,f15.7,a,1f15.7)') 'vdfmain 1.0:  jk=',jk,&
      !         & ' zexni:' , zexni (jl,jk), &
      !         & ' zthlm1:', zthlm1(jl,jk)
    endif

  enddo
enddo



!     ------------------------------------------------------------------

!*         3.  compute all surface related quantities
!          ------------------------------------------

 do jl=kidia,kfdia
   do jk = 1,ktiles
     ztskintiold(jl,jk) = ptskti(jl,jk)
   enddo   
 enddo   
 
 call vdfsurfexcdriver(&
  & cdconf, &
  & kidia, kfdia, klon, klevs, ktiles, kstep, &
  & klevsn, klevi, kdhvtls, kdhftls, &
  & kdhvtss, kdhftss, kdhvtts, kdhftts, &
  & kdhvtis, kdhftis, n_vmass, &
  & ptsphy, rvdifts, &
 ! input data, non-tiled
  & ktvl, ktvh, pcvl, pcvh, &
  & pum1(:,klev), pvm1(:,klev), ptm1(:,klev), &
  & pqm1(:,klev), paphm1(:,klev), pgeom1(:,klev), &
  & zcptgz(:,klev), psst, ptskm1m, pchar, &
  & pssrfl, pslrfl, pemis, ptice, ptsnow, &
  & pwlmx, pucurr, pvcurr, &
 ! input data, soil
  & ptsam1m, pwsam1m, &
 ! input data, tiled
  & pfrti, palbti, &
 ! updated data, tiled
  & pustrti, pvstrti, pahfsti, pevapti, &
  & ptskti, &
 ! updated data, non-tiled
  & pz0m, pz0h, &
 ! output data, tiled
  & pssrflti, zqsti, zdqsti, zcptsti, &
  & zcfhti, zcfqti, zcsatti, zcairti, &
 ! output data, non-tiled
  & pkh(:,klev), zcfm(:,klev), zkmfl, zkhfl, &
  & zkqfl, pevapsnw, zz0mw, zz0hw, zz0qw, &
  & zblend, zzcpts, zzqsa, zzbuom, &
  & zzzdl &
  & )

!--- temporary: set diff coefs at lowest level to zero ---
!zcfm(:,klev) = 0.   !see vdfdifm.F90
!zcfh(:,klev) = 0.   !not used, see vdfdifh.F90: Effective C_h at surface is composed of C_h's of all tiles (zcfhti)


!-- convert sfc fluxes from energy (W/m2) to kinematic (.. m/s) units --

!jl = kfdia
!jk = klev
!  write(0,'(a,i3,a,f15.7,a,f15.7)') 'vdfmain: ',jk, &
!    & ' pkh=',pkh(jl,jk), ' zcfm=',zcfm(jl,jk)

do jl=kidia,kfdia
  zrhoh(jl,klev) = paphm1(jl,klev) / ( rd*ptm1(jl,klev)*(1.0_jprb+retv*pqm1(jl,klev)) )
  zkhfl(jl) = zextshf(jl) / ( rcpd  * (1.0_jprb+rvtmp2*pqm1(jl,klev)) * zrhoh(jl,klev) )
  zkqfl(jl) = zextlhf(jl) / ( rlvtt *                                   zrhoh(jl,klev) )
  zkmfl(jl) = zextmomf(jl)
enddo
 
                        
!rn --- surface heat flux check ---
if (.false.) then
do jl=kidia,kfdia
  if (zkhfl(jl).lt.-10._jprb.or.zkqfl(jl).lt.-10._jprb) then
         
!    write(0,'(a,i,f10.6)') 'vdfmain:',jl,ptskm1m(jl)
    write(0,'(a,2f10.6)')   '  vdfmain: sfluxes:',zkhfl(jl),zkqfl(jl)
     
!    do jk = 1,ktiles
!     write(0,'(a,i,2f10.6)')'  vdfmain: tskin:',jk,pfrti(jl,jk),ptskti(jl,jk)
!    enddo   
     
  endif
enddo
endif

 
!rn --- surface heat flux limiters ---
do jl=kidia,kfdia
  zkhfl(jl) = max( zkhfl(jl), -1000.0_jprb / rcpd)
  zkqfl(jl) = max( zkqfl(jl), -2000.0_jprb / rlvtt )
enddo
 
 

!     ------------------------------------------------------------------

!*         4.     exchange coefficients
!                 ---------------------

!*         4.4  computation of the pbl extension


!          set pbl height-index to 1

do jl=kidia,kfdia
  khpbln(jl)=1
  kvartop(jl)  = 0
enddo
itop=1  !itop is used in some solvers. itop=1 means: always integrate over whole atmosphere depth


!*         4.5  boundary layer height for dianostics only
! 
! call vdfdpbl(kidia,kfdia,klon,klev,&
!  & pum1,pvm1,ptm1,pqm1,pgeom1,&
!  & zkmfl,zkhfl,zkqfl,pblh)  


!*         4.6  call the plume model interface

call vdfhghtn (kidia   , kfdia   , klon    , klev    , idraft   , ztmst , kstep, &
             & pum1    , pvm1    , pwm1    , ptm1    , pqm1    , plm1    , pim1    , pam1     , &
             & paphm1  , papm1   , pgeom1  , pgeoh   , pvervel , pqe     , pte      , &
             & zkmfl   , zkhfl   , zkqfl   , zmflx   , zfrac   , &
             & zr_max  , ztkeint , &
! diagnostic output
             & pextr2  , kfldx2  , pextra  , klevx   , kfldx    , &
!              
             & zuuh    , zvuh    , zslguh  , zqtuh   , zwuh  , &
             & zzptop  , iptop   , zzplcl  , iplcl   , iplzb    , &
             & pwuavg  , &
             & zfracc  , zqcc    , pfplvl  , pfplvn  , &
             & ldnodecp, llrundry, kpbltype, zwqt2 , zqt2uh)


!--- set some pbl heights based on i) pbl type and ii) various updraft properties ---
do jl=kidia,kfdia
    
  select case (kpbltype(jl))
    
      case(0)
          !stable pbl
          pzinv(jl)    = 0.0_jprb
          khpbln(jl)   = klev
          zcldbase(jl) = -100._jprb
          zcldtop(jl)  = -100._jprb
          kvartop(jl)  = 0
     
      case(1)
          !dry convective pbl
          pzinv(jl)    = zzptop(jl,1)
          khpbln(jl)   = iptop(jl,1)
          zcldbase(jl) = -100._jprb
          zcldtop(jl)  = -100._jprb
          kvartop(jl)  = iptop(jl,1)

      case(2)
          !stratocumulus
          pzinv(jl)    = zzptop(jl,1)
          khpbln(jl)   = iptop(jl,1)
          zcldbase(jl) = zzplcl(jl,1)
          zcldtop(jl)  = zzptop(jl,1)
          kvartop(jl)  = iptop(jl,1)
	  
      case(3)
          !shallow cumulus
          pzinv(jl)    = zzplcl(jl,1)
          khpbln(jl)   = iplcl(jl,1)+1   !equivalent to definition of top level: use first half level *below* boundary
          zcldbase(jl) = zzplcl(jl,1)
          zcldtop(jl)  = zzptop(jl,1)
          kvartop(jl)  = iptop(jl,1)
		
      case(4)
          !deep cumulus - only do a dry subcloud ml
          pzinv(jl)    = zzplcl(jl,1)
          khpbln(jl)   = iplcl(jl,1)+1   !equivalent to definition of top level: use first half level *below* boundary
          zcldbase(jl) = zzplcl(jl,1)
          zcldtop(jl)  = zzptop(jl,1)
          kvartop(jl)  = iptop(jl,1)
    
  end select 
  
 ! if (jl==kfdia) then
 !   write(0,'(a,i3,a,f10.3,a,i3,a,i3,a,f10.3,a,f10.3)') &
 !            & 'vdfmain 2.0:  kstep=',kstep, ' dt=',ptsphy, &
 !            & ' kvartop=', kvartop(jl), ' kpbltype=', kpbltype(jl), &
 !            & ' zcldbase=', zcldbase(jl), ' zcldtop=', zcldtop(jl)
 ! endif

enddo !jl


!--- set pbl indicator ---   
do jk=1,klev
  do jl=kidia,kfdia
    llpbl(jl,jk) = jk >= kvartop(jl) .and. kvartop(jl)>0
  enddo
enddo


!--- updraft precipitation tendencies ---
do jk=2,klev
  do jl=kidia,kfdia
      
      !calculate flux divergences
      zdzrhoi = rg/( paphm1(jl,jk)-paphm1(jl,jk-1) )
      zqtep(jl,jk)  = -( pfplvl(jl,jk) - pfplvl(jl,jk-1) ) * zdzrhoi &
                    & -( pfplvn(jl,jk) - pfplvn(jl,jk-1) ) * zdzrhoi
      zslgep(jl,jk) =   rlvtt * ( pfplvl(jl,jk) - pfplvl(jl,jk-1) ) * zdzrhoi &
                    & + rlstt * ( pfplvn(jl,jk) - pfplvn(jl,jk-1) ) * zdzrhoi

      !add contributions to total tendencies (to be accounted for in solver vdfdifh)
      zqte(jl,jk)  = zqte(jl,jk)  + zqtep(jl,jk)
      zslge(jl,jk) = zslge(jl,jk) + zslgep(jl,jk)

    !if (jl==kfdia ) then
    !  write(0,'(a,i3,a,f15.7,a,f15.7,a,f15.7,a,f15.7)') 'vdfmain 3.0:  jk=',jk,&
    !           & ' zslge:', zslge(jl,jk) * ztmst / rcpd, &
    !           & ' zslgep:', zslgep(jl,jk) * ztmst / rcpd, &
    !           & ' zqte:', 1000._jprb * zqte(jl,jk) * ztmst, &
    !           & ' zqtep:', 1000._jprb * zqtep(jl,jk) * ztmst
    !endif

  enddo
enddo


!*         4.7  exchange coefficients above the surface layer

call vdfexcu(kidia  , kfdia  , klon   , klev    , idraft  , ztmst  , pz0m   , &
           & phrlw  , phrsw  , &
           & pum1   , pvm1   , ptm1   , pqm1    , plm1    , pim1   , &
           & paphm1 , papm1  , pgeom1 , pgeoh   , zcptgz  , &
           & pextr2 , kfldx2 , pextra , klevx   , kfldx   , &
           & zkmfl  , zkhfl  , zkqfl  , zcfm    , zcfh    , ztauxcg, ztauycg, &
           & zmflx  , kvartop , &
           & pzinv  , khpbln , pkh    , zcldbase, zcldtop , kpbltype, pk_les)  

!do jk=1,klev
!  do jl=kidia,kfdia
!    if (jl==kfdia .and. llpbl(jl,jk)) then
!      write(0,'(a,i3,a,f15.7,a,2f15.7)') 'vdfmain:  jk=',jk,&
!             & ' C_m:', zcfm(jl,jk), ' C_h:', zcfh(jl,jk), pkh(jl,jk)
!    endif
!  enddo
!enddo


!*         4.8  mass flux modifications in bulk plume (for solver stability, in case llmultisolv=.F.)

do jd=1,1

  !-- remove single and/or double massflux layers --
  do jl=kidia,kfdia
    if ( zmflx(jl,klev-2,jd) < 1.e-40_jprb ) then
    !if ( zmflx(jl,klev-3,jd) < 1.e-40_jprb ) then
      zmflx(jl,klev-1,jd) = 0._jprb
      !zmflx(jl,klev-2,jd) = 0._jprb
    endif
  enddo
  
  !-- prune massflux-spikes in top pbl layer --
  do jk=2,klev-2
  do jl=kidia,kfdia
      if ( zmflx(jl,jk-1,jd) < 1.e-40_jprb .and. zmflx(jl,jk,jd) > 1.e-40_jprb ) then
        zmflx(jl,jk,jd) = min( zmflx(jl,jk,jd) , 1.5_jprb*zmflx(jl,jk+1,jd) )
      endif
  enddo
  enddo
  
enddo

! 
! !*         4.9     turbulent orographic drag coefficients 
! 
! call vdftofdc(kidia,kfdia,klon,klev,ztmst,&
!  & pum1,pvm1,pgeom1,psigflt,&
!  & ztofdc)  

ztofdc = 0._jprb


!     ------------------------------------------------------------------

!*         5.     solve advection-diffusion equation
!                 ----------------------------------

!*         5.1  momentum

zmflxm(:,:,:)=zmflx(:,:,:)
zucurr(:)    =0._jprb   ! ocean currents not yet active
zvcurr(:)    =0._jprb   ! ...
zhu1=rvdifts*ztmst
zsoc(kidia:kfdia,1:klev)=psobeta(kidia:kfdia,1:klev)!*zhu1

call vdfdifm (kidia, kfdia, klon , klev  , idraft , itop  , &
            & ztmst, zextmomf, llsfcflx,&
            & pum1 , pvm1 , paphm1, zcfm   , zmflxm, zuuh , zvuh ,&
            & ztofdc,psoteu,psotev,zsoc  , &
            & pvom , pvol , zucurr,zvcurr, zudif  , zvdif , ztaux, ztauy, &
            & llmultisolv )  


!*         5.2  generalized liquid water static energy and total water

call vdfdifh (kidia  , kfdia  , klon   , klev   , idraft , itop   , ktiles, &
            & ztmst  , zextshf, zextlhf, llsfcflx, &
            & pfrti  , pssrflti,pslrfl , pemis  , pevapsnw, &
            & zslgm1 , ptm1   , pqm1   , zqtm1  , paphm1 , &
            & zcfh   , zcfhti , zcfqti , zmflx  , zslguh , zqtuh  , &
            & zslgdif, zqtdif , zcptsti, zqsti  , zcairti, zcsatti, &
            & zdqsti , ptskti , ptskrad, ptsam1m(:,1)    , ptsnow , ptice  , psst, &
            & ztsktip1,zslge  , pte    , zqte   , &
            & pevapti, pahfsti, zahflti, zstr   , zg0, &
            & llmultisolv )  


!*         5.3  incrementation of u and v tendencies, storage of
!*              the dissipation, computation of multilevel fluxes.

call vdfincr (kidia  , kfdia  , klon   , klev   , itop   , ztmst  , &
            & pum1   , pvm1   , zslgm1 , ptm1   , zqtm1  , paphm1 , pgeom1 , &
            & zcfm   , ztofdc , psoteu , psotev , zsoc   ,&
            & zudif  , zvdif  , pucurr , pvcurr , zslgdif, zqtdif , &
            & pvom   , pvol   , zslge  , zqte   , zslgewodis, &
            & pvdis  , pvdisg , pstrtu , pstrtv , pstrsou, pstrsov , ptofdu , ptofdv)  
! 
! !          5.4  solve for tracers
! if (lvdftrac .and. ktrac > 0) then 
!   call vdfdifc(kidia,kfdia,klon,klev,itop,ktrac,&
!              & ztmst,pcm1,ptenc,paphm1,zcfh,pcflx)
! endif

!  do jk=1,klev
!  do jl=kidia,kfdia
!    if (jl==kidia ) then
!      write(0,'(a,i3,a,f15.7,a,f15.7)') 'vdfmain 3.5:  jk=',jk,&
!               & ' zslge:', 3600.*zslge(jl,jk) / rcpd, &
!               & ' zqte:', 3600.*1000._jprb * zqte(jl,jk)
!    endif
!  enddo
!  enddo




!     ------------------------------------------------------------------
!
!*        6.     time integration of conserved state variables qt and slg
!                 --------------------------------------------------------
!
!         compute the conserved state after rad+dyn *and* pbl conv+diff.
!         this will be used later to obtain the tendencies of the non-conserved 
!         prognostic variables t and qv.
!
  
  do jk=1,klev
    do jl=kidia,kfdia
      
      !--wipe out any tendency above pbl--
      if (lltropdamp .and. .not.llpbl(jl,jk)) then
        zqte(jl,jk)  = zqte(jl,jk) - zqtea(jl,jk)
        zslge(jl,jk) = zslge(jl,jk) - zslgea(jl,jk)
      endif
      
      !integrate in time
      zqtupd(jl,jk)  = zqtm1(jl,jk)  + zqte(jl,jk)  * ztmst
      zslgupd(jl,jk) = zslgm1(jl,jk) + zslge(jl,jk) * ztmst
      

      !if (jl==kfdia .and. llpbl(jl,jk)) then
      !  write(0,'(a,i3,a,f15.7,a,f15.7,a,f15.7,a,f15.7,a,f15.7)') 'vdfmain 4.0:  jk=',jk,&
      !         & ' zqtm1=',1000._jprb*zqtm1(jl,jk), &
      !         & ' zslgm1:', zslgm1(jl,jk)/rcpd, &
      !         & ' zqte:', 1000._jprb*zqte(jl,jk) * ztmst, &
      !         & ' zslge:', zslge(jl,jk) * ztmst / rcpd, &
      !         & ' zslgupd:', zslgupd(jl,jk)/ rcpd
      !endif

      !total specific humidity limiter
      zqtupd(jl,jk) = max( 0._jprb, zqtupd(jl,jk))
      zqte(jl,jk)   = (zqtupd(jl,jk) - zqtm1(jl,jk) ) * zrtmst
      
      if (lldiag) then
        !--output of rad+dyn and conv+diff tendencies of qt and thl --
        pextra(jl,jk,71) = zslgea(jl,jk) * 3600._jprb * 24._jprb / rcpd
        pextra(jl,jk,72) = (zslge(jl,jk)-zslgea(jl,jk)) * 3600._jprb * 24._jprb / rcpd

        pextra(jl,jk,73) = zqtea(jl,jk) * 3600._jprb * 24._jprb * 1000._jprb
        pextra(jl,jk,74) = (zqte(jl,jk)-zqtea(jl,jk)) * 3600._jprb * 24._jprb * 1000._jprb
      endif  

    enddo
  enddo



!     ------------------------------------------------------------------

!*         7.     surface fluxes - tiles
!                 ----------------------
!*         and    compute 2m temperature and humidity, 10m wind,
!*                  and gustiness

!  compute wind speed at blending height
! 
! call vdffblend(kidia,kfdia,klon,klev, &
!  & pum1, pvm1, pgeom1, pucurr, pvcurr, zblend, &
!  & zfblend)

! wrap-up computations for the surface and 2t/2d/10u/10v/gustiness computation

! call surfpp( kidia=kidia,kfdia=kfdia,klon=klon,ktiles=ktiles, &
!  & kdhvtls=kdhvtls,kdhftls=kdhftls, &
!  & ptstep=ptsphy, &
! ! input
!  & pfrti=pfrti, pahflti=zahflti, pg0ti=zg0, &
!  & pstrtulev=pstrtu(:,klev), pstrtvlev=pstrtv(:,klev), ptskm1m=ptskm1m, &
!  & pumlev=pum1(:,klev), pvmlev=pvm1(:,klev), pqmlev=pqm1(:,klev), &
!  & pgeomlev=pgeom1(:,klev), pcptspp=zzcpts, pcptgzlev=zcptgz(:,klev), &
!  & paphms=paphm1(:,klev), pz0mw=zz0mw, pz0hw=zz0hw, pz0qw=zz0qw, &
!  & pzdl=zzzdl, pqsapp=zzqsa, pblend=zblend, pfblend=zfblend, pbuom=zzbuom, &
!  & pz0m=pz0m, pevapsnw=pevapsnw,pssrflti=pssrflti, pslrfl=pslrfl, psst=psst, &
!  & pucurr=pucurr, pvcurr=pvcurr, &
! ! updated
!  & pahfsti=pahfsti, pevapti=pevapti, ptske1=ptske1,ptsktip1=ztsktip1, &
! ! output
!  & pdiftslev=pdifts(:,klev), pdiftqlev=pdiftq(:,klev), pustrti=pustrti, &
!  & pvstrti=pvstrti,  ptskti=ptskti, pahflev=pahflev, pahflsb=pahflsb, &
!  & pfwsb=pfwsb, pu10m=pu10m, pv10m=pv10m, pt2m=pt2m, pd2m=pd2m, pq2m=pq2m, &
!  & pgust=pgust, pzidlwv=pzidlwv, &
! ! output ddh
!  & pdhtls=pdhtls &
!  & )
! 


pdiftl  (kidia:kfdia,klev) = 0.0_jprb
pdifti  (kidia:kfdia,klev) = 0.0_jprb
zdiftqt (kidia:kfdia,klev) = 0.0_jprb
zdiftslg(kidia:kfdia,klev) = 0.0_jprb

! zdiftqt (kidia:kfdia,klev) = pdiftq(kidia:kfdia,klev)
! zdiftslg(kidia:kfdia,klev) = pdifts(kidia:kfdia,klev)
DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
    zdiftslg(JL,klev) = zdiftslg(JL,klev) + PFRTI(JL,JT)*PAHFSTI(JL,JT)
    zdiftqt (JL,klev) = zdiftqt (JL,klev) + PFRTI(JL,JT)*PEVAPTI(JL,JT)
  ENDDO
ENDDO


! 
! !     ------------------------------------------------------------------
! 
! !*         8.     slg, qt, u, v flux computations and t,skin tendency
! !                 ---------------------------------------------------
! 
 do jl=kidia,kfdia
   zdiftqt (jl,0) = 0.0_jprb
   pdiftq  (jl,0) = 0.0_jprb
   pdiftl  (jl,0) = 0.0_jprb
   pdifti  (jl,0) = 0.0_jprb
   pdifts  (jl,0) = 0.0_jprb
   zdiftslg(jl,0) = 0.0_jprb
   pstrtu  (jl,0) = 0.0_jprb
   pstrtv  (jl,0) = 0.0_jprb
 enddo
 
 do jk=klev-1,1,-1
   do jl=kidia,kfdia
     zgdph = - (paphm1(jl,jk)-paphm1(jl,jk+1)) * zrg
 !...change in slg,qt,u,v tendencies are converted to fluxes
     zdiftslg(jl,jk) = ( zslgewodis(jl,jk+1) - zslgea(jl,jk+1) ) * zgdph &
                     & + zdiftslg(jl,jk+1)  
     zdiftqt(jl,jk)  = (zqte (jl,jk+1)-zqtea(jl,jk+1))*zgdph + zdiftqt(jl,jk+1)
   enddo
 enddo
 
  
  
  
!     ------------------------------------------------------------------

!*         10.    qt variance budget
!                 ------------------
!
!          flux-gradient production, transport and dissipation
!

  !-- recall variance from previous timestep ----
  do jk=1,klev
    do jl=kidia,kfdia
      !pvar(jl,jk) = pextra(jl,jk,43)**2
      pvar(jl,jk) = 0._jprb
    enddo
  enddo
      
      
      
  do jk=klev-1,2,-1
    do jl=kidia,kfdia
      
      
      !-- recall and update the variance dissipation timescale ztau --
!      ztau(jl,jk)    = pextra(jl,jk,65)
!      if (llpbl(jl,jk)) then
!        ztau(jl,jk) =  pzinv(jl) / max(pwuavg(jl),0.01_jprb)
!        !ztau(jl,jk) =  2._jprb * pzinv(jl) / max(pwuavg(jl),0.01_jprb)
!        
!        !  implicit
!        !ztaunew = 2._jprb * pzinv(jl) / max(pwuavg(jl),0.01_jprb)
!        !ztau(jl,jk) =  ztaunew +  ( ztau(jl,jk) - ztaunew ) * exp( - ztmst / ztau(jl,jk) )
!        
!        ztau(jl,jk) = max(ztau(jl,jk), 100.0_jprb)
!      else
!        ztau(jl,jk) = ztau(jl,jk) + ztmst
!        !ztau(jl,jk) = 2.0e3_jprb
!      endif
!      ztau(jl,jk) = min( 1.0e3_jprb, ztau(jl,jk) )
!      !ztau(jl,jk) = min( 2.0e3_jprb, ztau(jl,jk) )
      
      ztau(jl,jk) =  pzinv(jl) / max(pwuavg(jl),0.01_jprb)

      
      !-- do the individual variance budget terms --
      zvardiss   = 0._jprb
      zvargen    = 0._jprb
      zvartrans  = 0._jprb
      zdqtdz  = 0._jprb
      zwqtf   = 0._jprb
      
      if (llpbl(jl,jk)) then
      
         !--- i  flux-gradient variance production at full level ---
        if (kpbltype(jl)==2 .and. jk<=kvartop(jl)+1 ) then
          !for stratocumulus, protect variance production in top pbl layer against strong capping gradient
          zdqtdz = (zqtupd(jl,jk)-zqtupd(jl,jk+1)) * rg / (pgeom1(jl,jk)-pgeom1(jl,jk+1))
        else
          !otherwise, do it truely centered
          zdqtdz = (zqtupd(jl,jk-1)-zqtupd(jl,jk+1)) * rg / (pgeom1(jl,jk-1)-pgeom1(jl,jk+1))
        endif
          
        !if ( jk==kvartop(jl) ) then
        !  zwqtf = 0._jprb
        !else
          zwqtf = (zdiftqt(jl,jk-1) + zdiftqt(jl,jk))/2._jprb  
        !endif
        zwqtf = -(rd * ptm1(jl,jk) / papm1(jl,jk)) * zwqtf
           
        zvargen = -2._jprb * zwqtf * zdqtdz
        zvargen = max(zvargen,0._jprb)   ! exclude countergradient flow
          
          
        !--- ii   variance transport at full level ---
        zvartrans = - (zwqt2(jl,jk-1) - zwqt2(jl,jk) ) * rg / (pgeoh(jl,jk-1)-pgeoh(jl,jk))


        !--- iii  variance dissipation at full level (for output only) ---
        zvardiss = min(0._jprb, -zvargen -zvartrans)       
           
!      else
!      
!        !--- iii  variance dissipation (for output only) ---
!        zvardiss = min(0._jprb, -pvar(jl,jk)/ztau(jl,jk) )
             
             
       !if ( jk==kvartop(jl)+1 ) then
       !   zvargen = zvargen * 4._jprb
       !endif
       
       
      endif
      
      
      !--- update the variance ---
!      !if (pvar(jl,jk)>0._jprb) then
!      if (pvar(jl,jk)>1.0e-10_jprb) then
!        !--- if variance exists at t-1, do it prognostically (implicit) ---
!        pvar(jl,jk) = (zvargen+zvartrans) * ztau(jl,jk) + &
!                    & (pvar(jl,jk) - (zvargen+zvartrans) * ztau(jl,jk) ) * &
!                    &     exp( - ztmst / ztau(jl,jk) )
!      else
        !--- if no variance exists at t-1, do it diagnostically ---
        pvar(jl,jk) = (zvargen+zvartrans) * ztau(jl,jk) 
!      endif              

      pvar(jl,jk) = max(pvar(jl,jk),0._jprb)
      
      !if (jl==kfdia .and. llpbl(jl,jk)) then
      !  write(0,'(a,i3,a,f10.3,a,f15.7,a,f15.7,a,f15.7)') 'vdfmain:  jk=', jk,  &
      !       & 'pvar='     , 1000000._jprb * pvar(jl,jk), &
      !       & ' zvargen='  , 1000000._jprb * zvargen, &
      !       & ' zvartrans=', 1000000._jprb * zvartrans, &
      !       & ' zqtdt=', 1000._jprb * zdqtdz
      !endif

      if (lldiag) then
        !pextra(jl,jk,57) = 1000000._jprb * (zvargen+zvartrans) * ztau(jl,jk)
        pextra(jl,jk,58) = 1000000._jprb * pvar(jl,jk) 
        pextra(jl,jk,59) = 1000000._jprb * zvartrans 
        pextra(jl,jk,60) = 1000000._jprb * zvargen           
        pextra(jl,jk,61) = 1000000._jprb * zvardiss
        pextra(jl,jk,62) = 1000._jprb * zdqtdz
        pextra(jl,jk,63) = 1000._jprb * zwqtf    
        pextra(jl,jk,64) = 1000000._jprb * zwqt2(jl,jk) 
        pextra(jl,jk,65) = ztau(jl,jk)
      endif  
      
    enddo
  enddo
  
  
  !-- update the prognostic variance and variance timescale --
!!  if (kcnt>0) then
!    do jk=2,klev-1
!    do jl=kidia,kfdia
!!      pextra(jl,jk,65) = ztau(jl,jk)
!      pextra(jl,jk,43) = sqrt( pvar(jl,jk) )
!    enddo
!    enddo
!!  endif



      
!     ------------------------------------------------------------------
!
!*         11.    vectorized bimodal cloud scheme
!                 -------------------------------
!          
!          the edmf decomposition is extended into the cloud scheme, by doing a bimodal
!          pdf: one diffusive, one updraft. each pdf is gaussian, their 1st and 2nd moments are
!          parameterized. the moments of all pdfs are related, see lewellen and yoh (jas, 1993).
!          the scheme is formulated in {thl,qt} space, using vector calculus.
! 
!          variance closures:
!            * the overall variance is done through the full budget (see section 10).
!            * the updraft pdf mean and variance are derived from the resolved size density.
!
!          see neggers (jas, 2009)
!

  !  use new {qt,slg} state, and bulk updraft
  do jk=1,klev
    do jl=kidia,kfdia
      
      if (kpbltype(jl).gt.0 .and. llpbl(jl,jk) .and. zfrac(jl,jk,1)>0._jprb .and. zfrac(jl,jk-1,1)>0._jprb ) then
      
        !-- interpolate the bulk updraft fields to full levels --
        zaup     (jl,jk) = ( zfrac(jl,jk,1)     + zfrac    (jl,jk-1,1) ) / 2._jprb
        zslgup   (jl,jk) = ( zslguh(jl,jk,1)    + zslguh   (jl,jk-1,1) ) / 2._jprb + zslge(jl,jk)
        zqtup    (jl,jk) = ( zqtuh(jl,jk,1)     + zqtuh    (jl,jk-1,1) ) / 2._jprb + zqte(jl,jk)
        zsigqt2up(jl,jk) = ( zqt2uh(jl,jk)      + zqt2uh   (jl,jk-1)   ) / 2._jprb
      
      else
          
        zaup     (jl,jk) = 0._jprb
        zslgup   (jl,jk) = zslgupd(jl,jk)
        zqtup    (jl,jk) = zqtupd(jl,jk)
        zsigqt2up(jl,jk) = 0._jprb
          
      endif
      
    enddo
  enddo
  
  
call vdfcloud ( kidia     , kfdia   , klon    , klev   , idraft , &
                & papm1     , pgeom1  , pgeoh   , &
                & zqtupd    , zslgupd , kvartop , &
                & zaup      , zqtup   , zslgup  , zsigqt2up , &
                & pvar      , &
! diagnostic output
                & pextr2    , kfldx2  , pextra  , klevx  , kfldx  , &
!              
                & zcldfrac  , zqlav   , pldiff)



!     ------------------------------------------------------------------
!
!*         12.    net tendencies
!                 --------------
!
!          calculate net tendencies of the prognostic model variables ql, qi, a, qv and t


  do jk=1,klev
    do jl=kidia,kfdia
      
      !rn --- testing: temporarily switch off cloudiness ---
      !zcldfrac(jl,jk) = 0._jprb
      !zqlav(jl,jk)    = 0._jprb
      !pldiff(jl,jk)   = 0._jprb
      
      !rn --- testing: set cloud fraction and condensate to what's derived from size pdf ---
      !pldiff(jl,jk)   = 0._jprb
      !zcldfrac(jl,jk) = zfracc(jl,jk)
      !zqlav(jl,jk)    = zfracc(jl,jk) * zqcc(jl,jk)
      
      
      !--- convert back from conserved variables qt and slg to non-conserved qv and t ---
      !
      !  new cloud variables (liquid, ice, fraction) for tendency calculation:
      !
      !    within pbl: use new values from vdfcloud
      !    above pbl:  use t-1 input values + tendencies
      !
      !    result: * above pbl, vdfmain does not change cloud variables:
      !                tendencies only contain contributions by rad+dyn.
      !
      if (llpbl(jl,jk)) then    !reset cloudiness within pbl
      !if (.true.) then          !reset cloudiness everywhere! (also above pbl)
      
        !-- safety: total condensate can not be larger than total specific humidity ---
        zqlav(jl,jk) = min( zqtupd(jl,jk) ,zqlav(jl,jk) )
        zaupd(jl,jk) = zcldfrac(jl,jk)
        
        !-- decomposition of total condensate into ice and liquid ---
        zalfaw(jl,jk) = foealfa(ptm1(jl,jk))
!        if (jl==kidia ) then
!          write(0,'(a,i3,a,f15.7,a,f15.7,a,f15.7,a,f15.7,a,f15.7,a,f15.7)') 'vdfmain after vdfcloud:  jk=',jk,&
!               & ' zqtupd:'   , zqtupd(jl,jk), &
!               & ' pqm1:'   , pqm1(jl,jk), &
!               & ' plm1:'   , plm1(jl,jk), &
!               & ' zqlav:'    , zqlav(jl,jk), &
!               & ' zcldfrac:' , zcldfrac(jl,jk), &
!               & ' zalfa:', zalfaw(jl,jk)
!        endif
        zlupd(jl,jk) = zqlav(jl,jk) * zalfaw(jl,jk) 
        ziupd(jl,jk) = zqlav(jl,jk) * ( 1.0_jprb - zalfaw(jl,jk))
        
      else
      
        !-- outside pbl, maintain tendencies from rad+dyn ---
        zaupd(jl,jk) = pam1(jl,jk) + pae(jl,jk)*ztmst
        zlupd(jl,jk) = plm1(jl,jk) + ple(jl,jk)*ztmst
        ziupd(jl,jk) = pim1(jl,jk) + pie(jl,jk)*ztmst
        
      endif  
      
      
      !--- derive non-conserved properties qv and t ---
      zqupd(jl,jk)  = zqtupd(jl,jk) - zlupd(jl,jk) - ziupd(jl,jk)
      ztupd(jl,jk)  = ( zslgupd(jl,jk) - pgeom1(jl,jk) &
        &     + rlvtt * zlupd(jl,jk) + rlstt * ziupd(jl,jk) &
        &   ) / ( rcpd * ( 1.0_jprb + rvtmp2 * zqupd(jl,jk) ) )   !compare to t->slg conversion in section 2.2

      zthlupd(jl,jk) = ( ztupd(jl,jk) - rlvtt * zlupd(jl,jk) / rcpd ) * zexni(jl,jk)    !new thl


      !--- calculate the final tendencies between state at t-1 and state after rad + dyn + pbl ---
      pqe(jl,jk) = ( zqupd(jl,jk) - pqm1(jl,jk) ) * zrtmst
      pte(jl,jk) = ( ztupd(jl,jk) - ptm1(jl,jk) ) * zrtmst
      ple(jl,jk) = ( zlupd(jl,jk) - plm1(jl,jk) ) * zrtmst
      pie(jl,jk) = ( ziupd(jl,jk) - pim1(jl,jk) ) * zrtmst
      pae(jl,jk) = ( zaupd(jl,jk) - pam1(jl,jk) ) * zrtmst

      pthle(jl,jk) = ( zthlupd(jl,jk) - zthlm1(jl,jk) ) * zrtmst   !thl-tendency for ucla-les


      !--- t-check ---
      if ( ztupd(jl,jk).gt.400._jprb.or.ztupd(jl,jk).lt.100._jprb) then
      
        write(0,'(a,3i8)')    'vdfmain t alarm:',kcnt,jl,jk
        write(0,'(a,3i8)')    '               : ', kpbltype(jl),kvartop(jl),khpbln(jl)
        write(0,'(a,2f15.7)')    '        sfluxes: ', zkhfl(jl),zkqfl(jl)
        write(0,'(a,2f15.7)')    '    cld heights: ', zcldbase(jl),zcldtop(jl)
        write(0,'(a,4f15.7)')    '             t: ',&
           & ptm1(jl,jk),ztupd(jl,jk),pte(jl,jk)*ztmst,ztea(jl,jk)*ztmst
        write(0,'(a,5f15.7)')    '           slg: ',&
           & zslgm1(jl,jk)/rcpd       , zslgupd(jl,jk)/rcpd      , &
           & zslge(jl,jk)*ztmst/rcpd  , zslgea(jl,jk)*ztmst/rcpd , &
           & zslgep(jl,jk)*ztmst/rcpd
        write(0,'(a,5f15.7)')    '            qt: ',&
           & zqtm1(jl,jk)*1000._jprb        , zqtupd(jl,jk)*1000._jprb      , &
           & zqte(jl,jk)*ztmst*1000._jprb   , zqtea(jl,jk)*ztmst*1000._jprb , &
           & zqtep(jl,jk)*ztmst*1000._jprb
           
        if (jk<=klev-1) then
          zmgeom  = pgeom1(jl,jk)-pgeom1(jl,jk+1)      
          zcfnc1  = rvdifts * ztmst * rg**2 * paphm1(jl,jk) &
                   & /( zmgeom * rd * 0.5_jprb &
                   & *( ptm1(jl,jk  )*(1.0_jprb+retv*pqm1(jl,jk  )) &
                   & +  ptm1(jl,jk+1)*(1.0_jprb+retv*pqm1(jl,jk+1))))  
          zdummy2 = zmflx(jl,jk,2) / ( zcfnc1 * zmgeom * zrg )
          zdummy3 = zmflx(jl,jk,3) / ( zcfnc1 * zmgeom * zrg )
        else
          zdummy2 = 0.
          zdummy3 = 0.
	endif
	
	write(0,'(a,6f15.7)')    '          updr: ',&
!           & zmflx(jl,jk,2), zslguh(jl,jk,2)/rcpd, zqtuh(jl,jk,2)*1000._jprb, &
!           & zmflx(jl,jk,3), zslguh(jl,jk,3)/rcpd, zqtuh(jl,jk,3)*1000._jprb
           & zdummy2, zslguh(jl,jk,2)/rcpd, zqtuh(jl,jk,2)*1000._jprb, &
           & zdummy3, zslguh(jl,jk,3)/rcpd, zqtuh(jl,jk,3)*1000._jprb

!         do jkk= klev-5,klev
!           write(0,'(a,i6,4f15.7)')    '      m-struct: ',jkk, zmflx(jl,jkk,2), zwuh(jl,jkk,2), zmflx(jl,jkk,3), zwuh(jl,jkk,3)
!         enddo
	
      endif
      
!       
!       !--- q-check ----
!       if (zqtupd(jl,jk).lt.0._jprb.or.zqtupd(jl,jk).gt.0.05_jprb) then
!         write(0,'(3i6,a,4f15.7,a,3i6,2f15.7)') kcnt,jl,jk,'  terror! q=',&
!            & zqtupd(jl,jk),zqtm1(jl,jk),zlupd(jl,jk),plm1(jl,jk),' - ',&
!            & kvartop(jl),khpbln(jl),kpbltype(jl),zcldbase(jl),zcldtop(jl)
!       endif

    enddo
  enddo
  
  
  
!     ----------- output --------------------------------------
!
!      changes required in vdfouter.f90, stat.f90 and ncio.f90

  if (lldiag) then
    do jl=kidia,kfdia

      !Note: 40-44 already used in vdfcloud

      zrhoh(jl,1:klev) = paphm1(jl,1:klev) / ( rd*ptm1(jl,:)*(1.0_jprb+retv*pqm1(jl,:)) )
      zrhoh(jl,0)      = 0._jprb 
      pextra(jl,:,40) = - rlvtt * zrhoh(jl,:) * zdiftqt (jl,:)    !from m/s kg/kg to W/m2, change sign (+ upwards)
      pextra(jl,:,41) = -         zrhoh(jl,:) * zdiftslg(jl,:)    !from m/s J/kg  to W/m2, change sign (+ upwards)

      pextra(jl,1:klev,42) = 100._jprb  * zcldfrac(jl,:)
      pextra(jl,1:klev,43) = 1000._jprb * zqlav(jl,:)

      !do jk=1,klev
      !  if (jl==kfdia) then
      !    write(0,'(a,i3,a,f10.3,a,f10.3)') 'vdfmain:  jk=',jk,&
      !       & ' zcldfrac=', 100._jprb * zcldfrac(jl,jk),' zqlav=', 1000._jprb * zqlav(jl,jk)
      !  endif
      !enddo

    enddo
  endif

      

!     ------------------------------------------------------------------

!*         13.    flux computations
!                 -----------------
!


!*         13.1    q, ql, qi and s flux

do jk=klev-1,1,-1
  do jl=kidia,kfdia
    zgdph = - (paphm1(jl,jk)-paphm1(jl,jk+1)) * zrg
!...changes in q,l,i tendencies are converted to fluxes
    pdiftq(jl,jk) = (pqe (jl,jk+1) - zqea(jl,jk+1)) * zgdph + pdiftq(jl,jk+1)
    pdiftl(jl,jk) = (ple (jl,jk+1) - zlea(jl,jk+1)) * zgdph + pdiftl(jl,jk+1)
    pdifti(jl,jk) = (pie (jl,jk+1) - ziea(jl,jk+1)) * zgdph + pdifti(jl,jk+1)
!...slg=s-lc*ql-ld*qi (same for fluxes)
    pdifts(jl,jk) = zdiftslg(jl,jk) &
                & + rlvtt * pdiftl(jl,jk) + rlstt * pdifti(jl,jk)  
  enddo
enddo

!*         13.2  pbl precipitation enthalpy fluxes

do jk=1,klev
  do jl=kidia,kfdia
    pfhpvl(jl,jk) = -rlvtt*pfplvl(jl,jk)
    pfhpvn(jl,jk) = -rlstt*pfplvn(jl,jk)
  enddo    
enddo    

!write(0,'(a)' ) 'vdfmain: end'


! if (lhook) call dr_hook('vdfmain',1,zhook_handle)
end subroutine vdfmain


end module vdf
