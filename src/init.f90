!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module init

  use grid

  integer, parameter    :: nns = 500
  integer               :: ns
  integer               :: iseed = 0
  integer               :: ipsflg = 1
  integer               :: itsflg = 1
!  integer, dimension(1) :: seed
  real, dimension(nns)  :: us,vs,ts,thds,ps,hs,rts,rss,tks,xs
  real                  :: zrand = 200.
  character  (len=80)   :: hfilin = 'test.'  



contains
  !
  ! ----------------------------------------------------------------------
  ! INITLZ:  this is the main driver for the model's initializ-
  ! ation routines.  it initializes the model according to runtype
  !
  subroutine initialize
!irina use lsvarflg
    use step, only : time, outflg,lsvarflg,case_name
    use stat, only : init_stat
    use mpi_interface, only : appl_abort, myid
    use thrm, only : thermo
!cgils
    use forc, only : lstendflg
    

    implicit none





    if (runtype == 'INITIAL') then
       time=0.
       call arrsnd
       call basic_state
       call fldinit
       dt  = dtlong
    else if (runtype == 'HISTORY') then
       call hstart
    else
       if (myid == 0) print *,'  ABORTING:  Invalid Runtype'
       call appl_abort(0)
    end if
    call sponge_init
    call init_stat(time+dt,filprf,expnm,nzp)
    !
    !irina
       if (lsvarflg) then
       call lsvar_init
       end if
    !cgils
       if (lstendflg) then
       call lstend_init
       end if
    !    
    ! write analysis and history files from restart if appropriate
    ! 
    if (outflg) then
       if (runtype == 'INITIAL') then
          call write_hist(1, time)
          call init_anal(time)
          call thermo(level)
          call write_anal(time)
       else
          call init_anal(time+dt)
          call write_hist(0, time)
       end if
    end if

    return
  end subroutine initialize
  !
  !----------------------------------------------------------------------
  ! FLDINIT: Initializeds 3D fields, mostly from 1D basic state
  ! 
  subroutine fldinit

    use defs, only : alvl, cpr, cp, p00
    use util, only : azero, atob
    use thrm, only : thermo, rslf
    use mpi_interface,  only : nyprocs, nxprocs
    use step, only : case_name

    implicit none

    integer :: i,j,k
    real    :: exner, pres, tk, rc, xran(nzp)
    real    :: xc,yc,zc,dist

    call htint(ns,ts,hs,nzp,th0,zt)

    do j=1,nyp
       do i=1,nxp
          a_ustar(i,j) = 0.
          do k=1,nzp
             a_up(k,i,j)    = u0(k)
             a_vp(k,i,j)    = v0(k)
             a_tp(k,i,j)    = (th0(k)-th00)
             if (associated (a_rp)) a_rp(k,i,j)   = rt0(k)
             a_theta(k,i,j) = th0(k)
             a_pexnr(k,i,j) = 0.
             if (case_name == 'bubble') then
                xc = (nxp-4)*deltax*nxprocs
                yc = (nyp-4)*deltay*nyprocs
                zc = 1400
                dist = sqrt((xt(i)-xc)**2+(yt(i)-yc)**2+(zt(i)-zc)**2)
                a_tp(k,i,j) = a_tp(k,i,j) + max(0.,2.*(1.-dist/1400))
             end if
          end do
       end do
    end do

    if ( allocated (vapor) ) vapor = a_rp

    if ( allocated (liquid) .and. itsflg == 0) then
       do j=1,nyp
          do i=1,nxp
             do k=1,nzp
                exner = (pi0(k)+pi1(k))/cp
                pres  = p00 * (exner)**cpr
                if (itsflg == 0) then
                   tk    = th0(k)*exner
                   rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                   a_tp(k,i,j) = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                   vapor(k,i,j) = a_rp(k,i,j)-rc
                end if
                if (itsflg == 2) then
                   tk    = th0(k)
                   a_theta(k,i,j) = tk/exner
                   rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                   a_tp(k,i,j) = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                   vapor(k,i,j) = a_rp(k,i,j)-rc
                end if
             end do
          end do
       end do
    end if

    k=1
    do while( zt(k+1) <= zrand .and. k < nzp)
       k=k+1
       xran(k) = 0.2*(zrand - zt(k))/zrand
    end do
    call random_pert(nzp,nxp,nyp,zt,a_tp,xran,k) 

    if (associated(a_rp)) then
       k=1
       do while( zt(k+1) <= zrand .and. k < nzp)
          k=k+1
          xran(k) = 5.0e-5*(zrand - zt(k))/zrand
       end do
       call random_pert(nzp,nxp,nyp,zt,a_rp,xran,k) 
    end if

    call azero(nxyzp,a_wp)
    !    
    ! initialize thermodynamic fields
    !
    call thermo (level)

    call atob(nxyzp,a_pexnr,press)

    return
  end subroutine fldinit
  !----------------------------------------------------------------------
  ! SPONGE_INIT: Initializes variables for sponge layer
  ! 
  subroutine sponge_init

    use mpi_interface, only: myid

    implicit none

    integer :: k,kk

    if (nfpt > 0) then
       allocate (spngt(max(1,nfpt)), spngm(max(1,nfpt)))

       do k=nzp-nfpt,nzp-1
          kk = k + 1 - (nzp-nfpt)
          spngt(kk)=max(0.,(zm(nzp)-zt(k))/((zm(nzp)-zm(nzp-nfpt))*distim))
          spngm(kk)=max(0.,(zm(nzp)-zm(k))/((zm(nzp)-zm(nzp-nfpt))*distim))
          spngt(kk) = max(0.,(1./distim - spngt(kk)))
          spngm(kk) = max(0.,(1./distim - spngm(kk)))
       end do

       if(myid == 0) then
          print "(//' ',49('-')/)"
          print '(2X,A17)', 'Sponge Layer Init '
          print '(3X,A12,F9.1,A1)', 'Starting at ', zt(nzp-nfpt), 'm'
          print '(3X,A18,F9.1,A1)', 'Minimum timescale ', 1/spngm(nfpt),'s'
       end if
    end if

    return
  end subroutine sponge_init

  !
  !
  ! ----------------------------------------------------------------------
  ! ARRSND: Arranges the sounding input into proper arrays
  !
  subroutine arrsnd

    use defs, only          : p00,p00i,cp,cpr,rcp,r,g,ep2,alvl,Rm,ep
    use thrm, only          : rslf
    use mpi_interface, only : appl_abort, myid
    use step, only          : case_name

    implicit none

    integer :: k, iterate
    real    :: tavg, zold2, zold1, x1, xx, yy, zz, til
    character (len=245) :: fm0 = &
         "(/,' -------------------------------------------------',/,"       //&
         "'  Sounding Input: ',//,7x,'ps',9x,'hs',7x,'ts',6x ,'thds',6x," // &
         "'us',7x,'vs',7x,'rts',5x,'rel hum',/,6x,'(Pa)',7X,'(m)',6X,'(K)'"// &
         ",6X,'(K)',6X,'(m/s)',4X,'(m/s)',3X,'(kg/kg)',5X,'(%)',/,1x/)"
    character (len=36) :: fm1 = "(f11.1,f10.1,2f9.2,2f9.2,f10.5,f9.1)"
    !
    ! arrange the input sounding
    !
    if (ps(1) == 0.) then
       open (1,file='sound_in',status='old',form='formatted')
       do ns=1,nns
          read (1,*,end=100) ps(ns),ts(ns),rts(ns),us(ns),vs(ns)
       end do
       close (1)
    end if
100 continue

    ns=1
    do while (ps(ns) /= 0. .and. ns <= nns)
       !
       ! filling relative humidity array only accepts sounding in mixing
       ! ratio (g/kg) converts to (kg/kg)
       !
       
       rts(ns)=rts(ns)*1.e-3
       !
       ! filling pressure array:
       ! ipsflg = 0 :pressure in millibars
       ! 1 :pressure array is height in meters (ps(1) is surface pressure)
       !
       select case (ipsflg)
       case (0)
          ps(ns)=ps(ns)*100.
       case default
          xs(ns)=(1.+ep2*rts(ns))
          if (case_name == 'bubble') then
            xs(ns) = rts(ns)*1e3
          end if
          if (ns == 1)then
             ps(ns)=ps(ns)*100.
             zold2=0.
             hs(1) = 0.
          else
             hs(ns) = ps(ns)
             zold1=zold2
             zold2=ps(ns)
             tavg=(ts(ns)*xs(ns)+ts(ns-1)*xs(ns-1)*(p00**rcp)             &
                  /ps(ns-1)**rcp)*.5
             ps(ns)=(ps(ns-1)**rcp-g*(zold2-zold1)*(p00**rcp)/(cp*tavg))**cpr
          end if
       end select
       !
       ! filling temperature array:
       ! itsflg = 0 :potential temperature in kelvin
       !          1 :liquid water potential temperature in kelvin
       !          2 :temperature
       !
       select case (itsflg)
       case (0)
          tks(ns)=ts(ns)*(ps(ns)*p00i)**rcp
       case (1)
          til=ts(ns)*(ps(ns)*p00i)**rcp
          xx=til
          yy=rslf(ps(ns),xx)
          zz=max(rts(ns)-yy,0.)
          if (zz > 0.) then
             do iterate=1,3
                x1=alvl/(cp*xx)
                xx=xx - (xx - til*(1.+x1*zz))/(1. + x1*til                &
                     *(zz/xx+(1.+yy*ep)*yy*alvl/(Rm*xx*xx)))
                yy=rslf(ps(ns),xx)
                zz=max(rts(ns)-yy,0.)
             enddo
          endif
          tks(ns)=xx
       case (2)
          tks(ns) = ts(ns) ! a long way of saying do nothing
       case default
          if (myid == 0) print *, '  ABORTING: itsflg not supported'
          call appl_abort(0)
       end select
       if (case_name == 'bubble') then
         rts(ns) = xs(ns)*rslf(ps(ns),tks(ns))
       end if
       ns = ns+1
    end do
    ns=ns-1


    !                                  
    ! compute height levels of input sounding.
    !
    if (ipsflg == 0) then
       do k=2,ns
          hs(k)=hs(k-1)-r*.5 *(tks(k)*(1.+ep2*rts(k))                      &
               +tks(k-1)*(1.+ep2*rts(k-1)))*(log(ps(k))-log(ps(k-1)))/g
       end do
    end if

    if (hs(ns) < zt(nzp)) then
       if (myid == 0) print *, '  ABORTING: Model top above sounding top'
       if (myid == 0) print '(2F12.2)', hs(ns), zt(nzp)
       call appl_abort(0)
    end if

    do k=1,ns
       thds(k)=tks(k)*(p00/ps(k))**rcp
    end do

    do k=1,ns
       xs(k)=100.*rts(k)/rslf(ps(k),tks(k))
    end do

    if(myid == 0) then
       write(6,fm0)
       write(6,fm1)(ps(k),hs(k),tks(k),thds(k),us(k),vs(k),rts(k),xs(k),k=1,ns)
    endif

604 format('    input sounding needs to go higher ! !', /,                &
         '      sounding top (m) = ',f12.2,'  model top (m) = ',f12.2)
    return
  end subroutine arrsnd
  !
  !----------------------------------------------------------------------
  ! BASIC_STATE: This routine computes the basic state values
  ! of pressure, density, moisture and temperature.  The basi!state 
  ! temperature is assumed to be a the volume weighted average value of
  ! the sounding
  !
  subroutine basic_state

    use defs, only : cp, rcp, cpr, r, g, p00, p00i, ep2
    use mpi_interface, only : myid

    implicit none

    integer k
    real :: v1da(nzp), v1db(nzp), v1dc(nzp), exner

    character (len=305) :: fmt =  &
         "(/,' -------------------------------------------------',/,"     //&
         "'  Basic State: ',//,4X,'Z',6X,'U0',6X,'V0',6X,'DN0',6X,' P0'"   //&
         ",6X,'PRESS',4X,'TH0',6X,'THV',5X,'RT0',/,3X,'(m)',5X,'(m/s)'"     //&
         ",3X,'(m/s)',2X,'(kg/m3)',2X,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X"      //&
         ",'(K)',4X,'(g/kg)',//,(1X,F7.1,2F8.2,F8.3,2F10.2,2F8.2,F7.2))"

    !

    call htint(ns,thds,hs,nzp,th0,zt)
    call htint(ns,us,hs,nzp,u0,zt)
    call htint(ns,vs,hs,nzp,v0,zt)

    if (level >= 1) then
       call htint(ns,rts,hs,nzp,rt0,zt)
       rt0(1)=rt0(2)
    else
       do k=1,nzp
          rt0(k)=0.
       end do
    end if
    !
    ! calculate theta_v for an unsaturated layer, neglecting condensate here is
    ! okay as this is only used for the first estimate of pi1, which will be 
    ! updated in a consistent manner on the first dynamic timestep
    !
    do k=1,nzp
       v1dc(k)=th0(k) * (1.+ep2*rt0(k)) ! theta_v assuming unsaturated
    end do
    !
    ! calculate pressure for actual initial state
    !
    pi1(1)=cp*(ps(1)*p00i)**rcp+g*(hs(1)-zt(1))/v1dc(1)
    do k=2,nzp
       pi1(k) = pi1(k-1)-g/(dzm(k-1)*0.5*(v1dc(k)+v1dc(k-1)))
    end do
    !
    ! calculate hydrostatic exner function associated with th00 constant along
    ! with associated basic state density
    !
    pi0(1)=cp*(ps(1)*p00i)**rcp + g*(hs(1)-zt(1))/th00
    dn0(1)=((cp**(1.-cpr))*p00)/(r*th00*pi0(1)**(1.-cpr))
    do k=2,nzp
       pi0(k)=pi0(1) + g*(zt(1) - zt(k))/th00
       dn0(k)=((cp**(1.-cpr))*p00)/(r*th00*pi0(k)**(1.-cpr))
       u0(k)=u0(k)-umean
       v0(k)=v0(k)-vmean
    end do
    !
    ! define pi1 as the difference between pi associated with th0 and pi 
    ! associated with th00, thus satisfying pi1+pi0 = pi = cp*(p/p00)**(R/cp)
    !
    do k=1,nzp
       pi1(k) = pi1(k)-pi0(k)
    end do
    !
    do k=1,nzp
       exner = (pi0(k) + pi1(k))/cp
       v1db(k)=p00*(exner)**cpr      ! pressure
       v1da(k)=p00*(pi0(k)/cp)**cpr  ! pressure associated with pi0
    end do

    u0(1) = u0(2)
    v0(1) = v0(2)
    psrf  = ps(1)

    if(myid == 0) write (*,fmt) (zt(k),u0(k),v0(k),dn0(k),v1da(k),v1db(k), &
         th0(k),v1dc(k),rt0(k)*1000.,k=1,nzp)

    return
  end subroutine basic_state
  !
  !---------------------------------------------------------------------
  ! HTINT: Height interpolation of field on one grid, to field on another
  !
  subroutine htint(na,xa,za,nb,xb,zb)

    implicit none
    integer, intent (in) :: na, nb
    real, intent (in)    :: xa(na),za(na),zb(nb)
    real, intent (out)   :: xb(nb)

    integer :: l, k
    real    :: wt

    l = 1
    do k=1,nb
       if (zb(k) <= za(na)) then
          do while ( zb(k) > za(l+1) .and. l < na)
             l=l+1
          end do
          wt=(zb(k)-za(l))/(za(l+1)-za(l))
          xb(k)=xa(l)+(xa(l+1)-xa(l))*wt    
       else
          wt=(zb(k)-za(na))/(za(na-1)-za(na))
          xb(k)=xa(na)+(xa(na-1)-xa(na))*wt
       end if
    end do

    return
  end subroutine htint
  !
  !----------------------------------------------------------------------
  ! HSTART:  This subroutine reads a history file and does
  ! a history start
  !
  subroutine hstart

    use step, only : time
    use mpi_interface, only : myid

    implicit none

    call read_hist(time, hfilin)

    if(myid == 0) &
         print "(//' ',49('-')/,' ',/,' History read from: ',A60)",hfilin

    return
  end subroutine hstart
  !
  !----------------------------------------------------------------------
  ! RANDOM_PERT: initialize field between k=2 and kmx with a
  ! random perturbation of specified magnitude
  !
  subroutine random_pert(n1,n2,n3,zt,fld,xmag,kmx)

    use mpi_interface, only :  nypg,nxpg,myid,wrxid,wryid,xoffset,yoffset, &
         double_scalar_par_sum

    use util, only : sclrset
    implicit none
    integer, allocatable :: seed(:)

    integer, intent(in) :: n1,n2,n3,kmx
    real, intent(inout) :: fld(n1,n2,n3)
    real, intent(in)    :: zt(n1),xmag(n1)

    real (kind=8) :: rand(3:n2-2,3:n3-2),  xx, xxl
    real (kind=8), allocatable :: rand_temp(:,:)
    integer :: i,j,k,n,n2g,n3g

    rand=0.0

    call random_seed(size=n)
    allocate (seed(n))
    seed = iseed * (/ (i, i = 1, n) /)
    call random_seed(put=seed)
    deallocate (seed)
    ! seed must be a double precision odd whole number greater than 
    ! or equal to 1.0 and less than 2**48.
    !seed(1) = iseed        
    !call random_seed(put=seed)
    n2g = nxpg
    n3g = nypg

    do k=2,kmx
       allocate (rand_temp(3:n2g-2,3:n3g-2))
       call random_number(rand_temp)
       rand(3:n2-2, 3:n3-2)=rand_temp(3+xoffset(wrxid):n2+xoffset(wrxid)-2, &
            3+yoffset(wryid):n3+yoffset(wryid)-2)
       deallocate (rand_temp)

       xx = 0.
       do j=3,n3-2
          do i=3,n2-2
             fld(k,i,j) = fld(k,i,j) + rand(i,j)*xmag(k)
          end do
       end do

       xxl = xx
       call double_scalar_par_sum(xxl,xx)
       xx = xx/real((n2g-4)*(n3g-4))
       fld(k,:,:)= fld(k,:,:) - xx
    end do

    if(myid == 0) then
       print *
       print *,'-------------------------------------------------'
       print 600,zt(kmx),rand(3,3),xx
       print *,'-------------------------------------------------'
    endif

    call sclrset('cnst',n1,n2,n3,fld)

    return

600 format(2x,'Inserting random temperature perturbations',    &
         /3x,'Below: ',F7.2,' meters;',                        &
         /3x,'with test value of: ',E12.5,                     &
         /3x,'and a magnitude of: ',E12.5)
  end subroutine random_pert

!irina
  !----------------------------------------------------------------------
  ! Lsvar_init if lsvarflg is true reads the lsvar forcing from the respective
  ! file lscale_in
  ! 
  subroutine lsvar_init

   use forc,only   : t_ls,div_ls,sst_ls,ugeo_ls,vgeo_ls

    implicit none
    
    ! reads the time varying lscale forcings
    !
    if (t_ls(2) == 0.) then
       open (1,file='lscale_in',status='old',form='formatted')
         ! print *, 'lsvar_init read'                 
       do ns=1,nns
          read (1,*,end=100) t_ls(ns),div_ls(ns),sst_ls(ns),&
                             ugeo_ls(ns),vgeo_ls(ns)
          !print *, t_ls(ns),div_ls(ns),sst_ls(ns), ugeo_ls(ns),vgeo_ls(ns)                 
       end do
       close (1)
    end if
100 continue
 
    return
  end subroutine lsvar_init

!cgils
  !----------------------------------------------------------------------
  ! Lstend_init if lstendflg is true reads the lstend  from the respective
  ! file lstend_in
  ! 
  subroutine lstend_init

   use grid,only   : wfls,dqtdtls,dthldtls
    use mpi_interface, only : myid
   

   implicit none
   
   real     :: tmp1

    ! reads the time varying lscale forcings
    !
 !   print *, wfls
    if (wfls(2) == 0.) then
 !         print *, 'lstend_init '                 
       open (1,file='lstend_in',status='old',form='formatted')
        if(myid == 0)  print *, 'lstend_init read'                 
       do ns=1,nzp-1
          read (1,'(F10.3,2X, E10.3,2X, E10.3,2X, E10.3)',end=100) tmp1,wfls(ns),dqtdtls(ns),dthldtls(ns)
        if(myid == 0)  print *, ns, tmp1,wfls(ns),dqtdtls(ns),dthldtls(ns)
       end do
       close (1)
    end if
100 continue
 
    return
  end subroutine lstend_init


  !


end module init
