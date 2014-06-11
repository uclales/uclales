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
  use ncio

  implicit none

  integer, parameter    :: nns = 500
  integer               :: ns
  integer               :: iseed = 0
  integer               :: ipsflg = 1
  integer               :: itsflg = 1
  integer               :: irsflg = 1
!  integer, dimension(1) :: seed
  real, dimension(nns)  :: us,vs,ts,thds,ps,hs,rts,rss,tks,xs,xsi,thl
  real                  :: zrand = 200.
  character  (len=80)   :: hfilin = 'test.'
  logical               :: lhomrestart = .false.
  real                  :: mag_pert_q = 5.0e-5
  real                  :: mag_pert_t = 0.2

contains
  !
  ! ----------------------------------------------------------------------
  ! INITLZ:  this is the main driver for the model's initializ-
  ! ation routines.  it initializes the model according to runtype
  !
  subroutine initialize
!irina use lsvarflg
    use step, only : time, outflg,lsvarflg, &
! LINDA, b
    lanom
!LINDA, e
    use stat, only : init_stat, write_ps, statistics
    use grid, only : isfctyp,level
    use mpi_interface, only : appl_abort, myid
    use thrm, only : thermo
!cgils
    use forc, only : lstendflg

    use mcrp, only : initmcrp
    use modcross, only : initcross, triggercross
    use grid, only : nzp, dn0, u0, v0, zm, zt, isfctyp
    use modparticles, only: init_particles, lpartic, lpartdump, lpartstat, initparticledump, initparticlestat, write_particle_hist, particlestat
    use lsmdata, only : initlsm_simple, initlsm
    use step, only : sst

    implicit none

    real ::   &
         t_ano(nzp,nxp-4,nyp-4), &
         q_ano(nzp,nxp-4,nyp-4), &
         u_ano(nzp,nxp-4,nyp-4), &
         v_ano(nzp,nxp-4,nyp-4), &
         w_ano(nzp,nxp-4,nyp-4)   
    integer :: k, i, j

    if (runtype == 'INITIAL') then
       time = 0.
       call random_init
       call arrsnd
       call basic_state
       call fldinit

       if (lanom) then
          call larm_init_anom (nzp,nxp-4,nyp-4,t_ano,q_ano,u_ano,v_ano,w_ano)
          do k=1,nzp
            do i=3,nxp-2
              do j=3,nyp-2
                a_tp(k,i,j) = a_tp(k,i,j) + t_ano(k,i-2,j-2)
                if(level>0) a_rp(k,i,j) = max(0.0,a_rp(k,i,j) + q_ano(k,i-2,j-2))
                a_up(k,i,j) = a_up(k,i,j) + u_ano(k,i-2,j-2)
                a_vp(k,i,j) = a_vp(k,i,j) + v_ano(k,i-2,j-2)
                a_wp(k,i,j) = a_wp(k,i,j) + w_ano(k,i-2,j-2)
              end do
            end do
          end do
       end if
       dt  = dtlong
    else if (runtype == 'HISTORY') then
       call hstart
       if (lhomrestart) then
          call homogenize
       end if
    else
       if (myid == 0) print *,'  ABORTING:  Invalid Runtype'
       call appl_abort(0)
    end if
    call sponge_init

    call initmcrp(level)

    if(isfctyp==5)  call initlsm(sst,time)
    if(isfctyp==55) call initlsm_simple

    call init_stat(time+dt,filprf,expnme,nzp)
     
    if (lsvarflg) then
       call lsvar_init
    end if

    if (lstendflg) then
       call lstend_init
    end if

    if (lpartic) then
      if(runtype == 'INITIAL') then
        call init_particles(.false.)
      else
        call init_particles(.true.,hfilin)
      end if
      if(lpartdump) call initparticledump(time)
      if(lpartstat) call initparticlestat(time)
    end if

    ! write analysis and history files from restart if appropriate
    !
    if (outflg) then
       if (runtype == 'INITIAL') then
          call write_hist(1, time)
          if(lpartic) call write_particle_hist(1,time)
          call init_anal(time)
          call thermo(level)
          call write_anal(time)
          call initcross(time, filprf)
          call triggercross(time)
          call statistics (time)
          call write_ps(nzp,dn0,u0,v0,zm,zt,time)
          if(lpartic) call particlestat(.false.,time)
          if(lpartic) call particlestat(.true.,time)
       else
          call init_anal(time+dt)
          call initcross(time, filprf)
          call thermo(level)
          call write_hist(0, time)
          if(lpartic) call write_particle_hist(0,time)
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
    use step, only : case_name, lanom
    use mpi_interface, only : myid

    implicit none

    integer :: i,j,k
    real    :: exner, pres, tk, rc, xran(nzp), zc, dist, xc
    real, dimension(nzp)  :: thli
! LINDA, b
    real    :: qv, rh
    real, allocatable :: f_xyz_3d(:,:,:)
    allocate(f_xyz_3d(nzp,nxp,nyp))
    f_xyz_3d = 0.0
! LINDA, e

    call htint(ns,ts,hs,nzp,th0,zt)
    call htint(ns,thl,hs,nzp,thli,zt)

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
          end do
       end do
    end do

    if ( allocated (vapor) ) vapor = a_rp

    if ( allocated (liquid)) then
       do j=1,nyp
          do i=1,nxp
             do k=1,nzp
                exner = (pi0(k)+pi1(k))/cp
                pres  = p00 * (exner)**cpr
                tk    = th0(k)*exner
                rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                a_tp(k,i,j) = thli(k) - th00
                vapor(k,i,j) = a_rp(k,i,j)-rc
             end do
          end do
       end do
    end if
    if (case_name == 'bubble') then
      xc = 1e4
      zc = 1400
      do j=1,nyp
        do i=1,nxp
          do k=1,nzp
            if (zt(k)>0 .and. zt(k)< 2 * zc) then
              dist = (xt(i)**2 + yt(j)**2)/xc**2 + (zt(k) - zc)**2/zc**2
              a_tp(k,i,j) = a_tp(k,i,j) + 2. *max(0.,(1-dist))
            end if
          end do
        end do
      end do
! LINDA, b
    elseif (case_name == 'squall') then

      CALL  squall3d_Morrison(f_xyz_3d)

      do j=1,nyp
        do i=1,nxp
          do k=1,nzp
             exner = (pi0(k)+pi1(k))/cp
             pres  = p00 * (exner)**cpr
             tk    = a_theta(k,i,j)*exner
             rc    = max(0.,a_rp(k,i,j)-rslf(pres,tk))
             rh    = vapor(k,i,j)/rslf(pres,tk)

             ! add perturbation to theta
             a_theta(k,i,j) = a_theta(k,i,j) + f_xyz_3d(k,i,j)

             tk    = a_theta(k,i,j)*exner
             qv    = rh*rslf(pres,tk)
             a_tp(k,i,j)  = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
             vapor(k,i,j) = qv
             a_rp(k,i,j)  = qv+rc

          end do
        end do
      end do

      WRITE (*,*) '=== SQUALL3D-Testcase: dTmax = ', &
           MINVAL( f_xyz_3d(:,:,:) ),'myid=',myid

    end if
! LINDA, e

    k=1
    do while( zt(k+1) <= zrand .and. k < nzp)
       k=k+1
       xran(k) = mag_pert_t*(zrand - zt(k))/zrand
    end do
! LINDA, b
    if (.not.lanom) call random_pert(nzp,nxp,nyp,zt,a_tp,xran,k) 
! LINDA, e

    if (associated(a_rp)) then
       k=1
       do while( zt(k+1) <= zrand .and. k < nzp)
          k=k+1
          xran(k) = mag_pert_q*(zrand - zt(k))/zrand
       end do
! LINDA, b
    if (.not.lanom) call random_pert(nzp,nxp,nyp,zt,a_rp,xran,k) 
! LINDA, e
    end if
    call azero(nxyzp,a_wp)
    !
    ! initialize thermodynamic fields
    !
    call thermo (level)
    call atob(nxyzp,a_pexnr,press)
    deallocate (f_xyz_3d)
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

          print*,spngt(kk),spngm(kk)
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
    use thrm, only          : rslf,rsif
    use mpi_interface, only : appl_abort, myid

    implicit none

    integer :: k, iterate, iterate1
    real    :: tavg, zold2, zold1, x1, xx, yy, zz, til
    character (len=260) :: fm0 = &
         "(/,' -------------------------------------------------',/,"       //&
         "'  Sounding Input: ',//,7x,'ps',9x,'hs',7x,'ts',6x ,'thds',6x," // &
         "'us',7x,'vs',7x,'rts',5x,'rel hum',5x,'rhi'/,6x,'(Pa)',7X,'(m)',6X,'(K)'"// &
         ",6X,'(K)',6X,'(m/s)',4X,'(m/s)',3X,'(kg/kg)',5X,'(%)',5X,'(%)'/,1x/)"
    character (len=37) :: fm1 = "(f11.2,f10.2,2f9.2,2f9.2,f10.5,2f9.1)"
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
       ! irsflg = 1:
       ! filling relative humidity array only accepts sounding in mixing
       ! ratio (g/kg) converts to (kg/kg)
       ! irsflg = 0:
       ! use relative humidity sounding
          if (irsflg == 0) then
            xs(ns)  = rts(ns)
            if (ns > 1) then
              rts(ns) = xs(ns)*rslf(ps(ns-1),tks(ns-1))  !cheat a tiny bit - fix it later.
            else
              rts(ns)  = xs(ns)*rslf(ps(ns)*100.,ts(ns))
            end if
          else
            rts(ns) = rts(ns)*1.e-3
          end if


       !
       ! filling pressure array:
       ! ipsflg = 0 :pressure in millibars
       ! 1 :pressure array is height in meters (ps(1) is surface pressure)
       !
       select case (ipsflg)
       case (0)
          ps(ns)=ps(ns)*100.
       case default
          if (ns == 1)then
             ps(ns)=ps(ns)*100.
             zold2=0.
             hs(1) = 0.
          else
             hs(ns) = ps(ns)
             zold1=zold2
             zold2=ps(ns)
             if ((itsflg==0)) then
                tavg=(ts(ns)*(1.+ep2*rts(ns))+ts(ns-1)*(1.+ep2*rts(ns-1))*(p00**rcp)             &
                  /ps(ns-1)**rcp)*.5
                ps(ns)=(ps(ns-1)**rcp-g*(zold2-zold1)*(p00**rcp)/(cp*tavg))**cpr
             elseif ((itsflg==1)) then
                tavg=(ts(ns)*(1.+ep2*rts(ns))+ts(ns-1)*(1.+ep2*rts(ns-1))*(p00**rcp)             &
                  /ps(ns-1)**rcp)*.5
                ps(ns)=(ps(ns-1)**rcp-g*(zold2-zold1)*(p00**rcp)/(cp*tavg))**cpr
             else !itsflg=2
                tavg=(ts(ns)*(1.+ep2*rts(ns))+tks(ns-1)*(1.+ep2*rts(ns-1)))*.5
                ps(ns)=ps(ns-1)*exp(-g*(zold2-zold1)/(r*tavg))
             end if
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
          thl(ns) = ts(ns)
          til=ts(ns)*(ps(ns)*p00i)**rcp
          xx=til
          yy=rslf(ps(ns),xx)
          zz=max(rts(ns)-yy,0.)
          if (zz > 0.) then
             do iterate=1,1
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

       ! at this point rts, ps and tks are approximately known for all cases

       do iterate1 = 1,5
         if (irsflg == 0) then
            rts(ns) = xs(ns)*rslf(ps(ns),tks(ns))
         end if
         if (ipsflg == 1 .and. ns > 1) then
            tavg=(tks(ns)*(1.+ep2*rts(ns))+tks(ns-1)*(1.+ep2*rts(ns-1)))*.5
            ps(ns)=ps(ns-1)*exp(-g*(hs(ns)-hs(ns-1))/(r*tavg))
         end if
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
            ts(ns)=tks(ns)*(p00/ps(ns))**rcp     ! update ts for fldinit
         case default
         end select
         
       end do
       if (itsflg==1) then
          ts(ns)=tks(ns)*(p00/ps(ns))**rcp     ! update ts for fldinit
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

    ! check if model top is below sounding top.

    if (hs(ns) < zt(nzp)) then
       if (myid == 0) print *, '  ABORTING: Model top above sounding top'
       if (myid == 0) print '(2F12.2)', hs(ns), zt(nzp)
       call appl_abort(0)
    end if

    do k=1,ns
       thds(k)=tks(k)*(p00/ps(k))**rcp ! thds goes into the basic state
    end do

    do k=1,ns
       xs(k)=100.*rts(k)/rslf(ps(k),tks(k))
    end do

    do k=1,ns
       xsi(k)=100.*rts(k)/rsif(ps(k),tks(k))
    end do
    ! calculate thl for fldinit for itsflg=0,2
    if ((itsflg==0).or.(itsflg==2)) then
       do k=1,ns
            yy=rslf(ps(k),tks(k))
            zz=max(rts(k)-yy,0.)
            if (zz > 0.) then
               thl(k)= ts(k)*exp(-(alvl/cp)*zz/tks(k))
            else
               thl(k)= ts(k)
            end if
         end do
      end if

    if(myid == 0) then
       write(6,fm0)
       write(6,fm1)(ps(k),hs(k),tks(k),thds(k),us(k),vs(k),rts(k),xs(k),xsi(k),k=1,ns)
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
    use thrm, only :rslf

    implicit none

    integer k
    real :: v1da(nzp), v1db(nzp), v1dc(nzp),x0(nzp), exner

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
       pi1(k) = pi1(k-1)-g/(dzi_m(k-1)*0.5*(v1dc(k)+v1dc(k-1)))
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
       if(irsflg==0) then
         call htint(ns,xs,hs,nzp,x0,zt)
         rt0(k) = 1e-2*x0(k)*rslf(v1db(k),exner*th0(k))
       end if
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

    integer, intent(in) :: n1,n2,n3,kmx
    real, intent(inout) :: fld(n1,n2,n3)
    real, intent(in)    :: zt(n1),xmag(n1)

    real (kind=8) :: rand(3:n2-2,3:n3-2),  xx, xxl
    real (kind=8), allocatable :: rand_temp(:,:)
    integer :: i,j,k,n2g,n3g

    rand=0.0

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

  subroutine random_init
    use mpi_interface, only: myid
    integer :: i, n
    integer, allocatable, dimension(:) :: seed
    call random_seed(size=n)
    allocate (seed(n))
    seed = iseed * (/ (i, i = 1, n) /) + myid
    call random_seed(put=seed)
    deallocate (seed)
  end subroutine random_init
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

   real     :: lowdthldtls,highdthldtls,lowdqtdtls,highdqtdtls,lowwfls,highwfls,highheight,lowheight,fac
   integer :: k

    ! reads the time varying lscale forcings
    !
 !   print *, wfls
    if (wfls(2) == 0.) then
 !         print *, 'lstend_init '
        open (1,file='lstend_in',status='old',form='formatted')
        read (1,*,end=100) lowheight,lowwfls,lowdqtdtls,lowdthldtls
        read (1,*,end=100) highheight,highwfls,highdqtdtls,highdthldtls
        if(myid == 0)  print *, 'lstend_init read'
        do  k=2,nzp-1
          if (highheight<zt(k)) then
            lowheight = highheight
            lowwfls = highwfls
            lowdqtdtls = highdqtdtls
            lowdthldtls = highdthldtls
            read (1,*) highheight,highwfls,highdqtdtls,highdthldtls
          end if
          fac = (highheight-zt(k))/(highheight - lowheight)
          wfls(k) = fac*lowwfls + (1-fac)*highwfls
          dqtdtls(k) = fac*lowdqtdtls + (1-fac)*highdqtdtls
          dthldtls(k) = fac*lowdthldtls + (1-fac)*highdthldtls
        end do
       close (1)
    end if
100 continue

    return
  end subroutine lstend_init

  subroutine homogenize
    use util, only : get_avg3,azero, atob
    use thrm, only : thermo
    use grid, only : a_ustar, a_tstar, a_rstar
    implicit none

    integer :: i,j,k,n
    real    :: exner, pres, tk, rc, xran(nzp), zc, dist, xc
    real :: ust(3,nxp,nyp), ustmn(3), prof(nzp)
    ust(1,:,:) = a_ustar
    ust(2,:,:) = a_tstar
    ust(3,:,:) = a_rstar
    call get_avg3(3,nxp,nyp,ust,ustmn)
    a_ustar = ustmn(1)
    a_tstar = ustmn(2)
    a_rstar = ustmn(3)
    do n = 1,nscl
      call get_avg3(nzp,nxp,nyp,a_xp(:,:,:,n),prof)
      do k = 1,nzp
        a_xp(k,:,:,n) = prof(k)
      end do
    end do

    k=1
    do while( zt(k+1) <= zrand .and. k < nzp)
       k=k+1
       xran(k) = 0.02*(zrand - zt(k))/zrand
       !xran(k) = 0.05*(zrand - zt(k))/zrand
    end do
    call random_pert(nzp,nxp,nyp,zt,a_tp,xran,k)

    if (associated(a_rp)) then
       k=1
       do while( zt(k+1) <= zrand .and. k < nzp)
          k=k+1
          xran(k) = 1.0e-5*(zrand - zt(k))/zrand
          !xran(k) = 1.0e-5*(zrand - zt(k))/zrand
       end do
       call random_pert(nzp,nxp,nyp,zt,a_rp,xran,k)
    end if
    call azero(nxyzp,a_wp)
    !
    ! initialize thermodynamic fields
    !
    call thermo (level)
    call atob(nxyzp,a_pexnr,press)

  end subroutine homogenize 

!--------------------------------------------------------------------------!
! routine to start the model with noise resulting from previous simulation !
! use cdo script cdo_anomaly in misc/scripts/                              !
!--------------------------------------------------------------------------!
 subroutine larm_init_anom (n1,n2,n3,t_ano,q_ano,u_ano,v_ano,w_ano)

   use netcdf
   use mpi_interface, only:myid,pecount, wrxid, wryid

   implicit none

   integer, intent(in) :: n1,n2,n3
   integer :: nx,ny,nz
   integer :: k, i, j
   integer             ::ncid,status
   real, intent(inout) ::   &
         t_ano(n1,n2,n3), &
         q_ano(n1,n2,n3), &
         u_ano(n1,n2,n3), &
         v_ano(n1,n2,n3), &
         w_ano(n1,n2,n3)   

   character (len=88) :: lfname
   character (len=80) :: fname
   integer            :: varid, dimid

    fname =  trim(filprf)
    if (pecount > 1) then
       write(lfname,'(a,a6,i4.4,i4.4,a3)') trim(fname),'.anom.',wrxid,wryid,'.nc'
    else
       write(lfname,'(a,a8)') trim(fname),'.anom.nc'
    end if
    print*, 'opening file: ', lfname

!*  Open
    status=nf90_open(lfname,nf90_nowrite,ncid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    if (myid==0) print*,'opened netcdf file'

    status = nf90_inq_dimid(ncid, "xt", DimID)
    status=nf90_inquire_dimension(ncid,dimid,len=nx)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status = nf90_inq_dimid(ncid, "yt", DimID)
     status=nf90_inquire_dimension(ncid,dimid,len=ny)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status = nf90_inq_dimid(ncid, "zt", DimID)
    status=nf90_inquire_dimension(ncid,dimid,len=nz)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    print*,'nx=',nx,'ny=',ny,'nz=',nz
!* Read
    status=nf90_inq_varid(ncid,"t",varid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status=nf90_get_var(ncid,varid,t_ano)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    if (myid == 0) print*,'read in t'
    status=nf90_inq_varid(ncid,"q",varid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status=nf90_get_var(ncid,varid,q_ano)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    if (myid == 0) print*,'read in q'
    status=nf90_inq_varid(ncid,"u",varid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status=nf90_get_var(ncid,varid,u_ano)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    if (myid == 0) print*,'read in u'
    status=nf90_inq_varid(ncid,"v",varid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status=nf90_get_var(ncid,varid,v_ano)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    if (myid == 0) print*,'read in v'
    status=nf90_inq_varid(ncid,"w",varid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    status=nf90_get_var(ncid,varid,w_ano)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)
    if (myid == 0) print*,'read in w'
!* Close
    status=nf90_close(ncid)
    if (status.ne.nf90_noerr) print*,nf90_strerror(status)   

 end subroutine larm_init_anom

! linda, e
  !----------------------------------------------------------------------------
  ! Shape function for the "Squall3D"-testcase temperature disturbance
  ! of Hugh Morrison (NCAR) for the 8th WMO Cloud Modeling Workshop case #2
  ! including updates from September, 2012
  ! Linda Schlemmer, September 2012
  !----------------------------------------------------------------------------

  subroutine  squall3d_Morrison(f_xyz)

    use mpi_interface, only: myid,nxpg,nypg
    use defs, only : pi, cp, rcp, cpr, r, g, p00, p00i, ep2

    implicit none

    ! spatial shape function of the bubble:
    real,       intent(out) :: f_xyz(nzp,nxp,nyp)

    integer      :: i, j, k, n

    ! variables for squallline initialization:
    real   :: &
         v_mag               ,&! max magnitude (K) of random variation
         rndm_nmbr(10000000)       ! just some random numbers

    real, allocatable :: noisedummy(:)

    real::bub_radx,bub_radz,bub_centi,bub_centj,bub_centz


    f_xyz(:,:,:) = 0.0
    v_mag             =      0.04                  ! max. rel. magnitude of random variation.

    CALL seed_random_number( 404 +myid)
    ! on some compilers, the random series shows problems within the first few
    ! hundred random numbers (the numbers are not really random, but
    ! can be monotonic and reproducibly the same on all processors).
    ! Only after some numbers, the series gets more random type.
    ! Therefore, fetch 5000 dummy random numbers, before doing the
    ! really needed random numbers:
    ALLOCATE(noisedummy(5000))
    CALL RANDOM_NUMBER(noisedummy)
    DEALLOCATE(noisedummy)
      
    CALL RANDOM_NUMBER(rndm_nmbr)
    
    ! in order to initialize the squall line a cold pool with a horizontal extent of 200km 
    ! with a maximum temperature disturbance of -5K at the surface, decreasing linearly 
    ! to 0 at a height of 4.5 km is applied to theta.
    ! since uclales does not have open boundaries and a sponge is used instead, I decided to
    ! add 100km to the domain in the zonal direction. The cold pool starts at x=100km and stops
    ! at x=300km. Between x=275km and x=300km random temperature perturbations are added up to
    ! a height of 4.5 km

    n=0
    DO k = 1, nzp
       DO i = 2, nxp-1
          ! distance from the left boundary of the domain
          IF (( xt(i)>-12500.).and.(xt(i)<12500.).and.(zt(k).le.4500.)) THEN
             DO j = 2, nyp-1
                n=n+1
                f_xyz(k,i,j) = v_mag*(rndm_nmbr(n)-0.5)
             ENDDO
          ELSE
             f_xyz(k,i,:) = 0.
          END IF
       ENDDO
    ENDDO
    
    IF (myid == 0) THEN
      WRITE (*,*) '=== Morrison SQUALL3D-Testcase:'
    ENDIF
    
  END SUBROUTINE squall3d_Morrison
  !==============================================================================
  !==============================================================================
  !
  ! Initialisation of the random number generator on a parallel machine
  ! with the following properties:
  !
  ! - if initialized with no "iseed_in" (optional integer parameter), then
  !   the resulting random number series (RNS) will be different for each
  !   model run, because seed is determined from the microseconds part of the actual date_and_time().
  !   Additionally, the PE-number is blended into the seed to ensure
  !   that the RNS is also different on each processor, even if called at the exact
  !   same time.
  !
  ! - if "iseed_in" is provided, this is used to generate the same seed on 
  !   a processor with a given PE-number each time the program runs, but
  !   again different seeds on PEs with different PE-numbers.
  !   This enables parallel generation of different random number series
  !   on each processor, which are however the same at each successive program run on
  !   each corresponding PE.
  !
  ! - NOTE: If you would like to, e.g., impose a random but reproducible noise on a 
  !   model field (i.e., INDEPENDENT of the number of PEs) which stays the same for each
  !   successive model run, then it is proposed that you calculate the noisy field
  !   globally on one processor and distribute it afterwards to the single nodes using
  !   the SR distribute_field() from parallel_utilities.f90.
  !   THIS IS DONE IN SR gen_bubnoise() BELOW !
  !
  !==============================================================================
  !==============================================================================

  !.. PGI-friendly version:
  subroutine seed_random_number(iseed_in) 
    use mpi_interface, only: myid,nxprocs,nyprocs
    
    implicit none 
    
    !.. local vars
    integer, optional, intent(in) :: iseed_in

! LS2011b, for cscs, we do need the detailed type specifications, pgi will terminate
! with floating point exceptions.

!!! UB>> Older settings with detailed type specifications seem not
!!!      to be necessary any more (aside from pgi-compiler, which
!!!      cannot be tested at DWD!
! LS, for cscs, we do need the detailed type specifications, pgi will terminate
! with floating point exceptions.
!    INTEGER*4 :: i
!    INTEGER*4 :: k
!    INTEGER*4 :: i1
!    INTEGER :: zeit(8), i2, iseed,num_compute
!    INTEGER*4, ALLOCATABLE :: seed(:)


    INTEGER :: i
    INTEGER :: k
    INTEGER :: i1
    INTEGER :: zeit(8), i2, iseed,num_compute
    INTEGER, ALLOCATABLE :: seed(:)

    num_compute=nxprocs*nyprocs
!LS2011e
    
    ! for pgi-compiler, the system_clock starts at 0 when system_clock is first called during a program.
    ! So by default, it measures a time *difference*, anticipating that the user
    ! only wants to time his program. So, at the first call, the returned time is always 0.
    ! Only for the subsequent calls, the time increases. What a nonsense!

    ! Unfortunately, this is unusable for the purpose of initializing the random number generator with
    ! a different seed for every program run, since this will all times lead
    ! to the same result, independent of a certain
    ! random or varying component.
    ! This is different from other compilers, where system_clock() delivers the elapsed time
    ! since 1.1.1970 in milliseconds, modulo HUGE(int).
    !
    ! So, we do it differently:
    ! First, we use HUGE() to determine the maxint value i2:
    i2 = HUGE(i1)
    IF (.NOT.PRESENT(iseed_in)) THEN
      ! and then we use DATE_AND_TIME() to get the milliseconds part of the actual time,
      ! which later will serve as the varying component from program run to program run:
      !CALL DATE_AND_TIME(values=zeit)
      iseed=17!zeit(8)
    ELSE
      ! or, if it is desired, we use just iseed_in, which leeds to the same random numbers 
      ! everytime:
      iseed = iseed_in
    END IF

    ! get length of seed vector:
    CALL RANDOM_SEED(SIZE=k)
    ALLOCATE( seed(k) )
    seed = 0

    ! However, in any case we want to have a different random number series on each task,
    ! so this is achieved by merging in my_cart_id into the seed.
    ! The seed itself is constructed in a way that it is (multiply) folded into
    ! the number range of integer*4 data type, in order to break somehow the monotonicity
    ! in the seed vector. Monotonicity in the seed leads to a number series, from
    ! which the first 100 elements or so are not random but very close to 0, and only
    ! afterwards convert to more random behaviour.
    DO i=1,k
      seed(i) = myid+MOD(INT(i2/11*13*((MOD(i,5)+i)*i) + i2*int(real(iseed)/1000.0) + &
          i2*int(0.95/real(myid+1)), kind=KIND(i2)), i2)
    END DO
    CALL RANDOM_SEED(PUT=seed)

    IF (k >= 4) THEN
      WRITE(*,'(a,i3,a,4(x,i14))') '    SEED_RANDOM_NUMBER (first 4 of ',k,'):', seed(1:4)
    ELSE
      WRITE(*,*) '    SEED_RANDOM_NUMBER : ', seed
    END IF

    DEALLOCATE( seed )    
    RETURN
  END SUBROUTINE seed_random_number
! LINDA, e
end module init
