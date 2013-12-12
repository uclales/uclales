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
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module step

  implicit none

  integer :: istpfl = 1
  real    :: timmax = 18000.
  real    :: timrsm = 86400.
  real    :: wctime = 1.e10
  logical :: corflg = .false.
  logical :: rylflg = .true.

  real    :: frqhis =  9000.
  real    :: frqanl =  3600.
  real    :: frqcross =  3600.
  real    :: radfrq =  0.

  real    :: time   =  0.
  real    :: strtim =  0.0
  real    :: cntlat =  31.5 ! 30.0
  logical :: outflg = .true.
  logical :: statflg = .false.
  real    :: tau = 900.
!irina
  real    :: sst=292.
  real    :: div = 0.0
  logical :: lsvarflg = .false.
  character (len=8) :: case_name = 'astex'

  integer :: istp
! linda,b
  logical ::lanom=.false.
!linda,e

  ! Flags for sampling, statistics output, etc.
  logical :: savgflg=.false.,anlflg=.false.,hisflg=.false.,crossflg=.false.,lpdumpflg=.false.

contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine model:  This is the main driver for the model's time
  ! integration.  It calls the routine tstep, which steps through the
  ! physical processes active on a time-step and updates variables.  It
  ! then checks to see whether or not different output options are
  ! satisfied.
  !
  subroutine stepper

    use mpi_interface, only : myid, broadcast_dbl, double_scalar_par_max,mpi_get_time
    use grid, only : dt, dtlong, zt, zm, nzp, dn0, u0, v0, level, &
         write_hist
    use ncio, only : write_anal, close_anal
    use modcross, only : triggercross, exitcross, lcross
    use stat, only : savg_intvl, ssam_intvl, write_ps, close_stat
    use thrm, only : thermo
    use modparticles, only : lpartic, exit_particles, lpartdump, exitparticledump, &
         lpartstat, exitparticlestat, write_particle_hist, particlestat, &
	 balanced_particledump,frqpartdump, ldropstart


    real, parameter    :: peak_cfl = 0.5, peak_peclet = 0.5

    real    :: tplsdt,begtime,cflmax,gcflmax,pecletmax,gpecletmax
    double precision    :: t0,t1,t2
    integer :: iret

    real    :: dt_prev

    !
    ! Timestep loop for program
    !
    begtime = time
    istp    = 0
    t2      = 0.
    call mpi_get_time(t0)
    call broadcast_dbl(t0, 0)
    do while (time < timmax .and. (t2-t0) < wctime)
       call mpi_get_time(t1)
       istp = istp + 1

       call stathandling
       if(myid .eq. 0 .and. statflg) print*,'     sampling stat at t=',time+dt

       call t_step
       time = time + dt

       call cfl(cflmax)
       call double_scalar_par_max(cflmax,gcflmax)
       cflmax = gcflmax
       call peclet(pecletmax)
       call double_scalar_par_max(pecletmax,gpecletmax)
       pecletmax = gpecletmax
       dt_prev = dt
       dt = min(dtlong,dt*peak_cfl/(cflmax+epsilon(1.)))

       !
       ! output control
       !
       !! Sample particles; automatically samples when savgflg=.true., don't sample double...
       !if(lpartic .and. lpartstat .and. statflg .and. (savgflg .eqv. .false.)) call particlestat(.false.,time+dt)

       if(savgflg) then
         if(myid==0) print*,'     profiles at time=',time
         call write_ps(nzp,dn0,u0,v0,zm,zt,time)
         if (lpartic .and. lpartstat) call particlestat(.true.,time)
       end if

       if (hisflg) then
         if(myid==0) print*,'     history at time=',time
         call write_hist(2, time)
         if(lpartic) call write_particle_hist(2,time)
       end if

       if (anlflg) then
         if(myid==0) print*,'     analysis at time=',time
         call thermo(level)
         call write_anal(time)
       end if

       if (crossflg) then
         if(myid==0) print*,'     cross at time=',time
         call thermo(level)
         call triggercross(time)
       end if

       if (lpdumpflg.and.(time.ge.ldropstart)) then
         if(myid==0) print*,'     particle dump at time=',time
         call balanced_particledump(time)
       end if

       statflg   = .false.
       savgflg   = .false.
       hisflg    = .false.
       anlflg    = .false.
       crossflg  = .false.
       lpdumpflg = .false.
       
       
       ! REMOVE THIS?
       !irina
       !if (mod(tplsdt,savg_intvl)<dt .or. time>=timmax .or. time>=timrsm .or. time==dt)   &
       !if (mod(tplsdt,savg_intvl)<dt .or. time>=timmax .or. time==dt) then
       !  call write_hist(1, time)
       !  if(lpartic) call write_particle_hist(1,time)
       !end if

       if(myid == 0) then
          call mpi_get_time(t2)
          if (mod(istp,istpfl) == 0 ) then
              if (wctime.gt.1e9) then
                print "('   Timestep # ',i6," //     &
                       "'   Model time(sec)=',f12.2,3x,'dt(sec)=',f8.4,'   CPU time(sec)=',f8.3)",     &
                       istp, time, dt_prev, t2-t1
              else
                print "('   Timestep # ',i6," //     &
                       "'   Model time(sec)=',f12.2,3x,'dt(sec)=',f8.4,'   CPU time(sec)=',f8.3'  WC Time left(sec) = ',f10.2)",     &
                       istp, time, dt_prev, t2-t1, wctime-t2+t0
              end if
          end if
       endif

       call broadcast_dbl(t2, 0)
    enddo

    call write_hist(1, time)

    iret = close_anal()

    if (lpartic) then
      call write_particle_hist(1, time)
      call exit_particles
      if(lpartdump) call exitparticledump
      if(lpartstat) call exitparticlestat
    end if

    if (lcross) call exitcross

    iret = close_stat()

    if ((t2-t0) .ge. wctime .and. myid == 0) write(*,*) '  Wall clock limit wctime reached, stopped simulation for restart'
    if (time.ge.timmax .and. myid == 0) write(*,*) '  Max simulation time timmax reached. Finished simulation successfully'

  end subroutine stepper

  !
  ! -------------------------------------------
  ! subroutine stathandling: Event handling statistics and output.
  ! Changed dt when nextevent-time < dt to nextevent-time to end
  ! up exactly at the required sampling or output time
  ! EXPERIMENTAL!! BvS, Sep2012
  !
  subroutine stathandling
    use grid, only          : dt
    use defs, only          : long
    use stat, only          : savg_intvl, ssam_intvl
    use modcross, only      : lcross
    use mpi_interface, only : myid
    use modparticles, only  : lpartic,lpartdump,frqpartdump

    integer(kind=long)  :: itnssam=1e16,itnsavg=1e16,itnanl=1e16,itnhist=1e16,itncross=1e16,itnlpdump=1e16
    integer(kind=long)  :: issam_intvl=1e16,isavg_intvl=1e16,ifrqanl=1e16,ifrqhis=1e16,ifrqcross=1e16,ifrqlpdump=1e16
    integer(kind=long)  :: itime,idt
    real, parameter     :: tres = 1e9

    ! Switch to integer microseconds
    issam_intvl  = int(ssam_intvl   *tres,long)
    isavg_intvl  = int(savg_intvl   *tres,long)
    ifrqanl      = int(frqanl       *tres,long)
    ifrqhis      = int(frqhis       *tres,long)
    ifrqcross    = int(frqcross     *tres,long)
    ifrqlpdump   = int(frqpartdump  *tres,long)

    itime        = int(time * tres,long)
    idt          = int(dt   * tres,long)

    ! Time next events
    itnssam      = (floor(itime/(ssam_intvl*tres))  + 1) * issam_intvl
    itnsavg      = (floor(itime/(savg_intvl*tres))  + 1) * isavg_intvl
    itnanl       = (floor(itime/(frqanl*tres))      + 1) * ifrqanl
    itnhist      = (floor(itime/(frqhis*tres))      + 1) * ifrqhis
    if(lcross) &
      itncross   = (floor(itime/(frqcross*tres))    + 1) * ifrqcross
    if(lpartic .and. lpartdump) &
      itnlpdump  = (floor(itime/(frqpartdump*tres)) + 1) * ifrqlpdump

    ! Limit time step to first event
    idt = min(itnssam-itime,itnsavg-itime,&                    ! Profile sampling and writing
               itnanl-itime,itnhist-itime,itncross-itime,&     ! Analysis, history and cross-sections
               itnlpdump-itime,&                               ! Lagrangian particles
               int(timmax*tres,long)-itime,idt)                ! End of simulation, current dt

    ! And back to normal seconds
    dt   = idt  / tres
    time = itime / tres        ! This can be tricky.................

    ! Set flags
    statflg  = .false.
    savgflg  = .false.
    hisflg   = .false.
    anlflg   = .false.
    crossflg = .false.
    lpdumpflg= .false.

    if(mod(itime+idt,issam_intvl) .eq. 0) then
      statflg    = .true.
    end if
    if(mod(itime+idt,isavg_intvl) .eq. 0) then
      statflg    = .true.
      savgflg    = .true.
    end if
    if(mod(itime+idt,ifrqanl)     .eq. 0) then
      anlflg     = .true.
    end if
    if(mod(itime+idt,ifrqhis)     .eq. 0) then
      hisflg     = .true.
    end if
    if((mod(itime+idt,ifrqcross)   .eq. 0) .and. lcross) then
      crossflg   = .true.
    end if
    if((mod(itime+idt,ifrqlpdump)  .eq. 0) .and. lpartic .and. lpartdump) then
      lpdumpflg  = .true.
    end if

    !if(myid.eq.0) print*,statflg,savgflg,anlflg,hisflg,crossflg,lpdumpflg
    !stop
    !if(myid.eq.0) print*,'dt old:new',itime/tres,idtt/tres,idt/tres,(itnssam-itime)/tres,statflg

  end subroutine stathandling

  !
  !----------------------------------------------------------------------
  ! subroutine t_step: Called by driver to timestep through the LES
  ! routines.  Within many subroutines, data is accumulated during
  ! the course of a timestep for the purposes of statistical analysis.
  !
  subroutine t_step

    use mpi_interface, only : myid, appl_abort
    use grid, only : level, dt, nstep, a_tt, a_up, a_vp, a_wp, dxi, dyi, dzi_t, &
         nxp, nyp, nzp, dn0,a_scr1, u0, v0, a_ut, a_vt, a_wt, zt, a_ricep, a_rct, a_rpt, &
         lwaterbudget, a_xt2
    use stat, only : sflg, statistics
    use sgsm, only : diffuse
    !use sgsm_dyn, only : calc_cs
    use srfc, only : surface
    use thrm, only : thermo
    use mcrp, only : micro, lpartdrop
    use prss, only : poisson
    use advf, only : fadvect
    use advl, only : ladvect
    use forc, only : forcings
    use lsvar, only : varlscale
    use util, only : velset,get_avg
    use centered, only:advection_scalars
    use modtimedep, only : timedep
    use modparticles, only : particles, lpartic, particlestat,lpartstat, &
         deactivate_drops, activate_drops, lpartmass, grow_drops


    logical, parameter :: debug = .false.
    real :: xtime
    character (len=8) :: adv='monotone'

    xtime = time/86400. + strtim
    call timedep(time,timmax, sst)

!    !new for mass budget
!    if(lpartic .and. lpartdrop .and. lpartmass) call grow_drops
!    if(lpartic .and. lpartdrop) call deactivate_drops(time+dt)
!    if(lpartic .and. lpartdrop) call activate_drops(time+dt)

    do nstep = 1,3

       ! Add additional criteria to ensure that some profile statistics that are
       ! updated every 'ssam_intvl' outside the main statistics module
       ! are not updated (summed) in all three RK substeps.

       if (statflg .and. nstep.eq.3) then
          sflg = .True.
       end if

       ! old position not suitable for mass budget
       if (lpartic) then
         call particles(time,timmax)
       end if

       call tendencies(nstep)
       call thermo(level)

       if (lsvarflg) then
         call varlscale(time,case_name,sst,div,u0,v0)
       end if

       xtime = xtime - strtim
       call surface(sst,xtime)
       xtime = xtime + strtim

       call diffuse(time)

       if (adv=='monotone') then
          call fadvect
       elseif ((adv=='second').or.(adv=='third').or.(adv=='fourth')) then
          call advection_scalars(adv)
       else 
          print *, 'wrong specification for advection scheme'
          call appl_abort(0)
       endif

       call ladvect

       if (level >= 1) then
          if (lwaterbudget) then
             call thermo(level,1)
          else
             call thermo(level)
          end if
          call forcings(xtime,cntlat,sst,div,case_name,time)
          call micro(level,istp)
       end if

       call corlos
       call buoyancy
       call sponge
       call decay
       call update (nstep)
       call poisson
       call velset(nzp,nxp,nyp,a_up,a_vp,a_wp)
       
!       !new position for mass budget
!       if (lpartic) then
!         call particles(time,timmax)
!       end if

    end do

    !old not suitable for mass budget
    if(lpartic .and. lpartdrop .and. lpartmass) call grow_drops
    if(lpartic .and. lpartdrop) call deactivate_drops(time+dt)
    if(lpartic .and. lpartdrop) call activate_drops(time+dt)

    if (statflg) then
       if (debug) WRITE (0,*) 't_step statflg thermo, myid=',myid
       call thermo (level)
       if (debug) WRITE (0,*) 't_step statflg statistics, myid=',myid
       call statistics (time+dt)
       if(lpartic .and. lpartstat) call particlestat(.false.,time+dt)
       sflg = .False.
    end if

  end subroutine t_step
  !
  !----------------------------------------------------------------------
  ! subroutine tend0: sets all tendency arrays to zero
  !
  subroutine tendencies(nstep)

    use grid, only : a_ut, a_vt, a_wt, a_tt, a_rt, a_rpt, a_npt, &
                     a_ricet,a_nicet,a_rsnowt, a_rgrt,&
                     a_rhailt,a_nhailt,a_nsnowt, a_ngrt,&
                     a_xt1, a_xt2, nscl, nxyzp, level, &
                     lwaterbudget, a_rct, ncld, &
                     lcouvreux, a_cvrxt, ncvrx
    use util, only : azero

    integer, intent (in) :: nstep

    select case(nstep)
    case default

       call azero(nxyzp*nscl,a_xt1)
       a_ut => a_xt1(:,:,:,1)
       a_vt => a_xt1(:,:,:,2)
       a_wt => a_xt1(:,:,:,3)
       a_tt => a_xt1(:,:,:,4)
       if (level > 0) a_rt  =>a_xt1(:,:,:,5)
       if (level >= 3) then
          a_rpt =>a_xt1(:,:,:,6)
          a_npt =>a_xt1(:,:,:,7)
       end if
       if (lwaterbudget) a_rct =>a_xt1(:,:,:,ncld)
       if (lcouvreux)    a_cvrxt =>a_xt1(:,:,:,ncvrx)
       if (level >= 4) then
          a_ricet  =>a_xt1(:,:,:, 8)
          a_nicet  =>a_xt1(:,:,:, 9)
          a_rsnowt =>a_xt1(:,:,:,10)
          a_rgrt   =>a_xt1(:,:,:,11)
       end if
       if (level >= 5) then
          a_nsnowt =>a_xt1(:,:,:,12)
          a_ngrt   =>a_xt1(:,:,:,13)
          a_rhailt =>a_xt1(:,:,:,14)
          a_nhailt =>a_xt1(:,:,:,15)
       end if

    case(2)
       call azero(nxyzp*nscl,a_xt2)
       a_ut => a_xt2(:,:,:,1)
       a_vt => a_xt2(:,:,:,2)
       a_wt => a_xt2(:,:,:,3)
       a_tt => a_xt2(:,:,:,4)
       if (level > 0) a_rt  =>a_xt2(:,:,:,5)
       if (level >= 3) then
          a_rpt =>a_xt2(:,:,:,6)
          a_npt =>a_xt2(:,:,:,7)
       end if
       if (lwaterbudget) a_rct =>a_xt2(:,:,:,ncld)
       if (lcouvreux)    a_cvrxt =>a_xt2(:,:,:,ncvrx)
       if (level >= 4) then
          a_ricet  =>a_xt2(:,:,:, 8)
          a_nicet  =>a_xt2(:,:,:, 9)
          a_rsnowt =>a_xt2(:,:,:,10)
          a_rgrt   =>a_xt2(:,:,:,11)
       end if
       if (level >= 5) then
          a_nsnowt =>a_xt2(:,:,:,12)
          a_ngrt   =>a_xt2(:,:,:,13)
          a_rhailt =>a_xt2(:,:,:,14)
          a_nhailt =>a_xt2(:,:,:,15)
       end if

    end select

  end subroutine tendencies
  !
  ! ----------------------------------------------------------------------
  ! subroutine update:
  !
  subroutine update(nstep)

!irina
    use grid, only : a_xp, a_xt1, a_xt2, a_up, a_vp, a_wp, a_sp, dzi_t, dt,  &
         nscl, nxp, nyp, nzp, newvar,level, a_rpp,a_ricep,a_nicep,a_rsnowp,a_rgrp,a_npp,rkalpha,rkbeta, &
         a_nsnowp,a_ngrp,a_rhailp,a_nhailp,a_rp,liquid
    use util, only : sclrset,velset

    integer, intent (in) :: nstep

    integer :: n

    a_xp = a_xp + dt *(rkalpha(nstep)*a_xt1 + rkbeta(nstep)*a_xt2)

    call velset(nzp,nxp,nyp,a_up,a_vp,a_wp)

    do n=4,nscl
       call newvar(n,istep=nstep)
       call sclrset('mixd',nzp,nxp,nyp,a_sp,dzi_t)
       ! for debugging
       !call sclrset('mixd',nzp,nxp,nyp,a_sp,dzi_t,n)
    end do

    if (level >= 1) then
       where (a_rp < 0.) 
          a_rp=0.
       end where
    end if
    if (level >= 3) then
       a_rpp(1,:,:) = 0.
       a_npp(1,:,:) = 0.
       where (a_rpp < 0.)
          a_rpp=0.
       end where
       where (a_npp < 0.)
          a_npp=0.
       end where
       where (liquid < 0.)
          liquid=0.
       end where
    end if
    if (level >= 4) then
       a_ricep(1,:,:) = 0.
       a_nicep(1,:,:) = 0.
       a_rsnowp(1,:,:) = 0.
       a_rgrp(1,:,:) = 0.
       where (a_ricep < 0.)
          a_ricep=0.
       end where
       where (a_ricep < 0.)
          a_ricep=0.
       end where
       where (a_ricep <= 0.)
          a_nicep=0.
       end where
       where (a_nicep < 0.)
          a_nicep=0.
       end where
       where (a_rsnowp < 0.)
          a_rsnowp=0.
       end where
       where (a_rgrp < 0.)
          a_rgrp=0.
       end where
    end if
    if (level >= 5) then
       a_nsnowp(1,:,:) = 0.
       a_ngrp(1,:,:)   = 0.
       a_rhailp(1,:,:) = 0.
       a_nhailp(1,:,:) = 0.
       where (a_rhailp < 0.)
          a_rhailp=0.
       end where
       where (a_rhailp <= 0.)
          a_nhailp=0.
       end where
       where (a_rsnowp <= 0.)
          a_nsnowp=0.
       end where
       where (a_rgrp <= 0.)
          a_ngrp=0.
       end where
       where (a_nhailp < 0.)
          a_nhailp=0.
       end where
       where (a_nsnowp < 0.)
          a_nsnowp=0.
       end where
       where (a_ngrp < 0.)
          a_ngrp=0.
       end where
    end if

  end subroutine update
  !
  !----------------------------------------------------------------------
  ! Subroutine cfl: Driver for calling CFL computation subroutine
  !
  subroutine cfl(cflmax)

    use grid, only : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzi_t,dt
    use stat, only : fill_scalar

    real, intent (out)   :: cflmax

    cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzi_t,dt)
    call fill_scalar(1,cflmax)

  end subroutine cfl
  !
  !----------------------------------------------------------------------
  ! Subroutine cfll: Gets the peak CFL number
  !
  real function cfll(n1,n2,n3,u,v,w,dxi,dyi,dzi_t,dt)

    integer, intent (in) :: n1, n2, n3
    real, dimension (n1,n2,n3), intent (in) :: u, v, w
    real, intent (in)    :: dxi,dyi,dzi_t(n1),dt

    integer :: i, j, k
    cfll=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             cfll=max(cfll, dt * max(abs(u(k,i,j)*dxi),                    &
                  abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzi_t(k))))
          end do
       end do
    end do

  end function cfll
  !
  !----------------------------------------------------------------------
  ! Subroutine peclet: Driver for calling peclet computation subroutine
  !
  subroutine peclet(pecletmax)

    use grid, only : a_km,nxp,nyp,nzp,dxi,dyi,dzi_t,dt
    use stat, only : fill_scalar

    real, intent (out)   :: pecletmax

    pecletmax =  pecletl(nzp,nxp,nyp,a_km,dxi,dyi,dzi_t,dt)
    call fill_scalar(1,pecletmax)

  end subroutine peclet
  !
  !----------------------------------------------------------------------
  ! Subroutine pecletl: Gets the peak peclet number
  !
  real function pecletl(n1,n2,n3,km,dxi,dyi,dzi_t,dt)

    integer, intent (in) :: n1, n2, n3
    real, dimension (n1,n2,n3), intent (in) :: km
    real, intent (in)    :: dxi,dyi,dzi_t(n1),dt

    integer :: i, j, k
    pecletl=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             pecletl=max(pecletl, dt * km(k,i,j)*max(dxi,dyi,dzi_t(k)))
          end do
       end do
    end do

  end function pecletl
  !
  ! ----------------------------------------------------------------------
  ! subroutine buoyancy:
  !
  subroutine buoyancy

    use grid, only : a_up, a_vp, a_wp, a_wt, vapor, a_theta, a_scr1, a_scr3,liquid,&
         a_rp,a_rpp,a_ricep, a_rsnowp, a_rgrp, a_rhailp, nxp, nyp, nzp, dzi_m, th00, level, pi1
    use stat, only : sflg, comp_tke
    use util, only : ae1mm
    use thrm, only : update_pi1
    real, dimension(nzp,nxp,nyp) :: rl
    real, dimension (nzp) :: awtbar
    rl = 0.
    if (level>1) rl = liquid
    if (level>2) rl = rl + a_rpp
    if (level>3) rl = rl + a_ricep + a_rsnowp + a_rgrp
    if (level>4) rl = rl + a_rhailp

    if(level>0) then
      call boyanc(nzp,nxp,nyp,level,a_wt,a_theta,th00,a_scr1,vapor,rl)
    else
      call boyanc(nzp,nxp,nyp,level,a_wt,a_theta,th00,a_scr1)
    end if

    call ae1mm(nzp,nxp,nyp,a_wt,awtbar)
    call update_pi1(nzp,awtbar,pi1)

    if (sflg)  call comp_tke(nzp,nxp,nyp,dzi_m,th00,a_up,a_vp,a_wp,a_scr1,a_scr3)

  end subroutine buoyancy
  !
  ! ----------------------------------------------------------------------
  ! subroutine boyanc:
  !
  subroutine boyanc(n1,n2,n3,level,wt,th,th00,scr,rv,rl)

    use defs, only: g, ep2

    integer, intent(in) :: n1,n2,n3,level
    real, intent(in)    :: th00,th(n1,n2,n3)
    real, intent(inout) :: wt(n1,n2,n3)
    real, intent(out)   :: scr(n1,n2,n3)
    real, intent(in), optional :: rv(n1,n2,n3),rl(n1,n2,n3)

    integer :: k, i, j
    real :: gover2

    gover2  = 0.5*g

    do j=3,n3-2
       do i=3,n2-2
          if (level >= 2) then
             do k=1,n1
                scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rv(k,i,j))-th00)       &
                     /th00-(rl(k,i,j)))
             end do
          else
             do k=1,n1
                scr(k,i,j)=gover2*(th(k,i,j)/th00-1.)
             end do
          end if

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)+scr(k,i,j)+scr(k+1,i,j)
          end do
       end do
    end do

  end subroutine boyanc
  !
  ! ----------------------------------------------------------------------
  ! subroutine corlos:  This is the coriolis driver, its purpose is to
  ! from the coriolis accelerations for u and v and add them into the
  ! accumulated tendency arrays of ut and vt.
  !
  subroutine corlos

    use defs, only : omega
    use grid, only : a_up, a_vp, a_ut, a_vt, nxp, nyp, nzp, u0, v0

    logical, save :: initialized = .False.
    real, save    :: fcor

    integer :: i, j, k

    if (corflg) then
       if (.not.initialized) fcor=2.*omega*sin(cntlat*0.01745329)
       do j=3,nyp-2
          do i=3,nxp-2
             do k=2,nzp
                a_ut(k,i,j)=a_ut(k,i,j) - fcor*(v0(k)-0.25*                  &
                     (a_vp(k,i,j)+a_vp(k,i+1,j)+a_vp(k,i,j-1)+a_vp(k,i+1,j-1)))
                a_vt(k,i,j)=a_vt(k,i,j) + fcor*(u0(k)-0.25*                  &
                     (a_up(k,i,j)+a_up(k,i-1,j)+a_up(k,i,j+1)+a_up(k,i-1,j+1)))
             end do
          end do
       end do
       initialized = .True.
    end if

  end subroutine corlos
!
! ----------------------------------------------------------------------
! subroutine sponge: does the rayleigh friction for the momentum terms,
! and newtonian damping of thermal term the damping is accumulated with the
! other tendencies
!
  subroutine sponge

    use grid, only : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
         nfpt, spngt, spngm, nzp, nxp, nyp, th0, th00

    integer :: i, j, k, kk

    if (maxval(spngt) > epsilon(1.) .and. nfpt > 1) then
       do j=3,nyp-2
          do i=3,nxp-2
             do k=nzp-nfpt,nzp-1
                kk = k+1-(nzp-nfpt)
                a_tt(k,i,j)=a_tt(k,i,j)-spngt(kk)*(a_tp(k,i,j)-th0(k)+th00)
                a_ut(k,i,j)=a_ut(k,i,j)-spngt(kk)*(a_up(k,i,j)-u0(k))
                a_vt(k,i,j)=a_vt(k,i,j)-spngt(kk)*(a_vp(k,i,j)-v0(k))
                a_wt(k,i,j)=a_wt(k,i,j)-spngm(kk)*(a_wp(k,i,j))
             end do
          end do
       end do
    end if

  end subroutine sponge

  subroutine decay
    use grid, only : lcouvreux, a_cvrxp, a_cvrxt, nxp, nyp, nzp, dt
    integer :: i, j, k
    real    :: rate
    if (lcouvreux) then
      rate = 1./(max(tau, dt))
      do j = 3, nyp - 2
        do i = 3, nxp - 2
          do k = 2, nzp
            a_cvrxt(k,i,j) = a_cvrxt(k,i,j) - a_cvrxp(k,i,j) * rate
          end do
        end do
      end do
    end if

  end subroutine

  !
  ! --------------------------------------------------------------------
  ! subroutine get_diverg: gets velocity tendency divergence and puts it
  ! into a complex value array for use in pressure calculation
  !
  real function divergence(n1,n2,n3,u,v,w,dn0,dz,dx,dy)

    integer, intent (in)  :: n1,n2,n3
    real, intent (in)     :: dz(n1),dn0(n1),dx,dy
    real, dimension (n1,n2,n3), intent (in) :: u, v, w

    real :: s1(n2,n3,n1)

    integer :: k,i,j,l,m
    real    :: xf,yf,zf,wf1,wf2

    s1(:,:,:) = 0.0
    m=0
    do j=3,n3-2
       m=m+1
       l=0
       do i=3,n2-2
          l=l+1
          do k=2,n1-1
             wf1=0.5*(dn0(k+1)+dn0(k))
             wf2=0.5*(dn0(k)+dn0(k-1))
             if (k == 2 )   wf2=0.
             if (k == n1-1) wf1=0.
             xf=dn0(k)*dx
             yf=dn0(k)*dy
             zf=dz(k)
             s1(i,j,k)= (wf1*w(k,i,j)-wf2*w(k-1,i,j))*zf                     &
                  +(v(k,i,j)-v(k,i,j-1))*yf+(u(k,i,j)-u(k,i-1,j))*xf
          enddo
       enddo
    enddo
!
! save mxdiv to a statistics array, no reduction necessary as this is done
! in post processing
!
    divergence = maxval(s1)

  end function divergence

end module step
