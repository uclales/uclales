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
  logical :: corflg = .false.
  logical :: rylflg = .true.

  real    :: frqhis =  9000.
  real    :: frqanl =  3600.
  real    :: radfrq =  0.

  real    :: time   =  0.
  real    :: strtim =  0.0
  real    :: cntlat =  31.5 ! 30.0
  logical :: outflg = .true.
  logical :: statflg = .false.
!irina  
  real    :: sst=292.
  real    :: div = 3.75e-6
  logical :: lsvarflg = .false.
  character (len=8) :: case_name = 'astex'

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

    use mpi_interface, only : myid, double_scalar_par_max

    use grid, only : dt, dtlong, zt, zm, nzp, dn0, u0, v0, level, &
         write_hist, write_anal, close_anal 
         
    use stat, only : savg_intvl, ssam_intvl, write_ps, close_stat
    use thrm, only : thermo

    real, parameter    :: peak_cfl = 0.5

    real    :: t1,t2,tplsdt,begtime,cflmax,gcflmax
    integer :: istp, iret

    !
    ! Timestep loop for program
    !
    begtime = time
    istp = 0
  !irina  
  !  timrsm = timrsm + begtime
   
  ! print *, 'time timrsm',time, timrsm
  !  timmax = min(timmax,timrsm)
  ! print *, 'timmax',timmax

    do while (time + 0.1*dt < timmax)

       call cpu_time(t1)           !t1=timing()

       istp = istp+1
       tplsdt = time + dt + 0.1*dt
       statflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dt  &
            .or. tplsdt >= timmax .or. tplsdt < 2.*dt) 

       call t_step
       time  = time + dt

       call cfl(cflmax)
       call double_scalar_par_max(cflmax,gcflmax)
       cflmax = gcflmax
       dt = min(dtlong,dt*peak_cfl/(cflmax+epsilon(1.)))
       !
       ! output control
       !
       if (mod(tplsdt,savg_intvl)<dt .or. time>=timmax .or. time==dt) &
       call write_ps(nzp,dn0,u0,v0,zm,zt,time)

       if ((mod(tplsdt,frqhis) < dt .or. time >= timmax) .and. outflg)   &
            call write_hist(2, time)
       !irina     
       !if (mod(tplsdt,savg_intvl)<dt .or. time>=timmax .or. time>=timrsm .or. time==dt)   &
       if (mod(tplsdt,savg_intvl)<dt .or. time>=timmax .or. time==dt)   &
            call write_hist(1, time)

!irina more frequent outputs for certain hours in astex

!       frqanl=3600.

!       if (time>=7200 .and. time<=10800) then
!       frqanl=300.
!       end if
 
!       if (time>=25200 .and. time<=28800) then
!       frqanl=300.
!       end if

!       if (time>=39600 .and. time<=43200) then
!       frqanl=300.
!       end if

!       if (time>=68400 .and. time<=72000) then
!       frqanl=300.
!       end if

!       if (time>=82800 .and. time<=86400) then
!       frqanl=300.
!       end if

!       if (time>=126000 .and. time<=129600) then
!       frqanl=300.
!       end if

       if ((mod(tplsdt,frqanl) < dt .or. time >= timmax) .and. outflg) then
          call thermo(level)
          call write_anal(time)
       end if

       if(myid == 0) then
          call cpu_time(t2)           !t1=timing()
          if (mod(istp,istpfl) == 0 ) print "('   Timestep # ',i5," //     &
              "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3,'   Est. CPU Time left(sec) = ',f10.2)",     &
              istp, time, t2-t1, t2*(timmax/time-1)
       endif

    enddo

    call write_hist(1, time)
    iret = close_anal()
    iret = close_stat()

  end subroutine stepper
  ! 
  !----------------------------------------------------------------------
  ! subroutine t_step: Called by driver to timestep through the LES
  ! routines.  Within many subroutines, data is accumulated during
  ! the course of a timestep for the purposes of statistical analysis.
  ! 
  subroutine t_step

    use grid, only : level, dt, nstep, a_tt, a_up, a_vp, a_wp, dxi, dyi, dzt, &
         nxp, nyp, nzp, dn0,a_scr1, u0, v0, a_ut, a_vt, a_wt, zt
    use stat, only : sflg, statistics
    use sgsm, only : diffuse
    !irina
    use srfc, only : surface
    !use srfc, only : surface,sst
    use thrm, only : thermo
    use mcrp, only : micro
    use prss, only : poisson
    use advf, only : fadvect
    use advl, only : ladvect
    use forc, only : forcings
    use lsvar, only : varlscale
    use util, only : velset,get_avg

    integer :: k
    real :: xtime
    character (len=11)    :: fname = 'debugXX.dat'

    xtime = time/86400. + strtim

    do nstep = 1,3

       ! Add additional criteria to ensure that some profile statistics that are  
       ! updated every 'ssam_intvl' outside the main statistics module
       ! are not updated (summed) in all three RK substeps.

       if (statflg .and. nstep.eq.3) then
          sflg = .True.
       end if

       call tendencies(nstep)
       call thermo(level)
!irina : should I called only at nstep=3 ?
       if (lsvarflg) then
      ! print *, 'step', lsvarflg
       call varlscale(time,case_name,sst,div,u0,v0)
       end if
       call surface(sst)      
!       
       call diffuse
       call fadvect
       call ladvect
       if (level >= 1) then
          call thermo(level)
!irina   
        ! print *, 'sst bef forcing',sst     
         !print *, "tt",a_tt(:,5,5)
          call forcings(xtime,cntlat,sst,div,case_name)
         !print *, "tt",a_tt(:,5,5)
          call micro(level)
       end if
       call corlos 
       call buoyancy
       call sponge
       call update (nstep)
       call poisson 
       call velset(nzp,nxp,nyp,a_up,a_vp,a_wp)

    end do
 
    if (statflg) then 
       call thermo (level)
       call statistics (time+dt)
       sflg = .False.
    end if

  end subroutine t_step
  ! 
  !----------------------------------------------------------------------
  ! subroutine tend0: sets all tendency arrays to zero
  ! 
  subroutine tendencies(nstep)

    use grid, only : a_ut, a_vt, a_wt, a_tt, a_rt, a_rpt, a_npt, a_ninuct, &
                     a_micet,a_nicet,a_msnowt,a_nsnowt, a_mgrt, a_ngrt,&
                     a_xt1, a_xt2, nscl, nxyzp, level
    use util, only : azero

    integer, intent (in) :: nstep
  
    select case(nstep)
    case default
  
       call azero(nxyzp*nscl,a_xt1)
       a_ut => a_xt1(:,:,:,1)
       a_vt => a_xt1(:,:,:,2)
       a_wt => a_xt1(:,:,:,3)
       a_tt => a_xt1(:,:,:,4)
       if (level >= 0) a_rt  =>a_xt1(:,:,:,5)
       if (level >= 3) then
          a_rpt =>a_xt1(:,:,:,6)
          a_npt =>a_xt1(:,:,:,7)
       end if
       if (level >= 4) then
          a_ninuct =>a_xt1(:,:,:, 8)
          a_micet  =>a_xt1(:,:,:, 9)
          a_nicet  =>a_xt1(:,:,:,10)
          a_msnowt =>a_xt1(:,:,:,11)
          a_nsnowt =>a_xt1(:,:,:,12)
          a_mgrt   =>a_xt1(:,:,:,13)
          a_ngrt   =>a_xt1(:,:,:,14)
       end if

    case(2)
       call azero(nxyzp*nscl,a_xt2)
       a_ut => a_xt2(:,:,:,1)
       a_vt => a_xt2(:,:,:,2)
       a_wt => a_xt2(:,:,:,3)
       a_tt => a_xt2(:,:,:,4)
       if (level >= 0) a_rt  =>a_xt2(:,:,:,5)
       if (level >= 3) then
          a_rpt =>a_xt2(:,:,:,6)
          a_npt =>a_xt2(:,:,:,7)
       end if
       if (level >= 4) then
          a_ninuct =>a_xt2(:,:,:, 8)
          a_micet  =>a_xt2(:,:,:, 9)
          a_nicet  =>a_xt2(:,:,:,10)
          a_msnowt =>a_xt2(:,:,:,11)
          a_nsnowt =>a_xt2(:,:,:,12)
          a_mgrt   =>a_xt2(:,:,:,13)
          a_ngrt   =>a_xt2(:,:,:,14)
       end if
    end select

  end subroutine tendencies
  ! 
  ! ----------------------------------------------------------------------
  ! subroutine update:
  !
  subroutine update(nstep)

!irina
    use grid, only : a_xp, a_xt1, a_xt2, a_up, a_vp, a_wp, a_sp, dzt, dt,  &
         nscl, nxp, nyp, nzp, newvar,level, a_rpp,a_npp
    use util, only : sclrset,velset

    real, parameter ::  alpha(3) = (/ 8./15., -17./60.,  3./4. /), &
         beta(3)  = (/    0.0,   5./12., -5./12./)

    integer, intent (in) :: nstep

    integer :: n

    a_xp = a_xp + dt *(alpha(nstep)*a_xt1 + beta(nstep)*a_xt2)

    call velset(nzp,nxp,nyp,a_up,a_vp,a_wp)

    do n=4,nscl
       call newvar(n,istep=nstep)
       call sclrset('mixd',nzp,nxp,nyp,a_sp,dzt)
    end do

!irina
   if (level >= 3) then
   where (a_rpp < 0.) 
    a_rpp=0.
    end where
   where (a_npp < 0.) 
    a_npp=0.
    end where
   end if


  end subroutine update
  ! 
  !----------------------------------------------------------------------
  ! Subroutine cfl: Driver for calling CFL computation subroutine
  ! 
  subroutine cfl(cflmax)

    use grid, only : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dt
    use stat, only : fill_scalar

    real, intent (out)   :: cflmax
      
    cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dt)
    call fill_scalar(1,cflmax)

  end subroutine cfl
  ! 
  !----------------------------------------------------------------------
  ! Subroutine cfll: Gets the peak CFL number
  ! 
  real function cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dt)

    integer, intent (in) :: n1, n2, n3
    real, dimension (n1,n2,n3), intent (in) :: u, v, w
    real, intent (in)    :: dxi,dyi,dzt(n1),dt

    integer :: i, j, k
    cfll=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             cfll=max(cfll, dt * max(abs(u(k,i,j)*dxi),                    &
                  abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
          end do
       end do
    end do

  end function cfll
  ! 
  ! ----------------------------------------------------------------------
  ! subroutine buoyancy:
  !
  subroutine buoyancy

    use grid, only : a_up, a_vp, a_wp, a_wt, vapor, a_theta, a_scr1, a_scr3,&
         a_rp, nxp, nyp, nzp, dzm, th00, level, pi1
    use stat, only : sflg, comp_tke
    use util, only : ae1mm
    use thrm, only : update_pi1

    real, dimension (nzp) :: awtbar

    call boyanc(nzp,nxp,nyp,level,a_wt,a_theta,a_rp,vapor,th00,a_scr1)
    call ae1mm(nzp,nxp,nyp,a_wt,awtbar)
    call update_pi1(nzp,awtbar,pi1)

    if (sflg)  call comp_tke(nzp,nxp,nyp,dzm,th00,a_up,a_vp,a_wp,a_scr1,a_scr3)

  end subroutine buoyancy
  ! 
  ! ----------------------------------------------------------------------
  ! subroutine boyanc:
  !
  subroutine boyanc(n1,n2,n3,level,wt,th,rt,rv,th00,scr)

    use defs, only: g, ep2

    integer, intent(in) :: n1,n2,n3,level
    real, intent(in)    :: th00,th(n1,n2,n3),rt(n1,n2,n3),rv(n1,n2,n3)
    real, intent(inout) :: wt(n1,n2,n3)
    real, intent(out)   :: scr(n1,n2,n3)

    integer :: k, i, j
    real :: gover2

    gover2  = 0.5*g
 
    do j=3,n3-2
       do i=3,n2-2
          if (level >= 2) then
             do k=1,n1
                scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rv(k,i,j))-th00)       &
                     /th00-(rt(k,i,j)-rv(k,i,j)))
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
