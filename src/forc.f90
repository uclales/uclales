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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module forc

  use defs, only      : cp, pi
  use radiation, only : d4stream,surfacerad
  !irina
  use rad_gcss, only  : gcss_rad
  !cgils
  use grid, only      : wfls, dthldtls, dqtdtls, sfc_albedo, lrad_ca
  use modnudge, only  : nudge, nudge_bound, lnudge_bound !LINDA 
  use stat, only : sflg
  implicit none

  !character (len=5), parameter :: case_name = 'xxxx'
  !irina
  !initializes the large scale forcing
  integer, parameter    :: nls = 500
  real, dimension(nls)  :: t_ls=0.
  real, dimension(nls)  :: div_ls=0.
  real, dimension(nls)  :: sst_ls=0.
  real, dimension(nls)  :: ugeo_ls=0.
  real, dimension(nls)  :: vgeo_ls=0.
  !cgils
  logical :: lstendflg=.false.

contains
  !
  ! -------------------------------------------------------------------
  ! subroutine forcings:  calls the appropriate large-scale forcings
  !irina
  subroutine forcings(time_in, cntlat, sst, div, case_name,time_in_2)

!irina
    use grid, only: nxp, nyp, nzp, zm, zt, dzi_t, dzi_m, dn0, iradtyp, isfctyp, liquid  &
         , a_rflx, a_sflx, albedo, a_tt, a_tp, a_rt, a_rp, a_pexnr, a_scr1 &
         , vapor, a_rpp,a_ricep,a_nicep,a_rgrp, CCN, pi0, pi1, level, a_ut, a_up, a_vt, a_vp,a_theta,&
          a_lflxu, a_lflxd, a_sflxu, a_sflxd,sflxu_toa,sflxd_toa,lflxu_toa,lflxd_toa,a_lflxu_ca, &
          a_lflxd_ca, a_sflxu_ca, a_sflxd_ca, lflxd_toa_ca, lflxu_toa_ca, sflxd_toa_ca, sflxu_toa_ca, a_wt, xt

    use mpi_interface, only : myid, appl_abort
    use util, only : get_avg

!irina
    real, optional, intent (in) :: time_in, cntlat, sst, div, time_in_2
    real, dimension (nzp):: um,vm

    character (len=8), intent (in) :: case_name
    integer :: i, j, k, kp1
    real, dimension(nzp,nxp,nyp) :: dum0, dum1, dum2, dum3, dum4

!irina
    real :: alpha, beta, gamma, delta, uinf !linda

    select case(iradtyp)
    case (1)
        call case_forcing(nzp,nxp,nyp,case_name,zt,dzi_t,dzi_m,a_tp,a_rp,a_tt,a_rt)
    case (2)
       select case(level)
       case(1)
          call smoke_rad(nzp, nxp, nyp, dn0, a_rflx, zm, dzi_t,a_tt,a_rp)
       case(2)
          call gcss_rad(nzp, nxp, nyp, cntlat, time_in, case_name, div, sst, liquid, dn0,   &
              a_rflx, a_sflx, zt, zm, dzi_t, a_tt, a_tp, a_rt, a_rp)
       case(3)
          call gcss_rad(nzp, nxp, nyp, cntlat, time_in, case_name, div, sst, liquid, dn0,   &
              a_rflx, a_sflx, zt, zm, dzi_t, a_tt, a_tp, a_rt, a_rp)
       end select
    case (3)
       call bellon(nzp, nxp, nyp, a_rflx, a_sflx, zt, dzi_t, dzi_m, a_tt, a_tp&
            ,a_rt, a_rp, a_ut, a_up, a_vt, a_vp)
    case (4)
       if (present(time_in) .and. present(cntlat) .and. present(sst)) then
          !irina
          !a_scr1 = a_theta/a_pexnr
          select case (level)
          case(3)
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, CCN,   &
                  dn0, pi0, pi1, dzi_t, a_pexnr, a_theta, vapor, liquid, a_tt,&
                  a_rflx, a_sflx, a_lflxu, a_lflxd,a_sflxu,a_sflxd, albedo, &
                  rr=a_rpp,sflxu_toa=sflxu_toa,sflxd_toa=sflxd_toa,&
                  lflxu_toa=lflxu_toa,lflxd_toa=lflxd_toa)
          case(4,5)
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, CCN,   &
                  dn0, pi0, pi1, dzi_t, a_pexnr, a_theta, vapor, liquid, a_tt,&
                  a_rflx, a_sflx, a_lflxu, a_lflxd,a_sflxu,a_sflxd, albedo, &
                  rr=a_rpp,sflxu_toa=sflxu_toa,sflxd_toa=sflxd_toa,&
                  lflxu_toa=lflxu_toa,lflxd_toa=lflxd_toa,ice=a_ricep,nice=a_nicep,grp=a_rgrp)
          case default
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, CCN,    &
                  dn0, pi0, pi1, dzi_t, a_pexnr, a_theta, vapor, liquid, a_tt, &
                  a_rflx, a_sflx,a_lflxu, a_lflxd,a_sflxu,a_sflxd,albedo, &
                  sflxu_toa=sflxu_toa,sflxd_toa=sflxd_toa,&
                  lflxu_toa=lflxu_toa,lflxd_toa=lflxd_toa)
          end select
          if (sflg .and. lrad_ca) then
            dum0 = 0.
            dum1 = 0.
            dum2 = 0.
            dum3 = 0.
            dum4 = 0.
            call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,   &
               dn0, pi0, pi1, dzi_t, a_pexnr, a_theta, vapor, dum0, dum1,&
               dum2, dum3, a_lflxu_ca, a_lflxd_ca,a_sflxu_ca,a_sflxd_ca, albedo, &
               rr=dum4,sflxu_toa=sflxu_toa_ca,sflxd_toa=sflxd_toa_ca,&
               lflxu_toa=lflxu_toa_ca,lflxd_toa=lflxd_toa_ca)
          end if
        else
          if (myid == 0) print *, '  ABORTING: inproper call to radiation'
          call appl_abort(0)
       end if

    ! BvS: Simple parameterized surface radiation
    ! Sw-up/down = f(lat,lon,doy,tUTC,albedo)
    ! Lw-up/down = Boltzman
    case (5)
      if(isfctyp .ne. 5) then
        print*,'surface rad only works with interactive land surface'
        stop
      end if
      if (present(time_in) .and. present(cntlat)) then
        call surfacerad(cntlat,time_in)
      else
        print*,'improper call surfacrad, stopping'
        stop
      end if


    end select
!irina
          !
          ! subsidence
          !
!cgils
! LINDA, b
    if (lnudge_bound) call nudge_bound
! LINDA, e
    if (lstendflg) then

      do j=3,nyp-2
          do i=3,nxp-2
             do k=2,nzp-2
                kp1 = k+1
               a_tt(k,i,j) = a_tt(k,i,j) - &
                       wfls(k)*(a_tp(kp1,i,j)-a_tp(k,i,j))*dzi_t(k)+dthldtls(k)
               a_rt(k,i,j)=a_rt(k,i,j) - &
                       wfls(k)*(a_rp(kp1,i,j)-a_rp(k,i,j))*dzi_t(k)+dqtdtls(k)
             end do
          enddo
       enddo
    end if

!cgils: Nudging
    call nudge(time_in_2)

    if (case_name == 'squall') then

       ! add u forcing to initialize squall line. Proposed change from George Bryan, 12-12-12

      alpha = 0.1

      IF(time_in_2.le.3600.0)THEN
        gamma = 1.0
        if(time_in_2.ge.3300.0)THEN
           gamma = 1.0-(time_in_2-3300.0)/(3600.0-3300.0)
        endif
        do k=1,nzp
           do j=3,nyp
              do i=3,nxp
                 beta  = (xt(i)-0.)/10000.0
                 delta = (zt(k)-0.)/10000.0
                 if((abs(beta).lt.1.0).and.(abs(delta).lt.1.0)) then
                    a_ut(k,i,j)=a_ut(k,i,j)+alpha*gamma*cos(0.5*pi*beta)*(cosh(2.5*delta))**(-2.)
                 endif
              enddo
           enddo
        enddo
      ENDIF


    endif


  end subroutine forcings
  !
  ! -------------------------------------------------------------------
  ! subroutine smoke_rad:  call simple radiative parameterization for
  ! the smoke cloud
  !
  subroutine smoke_rad(n1,n2,n3,dn0,flx,zm,dzi_t,tt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: zm(n1),dzi_t(n1),dn0(n1),rt(n1,n2,n3)
    real, intent (inout):: tt(n1,n2,n3)
    real, intent (out)  :: flx(n1,n2,n3)
    real, parameter     :: xka= 50.0, fr0=60.0

    integer :: i,j,k, km1, ki
    real    :: smoke(n2,n3)

    smoke=0.
    do j=3,n3-2
       do i=3,n2-2
          ki = n1
          do k=1,n1
             km1=max(1,k-1)
             smoke(i,j)=smoke(i,j)+max(0.,rt(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
          enddo

          do k=2,n1
             km1=max(2,k-1)
             smoke(i,j)=smoke(i,j)-max(0.,rt(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=fr0*exp(-1.*xka*smoke(i,j))
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzi_t(k)/(dn0(k)*cp)
          enddo
       enddo
    enddo

  end subroutine smoke_rad
  !
  ! -------------------------------------------------------------------
  ! subroutine case_forcing: adjusts tendencies according to a specified
  ! large scale forcing.  Normally case (run) specific.
  !
  subroutine case_forcing(n1,n2,n3,case_name,zt,dzi_t,dzi_m,tl,rt,tt,rtt)

    use mpi_interface, only : pecount, double_scalar_par_sum,myid, appl_abort
    use stat, only : get_zi

    integer, intent (in):: n1,n2, n3
    real, dimension (n1), intent (in)          :: zt, dzi_t, dzi_m
    real, dimension (n1,n2,n3), intent (in)    :: tl, rt
    real, dimension (n1,n2,n3), intent (inout) :: tt, rtt
    character (len=5), intent (in) :: case_name

    integer :: i,j,k,kp1
    real, dimension (n1) :: sf
    real, parameter :: zmx_sub = 2260. ! originally 2260.

    real (kind=8) :: zig, zil
    real          :: zibar

    select case (trim(case_name))
    case('rico')

       !
       ! calculate subsidence factor (wsub / dz)
       !
       do k=2,n1-2
          if (zt(k) < zmx_sub) then
             sf(k) =  -0.005*zt(k)/zmx_sub
          else
             sf(k) =  -0.005
          end if
          sf(k) = sf(k)*dzi_t(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! subsidence
                !
                kp1 = k+1
                tt(k,i,j)  =  tt(k,i,j) - ( tl(kp1,i,j) - tl(k,i,j) )*sf(k)
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                !
                ! temperature advection and radiative cooling
                !
                tt(k,i,j) = tt(k,i,j)  - 2.5/86400.
                !
                ! moisture advection
                !
                if (zt(k) <= 2980.) then
                   rtt(k,i,j) = rtt(k,i,j)  - (1. -  1.3456*zt(k)/2980.)/8.64e7
                else
                   rtt(k,i,j) = rtt(k,i,j)  + .3456/8.64e7
                end if
             enddo
          enddo
       enddo

    case ('bomex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       do k=2,n1-2
          if (zt(k) < 1500.) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065  + 0.0065*(zt(k)-1500.)/600.)
          end if
          sf(k) = sf(k)*dzi_t(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                if (zt(k) < 1500.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5
                else if (zt(k) < 2000.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1.- (zt(k)-1500.)*1.e-3)
                end if
                !
                ! moisture advection
                !
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                if (zt(k) < 300.) then
                   rtt(k,i,j) = rtt(k,i,j)  - 1.2e-8
                elseif (zt(k) < 500.) then
                   rtt(k,i,j) = rtt(k,i,j)  - 1.2e-8*(1.- (zt(k)-300.)/200.)
                end if
             enddo
          enddo
       enddo
    case ('atex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       zil = get_zi (n1, n2, n3, 2, rt, dzi_m, zt, 6.5e-3)
       call double_scalar_par_sum(zil,zig)
       zibar = real(zig/pecount)

       do k=2,n1-2
          if (zt(k) < zibar) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065*(1 - (zt(k)-zibar)/300.))
          end if
          sf(k) = sf(k)*dzi_t(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                if (zt(k) < zibar) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1. + (1.- zt(k)/zibar)/2.)
                else if (zt(k) < zibar+300.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1.- (zt(k)-zibar)/300.)
                end if
                !
                ! moisture advection
                !
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                if (zt(k) < zibar) rtt(k,i,j) = rtt(k,i,j)  - 1.5e-8
             enddo
          enddo
       enddo
    case default
       if (myid == 0) print *, '  ABORTING: inproper call to radiation'
       call appl_abort(0)
    end select

  end subroutine case_forcing
  !
  ! -------------------------------------------------------------------
  ! subroutine bellon_rad:  call simple radiative parameterization
  !
  subroutine bellon(n1,n2,n3,flx,sflx,zt,dzi_t,dzi_m,tt,tl,rtt,rt, ut,u,vt,v)

    integer, intent (in) :: n1,n2, n3

    real, dimension (n1), intent (in)            :: zt, dzi_t, dzi_m
    real, dimension (n1, n2, n3), intent (inout) :: tt, tl, rtt, rt, ut,u,vt,v
    real,  dimension (n1, n2, n3), intent (out)  :: flx, sflx
    real, parameter      :: w0= 7.5e-3, H=1000., Qrate = 2.5/86400.

    integer :: i,j,k,kp1
    real    :: grad,wk

    do j=3,n3-2
       do i=3,n2-2
          !
          ! subsidence
          !
          flx(1,i,j)  = 0.
          sflx(1,i,j) = 0.
          do k=2,n1-2
             kp1 = k+1
             wk = w0*(1.-exp(-zt(k)/H))
             grad = Qrate/wk
             flx(k,i,j)  = wk*((tl(kp1,i,j)-tl(k,i,j))*dzi_t(k)-grad)
             sflx(k,i,j) = wk*((rt(kp1,i,j)-rt(k,i,j))*dzi_t(k)-grad)
             tt(k,i,j) = tt(k,i,j) + flx(k,i,j)
             rtt(k,i,j)=rtt(k,i,j) + &
                  wk*(rt(kp1,i,j)-rt(k,i,j))*dzi_t(k)
             ut(k,i,j) =  ut(k,i,j) + &
                  wk*(u(kp1,i,j)-u(k,i,j))*dzi_m(k)
             vt(k,i,j) =  vt(k,i,j) + &
                  wk*(v(kp1,i,j)-v(k,i,j))*dzi_m(k)
          end do
          flx(n1,  i,j)  = 0.
          flx(n1-1,i,j)  = 0.
          sflx(n1,  i,j) = 0.
          sflx(n1-1,i,j) = 0.
       enddo
    enddo

  end subroutine bellon


end module forc
