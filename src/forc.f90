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

  use defs, only      : cp
  use radiation, only : d4stream
  implicit none

  character (len=5), parameter :: case_name = 'xxxx'
    
contains
  !
  ! -------------------------------------------------------------------
  ! subroutine forcings:  calls the appropriate large-scale forcings
  !
  subroutine forcings(time_in, cntlat, sst)

    use grid, only: nxp, nyp, nzp, zm, zt, dzt, dzm, dn0, iradtyp, liquid  &
         , a_rflx, a_sflx, albedo, a_tt, a_tp, a_rt, a_rp, a_pexnr, a_scr1 &
         , vapor, a_rpp, CCN, pi0, pi1, level, a_ut, a_up, a_vt, a_vp,a_theta

    use mpi_interface, only : myid, appl_abort

    real, optional, intent (in) :: time_in, cntlat, sst

    real :: xka, fr0, fr1, div, xref1, xref2

    if (trim(case_name) == 'atex') then
       xka = 130.
       fr0 = 74.
       fr1 = 0.
       div = 0.
    else
       xka = 85.
       fr0 = 70.
       fr1 = 22.
       div = 3.75e-6
    end if

    select case(iradtyp)
    case (1)
        call case_forcing(nzp,nxp,nyp,zt,dzt,dzm,a_tp,a_rp,a_tt,a_rt)
    case (2)
       select case(level)
       case(1) 
          call smoke_rad(nzp, nxp, nyp, dn0, a_rflx, zm, dzt,a_tt,a_rp)
       case(2)
          call gcss_rad(nzp, nxp, nyp, xka, fr0, fr1, div, liquid, dn0,   &
               a_rflx, zt, zm, dzt, a_tt, a_tp, a_rt, a_rp)
!irina add case 3
       case(3)
          call gcss_rad(nzp, nxp, nyp, xka, fr0, fr1, div, liquid, dn0,   &
               a_rflx, zt, zm, dzt, a_tt, a_tp, a_rt, a_rp)
       end select
!
       if (trim(case_name) == 'atex') call case_forcing(nzp, nxp, nyp,    &
            zt, dzt, dzm, a_tp, a_rp, a_tt, a_rt)
    case (3)
       if (present(time_in) .and. present(cntlat) .and. present(sst)) then
          a_scr1 = a_theta/a_pexnr 
          if (level == 3) then
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,   &
                  dn0, pi0, pi1, dzt, a_pexnr, a_scr1, vapor, liquid, a_tt,&
                  a_rflx, a_sflx, albedo, rr=a_rpp)
          else
             xref1 = 0.
             xref2 = 0.
             call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,    &
                  dn0, pi0, pi1, dzt, a_pexnr, a_scr1, vapor, liquid, a_tt, &
                  a_rflx, a_sflx, albedo)
             xref1 = xref1 + a_sflx(nzp,3,3)/albedo(3,3)
             xref2 = xref2 + a_sflx(nzp,3,3)
             albedo(3,3) = xref2/xref1
          end if
       else
          if (myid == 0) print *, '  ABORTING: inproper call to radiation'
          call appl_abort(0)
       end if
    case (4)
       call bellon(nzp, nxp, nyp, a_rflx, a_sflx, zt, dzt, dzm, a_tt, a_tp&
            ,a_rt, a_rp, a_ut, a_up, a_vt, a_vp)
    end select 

  end subroutine forcings
  !
  ! -------------------------------------------------------------------
  ! subroutine gcss_rad:  call simple radiative parameterization and 
  ! simultaneously update fields due to vertical motion as given by div
  !
  subroutine gcss_rad(n1,n2,n3,xka,fr0,fr1,div,rc,dn0,flx,zt,zm,dzt,   &
       tt,tl,rtt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: xka, fr0, fr1, div
    real, intent (in)   :: zt(n1),zm(n1),dzt(n1),dn0(n1),rc(n1,n2,n3),   &
         tl(n1,n2,n3),rt(n1,n2,n3)
    real, intent (inout):: tt(n1,n2,n3),rtt(n1,n2,n3)
    real, intent (out)  :: flx(n1,n2,n3)

    integer :: i, j, k, km1, kp1, ki
    real    :: lwp(n2,n3), fact

    lwp=0.
    do j=3,n3-2
       do i=3,n2-2
          ki = n1
          do k=1,n1
             km1=max(1,k-1)
             lwp(i,j)=lwp(i,j)+max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
             flx(k,i,j)=fr1*exp(-1.*xka*lwp(i,j))
             if ( (rc(k,i,j) > 0.01e-3) .and. (rt(k,i,j) >= 0.008) ) ki=k
          enddo

          fact = div*cp*dn0(ki)
          do k=2,n1
             km1=max(2,k-1)
             lwp(i,j)=lwp(i,j)-max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=flx(k,i,j)+fr0*exp(-1.*xka*lwp(i,j))
             if (zm(k) > zm(ki) .and. ki > 1 .and. fact > 0.) then
                flx(k,i,j)=flx(k,i,j) + fact*(0.25*(zm(k)-zm(ki))**1.333 + &
                  zm(ki)*(zm(k)-zm(ki))**0.333333)
             end if
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzt(k)/(dn0(k)*cp)
          enddo
          !
          ! subsidence
          !
          if (div /= 0.) then
             do k=2,n1-2
                kp1 = k+1
                tt(k,i,j) = tt(k,i,j) + &
                        div*zt(k)*(tl(kp1,i,j)-tl(k,i,j))*dzt(k)
                rtt(k,i,j)=rtt(k,i,j) + &
                        div*zt(k)*(rt(kp1,i,j)-rt(k,i,j))*dzt(k)
             end do
          end if
       enddo
    enddo

  end subroutine gcss_rad
  !
  ! -------------------------------------------------------------------
  ! subroutine smoke_rad:  call simple radiative parameterization for 
  ! the smoke cloud
  !
  subroutine smoke_rad(n1,n2,n3,dn0,flx,zm,dzt,tt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: zm(n1),dzt(n1),dn0(n1),rt(n1,n2,n3)
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
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzt(k)/(dn0(k)*cp)
          enddo
       enddo
    enddo

  end subroutine smoke_rad
  !
  ! -------------------------------------------------------------------
  ! subroutine case_forcing: adjusts tendencies according to a specified
  ! large scale forcing.  Normally case (run) specific.
  !
  subroutine case_forcing(n1,n2,n3,zt,dzt,dzm,tl,rt,tt,rtt)

    use mpi_interface, only : pecount, double_scalar_par_sum,myid, appl_abort
    use stat, only : get_zi

    integer, intent (in):: n1,n2, n3
    real, dimension (n1), intent (in)          :: zt, dzt, dzm
    real, dimension (n1,n2,n3), intent (in)    :: tl, rt
    real, dimension (n1,n2,n3), intent (inout) :: tt, rtt

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
          sf(k) = sf(k)*dzt(k)
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
          sf(k) = sf(k)*dzt(k)
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
       zil = get_zi (n1, n2, n3, 2, rt, dzm, zt, 6.5e-3)
       call double_scalar_par_sum(zil,zig)
       zibar = real(zig/pecount)

       do k=2,n1-2
          if (zt(k) < zibar) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065*(1 - (zt(k)-zibar)/300.))
          end if
          sf(k) = sf(k)*dzt(k)
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
  subroutine bellon(n1,n2,n3,flx,sflx,zt,dzt,dzm,tt,tl,rtt,rt, ut,u,vt,v)

    integer, intent (in) :: n1,n2, n3

    real, dimension (n1), intent (in)            :: zt, dzt, dzm
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
             flx(k,i,j)  = wk*((tl(kp1,i,j)-tl(k,i,j))*dzt(k)-grad)
             sflx(k,i,j) = wk*((rt(kp1,i,j)-rt(k,i,j))*dzt(k)-grad)
             tt(k,i,j) = tt(k,i,j) + flx(k,i,j)
             rtt(k,i,j)=rtt(k,i,j) + &
                  wk*(rt(kp1,i,j)-rt(k,i,j))*dzt(k)
             ut(k,i,j) =  ut(k,i,j) + &
                  wk*(u(kp1,i,j)-u(k,i,j))*dzm(k)
             vt(k,i,j) =  vt(k,i,j) + &
                  wk*(v(kp1,i,j)-v(k,i,j))*dzm(k)
          end do
          flx(n1,  i,j)  = 0.
          flx(n1-1,i,j)  = 0.
          sflx(n1,  i,j) = 0.
          sflx(n1-1,i,j) = 0.
       enddo
    enddo

  end subroutine bellon
 
end module forc
