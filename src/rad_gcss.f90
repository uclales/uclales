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
module rad_gcss

  use defs, only      : cp,pi
  implicit none
!
contains
  !
  ! -------------------------------------------------------------------
  ! subroutine gcss_rad:  call simple radiative parameterization for the LW and
  ! SW fluxes (LW similar to DYCOMS, but without the 3rd term, SW Delta Eddington EUROCS)
  ! case and simultaneously update fields due to vertical motion as given by div

  subroutine gcss_rad(n1,n2,n3,alat,time,case_name,div,sst, rc,dn0,flx,swn,zt,zm,dzi_t,   &
       tt,tl,rtt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: div, sst, alat, time
    real, intent (in)   :: zt(n1),zm(n1),dzi_t(n1),dn0(n1),rc(n1,n2,n3),   &
         tl(n1,n2,n3),rt(n1,n2,n3)
    real, intent (inout):: tt(n1,n2,n3),rtt(n1,n2,n3)
    real, intent (inout)  :: flx(n1,n2,n3),swn(n1,n2,n3)
    character (len=5), intent (in) :: case_name

    integer :: i, j, k, km1, kp1,ki
    real    :: lwp(n2,n3), mu,tauc, tau(n1), xka,fr0,fr1,xkb,fact
!     real,parameter :: dens_air=1.14
    real,parameter :: rho_l=1000.
    real,parameter :: reff=1.e-05


!irina
!default values
       xka = 85.
       fr0 = 70.
       fr1 = 22.

    if (trim(case_name) == 'atex') then
       xka = 130.
       fr0 = 74.
       fr1 = 0.
    else if (trim(case_name) == 'astex' .or. trim(case_name) == 'trans') then
       xka = 130.
       xkb = 80.
       fr0 = 70.
       fr1 = 10.
    end if

!print *, 'uses astex rad'

! determine the solar geometery, as measured by mu, the cosine of the 
! solar zenith angle
!
mu = zenith(alat,time)

!print *, 'mu', mu

   lwp=0.
!print *, 'print 1'
   swn(:,:,:)=0.
!print *, 'print 2'

    do j=3,n3-2
       do i=3,n2-2
        ki = n1
          do k=1,n1
             km1=max(1,k-1)
          if (trim(case_name) == 'astex' .or. trim(case_name) == 'trans') then
             lwp(i,j)=lwp(i,j)+max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
             flx(k,i,j)=fr1*exp(-1.*xka*lwp(i,j))
          else
          !this is for dycoms in fact
             lwp(i,j)=lwp(i,j)+max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
             flx(k,i,j)=fr1*exp(-1.*xka*lwp(i,j))
             if ( (rc(k,i,j) > 0.01e-3) .and. (rt(k,i,j) >= 0.008) ) ki=k
          end if
          enddo

!print *, 'astex rad after lw'

         !SW
         if (mu > 0.035) then  !factor needed for security
          tauc = 0.           ! tau cloud
          do k=1,n1
             km1=max(1,k-1)
             tau(k) = 0.      ! tau dz
             if (rc(k,i,j) > 1e-5) then
             tau(k)=max(0.,1.5*rc(k,i,j)*dn0(k)*(zm(k)-zm(km1))/reff/rho_l)
             tauc=tauc+tau(k)
             end if
           end do  
!print *, 'astex rad before sunray',i,j
          call sunray(mu,tau,tauc,i,j,swn)
!print *, 'astex rad after sunray'
          end if

  fact = div*cp*dn0(ki)

!print *, 'astex rad before lw2'
		flx(1,i,j)=flx(1,i,j)+fr0*exp(-1.*xka*lwp(i,j))
          do k=2,n1
             km1=max(2,k-1)
            lwp(i,j)=lwp(i,j)-max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=flx(k,i,j)+fr0*exp(-1.*xka*lwp(i,j))
          if (trim(case_name) .ne. 'astex' .and. trim(case_name) .ne. 'trans') then
             if (zm(k) > zm(ki) .and. ki > 1 .and. fact > 0.) then
                flx(k,i,j)=flx(k,i,j) + fact*(0.25*(zm(k)-zm(ki))**1.333 + &
                  zm(ki)*(zm(k)-zm(ki))**0.333333)
             end if
             end if
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzi_t(k)/(dn0(k)*cp)
             tt(k,i,j) =tt(k,i,j)+(swn(k,i,j)-swn(km1,i,j))*dzi_t(k)/(dn0(k)*cp)
          enddo

!print *, 'astex rad after lw2'

          !
          ! subsidence
          !
          if (div /= 0.) then
             do k=2,n1-2
                kp1 = k+1
                tt(k,i,j) = tt(k,i,j) + &
                        div*zt(k)*(tl(kp1,i,j)-tl(k,i,j))*dzi_t(k)
                rtt(k,i,j)=rtt(k,i,j) + &
                        div*zt(k)*(rt(kp1,i,j)-rt(k,i,j))*dzi_t(k)
             end do
          end if
       enddo
    enddo

  ! ---------------------------------------------------------------------------
  ! Return the cosine of the solar zenith angle give the decimal day and 
  ! the latitude
  ! 
!print *, 'end astex rad'

  end subroutine gcss_rad
  !

  !                                                                      c
  subroutine sunray(mu,tau,tauc,i,j,swn)

!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   the sunray model is described by fouquart and  bonnel                  c
!   (1980, contr. atmos. phys.).                                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   use grid, only: nxp, nyp, nzp, sfc_albedo

  real, intent(inout), dimension (nzp) :: tau
  real, intent(out), dimension (nzp,nxp,nyp) :: swn
  integer, intent(in) :: i,j
  real, intent(in) :: mu
  real,allocatable, dimension (:) :: taude
  real :: gcde,tauc,taucde &
           ,taupath,t1,t2,t3,c1,c2 &
           ,omega,omegade,ff,x1,x2,x3,rk,mu2,rp,alpha,beta,rtt &
           ,exmu0,expk,exmk,xp23p,xm23p,ap23b
  integer :: k

  real :: sw0        = 1100.0 !direct component at top of the cloud (W/m^2), diffuse not possible
  real :: gc         = 0.85   !asymmetry factor of droplet scattering angle distribution
  !real :: sfc_albedo = 0.05   !ground surface albedo

  allocate(taude(nzp))

  taucde = 0.         ! tau' cloud
  taupath = 0.
  do k=1,nzp
     taude(k) =0.     ! tau' laagje dz
  end do

! omega=0.9989-4.e-3*exp(-0.15*tauc)  !   fouquart and bonnel (1980)
  omega=1.-1.e-3*(0.9+2.75*(mu+1.)*exp(-0.09*tauc)) !fouquart

! the equations for the delta-eddington approximation are equal to those
! for the eddington approximation with transformed parameters g, omega
! and tau (joseph, wiscomb and weinman, 1976, j.a.s.).
! parameternames: x -> xde (delta-eddington)

  ff=gc*gc
  gcde=gc/(1.+gc)
  taucde=(1.0-omega*ff)*tauc

  do k =1,nzp
    taude(k)=(1.e0-omega*ff)*tau(k)
  end do

  omegade=(1.0-ff)*omega/(1.e0-omega*ff)

! the solution of the eddington equations are given by shettle and weinman
! (1970, j.a.s.).

  x1=1.0-omegade*gcde
  x2=1.0-omegade
  rk=sqrt(3.0*x2*x1)
  mu2=mu*mu
  x3=4.0*(1.0-rk*rk*mu2)
  rp=sqrt(3.0*x2/x1)
  alpha=3.e0*omegade*mu2*(1.0+gcde*x2)/x3
  beta=3.*omegade*mu*(1.0+3.0*gcde*mu2*x2)/x3

  rtt=2.0/3.0
  exmu0= exp(-taucde/mu)
  expk=  exp(rk*taucde)
  exmk=1.0/expk
  xp23p=1.0+rtt*rp
  xm23p=1.0-rtt*rp
  ap23b=alpha+rtt*beta

  t1=1-sfc_albedo-rtt*(1.+sfc_albedo)*rp
  t2=1-sfc_albedo+rtt*(1.+sfc_albedo)*rp
  t3=(1-sfc_albedo)*alpha-rtt*(1+sfc_albedo)*beta+sfc_albedo*mu
  c2=(xp23p*t3*exmu0-t1*ap23b*exmk)/(xp23p*t2*expk-xm23p*t1*exmk)
  c1=(ap23b-c2*xm23p)/xp23p

  do k = nzp,1,-1
      taupath = taupath + taude(k)
      swn(k,i,j)=sw0*(4./3.)*(rp*(c1*exp(-rk*taupath) &
                 -c2*exp(rk*taupath)) &
                 -beta*exp(-taupath/mu)) &
                 +mu*sw0*exp(-taupath/mu)
  end do

  deallocate(taude)

!  return
  end subroutine sunray

!
  real function zenith(alat,time)

    real, intent (in)  :: alat, time

    real :: lamda, d, sig, del, h, day

    day    = floor(time)
    lamda  = alat*pi/180.
    d      = 2.*pi*int(time)/365.
    sig    = d + pi/180.*(279.9340 + 1.914827*sin(d) - 0.7952*cos(d) &
         &                      + 0.019938*sin(2.*d) - 0.00162*cos(2.*d))
    del    = asin(sin(23.4439*pi/180.)*sin(sig))
    h      = 2.*pi*((time-day)-0.5)
    zenith = sin(lamda)*sin(del) + cos(lamda)*cos(del)*cos(h)

  end function zenith

end module rad_gcss
