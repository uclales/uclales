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
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!

module thrm

  implicit none

contains
!
! -------------------------------------------------------------------------
! THERMO: calculates thermodynamics quantities according to level.  Level
! is passed in to allow level of diagnosis to be determined by call rather
! than by runtype
!
  subroutine thermo (level,opt)

    use grid, only : liquid, vapor, a_theta, a_pexnr, press, a_scr1,  &
         a_scr2, a_rp, a_tp, nxp, nyp, nzp, th00, pi0, pi1,a_rpp,rsi,a_ricep,a_rsnowp,a_rgrp, &
         cnd_acc, cev_acc, zm, dn0, nstep, a_cld
    integer, intent (in) :: level
    integer, intent (in), optional :: opt

    select case (level)
      case default
         call drythrm(nzp,nxp,nyp,a_pexnr,press,a_tp,a_theta,a_scr1,pi0,   &
              pi1,th00)
      case (1)
         call drythrm(nzp,nxp,nyp,a_pexnr,press,a_tp,a_theta,a_scr1,pi0,   &
              pi1,th00,a_rp,vapor)
      case (2,3)
         if (present(opt)) then ! for lwaterbudget = .true.
            call satadjst(level,nzp,nxp,nyp,a_pexnr,press,a_tp,a_theta,a_scr1,pi0,  &
                 pi1,th00,a_rp,vapor,liquid,a_scr2, &
                 cloud=a_cld,tcond=cnd_acc,tevap=cev_acc,zm=zm,rho=dn0,nstep=nstep)
         else
            call satadjst(level,nzp,nxp,nyp,a_pexnr,press,a_tp,a_theta,a_scr1,pi0,  &
                 pi1,th00,a_rp,vapor,liquid,a_scr2)
         end if
      case (4,5)
         call satadjst(level,nzp,nxp,nyp,a_pexnr,press,a_tp,a_theta,a_scr1,pi0,  &
              pi1,th00,a_rp,vapor,liquid,a_scr2,rsi=rsi)
    end select

  end subroutine thermo

!
! -------------------------------------------------------------------------
! update_pi1:  this routine updates a pressure associated with the 
! subtraction of a mean acceleration, only incrementing it for dynamic and 
! thermal effects for layers above the surface
! 
  subroutine update_pi1(n1,awtbar,pi1)

    use grid, only : th00, zt

    integer, intent (in) :: n1
    real, intent (in) , dimension (n1) :: awtbar
    real, intent (inout), dimension (n1) :: pi1

    integer :: k

    do k=2,n1
       pi1(k) = pi1(k-1) + awtbar(k-1)*(zt(k)-zt(k-1))/th00 
    end do

  end subroutine update_pi1

!
! -------------------------------------------------------------------------
! DRYTHRM:  this routine calculates theta, pressure and optionally
! moisture (vapor phase only)
! 
  subroutine drythrm(n1,n2,n3,pp,p,thil,theta,t,pi0,pi1,th00,rt,rv)

  use defs, only : cp, cpr, p00

  integer, intent (in) :: n1,n2,n3
  real, intent (in)    :: pi0(n1),pi1(n1),th00
  real, intent (in)    :: pp(n1,n2,n3),thil(n1,n2,n3)
  real, intent (out)   :: p(n1,n2,n3),theta(n1,n2,n3),t(n1,n2,n3)
  real, intent (in), optional  :: rt(n1,n2,n3)
  real, intent (out), optional :: rv(n1,n2,n3)
  !logical, intent(in)  :: moist

  integer :: k,i,j
  real    :: exner

  do j=3,n3-2
    do i=3,n2-2
      do k=1,n1
        exner        = (pi0(k)+pi1(k)+pp(k,i,j))/cp
        p(k,i,j)     = p00 * (exner)**cpr
        theta(k,i,j) = thil(k,i,j)+th00
        t(k,i,j)     = theta(k,i,j)*exner
        if(present(rv)) rv(k,i,j) = rt(k,i,j)
      enddo
    enddo
  enddo

  end subroutine drythrm



! 
! -------------------------------------------------------------------------
! SATADJST:  this routine calculates theta, and pressure and diagnoses
! liquid water using a saturation adjustment for warm-phase systems
!
  subroutine satadjst(level,n1,n2,n3,pp,p,tl,th,tk,pi0,pi1,th00,rt,rv,rc,rs,rsi,&
       cloud,tcond,tevap,zm,rho,nstep)

    use defs, only : cp, cpr, alvl, ep, Rm, p00,t_hn,tmelt
    use mpi_interface, only : appl_abort

    integer, intent (in) ::  n1,n2,n3,level

    real, intent (in), dimension (n1,n2,n3)    :: pp, tl, rt
    real, intent (in), dimension (n1)          :: pi0, pi1
    real, intent (in)                          :: th00
    real, intent (out), dimension (n1,n2,n3)   :: rv,rs,th,tk,p
    real, intent (out), dimension (n1,n2,n3) :: rc
    real, intent (out), optional, dimension (n1,n2,n3) :: rsi
    real, intent (inout), optional, dimension (n1,n2,n3) :: cloud
    real, intent (inout), optional, dimension (n2,n3)    :: tcond,tevap
    real, intent (in), optional, dimension (n1)  :: zm,rho
    integer, intent (in), optional :: nstep

    integer :: k, i, j, iterate
    real    :: exner,tli,txi,tx1,tx,rsx,rix,rcx,ravail,dtx,part
    real, parameter :: epsln = 1.e-4

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             dtx = 2*epsln
             iterate = 1
             exner = (pi0(k)+pi1(k)+pp(k,i,j))/cp
             p(k,i,j) = p00 * (exner)**cpr
             tli=(tl(k,i,j)+th00)*exner
             tx=tli
             rsx=rslf(p(k,i,j),tx)
             if (level>3) rix=rsif(p(k,i,j),tx)
             rcx=max(rt(k,i,j)-rsx,0.)
             ravail = rt(k,i,j)
             part = 1.
             if (rcx > 0.) then
                do while (dtx > epsln .and. iterate < 20)
                   txi=alvl/(cp*tx)
                   tx1=tx - (tx - tli*(1.+txi*rcx))/(1. + txi*tli                &
                        *(rcx/tx+(1.+rsx*ep)*rsx*alvl/(Rm*tx*tx)))
                   dtx = abs(tx1-tx)
                   tx  = tx1
                   iterate = iterate+1
                   rsx=rslf(p(k,i,j),tx)
                   if (level>3) then
                     rix=rsif(p(k,i,j),tx)
! !                     part = max(0.,min(1.,(tx-t_hn)/(tmelt-t_hn)))
                   end if
                   rcx = part*max(ravail-rsx,0.)
                end do

                ! Test BvS for calcrl func
                !print*,'model=',rcx,'single=',calcrl(tl(k,i,j)+th00,rt(k,i,j),p(k,i,j),exner)

                if (dtx > epsln) then
                    print *, k,i,j,p(k,i,j), tli,rcx  ,exner
                    print *, '  ABORTING: thrm', dtx, epsln
                   call appl_abort(0)
                endif
             endif
             if (present(tcond)) then
                if (k.gt.1) then
                   ! diagnostic to get precipitation efficiency
                   tcond(i,j) = tcond(i,j) + max(rcx-cloud(k,i,j),0.)*(zm(k)-zm(k-1))*rho(k)
                   tevap(i,j) = tevap(i,j) + max(cloud(k,i,j)-rcx,0.)*(zm(k)-zm(k-1))*rho(k)
                end if
                cloud(k,i,j) = rcx
             end if
             rc(k,i,j)=rcx
             rv(k,i,j)=rt(k,i,j)-rc(k,i,j)
             rs(k,i,j)  = rsx
             if (level>3) rsi(k,i,j) = rix
             tk(k,i,j)=tx
             th(k,i,j)=tk(k,i,j)/exner
          enddo
       enddo
    enddo


  end subroutine satadjst
! 
! 
!   subroutine satpart(n1,n2,n3,rv,rc,rsup,tk)
!     use defs, only : t_hn, tmelt
!     integer, intent(in) :: n1,n2,n3
!     real,dimension(n1,n2,n3), intent(inout) :: rv,rc
!     real,dimension(n1,n2,n3), intent(out)   :: rsup
!     real,dimension(n1,n2,n3), intent(in)    :: tk
!     integer :: i,j,k
!     real :: part
!     
!      do j=3,n3-2
!        do i=3,n2-2
!           do k=1,n1
!             part = max(0.,min(1.,(tk(k,i,j)-t_hn)/(tmelt-t_hn)))
!             rsup(k,i,j) = rc(k,i,j)*(1-part)
!             rc(k,i,j)   = part*rc(k,i,j)
!             rv(k,i,j)   = rsup(k,i,j)+rv(k,i,j)
!         end do
!       end do
!     end do
! 
!   end subroutine satpart
! 

! ---------------------------------------------------------------------
! This function calculates the water saturation vapor mixing ratio as a
! function of temperature and pressure
!
  real elemental function rslf(p,t)

    real, intent (in) :: p, t
    rslf=.622*esl(t)/(p-esl(t))

  end function rslf
!
! ---------------------------------------------------------------------
!
  real elemental function esl(t)

    use defs, only : tmelt
    real, intent (in) :: t
    real, parameter :: c0=0.6105851e+03, c1=0.4440316e+02,    &
                       c2=0.1430341e+01, c3=0.2641412e-01,    &
                       c4=0.2995057e-03, c5=0.2031998e-05,    &
                       c6=0.6936113e-08, c7=0.2564861e-11,    &
                       c8=-.3704404e-13
    real :: x
    x   = max(-80.,t-tmelt)
    esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

  end function esl
!
! ---------------------------------------------------------------------
! This function calculates the ice saturation vapor mixing ratio as a
! function of temperature and pressure
!
  real elemental function rsif(p,t)

    real, intent (in) :: p, t
    rsif = .622*esi(t)/(p-esi(t))

  end function rsif
!
! ---------------------------------------------------------------------
!
  real elemental function esi(t)

    use defs, only : tmelt
    real, intent (in) :: t
    real, parameter   :: c0_i=0.6114327e+03, c1_i=0.5027041e+02,    &
                         c2_i=0.1875982e+01, c3_i=0.4158303e-01,    &
                         c4_i=0.5992408e-03, c5_i=0.5743775e-05,    &
                         c6_i=0.3566847e-07, c7_i=0.1306802e-09,    &
                         c8_i=0.2152144e-12
    real  :: x
    x   = max(-80.,t-tmelt)
    esi = c0_i+x*(c1_i+x*(c2_i+x*(c3_i+x*(c4_i+x*(c5_i+x*(c6_i+x*(c7_i+x*c8_i)))))))

  end function esi
!
! -------------------------------------------------------------------------
! Calculate liquid water (satadjust) for single grid point
!
  real elemental function calcrl(thl,rt,p,exner)
    
    use defs, only    : alvl,cp,Rm,ep

    real, intent(in) :: thl, rt, p, exner
    real             :: rl, tl, tx_old, tx, rs, txi 
    integer          :: niter 

    tl     = thl * exner   ! Liquid water temperature
    tx     = thl           ! Absolute temp to solve for
    niter  = 0         
    tx_old = 1.e9
   
    do while ((abs(tx-tx_old) .gt. 1e-4) .and. niter < 20)
      tx_old = tx
      rs     = rslf(p,tx)        ! Saturation mixing ratio
      rl     = max(0.,rt-rs)     ! Liquid water mixing ratio
      txi    = alvl / (cp * tx)
      tx     = tx-(tx-tl*(1.+txi*rl))/(1.+txi*tl*(rl/tx+(1.+rs*ep)*rs*alvl/(Rm*tx*tx)))
      niter  = niter + 1  
    end do

    calcrl = max(0.,rt-rslf(p,tx))

  end function calcrl
 
! -------------------------------------------------------------------------
! FLL_TKRS: Updates scratch arrays with temperature and saturation mixing
! ratio
! 
  subroutine fll_tkrs(n1,n2,n3,th,pp,pi0,pi1,tk,rs)

  use defs, only : cp, R,cpr,p00

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: th(n1,n2,n3), pp(n1,n2,n3)
    real, intent (in)    :: pi0(n1), pi1(n1)
    real, intent (out)   :: tk(n1,n2,n3)
    real, optional, intent (out)   :: rs(n1,n2,n3)

    integer :: i, j, k
    real    :: exner

    do j=3,n3-2
      do i=3,n2-2
        do k=1,n1
          exner=(pi0(k)+pi1(k)+pp(k,i,j))/cp
          tk(k,i,j)=th(k,i,j)*exner
          if (present(rs)) rs(k,i,j)=rslf(p00*(exner)**cpr,tk(k,i,j))
        end do
      end do
    end do

  end subroutine fll_tkrs
! 
! -------------------------------------------------------------------------
! BRUVAIS:  Cacluates the brunt-vaisaila frequency in accordance with the
! thermodynamic level
! 
  subroutine bruvais(n1,n2,n3,level,th,tl,en2,dzi_m,th00,rt,rs)

    use defs, only : g, R, cp, alvl, ep, ep2

    integer, intent (in) ::  n1, n2, n3, level
    real, intent (in)    ::  th(n1,n2,n3), tl(n1,n2,n3),dzi_m(n1), th00
    real, intent (in), optional :: rt(n1,n2,n3), rs(n1,n2,n3)
    real, intent (out)   ::  en2(n1,n2,n3)

    integer :: i, k, j, kp1
    real    :: c1, c2, c3, tvk, tvkp1, rtbar, rsbar, aa, bb

    select case(level)
    case(0)
      do j=3,n3-2
         do i=3,n2-2
            ! Inner loop vectorized
            do k=1,n1-1
               c1=(1.+ep*alvl/R/th(k,i,j))/ep
               c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
               c3=alvl/(cp*th(k,i,j))
               en2(k,i,j)=g*dzi_m(k)*((th(k+1,i,j)-th(k,i,j))/th00)
            end do
            en2(n1,i,j)=en2(n1-1,i,j)
         end do
      end do
    case(1)
      do j=3,n3-2
         do i=3,n2-2
            ! Inner loop vectorized
            do k=1,n1-1
               c1=(1.+ep*alvl/R/th(k,i,j))/ep
               c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
               c3=alvl/(cp*th(k,i,j))
               tvk=th(k,i,j)*(1.+ep2*rt(k,i,j))
               tvkp1=th(k+1,i,j)*(1.+ep2*rt(k+1,i,j))
               en2(k,i,j)=g*dzi_m(k)*(tvkp1-tvk)/th00
            end do
            en2(n1,i,j)=en2(n1-1,i,j)
         end do
      end do
    case(2)
      do j=3,n3-2
         do i=3,n2-2
            ! Not vectorized -> if statement...
            do k=1,n1-1
               c1=(1.+ep*alvl/R/th(k,i,j))/ep
               c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
               c3=alvl/(cp*th(k,i,j))
               rtbar=0.5*(rt(k,i,j)+rt(k+1,i,j))
               rsbar=0.5*(rs(k,i,j)+rs(k+1,i,j))
               kp1=min(n1-1,k+1)
               if (rt(k,i,j) > rs(k,i,j) .and. rt(kp1,i,j) > rs(kp1,i,j)) then
                  aa=(1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
                  bb=(c3*aa - 1.)
               else
                  aa=(1.00 + ep2*rtbar)
                  bb=ep2
               end if
               en2(k,i,j)=g*dzi_m(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
                    + bb*(rt(k+1,i,j)-rt(k,i,j)))
            end do
            en2(n1,i,j)=en2(n1-1,i,j)
         end do
      end do
    case(3,4,5)
      do j=3,n3-2
         do i=3,n2-2
            ! Not vectorized -> if statement...
            do k=1,n1-1
               c1=(1.+ep*alvl/R/th(k,i,j))/ep
               c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
               c3=alvl/(cp*th(k,i,j))
               rtbar=0.5*(rt(k,i,j)+rt(k+1,i,j))
               rsbar=0.5*(rs(k,i,j)+rs(k+1,i,j))
               kp1=min(n1-1,k+2)
               if (rt(k,i,j) > rs(k,i,j) .and. rt(kp1,i,j) > rs(kp1,i,j)) then
                  aa=(1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
                  bb=(c3*aa - 1.)
               else
                  aa=(1.00 + ep2*rtbar)
                  bb=ep2
               end if
               en2(k,i,j)=g*dzi_m(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
                    + bb*(rt(k+1,i,j)-rt(k,i,j)))
            end do
            en2(n1,i,j)=en2(n1-1,i,j)
         end do
      end do
    case default
      WRITE (*,*) 'level=',level,', not supported in bruvais'
      stop 
    end select

    ! BvS: code doesn't vectorize with select statement in loop
    ! With vectorization, level 0,1 factor 15 faster 
    !do j=3,n3-2
    !   do i=3,n2-2
    !      do k=1,n1-1
    !         c1=(1.+ep*alvl/R/th(k,i,j))/ep
    !         c2=ep*alvl*alvl/(R*cp*th(k,i,j)*th(k,i,j))
    !         c3=alvl/(cp*th(k,i,j))
    !         select case(level) 
    !         case (0)
    !            en2(k,i,j)=g*dzi_m(k)*((th(k+1,i,j)-th(k,i,j))/th00)
    !         case (1)
    !            tvk=th(k,i,j)*(1.+ep2*rt(k,i,j))
    !            tvkp1=th(k+1,i,j)*(1.+ep2*rt(k+1,i,j))
    !            en2(k,i,j)=g*dzi_m(k)*(tvkp1-tvk)/th00
    !         case (2)
    !            rtbar=0.5*(rt(k,i,j)+rt(k+1,i,j))
    !            rsbar=0.5*(rs(k,i,j)+rs(k+1,i,j))
    !            kp1=min(n1-1,k+1)
    !            if (rt(k,i,j) > rs(k,i,j) .and. rt(kp1,i,j) > rs(kp1,i,j)) then
    !               aa=(1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
    !               bb=(c3*aa - 1.)
    !            else
    !               aa=(1.00 + ep2*rtbar)
    !               bb=ep2
    !            end if
    !            en2(k,i,j)=g*dzi_m(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
    !                 + bb*(rt(k+1,i,j)-rt(k,i,j)))
    !         case (3,4,5)
    !            rtbar=0.5*(rt(k,i,j)+rt(k+1,i,j))
    !            rsbar=0.5*(rs(k,i,j)+rs(k+1,i,j))
    !            kp1=min(n1-1,k+2)
    !            if (rt(k,i,j) > rs(k,i,j) .and. rt(kp1,i,j) > rs(kp1,i,j)) then
    !               aa=(1. - rtbar + rsbar*c1)/(1. + c2*rsbar)
    !               bb=(c3*aa - 1.)
    !            else
    !               aa=(1.00 + ep2*rtbar)
    !               bb=ep2
    !            end if
    !            en2(k,i,j)=g*dzi_m(k)*(aa*(tl(k+1,i,j)-tl(k,i,j))/th00        &
    !                 + bb*(rt(k+1,i,j)-rt(k,i,j)))
    !         case default 
    !            WRITE (*,*) 'level=',level,', not supported in bruvais'
    !            stop 
    !         end select
    !      end do
    !      en2(n1,i,j)=en2(n1-1,i,j)
    !   end do
    !end do

  end subroutine bruvais

end module thrm
