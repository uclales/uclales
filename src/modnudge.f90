!> \file modnudge.f90
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudgeT
!>

!>
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudgeT
!>
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modnudge


implicit none
PRIVATE
PUBLIC :: nudge,lnudge,tnudgefac, qfloor, zfloor, znudgemin, znudgeplus, nudge_bound, lnudge_bound
SAVE
  real, dimension(:,:), allocatable :: tnudge,unudge,vnudge,wnudge,thlnudge,qtnudge
  real, dimension(:)  , allocatable :: timenudge
  real :: tnudgefac = 1., qfloor = -1.,  zfloor = 1200., znudgemin = -1., znudgeplus = -1.
  logical :: lnudge,lunudge,lvnudge,lwnudge,lthlnudge,lqtnudge
  integer :: ntnudge = 100
  logical :: firsttime = .true.
! LINDA, b
  ! arrays for nuding
  logical :: lnudge_bound = .false.

  ! arrays for initial values
  real, dimension(:),allocatable::tn,rn,un,vn,wn
!RV
  real, allocatable :: rlz(:) 
  integer kstrt
!rv
  real:: coef=1./60.
  real, allocatable :: rlx(:,:)
! LINDA, e
contains
  subroutine initnudge(time)
    use grid, only : nzp,zt,th00,umean,vmean
    use mpi_interface, only : myid
    use defs, only : pi
    implicit none

    integer :: ierr,k,t,ifinput = 19
    real,allocatable,dimension(:) :: height
    real, intent(in) :: time
    character(1) :: chmess1
    real :: highheight,highqtnudge,highthlnudge,highunudge,highvnudge,highwnudge,hightnudge
    real :: lowheight,lowqtnudge,lowthlnudge,lowunudge,lowvnudge,lowwnudge,lowtnudge
    real :: fac
    allocate(tnudge(nzp,ntnudge),unudge(nzp,ntnudge),vnudge(nzp,ntnudge),wnudge(nzp,ntnudge),thlnudge(nzp,ntnudge),qtnudge(nzp,ntnudge))
    allocate(timenudge(0:ntnudge), height(nzp))
    tnudge = 0
    unudge=0
    vnudge=0
    wnudge=0
    thlnudge=0
    qtnudge=0
    timenudge=0
    height = 0.


    if (.not. lnudge) return
      t = 0
      open (ifinput,file='nudge_in')
      ierr = 0
      readloop: do
        t = t + 1
        chmess1 = "#"
        ierr = 1 ! not zero
        do while (.not.(chmess1 == "#" .and. ierr ==0)) !search for the next line consisting of "# time", from there onwards the profiles will be read
          read(ifinput,*,iostat=ierr) chmess1,timenudge(t)
          if (ierr < 0) exit readloop

        end do
        write(6,*) ' height    t_nudge    u_nudge    v_nudge    w_nudge    thl_nudge    qt_nudge'
        read (ifinput,*)  lowheight , lowtnudge ,  lowunudge , lowvnudge , lowwnudge , lowthlnudge, lowqtnudge
        read (ifinput,*)  highheight , hightnudge ,  highunudge , highvnudge , highwnudge , highthlnudge, highqtnudge
        do  k=2,nzp-1
          if (highheight<zt(k)) then
            lowheight = highheight
            lowtnudge = hightnudge
            lowunudge = highunudge
            lowvnudge = highvnudge
            lowwnudge = highwnudge
            lowthlnudge= highthlnudge
            lowqtnudge=highqtnudge
            read (ifinput,*)  highheight , hightnudge ,  highunudge , highvnudge , highwnudge , highthlnudge, highqtnudge
          end if
          fac = (highheight-zt(k))/(highheight - lowheight)
          tnudge(k,t) = fac*lowtnudge + (1-fac)*hightnudge
          unudge(k,t) = fac*lowunudge + (1-fac)*highunudge
          vnudge(k,t) = fac*lowvnudge + (1-fac)*highvnudge
          wnudge(k,t) = fac*lowwnudge + (1-fac)*highwnudge
          thlnudge(k,t) = fac*lowthlnudge + (1-fac)*highthlnudge
          qtnudge(k,t) = fac*lowqtnudge + (1-fac)*highqtnudge
        end do
        if (myid == 0) then
          do k=nzp-1,1,-1
            write (6,'(2f10.1,6e12.4)') &
                  zt (k), &
                  height (k), &
                  tnudge (k,t), &
                  unudge (k,t), &
                  vnudge (k,t), &
                  wnudge (k,t), &
                  thlnudge(k,t), &
                  qtnudge(k,t)
          end do
        end if
      end do readloop
      close(ifinput)
      if (znudgemin>0) then
        do k = 1,nzp-1
          if (zt(k)<=znudgemin) then
            tnudge(k,:) = 1e10
          else if (zt(k)<=znudgeplus) then
            tnudge(k,:)  = 2.*tnudgefac/(1-cos(pi*(zt(k)-znudgemin)/(znudgeplus-znudgemin)))
          else
            tnudge(k,:) = tnudgefac
          end if
        end do
      else
        tnudge  = tnudgefac*tnudge
      end if
      thlnudge = thlnudge - th00
      unudge = unudge - umean
      vnudge = vnudge - vmean
    lunudge = any(abs(unudge)>1e-8)
    lvnudge = any(abs(vnudge)>1e-8)
    lwnudge = any(abs(wnudge)>1e-8)
    lthlnudge = any(abs(thlnudge)>1e-8)
    lqtnudge  = any(abs(qtnudge)>1e-8)

  end subroutine initnudge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nudge(timein)
    use grid, only : dt, nxp, nyp, nzp, a_ut, a_vt, a_wt, a_tt, a_rt, a_up,a_vp,a_tp,a_rp,zt,a_wp
    use util, only : get_avg3
    implicit none
    real, intent (in) :: timein
    integer k,t,i,j
    real :: dtm,dtp,currtnudge, nudgefac
    real, dimension(nzp) :: uav, vav, tav, qav

    if (firsttime) then
      firsttime = .false.
      call initnudge(timein)
    end if
    if (.not.(lnudge)) return

    t=1
    do while(timein>timenudge(t))
      t=t+1
    end do
    if (timein/=timenudge(1)) then
      t=t-1
    end if

    dtm = ( timein-timenudge(t) ) / ( timenudge(t+1)-timenudge(t) )
    dtp = ( timenudge(t+1)-timein)/ ( timenudge(t+1)-timenudge(t) )

    call get_avg3(nzp, nxp, nyp,a_up,uav)
    call get_avg3(nzp, nxp, nyp,a_vp,vav)
    call get_avg3(nzp, nxp, nyp,a_tp,tav)
    call get_avg3(nzp, nxp, nyp,a_rp,qav)

    do j=3,nyp-2
    do i=3,nxp-2
    do k=2,nzp-1
      currtnudge = max(dt,tnudge(k,t)*dtp+tnudge(k,t+1)*dtm)
      if(lunudge  ) a_ut(k,i,j)=a_ut(k,i,j)-&
          (uav(k)-(unudge  (k,t)*dtp+unudge  (k,t+1)*dtm))/currtnudge
      if(lvnudge  ) a_vt(k,i,j)=a_vt(k,i,j)-&
          (vav(k)-(vnudge  (k,t)*dtp+vnudge  (k,t+1)*dtm))/currtnudge
    end do
    end do
    end do
    do j=3,nyp-2
    do i=3,nxp-2
    do k=2,nzp-1
      currtnudge = max(dt,tnudge(k,t)*dtp+tnudge(k,t+1)*dtm)
      if(lthlnudge  ) a_tt(k,i,j)=a_tt(k,i,j)-&
          (tav(k) - (thlnudge  (k,t)*dtp+thlnudge  (k,t+1)*dtm))/currtnudge
      if(lqtnudge  ) then
        if (qav(k) < qfloor .and. zt(k) < zfloor) then
          nudgefac = qfloor
!           currtnudge = 3600.
        else
          nudgefac = qtnudge  (k,t)*dtp+qtnudge  (k,t+1)*dtm
        end if
        a_rt(k,i,j)= a_rt(k,i,j) - (qav(k)-nudgefac)/currtnudge
      end if
    end do
    end do
    end do

  end subroutine nudge

! LINDA, b
!--------------------------------------------------------------------------!
!This routine nudges simulated values of temperature, humidity             !
!and horizontal and vertical wind speed back to their initial state at a   !
!time scale tau. The nudging is only done in a small zone at the Eastern   !
!and Western edge of the domain.                                           !
!The effect of the nudging is to obtain relaxation boundary conditions     !
!which can e.g. be used to model the development of a squall line          !
!                                                                          !
! Linda Schlemmer, December 2011                                           !
!--------------------------------------------------------------------------!
  subroutine nudge_bound

  use grid, only: nxp,nyp,nzp,a_tt,a_tp,a_rt,a_rp,&
                  a_ut,a_up,a_vt,a_vp,a_wt,a_wp,liquid,outtend,wtendr, nstep,rkalpha,rkbeta, dt
  use mpi_interface, only : myid

  IMPLICIT NONE
  integer:: k,t,i,j
  real    :: rk !rv


  if (firsttime) then
     firsttime = .false.
     !allocate(rlx(nxp,nyp)) 
     !allocate(tn(nzp),rn(nzp),un(nzp),vn(nzp),wn(nzp))
     allocate(rlz(nzp))
     allocate(rn(nzp))
     rlz(:)=0.0
     call initnudge_bound
     
     !test values:
     if(myid==0) print*,' ------------------------------------------------------------ '
     if(myid==0) print*,' ------------ in nudge_bound ------ rlz(:)=',rlz(kstrt-1:nzp)
     if(myid==0) print*,' kstrt=',kstrt
     if(myid==0) print*,' rn(kstrt:nzp)=',rn(kstrt:nzp)
     if(myid==0) print*,' coef=',coef
     if(myid==0) print*,' ------------------------------------------------------------ '
  end if

! relax towards initial conditions
  if(outtend) then !RV
     if(nstep==1) rk = rkalpha(1)+rkalpha(2)!then 
     if(nstep==2) rk = rkbeta(2)+rkbeta(3)
     if(nstep==3) rk = rkalpha(3)
  end if   !rv
  
  !do k=1,nzp
  do k=kstrt,nzp !rv
     do i=1,nxp
        do j=1,nyp
           !a_tt(k,i,j)=a_tt(k,i,j)-(a_tp(k,i,j)-tn(k))*coef*rlx(i,j)
           a_rt(k,i,j)=a_rt(k,i,j)-(a_rp(k,i,j)-rn(k))*coef*rlz(k)
           if(outtend) then !RV
              wtendr(k,i,j)=wtendr(k,i,j)-(a_rp(k,i,j)-rn(k))*coef*rlz(k)*rk*dt
           end if !rv
           !a_ut(k,i,j)=a_ut(k,i,j)-(a_up(k,i,j)-un(k))*coef*rlx(i,j)
           !a_vt(k,i,j)=a_vt(k,i,j)-(a_vp(k,i,j)-vn(k))*coef*rlx(i,j)
           !a_wt(k,i,j)=a_wt(k,i,j)-(a_wp(k,i,j)-wn(k))*coef*rlx(i,j)
        enddo
     enddo
  enddo
  
  end subroutine nudge_bound
!--------------------------------------------------------------------------


  subroutine initnudge_bound

  use mpi_interface, only : myid,nxprocs,nyprocs
  use grid, only: nxp,nyp,nzp,a_tp,a_rp,a_up,a_vp,a_wp,deltax,runtype,iradtyp,zt
  use defs, only: pi

  IMPLICIT NONE

  integer k,t,i,j
  real::eps=0.01
  integer:: nnudge ! number of points where relaxation is done
  real:: xnudge    ! relaxation zone [km], to be tested
  real:: rnudge    ! relaxation zone [points]
  logical::flg


  !RV, simple nudging for z>8.8km to get rid of liquid in large domain runs
  if(iradtyp == 3) then
     
     ! prepare rlz
     flg=.false.
     do k=1,nzp
        if(zt(k)<8800.) then
           rlz(k)=0.
        else
           rlz(k)=(zt(k)/zt(nzp))**10
           if(.not. flg) then
              flg=.true.
              kstrt=k
           end if
        end if
     enddo

     ! read in/save initial conditions
     if(runtype == 'INITIAL') then
        do k=1,nzp
           rn(k)=a_rp(k,3,3)
        enddo
        
     else !'HISTORY'; read from nudge_in file!
        open (21,file='nudge_in')
        !write(6,*) 'rn(k)' 
        do k=1,nzp
           read(21,*,end=100) rn(k)
           rn(k)=rn(k)/1000. !needs to be in kg/kg, nudge_in provides g/kg
           !write(6,*) rn(k)
        end do
        close(21)
     end if
100  continue!rv
     
  else     
     
     flg=.false.
     
     xnudge=10.0
     rnudge=xnudge*1000.0/deltax
     nnudge=int(rnudge)
     if (nnudge>nxp) flg=.true.
     
     ! western boundary
     if (mod(myid,nxprocs)<eps) then
        print*,'western boundary, myid=',myid
        if (flg) nnudge=nxp
        do i=1,nnudge
           !       rlx(i,:)=1.0 ! step function
           rlx(i,:)=(cos(real(i)*pi/(rnudge*2)))**2!cos2-function
        enddo
     endif
     ! relaxation zone extends over more than one processor
     if ((myid>0).and.(mod(myid-1,nxprocs)<eps).and.(flg)) then
        print*,'western boundary, 2nd proc, myid=',myid
        do i=1,nnudge-nxp+4
           !       rlx(i,:)=1.0 ! step function
           rlx(i,:)=(cos(real(i+nxp-4)*pi/(rnudge*2)))**2!cos2-function
        enddo
     endif
     
     ! eastern boundary
     if (mod(myid+1,nxprocs)<eps) then
        print*,'eastern boundary, myid=',myid
        if (flg) nnudge=nxp
        do i=nxp-nnudge+1,nxp
           !       rlx(i,:)=1.0! step function
           rlx(i,:)=(cos(real(nxp-i+1)*pi/(rnudge*2)))**2!cos2-function
        enddo
     endif
     ! relaxation zone extends over more than one processor
     if ((mod(myid+2,nxprocs)<eps).and.(flg)) then
        print*,'eastern boundary, 2nd proc, myid=',myid
        do i=2*nxp-nnudge+1-4,nxp
           !       rlx(i,:)=1.0 ! step function
           rlx(i,:)=(cos(real(2*nxp-i+1-4)*pi/(rnudge*2)))**2!cos2-function
        enddo
     endif
     
     ! read in/save initial conditions
     
     do k=1,nzp
        tn(k)=a_tp(k,3,3)
        rn(k)=a_rp(k,3,3)
        un(k)=a_up(k,3,3)
        vn(k)=a_vp(k,3,3)
        wn(k)=0.0
     enddo
     
  end if
  
  
end subroutine initnudge_bound
! LINDA, e

end module modnudge
