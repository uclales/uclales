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
PUBLIC :: nudge,lnudge,tnudgefac
SAVE
  real, dimension(:,:), allocatable :: tnudge,unudge,vnudge,wnudge,thlnudge,qtnudge
  real, dimension(:)  , allocatable :: timenudge
  real :: tnudgefac = 1.
  logical :: lnudge,lunudge,lvnudge,lwnudge,lthlnudge,lqtnudge
  integer :: ntnudge = 100
  logical :: firsttime = .true.
contains
  subroutine initnudge(time)
    use grid, only : nzp,zt,th00
    use mpi_interface, only : myid
    implicit none

    integer :: ierr,k,t,ifinput = 19
    real,allocatable,dimension(:) :: height
    real, intent(in) :: time
    character(1) :: chmess1
    allocate(tnudge(nzp,ntnudge),unudge(nzp,ntnudge),vnudge(nzp,ntnudge),wnudge(nzp,ntnudge),thlnudge(nzp,ntnudge),qtnudge(nzp,ntnudge))
    allocate(timenudge(0:ntnudge), height(nzp))
    tnudge = 0
    unudge=0
    vnudge=0
    wnudge=0
    thlnudge=0
    qtnudge=0
    timenudge=0



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
        do  k=2,nzp-1
          read (ifinput,*) &
                height (k), &
                tnudge (k,t), &
                unudge (k,t), &
                vnudge (k,t), &
                wnudge (k,t), &
                thlnudge(k,t), &
                qtnudge(k,t)
        end do
        if (myid == 0) then
          do k=nzp-1,1,-1
            write (6,'(2f7.1,6e12.4)') &
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
      timenudge = timenudge/86400+time
      tnudge  = tnudgefac*tnudge
      thlnudge = thlnudge - th00
    lunudge = any(abs(unudge)>1e-8)
    lvnudge = any(abs(vnudge)>1e-8)
    lwnudge = any(abs(wnudge)>1e-8)
    lthlnudge = any(abs(thlnudge)>1e-8)
    lqtnudge  = any(abs(qtnudge)>1e-8)

  end subroutine initnudge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nudge(time)
!     use modglobal, only : rdt
!     use modfields, only : u0av,v0av,qt0av,thl0av
    use grid, only : dt, nxp, nyp, nzp, a_ut, a_vt, a_wt, a_tt, a_rt, a_up,a_vp,a_tp,a_rp,zt,a_wp
    implicit none
    real, intent (in) :: time
    integer k,t,i,j
    real :: dtm,dtp,currtnudge
    if (firsttime) then
      firsttime = .false.
      call initnudge(time)
    end if
    if (.not.(lnudge)) return
!     if (rk3step/=3) return

    t=1
    do while(time>timenudge(t))
      t=t+1
    end do
    if (time/=timenudge(1)) then
      t=t-1
    end if

    dtm = ( time-timenudge(t) ) / ( timenudge(t+1)-timenudge(t) )
    dtp = ( timenudge(t+1)-time)/ ( timenudge(t+1)-timenudge(t) )

    do j=3,nyp-2
    do i=3,nxp-2
    do k=2,nzp-1
      currtnudge = max(dt,tnudge(k,t)*dtp+tnudge(k,t+1)*dtm)
      if(lunudge  ) a_ut(k,i,j)=a_ut(k,i,j)-&
          (a_up(k,i,j)-(unudge  (k,t)*dtp+unudge  (k,t+1)*dtm))/currtnudge
      if(lvnudge  ) a_vt(k,i,j)=a_vt(k,i,j)-&
          (a_vp(k,i,j)-(vnudge  (k,t)*dtp+vnudge  (k,t+1)*dtm))/currtnudge
      if(lthlnudge  ) a_tt(k,i,j)=a_tt(k,i,j)-&
          (a_tp(k,i,j)-(thlnudge  (k,t)*dtp+thlnudge  (k,t+1)*dtm))/currtnudge
      if(lqtnudge  ) a_rt(k,i,j)=a_rt(k,i,j)-&
          (a_rp(k,i,j)-(qtnudge  (k,t)*dtp+qtnudge  (k,t+1)*dtm))/currtnudge
    end do
    end do
    end do
  end subroutine nudge


end module
