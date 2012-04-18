!> \file modtimedepsv.f90
!!  Prescribes surface values, fluxes and LS forcings at certain times for scalars

!>
!!  Prescribes surface values, fluxes and LS forcings at certain times for scalars
!>
!!  \author Roel Neggers, KNMI
!!  \author Thijs Heus,MPI-M
!!  \author Stephan de Roode, TU Delft
!!  \author Simon Axelsen, UU
!!  \par Revision list
!! \todo documentation
!  This file is part of DALES.
!
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

module modtimedep


implicit none

save
! switches for timedependent surface fluxes and large scale forcings
  logical       :: ltimedep     = .false. !< Overall switch, input in namoptions
  logical       :: ltimedepz    = .true.  !< Switch for large scale forcings
  logical       :: ltimedepsurf = .true.  !< Switch for surface fluxes
  logical       :: firsttime    = .true.
  integer, parameter    :: kflux = 200
  integer, parameter    :: kls   = 200
  real, allocatable     :: timeflux (:)
  real, allocatable     :: wqsurft  (:)
  real, allocatable     :: wtsurft  (:)
  real, allocatable     :: thlst    (:)
  real, allocatable     :: qtst     (:)
  real, allocatable     :: pst      (:)

  real, allocatable     :: timels  (:)
  real, allocatable     :: wflst   (:,:)
  real, allocatable     :: dqtdtlst(:,:)
  real, allocatable     :: dthldtlst(:,:)



contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedep(time,time_end)
!     use modglobal, only :ifnamopt,fname_options,dtmax, btime,cexpnr,nzp,nzp,ifinput,runtime
!     use modsurfdata,only :ps,qts,wqsurf,wtsurf,thls
!     use modtimedepsv, only : inittimedepsv
    use grid, only : nzp, zt
    real, intent(in) :: time,time_end
    character (80):: chmess
    character (1) :: chmess1
    integer :: k,t, ierr, ifinput = 19
    real :: dummyr, fac
!     real, allocatable, dimension (:) :: height
    real :: highheight,highwflst, highdthldt, highdqtdt
    real :: lowheight,  lowwflst,  lowdthldt,  lowdqtdt

    firsttime = .false.

    allocate(timeflux (0:kflux))
    allocate(wqsurft  (kflux))
    allocate(wtsurft  (kflux))
    allocate(thlst   (0:kls))
    allocate(qtst    (0:kls))
    allocate(pst    (0:kls))

    allocate(timels  (0:kls))
    allocate(wflst   (nzp,kls))
    allocate(dqtdtlst (nzp,kls))
    allocate(dthldtlst(nzp,kls))

    timels    = 0
    wflst     = 0
    dqtdtlst  = 0
    dthldtlst = 0


!    --- load lsforcings---


    open(ifinput,file='ls_flux_in')
    read(ifinput,'(a80)') chmess
    write(6,*) chmess
    read(ifinput,'(a80)') chmess
    write(6,*) chmess
    read(ifinput,'(a80)') chmess
    write(6,*) chmess

    timeflux = 0
    timels   = 0


!      --- load fluxes---

    t    = 0
    ierr = 0
    do while (timeflux(t) < (time_end))
      t=t+1
      read(ifinput,*, iostat = ierr) timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
      write(*,'(i8,6e12.4)') t,timeflux(t), wtsurft(t), wqsurft(t),thlst(t),qtst(t),pst(t)
      if (ierr < 0) then
          stop 'STOP: No time dependend data for end of run (surface fluxes)'
      end if
    end do
    if(timeflux(1)>(time_end)) then
        write(6,*) 'Time dependent surface variables do not change before end of'
        write(6,*) 'simulation. --> only large scale forcings'
        ltimedepsurf=.false.
    endif
! flush to the end of fluxlist
    do while (ierr ==0)
      read (ifinput,*,iostat=ierr) dummyr
    end do
    backspace (ifinput)
!     ---load large scale forcings----

    t = 0

    do while (timels(t) < time_end)
      t = t + 1
      chmess1 = "#"
      ierr = 1 ! not zero
      do while (.not.(chmess1 == "#" .and. ierr ==0)) !search for the next line consisting of "# time", from there onwards the profiles will be read
        read(ifinput,*,iostat=ierr) chmess1,timels(t)
        if (ierr < 0) then
          stop 'STOP: No time dependend data for end of run'
        end if
      end do
      write (*,*) 'timels = ',timels(t)
      read (ifinput,*)  lowheight , lowwflst,lowdthldt,lowdqtdt
      read (ifinput,*)  highheight , highwflst,highdthldt,highdqtdt
      do k=2,nzp
        if (highheight<zt(k)) then
          lowheight = highheight
          lowwflst  = highwflst
          lowdthldt = highdthldt
          lowdqtdt  = highdqtdt
          
          read (ifinput,*) highheight, highwflst, highdthldt, highdqtdt
        end if
        fac = (highheight-zt(k))/(highheight - lowheight)
        wflst(k,t)     = fac*lowwflst  + (1-fac)*highwflst
        dthldtlst(k,t) = fac*lowdthldt + (1-fac)*highdthldt
        dqtdtlst(k,t)  = fac*lowdqtdt  + (1-fac)*highdqtdt
      end do
    end do

    if ((timels(1) > time_end) .or. (timeflux(1) > time_end)) then
      write(6,*) 'Time dependent large scale forcings sets in after end of simulation -->'
      write(6,*) '--> only time dependent surface variables'
      ltimedepz=.false.
    end if

    close(ifinput)


  end subroutine inittimedep

  subroutine timedep(time, time_end, sst)

!-----------------------------------------------------------------|
!                                                                 |
!*** *timedep*  calculates ls forcings and surface forcings       |
!               case as a funtion of timee                        |
!                                                                 |
!      Roel Neggers    K.N.M.I.     01/05/2001                    |
!                                                                 |
!                                                                 |
!    calls                                                        |
!    * timedepz                                                   |
!      calculation of large scale advection, radiation and        |
!      surface fluxes by interpolation between prescribed         |
!      values at certain times                                    |
!                                                                 |
!    * timedepsurf                                                |
!      calculation  surface fluxes by interpolation               |
!      between prescribed values at certain times                 |
!                                                                 |
!                                                                 |
!-----------------------------------------------------------------|
    
    real, intent(in) :: time, time_end
    real, intent(inout) :: sst
    
    if (.not. ltimedep) return
    if (firsttime) call inittimedep(time,time_end)
    call timedepz(time)
    call timedepsurf(time, sst)
  end subroutine timedep

  subroutine timedepz(time)
    use grid, only : wfls,dthldtls,dqtdtls
!      use modfields,   only : ug, vg, dqtdtls,dqtdxls,dqtdyls, wfls,whls,thlprof,qtprof,thlpcar,dthldxls,dthldyls,dudxls,dudyls,dvdxls,dvdyls,dpdxl,dpdyl
!     use modglobal,   only : rtimee,om23_gs,zf,dzf,dzh,nzp,nzp,grav,llsadv
    real, intent(in) :: time

    integer t,k
    real fac

    if(.not.(ltimedepz)) return

    !---- interpolate ----
    t=1
    do while(time>timels(t))
      t=t+1
    end do
    if (time/=timels(1)) then
      t=t-1
    end if

    fac = ( time-timels(t) ) / ( timels(t+1)-timels(t) )
    wfls     = wflst    (:,t) + fac * ( wflst    (:,t+1) - wflst    (:,t) )
    dqtdtls  = dqtdtlst (:,t) + fac * ( dqtdtlst (:,t+1) - dqtdtlst (:,t) )
    dthldtls = dthldtlst(:,t) + fac * ( dthldtlst(:,t+1) - dthldtlst(:,t) )

  end subroutine timedepz

  subroutine timedepsurf(time, sst)
    use grid, only : psrf
    use srfc, only : dthcon, drtcon
    
!     use modglobal,   only : rtimee, lmoist
!     use modsurfdata, only : wtsurf,wqsurf,thls,qts,ps
!     use modsurface,  only : qtsurf

    implicit none
    real, intent(in) :: time
    real, intent(inout) :: sst

    integer t
    real fac


    if(.not.(ltimedepsurf)) return
  !     --- interpolate! ----
    t=1
    do while(time>timeflux(t))
      t=t+1
    end do
    if (time/=timeflux(t)) then
      t=t-1
    end if

    fac = ( time -timeflux(t) ) / ( timeflux(t+1)-timeflux(t))
    drtcon   = wqsurft(t) + fac * ( wqsurft(t+1) - wqsurft(t)  )
    dthcon   = wtsurft(t) + fac * ( wtsurft(t+1) - wtsurft(t)  )
    sst      = thlst(t)   + fac * ( thlst(t+1)   - thlst(t)    )
    psrf     = pst(t)     + fac * ( pst(t+1)   - pst(t)    )
  end subroutine timedepsurf



end module modtimedep
