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

module modparticles
  !--------------------------------------------------------------------------
  ! module modparticles: Langrangian particle tracking, ad(o/a)pted from DALES
  !--------------------------------------------------------------------------
  implicit none
  PUBLIC :: init_particles, particles, exit_particles

  ! For/from namelist  
  logical            :: lpartic     = .false.
  ! ..................

  character(30)      :: startfile
  integer            :: ifinput     = 1
  integer            :: np

  ! Particle structure
  TYPE :: particle_record
    real     :: unique, tstart
    integer  :: partstep
    real     :: x, x_prev, ures_prev, xstart, ures, usgs, usgs_prev
    real     :: y, y_prev, vres_prev, ystart, vres, vsgs, vsgs_prev
    real     :: z, z_prev, wres_prev, zstart, wres, wsgs, wsgs_prev
    !real     :: sigma2_sgs

    TYPE (particle_record), POINTER :: next,prev
  end TYPE

  integer            :: nplisted
  TYPE (particle_record), POINTER :: head, tail

contains
  !--------------------------------------------------------------------------
  ! subroutine particles: Main routine, called every .. step
  !--------------------------------------------------------------------------
  subroutine particles
    use grid, only : deltax, deltay
    !use step, only : time
    implicit none
    type (particle_record), pointer:: particle

    if ( np < 1 ) return      ! Just to be sure..
    !if (lpartsgs) call sgsinit

    particle => head
    do while( associated(particle) )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! How to deal with circular references......?
      !if (  time - particle%tstart >= 0 ) then
  
        particle%partstep = particle%partstep + 1

        !interpolation of the velocity field
        !if (  rtimee - particle%tstart >= 0 ) then
          !particle%ures = velocity_ures(particle%x,particle%y,particle%z) / dx
          !particle%vres = velocity_vres(particle%x,particle%y,particle%z) / dy
          !particle%wres = velocity_wres(particle%x,particle%y,particle%z) / dzf(floor(particle%z))

          particle%ures = 1. / deltax

        !  if (lpartsgs) then
        !    if (rk3step==1) then
        !      particle%usgs_prev = particle%usgs
        !      particle%vsgs_prev = particle%vsgs
        !      particle%wsgs_prev = particle%wsgs
        !    end if
        !    call sgshelpvar(particle)
        !    particle%usgs = velocity_usgs(particle) / dx
        !    particle%vsgs = velocity_vsgs(particle) / dy
        !    particle%wsgs = velocity_wsgs(particle) / dzf(floor(particle%z))
        !  end if
        !end if
      !end if


    particle => particle%next
    end do
    
    !****Statistics******
    !if (rk3step==3) then
    !  call statistics
    !  call writeparticles
    !end if

    !Time integration
    particle => head
    do while( associated(particle) )
      !if ( time - particle%tstart >= 0 ) then
    !    select case(intmeth)
    !    case(inomove)
    !            ! no movement
    !    case(irk3)
           call rk3(particle)
    !    case default
    !            stop 'PARTICLES ERROR: incorrect integration scheme'
    !    end select
      !end if
    particle => particle%next
    end do

    !Exchange particle to other processors
    !call partcommunicate

  end subroutine particles

  !--------------------------------------------------------------------------
  ! subroutine rk3
  !--------------------------------------------------------------------------
  subroutine rk3(particle)
    use grid, only : rkalpha, rkbeta, nstep, dt, deltax
    !use step, only : time
    implicit none
    real :: rk3coef

    TYPE (particle_record), POINTER:: particle

    rk3coef = dt / (4. - dble(nstep))
    particle%x   = particle%x + rkalpha(nstep) * (particle%ures+particle%usgs) * dt + rkbeta(nstep) * (particle%ures_prev + particle%usgs_prev) * dt
    particle%y   = particle%y + rkalpha(nstep) * (particle%vres+particle%vsgs) * dt + rkbeta(nstep) * (particle%vres_prev + particle%vsgs_prev) * dt
    particle%z   = particle%z + rkalpha(nstep) * (particle%wres+particle%wsgs) * dt + rkbeta(nstep) * (particle%wres_prev + particle%wsgs_prev) * dt

    particle%ures_prev = particle%ures
    particle%vres_prev = particle%vres
    particle%wres_prev = particle%wres

    !  call checkbound(particle)
    !  if (floor(particle%z)/=floor(particle%z_prev)) then
    !    particle%z = floor(particle%z) + (particle%z -floor(particle%z))*dzf(floor(particle%z_prev))/dzf(floor(particle%z))
    !  end if

   if ( nstep==3 ) then
      particle%ures_prev = 0.
      particle%vres_prev = 0.
      particle%wres_prev = 0.
    end if

  end subroutine rk3

  !--------------------------------------------------------------------------
  ! subroutine init_particles: initialize particles, reading initial position, 
  ! etc. Called from subroutine initialize (init.f90)
  !--------------------------------------------------------------------------
  subroutine init_particles
    use mpi_interface, only : myid
    use grid, only : zm, dzi_m, deltax, deltay, zt,dzi_t

    integer :: k, n, ierr, kmax
    real :: tstart, xstart, ystart, zstart
    type (particle_record), pointer:: particle

    kmax = size(zm)
    
    np = 0
    startfile = 'partstartpos'
    open(ifinput,file=startfile,status='old',position='rewind',action='read')
    read(ifinput,'(I10.1)') np
    if ( np < 1 ) return
   
    ! clear pointers to head and tail, sets nplisted=0
    call particle_initlist()

    ! read particles from partstartpos, create linked list
    do n = 1, np
      read(ifinput,*) tstart, xstart, ystart, zstart

    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! NEEDS TO BE CHANGED FOR 2D PARALLEL
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !if (floor(ystart/ysizelocal) == nprocs) ystart = 0.
      !if (floor(ystart/ysizelocal) == myid) then
        call add_particle(particle)

        particle%unique = n + myid/1000.0
        particle%x = xstart/deltax !+2
        particle%y = ystart/deltay !+2
        do k=kmax,1,-1
          if ( zm(k)<zstart ) exit
        end do
        particle%z = k + (zstart-zm(k))*dzi_m(k)
        particle%xstart = xstart
        particle%ystart = ystart
        particle%zstart = zstart
        particle%tstart = tstart
        particle%ures = 0.
        particle%vres = 0.
        particle%wres = 0.
        particle%usgs = 0.
        particle%vsgs = 0.
        particle%wsgs = 0.
        particle%usgs_prev = 0.
        particle%vsgs_prev = 0.
        particle%wsgs_prev = 0.
        particle%x_prev = particle%x
        particle%y_prev = particle%y
        particle%z_prev = particle%z
        particle%partstep = 0
    !    particle%sigma2_sgs = epsilon(particle%sigma2_sgs)
    !  !end if
    end do

    close(ifinput)
    !write(6,*) 'timee :',rtimee,': proc ',myid,': #particles: ',nplisted

  end subroutine init_particles

  !--------------------------------------------------------------------------
  ! subroutine exit_particles
  !--------------------------------------------------------------------------
  subroutine exit_particles
    implicit none

    do while( associated(tail) )
      call delete_particle(tail)
    end do

  print "(//' ',49('-')/,' ',/,'  Lagrangian particles removed.')"
  end subroutine exit_particles

  !--------------------------------------------------------------------------
  ! subroutine init_particles: descr
  !--------------------------------------------------------------------------
  subroutine add_particle(ptr)
    implicit none
    TYPE (particle_record), POINTER:: ptr
    TYPE (particle_record), POINTER:: new_p

    if( .not. associated(head) ) then 
      allocate(head)
      tail => head
      nullify(head%prev)
    else
      allocate(tail%next)
      new_p => tail%next
      new_p%prev => tail
      tail => new_p
    end if

    nplisted = nplisted + 1
    nullify(tail%next)
    ptr => tail

  end SUBROUTINE add_particle

  !--------------------------------------------------------------------------
  ! subroutine init_particles: clears pointers to start and end of linked
  ! particle list
  !--------------------------------------------------------------------------
  subroutine particle_initlist()
    implicit none

    nullify(head)
    nullify(tail)
    nplisted = 0
  end subroutine particle_initlist

  !--------------------------------------------------------------------------
  ! subroutine delete_particle
  !--------------------------------------------------------------------------
  subroutine delete_particle(ptr)
    implicit none
    TYPE (particle_record), POINTER:: ptr
    TYPE (particle_record), POINTER:: next_p,prev_p,cur_p

    cur_p => ptr
    if( .not. associated(cur_p) ) then         !error in calling ptr
      write(6,*) 'WARNING: cannot delete empty pointer'
      return
    end if

    if( .not. associated(head) ) then         !empty list
      write(6,*) 'WARNING: cannot delete elements in an empty list'
    else
      if( .not. associated(cur_p%next) ) then   ! last in list
        if( .not. associated(cur_p%prev) ) then ! last element
          nullify(head)
          nullify(tail)
        else
          tail => cur_p%prev
          nullify(tail%next)
        end if
      else
        if( .not. associated(cur_p%prev) ) then ! first in list
          if( .not. associated(cur_p%next) ) then !last element
            nullify(head)
            nullify(tail)
          else
            head => cur_p%next
            nullify(head%prev)
          end if
        else
          next_p => cur_p%next
          prev_p => cur_p%prev
          next_p%prev => prev_p
          prev_p%next => next_p
        end if
      end if

      nplisted = nplisted - 1
      deallocate(cur_p)
    end if

  end subroutine delete_particle

end module modparticles
