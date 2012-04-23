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
  logical            :: lpartdump   = .true.
  real               :: frqpartdump = 5

  character(30)      :: startfile
  integer            :: ifinput     = 1
  integer            :: np
  integer            :: tnextdump

  ! Particle structure
  type :: particle_record
    real     :: unique, tstart
    integer  :: partstep
    real     :: x, x_prev, ures_prev, xstart, ures, usgs, usgs_prev
    real     :: y, y_prev, vres_prev, ystart, vres, vsgs, vsgs_prev
    real     :: z, z_prev, wres_prev, zstart, wres, wsgs, wsgs_prev
    !real     :: sigma2_sgs

    type (particle_record), pointer :: next,prev
  end type

  integer            :: nplisted
  type (particle_record), pointer :: head, tail

  integer            :: ipunique, ipx, ipy, ipz, ipxstart, ipystart, ipzstart, iptsart, ipxprev, ipyprev, ipzprev
  integer            :: ipures, ipvres, ipwres, ipusgs, ipvsgs, ipwsgs, ipusgs_prev, ipvsgs_prev, ipwsgs_prev 
  integer            :: ipures_prev, ipvres_prev, ipwres_prev, ipartstep, nrpartvar

contains
  !
  !--------------------------------------------------------------------------
  ! subroutine particles: Main routine, called every .. step
  !--------------------------------------------------------------------------
  !
  subroutine particles(time)
    use grid, only : deltax, deltay, nstep, dzi_t, dt
    implicit none
    real, intent(in)               :: time
    type (particle_record), pointer:: particle

    if ( np < 1 ) return      ! Just to be sure..
    !if (lpartsgs) call sgsinit

    particle => head
    do while( associated(particle) )
      if (  time - particle%tstart >= 0 ) then
        particle%partstep = particle%partstep + 1

        ! Interpolation of the velocity field
        particle%ures = velocity_ures(particle%x,particle%y,particle%z) / deltax
        particle%vres = velocity_vres(particle%x,particle%y,particle%z) / deltay
        particle%wres = velocity_wres(particle%x,particle%y,particle%z) * dzi_t(floor(particle%z))

        !if (lpartsgs) then
        !  if (rk3step==1) then
        !    particle%usgs_prev = particle%usgs
        !    particle%vsgs_prev = particle%vsgs
        !    particle%wsgs_prev = particle%wsgs
        !  end if
        !  call sgshelpvar(particle)
        !  particle%usgs = velocity_usgs(particle) / dx
        !  particle%vsgs = velocity_vsgs(particle) / dy
        !  particle%wsgs = velocity_wsgs(particle) / dzf(floor(particle%z))
        !end if
      end if

    particle => particle%next
    end do
    
    !Time integration
    particle => head
    do while( associated(particle) )
      if ( time - particle%tstart >= 0 ) then
        call rk3(particle)
      end if

    if (nstep==3) then
      if(time + dt >= tnextdump) then
        call writeparticles   ! OLD!
        call particledump


        tnextdump = tnextdump + frqpartdump
      end if

      !call statistics
    end if

    particle => particle%next
    end do

    !Exchange particle to other processors
    call partcomm

  end subroutine particles


  !
  !--------------------------------------------------------------------------
  ! subroutine initparticledump
  !--------------------------------------------------------------------------
  !
  subroutine initparticledump
    use modnetcdf,       only : open_nc, addvar_nc
    use grid,            only : nzp, nxp, tname, tlongname, tunit, filprf
    implicit none
    
  end subroutine initparticledump

  !
  !--------------------------------------------------------------------------
  ! subroutine particledump
  !--------------------------------------------------------------------------
  !
  subroutine particledump
    use mpi_interface, only : mpi_comm_world, myid, mpi_integer, mpi_double_precision, ierror, wrxid, wryid, nxprocs, nyprocs,ranktable, mpi_status_size
    implicit none
    type (particle_record), pointer:: particle
    integer                :: nlocal, ii, i, start
    integer, allocatable   :: nremote(:)
    integer                :: status(mpi_status_size)
    integer                :: nvar = 4       ! id,x,y,z,u,v,w,..,..,..,
    real, allocatable      :: sendbuff(:), recvbuff(:)

    ! Count local particles
    nlocal = 0
    particle => head
    do while( associated(particle) )
      nlocal = nlocal + 1
      particle => particle%next
    end do

    ! Communicate number of local particles to main proces (0)
    allocate(nremote(0:(nxprocs*nyprocs)-1))
    call mpi_gather(nlocal,1,mpi_integer,nremote,1,mpi_integer,0,mpi_comm_world,ierror)

    ! Create buffer
    allocate(sendbuff(nvar * nlocal))
    ii = 1
    particle => head
    do while( associated(particle) )
      sendbuff(ii)   = particle%unique
      sendbuff(ii+1) = particle%x
      sendbuff(ii+2) = particle%y
      sendbuff(ii+3) = particle%z

      !print*,sendbuff(ii),sendbuff(ii+1),sendbuff(ii+2),sendbuff(ii+3)

      ii = ii + nvar
      particle => particle%next
    end do

    print*,'about to send...'

    call mpi_send(sendbuff,1,mpi_double_precision,0,myid,mpi_comm_world,ierror)

    print*,'did send..'    

    if(myid == 0) then
      do i = 0,(nxprocs*nyprocs)-1
        allocate(recvbuff(nremote(i)*nvar))
        call mpi_recv(recvbuff,1,mpi_double_precision,i,i,mpi_comm_world,status,ierror)   
        print*,i,recvbuff   
        deallocate(recvbuff)
      end do
    end if


    deallocate(nremote)
 
  end subroutine particledump


  !
  !--------------------------------------------------------------------------
  ! subroutine partcomm
  !--------------------------------------------------------------------------
  !
  subroutine partcomm
    use mpi_interface, only : wrxid, wryid, ranktable, nxg, nyg, xcomm, ycomm, ierror, mpi_status_size, mpi_integer, mpi_double_precision, mpi_comm_world, nyprocs, nxprocs
    implicit none

    type (particle_record), pointer:: particle,ptr
    real, allocatable, dimension(:) :: buffsend, buffrecv
    integer :: status(mpi_status_size)
    integer :: ii, n
    ! Number of particles to ('to') and from ('fr') N,E,S,W
    integer :: nton,ntos,ntoe,ntow
    integer :: nfrn,nfrs,nfre,nfrw
    integer :: nyloc, nxloc 

    nton = 0
    ntos = 0
    ntoe = 0
    ntow = 0

    nyloc = nyg / nyprocs
    nxloc = nxg / nxprocs

    ! --------------------------------------------
    ! First: all north to south (j) and vice versa
    ! --------------------------------------------
    particle => head
    do while(associated(particle) )
      if( particle%y >= nyloc + 3 ) nton = nton + 1
      if( particle%y < 3          ) ntos = ntos + 1
      particle => particle%next
    end do 

    call mpi_sendrecv(nton,1,mpi_integer,ranktable(wrxid,wryid+1),4, &
                      nfrs,1,mpi_integer,ranktable(wrxid,wryid-1),4, &
                      mpi_comm_world,status,ierror) 

    call mpi_sendrecv(ntos,1,mpi_integer,ranktable(wrxid,wryid-1),5, &
                      nfrn,1,mpi_integer,ranktable(wrxid,wryid+1),5, &
                      mpi_comm_world,status,ierror) 

    !if( nton > 0 ) allocate(buffsend(nrpartvar * nton))
    !if( ntos > 0 ) allocate(buffsend(nrpartvar * ntos))
    allocate(buffsend(nrpartvar * nton))
    allocate(buffrecv(nrpartvar * nfrs))

    if( nton > 0 ) then
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%y >= nyloc + 3 ) then
          particle%y      = particle%y      - nyloc
          particle%y_prev = particle%y_prev - nyloc

          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do
    end if

    call mpi_sendrecv(buffsend,nrpartvar*nton,mpi_double_precision,ranktable(wrxid,wryid+1),6, &
                      buffrecv,nrpartvar*nfrs,mpi_double_precision,ranktable(wrxid,wryid-1),6, &
                      mpi_comm_world, status, ierror)

    ii = 0
    do n = 1,nfrs
      call add_particle(particle)
      call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
      ii=ii+nrpartvar
    end do

    !if( nton > 0 ) deallocate(buffsend)
    !if( nfrs > 0 ) deallocate(buffrecv)
    !if( ntos > 0 ) allocate(buffsend(nrpartvar*ntos))
    !if( nfrn > 0 ) allocate(buffrecv(nrpartvar*nfrn))
    deallocate(buffsend)
    deallocate(buffrecv)
    allocate(buffsend(nrpartvar*ntos))
    allocate(buffrecv(nrpartvar*nfrn))

    if( ntos > 0 ) then
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%y < 3 ) then
          particle%y      = particle%y      + nyloc
          particle%y_prev = particle%y_prev + nyloc

          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)

          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do

    end if

    call mpi_sendrecv(buffsend,nrpartvar*ntos,mpi_double_precision,ranktable(wrxid,wryid-1),7, &
                      buffrecv,nrpartvar*nfrn,mpi_double_precision,ranktable(wrxid,wryid+1),7, &
                      mpi_comm_world, status, ierror)

    ii = 0
    do n = 1,nfrn
      particle => head
      call add_particle(particle)
      call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
      ii=ii+nrpartvar
    end do

    !if( ntos > 0 ) deallocate(buffsend)
    !if( nfrn > 0 ) deallocate(buffrecv)
    deallocate(buffsend)
    deallocate(buffrecv)

    ! --------------------------------------------
    ! Second: all east to west (i) and vice versa
    ! --------------------------------------------
    particle => head
    do while(associated(particle) )
      if( particle%x >= nxloc + 3 ) ntoe = ntoe + 1
      if( particle%x < 3          ) ntow = ntow + 1
      particle => particle%next
    end do 

    call mpi_sendrecv(ntoe,1,mpi_integer,ranktable(wrxid+1,wryid),8, &
                      nfrw,1,mpi_integer,ranktable(wrxid-1,wryid),8, &
                      mpi_comm_world,status,ierror) 

    call mpi_sendrecv(ntow,1,mpi_integer,ranktable(wrxid-1,wryid),9, &
                      nfre,1,mpi_integer,ranktable(wrxid+1,wryid),9, &
                      mpi_comm_world,status,ierror) 

    !if(ntoe > 0) allocate(buffsend(nrpartvar * ntoe))
    !if(ntow > 0) allocate(buffsend(nrpartvar * ntow))
    allocate(buffsend(nrpartvar * ntoe))
    allocate(buffrecv(nrpartvar * nfrw))

    if( ntoe > 0 ) then
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%x >= nxloc + 3 ) then
          particle%x      = particle%x      - nxloc
          particle%x_prev = particle%x_prev - nxloc

          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do
    end if

    call mpi_sendrecv(buffsend,nrpartvar*ntoe,mpi_double_precision,ranktable(wrxid+1,wryid),10, &
                      buffrecv,nrpartvar*nfrw,mpi_double_precision,ranktable(wrxid-1,wryid),10, &
                      mpi_comm_world, status, ierror)

    ii = 0
    do n = 1,nfrw
      call add_particle(particle)
      call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
      ii=ii+nrpartvar
    end do

    !if( ntoe > 0 ) deallocate(buffsend)
    !if( nfrw > 0 ) deallocate(buffrecv)
    !if( ntow > 0 ) allocate(buffsend(nrpartvar*ntow))
    !if( nfre > 0 ) allocate(buffrecv(nrpartvar*nfre))
    deallocate(buffsend)
    deallocate(buffrecv)
    allocate(buffsend(nrpartvar*ntow))
    allocate(buffrecv(nrpartvar*nfre))

    if( ntow > 0 ) then
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%x < 3 ) then
          particle%x      = particle%x      + nxloc
          particle%x_prev = particle%x_prev + nxloc

          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)

          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do

    end if

    call mpi_sendrecv(buffsend,nrpartvar*ntow,mpi_double_precision,ranktable(wrxid-1,wryid),11, &
                      buffrecv,nrpartvar*nfre,mpi_double_precision,ranktable(wrxid+1,wryid),11, &
                      mpi_comm_world, status, ierror)

    ii = 0
    do n = 1,nfre
      particle => head
      call add_particle(particle)
      call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
      ii=ii+nrpartvar
    end do

    !if( ntow > 0 ) deallocate(buffsend)
    !if( nfre > 0 ) deallocate(buffrecv)
    deallocate(buffsend)
    deallocate(buffrecv)

  end subroutine partcomm

  !
  !--------------------------------------------------------------------------
  ! subroutine partbuffer
  !--------------------------------------------------------------------------
  !
  subroutine partbuffer(particle, buffer, n, send)
    implicit none

    logical,intent(in)                :: send
    integer,intent(in)                :: n
    real,dimension(n+1:n+nrpartvar)   :: buffer
    TYPE (particle_record), POINTER:: particle

    if (send) then
      buffer(n+ipunique)    = particle%unique
      buffer(n+ipx)         = particle%x
      buffer(n+ipy)         = particle%y
      buffer(n+ipz)         = particle%z
      buffer(n+ipures)      = particle%ures
      buffer(n+ipvres)      = particle%vres
      buffer(n+ipwres)      = particle%wres
      buffer(n+ipures_prev) = particle%ures_prev
      buffer(n+ipvres_prev) = particle%vres_prev
      buffer(n+ipwres_prev) = particle%wres_prev
      buffer(n+ipusgs)      = particle%usgs
      buffer(n+ipvsgs)      = particle%vsgs
      buffer(n+ipwsgs)      = particle%wsgs
      buffer(n+ipusgs_prev) = particle%usgs_prev
      buffer(n+ipvsgs_prev) = particle%vsgs_prev
      buffer(n+ipwsgs_prev) = particle%wsgs_prev
      buffer(n+ipxstart)    = particle%xstart
      buffer(n+ipystart)    = particle%ystart
      buffer(n+ipzstart)    = particle%zstart
      buffer(n+iptsart)     = particle%tstart
      buffer(n+ipartstep)   = particle%partstep
      buffer(n+ipxprev)     = particle%x_prev
      buffer(n+ipyprev)     = particle%y_prev
      buffer(n+ipzprev)     = particle%z_prev
      
      !if (intmeth == irk3 .or. lpartsgs) then
      !   buffer(n+ipxprev) = particle%x_prev
      !   buffer(n+ipyprev) = particle%y_prev
      !   buffer(n+ipzprev) = particle%z_prev
      ! end if
      ! if (lpartsgs) then
      !   buffer(n+ipsigma2_sgs)=particle%sigma2_sgs
      ! end if
    else
      particle%unique       = buffer(n+ipunique)
      particle%x            = buffer(n+ipx)
      particle%y            = buffer(n+ipy)
      particle%z            = buffer(n+ipz)
      particle%ures         = buffer(n+ipures)
      particle%vres         = buffer(n+ipvres)
      particle%wres         = buffer(n+ipwres)
      particle%ures_prev    = buffer(n+ipures_prev)
      particle%vres_prev    = buffer(n+ipvres_prev)
      particle%wres_prev    = buffer(n+ipwres_prev)
      particle%usgs         = buffer(n+ipusgs)
      particle%vsgs         = buffer(n+ipvsgs)
      particle%wsgs         = buffer(n+ipwsgs)
      particle%usgs_prev    = buffer(n+ipusgs_prev)
      particle%vsgs_prev    = buffer(n+ipvsgs_prev)
      particle%wsgs_prev    = buffer(n+ipwsgs_prev)
      particle%xstart       = buffer(n+ipxstart)
      particle%ystart       = buffer(n+ipystart)
      particle%zstart       = buffer(n+ipzstart)
      particle%tstart       = buffer(n+iptsart)
      particle%partstep     = buffer(n+ipartstep)
      particle%x_prev       = buffer(n+ipxprev)
      particle%y_prev       = buffer(n+ipyprev)
      particle%z_prev       = buffer(n+ipzprev)
      
      !if (intmeth == irk3 .or. lpartsgs) then
      !  particle%x_prev =  buffer(n+ipxprev)
      !  particle%y_prev =  buffer(n+ipyprev)
      !  particle%z_prev =  buffer(n+ipzprev)
      !end if
      !if (lpartsgs) then
      !  particle%sigma2_sgs=buffer(n+ipsigma2_sgs)
      !end if
    end if

  end subroutine partbuffer

  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_ures: trilinear interpolation of u-component
  !--------------------------------------------------------------------------
  !
  function velocity_ures(x,y,z)
    use grid, only : a_up, dzi_m, dzi_t, zt, zm
    implicit none
    real, intent(in) :: x, y, z
    integer          :: xbottom, ybottom, zbottom
    real             :: velocity_ures, deltax, deltay, deltaz, sign

    xbottom = floor(x) - 1
    ybottom = floor(y - 0.5)
    zbottom = floor(z + 0.5)
    deltax = x - 1   - xbottom
    deltay = y - 0.5 - ybottom

    ! u(1,:,:) == u(2,:,:) with zt(1) = - zt(2). By multiplying u(1,:,:) with -1, 
    ! the velocity interpolates to 0 at the surface.  
    if (zbottom==1)  then
      sign = -1
    else
      sign = 1
    end if      

    deltaz = ((zm(floor(z)) + (z - floor(z)) / dzi_t(floor(z))) - zt(floor(z))) * dzi_m(floor(z))
    velocity_ures =  (1-deltaz) * (1-deltay) * (1-deltax) * sign * a_up(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) * sign * a_up(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) * sign * a_up(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) * sign * a_up(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *        a_up(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *        a_up(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *        a_up(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *        a_up(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function velocity_ures

  !--------------------------------------------------------------------------
  ! subroutine velocity_vres: trilinear interpolation of v-component
  !--------------------------------------------------------------------------
  function velocity_vres(x,y,z)
    use grid, only : a_vp, dzi_m, dzi_t, zt, zm
    implicit none
    real, intent(in) :: x, y, z
    integer          :: xbottom, ybottom, zbottom
    real             :: velocity_vres, deltax, deltay, deltaz, sign

    xbottom = floor(x - 0.5)
    ybottom = floor(y) - 1
    zbottom = floor(z + 0.5)
    deltax = x - 0.5 - xbottom
    deltay = y - 1   - ybottom

    ! v(1,:,:) == v(2,:,:) with zt(1) = - zt(2). By multiplying v(1,:,:) with -1, 
    ! the velocity interpolates to 0 at the surface.  
    if (zbottom==1)  then
      sign = -1
    else
      sign = 1
    end if      

    deltaz = ((zm(floor(z)) + (z - floor(z)) / dzi_t(floor(z))) - zt(floor(z))) * dzi_m(floor(z))
    velocity_vres =  (1-deltaz) * (1-deltay) * (1-deltax) * sign * a_vp(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) * sign * a_vp(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) * sign * a_vp(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) * sign * a_vp(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *        a_vp(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *        a_vp(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *        a_vp(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *        a_vp(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function velocity_vres

  !--------------------------------------------------------------------------
  ! subroutine velocity_wres: trilinear interpolation of w-component
  !--------------------------------------------------------------------------
  function velocity_wres(x,y,z)
    use grid, only : a_wp, dzi_m, dzi_t, zt, zm
    implicit none
    real, intent(in) :: x, y, z
    integer          :: xbottom, ybottom, zbottom
    real             :: velocity_wres, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z)
    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z - zbottom

    velocity_wres =  (1-deltaz) * (1-deltay) * (1-deltax) *  a_wp(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) *  a_wp(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) *  a_wp(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) *  a_wp(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *  a_wp(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *  a_wp(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *  a_wp(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *  a_wp(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1
   
  end function velocity_wres

  !--------------------------------------------------------------------------
  ! subroutine rk3
  !--------------------------------------------------------------------------
  subroutine rk3(particle)
    use grid, only : rkalpha, rkbeta, nstep, dt, deltax
    implicit none
    TYPE (particle_record), POINTER:: particle

    particle%x   = particle%x + rkalpha(nstep) * (particle%ures+particle%usgs) * dt + rkbeta(nstep) * (particle%ures_prev + particle%usgs_prev) * dt
    particle%y   = particle%y + rkalpha(nstep) * (particle%vres+particle%vsgs) * dt + rkbeta(nstep) * (particle%vres_prev + particle%vsgs_prev) * dt
    particle%z   = particle%z + rkalpha(nstep) * (particle%wres+particle%wsgs) * dt + rkbeta(nstep) * (particle%wres_prev + particle%wsgs_prev) * dt

    particle%ures_prev = particle%ures
    particle%vres_prev = particle%vres
    particle%wres_prev = particle%wres

    call checkbound(particle)
    
    !  if (floor(particle%z)/=floor(particle%z_prev)) then
    !    particle%z = floor(particle%z) + (particle%z -floor(particle%z))*dzf(floor(particle%z_prev))/dzf(floor(particle%z))
    !  end if

   if ( nstep==3 ) then
      particle%ures_prev   = 0.
      particle%vres_prev   = 0.
      particle%wres_prev   = 0.
      particle%x_prev      = particle%x
      particle%y_prev      = particle%y
      particle%z_prev      = particle%z
    end if

  end subroutine rk3

  !--------------------------------------------------------------------------
  ! subroutine checkbound
  !--------------------------------------------------------------------------
  subroutine checkbound(particle)
    use grid, only : nxp, nyp, zm
    implicit none

    type (particle_record), pointer:: particle

    ! Cyclic boundaries -> not needed with 2d parallel grid
    !particle%x      = modulo(particle%x-3,real(nxp-4))+3
    !particle%y      = modulo(particle%y-3,real(nyp-4))+3
    !particle%x_prev = modulo(particle%x_prev-3,real(nxp-4))+3
    !particle%y_prev = modulo(particle%y_prev-3,real(nyp-4))+3

    ! Reflect particles of surface and model top
    if (particle%z >= size(zm)) then
      particle%z = size(zm)-0.0001
      particle%wres = -abs(particle%wres)
    elseif (particle%z < 1.01) then
      particle%z = abs(particle%z-1.01)+1.01
      particle%wres =  abs(particle%wres)
    end if

    if (particle%z_prev >= size(zm)) then
      particle%z_prev = size(zm)-0.0001
    elseif (particle%z_prev < 1.01) then
      particle%z_prev = abs(particle%z_prev-1.01)+1.01
    end if

  end subroutine checkbound

  !--------------------------------------------------------------------------
  ! subroutine writeparticles -> raw dump of particle field
  !--------------------------------------------------------------------------
  subroutine writeparticles
    use mpi_interface, only : wrxid, wryid, nxg, nyg, myid, nxprocs, nyprocs 
    use grid, only : deltax, deltay, dzi_t, zm
    implicit none
    
    real, allocatable,dimension (:,:) :: partdata
    integer(KIND=selected_int_kind(11)), allocatable, dimension (:) :: partids
    type (particle_record), pointer:: particle
    integer :: ndata = 3         ! Hardcoded for now
    integer :: n, m
    integer :: ifoutput = 99
    character(3) :: cmyid
    write(cmyid,'(i3.3)') myid  
 
    n = 0

    allocate (partdata(ndata,nplisted))
    allocate (partids(nplisted))

    particle => head
    do while( associated(particle) )
      n = n + 1
      partids(n) = particle%unique

      if (ndata > 0) then
         partdata(1,n) = (wrxid * (nxg / nxprocs) + particle%x - 3) * deltax
      endif
      if (ndata > 1) then
         partdata(2,n) = (wryid * (nyg / nyprocs) + particle%y - 3) * deltay
      endif
      if (ndata > 2) then
         partdata(3,n) = zm(floor(particle%z)) + (particle%z-floor(particle%z)) / dzi_t(floor(particle%z))
      endif
      !if (ndata > 3) then
      !   partdata(4,n) = (particle%ures+particle%usgs)*dx
      !endif
      !if (ndata > 4) then
      !   partdata(5,n) = (particle%vres+particle%vsgs)*dy
      !endif
      !if (ndata > 5) then
      !   partdata(6,n) = (particle%wres+particle%wsgs)*dzf(floor(particle%z))
      !endif
      !if (ndata > 6) then
      !   partdata(7,n) = thlpart
      !endif
      !if (ndata > 7) then
      !   partdata(8,n) = thvpart
      !endif
      !if (ndata > 8) then
      !   partdata(9,n) = qtpart * 1000.
      !endif
      !if (ndata > 9) then
      !   partdata(10,n)= qlpart * 1000.
      !endif
      particle => particle%next
    end do



    ! Shouldn't be here......
    open(ifoutput,file='particles.'//cmyid,position='append',action='write')
    write(ifoutput,'(A2,I10,I10)') '# ',tnextdump,nplisted
    write(ifoutput,'(I10,3F12.5)') (partids(n),(partdata(m,n),m=1,ndata),n=1,nplisted)
    close(ifoutput) 
    
    deallocate (partdata)
    deallocate (partids)

  end subroutine writeparticles

  !--------------------------------------------------------------------------
  !
  ! BELOW: ONLY INIT / EXIT PARTICLES
  !
  !--------------------------------------------------------------------------
  ! subroutine init_particles: initialize particles, reading initial position, 
  ! etc. Called from subroutine initialize (init.f90)
  !--------------------------------------------------------------------------
  subroutine init_particles
    use mpi_interface, only : wrxid, wryid, nxg, nyg, myid, nxprocs, nyprocs
    use grid, only : zm, deltax, deltay, zt,dzi_t
    use grid, only : a_up, a_vp, a_wp

    integer :: k, n, ierr, kmax
    real :: tstart, xstart, ystart, zstart, ysizelocal, xsizelocal
    type (particle_record), pointer:: particle

    xsizelocal = (nxg / nxprocs) * deltax
    ysizelocal = (nyg / nyprocs) * deltay
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
      if(floor(xstart / xsizelocal) == wrxid) then
        if(floor(ystart / ysizelocal) == wryid) then
          call add_particle(particle)
          particle%unique      = n !+ myid/1000.0
          particle%x           = (xstart - (float(wrxid) * xsizelocal)) / deltax + 3.  ! +3 here for ghost cells.
          particle%y           = (ystart - (float(wryid) * ysizelocal)) / deltay + 3.  ! +3 here for ghost cells.
          do k=kmax,1,1
            if ( zm(k)<zstart ) exit
          end do
          particle%z           = k + (zstart-zm(k))*dzi_t(k)
          particle%xstart      = xstart
          particle%ystart      = ystart
          particle%zstart      = zstart
          particle%tstart      = tstart
          particle%ures        = 0.
          particle%vres        = 0.
          particle%wres        = 0.
          particle%ures_prev   = 0.
          particle%vres_prev   = 0.
          particle%wres_prev   = 0.
          particle%usgs        = 0.
          particle%vsgs        = 0.
          particle%wsgs        = 0.
          particle%usgs_prev   = 0.
          particle%vsgs_prev   = 0.
          particle%wsgs_prev   = 0.
          particle%x_prev      = particle%x
          particle%y_prev      = particle%y
          particle%z_prev      = particle%z
          particle%partstep    = 0
          !particle%sigma2_sgs = epsilon(particle%sigma2_sgs)
        end if
      end if
    end do

    ipunique       = 1
    ipx            = 2
    ipy            = 3
    ipz            = 4
    ipxstart       = 5
    ipystart       = 6
    ipzstart       = 7
    iptsart        = 8
    ipures         = 9
    ipvres         = 10
    ipwres         = 11
    ipures_prev    = 12
    ipvres_prev    = 13
    ipwres_prev    = 14
    ipusgs         = 15
    ipvsgs         = 16
    ipwsgs         = 17
    ipusgs_prev    = 18
    ipvsgs_prev    = 19
    ipwsgs_prev    = 20
    ipartstep      = 21
    ipxprev        = 22
    ipyprev        = 23
    ipzprev        = 24
    nrpartvar      = ipzprev

    !if (lpartsgs) then
    !  ipuresprev = nrpartvar+1
    !  ipvresprev = nrpartvar+2
    !  ipwresprev = nrpartvar+3
    !  ipxprev = nrpartvar+4
    !  ipyprev = nrpartvar+5
    !  ipzprev = nrpartvar+6
    !  nrpartvar=nrpartvar+6
    !end if

    !if (lpartsgs) then
    !  ipsigma2_sgs = nrpartvar + 1
    !  nrpartvar = nrpartvar + 1
    !end if

    ! Initialize particle dump to NetCDF
    if(lpartdump) call initparticledump

    ! Set first dump times
    tnextdump = frqpartdump

    close(ifinput)

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
