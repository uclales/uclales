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

!----------------------------------------------------------------------------
! To-do:
!   - merge interpolation functions -> are all pretty similar
!   - put lpartsgs switch in particle communication
!   - fix statistics -> possible error when particles move from proc2proc
!----------------------------------------------------------------------------


module modparticles
  !--------------------------------------------------------------------------
  ! module modparticles: Langrangian particle tracking, ad(o/a)pted from DALES
  !--------------------------------------------------------------------------
  implicit none
  PUBLIC :: init_particles, particles, exit_particles, initparticledump, initparticlestat

  ! For/from namelist  
  logical            :: lpartic      = .false.        ! Switch for particles
  logical            :: lpartsgs     = .false.        ! Switch for particle subgrid scheme
  logical            :: lpartstat    = .false.        ! Switch for particle statistics
  real               :: frqpartstat  =  3600.         ! Time interval for statistics writing
  real               :: avpartstat   =  60.           ! Averaging time before writing stats
  logical            :: lpartdump    = .false.        ! Switch for particle dump
  logical            :: lpartdumpui  = .false.        ! Switch for writing velocities to dump
  logical            :: lpartdumpsgs = .false.        ! Switch for writing subgrid velocities to dump
  real               :: frqpartdump  =  3600          ! Time interval for particle dump

  character(30)      :: startfile
  integer            :: ifinput     = 1
  integer            :: np
  integer            :: tnextdump, tnextstat

  ! Particle structure
  type :: particle_record
    real             :: unique, tstart
    integer          :: partstep
    real             :: x, xstart, x_prev, ures, ures_prev, usgs, usgs_prev
    real             :: y, ystart, y_prev, vres, vres_prev, vsgs, vsgs_prev
    real             :: z, zstart, z_prev, wres, wres_prev, wsgs, wsgs_prev
    real             :: sigma2_sgs
    type (particle_record), pointer :: next,prev
  end type

  integer            :: nplisted
  type (particle_record), pointer :: head, tail

  integer            :: ipunique, ipx, ipy, ipz, ipxstart, ipystart, ipzstart, iptsart, ipxprev, ipyprev, ipzprev
  integer            :: ipures, ipvres, ipwres, ipusgs, ipvsgs, ipwsgs, ipusgs_prev, ipvsgs_prev, ipwsgs_prev
  integer            :: ipures_prev, ipvres_prev, ipwres_prev, ipartstep, nrpartvar, isigma2_sgs

  ! Statistics and particle dump
  integer            :: ncpartid, ncpartrec             ! Particle dump
  integer            :: ncpartstatid, ncpartstatrec     ! Particle statistics
  integer            :: nstatsamp
  
  ! Arrays for local and domain averaged values
  real, allocatable, dimension(:)     :: npartprof,npartprofl, &
                                         uprof, uprofl, &
                                         vprof, vprofl, &
                                         wprof, wprofl, &
                                         fsprof, fsprofl, &
                                         eprof, eprofl

  real, allocatable, dimension(:,:,:) :: sgse,dthvdz
  real, parameter                     :: minsgse = 5e-5
  real, allocatable, dimension(:)     :: fs
  integer (KIND=selected_int_kind(10)):: idum = -12345
  real, parameter                     :: C0   = 6.
  real                                :: ceps, labda

  real                                :: dsigma2dx = 0, dsigma2dy = 0, dsigma2dz = 0, &
                                         dsigma2dt = 0, sigma2l = 0, epsl = 0, fsl = 0

contains
  !
  !--------------------------------------------------------------------------
  ! subroutine particles: Main routine, called every RK3-step
  !--------------------------------------------------------------------------
  !
  subroutine particles(time)
    use grid, only : dxi, dyi, nstep, dzi_t, dt, nzp
    !use mpi_interface, only : nxg, nyg
    implicit none
    real, intent(in)               :: time
    type (particle_record), pointer:: particle
    integer                        :: muz,mvz,mwz,k

    if ( np < 1 ) return   ! Just to be sure..

    if (lpartsgs .and. nstep == 1) then
      call calc_fields      ! Estimates SGS-TKE
      call fsubgrid(time)   ! Calculates fraction SGS-TKE / TOTAL-TKE
    end if

    particle => head
    do while( associated(particle) )
      if (  time - particle%tstart >= 0 ) then
        particle%partstep = particle%partstep + 1

        ! Interpolation of the velocity field
        particle%ures = velocity_ures(particle%x,particle%y,particle%z) * dxi
        particle%vres = velocity_vres(particle%x,particle%y,particle%z) * dyi
        particle%wres = velocity_wres(particle%x,particle%y,particle%z) * dzi_t(floor(particle%z))

        ! Subgrid velocities 
        if (lpartsgs .and. nstep == 1) then
          call prep_ui_sgs(particle)
          particle%usgs   = tend_usgs(particle) * dxi 
          particle%vsgs   = tend_vsgs(particle) * dyi
          particle%wsgs   = tend_wsgs(particle) * dzi_t(floor(particle%z))

          !if(abs(particle%usgs / dxi) > maxsgsu) then
          !  maxsgsu = abs(particle%usgs) / dxi
          !  muz     = particle%z
          !end if
          !if(abs(particle%vsgs / dyi) > maxsgsv) then
          !  maxsgsv = abs(particle%vsgs) / dyi
          !  mvz     = particle%z
          !end if
          !if(abs(particle%wsgs / dzi_t(floor(particle%z))) > maxsgsw) then
          !  maxsgsw = abs(particle%wsgs / dzi_t(floor(particle%z)))
          !  mwz     = particle%z
          !end if

        end if
      end if

    particle => particle%next
    end do
   
    ! Time integration
    particle => head
    do while( associated(particle) )
      if ( time - particle%tstart >= 0 ) then
        call rk3(particle)
      end if
    particle => particle%next
    end do

    ! Statistics
    if (nstep==3) then
      call checkdiv
      
      ! Particle dump
      if((time + dt > tnextdump) .and. lpartdump) then
        call particledump(time)
        tnextdump = tnextdump + frqpartdump
      end if

      ! Average statistics
      if((time + dt > tnextstat - avpartstat) .and. lpartstat) then
        call particlestat(.false.,time)
      end if

      ! Write statistics
      if((time + dt > tnextstat) .and. lpartstat) then
        call particlestat(.false.,time)
        call particlestat(.true.,time)
      end if
    end if

    !Exchange particle to other processors
    call partcomm

  end subroutine particles

  !
  !--------------------------------------------------------------------------
  ! subroutine fsubgrid : calculates fs (contribution sgs turbulence to
  !   total turbulence)
  !--------------------------------------------------------------------------
  !
  subroutine fsubgrid(time)
    use grid, only : a_up, a_vp, a_wp, nzp, nxp, nyp, dt, nstep, zt
    use mpi_interface, only : ierror, mpi_double_precision, mpi_sum, mpi_comm_world, nxg, nyg
    implicit none
  
    real, intent(in)               :: time
    integer    :: k
    real, allocatable, dimension(:)   :: &
       u_avl, v_avl, u2_avl, v2_avl, w2_avl, sgse_avl,     &
       u_av,  v_av,  u2_av,  v2_av,  w2_av,  sgse_av, e_res

    !integer :: ifoutput = 1

    allocate(u_avl(nzp), v_avl(nzp), u2_avl(nzp), v2_avl(nzp), w2_avl(nzp), sgse_avl(nzp),   &
             u_av(nzp),  v_av(nzp),  u2_av(nzp),  v2_av(nzp),  w2_av(nzp),  sgse_av(nzp),    &
             e_res(nzp))

    do k=1,nzp
      u_avl(k)    = sum(a_up(k,3:nxp-2,3:nyp-2))
      v_avl(k)    = sum(a_vp(k,3:nxp-2,3:nyp-2))
      w2_avl(k)   = sum(a_wp(k,3:nxp-2,3:nyp-2)**2.)
      sgse_avl(k) = sum(sgse(k,3:nxp-2,3:nyp-2))
    end do 

    call mpi_allreduce(u_avl,u_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(v_avl,v_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(w2_avl,w2_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(sgse_avl,sgse_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
  
    u_av        = u_av    / (nxg * nyg)  
    v_av        = v_av    / (nxg * nyg)     
 
    do k=1,nzp
      u2_avl(k) = sum((a_up(k,3:nxp-2,3:nyp-2) - u_av(k))**2.)
      v2_avl(k) = sum((a_vp(k,3:nxp-2,3:nyp-2) - v_av(k))**2.)
    end do 

    call mpi_allreduce(u2_avl,u2_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(v2_avl,v2_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

    u2_av       = u2_av   / (nxg * nyg)    
    v2_av       = v2_av   / (nxg * nyg)    
    w2_av       = w2_av   / (nxg * nyg)    
    sgse_av     = sgse_av / (nxg * nyg) 
 
    do k=2, nzp
      e_res(k)  = 0.5 * (u2_av(k) + v2_av(k) + 0.5*(w2_av(k) + w2_av(k-1)))
      fs(k)     = 1. !sgse_av(k) / (sgse_av(k) + e_res(k))
    end do

    e_res(1)    = -e_res(2)         ! e_res -> 0 at surface
    fs(1)       = 1. + (1.-fs(2))   ! fs    -> 1 at surface

    ! Raw statistics dump
    !if((time + dt > tnextstat) .and. lpartstat .and. nstep==1) then
    !  open(ifoutput,file='rawstat',position='append',action='write')
    !  write(ifoutput,'(A2,F10.2)') '# ',time
    !
    !  do k=1, nzp
    !    write(ifoutput,'(I10,9E15.6)') k,zt(k),u_av(k),u2_av(k),v_av(k),v2_av(k),w2_av(k),e_res(k),sgse_av(k),fs(k)
    !  end do

    !  close(ifoutput) 
    !end if

    deallocate(u_avl,v_avl,u2_avl,v2_avl,w2_avl,sgse_avl,u_av,v_av,u2_av,v2_av,w2_av,sgse_av,e_res)

  end subroutine fsubgrid

  !
  !--------------------------------------------------------------------------
  ! subroutine calc_fields : calculated additional 3d fields needed for 
  !   each particle: sgs-tke, gradient thetav, ...
  !--------------------------------------------------------------------------
  !
  subroutine calc_fields
    use grid, only             : a_km, nzp, zt, dxi, dyi, nxp, nyp, zm, a_rp, a_tp, dzi_t
    use defs, only             : pi, vonk
    use mpi_interface, only    : nxg, nyg
    implicit none
    
    real                       :: labda0, sz1,thvp,thvm
    integer                    :: i,j,k
    real, parameter            :: alpha = 1.5
    real, parameter            :: cf    = 2.5
    real, parameter            :: csx   = 0.23

    labda0 = (zm(2)/dxi/dyi)**(1./3.)

    !do j=1,nyp
    !   do i=1,nxp
    !      do k=2,nzp-1
    !        labda            = (1. / ((1. / labda0**2.) + (1. / (0.4 * (zt(k) + 0.001)**2.))))**0.5
    !        ceps             = 0.19 + 0.51 * (labda / labda0)
    !        sgse(k,i,j)      = (0.5*(a_km(k,i,j) + a_km(k-1,i,j)) / ((labda * (cf / (2. * pi)) * (1.5 * alpha)**(-1.5))))**2.  
    !      end do
    !      sgse(1,i,j)     = sgse(2,i,j)
    !      sgse(nzp,i,j)   = sgse(nzp-1,i,j)
    !   end do
    !end do

    do j=1,nyp
       do i=1,nxp
          do k=2,nzp-1
            ! 1. SGS-TKE
            ! ---------------------------
            labda             = labda0 !(1. / ((1. / labda0**2.) + (1. / (0.4 * (zt(k) + 0.001)**2.))))**0.5
            sz1               = (0.5*(a_km(k,i,j) + a_km(k-1,i,j)))**2.
            sgse(k,i,j)       = max(sz1 / (labda*pi*(csx**2))**2,minsgse)
 
            !sz1              = 1./sqrt(1./(labda0 * csx)**2 + 1./(zt(k) * vonk+0.001)**2)
            !ceps             = 0.19 + 0.51 * (labda / labda0)
            !sgse(k,i,j)      = (0.5*(a_km(k,i,j) + a_km(k-1,i,j)) / ((labda * (cf / (2. * pi)) * (1.5 * alpha)**(-1.5))))**2.  

            ! 2. Vertical gradient thetav
            ! ---------------------------
            thvp = (0.5*(a_tp(k,i,j) + a_tp(k+1,i,j))) * (1. + 0.61 * (0.5*(a_rp(k,i,j) + a_rp(k+1,i,j))))
            thvm = (0.5*(a_tp(k,i,j) + a_tp(k-1,i,j))) * (1. + 0.61 * (0.5*(a_rp(k,i,j) + a_rp(k-1,i,j))))
            dthvdz(k,i,j) = (thvp - thvm) * dzi_t(k)
          end do
          sgse(1,i,j)     = - sgse(2,i,j)
          sgse(nzp,i,j)   = sgse(nzp-1,i,j)
          dthvdz(1,i,j)   = dthvdz(2,i,j)
          dthvdz(nzp,i,j) = dthvdz(nzp-1,i,j)
       end do
    end do

  end subroutine calc_fields

  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_ures : trilinear interpolation of u-component
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

    deltaz = ((zm(floor(z)) + (z - floor(z)) / dzi_t(floor(z))) - zt(zbottom)) * dzi_m(zbottom)
    velocity_ures =  (1-deltaz) * (1-deltay) * (1-deltax) * sign * a_up(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) * sign * a_up(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) * sign * a_up(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) * sign * a_up(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *        a_up(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *        a_up(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *        a_up(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *        a_up(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function velocity_ures

  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_vres : trilinear interpolation of v-component
  !--------------------------------------------------------------------------
  !
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

    deltaz = ((zm(floor(z)) + (z - floor(z)) / dzi_t(floor(z))) - zt(zbottom)) * dzi_m(zbottom)
    velocity_vres =  (1-deltaz) * (1-deltay) * (1-deltax) * sign * a_vp(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) * sign * a_vp(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) * sign * a_vp(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) * sign * a_vp(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *        a_vp(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *        a_vp(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *        a_vp(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *        a_vp(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function velocity_vres
  
  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_wres : trilinear interpolation of w-component
  !--------------------------------------------------------------------------
  !
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

  !
  !--------------------------------------------------------------------------
  ! functions sca2part : trilinear interpolation of properties at scalar
  !   point (center of grid cell: T,q,tke,eps,...) to particle
  !--------------------------------------------------------------------------
  !
  function sca2part(x,y,z,input,lim)
    implicit none
    real, intent(in)                    :: x, y, z
    real, dimension(:,:,:), intent(in)  :: input
    logical, intent(in)                 :: lim

    integer  :: xbottom, ybottom, zbottom
    real     :: sca2part, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z + 0.5)
    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z + 0.5 - zbottom
    
    sca2part      =  (1-deltaz) * (1-deltay) * (1-deltax) *  input(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) *  input(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) *  input(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) *  input(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *  input(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *  input(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *  input(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *  input(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1
    
    if(lim) sca2part = max(sca2part,minsgse) 

  end function sca2part

  !
  !--------------------------------------------------------------------------
  ! subroutine prep_ui_sgs : prepares some variables for subgrid velocity
  !--------------------------------------------------------------------------
  ! 
  subroutine prep_ui_sgs(particle)
    use grid, only : dxi, dyi, dt, a_scr7, zm, zt, dzi_t, dzi_m, a_tp, a_rp, nxp,nyp, th00
    use defs, only : g
    implicit none

    real     :: deltaz
    integer  :: zbottom
    TYPE (particle_record), POINTER:: particle

    real     :: el,N,vgradthv

    zbottom      = floor(particle%z + 0.5)
    deltaz       = ((zm(floor(particle%z)) + (particle%z - floor(particle%z)) / dzi_t(floor(particle%z))) - zt(zbottom)) * dzi_m(zbottom)
    fsl          = (1-deltaz) * fs(zbottom) + deltaz * fs(zbottom+1)
    
    sigma2l      = sca2part(particle%x,particle%y,particle%z,sgse,.true.) * (2./3.)
    el           = sca2part(particle%x,particle%y,particle%z,sgse,.true.)
    !epsl         = sca2part(particle%x,particle%y,particle%z,a_scr7,.false.)
    
    dsigma2dx    = (2./3.) * (sca2part(particle%x+0.5,particle%y,particle%z,sgse,.true.) - &
                              sca2part(particle%x-0.5,particle%y,particle%z,sgse,.true.)) * dxi
    dsigma2dy    = (2./3.) * (sca2part(particle%x,particle%y+0.5,particle%z,sgse,.true.) - &
                              sca2part(particle%x,particle%y-0.5,particle%z,sgse,.true.)) * dyi
    dsigma2dz    = (2./3.) * (sca2part(particle%x,particle%y,particle%z+0.5,sgse,.true.) - &
                              sca2part(particle%x,particle%y,particle%z-0.5,sgse,.true.)) * dzi_t(floor(particle%z))
    dsigma2dt    = (sigma2l - particle%sigma2_sgs) / dt
    !dsigma2dt    = sign(1.,dsigma2dt) * min(abs(dsigma2dt),1. / dt)     ! Limit dsigma2dt term
    particle%sigma2_sgs = sigma2l

    N            = (g/th00) * sca2part(particle%x,particle%y,particle%z,dthvdz,.false.)
    labda        = min(labda,(0.76 * (sqrt(1.5 * sigma2l)/N)))
    ceps         = 0.19 + 0.51 * labda  

    !if(N .gt. 0.0001) then
    !  fsl = 1.
    !end if

  end subroutine prep_ui_sgs

  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_usgs : subgrid velocity particle
  !   NOTE: a_scr7 contains the dissipation rate (eps)
  !--------------------------------------------------------------------------
  !
  function tend_usgs(particle)
    use grid, only : dxi, dt
    implicit none

    real :: t1, t2, t3, tend_usgs
    TYPE (particle_record), POINTER:: particle

    t1        = (-0.75 * fsl * C0 * (particle%usgs / dxi) * (ceps/labda) * (1.5 * sigma2l)**0.5) * dt
    t2        = ( 0.5 * ((1. / sigma2l) * dsigma2dt * (particle%usgs / dxi) + dsigma2dx)) * dt
    t3        = ((fsl * C0 * (ceps/labda) * (1.5*sigma2l)**(1.5))**0.5 * xi(idum))

    ! 1. dissipation subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t1 + particle%usgs / dxi) .ne. sign(1.,particle%usgs / dxi)) t1 = - particle%usgs / dxi
    ! 2. decrease subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t2 + particle%usgs / dxi) .ne. sign(1.,particle%usgs / dxi)) t2 = - particle%usgs / dxi
    ! 3. Terms 1. and 2. combined may not exceed sugbrid velocity itself
    if(sign(1.,t1 + t2 + particle%usgs / dxi) .ne. sign(1.,particle%usgs / dxi)) then
      t1 = -0.5 * particle%usgs / dxi
      t2 = -0.5 * particle%usgs / dxi
    end if
   
    tend_usgs = (particle%usgs / dxi) + (t1 + t2 + t3)  

  end function tend_usgs

  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_vsgs : subgrid velocity particle
  !   NOTE: a_scr7 contains the dissipation rate (eps)
  !--------------------------------------------------------------------------
  !
  function tend_vsgs(particle)
    use grid, only : dyi, dt
    implicit none

    real :: t1, t2, t3, tend_vsgs
    TYPE (particle_record), POINTER:: particle

    t1        = (-0.75 * fsl * C0 * (particle%vsgs / dyi) * (ceps/labda) * (1.5 * sigma2l)**0.5) * dt
    t2        = ( 0.5 * ((1. / sigma2l) * dsigma2dt * (particle%vsgs / dyi) + dsigma2dy)) * dt
    t3        = ((fsl * C0 * (ceps/labda) * (1.5*sigma2l)**(1.5))**0.5 * xi(idum))

    ! 1. dissipation subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t1 + particle%vsgs / dyi) .ne. sign(1.,particle%vsgs / dyi)) t1 = - particle%vsgs / dyi
    ! 2. decrease subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t2 + particle%vsgs / dyi) .ne. sign(1.,particle%vsgs / dyi)) t2 = - particle%vsgs / dyi
    ! 3. Terms 1. and 2. combined may not exceed sugbrid velocity
    if(sign(1.,t1 + t2 + particle%vsgs / dyi) .ne. sign(1.,particle%vsgs / dyi)) then
      t1 = -0.5 * particle%vsgs / dyi
      t2 = -0.5 * particle%vsgs / dyi
    end if
 
    tend_vsgs = (particle%vsgs / dyi) + (t1 + t2 + t3) 

  end function tend_vsgs 

  !
  !--------------------------------------------------------------------------
  ! subroutine velocity_wsgs : subgrid velocity particle
  !   NOTE: a_scr7 contains the dissipation rate (eps)
  !--------------------------------------------------------------------------
  !
  function tend_wsgs(particle)
    use grid, only : dzi_t, dt
    implicit none

    real :: t1, t2, t3, tend_wsgs, dzi
    TYPE (particle_record), POINTER:: particle

    dzi        = dzi_t(floor(particle%z))
    t1        = (-0.75 * fsl * C0 * (particle%wsgs / dzi) * (ceps/labda) * (1.5 * sigma2l)**0.5) * dt
    t2        = ( 0.5 * ((1. / sigma2l) * dsigma2dt * (particle%wsgs / dzi) + dsigma2dz)) * dt
    t3        = ((fsl * C0 * (ceps/labda) * (1.5*sigma2l)**(1.5))**0.5 * xi(idum))

    ! 1. dissipation subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t1 + particle%wsgs / dzi) .ne. sign(1.,particle%wsgs / dzi)) t1 = - particle%wsgs / dzi
    ! 2. decrease subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t2 + particle%wsgs / dzi) .ne. sign(1.,particle%wsgs / dzi)) t2 = - particle%wsgs / dzi
    ! 3. Terms 1. and 2. combined may not exceed sugbrid velocity
    if(sign(1.,t1 + t2 + particle%wsgs / dzi) .ne. sign(1.,particle%wsgs / dzi)) then
      t1 = -0.5 * particle%wsgs / dzi
      t2 = -0.5 * particle%wsgs / dzi
    end if
 
    tend_wsgs = (particle%wsgs / dzi) + (t1 + t2 + t3) 

  end function tend_wsgs

  !
  !--------------------------------------------------------------------------
  ! subroutine rk3 : Third order Runge-Kutta scheme for spatial integration
  !--------------------------------------------------------------------------
  !
  subroutine rk3(particle)
    use grid, only : rkalpha, rkbeta, nstep, dt, dxi, dyi, dzi_t
    implicit none
    TYPE (particle_record), POINTER:: particle

    particle%x   = particle%x + rkalpha(nstep) * (particle%ures+particle%usgs) * dt + rkbeta(nstep) * (particle%ures_prev + particle%usgs) * dt
    particle%y   = particle%y + rkalpha(nstep) * (particle%vres+particle%vsgs) * dt + rkbeta(nstep) * (particle%vres_prev + particle%vsgs) * dt
    particle%z   = particle%z + rkalpha(nstep) * (particle%wres+particle%wsgs) * dt + rkbeta(nstep) * (particle%wres_prev + particle%wsgs) * dt

    particle%ures_prev = particle%ures
    particle%vres_prev = particle%vres
    particle%wres_prev = particle%wres

    particle%usgs_prev = particle%usgs
    particle%vsgs_prev = particle%vsgs
    particle%wsgs_prev = particle%wsgs

    call checkbound(particle)
   
   if ( nstep==3 ) then
      particle%ures_prev   = 0.
      particle%vres_prev   = 0.
      particle%wres_prev   = 0.
      particle%usgs_prev   = 0.
      particle%vsgs_prev   = 0.
      particle%wsgs_prev   = 0.
    end if

  end subroutine rk3

  !
  !--------------------------------------------------------------------------
  ! subroutine partcomm : send and receives particles between different
  !   processes and handles cyclic boundaries
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
  ! subroutine partbuffer : used to create send/receive buffer for 
  !   routine partcomm()
  !--------------------------------------------------------------------------
  !
  subroutine partbuffer(particle, buffer, n, send)
    implicit none

    logical,intent(in)                :: send
    integer,intent(in)                :: n
    real,dimension(n+1:n+nrpartvar)   :: buffer
    TYPE (particle_record), POINTER:: particle

    if (send) then
      buffer(n+ipunique)        = particle%unique
      buffer(n+ipx)             = particle%x
      buffer(n+ipy)             = particle%y
      buffer(n+ipz)             = particle%z
      buffer(n+ipures)          = particle%ures
      buffer(n+ipvres)          = particle%vres
      buffer(n+ipwres)          = particle%wres
      buffer(n+ipures_prev)     = particle%ures_prev
      buffer(n+ipvres_prev)     = particle%vres_prev
      buffer(n+ipwres_prev)     = particle%wres_prev
      buffer(n+ipusgs)          = particle%usgs
      buffer(n+ipvsgs)          = particle%vsgs
      buffer(n+ipwsgs)          = particle%wsgs
      buffer(n+ipusgs_prev)     = particle%usgs_prev
      buffer(n+ipvsgs_prev)     = particle%vsgs_prev
      buffer(n+ipwsgs_prev)     = particle%wsgs_prev
      buffer(n+ipxstart)        = particle%xstart
      buffer(n+ipystart)        = particle%ystart
      buffer(n+ipzstart)        = particle%zstart
      buffer(n+iptsart)         = particle%tstart
      buffer(n+ipartstep)       = particle%partstep
      buffer(n+ipxprev)         = particle%x_prev
      buffer(n+ipyprev)         = particle%y_prev
      buffer(n+ipzprev)         = particle%z_prev
      buffer(n+isigma2_sgs)     = particle%sigma2_sgs
    else
      particle%unique           = buffer(n+ipunique)
      particle%x                = buffer(n+ipx)
      particle%y                = buffer(n+ipy)
      particle%z                = buffer(n+ipz)
      particle%ures             = buffer(n+ipures)
      particle%vres             = buffer(n+ipvres)
      particle%wres             = buffer(n+ipwres)
      particle%ures_prev        = buffer(n+ipures_prev)
      particle%vres_prev        = buffer(n+ipvres_prev)
      particle%wres_prev        = buffer(n+ipwres_prev)
      particle%usgs             = buffer(n+ipusgs)
      particle%vsgs             = buffer(n+ipvsgs)
      particle%wsgs             = buffer(n+ipwsgs)
      particle%usgs_prev        = buffer(n+ipusgs_prev)
      particle%vsgs_prev        = buffer(n+ipvsgs_prev)
      particle%wsgs_prev        = buffer(n+ipwsgs_prev)
      particle%xstart           = buffer(n+ipxstart)
      particle%ystart           = buffer(n+ipystart)
      particle%zstart           = buffer(n+ipzstart)
      particle%tstart           = buffer(n+iptsart)
      particle%partstep         = buffer(n+ipartstep)
      particle%x_prev           = buffer(n+ipxprev)
      particle%y_prev           = buffer(n+ipyprev)
      particle%z_prev           = buffer(n+ipzprev)
      particle%sigma2_sgs       = buffer(n+isigma2_sgs)
    end if

  end subroutine partbuffer

  !
  !--------------------------------------------------------------------------
  ! subroutine particlestat : Time averages slab-averaged statistics and
  !   writes them to NetCDF
  !--------------------------------------------------------------------------
  !
  subroutine particlestat(dowrite,time)
    use mpi_interface, only : mpi_comm_world, myid, mpi_double_precision, mpi_sum, ierror, nxprocs, nyprocs, nxg, nyg
    use modnetcdf,     only : writevar_nc, fillvalue_double
    use grid,          only : tname,nzp,dxi,dyi,dzi_t,nxp,nyp
    implicit none

    logical, intent(in)     :: dowrite
    real, intent(in)        :: time
    integer                 :: k
    type (particle_record), pointer:: particle

    ! Time averaging step
    if(.not. dowrite) then
      particle => head
      do while(associated(particle))
        k               = floor(particle%z) + 1
        npartprofl(k)   = npartprofl(k) + 1
        uprofl(k)       = uprofl(k)     + (particle%ures / dxi)
        vprofl(k)       = vprofl(k)     + (particle%vres / dyi)
        wprofl(k)       = wprofl(k)     + (particle%wres / dzi_t(floor(particle%z)))
        particle => particle%next
      end do 

      if(lpartsgs) then
        do k = 1,nzp 
          fsprofl(k)      = fsprofl(k)    + fs(k) / (nxprocs * nyprocs)
          eprofl(k)       = eprofl(k)     + sum(sgse(k,3:nxp-2,3:nyp-2)) / (nxg * nyg)
        end do     
      end if

      nstatsamp = nstatsamp + 1
    end if

    ! Write to NetCDF
    if(dowrite) then
      do k = 1,nzp
        if(npartprofl(k) > 0) then 
          npartprofl(k) = npartprofl(k) /  nstatsamp
          uprofl(k)     = uprofl(k)     / (nstatsamp * npartprofl(k))
          vprofl(k)     = vprofl(k)     / (nstatsamp * npartprofl(k))
          wprofl(k)     = wprofl(k)     / (nstatsamp * npartprofl(k))
        else
          npartprofl(k) = 0.
          uprofl(k)     = 1e9
          vprofl(k)     = 1e9
          wprofl(k)     = 1e9
        end if
      end do       

      if(lpartsgs) then
        fsprofl  = fsprofl    / nstatsamp
        eprofl   = eprofl     / nstatsamp     
      end if

      call mpi_allreduce(npartprofl,npartprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(uprofl,uprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(vprofl,vprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(wprofl,wprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      
      if(lpartsgs) then
        call mpi_allreduce(fsprofl,fsprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(eprofl,eprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      end if

      if(myid == 0) then
        call writevar_nc(ncpartstatid,tname,time,ncpartstatrec)
        call writevar_nc(ncpartstatid,'np',npartprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'u',uprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'v',vprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'w',wprof,ncpartstatrec)
        if(lpartsgs) then
          call writevar_nc(ncpartstatid,'fs',fsprof,ncpartstatrec)
          call writevar_nc(ncpartstatid,'e',eprof,ncpartstatrec)
        end if
      end if

      npartprof  = 0
      npartprofl = 0
      uprof      = 0
      uprofl     = 0
      vprof      = 0
      vprofl     = 0
      wprof      = 0
      wprofl     = 0
      fsprof     = 0
      fsprofl    = 0
      eprof      = 0
      eprofl     = 0
      nstatsamp  = 0
      tnextstat  = tnextstat + frqpartstat
    end if

  end subroutine particlestat

  !
  !--------------------------------------------------------------------------
  ! subroutine particledump : communicates all particles to process #0 and
  !   writes it to NetCDF(4)
  !--------------------------------------------------------------------------
  !
  subroutine particledump(time)
    use mpi_interface, only : mpi_comm_world, myid, mpi_integer, mpi_double_precision, ierror, nxprocs, nyprocs, mpi_status_size, wrxid, wryid, nxg, nyg
    use grid,          only : tname, deltax, deltay, dzi_t, zm
    use modnetcdf,     only : writevar_nc, fillvalue_double
    implicit none

    real, intent(in)                     :: time
    type (particle_record), pointer:: particle
    integer                              :: nlocal, ii, i, pid, partid
    integer, allocatable, dimension(:)   :: nremote
    integer                              :: status(mpi_status_size)
    integer                              :: nvar
    real, allocatable, dimension(:)      :: sendbuff, recvbuff
    real, allocatable, dimension(:,:)    :: particles_merged

    nvar = 4                            ! id,x,y,z
    if(lpartdumpui)  nvar = nvar + 3    ! u,v,w
    if(lpartdumpsgs) nvar = nvar + 3    ! us,vs,ws

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
      sendbuff(ii+1) = (wrxid * (nxg / nxprocs) + particle%x - 3) * deltax
      sendbuff(ii+2) = (wryid * (nyg / nyprocs) + particle%y - 3) * deltay
      sendbuff(ii+3) = zm(floor(particle%z)) + (particle%z-floor(particle%z)) / dzi_t(floor(particle%z))
      if(lpartdumpui) then
        sendbuff(ii+4) = particle%ures * deltax
        sendbuff(ii+5) = particle%vres * deltay
        sendbuff(ii+6) = particle%wres / dzi_t(floor(particle%z))
      end if
      if(lpartdumpsgs) then
        sendbuff(ii+7) = particle%usgs * deltax
        sendbuff(ii+8) = particle%vsgs * deltay
        sendbuff(ii+9) = particle%wsgs / dzi_t(floor(particle%z))
      end if
      ii = ii + nvar
      particle => particle%next
    end do

    ! Dont send when main process
    if(myid .ne. 0) call mpi_send(sendbuff,nvar*nlocal,mpi_double_precision,0,myid,mpi_comm_world,ierror)

    if(myid .eq. 0) then
      allocate(particles_merged(np,nvar-1))

      ! Add local particles
      ii = 1
      do i=1,nlocal
        partid = int(sendbuff(ii))
        particles_merged(partid,1) = sendbuff(ii+1)
        particles_merged(partid,2) = sendbuff(ii+2)
        particles_merged(partid,3) = sendbuff(ii+3)
        if(lpartdumpui) then
          particles_merged(partid,4) = sendbuff(ii+4)
          particles_merged(partid,5) = sendbuff(ii+5)
          particles_merged(partid,6) = sendbuff(ii+6)
        end if
        if(lpartdumpsgs) then
          particles_merged(partid,7) = sendbuff(ii+7)
          particles_merged(partid,8) = sendbuff(ii+8)
          particles_merged(partid,9) = sendbuff(ii+9)
        end if
        ii = ii + nvar 
      end do 

      ! Add remote particles
      do pid = 1,(nxprocs*nyprocs)-1
        allocate(recvbuff(nremote(pid)*nvar))
        call mpi_recv(recvbuff,nremote(pid)*nvar,mpi_double_precision,pid,pid,mpi_comm_world,status,ierror)   
        ii = 1
        do i=1,nremote(pid)
          partid = int(recvbuff(ii))
          particles_merged(partid,1) = recvbuff(ii+1)
          particles_merged(partid,2) = recvbuff(ii+2)
          particles_merged(partid,3) = recvbuff(ii+3)
          if(lpartdumpui) then
            particles_merged(partid,4) = recvbuff(ii+4)
            particles_merged(partid,5) = recvbuff(ii+5)
            particles_merged(partid,6) = recvbuff(ii+6)
          end if
          if(lpartdumpsgs) then
            particles_merged(partid,7) = recvbuff(ii+7)
            particles_merged(partid,8) = recvbuff(ii+8)
            particles_merged(partid,9) = recvbuff(ii+9)
          end if
          ii = ii + nvar 
        end do 
        deallocate(recvbuff)
      end do

      ! Write to NetCDF
      call writevar_nc(ncpartid,tname,time,ncpartrec)
      call writevar_nc(ncpartid,'x',particles_merged(:,1),ncpartrec)
      call writevar_nc(ncpartid,'y',particles_merged(:,2),ncpartrec)
      call writevar_nc(ncpartid,'z',particles_merged(:,3),ncpartrec)
      if(lpartdumpui) then
        call writevar_nc(ncpartid,'u',particles_merged(:,4),ncpartrec)
        call writevar_nc(ncpartid,'v',particles_merged(:,5),ncpartrec)
        call writevar_nc(ncpartid,'w',particles_merged(:,6),ncpartrec)
      end if
      if(lpartdumpsgs) then
        call writevar_nc(ncpartid,'us',particles_merged(:,7),ncpartrec)
        call writevar_nc(ncpartid,'vs',particles_merged(:,8),ncpartrec)
        call writevar_nc(ncpartid,'ws',particles_merged(:,9),ncpartrec)
      end if
    end if
    
    if(myid==0) deallocate(particles_merged)
    deallocate(nremote, sendbuff)
 
  end subroutine particledump


  !
  !--------------------------------------------------------------------------
  ! Quick and dirty (Eulerian) divergence check, only for CPU #0
  !--------------------------------------------------------------------------
  !
  subroutine checkdiv
  use grid, only : dxi, dyi, dzi_t, a_up, a_vp, a_wp, nxp, nyp, nzp, dn0
  use mpi_interface, only : myid
  implicit none
  integer :: i,j,k
  real :: dudx,dvdy,dwdz,div
  real :: divmax,divtot
  real :: dnp,dnm

  div = 0.
  divmax = 0.
  divtot = 0.

  do j=3,nyp-2
    do i=3,nxp-2
      do k=2,nzp-2
        dnp  = 0.5 * (dn0(k) + dn0(k+1))
        dnm  = 0.5 * (dn0(k) + dn0(k-1))
        dudx = (a_up(k,i,j) - a_up(k,i-1,j)) * dxi * dn0(k)
        dvdy = (a_vp(k,i,j) - a_vp(k,i,j-1)) * dyi * dn0(k)
        dwdz = ((a_wp(k,i,j) * dnp) - (a_wp(k-1,i,j) * dnm)) * dzi_t(k)
        div  = dudx + dvdy + dwdz
        divtot = divtot + div
        if(abs(div) > divmax) divmax = abs(div)
      end do
    end do
  end do
 
  !if(myid==0) print*,'   divergence; max=',divmax, ', total=',divtot 

  end subroutine checkdiv

  !
  !--------------------------------------------------------------------------
  ! subroutine checkbound : bounces particles of surface and model top
  !--------------------------------------------------------------------------
  !
  subroutine checkbound(particle)
    use grid, only : nxp, nyp, nzp, zm
    implicit none

    type (particle_record), pointer:: particle

    ! Reflect particles of surface and model top
    if (particle%z >= nzp-1) then
      particle%z = nzp-1-0.0001
      particle%wres = -abs(particle%wres)
      particle%wsgs = -abs(particle%wsgs)
    !elseif(particle%z < 2.0) then
    !  particle%z = 2.01
    !  particle%wres =  abs(particle%wres)
    !  particle%wsgs =  abs(particle%wsgs)
    elseif (particle%z < 1.01) then
      particle%z = abs(particle%z-1.01)+1.01
      particle%wres =  abs(particle%wres)
      particle%wsgs =  abs(particle%wsgs)
    end if

  end subroutine checkbound

  !
  !--------------------------------------------------------------------------
  ! function xi & random : creates component of white Gaussian noise  
  !--------------------------------------------------------------------------
  !
  function xi(idum)
    implicit none
    integer (KIND=selected_int_kind(10)):: idum
    integer (KIND=selected_int_kind(10)):: iset
    real :: xi, fac, gset, rsq, v1, v2
    save iset, gset
    data iset /0/
    
    rsq   = 0.
    v1    = 0.
    v2    = 0.

    if (iset == 0) then
      do while (rsq >= 1 .or. rsq == 0)
        v1        = 2. * random(idum)-1.
        v2        = 2. * random(idum)-1.
        rsq       = v1 * v1 + v2 * v2
      end do
      fac         = sqrt(-2. * log(rsq) / rsq)
      gset        = v1 * fac
      xi          = v2 * fac
      iset        = 1
    else
      xi          = gset
      iset        = 0
    end if
    return
  end function xi

  function random(idum)
    implicit none
    integer, parameter :: ntab = 32
    integer (KIND=selected_int_kind(10)):: idum 
    integer (KIND=selected_int_kind(10)):: ia, im, iq, ir, iv(ntab), iy, ndiv, threshold=1
    real :: random, am, eps1, rnmx

    integer :: j, k
    save iv, iy

    ia     = 16807
    im     = 2147483647
    am     = 1. / real(im)
    iq     = 127773
    ir     = 2836
    ndiv   = 1 +  (im-1)/real(ntab)
    eps1   = 1.2E-7
    rnmx   = 1. - eps1

    data iv /ntab*0/ , iy /0/

    if (idum <= 0 .or. iy == 0 ) then
      idum    = max(idum,threshold)
      do j    = ntab + 8, 1, -1
        k     = idum / real(iq)
        idum  = ia * (idum - k * iq) - ir * k
        if (idum < 0) idum = idum + im
        if (j <= ntab ) iv(j) = idum
      end do
      iy = iv(1)
    end if
    
    k      = idum / real(iq)
    idum   = ia * (idum - k * iq) - ir * k
    if (idum <= 0) idum = idum + im
    j      = 1 + iy / real(ndiv)
    iy     = iv(j)
    iv(j)  = idum
    random = min(am*iy,rnmx)
    return
  end function random


  !--------------------------------------------------------------------------
  !
  ! BELOW: ONLY INIT / EXIT PARTICLES
  !
  !--------------------------------------------------------------------------
  !
  ! subroutine init_particles: initialize particles, reading initial position, 
  ! etc. Called from subroutine initialize (init.f90)
  !--------------------------------------------------------------------------
  !
  subroutine init_particles
    use mpi_interface, only : wrxid, wryid, nxg, nyg, myid, nxprocs, nyprocs
    use grid, only : zm, deltax, deltay, zt,dzi_t, nzp, nxp, nyp
    use grid, only : a_up, a_vp, a_wp

    integer :: k, n, kmax
    real :: tstart, xstart, ystart, zstart, ysizelocal, xsizelocal, firststart
    type (particle_record), pointer:: particle

    xsizelocal = (nxg / nxprocs) * deltax
    ysizelocal = (nyg / nyprocs) * deltay
    kmax = size(zm)

    firststart = 1e9

    np = 0
    startfile = 'partstartpos'
    open(ifinput,file=startfile,status='old',position='rewind',action='read')
    read(ifinput,'(I10.1)') np
    if ( np < 1 ) return
   
    ! clear pointers to head and tail
    nullify(head)
    nullify(tail)
    nplisted = 0

    ! read particles from partstartpos, create linked list
    do n = 1, np
      read(ifinput,*) tstart, xstart, ystart, zstart
      if(floor(xstart / xsizelocal) == wrxid) then
        if(floor(ystart / ysizelocal) == wryid) then
          call add_particle(particle)
          particle%unique         = n !+ myid/1000.0
          particle%x              = (xstart - (float(wrxid) * xsizelocal)) / deltax + 3.  ! +3 here for ghost cells.
          particle%y              = (ystart - (float(wryid) * ysizelocal)) / deltay + 3.  ! +3 here for ghost cells.
          do k=kmax,1,1
            if ( zm(k)<zstart ) exit
          end do
          particle%z              = k + (zstart-zm(k))*dzi_t(k)
          particle%xstart         = xstart
          particle%ystart         = ystart
          particle%zstart         = zstart
          particle%tstart         = tstart
          particle%ures           = 0.
          particle%vres           = 0.
          particle%wres           = 0.
          particle%ures_prev      = 0.
          particle%vres_prev      = 0.
          particle%wres_prev      = 0.
          particle%usgs           = 0.
          particle%vsgs           = 0.
          particle%wsgs           = 0.
          particle%usgs_prev      = 0.
          particle%vsgs_prev      = 0.
          particle%wsgs_prev      = 0.
          particle%x_prev         = particle%x
          particle%y_prev         = particle%y
          particle%z_prev         = particle%z
          particle%partstep       = 0
          particle%sigma2_sgs     = epsilon(particle%sigma2_sgs)

          if(tstart < firststart) firststart = tstart

        end if
      end if
    end do

    ipunique        = 1
    ipx             = 2
    ipy             = 3
    ipz             = 4
    ipxstart        = 5
    ipystart        = 6
    ipzstart        = 7
    iptsart         = 8
    ipures          = 9
    ipvres          = 10
    ipwres          = 11
    ipures_prev     = 12
    ipvres_prev     = 13
    ipwres_prev     = 14
    ipusgs          = 15
    ipvsgs          = 16
    ipwsgs          = 17
    ipusgs_prev     = 18
    ipvsgs_prev     = 19
    ipwsgs_prev     = 20
    ipartstep       = 21
    ipxprev         = 22
    ipyprev         = 23
    ipzprev         = 24
    isigma2_sgs     = 25
    nrpartvar       = isigma2_sgs

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
    !if(lpartdump) call initparticledump

    ! Set first dump times
    tnextdump = firststart
    tnextstat = 0
    nstatsamp = 0
 
    if(lpartstat) allocate(npartprof(nzp),npartprofl(nzp),uprof(nzp),uprofl(nzp),vprof(nzp),vprofl(nzp),wprof(nzp),wprofl(nzp),fsprof(nzp),fsprofl(nzp),eprof(nzp), eprofl(nzp))
    if(lpartsgs)  allocate(sgse(nzp,nxp,nyp),dthvdz(nzp,nxp,nyp),fs(nzp))
    
    close(ifinput)

  end subroutine init_particles
  
  !
  !--------------------------------------------------------------------------
  ! subroutine exit_particles
  !--------------------------------------------------------------------------
  !
  subroutine exit_particles
    use mpi_interface, only : myid
    implicit none

    do while( associated(tail) )
      call delete_particle(tail)
    end do

    if(lpartsgs)  deallocate(sgse,dthvdz)
    if(lpartstat) deallocate(npartprof,npartprofl,uprof,uprofl,vprof,vprofl,wprof,wprofl,fs,fsprof,fsprofl,eprof,eprofl)

    if(myid == 0) print "(//' ',49('-')/,' ',/,'  Lagrangian particles removed.')"
  end subroutine exit_particles

  !
  !--------------------------------------------------------------------------
  ! subroutine initparticledump : creates NetCDF file for particle dump. 
  !   Called from: init.f90
  !--------------------------------------------------------------------------
  !
  subroutine initparticledump(time)
    use modnetcdf,       only : open_nc, addvar_nc
    use grid,            only : nzp, tname, tlongname, tunit, filprf
    use mpi_interface,   only : myid
    use grid,            only : tname, tlongname, tunit
    implicit none

    real, intent(in)                  :: time
    character (40), dimension(2)      :: dimname, dimlongname, dimunit
    real, allocatable, dimension(:,:) :: dimvalues
    integer, dimension(2)             :: dimsize
    integer                           :: k
    integer, parameter                :: precis = 0

    allocate(dimvalues(np,2))

    dimvalues      = 0
    do k=1,np
      dimvalues(k,1) = k
    end do

    dimname(1)     = 'particle'
    dimlongname(1) = 'ID of particle'
    dimunit(1)     = '-'    
    dimsize(1)     = np
    dimname(2)     = tname
    dimlongname(2) = tlongname
    dimunit(2)     = tunit
    dimsize(2)     = 0
 
    if(myid == 0) then
      call open_nc(trim(filprf)//'.particles.nc', ncpartid, ncpartrec, time, .true.)
      call addvar_nc(ncpartid,'x','x-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      call addvar_nc(ncpartid,'y','y-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      call addvar_nc(ncpartid,'z','z-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      if(lpartdumpui) then
        call addvar_nc(ncpartid,'u','resolved u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
        call addvar_nc(ncpartid,'v','resolved v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
        call addvar_nc(ncpartid,'w','resolved w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      end if
      if(lpartdumpsgs) then
        call addvar_nc(ncpartid,'us','subgrid u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
        call addvar_nc(ncpartid,'vs','subgrid v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
        call addvar_nc(ncpartid,'ws','subgrid w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      end if
    end if 
 
  end subroutine initparticledump

  !
  !--------------------------------------------------------------------------
  ! subroutine initparticlestat : creates NetCDF file for particle statistics. 
  !   Called from: init.f90
  !--------------------------------------------------------------------------
  !
  subroutine initparticlestat(time)
    use modnetcdf,       only : open_nc, addvar_nc
    use grid,            only : nzp, tname, tlongname, tunit, filprf
    use mpi_interface,   only : myid
    use grid,            only : tname, tlongname, tunit, zname, zlongname, zunit, zt
    implicit none

    real, intent(in)                  :: time
    character (40), dimension(2)      :: dimname, dimlongname, dimunit
    real, allocatable, dimension(:,:) :: dimvalues
    integer, dimension(2)             :: dimsize

    allocate(dimvalues(nzp,2))
    dimvalues = 0
    dimvalues(1:nzp,1)  = zt(1:nzp)

    dimname(1)     = zname
    dimlongname(1) = zlongname
    dimunit(1)     = zunit
    dimsize(1)     = nzp-2
    dimname(2)     = tname
    dimlongname(2) = tlongname
    dimunit(2)     = tunit
    dimsize(2)     = 0
 
    if(myid == 0) then
      call open_nc(trim(filprf)//'.particlestat.nc', ncpartstatid, ncpartstatrec, time, .true.)
      call addvar_nc(ncpartstatid,'np','Number of particles','-',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'u','resolved u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'v','resolved v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'w','resolved w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues)
      if(lpartsgs) then
        call addvar_nc(ncpartstatid,'fs','subgrid turbulence fraction','-',dimname,dimlongname,dimunit,dimsize,dimvalues)
        call addvar_nc(ncpartstatid,'e','subgrid TKE','m2/s2',dimname,dimlongname,dimunit,dimsize,dimvalues)
      end if
    end if 
 
  end subroutine initparticlestat

  !
  !--------------------------------------------------------------------------
  ! subroutine exitparticledump
  !   Called from: step.f90
  !--------------------------------------------------------------------------
  !
  subroutine exitparticledump
    use modnetcdf,       only : close_nc
    implicit none
    call close_nc(ncpartid)     
  end subroutine exitparticledump

  !
  !--------------------------------------------------------------------------
  ! subroutine exitparticlestat
  !   Called from: step.f90
  !--------------------------------------------------------------------------
  !
  subroutine exitparticlestat
    use modnetcdf,       only : close_nc
    implicit none
    call close_nc(ncpartstatid)     
  end subroutine exitparticlestat

  !
  !--------------------------------------------------------------------------
  ! subroutine init_particles
  !--------------------------------------------------------------------------
  !
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

  !
  !--------------------------------------------------------------------------
  ! subroutine delete_particle
  !--------------------------------------------------------------------------
  !
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
