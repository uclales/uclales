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
! Doxygen:
!> Lagrangian Particle Tracking Module (LPTM). 
!! Tracks massless particles driven by the interpolated velocities from the 
!! Eulerian LES grid.
!>
!! \author Bart van Stratum
!! \todo 


module modparticles
  !--------------------------------------------------------------------------
  ! module modparticles: Langrangian particle tracking, ad(o/a)pted from DALES
  !--------------------------------------------------------------------------
  implicit none
  PUBLIC :: init_particles, particles, exit_particles, initparticledump, initparticlestat, write_particle_hist, particlestat, balanced_particledump

  ! For/from namelist  
  logical            :: lpartic        = .false.        !< Switch for enabling particles
  logical            :: lpartsgs       = .false.        !< Switch for enabling particle subgrid scheme
  logical            :: lrandsurf      = .false.        !< Switch for randomizing lowest grid level(s)
  logical            :: lpartstat      = .false.        !< Switch for enabling particle statistics
  logical            :: lpartdump      = .false.        !< Switch for particle dump
  logical            :: lpartdumpui    = .false.        !< Switch for writing velocities to dump
  logical            :: lpartdumpth    = .false.        !< Switch for writing temperatures (liquid water / virtual potential T) to dump
  logical            :: lpartdumpmr    = .false.        !< Switch for writing moisture (total / liquid (+rain if level==3) water mixing ratio) to dump
  real               :: frqpartdump    =  3600          !< Time interval for particle dump

  character(30)      :: startfile
  integer            :: ifinput        = 1
  integer            :: np      
  real               :: tnextdump
  real               :: randint   = 20.
  real               :: tnextrand = 6e6

  ! Particle structure
  type :: particle_record
    real             :: unique, tstart
    integer          :: partstep
    real             :: x, xstart, ures, ures_prev, usgs, usgs_prev
    real             :: y, ystart, vres, vres_prev, vsgs, vsgs_prev
    real             :: z, zstart, wres, wres_prev, wsgs, wsgs_prev, zprev
    real             :: sigma2_sgs
    type (particle_record), pointer :: next,prev
  end type

  integer            :: nplisted
  type (particle_record), pointer :: head, tail

  integer            :: ipunique, ipx, ipy, ipz, ipzprev, ipxstart, ipystart, ipzstart, iptsart
  integer            :: ipures, ipvres, ipwres, ipures_prev, ipvres_prev, ipwres_prev, ipartstep, nrpartvar
  integer            :: ipusgs, ipvsgs, ipwsgs, ipusgs_prev, ipvsgs_prev, ipwsgs_prev, ipsigma2_sgs

  ! Statistics and particle dump
  integer            :: ncpartid, ncpartrec             ! Particle dump
  integer            :: ncpartstatid, ncpartstatrec     ! Particle statistics
  integer            :: nstatsamp
  
  ! Arrays for local and domain averaged values
  real, allocatable, dimension(:)     :: npartprof,    npartprofl,  &
                                         uprof,        uprofl,      &
                                         vprof,        vprofl,      &
                                         wprof,        wprofl,      &
                                         u2prof,       u2profl,     &
                                         v2prof,       v2profl,     &
                                         w2prof,       w2profl,     &
                                         tkeprof,      tkeprofl,    &
                                         tprof,        tprofl,      &
                                         tvprof,       tvprofl,     &
                                         rtprof,       rtprofl,     &
                                         rlprof,       rlprofl,     &
                                         ccprof,       ccprofl,     &
                                         sigma2prof,   sigma2profl, &
                                         fsprof,       fsprofl

  integer (KIND=selected_int_kind(10)):: idum = -12345

  ! Prognostic sgs-velocity-related variables
  real, allocatable, dimension(:,:,:) :: sgse
  real, allocatable, dimension(:,:,:) :: rese
  real, allocatable, dimension(:)     :: fs
  real, allocatable, dimension(:,:,:) :: fs_local
  real, parameter                     :: minsgse = 5e-5
  real, parameter                     :: C0      = 4.
  real                                :: dsigma2dx = 0, dsigma2dy = 0, dsigma2dz = 0, &
                                         dsigma2dt = 0, sigma2l = 0,   epsl = 0, fsl = 0
  real                                :: ceps, labda, spngl  

  ! Test switches...
  logical                             :: lfsloc = .false.  ! Local or global calculated fs()
  logical                             :: fixedfs = .true.  ! Fix fs at 1.

contains
  !
  !--------------------------------------------------------------------------
  ! Subroutine particles 
  !> Main driver of the LPTM, calls both the velocity 
  !> interpolation from the Eulerian grid and RK3 integration scheme. 
  !> called from: step.f90
  !--------------------------------------------------------------------------
  !
  subroutine particles(time,timmax) 
    use grid, only : dxi, dyi, nstep, dzi_t, dt, nzp, zm, a_km, nxp, nyp, nfpt
    use defs, only : pi
    use mpi_interface, only : myid 
    implicit none
    real, intent(in)               :: time          !< time of simulation, determines the timing of particle dumps to NetCDF
    real, intent(in)               :: timmax        !< end of simulation, required to write particle dump at last timestep
    type (particle_record), pointer:: particle

    if (lpartsgs .and. nstep == 1) then
      call calc_sgstke                ! Estimates SGS-TKE
      call fsubgrid                   ! Calculates bulk fraction SGS-TKE / TOTAL-TKE
      if(lfsloc) then   ! only needed when using local fs
        call calc_restke              ! Calculated Resolved TKE
        call fsubgrid_local           ! Calculated local   "     "     "      "
      end if
    end if

    ! Randomize particles lowest grid level
    if (lrandsurf .and. nstep==1 .and. time > tnextrand) then
      call globalrandomize()
      tnextrand = tnextrand + randint
    end if

    ! Interpolation
    if(np > 0 .and. nplisted > 0) then
      particle => head
      do while( associated(particle) )
        if ( time - particle%tstart >= 0 ) then
          particle%partstep = particle%partstep + 1
          ! Interpolation of the velocity field
          particle%ures = ui3d(particle%x,particle%y,particle%z) * dxi
          particle%vres = vi3d(particle%x,particle%y,particle%z) * dyi
          particle%wres = wi3d(particle%x,particle%y,particle%z) * dzi_t(floor(particle%z))
          if (lpartsgs .and. nstep == 1) then
            call prep_sgs(particle)
            particle%usgs   = usgs(particle) * dxi 
            particle%vsgs   = vsgs(particle) * dyi
            particle%wsgs   = wsgs(particle) * dzi_t(floor(particle%z))
          end if
        end if
      particle => particle%next
      end do
    
      ! Integration
      particle => head
      do while( associated(particle) )
        if ( time - particle%tstart >= 0 ) then
          call rk3(particle)
          call checkbound(particle)
        end if
      particle => particle%next
      end do
    end if

    ! Communicate particles to other procs 
    call partcomm

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Temporary hack, sync particle dump with statistics
    ! see also step.f90
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    ! Statistics
    !if (nstep==3) then
    !  ! Particle dump
    !  if((time + (0.5*dt) >= tnextdump .or. time + dt >= timmax) .and. lpartdump) then
    !    !call particledump(time)
    !    call balanced_particledump(time)
    !    !call rawparticledump(time)
    !    tnextdump = tnextdump + frqpartdump
    !  end if
    !  !call checkdiv
    !end if
  

  end subroutine particles

  !
  !--------------------------------------------------------------------------
  ! subroutine calc_sgstke
  !> calculates sgs-tke, based on Stevens et al. (1999)
  !--------------------------------------------------------------------------
  !
  subroutine calc_sgstke
    use grid, only             : nzp, zt, dxi, dyi, nxp, nyp, zm, a_rp, a_tp, dzi_t, dzi_m, a_up, a_vp, a_wp, th00
    use defs, only             : pi, vonk, g
    use mpi_interface, only    : nxg, nyg
    implicit none
    
    real                       :: thvp,thvm, ri
    real                       :: S2,N2,cm,ch,cl,l,ch1,ch2,ceps1,ceps2
    integer                    :: i,j,k,ip,im,jp,jm,kp,km
    real, parameter            :: alpha = 1.6
    real, parameter            :: gamma = 1.34
    real, parameter            :: ric   = 0.23
    !real, parameter            :: ris   = 0.06

    ! Uncorrected parameters 
    cl    = (2./3.)**0.5
    cm    = (1./pi) * (2./(3.*alpha))**(3./2.)
    ch    = (4./(3.*gamma)) * (1./pi) * (2./(3.*alpha))**0.5 
    ceps  = pi * (2./(3.*alpha))**(3./2.)
    ch1   = cm
    ch2   = ch - ch1
    ceps1 = cm - cl**2 * (ric**(-1) - 1) 
    ceps2 = ceps - ceps1
    labda = ((1/dzi_t(1))/dxi/dyi)**(1./3.)

    do j=2,nyp-1
       jp = j+1
       jm = j-1
       do i=2,nxp-1
          ip = i+1
          im = i-1
          do k=2,nzp-1
            kp = k+1
            km = k-1            

            ! 1. Brunt-Vaisala^2
            thvp = (0.5*(a_tp(k,i,j) + a_tp(kp,i,j))) * (1. + 0.61 * (0.5*(a_rp(k,i,j) + a_rp(kp,i,j))))
            thvm = (0.5*(a_tp(k,i,j) + a_tp(km,i,j))) * (1. + 0.61 * (0.5*(a_rp(k,i,j) + a_rp(km,i,j))))
            N2   = max((g/th00) * (thvp - thvm) * dzi_t(k),1e-15)

            ! 2. Calculate S^2 = 2xSijSij, prevents using scratch arrays in sgsm.f90
            ! and bypasses unnecessary double interpolation
            ! {11,22,33}^2 = dudx^2 + dvdy^2 + dwdz^2
            S2 = ( &
              2. * ((a_up(k,i,j)   -  a_up(k,im,j))   * dxi       )**2    + &
              2. * ((a_vp(k,i,j)   -  a_vp(k,i,jm))   * dyi       )**2    + &
              2. * ((a_wp(k,i,j)   -  a_wp(km,i,j))   * dzi_t(k)  )**2) 

            ! {13 = 31}^2 = dvdz^2 + dwdx^2
            S2 = S2 + 0.25 * ( &
              ((a_wp(k,i,j)   -  a_wp(k,im,j))   * dxi       + &
               (a_up(kp,im,j) -  a_up(k,im,j))   * dzi_m(k)  )**2   + &
              ((a_wp(km,i,j)  -  a_wp(km,im,j))  * dxi       + &
               (a_up(k,im,j)  -  a_up(km,im,j))  * dzi_m(km) )**2   + &
              ((a_wp(km,ip,j) -  a_wp(km,i,j))   * dxi       + &
               (a_up(k,i,j)   -  a_up(km,i,j))   * dzi_m(km) )**2   + &
              ((a_wp(k,ip,j)  -  a_wp(k,i,j))    * dxi       + &
               (a_up(kp,i,j)  -  a_up(k,i,j))    * dzi_m(k)  )**2)

            ! {23 = 32}^2 = dvdz^2 + dwdy^2
            S2 = S2 + 0.25 * ( &
              ((a_wp(k,i,j)   -  a_wp(k,i,jm))   * dyi       + &
               (a_vp(kp,i,jm) -  a_vp(k,i,jm))   * dzi_m(k)  )**2   + &
              ((a_wp(km,i,j)  -  a_wp(km,i,jm))  * dyi       + &
               (a_vp(k,i,jm)  -  a_vp(km,i,jm))  * dzi_m(km) )**2   + &
              ((a_wp(km,i,jp) -  a_wp(km,i,j))   * dyi       + &
               (a_vp(k,i,j)   -  a_vp(km,i,j))   * dzi_m(km) )**2   + &
              ((a_wp(k,i,jp)  -  a_wp(k,i,j))    * dyi       + &
               (a_vp(kp,i,j)  -  a_vp(k,i,j))    * dzi_m(k)  )**2)

            ! {12 = 21}^2 = dudy^2 + dvdx^2
            S2 = S2 + 0.25 * ( &
              ((a_up(k,im,jp) -  a_up(k,im,j))   * dyi      + &
               (a_vp(k,i,j)   -  a_vp(k,im,j))   * dxi      )**2   + &
              ((a_up(k,im,j)  -  a_up(k,im,jm))  * dyi      + &
               (a_vp(k,i,jm)  -  a_vp(k,im,jm))  * dxi      )**2   + &
              ((a_up(k,i,j)   -  a_up(k,i,jm))   * dyi      + &
               (a_vp(k,ip,jm) -  a_vp(k,i,jm))   * dxi      )**2   + &
              ((a_up(k,i,jp)  -  a_up(k,i,j))    * dyi      + &
               (a_vp(k,ip,j)  -  a_vp(k,i,j))    * dxi      )**2)

            ! No stability correction
            ri            = N2 / S2
            l             = labda
            sgse(k,i,j)   = max(minsgse,(cm/ceps) * (l**2) * S2 * (1. - ((ch/cm) * ri)))

          end do
          sgse(1,i,j)     = -sgse(2,i,j)
          sgse(nzp,i,j)   = sgse(nzp-1,i,j)
       end do
    end do

  end subroutine calc_sgstke

  !
  !--------------------------------------------------------------------------
  ! subroutine calc_restke
  !> calculates resolved tke at cell center
  !--------------------------------------------------------------------------
  !
  subroutine calc_restke
    use grid, only              : nzp, a_up, a_vp, a_wp, nxp, nyp
    use mpi_interface, only     : ierror, mpi_double_precision, mpi_sum, mpi_comm_world, nxg, nyg
    implicit none
    
    integer                           :: i,j,k,ip,im,jp,jm,kp,km
    real, allocatable, dimension(:)   :: u_avl, v_avl, u_av, v_av
   
    allocate(u_avl(nzp),v_avl(nzp),u_av(nzp),v_av(nzp))

    do k=1,nzp
      u_avl(k)    = sum(a_up(k,3:nxp-2,3:nyp-2))
      v_avl(k)    = sum(a_vp(k,3:nxp-2,3:nyp-2))
    end do 

    call mpi_allreduce(u_avl,u_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(v_avl,v_av,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

    u_av        = u_av    / (nxg * nyg)  
    v_av        = v_av    / (nxg * nyg)     

    do j=2,nyp-1
       jp = j+1
       jm = j-1
       do i=2,nxp-1
          ip = i+1
          im = i-1
          do k=2,nzp-1
            kp = k+1
            km = k-1            
            rese(k,i,j) = ((0.5 * (a_up(k,i,j) + a_up(k,im,j)) - u_av(k))**2.  + &
                           (0.5 * (a_vp(k,i,j) + a_vp(k,i,jm)) - v_av(k))**2.  + &
                           (0.5 * (a_wp(k,i,j) + a_wp(km,i,j))          )**2.) / 2. 
          end do
          rese(1,i,j)   = rese(2,i,j)
          rese(nzp,i,j) = rese(nzp-1,i,j)
       end do
    end do

    deallocate(u_avl,v_avl,u_av,v_av)

  end subroutine calc_restke

  !
  !--------------------------------------------------------------------------
  ! subroutine fsubgrid_local : 
  !> Calculates fs (contribution sgs turbulence to total turbulence)
  !> at grid center
  !--------------------------------------------------------------------------
  !
  subroutine fsubgrid_local
    use grid, only    : nzp, nxp, nyp, dzi_m, nfpt
    implicit none

    integer :: i,j,k,im,ip,jm,jp,km,kp  
    real    :: dfsdz 
 
    do j=2,nyp-1
      jp = j+1
      jm = j-1
      do i=2,nxp-1
        ip = i+1
        im = i-1
        do k=2,nzp
          kp = k+1
          km = k-1
          if(k > nzp-nfpt) then       ! fix fs in ghost cells & sponge layer
            fs_local(k,i,j) = fs_local(k-1,i,j)
          else                        ! calculate fs
            fs_local(k,i,j)  = sgse(k,i,j) / (sgse(k,i,j) + rese(k,i,j))
          end if
        end do
        ! Extrapolate for lowest half level
        dfsdz = (fs_local(3,i,j) - fs_local(2,i,j)) * dzi_m(2)       
        fs_local(1,i,j) = fs_local(2,i,j) - (dfsdz / dzi_m(1))
      end do
    end do

  end subroutine fsubgrid_local

  !
  !--------------------------------------------------------------------------
  ! subroutine fsubgrid : 
  !> Calculates fs (contribution sgs turbulence to total turbulence)
  !--------------------------------------------------------------------------
  !
  subroutine fsubgrid
    use grid, only : a_up, a_vp, a_wp, nzp, nxp, nyp, dt, nstep, zt, nfpt
    use mpi_interface, only : ierror, mpi_double_precision, mpi_sum, mpi_comm_world, nxg, nyg
    implicit none
  
    integer    :: k
    real, allocatable, dimension(:)   :: &
       u_avl, v_avl, u2_avl, v2_avl, w2_avl, sgse_avl,     &
       u_av,  v_av,  u2_av,  v2_av,  w2_av,  sgse_av, e_res

    if(fixedfs) then
      fs = 1.
    else
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
        if(k > nzp-nfpt) then       ! fix fs in ghost cells & sponge layer
          fs(k) = fs(k-1)
        else
          fs(k)     = sgse_av(k) / (sgse_av(k) + e_res(k))
        end if
      end do

      e_res(1)    = -e_res(2)                   ! e_res -> 0 at surface
      fs(1)       = fs(2) - (fs(3) - fs(2))     ! fs    -> extrapolate to surface 
      !fs(1)       = 1. + (1.-fs(2))            ! fs    -> 1 at surface

      !Raw statistics dump
      !CAUTION:	NEEDS TO BE FIXED, STAT TIMEKEEPING NO LONGER IN PARTICLES.F90
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
    end if

  end subroutine fsubgrid

  !
  !--------------------------------------------------------------------------
  ! subroutine prep_sgs : prepares some variables for subgrid velocity
  !--------------------------------------------------------------------------
  ! 
  subroutine prep_sgs(particle)
    use grid, only : dxi, dyi, dt, a_scr7, zm, zt, dzi_t, dzi_m, a_tp, a_rp, nxp, nyp, nzp, th00, nfpt, distim
    use defs, only : g,pi
    implicit none

    real     :: deltaz
    integer  :: zbottom
    TYPE (particle_record), POINTER:: particle

    zbottom      = floor(particle%z + 0.5)
    deltaz       = ((zm(floor(particle%z)) + (particle%z - floor(particle%z)) / dzi_t(floor(particle%z))) - zt(zbottom)) * dzi_m(zbottom)
    if(lfsloc) then
      fsl        = i3d(particle%x,particle%y,particle%z,fs_local)
    else
      fsl        = (1-deltaz) * fs(zbottom) + deltaz * fs(zbottom+1)
    end if    

    sigma2l      = i3d(particle%x,particle%y,particle%z,sgse) * (2./3.)
    
    dsigma2dx    = (2./3.) * (i3d(particle%x+0.5,particle%y,particle%z,sgse) - &
                              i3d(particle%x-0.5,particle%y,particle%z,sgse)) * dxi
    dsigma2dy    = (2./3.) * (i3d(particle%x,particle%y+0.5,particle%z,sgse) - &
                              i3d(particle%x,particle%y-0.5,particle%z,sgse)) * dyi
    dsigma2dz    = (2./3.) * (i3d(particle%x,particle%y,particle%z+0.5,sgse) - &
                              i3d(particle%x,particle%y,particle%z-0.5,sgse)) * dzi_t(floor(particle%z))
    dsigma2dt    = (sigma2l - particle%sigma2_sgs) / dt
    particle%sigma2_sgs = sigma2l

    ! Damping subgrid velocities in sponge layer
    ! spngl = damping time scale at particle height.
    !if(particle%z > (nzp - nfpt)) then
    !  zparticle = zm(floor(particle%z)) + (particle%z-floor(particle%z)) / dzi_t(floor(particle%z))
    !  spngl = max(0.,(zm(nzp) - zparticle) / ((zm(nzp) - zm(nzp-nfpt)) * distim))
    !  spngl = max(0.,(1. / distim - spngl))
    !else
    !  spngl = 0.
    !end if

  end subroutine prep_sgs

  !
  !--------------------------------------------------------------------------
  ! subroutine usgs 
  !> Calculation subgrid (u) velocity particle based on subgrid-TKE
  !--------------------------------------------------------------------------
  !
  function usgs(particle)
    use grid, only : dxi, dt
    implicit none

    real :: t1, t2, t3, usgs
    TYPE (particle_record), POINTER:: particle

    t1        = (-0.75 * fsl * C0 * (particle%usgs / dxi) * (ceps/labda) * (1.5 * sigma2l)**0.5) * dt
    t2        = ( 0.5 * ((1. / sigma2l) * dsigma2dt * (particle%usgs / dxi) + dsigma2dx)) * dt
    t3        = ((fsl * C0 * (ceps/labda) * (1.5*sigma2l)**(1.5))**0.5 * xi(idum))
    !ts        = -spngl * (particle%usgs / dxi)

    ! 1. dissipation subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t1 + particle%usgs / dxi) .ne. sign(1.,particle%usgs / dxi)) t1 = - particle%usgs / dxi
    ! 2. decrease subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t2 + particle%usgs / dxi) .ne. sign(1.,particle%usgs / dxi)) t2 = - particle%usgs / dxi
    ! 3. Terms 1. and 2. combined may not exceed sugbrid velocity itself
    if(sign(1.,t1 + t2 + particle%usgs / dxi) .ne. sign(1.,particle%usgs / dxi)) then
      t1 = -0.5 * particle%usgs / dxi
      t2 = -0.5 * particle%usgs / dxi
    end if
     
    usgs = (particle%usgs / dxi) + (t1 + t2 + t3)  

  end function usgs

  !
  !--------------------------------------------------------------------------
  ! subroutine vsgs
  !> Calculation subgrid (v) velocity particle based on subgrid-TKE
  !--------------------------------------------------------------------------
  !
  function vsgs(particle)
    use grid, only : dyi, dt
    implicit none

    real :: t1, t2, t3, vsgs
    TYPE (particle_record), POINTER:: particle

    t1        = (-0.75 * fsl * C0 * (particle%vsgs / dyi) * (ceps/labda) * (1.5 * sigma2l)**0.5) * dt
    t2        = ( 0.5 * ((1. / sigma2l) * dsigma2dt * (particle%vsgs / dyi) + dsigma2dy)) * dt
    t3        = ((fsl * C0 * (ceps/labda) * (1.5*sigma2l)**(1.5))**0.5 * xi(idum))
    !ts        = -spngl * (particle%vsgs / dyi)

    ! 1. dissipation subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t1 + particle%vsgs / dyi) .ne. sign(1.,particle%vsgs / dyi)) t1 = - particle%vsgs / dyi
    ! 2. decrease subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t2 + particle%vsgs / dyi) .ne. sign(1.,particle%vsgs / dyi)) t2 = - particle%vsgs / dyi
    ! 3. Terms 1. and 2. combined may not exceed sugbrid velocity
    if(sign(1.,t1 + t2 + particle%vsgs / dyi) .ne. sign(1.,particle%vsgs / dyi)) then
      t1 = -0.5 * particle%vsgs / dyi
      t2 = -0.5 * particle%vsgs / dyi
    end if
 
    vsgs = (particle%vsgs / dyi) + (t1 + t2 + t3) 

  end function vsgs 

  !
  !--------------------------------------------------------------------------
  ! subroutine wsgs
  !> Calculation subgrid (w) velocity particle based on subgrid-TKE
  !--------------------------------------------------------------------------
  !
  function wsgs(particle)
    use grid, only : dzi_t, dt
    implicit none

    real :: t1, t2, t3, wsgs, dzi
    TYPE (particle_record), POINTER:: particle

    dzi        = dzi_t(floor(particle%z))
    t1        = (-0.75 * fsl * C0 * (particle%wsgs / dzi) * (ceps/labda) * (1.5 * sigma2l)**0.5) * dt
    t2        = ( 0.5 * ((1. / sigma2l) * dsigma2dt * (particle%wsgs / dzi) + dsigma2dz)) * dt
    t3        = ((fsl * C0 * (ceps/labda) * (1.5*sigma2l)**(1.5))**0.5 * xi(idum))
    !ts        = -spngl * (particle%wsgs / dzi)

    ! 1. dissipation subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t1 + particle%wsgs / dzi) .ne. sign(1.,particle%wsgs / dzi)) t1 = - particle%wsgs / dzi
    ! 2. decrease subgrid velocity cannot exceed subgrid velocity itself 
    if(sign(1.,t2 + particle%wsgs / dzi) .ne. sign(1.,particle%wsgs / dzi)) t2 = - particle%wsgs / dzi
    ! 3. Terms 1. and 2. combined may not exceed sugbrid velocity
    if(sign(1.,t1 + t2 + particle%wsgs / dzi) .ne. sign(1.,particle%wsgs / dzi)) then
      t1 = -0.5 * particle%wsgs / dzi
      t2 = -0.5 * particle%wsgs / dzi
    end if
 
    wsgs = (particle%wsgs / dzi) + (t1 + t2 + t3) 

  end function wsgs

  !
  !--------------------------------------------------------------------------
  ! Subroutine globalrandomize 
  !> Randomizes the X,Y,Z positions of all particles in 
  !> the lowest grid level every RK3 cycle. Called from: particles()
  !--------------------------------------------------------------------------
  !
  subroutine globalrandomize()
    use mpi_interface, only : nxg, nyg, nyprocs, nxprocs, mpi_integer, mpi_double_precision, mpi_sum, mpi_comm_world, ierror, wrxid, wryid, ranktable, mpi_status_size, myid
    use grid, only : deltax, deltay
    implicit none
  
    type (particle_record), pointer:: particle,ptr
    real               :: zmax = 1.                  ! Max height in grid coordinates
    !integer            :: status(mpi_status_size)
    real, allocatable, dimension(:)    :: buffsend,buffrecv
    integer, allocatable, dimension(:) :: recvcount,displacements
    integer            :: nglobal, nlocal=0, ii, i
    real               :: xsizelocal,ysizelocal,tempx,tempy
    real               :: randnr(3)

    ! Count number of local particles < zmax
    particle => head
    do while(associated(particle) )
      !if( particle%z <= (1. + zmax) .and. particle%wres < 0. ) nlocal = nlocal + 1
      if( particle%z <= (1. + zmax) ) nlocal = nlocal + 1
      particle => particle%next
    end do 

    call mpi_allreduce(nlocal,nglobal,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)

    if(nlocal > 0) then            
      ! Give them a random location and place in send buffer
      allocate(buffsend(nrpartvar * nlocal))
      ii = 0
      particle => head
      do while(associated(particle) )
        !if( particle%z <= (1. + zmax) .and. particle%wres < 0. ) then
        if( particle%z <= (1. + zmax) ) then
          call random_number(randnr)          ! Random seed has been called from init_particles...
          particle%x    = randnr(1) * float(nxg) 
          particle%y    = randnr(2) * float(nyg)
          particle%z    = 1. + (zmax * randnr(3))
          particle%ures = 0. 
          particle%vres = 0. 
          particle%wres = 0. 
          particle%usgs = 0. 
          particle%vsgs = 0. 
          particle%wsgs = 0. 

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

    ! Communicate number of local particles to each proc
    allocate(recvcount(nxprocs*nyprocs))
    call mpi_allgather(nlocal*nrpartvar,1,mpi_integer,recvcount,1,mpi_integer,mpi_comm_world,ierror)

    ! Create array to receive particles from other procs
    allocate(displacements(nxprocs*nyprocs))
    displacements(1) = 0
    do i = 2,nxprocs*nyprocs
      displacements(i) = displacements(i-1) + recvcount(i-1) 
    end do
    allocate(buffrecv(sum(recvcount)))

    ! Send all particles to all procs
    call mpi_allgatherv(buffsend,nlocal*nrpartvar,mpi_double_precision,buffrecv,recvcount,displacements,mpi_double_precision,mpi_comm_world,ierror)

    ! Loop through particles, check if on this proc
    xsizelocal = nxg / nxprocs
    ysizelocal = nyg / nyprocs
    ii = 0
    do i=1,nglobal   
      tempx = buffrecv(ii+2)
      tempy = buffrecv(ii+3)
      ! If on proc: add particle
      if(floor(tempx/xsizelocal) == wrxid) then
        if(floor(tempy/ysizelocal) == wryid) then
          call add_particle(particle)
          call partbuffer(particle,buffrecv(ii+1:ii+nrpartvar),ii,.false.)
          particle%x = particle%x - (floor(wrxid * xsizelocal)) + 3
          particle%y = particle%y - (floor(wryid * ysizelocal)) + 3
          particle%z = particle%z !+ 1
        end if
      end if
      ii = ii + nrpartvar
    end do

    ! Cleanup
    if(nlocal>0) deallocate(buffsend,buffrecv)
    deallocate(recvcount,displacements)
    nlocal  = 0
    nglobal = 0

    !ii = 0
    !write(hname,'(i4.4,a4)') myid,'glob'
    !open(998,file=hname,position='append',action='write')
    !do k=1,nglobal
    !  write(998,'(4F8.2)') buffrecv(ii+1),buffrecv(ii+2),buffrecv(ii+3),buffrecv(ii+4)
    !  ii = ii + nrpartvar
    !end do
    !close(998)
    !ii = 0
    !write(hname,'(i4.4,a3)') myid,'loc'
    !open(999,file=hname,position='append',action='write')
    !do k=1,nlocal
    !  write(999,'(4F8.2)') buffsend(ii+1),buffsend(ii+2),buffsend(ii+3),buffsend(ii+4)
    !  ii = ii + nrpartvar
    !end do
    !close(999)
    
  end subroutine globalrandomize

  !
  !--------------------------------------------------------------------------
  ! Function ui3d
  !> Performs a trilinear interpolation from the Eulerian grid to the 
  !> particle position.
  !--------------------------------------------------------------------------
  !
  function ui3d(x,y,z)
    use grid, only : a_up, dzi_m, dzi_t, zt, zm
    implicit none
    real, intent(in) :: x, y, z                               !< local x,y,z position in grid coordinates
    integer          :: xbottom, ybottom, zbottom
    real             :: ui3d, deltax, deltay, deltaz, sign

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
    ui3d          =  (1-deltaz) * (1-deltay) * (1-deltax) * sign * a_up(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) * sign * a_up(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) * sign * a_up(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) * sign * a_up(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *        a_up(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *        a_up(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *        a_up(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *        a_up(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function ui3d

  !
  !--------------------------------------------------------------------------
  ! Function vi3d
  !> Performs a trilinear interpolation from the Eulerian grid to the 
  !> particle position.
  !--------------------------------------------------------------------------
  !
  function vi3d(x,y,z)
    use grid, only : a_vp, dzi_m, dzi_t, zt, zm
    implicit none
    real, intent(in) :: x, y, z                               !< local x,y,z position in grid coordinates
    integer          :: xbottom, ybottom, zbottom
    real             :: vi3d, deltax, deltay, deltaz, sign

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
    vi3d          =  (1-deltaz) * (1-deltay) * (1-deltax) * sign * a_vp(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) * sign * a_vp(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) * sign * a_vp(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) * sign * a_vp(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *        a_vp(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *        a_vp(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *        a_vp(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *        a_vp(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function vi3d
  
  !
  !--------------------------------------------------------------------------
  ! Function wi3d
  !> Performs a trilinear interpolation from the Eulerian grid to the 
  !> particle position.
  !--------------------------------------------------------------------------
  !
  function wi3d(x,y,z)
    use grid, only : a_wp, dzi_m, dzi_t, zt, zm
    implicit none
    real, intent(in) :: x, y, z                               !< local x,y,z position in grid coordinates
    integer          :: xbottom, ybottom, zbottom
    real             :: wi3d, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z)
    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z - zbottom

    wi3d          =  (1-deltaz) * (1-deltay) * (1-deltax) *  a_wp(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) *  a_wp(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) *  a_wp(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) *  a_wp(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *  a_wp(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *  a_wp(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *  a_wp(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *  a_wp(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1
  
  end function wi3d

  !
  !--------------------------------------------------------------------------
  ! Function i3d
  !> trilinear interpolation from grid center to particle position
  !> Requires a 3D field as fourth argument
  !--------------------------------------------------------------------------
  !
  function i3d(x,y,z,input)
    implicit none
    real, intent(in)                    :: x, y, z                               !< local x,y,z position in grid coordinates
    real, dimension(:,:,:), intent(in)  :: input                                 !< scalar field used to interpolate from
 
    integer  :: xbottom, ybottom, zbottom
    real     :: i3d, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z + 0.5)
    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z + 0.5 - zbottom
    
    i3d           =  (1-deltaz) * (1-deltay) * (1-deltax) *  input(zbottom    , xbottom    , ybottom    ) + &    !
    &                (1-deltaz) * (1-deltay) * (  deltax) *  input(zbottom    , xbottom + 1, ybottom    ) + &    ! x+1
    &                (1-deltaz) * (  deltay) * (1-deltax) *  input(zbottom    , xbottom    , ybottom + 1) + &    ! y+1
    &                (1-deltaz) * (  deltay) * (  deltax) *  input(zbottom    , xbottom + 1, ybottom + 1) + &    ! x+1,y+1
    &                (  deltaz) * (1-deltay) * (1-deltax) *  input(zbottom + 1, xbottom    , ybottom    ) + &    ! z+1
    &                (  deltaz) * (1-deltay) * (  deltax) *  input(zbottom + 1, xbottom + 1, ybottom    ) + &    ! x+1, z+1
    &                (  deltaz) * (  deltay) * (1-deltax) *  input(zbottom + 1, xbottom    , ybottom + 1) + &    ! y+1, z+1
    &                (  deltaz) * (  deltay) * (  deltax) *  input(zbottom + 1, xbottom + 1, ybottom + 1)        ! x+1,y+1,z+1

  end function i3d

  !
  !--------------------------------------------------------------------------
  ! Function i1d
  !> linear interpolation of grid center to particle position
  !> Requires a 1D profile as second argument
  !--------------------------------------------------------------------------
  !
  function i1d(z,input)
    implicit none
    real, intent(in)                    :: z                                     !< local z position in grid coordinates
    real, dimension(:), intent(in)      :: input                                 !< profile used to interpolate from
 
    integer  :: zbottom
    real     :: i1d, deltaz

    zbottom = floor(z + 0.5)
    deltaz = z + 0.5 - zbottom

    i1d    =  (1-deltaz) * input(zbottom) + deltaz * input(zbottom+1)

  end function i1d

  !
  !--------------------------------------------------------------------------
  ! Subroutine rk3 
  !> Third-order Runge-Kutta scheme for spatial integration of the particles.
  !--------------------------------------------------------------------------
  !
  subroutine rk3(particle)
    use grid, only : rkalpha, rkbeta, nstep, dt, dxi, dyi, dzi_t
    implicit none
    TYPE (particle_record), POINTER:: particle
      
    particle%zprev = particle%z

    particle%x     = particle%x + rkalpha(nstep) * (particle%ures + particle%usgs) * dt + rkbeta(nstep) * (particle%ures_prev + particle%usgs_prev) * dt
    particle%y     = particle%y + rkalpha(nstep) * (particle%vres + particle%vsgs) * dt + rkbeta(nstep) * (particle%vres_prev + particle%vsgs_prev) * dt
    particle%z     = particle%z + rkalpha(nstep) * (particle%wres + particle%wsgs) * dt + rkbeta(nstep) * (particle%wres_prev + particle%wsgs_prev) * dt

    particle%ures_prev = particle%ures
    particle%vres_prev = particle%vres
    particle%wres_prev = particle%wres

    particle%usgs_prev = particle%usgs
    particle%vsgs_prev = particle%vsgs
    particle%wsgs_prev = particle%wsgs
   
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
  ! Subroutine partcomm 
  !> Handles the cyclic boundary conditions (through MPI) and sends
  !> particles from processor to processor
  !--------------------------------------------------------------------------
  !
  subroutine partcomm
    use mpi_interface, only : wrxid, wryid, ranktable, nxg, nyg, xcomm, ycomm, ierror, mpi_status_size, mpi_integer, mpi_double_precision, mpi_comm_world, nyprocs, nxprocs,myid
    implicit none

    type (particle_record), pointer:: particle,ptr
    real, allocatable, dimension(:) :: buffsend, buffrecv
    integer :: status(mpi_status_size)
    integer :: request
    integer :: ii, n
    ! Number of particles to ('to') and from ('fr') N,E,S,W
    integer :: nton,ntos,ntoe,ntow
    integer :: nfrn,nfrs,nfre,nfrw
    integer :: nyloc, nxloc

    nton  = 0
    ntos  = 0
    ntoe  = 0
    ntow  = 0
    nfrn  = 0
    nfrs  = 0
    nfre  = 0
    nfrw  = 0

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
   
    ! ---------------------------
    if( nton > 0 ) then
      allocate(buffsend(nrpartvar * nton))

      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%y >= nyloc + 3 ) then
          particle%y      = particle%y      - nyloc
          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do
      call mpi_isend(buffsend,nrpartvar*nton,mpi_double_precision,ranktable(wrxid,wryid+1),6,mpi_comm_world,request,ierror)
    end if
   
    if(nfrs > 0) then
      allocate(buffrecv(nrpartvar * nfrs))
      call mpi_recv(buffrecv,nrpartvar*nfrs,mpi_double_precision,ranktable(wrxid,wryid-1),6,mpi_comm_world,status,ierror)
      ii = 0
      do n = 1,nfrs
        call add_particle(particle)
        call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
        ii=ii+nrpartvar
      end do
    end if

    ! Wait for comm to finish
    if(nton>0) then 
      call mpi_wait(request,status,ierror)
      deallocate(buffsend)
    end if
    if(nfrs>0) deallocate(buffrecv)

    ! ---------------------------
    if( ntos > 0 ) then
      allocate(buffsend(nrpartvar * ntos))
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%y < 3 ) then
          particle%y      = particle%y      + nyloc
          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do
      call mpi_isend(buffsend,nrpartvar*ntos,mpi_double_precision,ranktable(wrxid,wryid-1),7,mpi_comm_world,request,ierror)
    end if

    if(nfrn > 0) then
      allocate(buffrecv(nrpartvar * nfrn))
      call mpi_recv(buffrecv,nrpartvar*nfrn,mpi_double_precision,ranktable(wrxid,wryid+1),7,mpi_comm_world,status,ierror)
      ii = 0
      do n = 1,nfrn
        particle => head
        call add_particle(particle)
        call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
        ii=ii+nrpartvar
      end do
    end if

    ! Wait for comm to finish
    if(ntos>0) then 
      call mpi_wait(request,status,ierror)
      deallocate(buffsend)
    end if
    if(nfrn>0) deallocate(buffrecv)

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

    ! ---------------------------
    if( ntoe > 0 ) then
      allocate(buffsend(nrpartvar * ntoe))
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%x >= nxloc + 3 ) then
          particle%x      = particle%x      - nxloc
          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do

      call mpi_isend(buffsend,nrpartvar*ntoe,mpi_double_precision,ranktable(wrxid+1,wryid),10,mpi_comm_world,request,ierror)
    end if

    if(nfrw > 0) then
      allocate(buffrecv(nrpartvar * nfrw))
      call mpi_recv(buffrecv,nrpartvar*nfrw,mpi_double_precision,ranktable(wrxid-1,wryid),10,mpi_comm_world,status,ierror)
      ii = 0
      do n = 1,nfrw
        call add_particle(particle)
        call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
        ii=ii+nrpartvar
      end do
    end if

    ! Wait for comm to finish
    if(ntoe>0) then 
      call mpi_wait(request,status,ierror)
      deallocate(buffsend)
    end if
    if(nfrw>0) deallocate(buffrecv)

    ! ---------------------------
    if( ntow > 0 ) then
      allocate(buffsend(nrpartvar * ntow))
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%x < 3 ) then
          particle%x      = particle%x      + nxloc
          call partbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call delete_particle(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do
      call mpi_isend(buffsend,nrpartvar*ntow,mpi_double_precision,ranktable(wrxid-1,wryid),11,mpi_comm_world,request,ierror)
    end if

    if( nfre > 0) then
      allocate(buffrecv(nrpartvar * nfre))
      call mpi_recv(buffrecv,nrpartvar*nfre,mpi_double_precision,ranktable(wrxid+1,wryid),11,mpi_comm_world,status,ierror)
      ii = 0
      do n = 1,nfre
        particle => head
        call add_particle(particle)
        call partbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
        ii=ii+nrpartvar
      end do
    end if

    ! Wait for comm to finish
    if(ntow>0) then 
      call mpi_wait(request,status,ierror)
      deallocate(buffsend)
    end if
    if(nfre>0) deallocate(buffrecv)

  end subroutine partcomm

  !
  !--------------------------------------------------------------------------
  ! Subroutine partbuffer 
  !> Packs/receives particle records to/from an array, sendable over MPI 
  !--------------------------------------------------------------------------
  !
  subroutine partbuffer(particle, buffer, n, send)
    implicit none

    logical,intent(in)                :: send                               !< 
    integer,intent(in)                :: n
    real,dimension(n+1:n+nrpartvar)   :: buffer
    TYPE (particle_record), POINTER:: particle

    if (send) then
      buffer(n+ipunique)        = particle%unique
      buffer(n+ipx)             = particle%x
      buffer(n+ipy)             = particle%y
      buffer(n+ipz)             = particle%z
      buffer(n+ipzprev)         = particle%zprev
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
      buffer(n+ipsigma2_sgs)    = particle%sigma2_sgs
      buffer(n+ipxstart)        = particle%xstart
      buffer(n+ipystart)        = particle%ystart
      buffer(n+ipzstart)        = particle%zstart
      buffer(n+iptsart)         = particle%tstart
      buffer(n+ipartstep)       = particle%partstep
    else
      particle%unique           = buffer(n+ipunique)
      particle%x                = buffer(n+ipx)
      particle%y                = buffer(n+ipy)
      particle%z                = buffer(n+ipz)
      particle%zprev            = buffer(n+ipzprev)
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
      particle%sigma2_sgs       = buffer(n+ipsigma2_sgs)
      particle%xstart           = buffer(n+ipxstart)
      particle%ystart           = buffer(n+ipystart)
      particle%zstart           = buffer(n+ipzstart)
      particle%tstart           = buffer(n+iptsart)
      particle%partstep         = int(buffer(n+ipartstep))
    end if

  end subroutine partbuffer

  !
  !--------------------------------------------------------------------------
  ! Subroutine thermo 
  !> Calculates thermodynamic variables at particle position (thl, thv,
  !> qt, qs)  
  !--------------------------------------------------------------------------
  !
  subroutine thermo(px,py,pz,thl,thv,rt,rl)
    use thrm,         only : rslf
    use grid,         only : a_pexnr, a_rp, a_theta, a_tp, pi0, pi1,th00
    !use grid,         only : tname,nzp,dxi,dyi,dzi_t,nxp,nyp,umean,vmean, a_tp, a_rp, press, th00, a_pexnr, a_theta,pi0,pi1
    
    use defs,         only : p00,cp,R,Rm,tmelt,alvl,cpr,ep2,ep
    use thrm,         only : rslf
    implicit none
   
    real, intent(in)  :: px,py,pz
    real, intent(out) :: thl,thv,rt,rl 
    real, parameter   :: epsln = 1.e-4
    real              :: exner,ploc,tlloc,rsloc,dtx,tx,txi,tx1
    integer           :: iterate

    ! scalar interpolations and calculations
    exner   = (i1d(pz,pi0)+i1d(pz,pi1)+i3d(px,py,pz,a_pexnr)) / cp
    ploc    = p00 * exner**cpr               ! Pressure
    thl     = i3d(px,py,pz,a_tp) + th00      ! Liquid water potential T 
    tlloc   = thl * exner                    ! Liquid water T 
    rsloc   = rslf(ploc,tlloc)               ! Saturation vapor mixing ratio
    rt      = i3d(px,py,pz,a_rp)             ! Total water mixing ratio
    rl      = max(rt-rsloc,0.)               ! Liquid water mixing ratio

    if(rl > 0.) then
      dtx          = 2. * epsln
      iterate      = 1
      tx           = tlloc
      do while(dtx > epsln .and. iterate < 20)
        txi        = alvl / (cp * tx)
        tx1        = tx - (tx - tlloc * (1. + txi  * rl)) / &
                       (1. + txi * tlloc * (rl / tx + (1. + rsloc * ep) * rsloc * alvl / (Rm * tx * tx)))
        dtx        = abs(tx1 - tx)
        tx         = tx1
        iterate    = iterate + 1
        rsloc      = rslf(ploc,tx)
        rl         = max(rt-rsloc,0.)
      end do
    end if

    thv = i3d(px,py,pz,a_theta) * (1.+ep2*(rt-rl))  

  end subroutine thermo

  !
  !--------------------------------------------------------------------------
  ! Subroutine particlestat 
  !> Performs the sampling and saving of binned and slab averaged particle
  !> statistics. Output written to *.particlestat.nc  
  !--------------------------------------------------------------------------
  !
  subroutine particlestat(dowrite,time)
    use mpi_interface, only : mpi_comm_world, myid, mpi_double_precision, mpi_sum, ierror, nxprocs, nyprocs, nxg, nyg
    use modnetcdf,     only : writevar_nc, fillvalue_double
    use grid,          only : tname,dxi,dyi,dzi_t,nzp,umean,vmean,nzp
    implicit none

    logical, intent(in)     :: dowrite
    real, intent(in)        :: time
    integer                 :: k
    real                    :: thv,thl,rt,rl           ! From subroutine thermo
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
        u2profl(k)      = u2profl(k)    + (particle%ures / dxi)**2.
        v2profl(k)      = v2profl(k)    + (particle%vres / dyi)**2.
        w2profl(k)      = w2profl(k)    + (particle%wres / dzi_t(floor(particle%z)))**2.

        call thermo(particle%x,particle%y,particle%z,thl,thv,rt,rl)

        ! scalar profiles
        tprofl(k)       = tprofl(k)     + thl
        tvprofl(k)      = tvprofl(k)    + thv
        rtprofl(k)      = rtprofl(k)    + rt
        rlprofl(k)      = rlprofl(k)    + rl
        if(rl > 0.)     ccprofl(k)      = ccprofl(k)     + 1
        if(lpartsgs) then
          sigma2profl(k)        = sigma2profl(k) + particle%sigma2_sgs
          if(lfsloc) fsprofl(k) = fsprofl(k) + i3d(particle%x,particle%y,particle%z,fs_local)
        end if
        particle => particle%next
      end do 

      nstatsamp = nstatsamp + 1
    end if

    ! Write to NetCDF
    if(dowrite) then
      call mpi_allreduce(npartprofl,npartprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(uprofl,uprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(vprofl,vprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(wprofl,wprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(u2profl,u2prof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(v2profl,v2prof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(w2profl,w2prof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      ! scalars
      call mpi_allreduce(tprofl,tprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(tvprofl,tvprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(rtprofl,rtprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(rlprofl,rlprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(ccprofl,ccprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      if(lpartsgs) then
        call mpi_allreduce(sigma2profl,sigma2prof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        if(lfsloc) call mpi_allreduce(fsprofl,fsprof,nzp,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      end if 

      ! Divide summed values by ntime and nparticle samples and
      ! correct for Galilean transformation 
      do k = 1,nzp
        if(npartprof(k) > 0) then 
          npartprof(k) = npartprof(k) / (nstatsamp)
          uprof(k)     = uprof(k)     / (nstatsamp * npartprof(k)) + umean
          vprof(k)     = vprof(k)     / (nstatsamp * npartprof(k)) + vmean
          wprof(k)     = wprof(k)     / (nstatsamp * npartprof(k))
          u2prof(k)    = u2prof(k)    / (nstatsamp * npartprof(k)) - (uprof(k)-umean)**2.
          v2prof(k)    = v2prof(k)    / (nstatsamp * npartprof(k)) - (vprof(k)-vmean)**2.
          w2prof(k)    = w2prof(k)    / (nstatsamp * npartprof(k)) - wprof(k)**2. 
          tkeprof(k)   = 0.5 * (u2prof(k) + v2prof(k) + w2prof(k))
          ! scalars
          tprof(k)     = tprof(k)     / (nstatsamp * npartprof(k))
          tvprof(k)    = tvprof(k)    / (nstatsamp * npartprof(k))
          rtprof(k)    = rtprof(k)    / (nstatsamp * npartprof(k))
          rlprof(k)    = rlprof(k)    / (nstatsamp * npartprof(k))
          ccprof(k)    = ccprof(k)    / (nstatsamp * npartprof(k))
          if(lpartsgs) then
            sigma2prof(k)  = 1.5 * sigma2prof(k)  / (nstatsamp * npartprof(k))
            if(lfsloc) fsprof(k) = fsprof(k)   / (nstatsamp * npartprof(k))
          end if
        end if
      end do      

      if(myid==0) print*,'particles 1-2-3-4:',npartprof(2),npartprof(3),npartprof(4),npartprof(5)

      if(myid == 0) then
        call writevar_nc(ncpartstatid,tname,time,ncpartstatrec)
        call writevar_nc(ncpartstatid,'np',npartprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'u',uprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'v',vprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'w',wprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'u_2',u2prof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'v_2',v2prof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'w_2',w2prof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'tke',tkeprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'t',tprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'tv',tvprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'rt',rtprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'rl',rlprof,ncpartstatrec)
        call writevar_nc(ncpartstatid,'cc',ccprof,ncpartstatrec)
        if(lpartsgs) then
          if(lfsloc) then
            call writevar_nc(ncpartstatid,'fs',fsprof,ncpartstatrec)
          else
            call writevar_nc(ncpartstatid,'fs',fs,ncpartstatrec)
          end if
          call writevar_nc(ncpartstatid,'sgstke',sigma2prof,ncpartstatrec)
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
      u2prof     = 0
      u2profl    = 0
      v2prof     = 0
      v2profl    = 0
      w2prof     = 0
      w2profl    = 0
      tkeprof    = 0
      tkeprofl   = 0
      tprof      = 0
      tprofl     = 0
      tvprof     = 0
      tvprofl    = 0
      rtprof     = 0
      rtprofl    = 0
      rlprof     = 0
      rlprofl    = 0
      ccprof     = 0
      ccprofl    = 0
      if(lpartsgs) then
        fsprof      = 0
        fsprofl     = 0
        sigma2prof  = 0.
        sigma2profl = 0.
      end if
      nstatsamp  = 0
    end if

  end subroutine particlestat

  !
  !--------------------------------------------------------------------------
  ! subroutine rawparticledump : Quick'n'dirty binary particle dump per core 
  !   to test performance versus NetCDF
  !--------------------------------------------------------------------------
  !
  subroutine rawparticledump(time)
    use grid, only   : tname, deltax, deltay, dzi_t, zm
    use mpi_interface, only : nxprocs, nyprocs, wrxid, wryid, nxg, nyg

    implicit none

    real, intent(in)                     :: time
    real                                 :: px,py,pz,pu,pv,pw
    real                                 :: thl,thv,rt,rl
    integer                              :: ploc
    character (len=80)                   :: hname
    type (particle_record), pointer:: particle

    write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
    hname = 'rawparticles.'//trim(hname)

    open(123,file=hname,form='unformatted',position = 'append',action='write')
    write(123) time

    ! Loop through particles
    particle => head
    do while( associated(particle) )
      call thermo(particle%x,particle%y,particle%z,thl,thv,rt,rl)
      ploc = particle%unique
      px = (wrxid * (nxg / nxprocs) + particle%x - 3) * deltax
      py = (wryid * (nyg / nyprocs) + particle%y - 3) * deltay
      pz = zm(floor(particle%z)) + (particle%z-floor(particle%z)) / dzi_t(floor(particle%z))
      pu = particle%ures * deltax
      pv = particle%vres * deltay
      pw = particle%wres / dzi_t(floor(particle%z))

      write(123) ploc,px,py,pz,pu,pv,pw,thl,thv,rt,rl
      particle => particle%next
    end do

    close(123)

  end subroutine rawparticledump

  !
  !--------------------------------------------------------------------------
  ! subroutine balanced_particledump : More balanced communication, sends particles
  !   back to home processor before writing. Should decrease load on proc #0
  !--------------------------------------------------------------------------
  !
  subroutine balanced_particledump(time)
    use mpi_interface, only : mpi_comm_world, myid, mpi_integer, mpi_double_precision, ierror, nxprocs, nyprocs, mpi_status_size, wrxid, wryid, nxg, nyg, mpi_sum
    use grid,          only : tname, deltax, deltay, dzi_t, zm, umean, vmean
    use modnetcdf,     only : writevar_nc, fillvalue_double
    implicit none

    real, intent(in)                     :: time
    type (particle_record),       pointer:: particle
    integer                              :: status(mpi_status_size)
    integer                              :: nlocal,nprocs,p,nvar,start,end,nsr,isr,nvl,loc,bloc,ii
    integer, allocatable, dimension(:)   :: tosend,toreceive,base,sendbase,receivebase
    real, allocatable, dimension(:)      :: sendbuff,recvbuff
    real, allocatable, dimension(:,:)    :: sb_sorted
    integer, allocatable, dimension(:,:) :: status_array
    integer, allocatable, dimension(:)   :: req
    real                                 :: thl,thv,rt,rl

    if(time > tnextdump) then

      nvar = 4                             ! id,x,y,z
      if(lpartdumpui)  nvar = nvar + 3     ! u,v,w
      if(lpartsgs)     nvar = nvar + 3     ! us,vs,ws
      if(lpartdumpth)  nvar = nvar + 2     ! thl,tvh
      if(lpartdumpmr)  nvar = nvar + 2     ! rt,rl

      nprocs = nxprocs * nyprocs
      allocate(tosend(0:nprocs-1),toreceive(0:nprocs-1),base(0:nprocs-1),sendbase(0:nprocs-1),receivebase(0:nprocs-1))

      ! Find average number of particles per proc
      nlocal = floor(real(np) / nprocs)

      ! Determine how many particles to send to which proc
      tosend = 0
      particle => head
      do while( associated(particle) )
        p = floor((particle%unique-1) / nlocal)       ! Which proc to send to
        if(p .gt. nprocs-1) p = nprocs - 1            ! Last proc gets remaining particles
        tosend(p) = tosend(p) + 1
        particle => particle%next
      end do

      ! 1. Communicate nparticles to send/receive to/from each proc
      ! 2. Find start position for each proc in send/recv buff
      do p=0,nprocs-1
        call mpi_sendrecv(tosend(p),   1,mpi_integer,p,1, &
                          toreceive(p),1,mpi_integer,p,1, &
                          mpi_comm_world,status,ierror)
        if(p .eq. 0) then
          sendbase(p)    = 1
          receivebase(p) = 1
        else
          sendbase(p)    = sendbase(p-1)    + (tosend(p-1)    * nvar)
          receivebase(p) = receivebase(p-1) + (toreceive(p-1) * nvar)
        end if
      end do

      base = sendbase  ! will be changed during filling send buffer

      ! Allocate send/receive buffers
      allocate(sendbuff(sum(tosend)*nvar),recvbuff(sum(toreceive)*nvar))

      ! Fill send buffer
      particle => head
      do while( associated(particle) )
        p = floor((particle%unique-1) / nlocal)       ! Which proc to send to
        if(p .gt. nprocs-1) p = nprocs-1              ! Last proc gets remaining particles

        if(lpartdumpth .or. lpartdumpmr) call thermo(particle%x,particle%y,particle%z,thl,thv,rt,rl)

        sendbuff(base(p))           =  particle%unique
        sendbuff(base(p)+1)         = (wrxid * (nxg / nxprocs) + particle%x - 3) * deltax
        sendbuff(base(p)+2)         = (wryid * (nyg / nyprocs) + particle%y - 3) * deltay
        sendbuff(base(p)+3)         = zm(floor(particle%z)) + (particle%z-floor(particle%z)) / dzi_t(floor(particle%z))

        nvl = 3
        if(lpartdumpui) then
          sendbuff(base(p)+nvl+1)   = particle%ures * deltax
          sendbuff(base(p)+nvl+2)   = particle%vres * deltay
          sendbuff(base(p)+nvl+3)   = particle%wres / dzi_t(floor(particle%z))
          nvl = nvl + 3
          if(lpartsgs) then
            sendbuff(base(p)+nvl+1) = particle%usgs * deltax
            sendbuff(base(p)+nvl+2) = particle%vsgs * deltay
            sendbuff(base(p)+nvl+3) = particle%wsgs / dzi_t(floor(particle%z))
            nvl = nvl + 3
          end if
        end if
        if(lpartdumpth) then
          sendbuff(base(p)+nvl+1)   = thl
          sendbuff(base(p)+nvl+2)   = thv
          nvl = nvl + 2
        end if
        if(lpartdumpmr) then
          sendbuff(base(p)+nvl+1)   = rt
          sendbuff(base(p)+nvl+2)   = rl
        end if

        base(p)             = base(p) + nvar

        particle => particle%next
      end do

      ! Non=blocking send and receive from/to each other proc
      ! Find total #send/recv's for non-blocking send/recv request checking 
      nsr = 0
      do p = 0, nprocs-1
        if(tosend(p)    .gt. 0) nsr = nsr + 1
        if(toreceive(p) .gt. 0) nsr = nsr + 1
      end do
      ! Allocate status & request arrays
      allocate(status_array(mpi_status_size,nsr),req(nsr))
 
      ! Do send/receives 
      isr = 1 
      do p=0,nprocs-1
        if(tosend(p) .gt. 0) then
          start = sendbase(p)
          if(p .eq. nprocs-1) then
            end = size(sendbuff)
          else
            end = sendbase(p+1)-1
          end if
          call mpi_isend(sendbuff(start:end),tosend(p)*nvar,mpi_double_precision,p,(myid+1)*(p+nprocs),mpi_comm_world,req(isr),ierror)
          isr = isr + 1
        end if
        if(toreceive(p) .gt. 0) then
          start = receivebase(p)
          if(p .eq. nprocs-1) then
            end = size(recvbuff)
          else 
            end = receivebase(p+1)-1
          end if       
          call mpi_irecv(recvbuff(start:end),toreceive(p)*nvar,mpi_double_precision,p,(p+1)*(myid+nprocs),mpi_comm_world,req(isr),ierror)
          isr = isr + 1
        end if
      end do

      ! Wait for all communication to finish
      call mpi_waitall(nsr,req,status_array,ierror)

      ! Sort particles
      allocate(sb_sorted(size(recvbuff)/nvar,nvar-1))
      bloc = myid * nlocal + 1
      ii = 1
      do p = 1, size(recvbuff)/nvar

        loc = recvbuff(ii)-bloc+1
        sb_sorted(loc,1) = recvbuff(ii+1)
        sb_sorted(loc,2) = recvbuff(ii+2)
        sb_sorted(loc,3) = recvbuff(ii+3)
        nvl = 3 
        if(lpartdumpui) then
          sb_sorted(loc,nvl+1) = recvbuff(ii+nvl+1) + umean
          sb_sorted(loc,nvl+2) = recvbuff(ii+nvl+2) + vmean
          sb_sorted(loc,nvl+3) = recvbuff(ii+nvl+3)
          nvl = nvl + 3
          if(lpartsgs) then
            sb_sorted(loc,nvl+1) = recvbuff(ii+nvl+1)
            sb_sorted(loc,nvl+2) = recvbuff(ii+nvl+2)
            sb_sorted(loc,nvl+3) = recvbuff(ii+nvl+3)
            nvl = nvl + 3
          end if
        end if
        if(lpartdumpth) then
          sb_sorted(loc,nvl+1) = recvbuff(ii+nvl+1)
          sb_sorted(loc,nvl+2) = recvbuff(ii+nvl+2)
          nvl = nvl + 2
        end if
        if(lpartdumpmr) then
          sb_sorted(loc,nvl+1) = recvbuff(ii+nvl+1)
          sb_sorted(loc,nvl+2) = recvbuff(ii+nvl+2)
        end if

        ii = ii + nvar
      end do

      ! Write to NetCDF
      call writevar_nc(ncpartid,tname,time,ncpartrec)
      call writevar_nc(ncpartid,'x',sb_sorted(:,1),ncpartrec)
      call writevar_nc(ncpartid,'y',sb_sorted(:,2),ncpartrec)
      call writevar_nc(ncpartid,'z',sb_sorted(:,3),ncpartrec)
      nvl = 3
      if(lpartdumpui) then
        call writevar_nc(ncpartid,'u',sb_sorted(:,nvl+1),ncpartrec)
        call writevar_nc(ncpartid,'v',sb_sorted(:,nvl+2),ncpartrec)
        call writevar_nc(ncpartid,'w',sb_sorted(:,nvl+3),ncpartrec)
        nvl = nvl + 3
        if(lpartsgs) then
          call writevar_nc(ncpartid,'us',sb_sorted(:,nvl+1),ncpartrec)
          call writevar_nc(ncpartid,'vs',sb_sorted(:,nvl+2),ncpartrec)
          call writevar_nc(ncpartid,'ws',sb_sorted(:,nvl+3),ncpartrec)
          nvl = nvl + 3
        end if
      end if
      if(lpartdumpth) then
        call writevar_nc(ncpartid,'t', sb_sorted(:,nvl+1),ncpartrec)
        call writevar_nc(ncpartid,'tv',sb_sorted(:,nvl+2),ncpartrec)
        nvl = nvl + 2
      end if
      if(lpartdumpmr) then
        call writevar_nc(ncpartid,'rt',sb_sorted(:,nvl+1),ncpartrec)
        call writevar_nc(ncpartid,'rl',sb_sorted(:,nvl+2),ncpartrec)
      end if

      ! Cleanup!
      deallocate(tosend,toreceive,base,sendbase,receivebase)
      deallocate(sendbuff,recvbuff)
      deallocate(sb_sorted)
      deallocate(status_array,req)

    end if
 
  end subroutine balanced_particledump


  !
  !--------------------------------------------------------------------------
  ! Quick and dirty local divergence check
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
 
  if(myid==0) print*,'   divergence; max=',divmax, ', total=',divtot 

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
    elseif (particle%z < 1.01) then
      particle%z = abs(particle%z-1.01)+1.01
      particle%wres =  abs(particle%wres)
    end if

  end subroutine checkbound

  !
  !--------------------------------------------------------------------------
  ! function xi & random : creates component of white Gaussian noise (random
  !   number from Gaussian distr.) using the Box-Müller algorithm   
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
  subroutine init_particles(hot,hfilin)
    use mpi_interface, only : wrxid, wryid, nxg, nyg, myid, nxprocs, nyprocs, appl_abort, ierror,mpi_double_precision,mpi_comm_world,mpi_min
    use grid, only : zm, deltax, deltay, zt,dzi_t, nzp, nxp, nyp
    use grid, only : a_up, a_vp, a_wp

    logical, intent(in) :: hot
    character (len=80), intent(in), optional :: hfilin
    integer  :: k, n, kmax, io
    logical  :: exans
    real     :: tstart, xstart, ystart, zstart, ysizelocal, xsizelocal, firststartl, firststart
    real     :: pu,pts,px,py,pz,pzp,pxs,pys,pzs,pur,pvr,pwr,purp,pvrp,pwrp
    real     :: pus,pvs,pws,pusp,pvsp,pwsp,psg2
    integer  :: pstp,idot
    type (particle_record), pointer:: particle
    character (len=80) :: hname,prefix,suffix

    xsizelocal = (nxg / nxprocs) * deltax
    ysizelocal = (nyg / nyprocs) * deltay
    kmax = size(zm)

    firststartl = 1e9
    call init_random_seed()

    ! clear pointers to head and tail
    nullify(head)
    nullify(tail)
    nplisted = 0

    if(hot) then
    ! ------------------------------------------------------
    ! Warm start -> load restart file

      write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      idot = scan(hfilin,'.',.false.)
      prefix = hfilin(:idot-1)
      suffix = hfilin(idot+1:)
      hname = trim(hname)//'.'//trim(prefix)//'.particles.'//trim(suffix)
      inquire(file=trim(hname),exist=exans)
      if (.not.exans) then
         print *,'ABORTING: History file', trim(hname),' not found'
         call appl_abort(0)
      end if
      open (666,file=hname,status='old',form='unformatted')
      read (666,iostat=io) np,tnextdump
      do
        read (666,iostat=io) pu,pts,pstp,px,pxs,pur,purp,py,pys,pvr,pvrp,pz,pzs,pzp,pwr,pwrp,pus,pvs,pws,pusp,pvsp,pwsp,psg2
        if(io .ne. 0) exit
        call add_particle(particle)
        particle%unique         = pu
        particle%x              = px
        particle%y              = py
        particle%z              = pz
        particle%zprev          = pzp
        particle%xstart         = pxs
        particle%ystart         = pys
        particle%zstart         = pzs
        particle%tstart         = pts
        particle%ures           = pur
        particle%vres           = pvr
        particle%wres           = pwr
        particle%ures_prev      = purp
        particle%vres_prev      = pvrp
        particle%wres_prev      = pwrp
        particle%usgs           = pus
        particle%vsgs           = pvs
        particle%wsgs           = pws
        particle%usgs_prev      = pusp
        particle%vsgs_prev      = pvsp
        particle%wsgs_prev      = pwsp
        particle%sigma2_sgs     = psg2
        particle%partstep       = pstp
        if(pts < firststartl) firststartl = pts
      end do
      close(666)

      call mpi_allreduce(firststartl,firststart,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
      if(lrandsurf) then
        tnextrand = firststart
      else
        tnextrand = 9e9
      end if

    else                 
    ! ------------------------------------------------------
    ! Cold start -> load particle startpositions from txt

      np = 0
      startfile = 'partstartpos'
      open(ifinput,file=startfile,status='old',position='rewind',action='read')
      read(ifinput,'(I10.1)') np
      if ( np < 1 ) return
      ! read particles from partstartpos, create linked list
      do n = 1, np
        read(ifinput,*) tstart, xstart, ystart, zstart
        if(xstart < 0. .or. xstart > nxg*deltax .or. ystart < 0. .or. ystart > nyg*deltay .or. zstart < 0. .or. zstart > zm(nzp-1)) then
          if (myid == 0) print *, '  ABORTING: particle initialized outsize domain'
           call appl_abort(0)
        else 
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
              particle%zprev          = particle%z
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
              particle%sigma2_sgs     = 0.
              particle%partstep       = 0
              if(tstart < firststartl) firststartl = tstart
            end if
          end if
        end if
      end do

      ! Collect first start from other procs
      call mpi_allreduce(firststartl,firststart,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)

      ! Set first dump times
      tnextdump = firststart
      tnextrand = firststart
      !tnextstat = 0
      nstatsamp = 0
    end if

    ipunique        = 1
    ipx             = 2
    ipy             = 3
    ipz             = 4
    ipzprev         = 5
    ipxstart        = 6
    ipystart        = 7
    ipzstart        = 8
    iptsart         = 9
    ipures          = 10
    ipvres          = 11
    ipwres          = 12
    ipures_prev     = 13
    ipvres_prev     = 14
    ipwres_prev     = 15
    ipusgs          = 16
    ipvsgs          = 17
    ipwsgs          = 18
    ipusgs_prev     = 19
    ipvsgs_prev     = 20
    ipwsgs_prev     = 21
    ipsigma2_sgs    = 22
    ipartstep       = 23
    nrpartvar       = ipartstep

    ! 1D arrays for online statistics 
    if(lpartstat) then
      allocate(npartprof(nzp),npartprofl(nzp),     &
                   uprof(nzp),    uprofl(nzp),     &
                   vprof(nzp),    vprofl(nzp),     &
                   wprof(nzp),    wprofl(nzp),     &
                   u2prof(nzp),   u2profl(nzp),    &
                   v2prof(nzp),   v2profl(nzp),    &
                   w2prof(nzp),   w2profl(nzp),    &
                   tkeprof(nzp),  tkeprofl(nzp),   &
                   tprof(nzp),    tprofl(nzp),     &
                   tvprof(nzp),   tvprofl(nzp),    &
                   rtprof(nzp),   rtprofl(nzp),    &
                   rlprof(nzp),   rlprofl(nzp),    &
                   ccprof(nzp),   ccprofl(nzp))

      npartprof      = 0.
      npartprofl     = 0.
      uprof          = 0.
      vprof          = 0.
      wprof          = 0.
      u2prof         = 0.
      v2prof         = 0.
      w2prof         = 0.
      tkeprof        = 0.
      tprof          = 0. 
      tvprof         = 0. 
      rtprof         = 0. 
      rlprof         = 0. 
      ccprof         = 0. 
      uprofl         = 0. 
      vprofl         = 0.
      wprofl         = 0.
      u2profl        = 0.
      v2profl        = 0.
      w2profl        = 0.
      tkeprofl       = 0.
      tprofl         = 0.
      tvprofl        = 0.
      rtprofl        = 0.
      rlprofl        = 0.
      ccprofl        = 0.

      if(lpartsgs) then
        allocate(sigma2prof(nzp),sigma2profl(nzp), &
                 fsprof(nzp),    fsprofl(nzp))
        sigma2prof   = 0.
        sigma2profl  = 0.
        fsprof       = 0.
        fsprofl      = 0.
      end if
    end if  
    close(ifinput)

    if(lpartsgs)  allocate(sgse(nzp,nxp,nyp),rese(nzp,nxp,nyp),fs_local(nzp,nxp,nyp),fs(nzp))
    call init_random_seed()

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

    if(myid == 0) print "(//' ',49('-')/,' ',/,'  Lagrangian particles removed.')"
  end subroutine exit_particles

  !
  !--------------------------------------------------------------------------
  ! subroutine write_hist; writes history files for warm restart particles
  !   called from: init.f90, step.f90
  !--------------------------------------------------------------------------
  !
  subroutine write_particle_hist(htype, time)
    use mpi_interface, only : myid,wrxid,wryid
    use grid,          only : filprf
    implicit none
    integer, intent (in) :: htype
    real, intent (in)    :: time
    character (len=80)   :: hname
    type (particle_record), pointer:: particle
    integer              :: iblank

    write(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
    hname = trim(hname)//'.'//trim(filprf)//'.particles'

    select case(htype)
    case default
       hname = trim(hname)//'.iflg'
    case(0)
       hname = trim(hname)//'.R'
    case(1) 
       hname = trim(hname)//'.rst'
    case(2) 
       iblank=index(hname,' ')
       write (hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(time), 's'
    end select

    open(666,file=trim(hname), form='unformatted')
    write(666) np,tnextdump
    particle => head
    do while(associated(particle))
      write(666) particle%unique, particle%tstart, particle%partstep, & 
        particle%x, particle%xstart, particle%ures, particle%ures_prev, & 
        particle%y, particle%ystart, particle%vres, particle%vres_prev, & 
        particle%z, particle%zstart, particle%zprev, particle%wres, particle%wres_prev, &
        particle%usgs,      particle%vsgs,      particle%wsgs, &
        particle%usgs_prev, particle%vsgs_prev, particle%wsgs_prev, &
        particle%sigma2_sgs
      particle => particle%next
    end do 
    close(666)

  end subroutine write_particle_hist

  !
  !--------------------------------------------------------------------------
  ! subroutine initparticledump : creates NetCDF file for particle dump. 
  !   Called from: init.f90
  !--------------------------------------------------------------------------
  !
  subroutine initparticledump(time)
    use modnetcdf,       only : open_nc, addvar_nc
    use grid,            only : nzp, tname, tlongname, tunit, filprf
    use mpi_interface,   only : myid, nxprocs, nyprocs, wrxid, wryid
    use grid,            only : tname, tlongname, tunit, filprf
    implicit none

    real, intent(in)                  :: time
    character (40), dimension(2)      :: dimname, dimlongname, dimunit
    real, allocatable, dimension(:,:) :: dimvalues
    integer, dimension(2)             :: dimsize
    integer                           :: k,nlocal
    integer, parameter                :: precis = 0
    character (len=80)                :: hname

    nlocal = floor(real(np) / (nxprocs*nyprocs))
    if(myid .eq. nxprocs*nyprocs-1) then
      nlocal = np - (nxprocs*nyprocs-1) * nlocal
    end if 

    allocate(dimvalues(nlocal,2))

    dimvalues      = 0
    do k=1,nlocal
      dimvalues(k,1) = k
    end do

    dimname(1)     = 'particles'
    dimlongname(1) = 'ID of particle'
    dimunit(1)     = '-'    
    dimsize(1)     = nlocal
    dimname(2)     = tname
    dimlongname(2) = tlongname
    dimunit(2)     = tunit
    dimsize(2)     = 0

    write(hname,'(i4.4,i4.4)') wrxid,wryid
    hname = trim(filprf)//'.particles.'//trim(hname)//'.nc'
 
    call open_nc(hname, ncpartid, ncpartrec, time, .false.)
    call addvar_nc(ncpartid,'x','x-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
    call addvar_nc(ncpartid,'y','y-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
    call addvar_nc(ncpartid,'z','z-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
    if(lpartdumpui) then
      call addvar_nc(ncpartid,'u','resolved u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      call addvar_nc(ncpartid,'v','resolved v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      call addvar_nc(ncpartid,'w','resolved w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      if(lpartsgs) then
        call addvar_nc(ncpartid,'us','subgrid u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
        call addvar_nc(ncpartid,'vs','subgrid v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
        call addvar_nc(ncpartid,'ws','subgrid w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      end if
    end if
    if(lpartdumpth) then
      call addvar_nc(ncpartid,'t','liquid water potential temperature','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      call addvar_nc(ncpartid,'tv','virtual potential temperature','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
    end if
    if(lpartdumpmr) then
      call addvar_nc(ncpartid,'rt','total water mixing ratio','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
      call addvar_nc(ncpartid,'rl','liquid water mixing ratio','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
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
      call open_nc(trim(filprf)//'.particlestat.nc', ncpartstatid, ncpartstatrec, time, .false.)
      call addvar_nc(ncpartstatid,'np','Number of particles','-',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'u','resolved u-velocity of particle','m s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'v','resolved v-velocity of particle','m s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'w','resolved w-velocity of particle','m s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'u_2','resolved u-velocity variance of particle','m2 s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'v_2','resolved v-velocity variance of particle','m2 s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'w_2','resolved w-velocity variance of particle','m2 s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'tke','resolved TKE of particle','m s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'t','liquid water potential temperature','K',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'tv','virtual potential temperature','K',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'rt','total water mixing ratio','kg kg-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'rl','liquid water mixing ratio','kg kg-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
      call addvar_nc(ncpartstatid,'cc','cloud fraction','-',dimname,dimlongname,dimunit,dimsize,dimvalues)
      if(lpartsgs) then
        call addvar_nc(ncpartstatid,'fs','fraction subgrid TKE','-',dimname,dimlongname,dimunit,dimsize,dimvalues)
        call addvar_nc(ncpartstatid,'sgstke','subgrid TKE of particle','m s-1',dimname,dimlongname,dimunit,dimsize,dimvalues)
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

  end subroutine add_particle

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

  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
 
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
 
    deallocate(seed)
  end subroutine init_random_seed


  !-------- ARCHIVE ---------------------

  !
  !--------------------------------------------------------------------------
  ! subroutine initparticledump : creates NetCDF file for particle dump. 
  !   Called from: init.f90
  !--------------------------------------------------------------------------
  !
  !subroutine initparticledump(time)
  !  use modnetcdf,       only : open_nc, addvar_nc
  !  use grid,            only : nzp, tname, tlongname, tunit, filprf
  !  use mpi_interface,   only : myid
  !  use grid,            only : tname, tlongname, tunit
  !  implicit none

  !  real, intent(in)                  :: time
  !  character (40), dimension(2)      :: dimname, dimlongname, dimunit
  !  real, allocatable, dimension(:,:) :: dimvalues
  !  integer, dimension(2)             :: dimsize
  !  integer                           :: k
  !  integer, parameter                :: precis = 0

  !  allocate(dimvalues(np,2))

  !  dimvalues      = 0
  !  do k=1,np
  !    dimvalues(k,1) = k
  !  end do

  !  dimname(1)     = 'particles'
  !  dimlongname(1) = 'ID of particle'
  !  dimunit(1)     = '-'    
  !  dimsize(1)     = np
  !  dimname(2)     = tname
  !  dimlongname(2) = tlongname
  !  dimunit(2)     = tunit
  !  dimsize(2)     = 0
 
  !  if(myid == 0) then
  !    call open_nc(trim(filprf)//'.particles.nc', ncpartid, ncpartrec, time, .false.)
  !    call addvar_nc(ncpartid,'x','x-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !    call addvar_nc(ncpartid,'y','y-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !    call addvar_nc(ncpartid,'z','z-position of particle','m',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !    if(lpartdumpui) then
  !      call addvar_nc(ncpartid,'u','resolved u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !      call addvar_nc(ncpartid,'v','resolved v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !      call addvar_nc(ncpartid,'w','resolved w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !      if(lpartsgs) then
  !        call addvar_nc(ncpartid,'us','subgrid u-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !        call addvar_nc(ncpartid,'vs','subgrid v-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !        call addvar_nc(ncpartid,'ws','subgrid w-velocity of particle','m/s',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !      end if
  !    end if
  !    if(lpartdumpth) then
  !      call addvar_nc(ncpartid,'t','liquid water potential temperature','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !      call addvar_nc(ncpartid,'tv','virtual potential temperature','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !    end if
  !    if(lpartdumpmr) then
  !      call addvar_nc(ncpartid,'rt','total water mixing ratio','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !      call addvar_nc(ncpartid,'rl','liquid water mixing ratio','K',dimname,dimlongname,dimunit,dimsize,dimvalues,precis)
  !    end if
  !  end if 
 
  !end subroutine initparticledump

  !
  !--------------------------------------------------------------------------
  ! subroutine particledump : communicates all particles to process #0 and
  !   writes it to NetCDF(4)
  !--------------------------------------------------------------------------
  !
  !subroutine particledump(time)
  !  use mpi_interface, only : mpi_comm_world, myid, mpi_integer, mpi_double_precision, ierror, nxprocs, nyprocs, mpi_status_size, wrxid, wryid, nxg, nyg
  !  use grid,          only : tname, deltax, deltay, dzi_t, zm, umean, vmean
  !  use modnetcdf,     only : writevar_nc, fillvalue_double
  !  implicit none

  !  real, intent(in)                     :: time
  !  type (particle_record), pointer:: particle
  !  integer                              :: nlocal, ii, i, pid, partid
  !  integer, allocatable, dimension(:)   :: nremote
  !  integer                              :: status(mpi_status_size)
  !  integer                              :: nvar,nvl
  !  real, allocatable, dimension(:)      :: sendbuff, recvbuff
  !  real, allocatable, dimension(:,:)    :: particles_merged
  !  real                                 :: thl,thv,rt,rl

  !  nvar = 4                            ! id,x,y,z
  !  if(lpartdumpui)  nvar = nvar + 3    ! u,v,w
  !  if(lpartsgs)     nvar = nvar + 3    ! us,vs,ws
  !  if(lpartdumpth)  nvar = nvar + 2    ! thl,tvh
  !  if(lpartdumpmr)  nvar = nvar + 2    ! rt,rl

  !  ! Count local particles
  !  nlocal = 0
  !  particle => head
  !  do while( associated(particle) )
  !    nlocal = nlocal + 1
  !    particle => particle%next
  !  end do

  !  ! Communicate number of local particles to main proces (0)
  !  allocate(nremote(0:(nxprocs*nyprocs)-1))
  !  nremote = 0
  !  call mpi_gather(nlocal,1,mpi_integer,nremote,1,mpi_integer,0,mpi_comm_world,ierror)

  !  !! Create buffer
  !  allocate(sendbuff(nvar * nlocal))
  !  ii = 1
  !  particle => head
  !  do while( associated(particle) )
  !    if(lpartdumpth .or. lpartdumpmr) call thermo(particle%x,particle%y,particle%z,thl,thv,rt,rl)

  !    sendbuff(ii)   = particle%unique
  !    sendbuff(ii+1) = (wrxid * (nxg / nxprocs) + particle%x - 3) * deltax
  !    sendbuff(ii+2) = (wryid * (nyg / nyprocs) + particle%y - 3) * deltay
  !    sendbuff(ii+3) = zm(floor(particle%z)) + (particle%z-floor(particle%z)) / dzi_t(floor(particle%z))
  !    nvl = 3
  !    if(lpartdumpui) then
  !      sendbuff(ii+nvl+1) = particle%ures * deltax
  !      sendbuff(ii+nvl+2) = particle%vres * deltay
  !      sendbuff(ii+nvl+3) = particle%wres / dzi_t(floor(particle%z))
  !      nvl = nvl + 3
  !      if(lpartsgs) then
  !        sendbuff(ii+nvl+1) = particle%usgs * deltax
  !        sendbuff(ii+nvl+2) = particle%vsgs * deltay
  !        sendbuff(ii+nvl+3) = particle%wsgs / dzi_t(floor(particle%z))
  !        nvl = nvl + 3
  !      end if
  !    end if
  !    if(lpartdumpth) then
  !      sendbuff(ii+nvl+1) = thl
  !      sendbuff(ii+nvl+2) = thv
  !      nvl = nvl + 2
  !    end if
  !    if(lpartdumpmr) then
  !      sendbuff(ii+nvl+1) = rt
  !      sendbuff(ii+nvl+2) = rl
  !    end if
  !    ii = ii + nvar
  !    particle => particle%next
  !  end do

  !  ! Dont send when main process
  !  if(myid .ne. 0) call mpi_send(sendbuff,nvar*nlocal,mpi_double_precision,0,myid,mpi_comm_world,ierror)

  !  if(myid .eq. 0) then
  !    allocate(particles_merged(np,nvar-1))

  !    ! Add local particles
  !    ii = 1
  !    do i=1,nlocal
  !      partid = int(sendbuff(ii))
  !      particles_merged(partid,1) = sendbuff(ii+1)
  !      particles_merged(partid,2) = sendbuff(ii+2)
  !      particles_merged(partid,3) = sendbuff(ii+3)
  !      nvl = 3
  !      if(lpartdumpui) then
  !        particles_merged(partid,nvl+1) = sendbuff(ii+nvl+1)
  !        particles_merged(partid,nvl+2) = sendbuff(ii+nvl+2)
  !        particles_merged(partid,nvl+3) = sendbuff(ii+nvl+3)
  !        nvl = nvl + 3
  !        if(lpartsgs) then
  !          particles_merged(partid,nvl+1) = sendbuff(ii+nvl+1)
  !          particles_merged(partid,nvl+2) = sendbuff(ii+nvl+2)
  !          particles_merged(partid,nvl+3) = sendbuff(ii+nvl+3)
  !          nvl = nvl + 3
  !        end if
  !      end if
  !      if(lpartdumpth) then
  !        particles_merged(partid,nvl+1) = thl
  !        particles_merged(partid,nvl+2) = thv
  !        nvl = nvl + 2
  !      end if
  !      if(lpartdumpmr) then
  !        particles_merged(partid,nvl+1) = rt
  !        particles_merged(partid,nvl+2) = rl
  !      end if
  !      ii = ii + nvar 
  !    end do 

  !    ! Add remote particles
  !    do pid = 1,(nxprocs*nyprocs)-1
  !      allocate(recvbuff(nremote(pid)*nvar))
  !      call mpi_recv(recvbuff,nremote(pid)*nvar,mpi_double_precision,pid,pid,mpi_comm_world,status,ierror)   
  !      ii = 1
  !      do i=1,nremote(pid)
  !        partid = int(recvbuff(ii))
  !        particles_merged(partid,1) = recvbuff(ii+1)
  !        particles_merged(partid,2) = recvbuff(ii+2)
  !        particles_merged(partid,3) = recvbuff(ii+3)
  !        nvl = 3
  !        if(lpartdumpui) then
  !          particles_merged(partid,nvl+1) = recvbuff(ii+nvl+1)
  !          particles_merged(partid,nvl+2) = recvbuff(ii+nvl+2)
  !          particles_merged(partid,nvl+3) = recvbuff(ii+nvl+3)
  !          nvl = nvl + 3
  !          if(lpartsgs) then
  !            particles_merged(partid,nvl+1) = recvbuff(ii+nvl+1)
  !            particles_merged(partid,nvl+2) = recvbuff(ii+nvl+2)
  !            particles_merged(partid,nvl+3) = recvbuff(ii+nvl+3)
  !            nvl = nvl + 3
  !          end if
  !        end if
  !        if(lpartdumpth) then
  !          particles_merged(partid,nvl+1) = recvbuff(ii+nvl+1)
  !          particles_merged(partid,nvl+2) = recvbuff(ii+nvl+2)
  !          nvl = nvl + 2
  !        end if
  !        if(lpartdumpmr) then
  !          particles_merged(partid,nvl+1) = recvbuff(ii+nvl+1)
  !          particles_merged(partid,nvl+2) = recvbuff(ii+nvl+2)
  !        end if
  !        ii = ii + nvar 
  !      end do 
  !      deallocate(recvbuff)
  !    end do

  !    ! Correct for Galilean transformation
  !    ! Subgrid motion is assumed to be centered around the resolved velocity -> no galilean tranformation
  !    if(lpartdumpui) then
  !      particles_merged(:,4) = particles_merged(:,4) + umean
  !      particles_merged(:,5) = particles_merged(:,5) + vmean
  !    end if

  !    ! Write to NetCDF
  !    call writevar_nc(ncpartid,tname,time,ncpartrec)
  !    call writevar_nc(ncpartid,'x',particles_merged(:,1),ncpartrec)
  !    call writevar_nc(ncpartid,'y',particles_merged(:,2),ncpartrec)
  !    call writevar_nc(ncpartid,'z',particles_merged(:,3),ncpartrec)
  !    nvl = 3
  !    if(lpartdumpui) then
  !      call writevar_nc(ncpartid,'u',particles_merged(:,nvl+1),ncpartrec)
  !      call writevar_nc(ncpartid,'v',particles_merged(:,nvl+2),ncpartrec)
  !      call writevar_nc(ncpartid,'w',particles_merged(:,nvl+3),ncpartrec)
  !      nvl = nvl + 3
  !      if(lpartsgs) then
  !        call writevar_nc(ncpartid,'us',particles_merged(:,nvl+1),ncpartrec)
  !        call writevar_nc(ncpartid,'vs',particles_merged(:,nvl+2),ncpartrec)
  !        call writevar_nc(ncpartid,'ws',particles_merged(:,nvl+3),ncpartrec)
  !        nvl = nvl + 3
  !      end if
  !    end if
  !    if(lpartdumpth) then
  !      call writevar_nc(ncpartid,'t', particles_merged(:,nvl+1),ncpartrec)
  !      call writevar_nc(ncpartid,'tv',particles_merged(:,nvl+2),ncpartrec)
  !      nvl = nvl + 2
  !    end if
  !    if(lpartdumpmr) then
  !      call writevar_nc(ncpartid,'rt',particles_merged(:,nvl+1),ncpartrec)
  !      call writevar_nc(ncpartid,'rl',particles_merged(:,nvl+2),ncpartrec)
  !    end if
  !  end if
  !  
  !  if(myid==0) deallocate(particles_merged)
  !  deallocate(nremote, sendbuff)
 
  !end subroutine particledump

  !
  !--------------------------------------------------------------------------
  ! Subroutine randomize 
  !
  !> !!!! DEPRECATED !!!!!
  !
  !> Randomizes the X,Y,Z positions of all particles in 
  !> the lowest grid level every RK3 cycle. Called from: particles()
  !--------------------------------------------------------------------------
  !
  !subroutine randomize(once)
  !  use mpi_interface, only : nxg, nyg, nyprocs, nxprocs
  !  implicit none

  !  logical, intent(in) :: once            !> flag: randomize all particles (false) or only the onces which sink into the surface layer (true)?  
  !  real                :: zmax = 1.       ! Max height in grid coordinates
  !  integer             :: nyloc, nxloc 
  !  type (particle_record), pointer:: particle
  !  real                :: randnr(3)

  !  nyloc   = nyg / nyprocs
  !  nxloc   = nxg / nxprocs

  !  if(once) then 
  !    particle => head
  !    do while(associated(particle) )
  !      if( particle%z <= (1. + zmax) .and. particle%zprev > (1. + zmax) ) then
  !        call random_number(randnr)
  !        particle%x = (randnr(1) * nyloc) + 3 
  !        particle%y = (randnr(2) * nxloc) + 3
  !        !particle%z = zmax * randnr(3)    + 1
  !        particle%ures_prev = 0.
  !        particle%vres_prev = 0.
  !        particle%wres_prev = 0.
  !      end if
  !      particle => particle%next
  !    end do 
  !  else
  !    particle => head
  !    do while(associated(particle) )
  !      if( particle%z <= (1. + zmax) ) then
  !        call random_number(randnr)
  !        particle%x = (randnr(1) * nyloc) + 3 
  !        particle%y = (randnr(2) * nxloc) + 3
  !        !particle%z = zmax * randnr(3)    + 1 
  !      end if
  !      particle => particle%next
  !    end do 
  !  end if 
 
  !end subroutine randomize

end module modparticles
