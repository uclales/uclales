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
! This is an adaptation from the original Independent Pixel Code.
! In the IPA Code the spatial dimensions are the outermost loops.
! To allow for 3D approximations and solvers, we need to rearrange the loop structure,
! so that the spectral integration is outermost.
! The 1D routines are reused as much as possible to guarantee conformance to earlier results.
! The 3D interface can be called with iradtype:
! 6 -- 1D original delta 4 stream solver
! 7 -- thermal and solar tenstream solver
! The copying of the optical properties fields results in a performance penalty 
! -- hence if you do not use a 3D solver, please use the 1D interface with iradtype=4

module radiation_3d

  use grid, only       : nzp,nxp,nyp, deltax,deltay,zm, a_rhl,a_rhs
  use fuliou, only     : minSolarZenithCosForVis
  use mpi_interface, only : myid,ierror
  use radiation, only  : zenith, setup, pp, pt, ph, po, pre, pde, plwc, piwc, prwc, pgwc, &
      u0, fixed_sun, radMcICA, rad_eff_radius !namelist parameters

#ifdef HAVE_TENSTREAM
  use mpi_interface, only: nxpa,nypa
  use m_tenstream, only: init_tenstream,set_optical_properties,solve_tenstream, destroy_tenstream, tenstream_get_result, need_new_solution
  use m_data_parameters, only : ireals,iintegers
  use grid, only : dt,nstep
#endif      

#ifdef _XLF
  use ieee_arithmetic 
#define isnan ieee_is_nan
#endif

  implicit none

  private 
  public :: rad_3d

  character (len=10), parameter :: background = 'backrad_in'

  !those global vars are in radiation module...  ::  pp, pt, ph, po, pre, pde, plwc, piwc, prwc, pgwc
  real,parameter :: zero=0, one=1

  integer :: npts
  real    :: ee, day, time, alat, zz
  integer :: is,ie,js,je

  logical,save :: linit=.False.
  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.

#ifdef HAVE_TENSTREAM
  integer(iintegers) :: solution_uid    ! is solution uid, each subband has one
  real(ireals)       :: solution_time   ! is set to the approximate time of each solve
#endif

contains

  subroutine rad_3d(alat, time, sknt, sfc_albedo, CCN, dn0, &
          pi0, pi1, pip, th, rv, rc, tt, rflx, sflx ,lflxu, lflxd,sflxu,sflxd, &
          albedo, lflxu_toa, lflxd_toa, sflxu_toa, sflxd_toa, rr,ice,nice,grp)

      real, intent (in)    :: alat, time, sknt, sfc_albedo, CCN
      real, dimension (nzp), intent (in)                   :: dn0, pi0, pi1
      real, dimension (nzp,nxp,nyp), intent (in)           :: pip, th, rv, rc
      real, optional, dimension (nzp,nxp,nyp), intent (in) :: rr,ice,nice,grp
      real, dimension (nzp,nxp,nyp), intent (inout)        :: tt, rflx, sflx, lflxu, lflxd, sflxu, sflxd
      real, dimension (nxp,nyp), intent (out),optional     :: albedo, lflxu_toa, lflxd_toa, sflxu_toa, sflxd_toa

      real,allocatable,dimension(:,:,:) :: fus,fds,fdiv_sol
      real,allocatable,dimension(:,:,:) :: fuir,fdir,fdiv_th

      real,allocatable,dimension(:,:,:) :: hr_factor, dz ! convert from flux divergence to heating rate

      real :: p0(nzp)

      integer i,j,k,kk,ierr

#ifdef HAVE_TENSTREAM
      !if(.not.radMcICA.and.nstep.ne.1) return
#endif

      call init_rad_3d()

      if(u0.ge.minSolarZenithCosForVis) &
          call solar_rad()

      call thermal_rad()

      if(ldebug.and.myid.le.0) then
        do k=1,ubound(fds,1)
          print *,k,'solar',fds(k,3,3),fus(k,3,3),fdiv_sol(k,3,3),':: thermal',fdir(k,3,3),fuir(k,3,3),fdiv_th(k,3,3)
        enddo
      endif

      !copy from radiation grid, to dynamics grid
      do k=1,nzp
        kk = ubound(fus,1) - (k-1)
        sflx(k,is:ie,js:je) = fus(kk,:,:)  - fds(kk,:,:)
        !irina
        sflxu(k,is:ie,js:je)=fus (kk,:,:)
        sflxd(k,is:ie,js:je)=fds (kk,:,:)
        lflxu(k,is:ie,js:je)=fuir(kk,:,:)
        lflxd(k,is:ie,js:je)=fdir(kk,:,:)

        rflx  (k,is:ie,js:je) = sflx    (k ,is:ie,js:je) + fuir(kk,:,:) - fdir(kk,:,:)
        a_rhl(k,is:ie,js:je) = fdiv_th (kk,is:ie,js:je) * hr_factor(k, :,:)
        a_rhs(k,is:ie,js:je) = fdiv_sol(kk,is:ie,js:je) * hr_factor(k, :,:)

      end do

      if (present(albedo)) then
        if (u0 > minSolarZenithCosForVis) then
          albedo(is:ie,js:je) = fus(1,:,:)/fds(1,:,:)
        else
          albedo = -999.
        end if
      end if

      if (present(sflxu_toa)) then
        if (u0 > minSolarZenithCosForVis) then
          sflxu_toa(is:ie,js:je) = fus(1,:,:)
        else
          sflxu_toa(is:ie,js:je) = -999.
        end if
      end if
      if (present(sflxd_toa)) then
        if (u0 > minSolarZenithCosForVis) then
          sflxd_toa(is:ie,js:je) = fds(1,:,:)
        else
          sflxd_toa(is:ie,js:je) = -999.
        end if
      end if
      if (present(lflxu_toa)) then
        lflxu_toa(is:ie,js:je) = fuir(1,:,:)
      end if
      if (present(lflxd_toa)) then
        lflxd_toa(is:ie,js:je) = fdir(1,:,:)
      end if
      !TODO here the loop was truncated to 'nzp-3' in the original code -- why not use heating rate in every layer?!?
      do k=2,nzp-3
        tt(k,is:ie, js:je) = tt(k, is:ie, js:je) + a_rhl(k,is:ie,js:je) + a_rhs(k,is:ie,js:je)
#ifndef _XLF
        if(ldebug) then
          if(any(isnan(rflx     (k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency rflx  ',k,rflx    (k,:,:), any(isnan(sflx)),any(isnan(sflxu)),any(isnan(sflxd)),any(isnan(lflxu)),any(isnan(lflxd))
          if(any(isnan(fdiv_sol (k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency divsol',k,fdiv_sol(k,:,:)
          if(any(isnan(fdiv_th  (k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency divth ',k,fdiv_th (k,:,:)
          if(any(isnan(hr_factor(k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency hrf',k,hr_factor(k,:,:)
          if(any(isnan(a_rhl    (k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency rhl',k,a_rhl(k,:,:)
          if(any(isnan(a_rhs    (k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency rhs',k,a_rhs(k,:,:)
          if(any(isnan(tt       (k,:,:)))) print *,myid,'rad3d :: nan in radiation tendency tt',k,tt(k,:,:)
          if(any(isnan([hr_factor(k,:,:),a_rhl(k,:,:),a_rhs(k,:,:),tt(k,:,:),fdiv_sol (k,:,:),fdiv_th  (k,:,:),rflx     (k,:,:)]))) call exit(1)
        endif
#endif
      end do

      if(myid.eq.0.and.ldebug) print *,'calculate radiation ... done'
    contains
      subroutine thermal_rad()
          use grid, only       : iradtyp, a_tskin
          use ckd   , only: llimit, rlimit
          use fuliou, only: computeIRBandWeights, planck, select_bandg
          use defs, only: nv,nv1,pi
          use solver, only : qft
          use ckd, only: ir_bands,gPointWeight,kg
          use mpi, only    : mpi_comm_world,mpi_integer

          real,dimension(nv1)  :: fu1,fd1
          real,dimension(nv1,  is:ie,js:je) :: fd3d,fu3d,fdiv3d
          real,dimension(nv ,  is:ie,js:je) :: tau,w0
          real,dimension(nv ,4,is:ie,js:je) :: phasefct

          real,dimension(nv1,  is:ie,js:je) :: bf

          real, dimension(:), allocatable, save :: bandWeights

          real :: fuq2,xir_norm
          logical, parameter :: irWeighted = .False.

          real    :: randomNumber

          integer :: iband, igpt,nrbands,nrgpts, ib, ig
          integer :: ibandloop(3:nxp-2, 3:nyp-2), ibandg(3:nxp-2, 3:nyp-2) ! for McICA, save wavelength band information for each pixel


#ifdef HAVE_TENSTREAM
          real(ireals),dimension(:,:,:),allocatable :: edn,eup
          real(ireals),dimension(:,:,:),allocatable :: abso
          integer :: ierr
#endif

          real,allocatable,dimension(:,:,:) :: fuir_nca,fdir_nca,fdiv_th_nca,tau_nca

          if(iradtyp.eq.8) then
            allocate(fdiv_th_nca(nv1,is:ie,js:je) )
            allocate(fdir_nca   (nv1,is:ie,js:je) )
            allocate(fuir_nca   (nv1,is:ie,js:je) )
            allocate(tau_nca    (nv1,is:ie,js:je) )
            fdiv_th_nca=0
            fdir_nca   =0
            fuir_nca   =0
            tau_nca    =0
          endif

          fdir    = 0 
          fuir    = 0 
          fdiv_th = 0 

          if (.not. allocated(bandweights)) then 
            allocate(bandweights(size(ir_bands)))
            call computeIRBandWeights(ir_bands, irWeighted, bandWeights)
          end if

          if (radMcICA) then
            do j=js,je
              do i=is,ie
                call random_number(randomNumber)
                !
                ! Select a single band and g-point (ib, ig1) and use these as the 
                ! limits in the loop through the spectrum below. 
                !
                call select_bandg(ir_bands, bandweights, randomNumber, ib, ig)
                ibandloop(i,j) = ib
                ibandg   (i,j) = ig
                nrbands = 1
                nrgpts  = 1
              enddo 
            enddo
          else
            nrbands = size(ir_bands)
          end if
          do iband = 1, nrbands
            if(.not.radMcICA) nrgpts = kg(ir_bands(iband))
            do igpt = 1, nrgpts

              if(radMcICA) then
                if(iradtyp.eq.7 .or. iradtyp.eq.8) then
                  ! tenstream -- only use one random band.... due to
                  ! horizontal correlations, it does not make sense to have
                  !different spectral intervals horizontally
                  call mpi_bcast(ibandloop(is,js),1,mpi_integer,0,mpi_comm_world,ierr)
                  call mpi_bcast(ibandg   (is,js),1,mpi_integer,0,mpi_comm_world,ierr)
                  ib = ibandloop(is,js)
                  ig = ibandg   (is,js)
                else
                  ! will be set below, dont need them for solution_uid...
                  ib = -1
                  ig = -1
                endif
              else
                ib = iband
                ig = igpt
              endif

#ifdef HAVE_TENSTREAM
              if(iradtyp.eq.7) then
                solution_uid=get_band_uid(.False., ib,ig)
                solution_time = time*3600._ireals*24._ireals + dt*(nstep-1._ireals)/3._ireals !time is given in days + approx. a third at each rungekutta step
                if(.not.need_new_solution(solution_uid,solution_time)) then
                  call load_tenstream_solution(.False.,dz,u0,solution_uid,fd3d,fu3d,fdiv3d)

                  if (radMcICA) then
                    xir_norm = 1./bandweights(ib)
                  else
                    xir_norm = gPointWeight(ir_bands(ib), ig)
                  end if

                  fdir = fdir       + fd3d  *xir_norm
                  fuir = fuir       + fu3d  *xir_norm
                  fdiv_th = fdiv_th + fdiv3d*xir_norm
                  cycle ! if we successfully loaded a solution, just cycle this spectral band
                endif ! dont need to calc new solution
              endif
#endif

              do j=js,je
                do i=is,ie

                  if(radMcICA) then
                    select case(iradtyp)
                    case (6) ! d4stream with 3d interface
                      ib = ibandloop(i,j)
                      ig = ibandg   (i,j)
                    case (7) ! tenstream -- only use one random band.... due to
                      ! horizontal correlations, it does not make sense to have
                      !different spectral intervals horizontally
                      ib = ibandloop(is,js)
                      ig = ibandg   (is,js)
                    end select
                  else
                    ib = ib
                    ig = ig
                  endif

                  if (present(ice).and.present(grp)) then
                    call setup_rad_atmosphere(CCN, dn0, pi0, pi1,    &
                        pip(:,i,j), th(:,i,j), rv(:,i,j), rc(:,i,j),     &
                        hr_factor(:,i,j), dz(:,i,j),         &
                        rr=rr(:,i,j),ice=ice(:,i,j),nice=nice(:,i,j),grp=grp(:,i,j))
                    call optprop_rad_ir( ib, ig, pp, pt, ph, po, &
                        tau (:,i,j), w0  (:,i,j), phasefct(:,:,i,j),dz(:,i,j),  & 
                        plwc=plwc, pre=pre, piwc=piwc, pde=pde, pgwc=pgwc)
                  else
                    call setup_rad_atmosphere(CCN, dn0, pi0, pi1,    &
                        pip(:,i,j), th(:,i,j), rv(:,i,j), rc(:,i,j),     &
                        hr_factor(:,i,j), dz(:,i,j)         )
                    call optprop_rad_ir( ib, ig, pp, pt, ph, po,&
                        tau (:,i,j), w0  (:,i,j), phasefct(:,:,i,j),dz(:,i,j), & 
                        plwc=plwc, pre=pre)
                  end if

                  if (allocated(a_tskin)) then ! have lsm to give surface temp
                    call planck(pt, a_tskin(i,j), llimit(ir_bands(ib)), rlimit(ir_bands(ib)), bf(:,i,j))
                  else
                    call planck(pt, sknt, llimit(ir_bands(ib)), rlimit(ir_bands(ib)), bf(:,i,j))
                  endif

                end do ! j
              end do ! i


              select case(iradtyp)

              case (6,8) ! d4stream with 3d interface

                if(iradtyp.eq.8) tau_nca = tau

                do j = js, je  
                  do i = is, ie  

                    ! Solver expects cumulative optical depth
                    do k = 2, nv

                      tau(k,i,j) = tau(k, i,j) + tau(k-1, i,j)
                    end do
                    call qft (.false., ee, zero, zero, bf(:,i,j), tau(:,i,j), w0(:,i,j),                    & ! thermal qft
                        phasefct(:, 1, i,j), phasefct(:, 2, i,j), phasefct(:, 3, i,j), phasefct(:, 4, i,j), &
                        fu1, fd1)

                    if (radMcICA) then
                      ib = ibandloop(i,j)
                      ig = ibandg   (i,j)
                      xir_norm = 1./bandweights(ib)
                    else
                      xir_norm = gPointWeight(ir_bands(ib), ig)
                    end if

                    if(iradtyp.eq.8) then
                      fdir_nca(:,i,j) =fd1(:)
                      fuir_nca(:,i,j) =fu1(:)
                    endif

                    fdir(:,i,j) = fdir(:,i,j) + fd1(:) * xir_norm 
                    fuir(:,i,j) = fuir(:,i,j) + fu1(:) * xir_norm 
                  end do ! i
                end do ! j

                fdiv_th(1:nv,:,:) = (fdir(1:nv,:,:) - fuir(1:nv,:,:)) + (fuir(2:nv1,:,:) - fdir(2:nv1,:,:))
                fdiv_th( nv1,:,:) = (fdir( nv1,:,:) - fuir(nv1,:,:))


                if(iradtyp.eq.8) then
                  call nca_wrapper(nxp, nyp, nv, deltax, deltay, dz, tau_nca, w0, bf, fdir_nca , fuir_nca, fdiv3d)

                  !give 3d heating to fdiv_th, layers above nzp remain 1d heating
                  fdiv_th_nca(nv+2-nzp:nv,:,:) =  fdiv_th_nca(nv+2-nzp:nv,:,:) + fdiv3d(nv+2-nzp:nv,:,:)*xir_norm 
                endif


              case (7) !tenstr
#ifdef HAVE_TENSTREAM
                call tenstream_wrapper(.False., nxp,nyp,nv,deltax,deltay,dz, -one,-one, one-ee,zero, tau, w0, phasefct,bf, fd3d,fu3d,fdiv3d)
#else
                print *,'This build does not support the tenstream solver ... exiting!'
                call exit(1)
#endif

                if (radMcICA) then
                  xir_norm = 1./bandweights(ibandloop(is,js))
                else
                  xir_norm = gPointWeight(ir_bands(ib), ig)
                end if

                fdir = fdir       + fd3d  *xir_norm
                fuir = fuir       + fu3d  *xir_norm
                fdiv_th = fdiv_th + fdiv3d*xir_norm
              end select

            enddo !igpt
          enddo !iband

          if(iradtyp.eq.8) then
            fdiv_th(nv+2-nzp:nv,:,:) = fdiv_th_nca(nv+2-nzp:nv,:,:) 
          endif

      end subroutine
      subroutine solar_rad()
          use grid, only   : iradtyp
          use defs, only   : nv,nv1,SolarConstant,totalpower
          use solver, only : qft
          use fuliou, only : computesolarbandweights, select_bandg
          use ckd, only    : solar_bands,kg,power,gPointWeight
          use mpi, only    : mpi_comm_world,mpi_integer

!          use m_twostream, only: delta_eddington_twostream

          real :: fuq1,xs_norm
          real,dimension(nv1)  :: fu1,fd1
          real,dimension(nv1,  is:ie,js:je) :: fu3d,fd3d,fdiv3d
          real,dimension(nv ,  is:ie,js:je) :: tau,w0
          real,dimension(nv1,  is:ie,js:je) :: bf
          real,dimension(nv ,4,is:ie,js:je) :: phasefct
          real, dimension(:), allocatable, save :: bandWeights

          real    :: randomNumber
          logical, parameter :: solarWeighted = .true. 

          integer :: iband, igpt,nrbands,nrgpts, ib, ig
          integer :: ibandloop(3:nxp-2, 3:nyp-2), ibandg(3:nxp-2, 3:nyp-2) ! for McICA, save wavelength band information for each pixel

#ifdef HAVE_TENSTREAM
          real(ireals),dimension(:,:,:),allocatable :: edir,edn,eup
          real(ireals),dimension(:,:,:),allocatable :: abso
          integer :: ierr
#endif

          fus      = zero
          fds      = zero
          fdiv_sol = zero

          if(u0.le.minSolarZenithCosForVis) then
            return
          endif

          if (.not. allocated(bandweights)) then 
            allocate(bandweights(size(solar_bands)))
            call computeSolarBandWeights(solar_bands, solarWeighted, bandWeights)
          end if

          if (radMcICA) then
            do j=js,je
              do i=is,ie
                call random_number(randomNumber)
                !
                ! Select a single band and g-point (ib, ig1) and use these as the 
                ! limits in the loop through the spectrum below. 
                !
                call select_bandg(solar_bands, bandweights, randomNumber, ib, ig)
                ibandloop(i,j) = ib
                ibandg   (i,j) = ig
              enddo 
            enddo
            nrbands = 1
            nrgpts  = 1
          else
            nrbands = size(solar_bands)
          end if
          do iband = 1, nrbands
            if(.not.radMcICA) nrgpts = kg(solar_bands(iband))
            do igpt = 1, nrgpts

              if(radMcICA) then
                if(iradtyp.eq.7) then
                  ! tenstream -- only use one random band.... due to
                  ! horizontal correlations, it does not make sense to have
                  !different spectral intervals horizontally
                  call mpi_bcast(ibandloop(is,js),1,mpi_integer,0,mpi_comm_world,ierr)
                  call mpi_bcast(ibandg   (is,js),1,mpi_integer,0,mpi_comm_world,ierr)
                  ib = ibandloop(is,js)
                  ig = ibandg   (is,js)
                else
                  ! will be set below, dont need them for solution_uid...
                  ib = -1
                  ig = -1
                endif
              else
                ib = iband
                ig = igpt
              endif

#ifdef HAVE_TENSTREAM
              if(iradtyp.eq.7) then
                solution_uid = get_band_uid( .True., ib, ig )
                solution_time = time*3600._ireals*24._ireals + dt*(nstep-1._ireals)/3._ireals !time is given in days + approx. a third at each rungekutta step
                if(.not.need_new_solution(solution_uid,solution_time)) then
                  call load_tenstream_solution(.True.,dz,u0,solution_uid,fd3d,fu3d,fdiv3d)

                  fds     = fds      + fd3d  
                  fus     = fus      + fu3d  
                  fdiv_sol= fdiv_sol + fdiv3d
                  cycle ! cycle this spectral band
                endif ! dont need to calc new solution
              endif
#endif

              do j=js,je
                do i=is,ie

                  if(radMcICA) then
                    select case(iradtyp)
                    case (6) ! d4stream with 3d interface
                      ib = ibandloop(i,j)
                      ig = ibandg   (i,j)
                    case (7) ! tenstream -- only use one random band.... due to
                      ! horizontal correlations, it does not make sense to have
                      !different spectral intervals horizontally
                      ib = ibandloop(is,js)
                      ig = ibandg   (is,js)
                    end select
                  endif

                  if (present(ice).and.present(grp)) then
                    call setup_rad_atmosphere(CCN, dn0, pi0, pi1,    &
                        pip(:,i,j), th(:,i,j), rv(:,i,j), rc(:,i,j),     &
                        hr_factor(:,i,j), dz(:,i,j),         &
                        rr=rr(:,i,j),ice=ice(:,i,j),nice=nice(:,i,j),grp=grp(:,i,j))
                    call optprop_rad_vis( ib, ig, pp, pt, ph, po, &
                        tau (:,i,j), w0  (:,i,j), phasefct(:,:,i,j),dz(:,i,j), & 
                        plwc=plwc, pre=pre, piwc=piwc, pde=pde, pgwc=pgwc)
                  else
                    call setup_rad_atmosphere(CCN, dn0, pi0, pi1,    &
                        pip(:,i,j), th(:,i,j), rv(:,i,j), rc(:,i,j),     &
                        hr_factor(:,i,j), dz(:,i,j)         )
                    call optprop_rad_vis( ib, ig, pp, pt, ph, po,&
                        tau (:,i,j), w0  (:,i,j), phasefct(:,:,i,j),dz(:,i,j), & 
                        plwc=plwc, pre=pre)
                  end if

                  bf(:,i,j) = 0 ! no thermal emission in solar spectral range
                end do ! j
              end do ! i

              select case(iradtyp)

              case (6) ! d4stream with 3d interface

                do j = js, je  
                  do i = is, ie  

                    if(radMcICA) then
                      ib = ibandloop(i,j)
                      ig = ibandg   (i,j)
                    endif

                    ! Solver expects cumulative optical depth
                    do k = 2, nv
                      tau(k,i,j) = tau(k, i,j) + tau(k-1, i,j)
                    end do
                    call qft (.true., zero, sfc_albedo, u0, bf(:,i,j), tau(:,i,j), w0(:,i,j), &   ! Solar qft
                        phasefct(:, 1, i,j), phasefct(:, 2, i,j),    &
                        phasefct(:, 3, i,j), phasefct(:, 4, i,j), fu1, fd1)

                    if (radMcICA) then
                      xs_norm = power(solar_bands(ib))/ bandweights(ib)
                    else
                      xs_norm = gPointWeight(solar_bands(ib), ig)*power(solar_bands(ib))
                    end if
                    fd3d(:,i,j) = fd1(:) * xs_norm 
                    fu3d(:,i,j) = fu1(:) * xs_norm 
                    fdiv3d(1:nv,i,j) = ( (fd1(1:nv) - fu1(1:nv)) + (fu1(2:nv1) - fd1(2:nv1))   ) *xs_norm
                    fdiv3d( nv1,i,j) = ( (fd1(nv1) - fu1(nv1))                                 ) *xs_norm

                  end do ! i
                end do ! j

              case (7,8) !tenstr
                if (radMcICA) then
                  ib = ibandloop(is,js)
                  ig = ibandg   (is,js)
                  xs_norm = power(solar_bands(ib))/ bandweights(ib)
                else
                  xs_norm = gPointWeight(solar_bands(ib), ig)*power(solar_bands(ib))
                end if

#ifdef HAVE_TENSTREAM
                call tenstream_wrapper(.True., nxp,nyp,nv,deltax,deltay,dz, 180.,u0, sfc_albedo,xs_norm, tau, w0, phasefct,bf, fd3d,fu3d,fdiv3d)
#else
                print *,'This build does not support the tenstream solver ... exiting!'
                call exit(1)
#endif

              end select ! iradtyp

              fds     = fds      + fd3d  
              fus     = fus      + fu3d  
              fdiv_sol= fdiv_sol + fdiv3d


#ifndef _XLF
              if(ldebug) then
                do k=1,nv1
                  if(any(isnan( fds     (k,:,:)))) print *,myid,'nan in radiation tendency fds     ',k,fds     (k,:,:)
                  if(any(isnan( fus     (k,:,:)))) print *,myid,'nan in radiation tendency fus     ',k,fus     (k,:,:)
                  if(any(isnan( fdiv3d  (k,:,:)))) print *,myid,'nan in radiation tendency fdiv3d  ',k,fdiv3d  (k,:,:)
                  if(any(isnan( fdiv_sol(k,:,:)))) print *,myid,'nan in radiation tendency fdiv_sol',k,fdiv_sol(k,:,:)
                  if(any(isnan([fds     (k,:,:),fus     (k,:,:),fdiv_sol(k,:,:),fdiv3d  (k,:,:)]))) call exit(1)
                enddo
              endif
#endif

            enddo !igpt
          enddo !iband

          !
          ! In this model, we used the solar spectral irradiance determined by
          ! Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained 
          ! in the spectral region 0.2 - 4.0 um., thus scale solar fluxes by
          ! fuq1
          !
          fuq1     = SolarConstant / totalpower
          fus      = fus *fuq1
          fds      = fds *fuq1
          fdiv_sol = fdiv_sol*fuq1
      end subroutine
      subroutine init_rad_3d
          use defs,only: p00,cp,cpr,nv,nv1
          use ckd,    only : init_ckd
          use cldwtr, only : init_cldwtr, init_cldice, init_cldgrp

          if(.not.linit) then
            is=3; ie=nxp-2
            js=3; je=nyp-2
            p0(nzp) = (p00*(pi0(nzp)/cp)**cpr) / 100.
            p0(nzp-1) = (p00*(pi0(nzp-1)/cp)**cpr) / 100.
            call setup(background,nzp,npts,nv1,nv,p0,pi0)
            linit = .True.
            if (allocated(pre))   pre(:) = 0.
            if (allocated(pde))   pde(:) = 0.
            if (allocated(piwc)) piwc(:) = 0.
            if (allocated(prwc)) prwc(:) = 0.
            if (allocated(plwc)) plwc(:) = 0.
            if (allocated(pgwc)) pgwc(:) = 0.

            !Initializations from file: rad_d4stream
            call init_ckd
            call init_cldwtr
            call init_cldice
            call init_cldgrp

          end if

          ! reset 3d heating rate output from grid module
          a_rhl = 0
          a_rhs = 0

          ! allocate work arrays
          allocate( fus     (nv1, is:ie, js:je) ) ; fus=0
          allocate( fds     (nv1, is:ie, js:je) ) ; fds=0
          allocate( fdiv_sol(nv1, is:ie, js:je) ) ; fdiv_sol=0

          allocate( fuir    (nv1, is:ie, js:je) ) ; fuir=0
          allocate( fdir    (nv1, is:ie, js:je) ) ; fdir=0
          allocate( fdiv_th (nv1, is:ie, js:je) ) ; fdiv_th=0


          allocate( hr_factor(nzp, is:ie, js:je) ); hr_factor=0
          allocate( dz(nv,is:ie, js:je) )

          ! initialize surface albedo, emissivity and skin temperature.
          ee = 1.0

          ! determine the solar geometery, as measured by u0, the cosine of the
          ! solar zenith angle
          if (.not. fixed_sun) u0 = zenith(alat,time)

          ! TODO we do calculate profiles here and then again in thermal AND solar
          ! TODO -- this is not elegent :( however not that straightforward to circumvent....
          ! TODO need it here in case radiation solution is loaded from tenstream
          ! TODO and need it at optprop calls because lwc is just 1D vector ... 
          ! TODO  should refactor these parts also into 3D arrays but then loose historic coupling to 1D routines.
          do j=js,je
            do i=is,ie
              if (present(ice).and.present(grp)) then
                call setup_rad_atmosphere(CCN, dn0, pi0, pi1,    &
                    pip(:,i,j), th(:,i,j), rv(:,i,j), rc(:,i,j),     &
                    hr_factor(:,i,j), dz(:,i,j),         &
                    rr=rr(:,i,j),ice=ice(:,i,j),nice=nice(:,i,j),grp=grp(:,i,j))
              else
                call setup_rad_atmosphere(CCN, dn0, pi0, pi1,    &
                    pip(:,i,j), th(:,i,j), rv(:,i,j), rc(:,i,j),     &
                    hr_factor(:,i,j), dz(:,i,j)         )
              end if
            enddo
          enddo

      end subroutine
  end subroutine rad_3d
  subroutine setup_rad_atmosphere(CCN, dn0, pi0, pi1, pip, th, rv, rc, hr_fac,dz,    rr,ice,nice,grp)
      use grid, only: dzi_m
      use defs, only: cp,cpr,nv,nv1,p00,pi,roice,rowt
      use fuliou, only: thicks
      real, intent(in) :: CCN
      real, dimension (nzp), intent (in)           :: dn0, pi0, pi1
      real, dimension (nzp), intent (in)           :: pip, th, rv, rc

      real,intent(out) :: hr_fac(nzp),dz(nv)
      !thos global vars are in radiation module...  ::  pp, pt, ph, po, pre, pde, plwc, piwc, prwc, pgwc

      real, optional, dimension (nzp), intent (in) :: rr,ice,nice,grp
      real    :: exner(nzp), pres(nzp) 

      real :: prw, pri
      integer :: kk,k
      prw = (4./3.)*pi*rowt
      pri = (3.*sqrt(3.)/8.)*roice

      do k=1,nzp
        exner(k)= (pi0(k)+pi1(k)+pip(k))/cp
        pres(k) = p00 * (exner(k))**cpr
      end do
      pp(nv1) = 0.5*(pres(1)+pres(2)) / 100.
      do k=2,nzp
        kk = nv-(k-2)
        pt(kk) = th(k)*exner(k)
        ph(kk) = max(0.,rv(k))
        plwc(kk) = 1000.*dn0(k)*max(0.,rc(k))
        pre(kk)  = rad_eff_radius*1.e6*(plwc(kk)/(1000.*prw*CCN*dn0(k)))**(1./3.)
        pre(kk)=min(max(pre(kk),4.18),31.23)
        if (plwc(kk).le.0.) pre(kk) = 0.
        if (present(rr)) then
          prwc(kk) = 1000.*dn0(k)*rr(k)
        else
          prwc(kk) = 0.
        end if
        if (present(ice)) then
          piwc(kk) = 1000.*dn0(k)*ice(k)
          if (nice(k).gt.0.0) then
            pde(kk)  = 1.e6*(piwc(kk)/(1000.*pri*nice(k)*dn0(k)))**(1./3.)
            pde(kk)=min(max(pde(kk),20.),180.)
          else
            pde(kk)  = 0.0
          endif
        else
          piwc(kk) = 0.
          pde(kk) = 0.0
        end if
        if (present(grp)) then
          pgwc(kk) = 1000.*dn0(k)*grp(k)
        else
          pgwc(kk) = 0.
        end if
        if (k < nzp) pp(kk) = 0.5*(pres(k)+pres(k+1)) / 100.
      end do
      pp(nv-nzp+2) = pres(nzp)/100. - 0.5*(pres(nzp-1)-pres(nzp)) / 100.

!      do k=1,nv
!        if(pp(k).gt.pp(k+1)) then
!          print *,'Error in Pressure computation, upper layer has more pressure than lower',pp(k).gt.pp(k+1)
!          print *,'pp',pp
!          print *,'pres',pres
!          print *,'exner',exner
!          stop 'pressure error'
!        endif
!      enddo

      call thicks(pp, pt, ph, dz) 

      do k=2,nzp
        hr_fac(k)  = dzi_m(k)/(cp*dn0(k)*exner(k))
      enddo
  end subroutine

  !calc of optprop extracted from rad_ir:
  subroutine optprop_rad_ir( ibandloop, ibandg, pp, pt, ph, po, tau, w, pf, dz, plwc, pre, piwc, pde, pgwc)
      use cldwtr, only : cloud_water, cloud_ice, cloud_grp
      use fuliou, only : gascon, combineopticalproperties,gases
      use ckd, only: ir_bands,solar_bands,center
      integer, intent(in) :: ibandloop,ibandg

      real, intent (in)  :: pp (:) ! pressure at interfaces

      real, dimension(:), intent (in)  :: &
          pt,   & ! temperature [K] at mid points
          ph,   & ! humidity mixing ratio in kg/kg
          po      ! ozone mixing ratio

      real, dimension(:)  , intent(out) :: tau,w  ! dim: (nv)
      real, dimension(:,:), intent(out) :: pf     ! dim: (nv,4)
      real, dimension(:)  , intent(in)  :: dz     ! dim: (nv)

      real, optional, dimension(:), intent (in)  :: & ! dim: (nv)
          plwc, & ! cloud liquid water content [g/m^3]
          pre,  & ! effective radius of cloud droplets [microns]
          piwc, & ! cloud ice water content [g/m^3]
          pde,  & ! effective diameter of ice particles [microns]
          pgwc    ! graupel water content

      ! ----------------------------------------

      real, dimension (size(tau))   :: tw,ww,tg, tauNoGas, wNoGas
      real, dimension (size(tau))   :: ti,wi
      real, dimension (size(tau))   :: tgr,wgr
      real, dimension (size(tau),4) :: www, pfNoGas
      real, dimension (size(tau),4) :: wwi
      real, dimension (size(tau),4) :: wwgr

      integer :: ib, ig,k
      ! ----------------------------------------
      tau=0
      w=0
      pf=0

      ib = ibandloop
      ig = ibandg

      ! Water vapor continuum optical depth
      !
      call gascon ( center(ir_bands(ib)), pp, pt, ph, TauNoGas )
      wNoGas = 0.; pfNoGas  = 0.
      if (present(plwc)) then
        call cloud_water(ib + size(solar_bands), pre, plwc, dz, tw, ww, www)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw, ww, www)
      end if
      if (present(piwc)) then
        call cloud_ice(ib + size(solar_bands), pde, piwc, dz, ti, wi, wwi)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, ti, wi, wwi)
      end if
      if (present(pgwc)) then
        call cloud_grp(ib + size(solar_bands), pgwc, dz, tgr, wgr, wwgr)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgr, wgr, wwgr)
      end if

      tau = TauNoGas; w = wNoGas; pf = pfNoGas
      call gases (ir_bands(ib), ig, pp, pt, ph, po, tg )
      call combineOpticalProperties(tau, w, pf, tg)

#ifndef _XLF
      if(ldebug) then
        if(any(isnan([tau, w, pf]))) then
          do k=1,size(pt)
            print *,'DEBUG',k,'pp',pp(k),'pt',pt(k),'ph',ph(k),'dz',dz(k),'opt',tau(k), w(k), pf(k,:)
          enddo
          call exit(-1)
        endif
      endif
#endif
  end subroutine
  !calc of optprop extracted from rad_vis:
  subroutine optprop_rad_vis( ibandloop, ibandg, pp, pt, ph, po, tau, w, pf, dz, plwc, pre, piwc, pde, pgwc)
      use cldwtr, only : cloud_water, cloud_ice, cloud_grp
      use fuliou, only : rayle, gascon, combineopticalproperties,gases
      use ckd, only: solar_bands,power,center
      integer, intent(in) :: ibandloop,ibandg

      real, intent (in)  :: pp (:) ! pressure at interfaces

      real, dimension(:), intent (in)  :: &
          pt,   & ! temperature [K] at mid points
          ph,   & ! humidity mixing ratio in kg/kg
          po      ! ozone mixing ratio

      real, dimension(:)  , intent(out) :: tau,w  ! dim: (nv)
      real, dimension(:,:), intent(out) :: pf     ! dim: (nv,4)
      real, dimension(:)  , intent(in)  :: dz     ! dim: (nv)

      real, optional, dimension(:), intent (in)  :: & ! dim: (nv)
          plwc, & ! cloud liquid water content [g/m^3]
          pre,  & ! effective radius of cloud droplets [microns]
          piwc, & ! cloud ice water content [g/m^3]
          pde,  & ! effective diameter of ice particles [microns]
          pgwc    ! graupel water content

      ! ----------------------------------------

      real, dimension (size(tau))   :: tw,ww,tg,tgm, tauNoGas, wNoGas
      real, dimension (size(tau))   :: ti,wi
      real, dimension (size(tau))   :: tgr,wgr
      real, dimension (size(tau),4) :: www, pfNoGas
      real, dimension (size(tau),4) :: wwi
      real, dimension (size(tau),4) :: wwgr

      integer :: ib, ig, k
      ! ----------------------------------------
      tau=0
      w  =0
      pf =0

      ib = ibandloop
      ig = ibandg

      !
      ! Rayleigh scattering
      !
      call rayle ( ib, u0, power(solar_bands(ib)), pp, pt, dz, tauNoGas, &
          wNoGas, pfNoGas)
      !
      ! Water vapor continuum
      !
      call gascon ( center(solar_bands(ib)), pp, pt, ph, tgm )
      if(any(tgm > 0.)) &
          call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgm)
      !
      ! Cloud water
      !
      if (present(plwc)) then
        call cloud_water(ib, pre, plwc, dz, tw, ww, www)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw,ww,www)
      end if
      if (present(piwc)) then
        call cloud_ice(ib, pde, piwc, dz, ti, wi, wwi)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, ti,wi,wwi)
      end if 
      if (present(pgwc)) then
        call cloud_grp(ib,pgwc, dz, tgr, wgr, wwgr)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgr, wgr,wwgr)
      end if 

      tau = tauNoGas; w = wNoGas; pf = pfNoGas
      call gases (solar_bands(ib), ig, pp, pt, ph, po, tg )
      call combineOpticalProperties(tau, w, pf, tg)
      if(ldebug) then
        if(any(isnan([tau, w, pf]))) then
          do k=1,size(pt)
            print *,'DEBUG',k,'pp',pp(k),'pt',pt(k),'ph',ph(k),'dz',dz(k),'lwc',plwc(k),'opt',tau(k), w(k), pf(k,:)
          enddo
          call exit(-1)
        endif
      endif
  end subroutine
  function get_band_uid(lsolar,iband,ibandg)
      use ckd, only: solar_bands,ir_bands,kg
      logical,intent(in) :: lsolar
      integer,intent(in) :: iband,ibandg
      integer :: get_band_uid

      integer :: ib,ig

      get_band_uid=-100
      if(lsolar) then

        get_band_uid=0
        outer_s: do ib = 1,size(solar_bands)
          do ig = 1,kg(solar_bands(ib))
            get_band_uid = get_band_uid+1
            if( ib.eq.iband .and. ig.eq.ibandg ) exit outer_s
          enddo
        enddo outer_s

      else

        get_band_uid=500
        outer_ir: do ib = 1,size(ir_bands)
          do ig = 1,kg(ir_bands(ib))
            get_band_uid = get_band_uid+1
            if( ib.eq.iband .and. ig.eq.ibandg ) exit outer_ir
          enddo
        enddo outer_ir

      endif

      return
  end function

#ifdef HAVE_TENSTREAM
  subroutine tenstream_wrapper(lsolar, in_nxp,in_nyp,in_nv,in_dx,in_dy,dz, in_phi0,in_u0, in_albedo,in_incSolar, tau, w0, pf, bf, fdn,fup,fdiv)
      use mpi, only: mpi_comm_world

      logical                ,intent(in) :: lsolar
      integer                ,intent(in) :: in_nxp,in_nyp,in_nv
      real                   ,intent(in) :: in_dx,in_dy,in_phi0,in_u0,in_albedo,in_incSolar
      real,dimension(:,:,:)  ,intent(in) :: tau,w0,dz,bf ! have dimensions (nv,nxp-4,nyp-4); bf with (nv1)
      real,dimension(:,:,:,:),intent(in) :: pf

      real,dimension(:,:,:),intent(out) :: fdn,fup,fdiv

      integer(iintegers) :: nxp,nyp,nv
      real(ireals) :: dx,dy,phi0,u0,albedo
      real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,deltaz ! dim(in_nv, 3:in_nxp-2,3:in_nyp-2)
      real(ireals),allocatable,dimension(:,:,:) :: planck
      integer(iintegers) :: nxproc(size(nxpa)), nyproc(size(nypa))

      integer(iintegers) :: k
      real(ireals)       :: theta0,incSolar


      nxp=in_nxp; nyp=in_nyp; nv=in_nv
      dx=in_dx; dy=in_dy
      phi0=in_phi0; u0=in_u0; albedo=in_albedo
      nxproc=nxpa; nyproc=nypa
      if(ldebug.and.myid.eq.0) print *,'Tenstrwrapper lsol',lsolar,'nx/y',nxp,nyp,nv,'uid',solution_uid,solution_time,' shapes::',shape(dz),shape(fdn),shape(fup),shape(fdiv)

      if(lsolar .and. u0.gt.minSolarZenithCosForVis) then
        theta0=acos(u0)*180._ireals/3.141592653589793_ireals !rad2deg
        incSolar = in_incSolar
      else
        theta0=0
        incSolar=0
      endif

#ifdef TENSTREAM_SINGLE
        if(.not.allocated(deltaz) ) allocate(deltaz(nv, 3:nxp-2,3:nyp-2))
        if(.not.allocated(kabs  ) ) allocate(kabs  (nv, 3:nxp-2,3:nyp-2))
        if(.not.allocated(ksca  ) ) allocate(ksca  (nv, 3:nxp-2,3:nyp-2))
        if(.not.allocated(g     ) ) allocate(g     (nv, 3:nxp-2,3:nyp-2))
        deltaz = dz
        kabs   = max(epsilon(tau), tau * (one - w0) / dz )
        ksca   = max(epsilon(tau), tau *        w0  / dz )
        g      = min( one, pf (:,1,:,:)/3._ireals ) !todo: optprop used with ice microphysics gives g>1 -- we should fix the issue instead of constraining the value!

        if(.not. lsolar) then
          if(.not.allocated(planck) ) allocate(planck(nv+1, 3:nxp-2,3:nyp-2))
          planck = bf 
        endif

        if(ldebug) then
          if(any(isnan(kabs))) print *,myid,'tenstream_wrapper :: corrupt kabs',kabs,'::',tau,'::',w0,'::',deltaz
          if(any(isnan(ksca))) print *,myid,'tenstream_wrapper :: corrupt ksca',ksca,'::',tau,'::',w0,'::',deltaz
          if(any(isnan(g   ))) print *,myid,'tenstream_wrapper :: corrupt g   ',g   ,'::',pf (:,1,:,:)                                      
        endif

        call init_tenstream(MPI_COMM_WORLD, nv, nxp-4,nyp-4, dx,dy,phi0, theta0, albedo, nxproc=nxproc, nyproc=nyproc,  dz3d=deltaz)
        if(lsolar) then
          call set_optical_properties( kabs, ksca, g )
        else
          call set_optical_properties( kabs, ksca, g, planck)
        endif

        call solve_tenstream(incSolar,solution_uid,solution_time)

        call load_tenstream_solution(lsolar,dz,in_u0,solution_uid,fdn,fup,fdiv)

#else 

        call init_tenstream(MPI_COMM_WORLD, nv, nxp-4,nyp-4, dx,dy,phi0, theta0, albedo, nxproc=nxproc, nyproc=nyproc,  dz3d=dz)
        if(lsolar) then
          call set_optical_properties( max(epsilon(tau), tau * (one - w0) / dz ), max(epsilon(tau), tau *       w0  / dz ), min(one, pf (:,1,:,:)/3._ireals) )
        else
          call set_optical_properties( max(epsilon(tau), tau * (one - w0) / dz ), max(epsilon(tau), tau *       w0  / dz ), min(one, pf (:,1,:,:)/3._ireals), bf )
        endif

        call solve_tenstream(incSolar,solution_uid,solution_time)

        call load_tenstream_solution(lsolar,dz,in_u0,solution_uid,fdn,fup,fdiv)

#endif        

        if(ldebug) then
          if(any(isnan([fdn,fup,fdiv]))) then
            do k=1,nv+1
              print *,myid,'DEBUG',k,phi0, theta0,albedo
              print *,myid,'edn ::',fdn (k,:,:)
              print *,myid,'eup ::',fup (k,:,:)
              print *,myid,'div ::',fdiv(k,:,:)
            enddo
            call exit(-1)
          endif
          if(lsolar.and.any([fdn,fup].lt.-1._ireals)) then
            do k=1,nv+1
              print *,myid,'DEBUG value less than zero in solar rad',k,phi0, theta0,albedo
              print *,myid,'edn ::',k,minval(fdn (k,:,:))
              print *,myid,'eup ::',k,minval(fup (k,:,:))
              print *,myid,'div ::',k,minval(fdiv(k,:,:))
            enddo
            call exit(-1)
          endif
        endif

  end subroutine
  subroutine load_tenstream_solution(lsolar,dz,u0,uid,fdn,fup,fdiv)
      logical ,intent(in) :: lsolar
      real,dimension(:,:,:),intent(in) :: dz ! dimensions (nv,nxp-4,nyp-4)
      real,intent(in) :: u0
      integer(iintegers) :: uid
      real,dimension(:,:,:),intent(out) :: fdn,fup,fdiv
      real(ireals),dimension(:,:,:),allocatable :: edir,edn,eup
      real(ireals),dimension(:,:,:),allocatable :: abso

      integer :: i,j,k
      integer :: is,ie,js,je,ks,ke

      if(ldebug.and.myid.eq.0) print *,'load_tenstream_solution',shape(dz),shape(fdn),shape(fup),shape(fdiv)

      is = lbound(fdn,2); ie = ubound(fdn,2)
      js = lbound(fdn,3); je = ubound(fdn,3)
      ks = lbound(fdn,1); ke = ubound(fdn,1)

      !      print *,myid,'load_tenstream_solution inp',lsolar,u0
      !      print *,myid,'load_tenstream_solution dz ',dz       
      !      print *,myid,'load_tenstream_solution dim1',is,ie
      !      print *,myid,'load_tenstream_solution dim1',js,je
      !      print *,myid,'load_tenstream_solution dim1',ks,ke

#ifdef TENSTREAM_SINGLE
        if(.not.allocated(edn ) ) allocate(edn (ks:ke  ,is:ie,js:je))
        if(.not.allocated(eup ) ) allocate(eup (ks:ke  ,is:ie,js:je))
        if(.not.allocated(abso) ) allocate(abso(ks:ke-1,is:ie,js:je))

        if(lsolar .and. u0.gt.minSolarZenithCosForVis) then
          if(.not.allocated(edir) ) allocate(edir(ks:ke  ,is:ie,js:je))
          call tenstream_get_result(edir,edn,eup,abso,uid)
          edn = edn+edir
        else
          call tenstream_get_result(edir,edn,eup,abso,uid)
        endif

        fdiv(ks:ke-1,:,:) = abso(:,:,:) *dz

        fdiv(ke,:,:) = edn(ke,:,:) - eup(ke,:,:)
        deallocate(abso)

        fdn = edn
        fup = eup
#else        
        if(.not.allocated(abso) ) allocate(abso(ks:ke-1,is:ie,js:je))
        if(lsolar .and. u0.gt.minSolarZenithCosForVis) then
          if(.not.allocated(edir) ) allocate(edir(ks:ke  ,is:ie,js:je))
          call tenstream_get_result(edir,fdn,fup,abso,uid)
          fdn = fdn+edir
        else
          call tenstream_get_result(edir,fdn,fup,abso,uid)
        endif

        fdiv(ks:ke-1,:,:) = abso(:,:,:) *dz

        fdiv(ke,:,:) = fdn(ke,:,:) - fup(ke,:,:)
#endif

!      if( lsolar .and. any(fdiv.lt.0.) ) then
!        print *,'Found values smaller than 0 in divergence:',minval(fdiv)
!        call exit()
!      endif
!      print *,myid,'load_tenstream_solution edn ::',edn (is,js,:)
!      print *,myid,'load_tenstream_solution eup ::',eup (is,js,:)
!      print *,myid,'load_tenstream_solution  dz ::',dz  (is,js,:)
!      print *,myid,'load_tenstream_solution fdn ::',fdn (:,is,js)
!      print *,myid,'load_tenstream_solution fup ::',fup (:,is,js)
  end subroutine
#endif

  ! Caro's stuff
  subroutine nca_wrapper(nxp, nyp, nv, dx, dy, dz, tau, w0, bf, fdir, fuir, fdiv)
    
    use m_nca, only: nca

    use mpi_interface, only : cyclics, cyclicc, cyclicsp1

    integer :: req(16)
    integer,intent(in) :: nxp,nyp,nv
    real,intent(in) :: dx,dy
    real,dimension(:,:,:),intent(in) :: tau,w0,dz     ! have dimensions (nv,nxp,nyp)
    real,dimension(:,:,:),intent(in)  :: fdir,fuir,bf ! have dimensions (nv1,nxp,nyp)
    real,dimension(:,:,:),intent(out) :: fdiv         ! have dimensions (nv,1:nxp-4,1:nyp-4)
    real,dimension(nzp,nxp,nyp)   :: kabs             ! (nv,nxp,nyp)
    real,dimension(nzp,nxp,nyp)   :: heat             ! (nv,nxp,nyp)

    real :: fdir_ghosted(nzp+1,nxp,nyp)
    real :: fuir_ghosted(nzp+1,nxp,nyp)
    real :: w0_ghosted(nzp,nxp,nyp)
    real :: dz_ghosted(nzp,nxp,nyp)
    real :: bf_ghosted(nzp+1,nxp,nyp)
    real :: tau_ghosted(nzp,nxp,nyp)

 !   fdir_ghosted=-1!DEBUG
 !   fuir_ghosted=-1!DEBUG
 !   w0_ghosted=-1!DEBUG
 !   bf_ghosted=-1!DEBUG
 !   tau_ghosted=-1!DEBUG
 !   dz_ghosted=-1!DEBUG

    fdir_ghosted(:,3:nxp-2,3:nyp-2) = fdir(nv+1-nzp:nv+1,:,:)
    fuir_ghosted(:,3:nxp-2,3:nyp-2) = fuir(nv+1-nzp:nv+1,:,:)
    w0_ghosted(:,3:nxp-2,3:nyp-2) = w0(nv+1-nzp:nv,:,:)
    tau_ghosted(:,3:nxp-2,3:nyp-2) = tau(nv+1-nzp:nv,:,:)
    bf_ghosted(:,3:nxp-2,3:nyp-2) = bf(nv+1-nzp:nv+1,:,:)
    dz_ghosted(:,3:nxp-2,3:nyp-2) = dz(nv+1-nzp:nv,:,:)

    call cyclicsp1(nzp+1,nxp,nyp,fdir_ghosted,req)
    call cyclicsp1(nzp+1,nxp,nyp,fuir_ghosted,req)
    call cyclicc(nzp+1,nxp,nyp,fdir_ghosted,req)
    call cyclicc(nzp+1,nxp,nyp,fuir_ghosted,req)
    call cyclicsp1(nzp+1,nxp,nyp,bf_ghosted,req)
    call cyclicc(nzp+1,nxp,nyp,bf_ghosted,req)

    call cyclics(nzp,nxp,nyp,w0_ghosted,req)
    call cyclicc(nzp,nxp,nyp,w0_ghosted,req)
    call cyclics(nzp,nxp,nyp,tau_ghosted,req)
    call cyclicc(nzp,nxp,nyp,tau_ghosted,req)
    call cyclics(nzp,nxp,nyp,dz_ghosted,req)
    call cyclicc(nzp,nxp,nyp,dz_ghosted,req)

    ! Get kabs from optical thicknes (tau)
    kabs  = tau_ghosted*(1.-w0_ghosted) / dz_ghosted
    
    heat=0
    call nca(nzp, dz_ghosted, bf_ghosted, kabs, heat, dx, dy, fdir_ghosted, fuir_ghosted) 
    
    fdiv=0
  
    fdiv(nv+1-nzp:nv,:,:) = heat(1:nzp,3:nxp-2,3:nyp-2)
    
#ifndef _XLF
    if(any(isnan(fdiv))) print *, 'fdiv 2 shows nan', fdiv
#endif
    
  end subroutine nca_wrapper

end module
