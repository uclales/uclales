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
! \todo CMake as in microhh
! \todo Separate Surface Namelist?
! \todo Separate Soil from 3D fields?
! \todo Top level radiation hopping
! \todo reduce nr of warnings when compiling with gfortran/NAG
! \todo Nicer version of accumelated precip
program ucla_les

  implicit none

  real :: t1, t2

  call cpu_time(t1)
  call driver
  call cpu_time(t2)

  print "(/,' ',49('-')/,' ',A16,F10.1,' s')", '  Execution time: ', t2-t1
  stop ' ..... Normal termination'

contains

  !----------------------------------------------------------------------
  ! Subroutine Driver:  This is the main program driver.  It calls routines
  ! to read the model initialization file, and configure memory and pointes.
  ! It also calls the routines which initialize the model and timestep it.
  !
  subroutine driver

    use grid, only          : define_grid, define_vars, nxp, nyp, nzp, nxpart
    use init, only          : initialize
    use step, only          : stepper
    use mpi_interface, only : init_mpi, define_decomp,                    &
         init_alltoall_reorder, appl_finalize
    implicit none
    integer ierror

    call init_mpi
    call define_parm
    call define_decomp(nxp, nyp, nxpart)
    call define_grid
    call init_alltoall_reorder(nxp, nyp, nzp)
    call define_vars
    call initialize
    call stepper
    call appl_finalize(ierror)

    return
  end subroutine driver

  !
  ! ----------------------------------------------------------------------
  ! Subroutine Read_nl: Driver for reading model namelist
  !
  subroutine define_parm

    use util, only : fftinix,fftiniy
    use defs, only : SolarConstant
    use sgsm, only : csx, prndtl, clouddiff
    use advf, only : lmtr !,advs
    use advl, only : advm
    use srfc, only : zrough, ubmin, dthcon, drtcon, rh_srf, drag, lhomflx
    use step, only : timmax, timrsm, istpfl, corflg, outflg, frqanl, frqhis,          &
         frqcross , strtim, radfrq, cntlat,&
         case_name,lsvarflg, sst, div, wctime , tau, &                !irina
         lanom ! linda
    use modnetcdf, only : lsync, deflate_level
    use ncio, only : deflev => deflate_level
    use modcross, only : lcross, lxy,lxz,lyz,xcross,ycross,zcross, crossvars
    use forc, only : lstendflg, sfc_albedo !, w0, Qrate
    use grid, only : deltaz, deltay, deltax, nzp, nyp, nxp, nxpart,           &
         dtlong, dzrat,dzmax, th00, umean, vmean, naddsc, level,              &
         filprf, expnme, iradtyp, igrdtyp, nfpt, distim, lspongeinit, runtype,             &
         CCN, lwaterbudget, lcouvreux, prc_lev, isfctyp, sfc_albedo, lrad_ca, &
	 w0, Qrate, iradbel, wtgbel !RV
    use init, only : us, vs, ts, rts, ps, hs, ipsflg, itsflg,irsflg, iseed, hfilin,   &
         zrand,lhomrestart,mag_pert_q,mag_pert_t
    use stat, only : ssam_intvl, savg_intvl
    use mpi_interface, only : myid, appl_abort
    use radiation, only : u0, fixed_sun, rad_eff_radius, radMcICA
    use modnudge, only : lnudge,tnudgefac, qfloor, zfloor, znudgemin, znudgeplus, &
         lnudge_bound
    use modtimedep, only : ltimedep
    use mcrp, only : microseq,lrandommicro,timenuc,nin_set,cloud_type, lpartdrop
    use modparticles, only : lpartic, lpartsgs, lrandsurf, lpartstat, lpartdump, &
         lpartdumpui, lpartdumpth, lpartdumpmr, frqpartdump, ldropstart

    implicit none

    namelist /model/  &
         expnme    ,                         & ! experiment name
         nxpart    ,                         & ! whether partition in x direction?
         naddsc    ,                         & ! Number of additional scalars
         savg_intvl,                         & ! output statistics frequency
         ssam_intvl,                         & ! integral accumulate/ts print frequency
         corflg , cntlat ,                   & ! coriolis flag
         nfpt   , distim ,                   & ! rayleigh friction points, dissipation time
         lspongeinit     , & ! Sponge back to initial profile or bulk values
         level  , CCN    ,                   & ! Microphysical model Number of CCN per kg of air
         iseed  , zrand  ,                   & ! random seed
         mag_pert_q , mag_pert_t ,            & ! Magnitude of pertubations
         nxp    , nyp    , nzp   ,           & ! number of x, y, z points
         deltax , deltay , deltaz ,          & ! delta x, y, z (meters)
         dzrat  , dzmax  , igrdtyp,          & ! stretched grid parameters
         timmax , dtlong , istpfl , timrsm, wctime, & ! timestep control
         runtype, hfilin , filprf ,          & ! type of run (INITIAL or HISTORY)
         frqhis , frqanl, frqcross, outflg , & ! freq of history/anal writes, output flg
         lsync, lcross, lxy,lxz,lyz,xcross,ycross,zcross, crossvars,prc_lev,&
                  iradtyp, radfrq , strtim , sfc_albedo, & ! radiation type flag
         isfctyp, ubmin  , zrough ,          & ! surface parameterization type
         sst    , dthcon , drtcon ,          & ! SSTs, surface flx parameters
         csx    , prndtl ,                   & ! SGS model type, parameters
         ipsflg , itsflg , irsflg,           & ! sounding flags
         hs     , ps     , ts    ,           & ! sounding heights, pressure, temperature
         us     , vs     , rts   ,           & ! sounding E/W winds, water vapor
         umean  , vmean  , th00  ,           & ! gallilean E/W wind, basic state
         lanom  ,                            & ! LINDA, start with anomalies
         case_name, lmtr,                    & ! irina:name of the case, i.e. astex, rico, etc
         advm,                               & ! Advection scheme scalars, momentum
         lsvarflg,                           & ! irina:flag for time bvarying large scale forcing
         lstendflg,                          & ! irina:flag for time large scale advective tendencies
         div,                                & ! irina: divergence
         lnudge, tnudgefac, ltimedep, qfloor, zfloor,znudgemin, znudgeplus,  &             !thijs: Nudging
         lnudge_bound, &               ! LINDA, relaxation boundaries
         rh_srf, drag, &
         SolarConstant,u0,fixed_sun, rad_eff_radius, radMcICA, & ! SolarConstant (In case of prescribed TOA radiation
         lrandommicro, microseq,timenuc ,nin_set,cloud_type, &  !thijs: sequence of variables for microphysics
         lwaterbudget, &                 ! axel: flag for liquid water budget diagnostics (only level=3)
         lcouvreux , tau , &                    ! The Couvreux 'radioactive' scalar
         lrad_ca, &                        ! Clear air radiation statistics
         deflate_level , lhomflx,lhomrestart, &                         !Compression of the crosssections
         clouddiff, &
         lpartic,lpartsgs,lrandsurf,lpartstat,lpartdump, &           ! Particles
         lpartdumpui,lpartdumpth,lpartdumpmr,frqpartdump,&           ! Particles
         lpartdrop, ldropstart, &                                    ! Particles
	 w0, Qrate, iradbel, wtgbel                                    !RV: parameters for large-scale forcings in bellon_rad


    deflev = deflate_level
    ps       = 0.
    ts       = th00
    !
    ! these are for initializing the temp variables used in ffts in x and y
    ! directions.
    !
      fftinix=1
      fftiniy=1
    !
    ! read namelist from specified file
    !
    open  (1,status='old',file='NAMELIST')
    read  (1, nml=model)
    close (1)
!     write (0,model)
    !
    ! write file variable control to standard output
    !
    if (myid == 0) then
       if (runtype == 'HISTORY') then
          write (*,601) expnme, hfilin, timmax
       else
          write (*,600) expnme, timmax
       end if
       if (outflg) write (*,602) filprf, frqhis, frqanl
       !
       ! do some cursory error checking in namelist variables
       !

       if (min(nxp,nyp) < 5) then
          if (myid == 0) print *, '  ABORTING: min(nxp,nyp) must be > 4.'
          call appl_abort(0)
       endif

       if (nzp < 3 ) then
          if (myid == 0) print *, '  ABORTING: nzp must be > 2 '
          call appl_abort(0)
       endif

       if (cntlat < -90. .or. cntlat > 90.) then
          if (myid == 0) print *, '  ABORTING: central latitude out of bounds.'
          call appl_abort(0)
       endif

       if(isfctyp == 4 .and. level == 0) then 
          if (myid == 0) print *, ' Fixed buoyancy flux (isfctyp=4, Stevens, 2007, JAS) without moisture (level=0) not supported'
          call appl_abort(0)
       endif

       if(isfctyp == 5 .and. level == 0) then 
          if (myid == 0) print *, ' Land surface scheme (isfctyp=5) without moisture (level=0) not supported'
          call appl_abort(0)
       endif


    end if

600 format(//' ',49('-')/,' ',/,'  Initial Experiment: ',A50 &
         /,'  Final Time:         ',F7.1,' s'              )
601 format(//' ',49('-')/,' ',/,'  Restart Experiment: ',A50 &
         /,'  Restart File: ',A30,                           &
         /,'  Final Time: ',F10.1,' s'              )
602 format('  Output File Stem:   ',A50                      &
         /,'  History Frequency:  ',F7.1,                    &
         /,'  Analysis Frequency: ',F7.1)

    return
  end subroutine define_parm

end program ucla_les

