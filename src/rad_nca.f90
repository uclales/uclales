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
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2015, **
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------
module m_nca
  implicit none

contains
  subroutine nca (nlyr, dz, B, kabs_3d, hr, dx, dy, Edn, Eup)

    integer, intent(in) :: nlyr
    real, intent(in) :: dz(:,:,:), kabs_3d(:,:,:), dx, dy   ! dimensions including ghost values (nlyr,Nx,Ny)
    real, intent(in) :: Edn(:,:,:), Eup(:,:,:), B(:,:,:)    ! dimensions including ghost values (nlyr+1,Nx,Ny)
    real, intent(out) :: hr(:,:,:)                          ! dimensions including ghost values (nlyr,Nx,Ny)

    ! ############## Definition of variables ##################

    integer        :: ilyr          ! counter z levels
    integer        :: ixx           ! counter x grid boxes
    integer        :: iyy           ! counter y grid boxes

    real           :: Absup         ! Upwelling, Absorption, bottom
    real           :: Absdn         ! Downwelling Absorption, top
    real           :: Absup_s       ! Upwelling, Absorption, side
    real           :: Absdn_s       ! Downwelling Absorption, side  
    real           :: Emup          ! Upwelling, Emission, bottom
    real           :: Emdn          ! Downwelling Emission, top
    real           :: Emup_s        ! Upwelling, Emission, side
    real           :: Emdn_s        ! Downwelling Emission, side  

    real           :: HR_T          ! Heating Rate Top Contribution  
    real           :: HR_s          ! Heating Rate Side Contribution 

    real           :: mu            ! zenith angle

    real           :: ax,ay         ! integration boarder for side contributions
    real           :: bx,by         ! integration boarder for side contributions 
    real           :: cx,cy         ! integration boarder for side contributions
    real           :: az1,az        ! integration boarder for side contributions
    real           :: bz1,bz        ! integration boarder for side contributions 
    real           :: cz1,cz        ! integration boarder for side contributions 


    real           :: l      ! Averaged Transmitted Flux, Top Contribution 
    real           :: Trans  ! Transmitted Flux, Top Contribution

    real           :: L_up_3d(lbound(Eup,1):ubound(Eup,1),lbound(Eup,2):ubound(Eup,2),lbound(Eup,3):ubound(Eup,3))
    real           :: L_dn_3d(lbound(Edn,1):ubound(Edn,1),lbound(Edn,2):ubound(Edn,2),lbound(Edn,3):ubound(Edn,3))
    real           :: B2   ! Planck General
    real           :: B2_1 ! Planck Top Contribution

    ! ## Fitting Parameters
    real           :: afit1 
    real           :: bfit1 
    real           :: bfit2 
    real           :: afit2 
    real           :: cfit2 
    real           :: factor2 
    real           :: factor1 

    real,parameter :: eps = 0.0001
    real,parameter :: pi=3.141592653589793
 
    integer :: is,ie,js,je

    is = lbound(hr,2)+2
    ie = ubound(hr,2)-2
    js = lbound(hr,3)+2
    je = ubound(hr,3)-2

    hr = -9999999
    stop 'NCA not freely available.... please consult Carolin Klinger for an implementation.'
  end subroutine


end module m_nca
