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

  private
  public :: nca

contains
  subroutine nca (nlyr, Nx, Ny, dz, B, kabs_3d, hr_nca_3d_tmp, dx, dy, Ldn, Lup)

    integer, intent(in) :: Nx, Ny, nlyr
    real, intent(in) :: dz(:,:,:), B(:,:,:), kabs_3d(:,:,:), dx, dy
    real, intent(inout) :: Ldn(:,:,:), Lup(:,:,:)
    real, intent(out) :: hr_nca_3d_tmp(:,:,:)

    ! ############## Definition of variables ##################

    integer        :: ilyr = 0         ! counter z levels
    integer        :: ixx = 0          ! counter x grid boxes
    integer        :: iyy = 0          ! counter y grid boxes
    integer        :: iface = 0        ! counter side faces/directions of grid box

    real           :: Absup = 0        ! Upwelling, Absorption, lower/upper face
    real           :: Absdn = 0        ! Downwelling Absorption, lower/upper face
    real           :: Absup1 = 0       ! Upwelling, Absorption, face 1
    real           :: Absup2 = 0       ! Upwelling, Absorption, face 2
    real           :: Absup3 = 0       ! Upwelling, Absorption, face 3
    real           :: Absup4 = 0       ! Upwelling, Absorption, face 4
    real           :: Absdn1 = 0       ! Downwelling Absorption, face 1  
    real           :: Absdn2 = 0       ! Downwelling Absorption, face 2
    real           :: Absdn3 = 0       ! Downwelling Absorption, face 3
    real           :: Absdn4 = 0       ! Downwelling Absorption, face 4

    real           :: Emup = 0         ! Upwelling, Emission, lower/upper face
    real           :: Emdn = 0         ! Downwelling Emission, lower/upper face
    real           :: Emup1 = 0        ! Upwelling, Emission, face 1
    real           :: Emup2 = 0        ! Upwelling, Emission, face 2
    real           :: Emup3 = 0        ! Upwelling, Emission, face 3
    real           :: Emup4 = 0        ! Upwelling, Emission, face 4
    real           :: Emdn1 = 0        ! Downwelling Emission, face 1  
    real           :: Emdn2 = 0        ! Downwelling Emission, face 2  
    real           :: Emdn3 = 0        ! Downwelling Emission, face 3  
    real           :: Emdn4 = 0        ! Downwelling Emission, face 4  
   
    real           :: HR_up = 0        ! Downwelling Emission, face 1  
    real           :: HR_dn = 0        ! Downwelling Emission, face 2  
    real           :: HR_up_s = 0      ! Downwelling Emission, face 3  
    real           :: HR_dn_s = 0      ! Downwelling Emission, face 4  

    real           :: mu = 0           ! zenith angle

    real           :: areaweight = 0   ! area weight for non-cubic grid boxes

    real           :: ax = 0           ! integration boarder for side contributions
    real           :: bx = 0           ! integration boarder for side contributions 
    real           :: cx = 0           ! integration boarder for side contributions

    real           :: az = 0           ! integration boarder for side contributions
    real           :: bz = 0           ! integration boarder for side contributions 
    real           :: cz = 0           ! integration boarder for side contributions 

    real           :: pi=3.141592653589793
 
    real           :: factor = 0
    real           :: l = 0.0
    real           :: Trans = 0.0

    real           :: L_up_3d(nlyr,Nx,Ny)
    real           :: L_dn_3d(nlyr,Nx,Ny)
    real           :: B1 = 0
    real           :: B2 = 0
    real           :: B1_1 = 0
    real           :: B2_1 = 0

    hr_nca_3d_tmp=-999999
   
  end subroutine nca


end module m_nca
