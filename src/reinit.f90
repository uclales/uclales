module reinit

use util, only : phiset
use grid, only : write_anal

implicit none
private ! By keeping 'private' empty all subroutines and variables will be
        ! private to the module by default.
public reinitialize
integer, parameter :: VERTICAL_METHOD = 0
integer, parameter :: SUSSMAN_METHOD = 1
integer, parameter :: RUSSO_METHOD = 2
integer, parameter :: nt = 20
real, parameter :: eps = tiny(1.0)
real :: sign_eps2 

contains
!
!===============================================================================
! reinitialize: Driver routine for level set reinitialization.
!===============================================================================
subroutine reinitialize(phi, z, nxp, nyp, nzp, dxi, dyi, dzi)
integer, intent(in) :: nzp, nxp, nyp
real, intent(in)    :: dxi, dyi, dzi(nzp), z(nzp)
real, intent(inout) :: phi(nzp, nxp, nyp)
integer :: method
real :: dtau

method = RUSSO_METHOD
dtau = 0.75*min(min(1./dxi, 1./dyi), 1./maxval(dzi))
sign_eps2 = (0.5 / min(min(dxi, dyi), minval(dzi)))**2

select case(method)

case(VERTICAL_METHOD)
    call vertical(phi, nxp, nyp, nzp, z)

case(SUSSMAN_METHOD)
    stop "Reinitialization method not implemented."

case(RUSSO_METHOD)
    call russo(phi, nxp, nyp, nzp, dxi, dyi, dzi, dtau)

case default
    stop "Reinitialization method not implemented."

end select

end subroutine reinitialize
!
!===============================================================================
! vertival: Directly reinitializes the level set into a vertical signed-distance
!           function.
!===============================================================================
subroutine vertical(phi, nxp, nyp, nzp, z)
integer, intent(in) :: nzp, nxp, nyp
real, intent(in)    :: z(nzp)
real, intent(inout) :: phi(nzp, nxp, nyp)
integer :: k, i, j
real    :: ht

do j = 1, nyp
do i = 1, nxp
    !
    ! find inversion hight starting from the top
    !
    ht = z(nzp)
    k = nzp
    do while (k > 1)
        k = k - 1
        if (phi(k+1,i,j)*phi(k,i,j)<=0.) then
            ht = z(k) - phi(k,i,j) * (z(k+1)-z(k)) / (phi(k+1,i,j)-phi(k,i,j))
            k = 1
        end if
    end do
    !
    ! 1D reinitialization
    !
    do k = 1, nzp
      phi(k,i,j) = z(k) - ht
    end do
end do
end do

end subroutine vertical
!
!===============================================================================
! sussman: Reinitialize the level set using Sussman, Smereka, and Osher's (1994)
!          method.
!===============================================================================
subroutine sussman()

end subroutine sussman
!
!===============================================================================
! russo: Reinitialize the level set using Russo and Smereka's (2000) 'subcell
!        fix' with Sussman, Smereka and Osher's (1994) original method.
!===============================================================================
subroutine russo(phi, nxp, nyp, nzp, dxi, dyi, dzi, dtau)
integer, intent(in) :: nzp, nxp, nyp
real, intent(in)    :: dtau, dxi, dyi, dzi(nzp)
real, intent(inout) :: phi(nzp, nxp, nyp)
real, dimension(nzp, nxp, nyp) :: phin, phi0, d, G, phix, phiy, phiz
integer :: n, k, i, j, mask(nzp, nxp, nyp)

phi0(:,:,:) = phi(:,:,:)
call dRusso(d, phi0, dxi, dyi, dzi, nxp, nyp, nzp)
call maskRusso(mask, phi0, nxp, nyp, nzp)

do n = 1, nt
    phin(:,:,:) = phi(:,:,:)

    call gradientBackwards(phix, phiy, phiz, phi, nxp, nyp, nzp, dxi, dyi, dzi)
    call godunovHamiltonian(G, phi0, phix, phiy, phiz, nxp, nyp, nzp)

    do j = 3, nyp-2
    do i = 3, nxp-2
    do k = 2, nzp-1
        phi(k,i,j) = phin(k,i,j) - dtau * ( & 
              (1.0 - mask(k,i,j)) * smoothSign(phi0(k,i,j)) * (G(k,i,j) - 1.0) &
            + mask(k,i,j) * dxi   * (sign(1.0, phi0(k,i,j)) * abs(phin(k,i,j)) &
                                                                  - d(k,i,j) ) &
        )
    end do
    end do
    end do

    call phiset(nzp, nxp, nyp, phi, dzi)
end do

end subroutine russo
!
!===============================================================================
! dRusso : Returns an approximation for signed-distance function in Gamma cells.
!===============================================================================
subroutine dRusso(d, phi0, dxi, dyi, dzi, nxp, nyp, nzp)
real, intent(inout) :: d(nzp, nxp, nyp)
real, intent(in) :: dxi, dyi, dzi(nzp), phi0(nzp, nxp, nyp)
integer, intent(in) :: nxp, nyp, nzp
real :: phix, phiy, phiz
integer :: k, i, j
logical :: cutxp, cutxm, cutyp, cutym, cutzp, cutzm

d(:,:,:) = 0.0

do j = 3, nyp-2
do i = 3, nxp-2
do k = 2, nzp-1
    phix = 0.0
    phiy = 0.0
    phiz = 0.0

    cutxp = phi0(k,i,j)*phi0(k,i+1,j) < 0.0
    cutxm = phi0(k,i,j)*phi0(k,i-1,j) < 0.0
    cutyp = phi0(k,i,j)*phi0(k,i,j+1) < 0.0
    cutym = phi0(k,i,j)*phi0(k,i,j-1) < 0.0
    cutzp = phi0(k,i,j)*phi0(k+1,i,j) < 0.0
    cutzm = phi0(k,i,j)*phi0(k-1,i,j) < 0.0

    if (cutxp) then
        phix = ( phi0(k,i+1,j)  - phi0(k,i,j) ) * dxi
    end if
    if (cutxm) then
        phix = ( phi0(k,i,j)  - phi0(k,i-1,j) ) * dxi
    end if
    if (cutyp) then
        phiy = ( phi0(k,i,j+1)  - phi0(k,i,j) ) * dyi
    end if
    if (cutym) then
        phiy = ( phi0(k,i,j)  - phi0(k,i,j-1) ) * dyi
    end if
    if (cutzp) then
        phiz = ( phi0(k+1,i,j)  - phi0(k,i,j) ) * dzi(k)
    end if
    if (cutzm) then
        phiz = ( phi0(k,i,j)  - phi0(k-1,i,j) ) * dzi(k)
    end if

    if (cutxp.or.cutxm.or.cutyp.or.cutym.or.cutzp.or.cutzm) then
        d(k,i,j) = phi0(k,i,j) / (phix**2 + phiy**2 + phiz**2 + eps)**0.5
    end if
end do
end do
end do

end subroutine dRusso
!
!===============================================================================
! maskRusso : Computes a cell mask with ones for cells in Gamma and zero for
!             all other cells.
!===============================================================================
subroutine maskRusso(mask, phi0, nxp, nyp, nzp)
integer, intent(inout) :: mask(nzp, nxp, nyp)
integer, intent(in) :: nxp, nyp, nzp
real, intent(in) :: phi0(nzp, nxp, nyp)
integer :: k, i, j

mask(:,:,:) = 0

do j = 3, nyp-2
do i = 3, nxp-2
do k = 2, nzp-1
    if ( (phi0(k,i,j)*phi0(k,i-1,j) < 0.0).or. & 
         (phi0(k,i,j)*phi0(k,i+1,j) < 0.0).or. &
         (phi0(k,i,j)*phi0(k,i,j-1) < 0.0).or. & 
         (phi0(k,i,j)*phi0(k,i,j+1) < 0.0).or. &
         (phi0(k,i,j)*phi0(k+1,i,j) < 0.0).or. &
         (phi0(k,i,j)*phi0(k-1,i,j) < 0.0) ) then
        mask(k,i,j) = 1
    end if
end do
end do
end do

end subroutine maskRusso
!
!===============================================================================
! smootSign: Returns the maximum timestep to ensure stability.
!===============================================================================
real function smoothSign(phi)
real, intent(in) :: phi
smoothSign = phi / (phi**2 + sign_eps2)**0.5
end function smoothSign
!
!===============================================================================
! gradientBackwards: Computes backwards gradients of phi and puts them into
!                    phix, phiy, and phiz.
!===============================================================================
subroutine gradientBackwards(phix, phiy, phiz, phi, nxp, nyp, nzp, dxi, dyi, dzi)
real, intent (in) :: phi(nzp, nxp, nyp), dxi, dyi, dzi(nzp)
real, dimension(nzp, nxp, nyp), intent(inout) :: phix, phiy, phiz
integer, intent(in) :: nxp, nyp, nzp
integer :: k, i, j

phix(:,:,:) = 0.0
phiy(:,:,:) = 0.0
phiz(:,:,:) = 0.0

! I need gradients in the first layer of ghost cells for the Godunov
! Hamiltonian, so I index from i/j = 2 to nxp-1/nyp-1.
do j = 2, nyp-1
do i = 2, nxp-1
do k = 2, nzp-1
    phix(k,i,j) = dxi * (phi(k,i,j) - phi(k,i-1,j))
    phiy(k,i,j) = dyi * (phi(k,i,j) - phi(k,i,j-1))
    phiz(k,i,j) = dzi(k) * (phi(k,i,j) - phi(k-1,i,j))
end do
end do
end do

end subroutine gradientBackwards
!
!===============================================================================
! godunovHamiltonian: Computes the Hamiltonian |∇φ| and write it to 'G'. This
!                     routine takes phi0--the phi before reinitialization--and
!                     the backwards-differenced x and z derivatives phix and
!                     phiz as inputs.
!===============================================================================
subroutine godunovHamiltonian(G, phi0, phix, phiy, phiz, nxp, nyp, nzp)
real, intent(inout) :: G(nzp, nxp, nyp)
real, dimension(nzp, nxp, nyp), intent(in) :: phi0, phix, phiy, phiz
integer, intent(in) :: nxp, nyp, nzp
integer :: k, i, j

G(:,:,:) = 0.0

do j = 3, nyp-2
do i = 3, nxp-2
do k = 2, nzp-1
    G(k,i,j) = max( sign(1.0, phi0(k,i,j)), 0.0 ) * &
        ( max( max(phix(k,i,j), 0.0)**2, min(phix(k,i+1,j), 0.0)**2 ) + &
          max( max(phiy(k,i,j), 0.0)**2, min(phiy(k,i,j+1), 0.0)**2 ) + &
          max( max(phiz(k,i,j), 0.0)**2, min(phiz(k+1,i,j), 0.0)**2 ) )**0.5 &
             - min( sign(1.0, phi0(k,i,j)), 0.0 ) * &
        ( max( min(phix(k,i,j), 0.0)**2, max(phix(k,i+1,j), 0.0)**2 ) + &
          max( min(phiy(k,i,j), 0.0)**2, max(phiy(k,i,j+1), 0.0)**2 ) + &
          max( min(phiz(k,i,j), 0.0)**2, max(phiz(k+1,i,j), 0.0)**2 ) )**0.5
end do
end do
end do

end subroutine godunovHamiltonian

end module reinit
