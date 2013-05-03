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
!> Dynamic subgrid model
!>
!! \author Bart van Stratum

module sgsm_dyn
  implicit none
  include 'fftw3.f'

contains

  ! ---------------------------------------------------------------------
  ! Subroutine dynsub: dynamic subgrid scheme, should provide miracles
  ! for stable boundary layers....
  !
  ! Not the most efficient implementation, just for testing if miracles happen
  ! For now: only serial version....
  !
  subroutine calc_cs(timein)
    use grid, only             : a_up, a_vp, a_wp, nzp, nxp, nyp, dzi_t, dzi_m, dxi, dyi
    use mpi_interface, only    : nxg, nyg
    implicit none    

    real, intent(in)       :: timein 
    real, allocatable, dimension (:,:) :: uf1,uf2 
    integer :: k,i,j,ip,im,km,kp,jm,jp
    real :: alpha,wn1,wn2,Csz 

    real, allocatable, dimension(:,:) :: u_bar, v_bar, w_bar, u_hat, v_hat, w_hat, &
        S11n,      S12n ,      S13n ,      S22n,      S23n,      S33n, &
        S11_bar,   S12_bar,    S13_bar,    S22_bar,   S23_bar,   S33_bar, &
        S11_hat,   S12_hat,    S13_hat,    S22_hat,   S23_hat,   S33_hat, &
        S_S11_bar, S_S12_bar,  S_S13_bar,  S_S22_bar, S_S23_bar, S_S33_bar, & 
        S_S11_hat, S_S12_hat,  S_S13_hat,  S_S22_hat, S_S23_hat, S_S33_hat, & 
        Sn,  S_bar,S_hat, & 
        L11, L12, L13, L22, L23, L33, & 
        Q11, Q12, Q13, Q22, Q23, Q33, & 
        M11, M12, M13, M22, M23, M33, & 
        N11, N12, N13, N22, N23, N33, & 
        LM,  MM,  QN, NN 

    ! Constants etc
    alpha = sqrt(2.)  ! Basu & Porte-Agel, JAS2006
    wn1   = 1.
    wn2   = 1.

    allocate(u_bar(nxg,nyg),     v_bar(nxg,nyg),     w_bar(nxg,nyg),     u_hat(nxg,nyg),     v_hat(nxg,nyg),     w_hat(nxg,nyg))
    allocate(S11n(nxg,nyg),      S12n(nxg,nyg) ,     S13n(nxg,nyg) ,     S22n(nxg,nyg),      S23n(nxg,nyg),      S33n(nxg,nyg))
    allocate(S11_bar(nxg,nyg),   S12_bar(nxg,nyg),   S13_bar(nxg,nyg),   S22_bar(nxg,nyg),   S23_bar(nxg,nyg),   S33_bar(nxg,nyg))
    allocate(S11_hat(nxg,nyg),   S12_hat(nxg,nyg),   S13_hat(nxg,nyg),   S22_hat(nxg,nyg),   S23_hat(nxg,nyg),   S33_hat(nxg,nyg))
    allocate(S_S11_bar(nxg,nyg), S_S12_bar(nxg,nyg), S_S13_bar(nxg,nyg), S_S22_bar(nxg,nyg), S_S23_bar(nxg,nyg), S_S33_bar(nxg,nyg))
    allocate(S_S11_hat(nxg,nyg), S_S12_hat(nxg,nyg), S_S13_hat(nxg,nyg), S_S22_hat(nxg,nyg), S_S23_hat(nxg,nyg), S_S33_hat(nxg,nyg))
    allocate(Sn(nxg,nyg),        S_bar(nxg,nyg),     S_hat(nxg,nyg))
    allocate(L11(nxg,nyg),       L12(nxg,nyg),       L13(nxg,nyg),       L22(nxg,nyg),       L23(nxg,nyg),       L33(nxg,nyg)) 
    allocate(Q11(nxg,nyg),       Q12(nxg,nyg),       Q13(nxg,nyg),       Q22(nxg,nyg),       Q23(nxg,nyg),       Q33(nxg,nyg)) 
    allocate(M11(nxg,nyg),       M12(nxg,nyg),       M13(nxg,nyg),       M22(nxg,nyg),       M23(nxg,nyg),       M33(nxg,nyg)) 
    allocate(N11(nxg,nyg),       N12(nxg,nyg),       N13(nxg,nyg),       N22(nxg,nyg),       N23(nxg,nyg),       N33(nxg,nyg)) 
    allocate(LM(nxg,nyg),        MM(nxg,nyg),        QN(nxg,nyg),        NN(nxg,nyg)) 

    ! Large loop over height..
    do k=1,nzp 
      kp = k+1
      km = k-1

      ! Reset all vars to zero
      Sn(:,:)        = 0.
      S11n(:,:)      = 0.
      S12n(:,:)      = 0.
      S13n(:,:)      = 0.
      S22n(:,:)      = 0.
      S23n(:,:)      = 0.
      S33n(:,:)      = 0. 
      u_bar(:,:)     = 0.
      v_bar(:,:)     = 0.
      w_bar(:,:)     = 0.
      u_hat(:,:)     = 0.
      v_hat(:,:)     = 0.
      w_hat(:,:)     = 0.
      S11_bar(:,:)   = 0.
      S12_bar(:,:)   = 0.
      S13_bar(:,:)   = 0.
      S22_bar(:,:)   = 0.
      S23_bar(:,:)   = 0.
      S33_bar(:,:)   = 0. 
      S11_hat(:,:)   = 0.
      S12_hat(:,:)   = 0.
      S13_hat(:,:)   = 0.
      S22_hat(:,:)   = 0.
      S23_hat(:,:)   = 0.
      S33_hat(:,:)   = 0.
      S_S11_bar(:,:) = 0.
      S_S12_bar(:,:) = 0.
      S_S13_bar(:,:) = 0.
      S_S22_bar(:,:) = 0.
      S_S23_bar(:,:) = 0.
      S_S33_bar(:,:) = 0.
      S_S11_hat(:,:) = 0.
      S_S12_hat(:,:) = 0.
      S_S13_hat(:,:) = 0.
      S_S22_hat(:,:) = 0.
      S_S23_hat(:,:) = 0.
      S_S33_hat(:,:) = 0.
      S_bar(:,:)     = 0.
      S_hat(:,:)     = 0.
      L11(:,:)       = 0.
      L12(:,:)       = 0.
      L13(:,:)       = 0.
      L22(:,:)       = 0.
      L23(:,:)       = 0.
      L33(:,:)       = 0.
      Q11(:,:)       = 0.
      Q12(:,:)       = 0.
      Q13(:,:)       = 0.
      Q22(:,:)       = 0.
      Q23(:,:)       = 0.
      Q33(:,:)       = 0. 
      M11(:,:)       = 0.
      M12(:,:)       = 0.
      M13(:,:)       = 0.
      M22(:,:)       = 0.
      M23(:,:)       = 0.
      M33(:,:)       = 0. 
      N11(:,:)       = 0.
      N12(:,:)       = 0.
      N13(:,:)       = 0.
      N22(:,:)       = 0.
      N23(:,:)       = 0.
      N33(:,:)       = 0. 
      LM(:,:)        = 0.
      MM(:,:)        = 0.
      QN(:,:)        = 0.
      NN(:,:)        = 0.

      do j=3,nyp-2
        jp = j+1
        jm = j-1
        do i=3,nxp-2
          ip=i+1
          im=i-1
  
            ! Components Sij=(0.5*(dui/dxj+duj/dxi))^2 = (du/dx)^2
            ! S11,S22,S33 
            S11n(i,j) = ((a_up(k,i,j) - a_up(k,im,j)) * dxi     )**2.
            S22n(i,j) = ((a_vp(k,i,j) - a_vp(k,i,jm)) * dyi     )**2.
            S33n(i,j) = ((a_wp(k,i,j) - a_wp(km,i,j)) * dzi_t(k))**2.

            ! S12 = (0.5(du/dy+dv/dx))^2
            ! Averaging over 4 points: factor 0.25
            ! factor 0.5 in S12: 0.125
            ! Appears twice in SijSij: final factor 0.25...
            S12n(i,j) = 0.25 * ( &
              ((a_up(k,i,jp)  -  a_up(k,i,j))   * dyi      + &
               (a_vp(k,ip,j)  -  a_vp(k,i,j))   * dxi       )**2.    + &
              ((a_up(k,i,j)   -  a_up(k,i,jm))  * dyi      + &
               (a_vp(k,ip,jm) -  a_vp(k,i,jm))  * dxi       )**2.    + &
              ((a_up(k,im,j)  -  a_up(k,im,jm)) * dyi      + &
               (a_vp(k,i,jm)  -  a_vp(k,im,jm)) * dxi       )**2.    + &
              ((a_up(k,im,jp) -  a_up(k,im,j))  * dyi      + &
               (a_vp(k,i,j)   -  a_vp(k,im,j))  * dxi       )**2.    )

            ! S13 = 0.5(dw/dx+du/dz)^2
            S13n(i,j) = 0.25 * ( &
              ((a_wp(k,ip,j)  -  a_wp(k,i,j))   * dxi      + &
               (a_up(kp,i,j)  -  a_up(k,i,j))   * dzi_m(k)  )**2.    + &
              ((a_wp(km,ip,j) -  a_wp(km,i,j))  * dxi      + &
               (a_up(k,i,j)   -  a_up(km,i,j))  * dzi_m(km) )**2.    + &
              ((a_wp(km,i,j)  -  a_wp(km,im,j)) * dxi      + &
               (a_up(k,im,j)  -  a_up(km,im,j)) * dzi_m(km) )**2.    + &
              ((a_wp(k,i,j)   -  a_wp(k,im,j))  * dxi      + &
               (a_up(kp,im,j) -  a_up(k,im,j))  * dzi_m(k)  )**2.    )

            ! S23 = (0.5(dv/dz+dw/dyi))^2
            S23n(i,j) = 0.25 * ( &
              ((a_vp(kp,i,j)  -  a_vp(k,i,j))   * dzi_m(k) + &
               (a_wp(k,i,jp)  -  a_wp(k,i,j))   * dyi       )**2.    + &
              ((a_vp(k,i,j)   -  a_vp(km,i,j))  * dzi_m(km)+ &
               (a_wp(km,i,jp) -  a_wp(km,i,j))  * dyi       )**2.    + &
              ((a_vp(k,i,jm)  -  a_vp(km,i,jm)) * dzi_m(km)+ &
               (a_wp(km,i,j)  -  a_wp(km,i,jm)) * dyi       )**2.    + &
              ((a_vp(kp,i,jm) -  a_vp(k,i,jm))  * dzi_m(k) + &
               (a_wp(k,i,j)   -  a_wp(k,i,jm))  * dyi       )**2.    )

            ! unstagger velocities to strain point 
            u_bar(i,j) = 0.5 * (a_up(k,i,j) + a_up(k,im,j))
            v_bar(i,j) = 0.5 * (a_vp(k,i,j) + a_vp(k,i,jm))
            w_bar(i,j) = 0.5 * (a_wp(k,i,j) + a_wp(km,i,j)) 

        end do
      end do
  
      ! S=sqrt(2SijSij) 
      Sn(:,:) = sqrt(2.*(S11n(:,:) + S22n(:,:) + S33n(:,:) + S12n(:,:) + S13n(:,:) + S23n(:,:)))

      ! Weekend :)

    end do




    !allocate(uf1(nxg,nyg),uf2(nxg,nyg))
    !uf1(:,:) = a_up(1,3:nxg+2,3:nyg+2)
    !uf2(:,:) = a_up(1,3:nxg+2,3:nyg+2)

    !if(timein>60) then
    !  call sfilter(uf1,32,nxg,nyg) 
    !  do i=1,nxg
    !    print*,uf2(i,10),uf1(i,10)
    !  end do 

    !  stop
    !end if

    
    deallocate(u_bar, v_bar, w_bar, u_hat, v_hat, w_hat, &
        S11n,      S12n ,      S13n ,      S22n,      S23n,      S33n, &
        S11_bar,   S12_bar,    S13_bar,    S22_bar,   S23_bar,   S33_bar, &
        S11_hat,   S12_hat,    S13_hat,    S22_hat,   S23_hat,   S33_hat, &
        S_S11_bar, S_S12_bar,  S_S13_bar,  S_S22_bar, S_S23_bar, S_S33_bar, & 
        S_S11_hat, S_S12_hat,  S_S13_hat,  S_S22_hat, S_S23_hat, S_S33_hat, & 
        Sn,  S_bar,S_hat, & 
        L11, L12, L13, L22, L23, L33, & 
        Q11, Q12, Q13, Q22, Q23, Q33, & 
        M11, M12, M13, M22, M23, M33, & 
        N11, N12, N13, N22, N23, N33, & 
        LM,  MM,  QN, NN) 


  end subroutine calc_cs 

  ! ---------------------------------------------------------------------
  ! Subroutine sfilter: spectral cut-off filter using 2d FFTW routines
  !
  subroutine sfilter(var,kf)
    use mpi_interface, only    : nxg, nyg
    implicit none

    integer, intent(in)                       :: kf         !< cutoff wavenumber
    real, dimension(:,:), intent(inout)       :: var        !< variable to filter
    !integer, intent(in)                       :: nxin,nyin  !< dimension of input

    double precision in
    dimension in(nxg,nyg)
    double complex out
    dimension out(nxg/2+1,nyg)
    integer*8 plan, plan2

    in(:,:) = var(:,:)
    ! Recycle FFT plan??
    call dfftw_plan_dft_r2c_2d(plan,nxg,nyg,in,out,FFTW_ESTIMATE)    ! Make plan forward dft
    call dfftw_execute_dft_r2c(plan, in, out)
    call dfftw_destroy_plan(plan)

    out(kf+1:,:)          = cmplx(0,0.0)  !-> OK
    out(:,kf+1:nyg-kf+1)  = cmplx(0,0.0)  !-> OK

    ! Recycle FFT plan??
    call dfftw_plan_dft_c2r_2d(plan2,nxg,nyg,out,in,FFTW_ESTIMATE)   ! Make plan backward dft
    call dfftw_execute_dft_c2r(plan2, out, in)
    call dfftw_destroy_plan(plan2)
    
    var(:,:) = in(:,:)/(nxg*nyg)  ! data returned scaled by n elements, see fftw3 docs   

  end subroutine sfilter

end module sgsm_dyn

