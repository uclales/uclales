program test_reinit
    
    use reinit, only : cond_B, russo, in_Gamma
    implicit none
    call test_cond_B()
    call test_russo()

    contains

    ! This routine tests the truth table of cond_B() with synthetic input for
    ! each case of the table.
    ! 
    subroutine test_cond_B()
        integer, parameter :: n_tests = 8
        real, dimension(n_tests) :: phi_p2, phi_p1, phi, phi_m1, phi_m2
        logical, dimension(n_tests) :: res
        integer :: i

        ! truth table
        res = (/ .false., .true., .true., .true., &
                 .true.,  .true., .true., .true./)

        ! initialize values
        phi_p2 = 1e10
        phi_p1 = 1e10
        phi    = 1e10
        phi_m1 = 1e10
        phi_m2 = 1e10

        ! 0, 0, 0 => 0
        phi_p2(1) = 3.
        phi_p1(1) = 2.
        phi   (1) = 1.
        phi_m1(1) = .5
        phi_m2(1) = .25

        ! 0, 0, 1 => 1
        phi_p2(2) = -1.
        phi_p1(2) = 2.
        phi   (2) = 1.
        phi_m1(2) = .5
        phi_m2(2) = .25

        ! 0, 1, 0 => 1
        phi_p2(3) = 3.
        phi_p1(3) = 2.
        phi   (3) = 1.
        phi_m1(3) = .5
        phi_m2(3) = -2.

        ! 0, 1, 1 => 1
        phi_p2(4) = -1.
        phi_p1(4) = 2.
        phi   (4) = 1.
        phi_m1(4) = .5
        phi_m2(4) = -2.

        ! 1, 0, 0 => 1
        phi_p2(5) = 3.
        phi_p1(5) = 2.
        phi   (5) = 1.
        phi_m1(5) = 1.5
        phi_m2(5) = .25

        ! 1, 0, 1 => 1
        phi_p2(6) = -1.
        phi_p1(6) = 2.
        phi   (6) = 1.
        phi_m1(6) = 1.5
        phi_m2(6) = .25

        ! 1, 1, 0 => 1
        phi_p2(7) = 3.
        phi_p1(7) = 2.
        phi   (7) = 1.
        phi_m1(7) = 1.5
        phi_m2(7) = -1.

        ! 1, 1, 1 => 1
        phi_p2(8) = -3.
        phi_p1(8) = 2.
        phi   (8) = 1.
        phi_m1(8) = 1.5
        phi_m2(8) = -1.

        print *, "=============================================="
        print *, " Test 1: cond_B:"
        print *, "=============================================="
        do i = 1, n_tests
            if (res(i) .eqv. cond_B(phi_p2(i), phi_p1(i), phi(i), phi_m1(i), &
                                  phi_m2(i))) then
                print '(A, I2, A)', '    Test 1.',i , ' passed.'
            else
                print *, '  Test 1.',i, 'failed.'
            end if
        end do
        print *

    end subroutine test_cond_B


    subroutine test_russo()
        integer, parameter :: ntests = 5
        integer, dimension(ntests), parameter :: np_list = (/17, 33, 65, 129, 257/)
        integer :: gkmin, gkmax, i, j, k, l, n, np
        real, dimension(ntests) :: err_norm, err_norm2, max_err
        real :: dxi, dyi, dtau, h, r, err
        real, allocatable, dimension(:,:,:) :: phi0, phi, d
        real, allocatable, dimension (:)    :: xm, ym, zm, dzi

        print *, "=============================================="
        print *, " Test 2: Russo and Smereka's reinitialization"
        print *, "=============================================="
        do l = 1, ntests
            np = np_list(l)
            allocate (xm(np), ym(np), zm(np), dzi(np))
            allocate (phi0(np,np,np), phi(np,np,np), d(np,np,np))
            gkmin = 1
            gkmax = np

            ! The test domain is a (10 x 10 x 10) cube.
            h = 10.0 / (np-1)
            dxi    = 1.0 / h
            dyi    = 1.0 / h
            dzi(:) = 1.0 / h
            dtau   = 0.75 * h

            do i=1, np
                xm(i) = -5.0 + i*h
                ym(i) = -5.0 + i*h
                zm(i) = -5.0 + i*h
            end do

            r = 3
            do k=1, np
            do i=1, np
            do j=1, np
                phi0(k,i,j) = (0.1 + (xm(i) - r)**2 + (ym(j) - r)**2 + (zm(k) - r )**2) &
                            * (r - ( xm(i)**2 + ym(j)**2 + zm(k)**2 )**0.5)
                d(k,i,j)    = r - (xm(i)**2 + ym(j)**2 + zm(k)**2)**0.5
            end do
            end do
            end do
            phi(:,:,:) = phi0(:,:,:)

            print "(' on [', i3, ' x ', i3, ' x ', i3, ' ] grid, h = ', e22.15)", np-1, np-1, np-1, h
            ! Reinitialize
            call russo(phi, xm, ym, zm, np, np, np, gkmin, gkmax, dxi, dyi, dzi, dtau)

            ! Test zero level set displacement
            n = 0
            max_err(l)  = 0.0
            err_norm(l) = 0.0
            err_norm2(l) = 0.0
            do k=2, np-1
            do i=2, np-1
            do j=2, np-1
                if (in_Gamma(phi0, k, i, j, np, np, np)) then
                    n   = n+1
                    err = abs(phi(k,i,j) - d(k,i,j))
                    max_err(l)  = max(err, max_err(l))
                    err_norm(l) = err_norm(l) + err
                    err_norm2(l) = err_norm2(l) + err**2
                end if
            end do
            end do
            end do
            err_norm(l) = err_norm(l) / n
            err_norm2(l) = (err_norm2(l) / n)**0.5

            print '(A, E22.15, A, F6.2, A)', '    L_1:   ', err_norm(l), ' (', err_norm(l)*100/h, ' % of h.)'
            print '(A, E22.15, A, F6.2, A)', '    L_2:   ', err_norm2(l), ' (', err_norm2(l)*100/h, ' % of h.)'
            print '(A, E22.15, A, F6.2, A)', '    L_max: ', max_err(l), ' (', max_err(l)*100/h, ' % of h.)'
            print *

            deallocate(phi0, phi, d, xm, ym, zm, dzi)
        end do
 
        print *, "Order of convergence p:"
        do l = 1, ntests - 1
            print *, "    p(L_1)   = ", log(err_norm(l)/err_norm(l+1))/log(2.0)
        end do
        do l = 1, ntests - 1
            print *, "    p(L_2)   = ", log(err_norm2(l)/err_norm2(l+1))/log(2.0)
        end do
        do l = 1, ntests - 1
            print *, "    p(L_max) = ", log(max_err(l)/max_err(l+1))/log(2.0)
        end do

    end subroutine test_russo

end program
