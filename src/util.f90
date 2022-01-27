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
module util

  use mpi_interface, only : cyclics, cyclicc
  implicit none

  integer, save :: fftinix=0, fftiniy=0

contains
  ! ----------------------------------------------------------------------
  ! Subroutine sclrset: Sets upper and lower boundaries to a constant
  ! gradient via extrapolation, or a zero-gradient condition depending
  ! on the flag, typically used to set boundary conditions for scalars
  !
  subroutine sclrset(type,n1,n2,n3,a,dz,n)

    use mpi_interface, only : myid, appl_abort

    integer, intent(in)           :: n1,n2,n3
    integer, intent(in), optional :: n
    real, intent(in), optional    :: dz(n1)
    real, intent(inout)           :: a(n1,n2,n3)
    character (len=4)             :: type

    integer :: i,j,req(16)
    real :: dzf1,dzf2
    select case (type)
    case ('grad')
       dzf1  = dz(2)/dz(1)
       dzf2  = dz(n1-2)/dz(n1-1)
    case ('cnst')
       dzf1  = 0.
       dzf2  = 0.
    case ('mixd')
       dzf1  = 0.
       dzf2  = dz(n1-2)/dz(n1-1)
    case default
       if (myid == 0) print *, '  ABORTING:  BCs not supported'
       call appl_abort(0)
    end select

    do j=1,n3
       do i=1,n2
          a(1,i,j)  = a(2,i,j)    - dzf1*(a(3,i,j)    - a(2,i,j))
          a(n1,i,j) = a(n1-1,i,j) + dzf2*(a(n1-1,i,j) - a(n1-2,i,j))
       end do
    end do
    
    call cyclics(n1,n2,n3,a,req)
    call cyclicc(n1,n2,n3,a,req)

  end subroutine sclrset
  !
  ! ----------------------------------------------------------------------
  ! VELSET:  Sets boundary conditions for velocity
  !
  subroutine velset(n1,n2,n3,u,v,w)

    logical, parameter  :: noslip = .false.
    integer, intent(in) :: n1,n2,n3
    real, intent(inout) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

    integer :: i, j,req(16)

    call cyclics (n1,n2,n3,u,req)
    call cyclicc (n1,n2,n3,u,req)
    call cyclics (n1,n2,n3,v,req)
    call cyclicc (n1,n2,n3,v,req)
    call cyclics (n1,n2,n3,w,req)
    call cyclicc (n1,n2,n3,w,req)

    do j=1,n3
       do i=1,n2
          w(n1,i,j)    = 0.
          w(n1-1,i,j)  = 0.
          w(1,i,j)     = 0.
       enddo
    enddo
    if (noslip) then
       do j=1,n3
          do i=1,n2
             u(1,i,j)  = -u(2,i,j)
             u(n1,i,j) = -u(n1-1,i,j)
             v(1,i,j)  = -v(2,i,j)
             v(n1,i,j) = -v(n1-1,i,j)
          end do
       end do
    else
       do j=1,n3
          do i=1,n2
             u(1,i,j)  =  u(2,i,j)
             u(n1,i,j) = u(n1-1,i,j)
             v(1,i,j)  = v(2,i,j)
             v(n1,i,j) = v(n1-1,i,j)
          end do
       end do
    end if

  end subroutine velset
  !
  !---------------------------------------------------------------------
  ! GET_AVG: gets average field value from k'th level
  !
  real function get_avg(n1,n2,n3,k,a)

    integer, intent (in):: n1, n2, n3, k
    real, intent (in)   :: a(n1,n2,n3)

    integer :: i,j

    get_avg=0.
    do j=3,n3-2
       do i=3,n2-2
          get_avg = get_avg + a(k,i,j)
       enddo
    enddo
    get_avg = get_avg/real((n3-4)*(n2-4))

  end function get_avg
  !
  !---------------------------------------------------------------------
  ! GET_AVG3: gets average across outer two dimensions at each
  ! point along inner dimension
  !
  subroutine get_avg3(n1,n2,n3,a,avg)

    use mpi_interface, only : nypg,nxpg,real_array_par_sum

    integer,intent(in) :: n1,n2,n3
    real,intent(in) :: a(n1,n2,n3)
    real,intent(out) :: avg(n1)

    integer      :: k,i,j
    real :: lavg(n1),gavg(n1),x

    x = 1./(real(nypg-4)*real(nxpg-4))
    gavg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             gavg(k)=gavg(k)+a(k,i,j)
          end do
       end do
    end do
    lavg = gavg
    call real_array_par_sum(lavg,gavg,n1)
    avg(:) = real(gavg(:) * x)

  end subroutine get_avg3
  !
  !---------------------------------------------------------------------
  ! function get_cavg: gets conditionally average field from kth level
  !
  real function get_cavg(n1,n2,n3,k,a,x)

    integer, intent (in) :: n1,n2,n3,k
    real, intent (in)    :: a(n1,n2,n3),x(n2,n3)

    integer :: i,j
    real    :: avg, cnt

    avg=0.
    cnt=0.
    do j=3,n3-2
       do i=3,n2-2
          avg=avg+a(k,i,j)*x(i,j)
          cnt=cnt+x(i,j)
       end do
    end do
    get_cavg=avg/max(1.,cnt)

  end function get_cavg
  !
  !---------------------------------------------------------------------
  ! function get_csum: conditionally sum field a over indicator x
  !
  real function get_csum(n1,n2,n3,k,a,x)

    integer, intent (in) :: n1,n2,n3,k
    real, intent (in)    :: a(n1,n2,n3),x(n2,n3)

    integer :: i,j
    real    :: csum

    csum=0.
    do j=3,n3-2
       do i=3,n2-2
          csum=csum+a(k,i,j)*x(i,j)
       enddo
    enddo
    get_csum=csum

  end function get_csum
  !
  !---------------------------------------------------------------------
  ! function get_cor: gets mean correlation between two fields at a 
  ! given level
  !
  real function get_cor(n1,n2,n3,k,a,b)

    integer, intent (in) :: n1,n2,n3,k
    real, intent (inout) :: a(n1,n2,n3),b(n1,n2,n3)

    integer :: i,j
    real    :: avg
    avg=0.

    do j=3,n3-2
       do i=3,n2-2
          avg=avg+a(k,i,j)*b(k,i,j)
       end do
    end do
    get_cor=avg/real((n3-4)*(n2-4))

  end function get_cor
  !
  !---------------------------------------------------------------------
  ! function get_cor3: gets mean correlation accross outer two dimensions
  ! at each point along inner dimension
  !
  subroutine get_cor3(n1,n2,n3,a,b,avg)

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: a(n1,n2,n3),b(n1,n2,n3)
    real, intent (out)   :: avg(n1)

    integer :: k,i,j
    real    :: x

    do k=1,n1
       avg(k) = 0.
    end do
    x = 1./real((n3-4)*(n2-4))

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             avg(k)=avg(k)+a(k,i,j)*b(k,i,j)*x
          end do
       enddo
    enddo

  end subroutine get_cor3
  !
  !---------------------------------------------------------------------
  ! function get_var3: gets variance for a field whose mean is known
  !
  subroutine get_var3(n1,n2,n3,a,b,avg)


    use mpi_interface, only : nypg,nxpg,real_array_par_sum

    integer,intent(in) :: n1,n2,n3
    real,intent(in) :: a(n1,n2,n3)
    real,intent(in) :: b(n1)
    real,intent(out) :: avg(n1)

    integer      :: k,i,j
    real :: lavg(n1),gavg(n1),x

    x = 1./(real(nypg-4)*real(nxpg-4))
    gavg(:) = 0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             gavg(k)=gavg(k)+(a(k,i,j)-b(k))**2
          end do
       end do
    end do
    lavg = gavg
    call real_array_par_sum(lavg,gavg,n1)
    avg(:) = real(gavg(:) * x)


  end subroutine get_var3
  ! 
  !-------------------------------------------------------------------
  ! function get_var: gets square of a field from k'th level 
  ! 
  real function get_sqr(n1,n2,n3,k,a)

    integer, intent (in):: n1, n2, n3, k
    real, intent (in)   :: a(n1,n2,n3)

    integer :: i,j
    real    :: avg

    avg=0.
    do j=3,n3-2
       do i=3,n2-2
          avg=avg+a(k,i,j)**2
       end do
    end do
    get_sqr = avg/real((n3-4)*(n2-4))

  end function get_sqr
  !
  ! ----------------------------------------------------------------------
  ! subroutine tridiff: standard tri-diagonal solver for nh columns
  ! uses the LU decomposition of the equation 
  !
  !     cin1(k)*x(k-1)+ci(k)*x(k)+cip1(k)*x(k+1) = b(k)
  !
  ! such that
  !
  !         |  1   0   0   0   0 ... 0 | |u11 a12  0   0   0   0  ...  0  |
  !     L = | l21  1   0   0   0 ... 0 | | 0  u22 a23  0   0   0  ...  0  |
  !         |  0  l32  1   0   0 ... 0 | | 0   0  u33 a34  0   0  ...  0  |
  !         |                          | |                                |
  !         |  0   0   0   0   0 ... 1 | | 0   0   0   0   0   0  ... unn |
  !
  ! where aik =cip1(k) (i=k-1), =cin1(k) i=k+1
  !    u(1,1)   = a(1,1)
  !    l(i,i-1) = a(i,i-1)/u(i-1,i-1)
  !    u(i,i)   = a(i,i1)-l(i,i-1)*a(i-1,i)
  !
  ! and solves for x(k) = y(k)/u(k,k) - a(k,k+1)*x(k+1)/u(k,k) where
  ! y(k) = b(k) -l(k)y(k-1)
  !
  subroutine tridiff(nh,n1,nhdo,cin1,ci,cip1,rhs,cj,cjp1)

    integer, intent(in) :: nh,n1,nhdo
    real, intent(in)    :: cin1(nh,n1),ci(nh,n1),cip1(nh,n1)
    real, intent(inout) :: rhs(nh,n1),cj(nh,n1),cjp1(nh,n1)

    integer k,i
    real eps

    eps=sqrt(tiny(1.))

    do i=1,nhdo
       cjp1(i,2)=cip1(i,2)/ci(i,2)
       rhs(i,2)=rhs(i,2)/ci(i,2)
       do k=3,n1
          cj(i,k)=ci(i,k)-cin1(i,k)*cjp1(i,k-1)+eps
          cjp1(i,k)=cip1(i,k)/cj(i,k)
          rhs(i,k)=(rhs(i,k)-cin1(i,k)*rhs(i,k-1))/cj(i,k)
       enddo
       !
       ! here rhs = y(k)/u(k,k), cjp1=a(k,k+1)/u(k,k)
       !
       cj(i,n1)=rhs(i,n1)
       do k=n1-1,2,-1
          cj(i,k)=rhs(i,k)-cjp1(i,k)*cj(i,k+1)
       enddo
    end do

  end subroutine tridiff
  !
  ! --------------------------------------------------------------------
  !
  subroutine azero(n1, a1, a2, a3)

    integer, intent (in) :: n1
    real, intent (inout) :: a1(n1)
    real, optional, intent (inout) :: a2(n1), a3(n1)

    integer :: k

    do k=1,n1
       a1(k)=0.
       if (present(a2)) a2(k)=0.
       if (present(a3)) a3(k)=0.
    enddo

  end subroutine azero
  !
  ! --------------------------------------------------------------------
  ! subroutine ae1mm: subtracts mean value from given field (a=a-a_bar)
  !
  subroutine ae1mm(n1,n2,n3,a,abar)

    integer n1,n2,n3
    real, intent (inout), dimension (n1,n2,n3) :: a(n1,n2,n3)
    real, intent (out), dimension (n1)         :: abar(n1)

    integer :: i,j,k

    call get_avg3(n1,n2,n3,a,abar)

    do j=1,n3
       do i=1,n2
          do k=1,n1
             a(k,i,j)=a(k,i,j)-abar(k)
          enddo
       enddo
    enddo

  end subroutine ae1mm
  !
  ! --------------------------------------------------------------------
  ! subroutine atob: copies array a to array b
  !
  subroutine atob(nn,a,b)

    integer, intent (in) :: nn
    real, intent(in)     :: a(nn)
    real, intent(out)    :: b(nn)

    integer :: j

    do j=1,nn
       b(j)=a(j)
    end do

  end subroutine atob
  !
  !---------------------------------------------------------------------
  ! CRAYFFTUSE:  Uses the cray routines to do a 2D transform
  !
  subroutine get_fft_twodim(nx,ny,nz,a,wsavex,wsavey,isgn)

    use mpi_interface, only : xshuffle,yshuffle,nxg,nyg,nynzp,nxnzp

    integer, intent(in)    :: nx,ny,nz,isgn
    complex, intent(inout) :: a(nx,ny,nz)
    real, intent(inout)    :: wsavex(4*nxg+100),wsavey(4*nyg+100)

    integer :: k, j, i
    complex :: atmp((nx+1)*(ny+1)*(nz+1)),btmp(ny,nx,nz)

    call xshuffle(a,atmp,nx,ny,nz,1)
    call fft1dc(nxg,nynzp,atmp,wsavex,isgn,fftinix)
    call xshuffle(a,atmp,nx,ny,nz,-1)

    do k=1,nz
       do j=1,ny
          do i=1,nx
             btmp(j,i,k)=A(i,j,k)
          enddo
       enddo
    enddo

    call yshuffle(btmp,atmp,nx,ny,nz,1)
    call fft1dc(nyg,nxnzp,atmp,wsavey,isgn,fftiniy)
    call yshuffle(btmp,atmp,nx,ny,nz,-1)

    do k=1,nz
       do j=1,ny
          do i=1,nx
             A(i,j,k)= btmp(j,i,k)
          enddo
       enddo
    enddo

  end subroutine get_fft_twodim
  !
     
  subroutine calclevel(varin,varout,location, threshold)
    use grid, only : nzp, nxp, nyp, zt
    use mpi_interface, only : real_scalar_par_max, real_scalar_par_min
    real, intent(in), dimension(:,:,:) :: varin
    real, intent(in), optional :: threshold
    integer, intent(out) :: varout
    integer :: klocal
    real :: rlocal, rglobal
    character(*), intent(in) :: location
    integer :: i, j, k, km1
    real :: thres
    

     if (present(threshold)) then
      thres = threshold
    else
      thres = 0.0
    end if
   
    select case(location)
    case ('top')
      klocal = 0
      do j = 3, nyp - 2
        do i = 3, nxp - 2
          top:do k = nzp - 1, 2, -1
            if (varin(k,i,j) > thres) then
              klocal = max(klocal, k)
              exit top
            end if
          end do top
        end do
      end do
      rlocal = klocal
      call real_scalar_par_max(rlocal,rglobal) 
      varout = rglobal
    case ('base')
      klocal = nzp 
      do j = 3, nyp - 2
        do i = 3, nxp - 2
          base:do k = 2, nzp - 1
            if (varin(k,i,j) > thres) then
              klocal = min(klocal, k)
              exit base
            end if
          end do base
        end do
      end do
      rlocal = klocal
      call real_scalar_par_min(rlocal,rglobal) 
      varout = rglobal

    end select
  end subroutine calclevel

end module util
