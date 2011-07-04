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
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA,
! and TV Singh, Academic and Technology Services
!----------------------------------------------------------------------------
!
module mpi_interface
  !
  !    nxg = nxpg-4
  !    nyg = nypg-4
  !    xcomm, commxid - communicator for x side processors and rank wrt it
  !    ycomm, commyid - communicator for y side processors and rank wrt it
  !    nxnzp = nx*nzp
  !    nynzp = ny*nzp
  !    wrxid, wryid, nxprocs,nyprocs:(wrxid,wryid)=myid 
  !       in ranktable (nxprocs,nyprocs)
  !    nxpa,nypa: arrays containing nxp and nyp for all nxprocs and nyprocs resp.
  !    nynza, nxnza: arrays containing nynzp and nxnzp on nxprocs and nyprocs 
  !    resp.
  !
  implicit none

  integer :: myid, pecount, nxpg, nypg, nxg, nyg, nbytes, intsize
  integer :: xcomm, ycomm,commxid,commyid, MY_CMPLX, MY_SIZE
  integer :: nxnzp,nynzp
  integer :: wrxid, wryid, nxprocs, nyprocs
  integer, allocatable, dimension(:) :: xoffset, yoffset, nxpa, nypa

  ! these are the parameters used in the alltoallw call in the fft

contains
  !
  !----------------------------------------------------------------------
  ! INIT_MP: Initializes MPI
  !
  subroutine init_mpi

    integer ierror
    character (len=8) date

    myid=0
    pecount=1

    select case (kind(0.0))
    case (4)
       nbytes = 4
       MY_SIZE = 4
       MY_CMPLX = 8
    case (8)
       nbytes = 8
       MY_SIZE = 8
       MY_CMPLX = 16
    case default
       stop "kind not supported"
    end select

    select case(kind(0))
    case (4)
       intsize=4
    case (8)
       intsize=8
    case default
       stop "int kind not supported"
    end select
    !
    call date_and_time(date)
    if (myid == 0) print "(/1x,75('-'),/2x,A13,A8,/2x,A15,I2,A15,I2,A14)", &
         'UCLA LES 2.0 ',date, 'Computing using',nbytes,' byte reals and', &
         intsize," byte integers"

  end subroutine init_mpi
  !
  !----------------------------------------------------------------------
  ! DEFINE_DECOMP: Defines MPI Decomposition
  !
  subroutine define_decomp(nxp, nyp, nxpart)


    integer, intent(inout) :: nxp, nyp
    logical, intent(in) ::  nxpart


    nxprocs=1
    nyprocs=1

    !
    !   ranktable is the matrix having ranks of processes in x-y domain
    !
    wrxid = 0
    wryid = 0
    commxid=0
    commyid=0

    !
    ! there are two boundary points in each direction
    !
    nxpg = nxp
    nypg = nyp

    nxg=nxpg-4
    nyg=nypg-4

    allocate (nxpa(0:nxprocs-1), nypa(0:nyprocs-1))
    nxpa(0)=nxg
    nypa(0)=nyg

    !
    !  offsets for ecah processor in x and y, for a given grid (nxp x nyp)
    !
    allocate(xoffset(0:nxprocs-1),yoffset(0:nyprocs-1))

    xoffset = 0
    yoffset = 0

    if(nxp.lt.5) then
       print *, 'ABORT: X Horizontal domain size too small for ',nxprocs,    &
            ' processors.'
       print *, '       Increase nyp to ',nxprocs*5, ' or run on ',nxpg/5,   &
            ' or fewer processors'
       call appl_abort(0)
    endif
    if(nyp.lt.5) then
       print *, 'ABORT: Y Horizontal domain size too small for ',nyprocs,    &
            ' processors.'
       print *, '       Increase nyp to ',nyprocs*5, ' or run on ',nypg/5,   &
            ' or fewer processors'
       call appl_abort(0)
    endif

    if (myid == 0) print 61,'Sequential (shared memory) simulation'

61 format (/1x,49('-')/2x,A37)

  end subroutine define_decomp
  !
  !----------------------------------------------------------------------
  ! INIT_ALLTOALL_REORDERXY: Defines the mpi derived types to do a data 
  ! movement of the form A(m,n/p,z) -> B(n,m/p,z) for data of type MY_CMPLX
  !
  subroutine init_alltoall_reorder(nxp,nyp,nzp)

    integer, intent(in) :: nxp,nyp,nzp

    nxg=nxp-4
    nyg=nyp-4
    nxnzp=nxg*nzp
    nynzp=nyg*nzp

  end subroutine init_alltoall_reorder
  ! ---------------------------------------------------------------------
  ! subroutine cyclics: commits exchange cyclic x boundary conditions
  !
  subroutine cyclics(n1,n2,n3,var,req)

    integer, intent(in) :: n1,n2,n3,req(16)
    real, intent(inout) :: var(n1,n2,n3)

    if (n3 == 5) then
       var(:,:,1) = var(:,:,3)
       var(:,:,2) = var(:,:,3)
       var(:,:,4) = var(:,:,3)
       var(:,:,5) = var(:,:,3)
    end if
    if (n2 == 5) then
       var(:,1,:) = var(:,3,:)
       var(:,2,:) = var(:,3,:)
       var(:,4,:) = var(:,3,:)
       var(:,5,:) = var(:,3,:)
    end if

    var(:,:2,:) = var(:,n2-3:n2-2,:)
    var(:,n2-1:,:) = var(:,3:4,:)

    var(:,:,:2) = var(:,:,n3-3:n3-2)
    var(:,:,n3-1:) = var(:,:,3:4)
    var(:,:2,:2) = var(:,n2-3:n2-2,n3-3:n3-2)
    var(:,n2-1:,n3-1:) = var(:,3:4,3:4)
    
  end subroutine cyclics
  !
  ! ---------------------------------------------------------------------
  ! subroutine cyclicc: comits excahnging cyclic boundary conditions
  subroutine cyclicc(n1,n2,n3,var,req)

    integer :: n1,n2,n3,req(16)
    real :: var(n1,n2,n3)

  end subroutine cyclicc
  !
  ! ---------------------------------------------------------------------
  subroutine appl_abort(ierr)

    integer ierr
    stop 'Program Aborted'

  end subroutine appl_abort
  !
  ! ---------------------------------------------------------------------
  subroutine appl_finalize(ierr)

    integer ierr

  end subroutine appl_finalize
  !
  !---------------------------------------------------------------------------
  subroutine xshuffle(a,atmp,nx,ny,nz,isign)

    integer, intent(in):: nx,ny,nz,isign
    complex, intent(inout):: a(nx,ny,nz),atmp((nx+1)*(ny+1)*(nz+1))
    integer ll,i,j,k

    if(isign .eq. 1) then
       ll=0
       do k=1,nz
          do j=1,ny
             do i=1,nx
                ll=ll+1
                atmp(ll)=a(i,j,k)
             enddo
          enddo
       enddo

    else
       ll=0
       do k=1,nz
          do j=1,ny
             do i=1,nx
                ll=ll+1
                a(i,j,k)=atmp(ll)
             enddo
          enddo
       enddo

    endif

  end subroutine xshuffle
  !
  !---------------------------------------------------------------------------
  subroutine yshuffle(a,atmp,nx,ny,nz,isign)

    integer, intent(in):: nx,ny,nz,isign
    complex, intent(inout):: a(ny,nx,nz),atmp((nx+1)*(ny+1)*(nz+1))
    integer ierr,ll,i,j,k

    if(isign .eq. 1) then
       ll=0
       do k=1,nz
          do j=1,nx
             do i=1,ny
                ll=ll+1
                atmp(ll)=a(i,j,k)
             enddo
          enddo
       enddo
    else
       ll=0
       do k=1,nz
          do j=1,nx
             do i=1,ny
                ll=ll+1
                a(i,j,k)=atmp(ll)
             enddo
          enddo
       enddo

    endif

  end subroutine yshuffle
  !
  !---------------------------------------------------------------------------
  ! get maximum across processors
  !
  subroutine double_scalar_par_max(xxl,xxg)

    real(kind=8), intent(out) :: xxg
    real(kind=8), intent(in) :: xxl

    xxg=xxl

  end subroutine double_scalar_par_max
  !
  !---------------------------------------------------------------------------
  ! get maximum across processors
  !
  subroutine double_scalar_par_min(xxl,xxg)

    real(kind=8), intent(out) :: xxg
    real(kind=8), intent(in) :: xxl

    xxg=xxl

  end subroutine double_scalar_par_min
  !
  !---------------------------------------------------------------------------
  subroutine double_scalar_par_sum(xxl,xxg)

    real(kind=8), intent(out) :: xxg
    real(kind=8), intent(in) :: xxl

    xxg=xxl

  end subroutine double_scalar_par_sum
  !
  !---------------------------------------------------------------------------
  subroutine double_array_par_sum(xxl,xxg,n)

    integer, intent(in)::n
    real(kind=8), intent(out) :: xxg(n)
    real(kind=8), intent(in) :: xxl(n)

    xxg=xxl

  end subroutine double_array_par_sum
  
  subroutine broadcast(val, n, procsend)
   integer, intent(in) :: n, procsend
   real(kind=8), intent(inout) :: val(n)
   
  end subroutine broadcast


end module mpi_interface
