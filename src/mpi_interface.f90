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

  use mpi
  implicit none
  !
  !    nxg = nxpg-4
  !    nyg = nypg-4
  !    xcomm, commxid - communicator for x side processors and rank wrt it
  !    ycomm, commyid - communicator for y side processors and rank wrt it
  !    nxnzp = nx*nzp
  !    nynzp = ny*nzp
  !    wrxid, wryid, nxprocs,nyprocs: (wrxid,wryid)=myid in 
  !        ranktable (nxprocs,nyprocs)
  !    nxpa,nypa: arrays containing nxp and nyp for all nxprocs and nyprocs 
  !         respectively
  !    nynza, nxnza: arrays containing nynzp and nxnzp on nxprocs and nyprocs 
  !         respectively
  !

  integer :: myid, pecount, nxpg, nypg, nxg, nyg, nbytes, intsize, &
       MY_SIZE, MY_CMPLX
  integer :: xcomm, ycomm,commxid,commyid
  integer :: nxnzp,nynzp,fftinix,fftiniy
  integer :: wrxid, wryid, nxprocs, nyprocs
  integer, allocatable, dimension(:) :: xoffset, yoffset, nxpa, nypa, &
       nynza, nxnza

  ! these are the parameters used in the alltoallw call in the fft

  integer, allocatable, dimension(:,:) :: ranktable,xtype,ytype,xdisp,&
       ydisp,xcount,ycount

  integer :: stridetype,xstride,ystride,xystride,xylarry,xyzlarry,&
       fxytype,fxyztype

contains
  !
  !----------------------------------------------------------------------
  ! INIT_MP: Initializes MPI
  !
  subroutine init_mpi

    integer ierror
    character (len=8) date

    call mpi_init(ierror)  
    call mpi_comm_size(MPI_COMM_WORLD, pecount, ierror)    
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierror)

    select case (kind(0.0))
    case (4)
       nbytes = 4
       MY_SIZE = MPI_REAL
       MY_CMPLX = MPI_COMPLEX
    case (8)
       nbytes = 8
       MY_SIZE = MPI_DOUBLE_PRECISION
       MY_CMPLX = MPI_doUBLE_COMPLEX
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

    integer :: ierror, i,j, modx,mody, nxpj, nypj, irank
    integer :: worldgroup, xgroup, ygroup

    nxprocs=1
    nyprocs=pecount
    if (nyp .gt. nxp) then
       nyprocs=pecount
       nxprocs=1
       if( (nyp-4)/nyprocs.lt.5 .or. nxpart) then
          nyprocs=int(sqrt(real(pecount)))
          nxprocs=pecount/nyprocs
          do while (nyprocs*nxprocs .ne. pecount)
             nyprocs=nyprocs+1
             nxprocs=pecount/nyprocs
          end do

          if(nxprocs.gt.nyprocs)then
             i=nxprocs
             nxprocs=nyprocs
             nyprocs=i
          endif

          if(nxp .lt. 5) then
             print *, 'ABORTING: NXP too small, increase to at least 5'
             call mpi_abort(MPI_COMM_WORLD,0,ierror)
          elseif( (nxp-4)/nxprocs .le. 5) then
             nxprocs=1
             nyprocs=pecount
          endif

          if ( (nyp-4)/nyprocs .lt. 5) then
             print *, '  ABORTING: NYP too small for ',nyprocs,' processors.'
             print *, '  Increase to ',nyprocs*9, ' or run on ',nypg/9,       &
                  ' or fewer processors'
             call mpi_abort(MPI_COMM_WORLD,0,ierror)
          endif
       endif

    else

       nxprocs=pecount
       nyprocs=1
       if( (nxp-4)/nxprocs.lt.5 .or. nxpart) then
          nxprocs=int(sqrt(real(pecount)))
          nyprocs=pecount/nxprocs
          do while (nyprocs*nxprocs .ne. pecount)
             nxprocs=nxprocs+1
             nyprocs=pecount/nxprocs
          end do

          if(nyprocs.gt.nxprocs)then
             i=nyprocs
             nyprocs=nxprocs
             nxprocs=i
          endif

          if(nyp .lt. 5) then
             print *, 'ABORTING: NYP too small, increase to at least 5'
             call mpi_abort(MPI_COMM_WORLD,0,ierror)
          elseif( (nyp-4)/nyprocs.le. 5) then
             nyprocs=1
             nxprocs=pecount
          endif

          if ( (nxp-4)/nxprocs .lt. 5) then
             print *, '  ABORTING: NXP too small for ',nxprocs,' processors.'
             print *, '  Increase to ',nxprocs*9, ' or run on ',nxpg/9,       &
                  ' or fewer processors'
             call mpi_abort(MPI_COMM_WORLD,0,ierror)
          endif
       endif
    endif

    if( (nyp-4)/nyprocs.lt.5 .or. nxpart) then
       nyprocs=int(sqrt(real(pecount)))
       nxprocs=pecount/nyprocs
       do while (nyprocs*nxprocs .ne. pecount)
          nyprocs=nyprocs+1
          nxprocs=pecount/nyprocs
       end do
       if(nxprocs.gt.nyprocs)then
          i=nxprocs
          nxprocs=nyprocs
          nyprocs=i
       endif
    endif

    !ctvs    nxprocs=1
    !ctvs    nyprocs=2
    !
    !   ranktable is the matrix having ranks of processes in x-y domain
    !

    allocate(ranktable(-1:nxprocs,-1:nyprocs))

    do i = -1, nxprocs
       do j = -1, nyprocs
          ranktable(i,j) = 0
       enddo
    enddo
    irank = 0
    do j = 0, nyprocs - 1
       do i = 0, nxprocs - 1
          ranktable(i,j) = irank
          if (myid .eq. irank) THEN
             wrxid = i
             wryid = j
          endif
          irank = irank + 1
       enddo
    enddo
    do i = 0, nxprocs - 1
       ranktable(i, -1) = ranktable(i, nyprocs - 1)
       ranktable(i, nyprocs) = ranktable(i, 0)
    enddo
    do j = 0, nyprocs-1 
       ranktable(-1, j) = ranktable(nxprocs - 1, j)
       ranktable(nxprocs, j) = ranktable(0, j)
    enddo
    ranktable(-1,-1)=ranktable(nxprocs-1,nyprocs-1)
    ranktable(nxprocs,nyprocs)=ranktable(0,0)
    ranktable(-1,nyprocs)=ranktable(nxprocs-1,0)
    ranktable(nxprocs,-1)=ranktable(0,nyprocs-1)

    call mpi_comm_group(mpi_comm_world, worldgroup, ierror)
    call mpi_group_incl(worldgroup, nxprocs,ranktable(0:nxprocs-1,wryid),&
         xgroup,ierror)
    call mpi_comm_create(mpi_comm_world, xgroup,xcomm,ierror)
    call mpi_group_incl(worldgroup, nyprocs,ranktable(wrxid,0:nyprocs-1),&
         ygroup,ierror)
    call mpi_comm_create(mpi_comm_world, ygroup,ycomm,ierror)

    call mpi_comm_rank(xcomm,commxid,ierror)
    call mpi_comm_rank(ycomm,commyid,ierror)

    !
    ! there are two boundary points in each direction
    !
    nxpg = nxp
    nypg = nyp
    nxp = (nxpg-4)/nxprocs + 4
    nyp = (nypg-4)/nyprocs + 4

    modx = modulo(nxpg-4,nxprocs)
    mody = modulo(nypg-4,nyprocs)
    !
    ! offsets for each processor in x and y directons(nxp x nyp)
    !
    allocate(xoffset(0:nxprocs-1),yoffset(0:nyprocs-1))

    xoffset = 0

    do j=1,nxprocs-1
       if(j <= modx) then
          nxpj = nxp+1
       else
          nxpj = nxp
       endif
       xoffset(j) = xoffset(j-1)+nxpj-4
    enddo

    yoffset = 0

    do j=1,nyprocs-1
       if(j <= mody) then
          nypj = nyp+1
       else
          nypj = nyp
       endif
       yoffset(j) = yoffset(j-1)+nypj-4
    enddo

    if (nxpg > 5 .and. nxp == 5) then
       print *, 'ABORTING: Subdomain too finely discretized in x', nxpg, nxp
       call mpi_abort(MPI_COMM_WORLD,0,ierror)
    end if

    if (nypg > 5 .and. nyp == 5) then
       print *, 'ABORTING: Subdomain too finely discretized in y', nypg, nyp
       call mpi_abort(MPI_COMM_WORLD,0,ierror)
    end if

    if(nxp.lt.5) then
       print *, 'ABORT: X Horizontal domain size too small for ',nxprocs,    &
            ' processors.'
       print *, '       Increase nyp to ',nxprocs*5, ' or run on ',nxpg/5, &
            ' or fewer processors'
       call mpi_abort(MPI_COMM_WORLD,0,ierror)
    endif
    if(nyp.lt.5) then
       print *, 'ABORT: Y Horizontal domain size too small for ',nyprocs,    &
            ' processors.'
       print *, '       Increase nyp to ',nyprocs*5, ' or run on ',nypg/5, &
            ' or fewer processors'
       call mpi_abort(MPI_COMM_WORLD,0,ierror)
    endif

    if (myid == 0) then
       print 61, 'Processor count', pecount,'nxpl =', nxp,' nypl = ',nyp
       do i=0,min(nxprocs,nyprocs)-1
          print "(2x,A13,2I5)", 'x/y offset = ', xoffset(i), yoffset(i)
       end do
       if (nxprocs>nyprocs) print "(15x,I5)", xoffset(nyprocs:nxprocs-1)
       if (nxprocs<nyprocs) print "(15x,I5)", yoffset(nxprocs:nyprocs-1)
    end if

61 format (/1x,49('-')/2x,A15,I5,2(A6,I5))

  end subroutine define_decomp
  !
  !----------------------------------------------------------------------
  ! INIT_ALLTOALL_REORDERXY: Defines the mpi derived types to do a data 
  ! movement of the form A(m,n/p,z) -> B(n,m/p,z) for data of type MY_CMPLX
  !
  subroutine init_alltoall_reorder(nxp,nyp,nzp)

    integer, intent(in) :: nxp,nyp,nzp

    integer :: nx, ny, i, j, k, ii, jj, ierr, cnt, typesize,nynzg, nxnzg


    nx = max(1,nxp-4)
    ny = max(1,nyp-4)
    nxg=nxpg-4
    nyg=nypg-4

    allocate (nxpa(0:nxprocs-1), nypa(0:nyprocs-1))
    allocate (nxnza(0:nyprocs-1), nynza(0:nxprocs-1))
    allocate(xcount(0:nxprocs-1,2),xtype(0:nxprocs-1,2),xdisp(0:nxprocs-1,2))
    allocate(ycount(0:nyprocs-1,2),ytype(0:nyprocs-1,2),ydisp(0:nyprocs-1,2))

    ii = nxg/nxprocs
    ii= nxg - nxprocs*ii

    do i=0,nxprocs-1
       nxpa(i)=nxg/nxprocs
       if(i .lt. ii) nxpa(i)=nxpa(i)+1
    enddo

    jj = nyg/nyprocs
    jj= nyg - nyprocs*jj

    do i=0,nyprocs-1
       nypa(i)=nyg/nyprocs
       if(i .lt. jj) nypa(i)=nypa(i)+1
    enddo

    nynzg=ny*nzp    
    ii = nynzg/nxprocs
    ii= nynzg - nxprocs*ii
    do i=0,nxprocs-1
       nynza(i) = nynzg/nxprocs
       if (i .lt. ii) nynza(i)=nynza(i)+1
    enddo

    nxnzg=nx*nzp    
    jj = nxnzg/nyprocs
    jj= nxnzg - nyprocs*jj
    do i=0,nyprocs-1
       nxnza(i) = nxnzg/nyprocs
       if (i .lt. jj) nxnza(i)=nxnza(i)+1
    enddo

    nxnzp=nxnza(wryid)
    nynzp=nynza(wrxid)


    call MPI_TYPE_SIZE(MY_CMPLX,typesize, ierr)

    xcount=1
    ycount=1
    xdisp=0
    ydisp=0

    do i=1,nxprocs-1
       xdisp(i,1)=xdisp(i-1,1)+nx*nynza(i-1)*typesize
       xdisp(i,2)=xdisp(i-1,2)+nxpa(i-1)*typesize
    enddo

    do i=1,nyprocs-1
       ydisp(i,1)=ydisp(i-1,1)+ny*nxnza(i-1)*typesize
       ydisp(i,2)=ydisp(i-1,2)+nypa(i-1)*typesize
    enddo

    do i=0,nxprocs-1
       call mpi_type_contiguous(nx*nynza(i),MY_CMPLX, xtype(i,1),ierr)
       call mpi_type_commit(xtype(i,1),ierr)
       call mpi_type_vector(nynzp, nxpa(i), nxg, MY_CMPLX, xtype(i,2),ierr)
       call mpi_type_commit(xtype(i,2),ierr)
    enddo

    do i=0,nyprocs-1
       call mpi_type_contiguous(ny*nxnza(i),MY_CMPLX, ytype(i,1),ierr)
       call mpi_type_commit(ytype(i,1),ierr)
       call mpi_type_vector(nxnzp, nypa(i), nyg, MY_CMPLX, ytype(i,2),ierr)
       call mpi_type_commit(ytype(i,2),ierr)
    enddo

    call MPI_TYPE_VECTOR(nyp-4,nzp*2,nxp*nzp,MY_SIZE,stridetype,ierr)
    call MPI_TYPE_COMMIT(stridetype,ierr)
    call MPI_TYPE_VECTOR(nyp-4,nzp*2,nxp*nzp,MY_SIZE,xstride,ierr)
    call MPI_TYPE_COMMIT(xstride,ierr)
    call MPI_TYPE_VECTOR(2,nzp*(nxp-4),nxp*nzp,MY_SIZE,ystride,ierr)
    call MPI_TYPE_COMMIT(ystride,ierr)
    call MPI_TYPE_VECTOR(2,2*nzp,nxp*nzp,MY_SIZE,xystride,ierr)
    call MPI_TYPE_COMMIT(xystride,ierr)

    call MPI_TYPE_VECTOR(nyp-4,nxp-4,nxpg-4,MY_SIZE,fxytype,ierr)
    call MPI_TYPE_COMMIT(fxytype,ierr)

    call MPI_TYPE_VECTOR(nyp-4,(nxp-4)*nzp,(nxpg-4)*nzp,MY_SIZE,fxyztype,ierr)
    call MPI_TYPE_COMMIT(fxyztype,ierr)

    call MPI_TYPE_VECTOR(nyp-4,nxp-4,nxp,MY_SIZE,xylarry,ierr)
    call MPI_TYPE_COMMIT(xylarry,ierr)
    call MPI_TYPE_VECTOR(nyp-4,(nxp-4)*nzp,nxp*nzp,MY_SIZE,xyzlarry,ierr)
    call MPI_TYPE_COMMIT(xyzlarry,ierr)



  end subroutine init_alltoall_reorder

  ! ---------------------------------------------------------------------
  ! subroutine cyclics: commits exchange cyclic x boundary conditions
  !
  subroutine cyclics(n1,n2,n3,var,req)

    integer, intent(in) :: n1,n2,n3
    real, intent(inout) :: var(n1,n2,n3)
    integer req(16)
    integer :: ierror, stats(MPI_STATUS_SIZE,16), pxfwd, pxback, pyfwd, pyback
    integer :: pxyne,pxyse,pxynw,pxysw

    if (nypg == 5) then
       var(:,:,1) = var(:,:,3)
       var(:,:,2) = var(:,:,3)
       var(:,:,4) = var(:,:,3)
       var(:,:,5) = var(:,:,3)
    end if
    if (nxpg == 5) then
       var(:,1,:) = var(:,3,:)
       var(:,2,:) = var(:,3,:)
       var(:,4,:) = var(:,3,:)
       var(:,5,:) = var(:,3,:)
    end if

    pxfwd = ranktable(wrxid+1,wryid)
    pxback =ranktable(wrxid-1,wryid)
    pyfwd = ranktable(wrxid,wryid+1)
    pyback =ranktable(wrxid,wryid-1)

    pxyne = ranktable(wrxid+1,wryid+1)
    pxyse = ranktable(wrxid+1,wryid-1)
    pxynw = ranktable(wrxid-1,wryid+1)
    pxysw = ranktable(wrxid-1,wryid-1)

    call mpi_isend(var(1,n2-3,3),1,xstride, pxfwd, 130, &
         MPI_COMM_WORLD, req(1), ierror)
    call mpi_isend(var(1,3,3), 1, xstride, pxback, 140, &
         MPI_COMM_WORLD, req(2), ierror)
    call mpi_irecv(var(1,n2-1,3), 1, xstride, pxfwd, 140, &
         MPI_COMM_WORLD, req(3), ierror)
    call mpi_irecv(var(1,1,3), 1, xstride, pxback, 130, &
         MPI_COMM_WORLD, req(4), ierror)

    call mpi_isend(var(1,3,3), 1, ystride, pyback, 110, &
         MPI_COMM_WORLD, req(5), ierror)
    call mpi_irecv(var(1,3,n3-1), 1, ystride, pyfwd, 110, &
         MPI_COMM_WORLD, req(6), ierror)
    call mpi_isend(var(1,3,n3-3), 1, ystride, pyfwd, 120, &
         MPI_COMM_WORLD, req(7), ierror)
    call mpi_irecv(var(1,3,1), 1, ystride, pyback, 120, &
         MPI_COMM_WORLD, req(8), ierror)

    call mpi_isend(var(1,n2-3,n3-3), 1, xystride, pxyne, 150, &
         MPI_COMM_WORLD, req(9), ierror)
    call mpi_irecv(var(1,1,1), 1, xystride, pxysw, 150, &
         MPI_COMM_WORLD, req(10), ierror)
    call mpi_isend(var(1,n2-3,3), 1, xystride, pxyse, 160, &
         MPI_COMM_WORLD, req(11), ierror)
    call mpi_irecv(var(1,1,n3-1), 1, xystride, pxynw, 160, &
         MPI_COMM_WORLD, req(12), ierror)

    call mpi_isend(var(1,3,n3-3), 1, xystride, pxynw, 170, &
         MPI_COMM_WORLD, req(13), ierror)
    call mpi_irecv(var(1,n2-1,1), 1, xystride, pxyse, 170, &
         MPI_COMM_WORLD, req(14), ierror)
    call mpi_isend(var(1,3,3), 1, xystride, pxysw, 180, &
         MPI_COMM_WORLD, req(15), ierror)
    call mpi_irecv(var(1,n2-1,n3-1), 1, xystride, pxyne, 180, &
         MPI_COMM_WORLD, req(16), ierror)

  end subroutine cyclics
  !
  !
  ! ---------------------------------------------------------------------
  ! subroutine cyclicc: comits excahnging cyclic boundary conditions
  subroutine cyclicc(n1,n2,n3,var,req)

    integer :: ierror, stats(MPI_STATUS_SIZE,16)
    integer :: req(16),n1,n2,n3
    real :: var(n1,n2,n3)

    call mpi_waitall(16,req,stats,ierror)

  end subroutine cyclicc

  subroutine appl_abort(apperr)

    integer apperr,ierr
    call mpi_abort(MPI_COMM_WORLD,apperr,ierr)

  end subroutine appl_abort

  subroutine appl_finalize(ierr)

    integer ierr
    call mpi_finalize(ierr)

  end subroutine appl_finalize

  subroutine xshuffle(a,atmp,nx,ny,nz,isign)

    integer, intent(in):: nx,ny,nz,isign
    complex, intent(inout):: a(nx,ny,nz),atmp((nx+1)*(ny+1)*(nz+1))
    integer ierr,ll,i,j,k

    if(isign .eq. 1) then
       if(nxprocs .ne. 1)then
          call mpi_alltoallw( a,xcount(0:,1) , xdisp(0:,1), xtype(0:,1), atmp, &
               xcount(0:,2), xdisp(0:,2),xtype(0:,2),xcomm,ierr)
       else
          ll=0
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   ll=ll+1
                   atmp(ll)=a(i,j,k)
                enddo
             enddo
          enddo

       endif
    else
       if(nxprocs .ne. 1)then
          call mpi_alltoallw(atmp,xcount(0:,2),xdisp(0:,2),xtype(0:,2),a, &
               xcount(0:,1), xdisp(0:,1),xtype(0:,1),xcomm,ierr)
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
    endif

  end subroutine xshuffle

  subroutine yshuffle(a,atmp,nx,ny,nz,isign)

    integer, intent(in):: nx,ny,nz,isign
    complex, intent(inout):: a(ny,nx,nz),atmp((nx+1)*(ny+1)*(nz+1))
    integer ierr,ll,i,j,k

    if(isign .eq. 1) then
       if(nyprocs .ne. 1)then
          call mpi_alltoallw( a,ycount(0:,1),ydisp(0:,1),ytype(0:,1),atmp, &
               ycount(0:,2),ydisp(0:,2),ytype(0:,2),ycomm,ierr)
       else
          ll=0
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   ll=ll+1
                   atmp(ll)=a(i,j,k)
                enddo
             enddo
          enddo
       endif
    else
       if(nyprocs .ne. 1)then
          call mpi_alltoallw(atmp,ycount(0:,2),ydisp(0:,2),ytype(0:,2),a, &
               ycount(0:,1),ydisp(0:,1),ytype(0:,1),ycomm,ierr)
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
    endif

  end subroutine yshuffle
  !
  !---------------------------------------------------------------------------
  ! get maximum across processors
  !
  subroutine double_scalar_par_max(xxl,xxg)

    real(kind=8), intent(out) :: xxg
    real(kind=8), intent(in) :: xxl
    integer:: mpiop,ierror


    call mpi_allreduce(xxl,xxg,1,MPI_DOUBLE_PRECISION, MPI_MAX, &
         MPI_COMM_WORLD, ierror)

  end subroutine double_scalar_par_max


  subroutine double_scalar_par_sum(xxl,xxg)

    real(kind=8), intent(out) :: xxg
    real(kind=8), intent(in) :: xxl
    integer:: mpiop,ierror


    call mpi_allreduce(xxl,xxg,1,MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierror)

  end subroutine double_scalar_par_sum

  subroutine double_array_par_sum(xxl,xxg,n)

    integer, intent(in)::n
    real(kind=8), intent(out) :: xxg(n)
    real(kind=8), intent(in) :: xxl(n)
    integer:: mpiop,ierror


    call mpi_allreduce(xxl,xxg,n,MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierror)

  end subroutine double_array_par_sum


end module mpi_interface
