module readhist_3
contains
  subroutine read_hist_3(hname, nxp1, nyp1, nx1, ny1, nxt, nyt, nz, level, &
       nv2, nsmp, svctr, prc_acc, rev_acc, cnd_acc, cev_acc)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp1, nyp1, nz, nscl, nv2, unit
    logical :: UNITOK, UNITOP

    character(len=40)  :: filename
    character (len=40) :: hname

    integer :: level, nsmp

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, optional, dimension (nxt,nyt) :: &
         cnd_acc, &  ! accumulated condensation  [kg/m2] (diagnostic for 2D output)
         cev_acc     ! accumulated evaporation of cloud water [kg/m2] (diagnostic for 2D output)

    real, dimension (nxt,nyt) :: &
         rev_acc     ! accumulated evaporation of rainwater   [kg/m2] (diagnostic for 2D output)

    real, dimension (nxt,nyt,0) :: &
         prc_acc

    real, dimension (:,:), allocatable :: &
         svctr
    
  ! list of variables of the small grid

    integer :: nx1, ny1

    real, dimension (nx1,ny1) :: &
       cnd_acc_l, &  ! accumulated condensation  [kg/m2] (diagnostic for 2D output)
       cev_acc_l     ! accumulated evaporation of cloud water [kg/m2] (diagnostic for 2D output)

    real, dimension (nx1,ny1) :: &
       rev_acc_l     ! accumulated evaporation of rainwater   [kg/m2] (diagnostic for 2D output)

    real, dimension (nx1,ny1,0) :: &
       prc_acc_l

    ! Loop over all processors

    do xid = 1, nxp1
       do yid = 1, nyp1

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          unit = 300 + wryid-wryid/nyp1 + nyp1*wrxid
          inquire (unit=unit,exist=UNITOK,opened=UNITOP)
          if (.not.UNITOK) then
             print*,'unit does not exist', unit
             stop
          end if
          if (.not.UNITOP) then
             print*,'unit not opened'
             stop
          end if
          !
          ! open input file.
          !

          write(filename,'(i4.4,a1,i4.4)') wrxid,'_',wryid
          filename = trim(filename)//'.'//trim(hname)

          if(level>=3) then
             read(unit) prc_acc_l, rev_acc_l
          end if
          if(present(cnd_acc)) then
             read(unit) cnd_acc_l, cev_acc_l
          end if
          read(unit) nv2, nsmp
          if (xid==1.and.yid==1) then
             allocate (svctr(nz,nv2))
          end if
          read(unit) svctr

          close(unit)

          print "('Closed file: ',A60)",filename

          ! put the data of the small local array into the correct position of the large array
          ! take into account the overlap of the subdomains

          if (xid == 1) then
             istart = 1 
          else
             istart = (xid-1)*(nx1-4) +2 +1
          end if
          if (xid == nxp1) then
             iend = nxt
          else
             iend  = (xid)*(nx1-4) +2
          end if
          if (yid == 1) then
             jstart = 1
          else
             jstart = (yid-1)*(ny1-4) + 2 +1
          end if
          if (yid == nyp1) then
             jend = nyt
          else
             jend   = (yid)*(ny1-4) +2 
          end if

          do i = istart, iend
             do j = jstart, jend
                if (xid == 1) then
                   ii = i-istart+1
                else
                   ii = i-istart+3
                end if
                if (yid == 1) then
                   jj = j-jstart+1
                else
                   jj = j-jstart+3
                end if

                if(level>=3) then
                   prc_acc(i,j,:) = prc_acc_l(ii,jj,:)
                   rev_acc(i,j) = rev_acc_l(ii,jj)
                end if
                if(present(cnd_acc)) then
                   cnd_acc(i,j) = cnd_acc_l(ii,jj)
                   cev_acc(i,j) = cev_acc_l(ii,jj)
                end if
             end do
          end do
       end do ! loop over processors
    end do ! loop over processors
  end subroutine read_hist_3
end module readhist_3
