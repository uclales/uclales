module readhist_srfc_rad
contains
  subroutine read_hist_srfc_rad(nxp1, nyp1, nx1, ny1, nxt, nyt, a_flx)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp1, nyp1
    integer :: unit
    logical :: UNITOK, UNITOP

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension (100,nxt,nyt) :: & 
          a_flx

  ! list of variables of the small grid

    integer :: nx1, ny1

    real, dimension (100,nx1,ny1) :: &
       a_flx_l

    ! Loop over all processors

    do xid = 1, nxp1
       do yid = 1, nyp1

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          !
          ! open input file.
          !
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
          read(unit) a_flx_l

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

                do k = 1, 100
                   a_flx(k,i,j) = a_flx_l(k,ii,jj)
                end do
             end do
          end do
       end do ! loop over processors
    end do ! loop over processors
  end subroutine read_hist_srfc_rad
end module readhist_srfc_rad
