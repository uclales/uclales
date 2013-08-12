module readhist_2
contains
  subroutine read_hist_2(nxp1, nyp1, nx1, ny1, nxt, nyt, nz, a_xp)
    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp1, nyp1, nz
    integer :: unit
    logical :: UNITOK, UNITOP

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension (nz,nxt,nyt):: &
         a_xp

  ! list of variables of the small grid

    integer :: nx1, ny1

    real, dimension (nz,nx1,ny1):: a_xp_l

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

          read (unit) a_xp_l


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
                
                a_xp(:,i,j) = a_xp_l(:,ii,jj)

             end do
          end do
       end do ! loop over processors
    end do ! loop over processors
  end subroutine read_hist_2
end module readhist_2
