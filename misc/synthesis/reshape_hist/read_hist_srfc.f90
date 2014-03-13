module readhist_srfc
contains
  subroutine read_hist_srfc(nxp1, nyp1, nx1, ny1, nxt, nyt,&
       a_tsoil, a_phiw, a_tskin, a_qskin, a_wl)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp1, nyp1
    integer :: unit
    logical :: UNITOK, UNITOP

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension (4,nxt,nyt) :: &
          a_tsoil, &
          a_phiw

    real, dimension (nxt,nyt) :: &
         a_tskin, &
         a_qskin, &
         a_Wl

  ! list of variables of the small grid

    integer :: nx1, ny1

    real, dimension (4,nx1,ny1) :: &
       a_tsoil_l, &
       a_phiw_l

    real, dimension (nx1,ny1) :: &
       a_tskin_l, &
       a_qskin_l, &
       a_Wl_l

    ! Loop over all processors

    do xid = 1, nxp1
       do yid = 1, nyp1

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          !
          ! open input file.
          !
          unit = 3000 + wryid-wryid/nyp1 + nyp1*wrxid

          inquire (unit=unit,exist=UNITOK,opened=UNITOP)
          if (.not.UNITOK) then
             print*,'unit does not exist', unit
             stop
          end if
          if (.not.UNITOP) then
             print*,'unit not opened'
             stop
          end if
          read(unit) a_tsoil_l
          read(unit) a_phiw_l
          read(unit) a_tskin_l
          read(unit) a_qskin_l
          read(unit) a_wl_l

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

                do k = 1, 4
                   a_tsoil(k,i,j) = a_tsoil_l(k,ii,jj)
                   a_phiw(k,i,j)  = a_phiw_l(k,ii,jj)
                end do
                a_tskin(i,j) = a_tskin_l(ii,jj)
                a_qskin(i,j) = a_qskin_l(ii,jj)
                a_wl(i,j)    = a_wl_l(ii,jj)
             end do
          end do
       end do ! loop over processors
    end do ! loop over processors
  end subroutine read_hist_srfc
end module readhist_srfc
