module writehist_srfc_rad
contains
  subroutine write_hist_srfc_rad(nxp2, nyp2, nx2, ny2, nxt, nyt, a_flx)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp2, nyp2, unit

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension (100,nxt,nyt), intent(in) :: & 
          a_flx

  ! list of variables of the small grid

    integer :: nx2, ny2

    real, dimension(nx2) ::  xt_l, xm_l
    real, dimension(ny2) ::  yt_l, ym_l

    real, dimension (100,nx2,ny2) :: &
       a_flx_l

    ! Loop over all the processors

    do xid = 1, nxp2
       do yid = 1, nyp2
          
          ! get the data from the correct position of the large array for the small local array 

          if (xid == 1) then
             istart = 1 
          else
             istart = (xid-1)*(nx2-4) +2 -2 +1
          end if
          if (xid == nxp2) then
             iend = nxt
          else
             iend  = xid*(nx2-4) +4
          end if
          if (yid == 1) then
             jstart = 1
          else
             jstart = (yid-1)*(ny2-4) +2 -2 +1
          end if
          if (yid == nyp2) then
             jend = nyt
          else
             jend  = yid*(ny2-4) +4
          end if

          do i = istart, iend
             do j = jstart, jend
                ii = i-istart+1
                jj = j-jstart+1

                do k = 1, 100
                   a_flx_l(k,ii,jj) = a_flx(k,i,j)
                end do
             end do
          end do
          !
          ! open output file.
          !

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          unit = 10 + wryid-wryid/nyp2 + nyp2*wrxid

          write(unit) a_flx_l

       end do
    end do

  end subroutine write_hist_srfc_rad
end module writehist_srfc_rad
