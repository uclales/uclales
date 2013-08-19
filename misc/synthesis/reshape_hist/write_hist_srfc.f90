module writehist_srfc
contains
  subroutine write_hist_srfc(nxp2, nyp2, nx2, ny2, nxt, nyt, hname,&
       a_tsoil, a_phiw, a_tskin, a_qskin, a_wl)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp2, nyp2, unit

    character(len=40)             :: filename
    character(len=40), intent(in) :: hname
    logical :: exans

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension (4,nxt,nyt), intent(in) :: &
          a_tsoil, &
          a_phiw

    real, dimension (nxt,nyt), intent(in) :: &
         a_tskin, &
         a_qskin, &
         a_Wl

  ! list of variables of the small grid

    integer :: nx2, ny2

    real, dimension(nx2) ::  xt_l, xm_l
    real, dimension(ny2) ::  yt_l, ym_l

    real, dimension (4,nx2,ny2) :: &
       a_tsoil_l, &
       a_phiw_l

    real, dimension (nx2,ny2) :: &
       a_tskin_l, &
       a_qskin_l, &
       a_Wl_l

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

                do k = 1, 4
                   a_tsoil_l(k,ii,jj) = a_tsoil(k,i,j)
                   a_phiw_l(k,ii,jj)  = a_phiw(k,i,j)
                end do
                a_tskin_l(ii,jj) = a_tskin(i,j)
                a_qskin_l(ii,jj) = a_qskin(i,j)
                a_wl_l(ii,jj)    = a_wl(i,j)
             end do
          end do
          !
          ! open output file.
          !

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          unit = 10 + wryid-wryid/nyp2 + nyp2*wrxid
          write(filename,'(i4.4,a1,i4.4)') wrxid,'_',wryid
          filename = './out/'//trim(filename)//'.'//trim(hname)
          inquire(file=trim(filename),exist=exans)
          if (.not.exans) then
             print *,'ABORTING: History file', trim(filename),' does not exist'
             stop
          else
             open (unit,file=trim(filename),status='old',form='unformatted',position='append')!,convert='BIG_ENDIAN')

             write(unit) a_tsoil_l
             write(unit) a_phiw_l
             write(unit) a_tskin_l
             write(unit) a_qskin_l
             write(unit) a_wl_l

             close(unit)
          end if
       end do
    end do

  end subroutine write_hist_srfc
end module writehist_srfc
