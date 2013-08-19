module writehist_1
contains
  subroutine write_hist_1(time, hname, nxp2, nyp2, nx2, ny2, nxt, nyt, nz, nscl,&
       umean, vmean, th00, level, isfctyp, lwaterbudget, iradtyp, xt, xm, yt,   &
       ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf, seed_arr, nseed, dt,      &
       a_ustar, a_tstar, a_rstar, a_pexnr, nxp1, nyp1)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp2, nyp2, nz, nscl, nv2, nxp1, nyp1

    character(len=40)             :: filename
    real, intent(in)              :: time
    integer                       :: unit, seedct

    character (len=40), intent(in) :: hname
    integer :: n, nseed, nxpx, nypx, nzpx, iradtyp
    integer, dimension(:,:), allocatable :: seed_arr
    logical :: exans, lwaterbudget
    real    :: umean, vmean, th00, dt, psrf
    integer :: isfctyp, level, nsmp

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension(nxt) ::  xt, xm
    real, dimension(nyt) ::  yt, ym
    real, dimension(nz)  ::  zt, zm, u0, v0, pi0, pi1, th0, dn0, rt0 

    real, dimension (nxt,nyt) :: &
         a_ustar, &
         a_tstar, &
         a_rstar

    real, dimension (nz,nxt,nyt), intent(in) :: &
          a_pexnr

  ! list of variables of the small grid

    integer :: nx2, ny2

    real, dimension(nx2) ::  xt_l, xm_l
    real, dimension(ny2) ::  yt_l, ym_l

    real, dimension (nx2,ny2) :: &
       a_ustar_l, &
       a_tstar_l, &
       a_rstar_l

    real, dimension (nz,nx2,ny2) :: & 
       a_pexnr_l

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
             ii = i-istart+1
             xt_l(ii) = xt(i)
             xm_l(ii) = xm(i)
          end do

          do j = jstart, jend
             jj = j-jstart+1
             yt_l(jj) = yt(j)
             ym_l(jj) = ym(j)
          end do

          do i = istart, iend
             do j = jstart, jend
                ii = i-istart+1
                jj = j-jstart+1

                a_ustar_l(ii,jj)=a_ustar(i,j)
                a_tstar_l(ii,jj)=a_tstar(i,j)
                a_rstar_l(ii,jj)=a_rstar(i,j)
                do k = 1, nz
                   a_pexnr_l(k,ii,jj) = a_pexnr(k,i,j)
                end do
             end do
          end do

          !
          ! open output file.
          !

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          unit = 10 + wryid-wryid/nyp2 + nyp2*wrxid
          seedct = int(real(nyp1*wryid)/real(nyp2)) + &
               int(real(nyp1*nyp2)/real(nyp2))*int(real(nxp1*wrxid)/real(nxp2)) +1

          write(filename,'(i4.4,a1,i4.4)') wrxid,'_',wryid
          filename = './out/'//trim(filename)//'.'//trim(hname)
          inquire(file=trim(filename),exist=exans)
          if (exans) then
             print *,'ABORTING: History file', trim(filename),' exists already'
             stop
          else
             open (unit,file=trim(filename),status='new',form='unformatted')
             write (unit) time,th00,umean,vmean,dt,level,iradtyp,nz,nx2,ny2,nscl
             write (unit) nseed
             write (unit) seed_arr(:,seedct)

             write (unit) xt_l, xm_l, yt_l, ym_l, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf
             write (unit) a_ustar_l, a_tstar_l, a_rstar_l
             write (unit) a_pexnr_l

             print "('Wrote history to: ',A60)",filename
             close(unit)
          end if
       end do
    end do

  end subroutine write_hist_1
end module writehist_1
