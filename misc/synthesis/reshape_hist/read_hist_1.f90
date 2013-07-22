module readhist_1
contains
  subroutine read_hist_1(time, hname, nxp1, nyp1, nx1, ny1, nxt, nyt, nz, nscl,&
       umean, vmean, th00, level, isfctyp, lwaterbudget, iradtyp, xt, xm, yt,  &
       ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf, seed_arr, nseed, dt,     &
       a_ustar, a_tstar, a_rstar, a_pexnr)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp1, nyp1, nz, nscl, nv2

    character(len=40) :: filename
    real              :: time
    integer           :: unit, seedct

    character (len=40) :: hname
    integer :: n, nseed, nxpx, nypx, nzpx, nsclx, iradtyp, lvlx
    integer, dimension(:,:), allocatable :: seed_arr
    logical :: exans
    real    :: umean, vmean, th00, dt, psrf
    integer :: isfctyp, level, nsmp
    logical :: lwaterbudget

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension(nxt) ::  xt, xm
    real, dimension(nyt) ::  yt, ym
    real, dimension(nz)   ::  zt, zm, u0, v0, pi0, pi1, th0, dn0, rt0 

    real, dimension (nxt,nyt) :: &
         a_ustar, &
         a_tstar, &
         a_rstar

    real, dimension (nz,nxt,nyt) :: &
         a_pexnr

  ! list of variables of the small grid

    integer :: nx1, ny1

    real, dimension(nx1) ::  xt_l, xm_l
    real, dimension(ny1) ::  yt_l, ym_l

    real, dimension (nx1,ny1) :: &
       a_ustar_l, &
       a_tstar_l, &
       a_rstar_l

    real, dimension (nz,nx1,ny1) :: & 
       a_pexnr_l

    ! Loop over all processors

    do xid = 1, nxp1
       do yid = 1, nyp1

          wrxid = xid -1 ! zero-based
          wryid = yid -1

          !
          ! open input file.
          !
          unit = 300 + wryid-wryid/nyp1 + nyp1*wrxid
          seedct = wryid + nyp1*wrxid +1

          write(filename,'(i4.4,a1,i4.4)') wrxid,'_',wryid
          filename = trim(filename)//'.'//trim(hname)
          inquire(file=trim(filename),exist=exans)
          if (.not.exans) then
             print *,'ABORTING: History file', trim(filename),' not found'
             stop
          else
             open (unit,file=trim(filename),status='old',form='unformatted')
             read (unit) time,th00,umean,vmean,dt,lvlx,iradtyp,nzpx,nxpx,nypx,nsclx
             read (unit) nseed
             ! the seed is redistributed to the new fields. If nxp2>nxp1 then the seed from 
             ! adjacent domains is taken, if nxp2<nyp1 some seeds inbetween are dropped. 
             if (xid==1.and.yid==1) then
                allocate(seed_arr(nseed, nxp1*nyp1))
             end if
             read(unit) seed_arr(:,seedct)
             if (nxpx /= nx1 .or. nypx /= ny1 .or. nzpx /= nz)  then
                print *, 'WRONG GRID !', nx1, ny1, nz, nxpx, nypx, nzpx
                stop
             end if
             if (lvlx /= level)  then
                print *, 'WRONG LEVEL GIVEN', level, lvlx
                stop
             end if
             if (nsclx /= nscl)  then
                print *, 'NUMBER OF SCALARS WRONG', nscl, nsclx
                stop
             end if

             read (unit) xt_l, xm_l, yt_l, ym_l, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf
             read (unit) a_ustar_l, a_tstar_l, a_rstar_l
             read (unit) a_pexnr_l

             print "('History read from: ',A60)",filename

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
                if (xid == 1) then
                   ii = i-istart+1
                else
                   ii = i-istart+3
                end if
                xt(i) = xt_l(ii)
                xm(i) = xm_l(ii)
             end do

             do j = jstart, jend
                if (yid == 1) then
                   jj = j-jstart+1
                else
                   jj = j-jstart+3
                end if
                yt(j) = yt_l(jj)
                ym(j) = ym_l(jj)
             end do

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

                   a_ustar(i,j)=a_ustar_l(ii,jj)
                   a_tstar(i,j)=a_tstar_l(ii,jj)
                   a_rstar(i,j)=a_rstar_l(ii,jj)
                   do k = 1, nz
                      a_pexnr(k,i,j) = a_pexnr_l(k,ii,jj)
                   end do
                end do
             end do
          end if ! file existst
       end do ! loop over processors
    end do ! loop over processors
  end subroutine read_hist_1
end module readhist_1
