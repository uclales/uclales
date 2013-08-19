module writehist_3
contains
  subroutine write_hist_3(hname, nxp2, nyp2, nx2, ny2, nxt, nyt, nz, level, &
       nv2, nsmp, svctr, prc_acc, rev_acc, cnd_acc, cev_acc)

    implicit none

    integer :: wrxid, wryid, xid, yid, i, j, ii, jj, k, istart, iend, jstart, jend
    integer :: nxp2, nyp2, nz, nscl, nv2, unit

    character(len=40)              :: filename
    character (len=40), intent(in) :: hname
    logical :: exans

    integer :: level, nsmp

    ! list of variables of the large grid

    integer :: nxt, nyt

    real, dimension (nxt,nyt), optional, intent(in) :: &
         cnd_acc, &  ! accumulated condensation  [kg/m2] (diagnostic for 2D output)
         cev_acc     ! accumulated evaporation of cloud water [kg/m2] (diagnostic for 2D output)

    real, dimension (nxt,nyt), intent(in) :: &
         rev_acc     ! accumulated evaporation of rainwater   [kg/m2] (diagnostic for 2D output)
    
    real, dimension (nxt,nyt,0) :: &
         prc_acc

    real, dimension (:,:), allocatable :: &
         svctr
    
  ! list of variables of the small grid

    integer :: nx2, ny2

    real, dimension (nx2,ny2) :: &
       cnd_acc_l, &  ! accumulated condensation  [kg/m2] (diagnostic for 2D output)
       cev_acc_l     ! accumulated evaporation of cloud water [kg/m2] (diagnostic for 2D output)

    real, dimension (nx2,ny2) :: &
       rev_acc_l     ! accumulated evaporation of rainwater   [kg/m2] (diagnostic for 2D output)

    real, dimension (nx2,ny2,0) :: &
       prc_acc_l

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

                if(level>=3) then
                   prc_acc_l(ii,jj,:) = prc_acc(i,j,:)
                   rev_acc_l(ii,jj) = rev_acc(i,j)
                end if
                if(present(cnd_acc)) then
                   cnd_acc_l(ii,jj) = cnd_acc(i,j)
                   cev_acc_l(ii,jj) = cev_acc(i,j)
                end if
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
             if(level>=3) then
                write(unit) prc_acc_l, rev_acc_l
             end if
             if(present(cnd_acc)) then
                write(unit) cnd_acc_l, cev_acc_l
             end if
             write(unit) nv2, nsmp
             write(unit) svctr

             close(unit)

             print "('Closed file: ',A60)",filename
          end if
       end do
    end do

  end subroutine write_hist_3
end module writehist_3
