module debugg
    character*4 str1,str2
    integer ierror,n3d,n2d,n1dx,n1dy,n1dz
contains
  subroutine debuggg
    use grid, only : level,  a_tt, dzt, nxp, nyp, &
         a_wp, dtl, zt, zm, nzp, dn0, u0, v0, th0, a_tp,   &
         a_up, a_vp, a_pexnr, a_theta, a_rp, a_uc, a_vc,   &
         a_wc, write_hist, write_anal, close_anal,         &
         write_slice, close_slice, &
         a_rpp,a_rpt,a_npp,a_npt,    &
         a_rc, a_rv,a_press, a_pexnr,a_ustar, a_tstar, a_rstar,   &
         uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc, precip,xt,yt,zt,xm,ym,zm,   &
         dn0, dzt,dzm,a_ut,a_vt,a_wt
    use mpi_interface, only : init_mpi, define_decomp, nxpg, nypg, pecount, &
                              xoffset, yoffset, myid, init_alltoall_reorder, &
                              wrxid, wryid, nxprocs,nyprocs
    implicit none


!    write(str1,'(i4)')wrxid
!    write(str2,'(i4)')wryid
!    if(wrxid.ge.0 .and. wrxid .le. 9 )then
!      str1(1:3)='000'
!    elseif(wrxid.gt.9 .and. wrxid .le. 99 )then
!      str1(1:2)='00'
!    elseif(wrxid.gt.99 .and. wrxid .le. 999 )then
!      str1(1:1)='0'
!    endif
!
!    if(wryid.ge.0 .and. wryid .le. 9 )then
!      str2(1:3)='000'
!    elseif(wryid.gt.9 .and. wryid .le. 99 )then
!      str2(1:2)='00'
!    elseif(wryid.gt.99 .and. wryid .le. 999 )then
!      str2(1:1)='0'
!    endif
!
!    open(121212, file='s2d'//str1//str2//'.dat',form='formatted')
!    open(131313, file='single2d'//str1//str2//'.dat',form='formatted')
!    open(141414, file='singlesimple2d'//str1//str2//'.dat',form='formatted')
!    open(151515, file='twodim2d'//str1//str2//'.dat',form='formatted')
!!ctvs    write(121212,*)xoffset(wrxid),yoffset(wryid)
!!    write(121212,1020)th0

!ctvs 3-d 11 variables
    n1dx=2
    n1dy=2
    n1dz=5
    n3d=0
    if(allocated(a_up))then
    n3d=n3d+1
    endif
    if(allocated(a_vp))then
    n3d=n3d+1
    endif
    if(allocated(a_wp))then
    n3d=n3d+1
    endif
    if(allocated(a_uc))then
    n3d=n3d+1
    endif
    if(allocated(a_vc))then
    n3d=n3d+1
    endif
    if(allocated(a_wc))then
    n3d=n3d+1
    endif
    if(allocated(a_ut))then
    n3d=n3d+1
    endif
    if(allocated(a_vt))then
    n3d=n3d+1
    endif
    if(allocated(a_wt))then
    n3d=n3d+1
    endif
    if(allocated(a_theta))then
    n3d=n3d+1
    endif
    if(allocated(a_rc))then
    n3d=n3d+1
    endif
    if(allocated(a_rv))then
    n3d=n3d+1
    endif
    if(allocated(a_press))then
    n3d=n3d+1
    endif
    if(allocated(a_pexnr))then
    n3d=n3d+1
    endif
    n2d=0
    if(allocated(a_ustar))then
    n2d=n2d+1
    endif
    if(allocated(a_tstar))then
    n2d=n2d+1
    endif
    if(allocated(a_rstar))then
    n2d=n2d+1
    endif
    if(allocated(uw_sfc))then
    n2d=n2d+1
    endif
    if(allocated(vw_sfc))then
    n2d=n2d+1
    endif
    if(allocated(ww_sfc))then
    n2d=n2d+1
    endif
    if(allocated(wt_sfc))then
    n2d=n2d+1
    endif
    if(allocated(wq_sfc))then
    n2d=n2d+1
    endif
    if(allocated(precip))then
    n2d=n2d+1
    endif
    
    write(121212,1212)n3d,n2d,n1dx,n1dy,n1dz
 1212 format (1x,5I8)

    if(allocated(a_up))then
    write(121212,*)'a_up4444'
    write(121212,1020)a_up(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_vp))then
    write(121212,*)'a_vp4444'
    write(121212,1020)a_vp(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_wp))then
    write(121212,*)'a_wp4444'
    write(121212,1020)a_wp(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_uc))then
    write(121212,*)'a_uc4444'
    write(121212,1020)a_uc(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_vc))then
    write(121212,*)'a_vc4444'
    write(121212,1020)a_vc(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_wc))then
    write(121212,*)'a_wc4444'
    write(121212,1020)a_wc(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_ut))then
    write(121212,*)'a_ut4444'
    write(121212,1020)a_ut(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_vt))then
    write(121212,*)'a_vt4444'
    write(121212,1020)a_vt(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_wt))then
    write(121212,*)'a_wt4444'
    write(121212,1020)a_wt(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_theta))then
    write(121212,*)'a_theta4444'
    write(121212,1020)a_theta(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_rc))then
    write(121212,*)'a_rc4444'
    write(121212,1020)a_rc(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_rv))then
    write(121212,*)'a_rv4444'
    write(121212,1020)a_rv(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_press))then
    write(121212,*)'a_press4444'
    write(121212,1020)a_press(1:nzp, 3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_pexnr))then
    write(121212,*)'a_pexnr4444'
    write(121212,1020)a_pexnr(1:nzp, 3:nxp-2, 3:nyp-2)
    endif

!ctvs 2-d 9 variables
    if(allocated(a_ustar))then
    write(121212,*)'a_ustar4444'
    write(121212,1020)a_ustar(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_tstar))then
    write(121212,*)'a_tstar4444'
    write(121212,1020)a_tstar(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(a_rstar))then
    write(121212,*)'a_rstar4444'
    write(121212,1020)a_rstar(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(uw_sfc))then
    write(121212,*)'uw_sfc4444'
    write(121212,1020)uw_sfc(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(vw_sfc))then
    write(121212,*)'vw_sfc4444'
    write(121212,1020)vw_sfc(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(ww_sfc))then
    write(121212,*)'ww_sfc4444'
    write(121212,1020)ww_sfc(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(wt_sfc))then
    write(121212,*)'wt_sfc4444'
    write(121212,1020)wt_sfc(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(wq_sfc))then
    write(121212,*)'wq_sfc4444'
    write(121212,1020)wq_sfc(3:nxp-2, 3:nyp-2)
    endif
    if(allocated(precip))then
    write(121212,*)'precip4444'
    write(121212,1020)precip(3:nxp-2, 3:nyp-2)
    endif

!ctvs 1-d nxp 2 variables
    write(121212,*)'xt444444'
    write(121212,1020)xt(3:nxp-2)
    write(121212,*)'xm444444'
    write(121212,1020)xm(3:nxp-2)

!ctvs 1-d nyp 2 variables
    write(121212,*)'yt444444'
    write(121212,1020)yt(3:nyp-2)
    write(121212,*)'ym444444'
    write(121212,1020)ym(3:nyp-2)

!ctvs 1-d nzp 5 variables
    write(121212,*)'dn0444444'
    write(121212,1020)dn0(1:nzp)
    write(121212,*)'zt444444'
    write(121212,1020)zt(1:nzp)
    write(121212,*)'zm444444'
    write(121212,1020)zm(1:nzp)
    write(121212,*)'dzt444444'
    write(121212,1020)dzt(1:nzp)
    write(121212,*)'dzm444444'
    write(121212,1020)dzm(1:nzp)
    
 1020 format(1x,5E14.4)
    return
!ctvs    call mpi_finalize(ierror)
!ctvs    stop
end subroutine debuggg

    subroutine openfiles
    use grid, only : level,  a_tt, dzt, nxp, nyp, &
         a_wp, dtl, zt, zm, nzp, dn0, u0, v0, th0, a_tp,   &
         a_up, a_vp, a_pexnr, a_theta, a_rp, a_uc, a_vc,   &
         a_wc, write_hist, write_anal, close_anal,         &
         write_slice, close_slice, &
         a_rpp,a_rpt,a_npp,a_npt,    &
         a_rc, a_rv,a_press, a_pexnr,a_ustar, a_tstar, a_rstar,   &
         uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc, precip,xt,yt,zt,xm,ym,zm,   &
         dn0, dzt,dzm,a_ut,a_vt,a_wt
    use mpi_interface, only : init_mpi, define_decomp, nxpg, nypg, pecount, &
                              xoffset, yoffset, myid, init_alltoall_reorder, &
                              wrxid, wryid, nxprocs,nyprocs
    implicit none


    write(str1,'(i4)')wrxid
    write(str2,'(i4)')wryid
    if(wrxid.ge.0 .and. wrxid .le. 9 )then
      str1(1:3)='000'
    elseif(wrxid.gt.9 .and. wrxid .le. 99 )then
      str1(1:2)='00'
    elseif(wrxid.gt.99 .and. wrxid .le. 999 )then
      str1(1:1)='0'
    endif

    if(wryid.ge.0 .and. wryid .le. 9 )then
      str2(1:3)='000'
    elseif(wryid.gt.9 .and. wryid .le. 99 )then
      str2(1:2)='00'
    elseif(wryid.gt.99 .and. wryid .le. 999 )then
      str2(1:1)='0'
    endif

    open(121212, file='s2d'//str1//str2//'.dat',form='formatted')
    open(131313, file='single2d'//str1//str2//'.dat',form='formatted')
    open(141414, file='singlesimple2d'//str1//str2//'.dat',form='formatted')
    open(151515, file='twodim2d'//str1//str2//'.dat',form='formatted')
!ctvs    write(121212,*)xoffset(wrxid),yoffset(wryid)
!    write(121212,1020)th0
    return
    end subroutine openfiles

    subroutine closefiles
    implicit none
!ctvs    call flush(121212)
!ctvs    call flush(131313)
!ctvs    call flush(141414)
!ctvs    call flush(151515)
    close(121212)
    close(131313)
    close(141414)
    close(151515)
    return
    end subroutine closefiles
end module debugg
