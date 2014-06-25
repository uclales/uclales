!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module sgsm

  use stat, only : sflg, updtst, acc_tend, sgsflxs, sgs_vel
  use util, only : tridiff
  implicit none
!
! setting the prandtl number to a value less than zero enforces an exponential
! decay in its value as a function of height from the surface with a scale
! height of 100m.  This is the default
!

  real, parameter     :: tkemin=1.e-20
  real :: csx = 0.23
  real :: prndtl = 0.3333333333
  real :: clouddiff = -1.
  real, allocatable, dimension (:,:) :: sxy1, sxy2, sxy3, sxz1, sxz2, sxz3    &
       , sxz4, sxz5, sxz6  ,szx1, szx2, szx3, szx4, szx5
  real, allocatable, dimension (:,:,:) :: szxy 
  real, allocatable, dimension (:)   :: sz1, sz2, sz3, sz4, sz5, sz6, sz7, sz8

  integer :: k, i, j, indh, req(16)
  real    :: dti, dfact
  logical, save :: Initialized = .false.

contains 

  subroutine diffuse_init(n1,n2,n3)

    integer, intent (in) :: n1, n2, n3

    allocate(sxy1(n2,n3), sxy2(n2,n3), sxy3(n2,n3), sxz1(n2,n1), sxz2(n2,n1))
    allocate(szxy(n1,n2,n3)) 
    allocate(sxz3(n2,n1), sxz4(n2,n1), sxz5(n2,n1), sxz6(n2,n1))
    allocate(szx1(n1,n2), szx2(n1,n2), szx3(n1,n2), szx4(n1,n2), szx5(n1,n2))
    allocate(sz1(n1),sz2(n1),sz3(n1),sz4(n1),sz5(n1),sz6(n1),sz7(n1),sz8(n1))

    initialized = .true.
  end subroutine

  !
  ! ---------------------------------------------------------------------
  ! SUBROUTINE DIFFUSE: Driver for calculating sub-grid fluxes (thus it
  ! includes call to surface routines) 
  !
  subroutine diffuse(timein)

    use grid, only : newvar, nstep, a_up, a_ut, a_vp, a_vt, a_wp, a_wt       &
         ,a_rp, a_tp, a_sp, a_st, vapor, a_pexnr, a_theta,a_km               &
         , a_scr1, a_scr2, a_scr3, a_scr4, a_scr5, a_scr6, a_scr7, nscl, nxp, nyp    &
         , nzp, nxyp, nxyzp, zm, dxi, dyi, dzi_t, dzi_m, dt, th00, dn0, pi0, pi1     &
         ,level, uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc,liquid, a_cvrxp, trac_sfc    & 
	 , sgtendt, sgtendr, outtend, rkalpha, rkbeta

    use util, only         : atob, azero, get_avg3
    use mpi_interface, only: cyclics, cyclicc
    use thrm, only         : bruvais, fll_tkrs

    real, intent(in)       :: timein 
    integer :: n
    real :: rk

    ! Hack BvS: slowly increase smago constant...
    !csx = min(timein*0.23/3600.,0.23)    


    if (.not.Initialized) call diffuse_init(nzp, nxp, nyp)

    !
    ! ----------
    ! Calculate Deformation and stability for SGS calculations
    !

    if(level>0) then
      call fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,a_scr1,rs=a_scr2)
      call bruvais(nzp,nxp,nyp,level,a_theta,a_tp,a_scr3,dzi_m,th00,a_rp,a_scr2)
    else
      call fll_tkrs(nzp,nxp,nyp,a_theta,a_pexnr,pi0,pi1,a_scr1)
      call bruvais(nzp,nxp,nyp,level,a_theta,a_tp,a_scr3,dzi_m,th00)
    end if

    !
    !
    call deform(nzp,nxp,nyp,dzi_m,dzi_t,dxi,dyi,a_up,a_vp,a_wp,a_scr5,a_scr6     &
         ,a_scr4,a_scr2)

    ! ----------
    ! Calculate Eddy Viscosity/Diffusivity 
    !
    call smagor(nzp,nxp,nyp,sflg,dxi,dyi,dn0,a_scr3,a_scr2,a_km,a_scr7,zm)
    !
    ! Diffuse momentum
    !
    if (sflg) call acc_tend(nzp,nxp,nyp,a_up,a_vp,a_wp,a_ut,a_vt,a_wt         &
         ,sz4,sz5,sz6,1,'sgs')

    call diff_prep(nzp,nxp,nyp,a_scr5,a_scr6,a_scr4,a_km)

    call azero(nxyp,sxy1,a2=sxy2)

    call diff_vpt(nzp,nxp,nyp,dn0,dzi_m,dzi_t,dyi,dt,vw_sfc,sxy2,a_scr6         &
         ,a_scr5,a_km,a_vp,a_wp,a_vt,sz2)

    call diff_upt(nzp,nxp,nyp,dn0,dzi_m,dzi_t,dxi,dt,uw_sfc,sxy1,a_scr5         &
         ,a_km,a_up,a_wp,a_ut,sz1)

    call diff_wpt(nzp,nxp,nyp,dn0,dzi_m,dzi_t,dyi,dxi,dt,ww_sfc,sxy1,a_scr4     &
         ,a_km,a_wp,a_up,a_wt,sz3)

    call cyclics(nzp,nxp,nyp,a_wt,req)
    call cyclicc(nzp,nxp,nyp,a_wt,req)
    call cyclics(nzp,nxp,nyp,a_vt,req)
    call cyclicc(nzp,nxp,nyp,a_vt,req)
    call cyclics(nzp,nxp,nyp,a_ut,req)
    call cyclicc(nzp,nxp,nyp,a_ut,req)

    if (sflg) then
       call sgs_vel(nzp,nxp,nyp,sz1,sz2,sz3)
       call acc_tend(nzp,nxp,nyp,a_up,a_vp,a_wp,a_ut,a_vt,a_wt,sz4,sz5,sz6    &
            ,2,'sgs')
    end if

    !
    ! Diffuse scalars
    !
    
    if(outtend) then !RV
       !sgtendt = 0.!new at if (sflg)
       !sgtendr = 0.
       if(nstep==1) rk = rkalpha(1)+rkalpha(2)
       if(nstep==2) rk = rkbeta(2)+rkbeta(3)
       if(nstep==3) rk = rkalpha(3)
    end if    !rv

    do n=4,nscl
       call newvar(n,istep=nstep)
       call azero(nxyp,sxy1)
       call azero(nxyp,sxy2)
       if ( associated(a_tp,a_sp) ) call atob(nxyp,wt_sfc,sxy1)
       if ( associated(a_rp,a_sp) ) call atob(nxyp,wq_sfc,sxy1)
       if ( associated(a_cvrxp,a_sp) ) call atob(nxyp,trac_sfc,sxy1)

       if (sflg) call azero(nxyzp,a_scr1)

       if(outtend .and. n<=5) then
          if(n==4) call diffsclr(nzp,nxp,nyp,dt,dxi,dyi,dzi_m,dzi_t,dn0,sxy1,sxy2   & !RV
                     ,a_sp,a_scr2,a_st,a_scr1,sgtendt,rk)
          if(n==5) call diffsclr(nzp,nxp,nyp,dt,dxi,dyi,dzi_m,dzi_t,dn0,sxy1,sxy2   & !RV
                     ,a_sp,a_scr2,a_st,a_scr1,sgtendr,rk)
       else
          call diffsclr(nzp,nxp,nyp,dt,dxi,dyi,dzi_m,dzi_t,dn0,sxy1,sxy2   &
               ,a_sp,a_scr2,a_st,a_scr1)
       endif

  !     if(outtend .and. n==4) print *, 'sgtendt(30:40,30,30) of nstep ',nstep,' is ',sgtendt(30:40,30,30), '.' !RV: Added new print statment for bug fix
    
       if (sflg) then

          call get_avg3(nzp,nxp,nyp,a_scr1,sz1)
          call updtst(nzp,'sgs',n-3,sz1,1)

          if (associated(a_sp,a_tp))                                          &
             call sgsflxs(nzp,nxp,nyp,level,liquid,vapor,a_theta,a_scr1,'tl')
          if (associated(a_sp,a_rp))                                          &
             call sgsflxs(nzp,nxp,nyp,level,liquid,vapor,a_theta,a_scr1,'rt')

          if (outtend .and. n<=5) then !RV
             if(n==4) then
		 call get_avg3(nzp,nxp,nyp,sgtendt,sz1)
		 sgtendt = 0.
!		 print *, 'updtst of sgtendt with average (sz1(30:40)) of ' ,sz1(30:40), '.'  !RV: Added new print
             else if(n==5) then
		 call get_avg3(nzp,nxp,nyp,sgtendr,sz1)
		 sgtendr = 0.
	     end if
             call updtst(nzp,'tnd',n-1,sz1,1)
          endif !rv

       endif
       call cyclics(nzp,nxp,nyp,a_st,req)
       call cyclicc(nzp,nxp,nyp,a_st,req)
    enddo

  end subroutine diffuse
  !
  ! ---------------------------------------------------------------------
  ! subroutine deform: computes the components of the deviatoric strain
  ! tensor, then computes s = du_i/dx_j(du_i/dx_j + du_j/dx_i) at a 
  ! thermo point, the dummy arrays are respectively:
  !
  ! sxz1 -> div; 
  ! szx2 -> s11;
  ! szx3 -> s33, 
  ! szx4 -> s13
  !
  subroutine deform(n1,n2,n3,dzi_m,dzi_t,dx,dy,u,v,w,s12,s22,s23,s)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
    real, intent(in)    :: dzi_m(n1),dzi_t(n1),dx,dy

    real, intent(out)   :: s(n1,n2,n3), s22(n1,n2,n3)
    real, intent(out)   :: s12(n1,n2,n3),s23(n1,n2,n3)

    integer :: ip,im,jp,jm,kp
    real    :: y1a,y2a,y3a,y1b,y2b,y3b
    real    :: s11_wpt,s22_wpt,s33_wpt,s12_wpt,s13_wpt,s23_wpt

    !
    ! calculate components of the stress tensor at their natural locations
    !
    do j=1,n3
       jm=max(j-1,1)
       jp=min(j+1,n3)
       do i=1,n2
          im=max(i-1,1)
          ip=min(i+1,n2)
          do k=1,n1
             szx2(k,i) = 2.*(u(k,i,j)-u(k,im,j))*dx
             s22(k,i,j)= 2.*(v(k,i,j)-v(k,i,jm))*dy
             szx3(k,i) = 2.*(w(k,i,j)-w(max(1,k-1),i,j))*dzi_t(k)
             s12(k,i,j)= (u(k,i,jp)-u(k,i,j))*dy + (v(k,ip,j)-v(k,i,j))*dx
             szx1(k,i) = 0.333333*(szx2(k,i)+s22(k,i,j)+szx3(k,i))
          enddo

          do k=1,n1-1
             szx4(k,i) =(u(k+1,i,j)-u(k,i,j))*dzi_m(k)+(w(k,ip,j)-w(k,i,j))*dx
             s23(k,i,j)=(v(k+1,i,j)-v(k,i,j))*dzi_m(k)+(w(k,i,jp)-w(k,i,j))*dy
          end do
       end do
       !
       ! average to a w-point
       !
       do i=1,n2
          im=max(i-1,1)
          do k=1,n1-1
             kp=k+1
             y1a=(szx2(k,i)-szx1(k,i))
             y2a=(s22(k,i,j)-szx1(k,i))
             y3a=(szx3(k,i)-szx1(k,i))
             y1b=(szx2(kp,i)-szx1(kp,i))
             y2b=(s22(kp,i,j)-szx1(kp,i))
             y3b=(szx3(kp,i)-szx1(kp,i))

             s11_wpt=0.5*(y1a*y1a+y1b*y1b)
             s22_wpt=0.5*(y2a*y2a+y2b*y2b)
             s33_wpt=0.5*(y3a*y3a+y3b*y3b)
             s12_wpt=0.125*(s12(k,i,j)*s12(k,i,j) + s12(kp,i,j)*s12(kp,i,j) &
                  +s12(k,im,j)*s12(k,im,j) + s12(k,i,jm)*s12(k,i,jm)        &
                  +s12(kp,im,j)*s12(kp,im,j)+s12(k,im,jm)*s12(k,im,jm)      &
                  +s12(kp,i,jm)*s12(kp,i,jm)+s12(kp,im,jm)*s12(kp,im,jm))
             s13_wpt=0.5*(szx4(k,i)*szx4(k,i)+szx4(k,im)*szx4(k,im))
             s23_wpt=0.5*(s23(k,i,j)*s23(k,i,j)+s23(k,i,jm)*s23(k,i,jm))
             s(k,i,j)= 0.5*(s11_wpt+s22_wpt+s33_wpt)                        &
                  + s12_wpt + s13_wpt + s23_wpt
          end do
       end do
    end do

  end subroutine deform
  !
  ! ----------------------------------------------------------------------
  ! Subroutine smagor:  computes visc/diff, upon entering the routine
  ! kh is filled with s2 and ri is filled with N^2.  On statistical 
  ! timsteps, SGS energy, dissipation, viscosity, diffusivity and 
  ! lengthscales are stored.
  !
  subroutine smagor(n1,n2,n3,sflg,dxi,dyi,dn0,ri,kh,km,szxy,zm)
      
    use defs, only          : pi, vonk
    use stat, only          : tke_sgs
    use util, only          : get_avg3, get_cor3, calclevel
    use mpi_interface, only : cyclics, cyclicc
    use grid, only          : liquid

    implicit none

    logical, intent(in) :: sflg
    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: dxi,dyi,zm(n1),dn0(n1)
    real, intent(inout) :: ri(n1,n2,n3),kh(n1,n2,n3)
    real, intent(out)   :: km(n1,n2,n3),szxy(n1,n2,n3)
    integer             :: cb, ct
    real    :: delta,pr,labn

    pr    = abs(prndtl)
    
    delta = 1./dxi
    delta = (zm(2)/dxi/dyi)**0.333333333

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-2
             ri(k,i,j) = max( -1., ri(k,i,j)/(kh(k,i,j) + 1.e-12) )
             !
             ! variable km represents what is commonly known as Km, the eddy viscosity
             ! variable kh represents strain rate factor S^2 (dummy variable)
             !
             
             ! Original Bjorn:
             !km(k,i,j) = sqrt(max(0.,kh(k,i,j)*(1.-ri(k,i,j)/pr))) &
             !     *0.5*(dn0(k)+dn0(k+1))/(1./(delta*csx)**2+1./(zm(k)*vonk)**2)

             ! BvS: split out wall damping and stability correction
             labn      = (1./(delta*csx)**2.+1./(zm(k)*vonk)**2.)
             km(k,i,j) =  (dn0(k)+dn0(k+1))/2. * sqrt(max(0.,kh(k,i,j))) * sqrt(max(0.,(1.-ri(k,i,j)/pr))) / labn 
             !
             ! after kh is multiplied with the factor (1-ri/pr), the product of kh 
             ! and km represents the dissipation rate epsilon 
             !
             kh(k,i,j) = kh(k,i,j) *(1.-(ri(k,i,j)/pr))             
          enddo
          kh(1,i,j)    = kh(2,i,j)
          kh(n1,i,j)   = kh(n1-2,i,j)
          km(1,i,j)    = km(2,i,j)
          km(n1,i,j)   = km(n1-2,i,j)
          km(n1-1,i,j) = km(n1-2,i,j)    
       enddo
    enddo

    call cyclics(n1,n2,n3,km,req)
    call cyclicc(n1,n2,n3,km,req)


    if (sflg) then
       call get_cor3(n1,n2,n3,km,km,sz1)
       ! 
       ! The product km and kh represent the local dissipation rate
       ! 
       call get_cor3(n1,n2,n3,km,kh,sz2)
       call updtst(n1,'sgs',-2,sz2,1)      ! dissipation averaged over domain
       do k=1,n1
          !
          ! the factor 1/pi^2 probably represents the ratio of the constants
          ! Cm/Ce that appears in the definition of TKE, the factor csx^2
          ! will cancel out with the csx^2 that appears in the numerator of  
          ! variable sz1 which corresponds to Km^2.
          !
          tke_sgs(k) = sz1(k)/(delta*pi*(csx**2))**2
          sz1(k) = 1./sqrt(1./(delta*csx)**2.+1./(zm(k)*vonk+0.001)**2.)
       end do
       call updtst(n1,'sgs',-1,tke_sgs,1) ! sgs tke
       call updtst(n1,'sgs',-5,sz1,1)      ! mixing length
       call updtst(n1,'sgs',-6,sz1,1)      ! dissipation lengthscale
    end if

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
            !
            ! The product km and kh represent the local dissipation rate
            !
            szxy(k,i,j) = km(k,i,j)*kh(k,i,j)
            !
            ! What is known as the 'physical' eddy diffusivity, Kh, is yet calculated from Km 
            !
            kh(k,i,j) = km(k,i,j)/pr
            if (prndtl < 0.) then
               kh(k,i,j) = kh(k,i,j) * exp(zm(k)/(-100.))
            end if
          enddo
       enddo
    enddo
    if (clouddiff> 0) then ! Additional diffusion outside of the clouds - but in the cloud layer
      call calclevel(liquid, cb, 'base')
      call calclevel(liquid, ct, 'top')    
      do j=3,n3-2
        do i=3,n2-2
            do k=cb,ct
              if (liquid(k,i,j) <1e-10) then
                kh(k,i,j) = clouddiff * kh(k,i,j) 
              end if
            enddo
        enddo
      enddo
        
    
    end if
    call cyclics(n1,n2,n3,kh,req)
    call cyclicc(n1,n2,n3,kh,req)

    if (sflg) then
       call get_avg3(n1,n2,n3,km,sz3)
       call updtst(n1,'sgs',-3,sz3,1)  ! eddy viscosity
       call get_avg3(n1,n2,n3,kh,sz2)
       call updtst(n1,'sgs',-4,sz2,1)  ! eddy diffusivity
    end if

  end subroutine smagor
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_prep: multiplies the strain components computed in 
  ! "deform" by the appropriately averaged value of km, in preperation
  ! for use by the diffusion routines
  !
  subroutine  diff_prep(n1,n2,n3,s12,s22,s23,km)

    integer, intent(in) :: n1,n2,n3
    real, intent(inout) :: s22(n1,n2,n3)
    real, intent(inout) :: s12(n1,n2,n3),s23(n1,n2,n3)
    real, intent(inout) :: km(n1,n2,n3)

    integer :: ip,jp

    do j=2,n3-1
       jp=min(j+1,n3)
       do i=2,n2-1
          ip=min(i+1,n2)
          do k=2,n1-1
             s22(k,i,j)=-s22(k,i,j)*0.5*(km(k,i,j)+km(k-1,i,j))
             s12(k,i,j)= -s12(k,i,j)*0.125*(km(k,i,j)+km(k,ip,j)+km(k,i,jp) &
                  +km(k,ip,jp)+km(k-1,i,j)+km(k-1,ip,j)+km(k-1,i,jp)        &
                  +km(k-1,ip,jp))
          enddo

          do k=1,n1-1
             s23(k,i,j)=-s23(k,i,j)*0.5*(km(k,i,j)+km(k,i,jp))
          enddo
       enddo
    enddo

  end subroutine diff_prep
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_upt: computes the diffusivity of velocities using
  ! a tri-diagnonal solver in the vertical at a u point, the deformation
  ! tensor component d31 is passed in via the tendency array
  !
  subroutine  diff_upt(n1,n2,n3,dn0,dzi_m,dzi_t,dxi,dt,sflx,tflx,sij,km,   &
       u,w,tnd,flx)

    logical, parameter  :: noslip = .false.
    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: sij(n1,n2,n3)
    real, intent(in)    :: km(n1,n2,n3),u(n1,n2,n3),w(n1,n2,n3)
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3)
    real, intent(in)    :: dn0(n1),dzi_m(n1),dzi_t(n1),dxi,dt

    real, intent(inout) :: flx(n1), tnd(n1,n2,n3)

    dti   = 1.0/dt
    do k=1,n1
       sz7(k) = 0.
       sz8(k) = 0.
       flx(k) = 0.
       szx5(k,n2) = 0.
    end do

    do j=3,n3-2
       indh=0
       do i=3,n2-2
          indh=indh+1
          sz8(1)=dzi_m(1)*(km(1,i,j)+km(1,i+1,j))
          sz7(n1-1)  =.5*tflx(i,j)*(dn0(n1)+dn0(n1-1))
          sz7(1)     =.5*sflx(i,j)*(dn0(1)+dn0(2))
          do k=2,n1-1
             sz8(k)=dzi_m(k)*(km(k,i,j)+km(k,i+1,j))
             sz7(k)= (-(w(k,i+1,j)-w(k,i,j))*dxi)*0.5*(km(k,i,j)+km(k,i+1,j))
             sxz4(indh,k)=u(k,i,j)*dn0(k) + dt*dzi_t(k)*(sz7(k-1)-sz7(k))
             sxz3(indh,k)=-0.5*dt*dzi_t(k)*sz8(k)
             sxz2(indh,k)=-0.5*dt*dzi_t(k)*sz8(k-1)
             sxz1(indh,k)=dn0(k)-sxz2(indh,k)-sxz3(indh,k)
          end do
          !
          ! Boundary conditions
          !
          if (noslip) then
             sxz1(indh,2)    = dn0(2)-2.*sxz2(indh,2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)-2.*sxz3(indh,n1-1)
          else
             sxz1(indh,2)    = dn0(2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)
          end if

          sxz3(indh,n1-1) = 0.
          sxz2(indh,2)    = 0.
          flx(1)=flx(1)+sz7(1)
          flx(n1-1)=flx(n1-1)+sz7(n1-1)

          do k=2,n1-1
             szx5(k,i) = (-2.*(u(k,i,j)-u(k,i-1,j))*dxi)*0.5*                 &
                  (km(k,i,j)+km(k-1,i,j))
          end do
       enddo

       do k=2,n1
          szx5(k,n2-1) = (-2.*(u(k,n2-1,j)-u(k,n2-2,j))*dxi)*0.5*             &
                  (km(k,n2-1,j)+km(k-1,n2-1,j))
       end do

       call tridiff(n2,n1-1,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and 
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !   
       indh=0
       do i=3,n2-2
          indh=indh+1
          do k=2,n1-1
             dfact = 1./dn0(k)
             tnd(k,i,j)=tnd(k,i,j) + dti*(sxz1(indh,k)-u(k,i,j)) - dfact *    &
                  ((szx5(k,i+1)-szx5(k,i))*dxi + (sij(k,i,j)-sij(k,i,j-1))*dxi)

             if (k < n1-1) flx(k)= flx(k)-dzi_m(k)*(km(k,i,j)+km(k,i+1,j))      &
                  *(sxz1(indh,k+1)-sxz1(indh,k))*.5
          enddo
          tnd(n1,i,j) = 0.
       enddo
    enddo
   
  end subroutine diff_upt
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_vpt: computes the diffusivity of velocities using
  ! a tri-diagnonal solver in the vertical at u or v pts depending on
  ! the values of ip and jp and the input arguments
  !
  subroutine  diff_vpt(n1,n2,n3,dn0,dzi_m,dzi_t,dyi,dt,sflx,tflx,sii,sij,  &
       km,v,w,tnd,flx)

    logical, parameter  :: noslip = .false.
    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: dn0(n1),dzi_m(n1),dzi_t(n1),dyi,dt
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3)
    real, intent(in)    :: sii(n1,n2,n3),sij(n1,n2,n3)
    real, intent(in)    :: km(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

    real, intent(inout) :: flx(n1)
    real, intent(inout) :: tnd(n1,n2,n3)

    dti = 1.0/dt
    do k=1,n1
       sz7(k)  = 0.
       sz8(k)  = 0.
       flx(k)  = 0.
    end do

    do j=3,n3-2
       indh  = 0
       do i=3,n2-2
          indh=indh+1
          sz8(1)=dzi_m(1)*(km(1,i,j)+km(1,i,j+1))
          sz7(n1-1)  =.5*(tflx(i,j))*(dn0(n1)+dn0(n1-1))
          sz7(1)     =.5*(sflx(i,j))*(dn0(1)+dn0(2))
          do k=2,n1-1
             sz8(k)=dzi_m(k)*(km(k,i,j)+km(k,i,j+1))
             sz7(k)= (-(w(k,i,j+1)-w(k,i,j))*dyi)*0.5*(km(k,i,j)+km(k,i,j+1))
             sxz4(indh,k)=v(k,i,j)*dn0(k) + dt*dzi_t(k)*(sz7(k-1)-sz7(k))
             sxz3(indh,k)=-0.5*dt*dzi_t(k)*sz8(k)
             sxz2(indh,k)=-0.5*dt*dzi_t(k)*sz8(k-1)
             sxz1(indh,k)=dn0(k)-sxz2(indh,k)-sxz3(indh,k)
          end do
          !
          ! Boundary conditions
          !
          if (noslip) then
             sxz1(indh,2)    = dn0(2)-2.*sxz2(indh,2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)-2.*sxz3(indh,n1-1)
          else
             sxz1(indh,2)    = dn0(2)-sxz3(indh,2)
             sxz1(indh,n1-1) = dn0(n1-1)-sxz2(indh,n1-1)
          end if

          sxz3(indh,n1-1) = 0.
          sxz2(indh,2)    = 0.
          flx(1)=flx(1)+sz7(1)
          flx(n1-1)=flx(n1-1)+sz7(n1-1)
       enddo

       call tridiff(n2,n1-1,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
       !
       ! Back out diffusive tendencies from vertical implicit solution and 
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh=0
       do i=3,n2-2
          indh=indh+1
          do k=2,n1-1
             dfact = 1./dn0(k) 
             tnd(k,i,j)=tnd(k,i,j) + dti*(sxz1(indh,k)-v(k,i,j)) - dfact *    &
                 ((sii(k,i,j+1)-sii(k,i,j))*dyi+(sij(k,i,j)-sij(k,i-1,j))*dyi)
             if (k < n1-1) flx(k)= flx(k)-dzi_m(k)*(km(k,i,j)+km(k,i,j+1))      &
                  *(sxz1(indh,k+1)-sxz1(indh,k))*.5
          enddo
       enddo
    enddo

  end subroutine diff_vpt
  !
  ! ----------------------------------------------------------------------
  ! subroutine diff_wpt: computes the diffusivity of velocities at a
  ! wpt
  !
  subroutine  diff_wpt(n1,n2,n3,dn0,dzi_m,dzi_t,dxi,dyi,dt,sflx,tflx,s23,km,w,u   &
       ,tnd,flx)

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: s23(n1,n2,n3)
    real, intent(in)    :: km(n1,n2,n3),w(n1,n2,n3),u(n1,n2,n3)
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3)
    real, intent(in)    :: dn0(n1),dzi_m(n1),dzi_t(n1),dxi,dyi,dt

    real, intent(inout) :: flx(n1), tnd(n1,n2,n3)

    integer :: kp1,im1,jm1

    dti = 1.0/dt
    do k=1,n1
       sz8(k)  = 0.
       flx(k)  = 0.
    end do

    do k=1,n1
       do i=1,n2
          sxz1(i,k) = 0.0
       end do
    end do

    do j=3,n3-2
       indh=0
       do i=3,n2-2
          indh=indh+1

          sz8(1)=dzi_t(2)*.5*(km(1,i,j)+km(2,i,j))
          do k=2,n1-2
             kp1 = k+1
             sz8(k)=dzi_t(kp1)*.5*(km(k,i,j)+km(kp1,i,j))
             sxz4(indh,k)=w(k,i,j)*(dn0(k)+dn0(kp1))*.5
             sxz3(indh,k)=-dt*dzi_m(k)*sz8(k)
             sxz2(indh,k)=-dt*dzi_m(k)*sz8(k-1)
             sxz1(indh,k)=(dn0(k)+dn0(kp1))*0.5 - sxz2(indh,k) - sxz3(indh,k)
          end do
          sxz2(indh,2)    = 0.
          sxz3(indh,n1-2) = 0.

          flx(1)=flx(1)+sflx(i,j)*dn0(2)
          flx(n1-1)=flx(n1-1)+tflx(i,j)*dn0(n1-1)

          do k=2,n1-1
             szx5(k,i) = ((u(k+1,i,j)-u(k,i,j))*dzi_m(k)+(w(k,i+1,j)-w(k,i,j))  &
                  *dxi)*(-0.5)*(km(k,i,j)+km(k,i+1,j))
          end do
          sxz4(indh,2)   =sxz4(indh,2) + dt*dzi_m(2)*sflx(i,j)*dn0(2)
          sxz4(indh,n1-2)=sxz4(indh,n1-2)-dt*dzi_m(n1-2)*tflx(i,j)*dn0(n1-1)
       end do

       do k=2,n1-1
          szx5(k,2) =  ((u(k+1,2,j)-u(k,2,j))*dzi_m(k) + (w(k,3,j)-w(k,2,j))    &
               *dxi)*(-0.5)*(km(k,2,j)+km(k,3,j))
       end do

       call tridiff(n2,n1-2,indh,sxz2,sxz1,sxz3,sxz4,sxz1,sxz3)
            

       !
       ! Back out diffusive tendencies from vertical implicit solution and 
       ! solve for horizontal tendencies.  Accumulate vertical flux array
       !
       indh= 0
       do i=3,n2-2
          indh=indh+1
          im1 = max(i-1,2)
          jm1 = max(j-1,2)
          do k=2,n1-2
             dfact = 1./((dn0(k)+dn0(k+1))*.5)
             tnd(k,i,j)=tnd(k,i,j) + dti*(sxz1(indh,k)-w(k,i,j))- dfact *     &
                  ((szx5(k,i)-szx5(k,im1))*dxi + (s23(k,i,j)-s23(k,i,jm1))*dyi)
             flx(k) = flx(k)-dzi_t(k)*(km(k,i,j)+km(k+1,i,j))*0.5               &
                  *(sxz1(indh,k)-sxz1(indh,k-1))
          enddo
       enddo
    enddo

  end subroutine diff_wpt
  !
  ! -----------------------------------------------------------------------
  ! subroutine diffsclr: computes the diffusivity of a scalar using
  ! a tri-diagnonal solver in the vertical
  !
  subroutine diffsclr(n1,n2,n3,dt,dxi,dyi,dzi_m,dzi_t,dn0,sflx,tflx,scp,xkh,sct &
       ,flx,flxtend,rk)

    use grid, only : outtend

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: xkh(n1,n2,n3),scp(n1,n2,n3)
    real, intent(in)    :: sflx(n2,n3),tflx(n2,n3),dn0(n1)
    real, intent(in)    :: dxi,dyi,dzi_m(n1),dzi_t(n1),dt

    real, optional, intent(inout)   :: flxtend(n1,n2,n3) !RV
    real, optional, intent(in)   :: rk !rv

    real, intent(inout) :: flx(n1,n2,n3),sct(n1,n2,n3)
    !
    ! compute vertical diffusion matrix coefficients for scalars,
    ! Coefficients need only be calculated once and can be used repeatedly
    ! for other scalars
    !

    dti       = 1.0/dt
    do k=1,n1
       sz7(k)   = 0.
    end do

    do j=3,n3-2
       do i=2,n2-2
          do k=2,n1-1
             szx1(k,i)=-(scp(k,i+1,j)-scp(k,i,j))*dxi*.25*(xkh(k,i,j)  +     &
                  xkh(k,i+1,j)+xkh(k-1,i,j)+xkh(k-1,i+1,j)) 
          enddo
       enddo
       !
       ! Set up Tri-diagonal Matrix
       !
       indh=0
       do i=3,n2-2
          indh=indh+1
          do k=2,n1-1
             if (k < n1-1) sz7(k)=dt*dzi_m(k)*xkh(k,i,j)
             sxz1(indh,k)=-dzi_t(k)*sz7(k-1)
             sxz2(indh,k)=-dzi_t(k)*sz7(k)
             sxz3(indh,k)=dn0(k)-sxz1(indh,k)-sxz2(indh,k)
             sxz4(indh,k)=scp(k,i,j)*dn0(k)
          enddo
          sxz4(indh,2)=scp(2,i,j)*dn0(2)                                     &
               + sflx(i,j)*(dn0(1)+dn0(2))     *.5 *dt*dzi_t(2)
          sxz4(indh,n1-1)=scp(n1-1,i,j)*dn0(n1-1)                            &
               - tflx(i,j)*(dn0(n1-1)+dn0(n1)) *.5 *dt*dzi_t(n1-1)
       enddo

       call tridiff(n2,n1-1,indh,sxz1,sxz3,sxz2,sxz4,sxz5,sxz6)
       !
       ! compute scalar tendency in addition to vertical flux
       !
       indh=0
       do i=3,n2-2
          flx(1,i,j)   =sflx(i,j)*(dn0(1)+dn0(2))*.5
          flx(n1-1,i,j)=tflx(i,j)*(dn0(n1)+dn0(n1-1))*.5
          flx(n1,i,j)  =0.
          indh=indh+1
          do k=2,n1-1
             sct(k,i,j)= sct(k,i,j) + dti*(sxz5(indh,k)-scp(k,i,j))           &
                  -((szx1(k,i)-szx1(k,i-1))                                   &
                  *dxi + (-(scp(k,i,j+1)-scp(k,i,j))*dyi*0.25*(xkh(k,i,j)     &
                  +xkh(k,i,j+1)+xkh(k-1,i,j)+xkh(k-1,i,j+1))+(scp(k,i,j)      &
                  -scp(k,i,j-1))*dyi*0.25*(xkh(k,i,j-1)+xkh(k,i,j)            &
                  +xkh(k-1,i,j-1)+xkh(k-1,i,j)))*dyi) /dn0(k)
             if (k<n1-1) flx(k,i,j)=-xkh(k,i,j)*(sxz5(indh,k+1)-sxz5(indh,k)) &
                  *dzi_m(k)
             if(outtend) then   !RV: t&r tendencies stored
                flxtend(k,i,j)= flxtend(k,i,j) + dti*rk*(sxz5(indh,k)-scp(k,i,j))   &
                     -rk*((szx1(k,i)-szx1(k,i-1))*dxi                               &
                     + (-(scp(k,i,j+1)-scp(k,i,j))*dyi*0.25*(xkh(k,i,j)          &
                     +xkh(k,i,j+1)+xkh(k-1,i,j)+xkh(k-1,i,j+1))+(scp(k,i,j)      &
                     -scp(k,i,j-1))*dyi*0.25*(xkh(k,i,j-1)+xkh(k,i,j)            &
                     +xkh(k-1,i,j-1)+xkh(k-1,i,j)))*dyi) /dn0(k)
             end if !rv
          end do
       enddo
    enddo

  end subroutine diffsclr

end module sgsm



