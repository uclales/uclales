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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
! PRSS: Pressure Solver:  Solves the anelastic or bousinessq system 
! for the pressure using a fractional step method, which is implemented
! using fft's and a tri-diagonal solver in the vertical
!
module prss

  use defs, only: pi
  implicit none

contains
!
!----------------------------------------------------------------------
! subroutine poisson: called by timesteping driver to invert the 
! poisson equation for pressure and apply the velocity tendencies.
!
  subroutine poisson

    use grid, only : nxp, nyp, nzp, dxi, dyi, dzi_m, dzi_t, dt, a_up, a_ut,      &
         a_vp, a_vt, a_wp, a_wt, press, a_pexnr, th00, dn0, wsavex, wsavey
    use stat, only : fill_scalar, sflg
    use util, only : ae1mm

    complex, allocatable     :: s1(:,:,:)
    real    :: mxdiv, awpbar(nzp)
    integer :: ix,iy

    ix=max(1,nxp-4)
    iy=max(1,nyp-4)

    allocate (s1(ix,iy,nzp))
    s1=0.0
    !
    ! ------
    ! Pressure Solve
    !

    call poiss(nzp,nxp,nyp,ix,iy,a_ut,a_vt,a_wt,a_up,a_vp,a_wp,a_pexnr,     &
         press,dn0,th00,dzi_t,dzi_m,dxi,dyi,dt,s1,wsavex,wsavey)
    call ae1mm(nzp,nxp,nyp,a_wp,awpbar)

    !
    ! -------

    if (sflg) then
  
       call get_diverg(nzp,nxp,nyp,ix,iy,s1,a_up,a_vp,a_wp, &
            dn0,dzi_t,dxi,dyi,dt,mxdiv)
       call fill_scalar(2,mxdiv)
       call prs_cor(nzp,nxp,nyp,a_pexnr,a_up,a_vp,a_wp,dzi_m,dxi,dyi,th00)
    end if
   
    deallocate (s1)

  end subroutine poisson
  !
  ! --------------------------------------------------------------------
  ! subroutine poiss: called each timestep to evaluate the pressure
  ! in accordance with the anelastic continuity equation, and then apply
  ! the pressure to the velocity terms for three dimensional flow, 
  ! cyclic in x and y.  pp and pc are used as scratch arrays in the
  ! call to trdprs.  pp is filled with its diagnostic value in fll_prs
  !
  subroutine poiss(n1,n2,n3,ix,iy,ut,vt,wt,u,v,w,pp,pc,dn0,th00,dzi_t,dzi_m, &
       dx,dy,dt,s1,wsvx,wsvy)

    use util, only  : get_fft_twodim, velset
    use grid, only  : rkalpha, rkbeta, nstep

    integer :: n1,n2,n3,ix,iy
    real    :: pp(n1,n2,n3),pc(n1,n2,n3),dmy
    real, dimension(n1,n2,n3)  :: ut, vt, wt, u, v, w
    real    :: wsvx(1:),wsvy(1:),dn0(n1),dzi_t(n1),dzi_m(n1),dx,dy,dt,th00
    complex :: s1(ix,iy,n1)

    call get_diverg(n1,n2,n3,ix,iy,s1,u,v,w,dn0,dzi_t,dx,dy,dt,dmy)
    call get_fft_twodim(ix,iy,n1,s1,wsvx,wsvy,-1)
    call trdprs(n1,ix,iy,s1,dn0,dzi_t,dzi_m,dx,dy)
    call get_fft_twodim(ix,iy,n1,s1,wsvx,wsvy,+1)

    ! The following routines calculate the pressure gradient contribution to the
    ! velocity tendencies, note that just the velocities and the divergence are updated, 
    ! and not the tendencies themselves (these are not used anymore after the call to the
    ! poisson solver within the main stepper in step.f90, in fact, they are put to zero at 
    ! the start of the new RK timestep in the call to 'tendencies'

    call prs_grad(n1,n2,n3,ix,iy,s1,pp,u,v,w,dzi_m,dx,dy,dt)   
    call get_diverg(n1,n2,n3,ix,iy,s1,u,v,w,dn0,dzi_t,dx,dy,dt,dmy)

    pp(:,:,:) = pp(:,:,:)/th00/(rkalpha(nstep)+rkbeta(nstep))
    pc(:,:,:) = pp(:,:,:)
    
  end subroutine poiss
  !
  ! --------------------------------------------------------------------
  ! subroutine get_diverg: gets velocity tendency divergence and puts it 
  ! into a complex value array for use in pressure calculation
  !
  subroutine get_diverg(n1,n2,n3,ix,iy,s1,u,v,w,dn0,dz,dx,dy,dt,mxdiv)

    integer, intent (in)  :: n1,n2,n3,ix,iy
    real, intent (in)     :: dz(n1),dn0(n1),dx,dy,dt
    real, dimension (n1,n2,n3), intent (in) :: u, v, w
    real, intent (out)    :: mxdiv
    complex, intent (out) :: s1(ix,iy,n1)

    integer :: k,i,j,l,m
    real    :: xf,yf,zf,wf1,wf2,dti

    s1(:,:,:) = (0.0,0.0)
    dti = 1./dt
    m=0
    do j=3,n3-2
       m=m+1
       l=0
       do i=3,n2-2
          l=l+1
          do k=2,n1-1
             wf1=0.5*(dn0(k+1)+dn0(k)) 
             wf2=0.5*(dn0(k)+dn0(k-1)) 
             if (k == 2 )   wf2=0.
             if (k == n1-1) wf1=0.
             xf=dn0(k)*dx 
             yf=dn0(k)*dy 
             zf=dz(k)
             s1(l,m,k)= &!(wf1*wt(k,i,j)-wf2*wt(k-1,i,j))*zf                  &
                  !+(vt(k,i,j)-vt(k,i,j-1))*yf +(ut(k,i,j)-ut(k,i-1,j))*xf  &
                  ((wf1*w(k,i,j)-wf2*w(k-1,i,j))*zf                     &
                  +(v(k,i,j)-v(k,i,j-1))*yf + (u(k,i,j)-u(k,i-1,j))*xf)
          enddo
       enddo
    enddo
           

!
! save mxdiv to a statistics array, no reduction necessary as this is done
! in post processing
!
    mxdiv = maxval(real(s1))
    !print *, mxdiv

  end subroutine get_diverg
  !
  !----------------------------------------------------------------------
  ! subroutine prs_grad: writes the pressure to the appropriate array and
  ! calculates the pressure gradient contribution to the tendency
  ! Note that not the tendencies, but the actual velocities are updated and output
  ! within this routine
  !
  subroutine prs_grad(n1,n2,n3,ix,iy,s1,p,u,v,w,dz,dx,dy,dt)

    use mpi_interface, only : cyclics,cyclicc

    integer, intent (in) :: n1,n2,n3,ix,iy
    real, intent (in)    :: dz(n1),dx,dy,dt
    real, intent (inout) :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
    real, intent (out)   :: p(n1,n2,n3)
    complex, intent (in) :: s1(ix,iy,n1)

    integer :: i,j,k,l,m,req(16)

    p(:,:,:)=0.0
    do k=2,n1-1
       l=0
       do i=3,n2-2
          l=l+1
          m=0
          do j=3,n3-2
             m=m+1
             p(k,i,j) =real(s1(l,m,k))
          enddo
       enddo
    enddo
    call cyclics(n1,n2,n3,p,req)
    call cyclicc(n1,n2,n3,p,req)

! Pressure in following calculation contains a timestep term

    do j=1,n3-1
       do i=1,n2-1
          do k=2,n1-1
             if (k /= n1-1) w(k,i,j) = w(k,i,j)-dz(k)*(p(k+1,i,j)-p(k,i,j))
             u(k,i,j) = u(k,i,j) - dx*(p(k,i+1,j)-p(k,i,j))
             v(k,i,j) = v(k,i,j) - dy*(p(k,i,j+1)-p(k,i,j))
          end do
       end do
    end do

    p = p/dt
  end subroutine prs_grad
  !
  !---------------------------------------------------------------------
  ! TRDPRS: solves for the wave number (l,m) component of 
  ! pressure in a vertical column using a tri-diagonal solver.
  !
  subroutine trdprs(n1,ix,iy,s1,dn0,dzi_t,dzi_m,dx,dy)  

    use mpi_interface, only : yoffset, nypg, xoffset, wrxid, wryid, nxpg
    use util, only          : tridiff

    integer, intent (in)    :: n1,ix,iy
    real, intent (in)       :: dn0(n1),dzi_t(n1),dzi_m(n1),dx,dy
    complex, intent (inout) :: s1(ix,iy,n1)

    real    :: ak(ix,n1),dk(ix,n1),bk(ix,n1),ck(ix,n1)
    real    :: xk(ix,n1),yk(ix,n1),wv(ix,iy)

    integer :: k,l,m
    real    :: fctl,fctm,xl,xm,af,cf
    integer :: xof, yof

    fctl=2.*pi/float(nxpg-4)
    fctm=2.*pi/float(nypg-4)

    xof=xoffset(wrxid)
    yof=yoffset(wryid)

    do l=1,ix
          if(l+xof .le. int((nxpg-4)/2)+1) then
            xl=float(l-1+xof)
          else
             xl=float(l-(nxpg-4)-1+xof)
          endif
      
       do m=1,iy
          if(m+yof .le. int((nypg-4)/2)+1) then
             xm=float(m-1+yof)
          else
             xm=float(m-(nypg-4)-1+yof)
          endif
          wv(l,m)=2.*((cos(fctl*xl)-1.)*dx*dx + (cos(fctm*xm)-1.)*dy*dy)
       enddo
    enddo

    if(wrxid.eq.0 .and. wryid .eq.0 ) then
       wv(1,1)=0.
    endif
    !
    ! configure vectors for tri-diagonal solver
    !
    do m=1,iy
       do k=2,n1-1 
          af=(dn0(k)+dn0(k-1))*.5
          cf=(dn0(k+1)+dn0(k))*.5
          if (k == 2   )af=0.
          if (k == n1-1)cf=0.
          do l=1,ix
             ak(l,k)=dzi_t(k)*dzi_m(k-1)*af
             bk(l,k)=s1(l,m,k)
             ck(l,k)=dzi_t(k)*dzi_m(k)*cf
             dk(l,k)=dn0(k)*wv(l,m)-(ak(l,k)+ck(l,k))
          enddo
       enddo
       !
       ! solve for fourier components, x_k, given a tri-diagonal matrix of the
       ! form a_k x_k-1 + d_k x_k + c_k x_k+1 = b_k.  y_k is a scratch array.
       !
       call tridiff(ix,n1-1,ix,ak,dk,ck,bk,xk,yk)

       do k=2,n1-1
          do l=1,ix
             if (m+yof+l+xof>2) bk(l,k)=aimag(s1(l,m,k))
             if (m+yof+l+xof>2) s1(l,m,k)=xk(l,k)
          enddo
       enddo
      
       call tridiff(ix,n1-1,ix,ak,dk,ck,bk,xk,yk)

       do k=2,n1-1
          do l=1,ix
             if (m+yof+l+xof > 2)                &
                  s1(l,m,k)=cmplx(real(s1(l,m,k)),xk(l,k))
          enddo
       enddo

    enddo

  end subroutine trdprs
  !
  !---------------------------------------------------------------------
  ! Subroutine Prs_cor: correlate the pressure tendency with velocity
  ! field for TKE budget
  !
  subroutine prs_cor(n1,n2,n3,p,u,v,w,dz,dx,dy,th00)

    use stat, only : updtst
    use util, only : get_cor

    integer, intent (in) :: n1,n2,n3
    real, intent (in)    :: p(n1,n2,n3),dz(n1),dx,dy,th00
    real, intent (in)    :: u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)

    real, dimension (n2,n3) :: pgx, pgy, pgz, ufld, vfld, wfld
    real    :: fx, fy, fz, v1da(n1), v1db(n1), v1dc(n1)
    integer :: i,j,k,ip1,jp1

    v1da = 0.0
    v1db = 0.0
    v1dc = 0.0

    do k=2,n1-1
       fx=dx*th00
       fy=dy*th00
       fz=dz(k)*th00
       do j=1,n3
          do i=1,n2
             ip1 = min(n2,i+1)
             jp1 = min(n3,j+1)
             pgx(i,j) = -fx*(p(k,ip1,j)-p(k,i,j))
             pgy(i,j) = -fy*(p(k,i,jp1)-p(k,i,j))
             pgz(i,j) = -fz*(p(k+1,i,j)-p(k,i,j))
             ufld(i,j) = u(k,i,j)
             vfld(i,j) = v(k,i,j)
             wfld(i,j) = w(k,i,j)
          end do
       end do
       v1da(k) = get_cor(1,n2,n3,1,ufld,pgx)
       v1db(k) = get_cor(1,n2,n3,1,vfld,pgy)
       v1dc(k) = get_cor(1,n2,n3,1,wfld,pgz)
    enddo
    call updtst(n1,'prs',1,v1da,1)
    call updtst(n1,'prs',2,v1db,1)
    call updtst(n1,'prs',3,v1dc,1)

  end subroutine prs_cor

end module prss
