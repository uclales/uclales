!----------------------------------------------------------------------------
! centered differences scheme for scalar advection. Original code provided by 
! Georgis Matheou (Matheou et al., 2011, MWR, 139, 2918-2939), 
!modified by Linda Schlemmer, March 2012
!----------------------------------------------------------------------------

module centered

  implicit none

  save 

  logical :: advRainMonFlag = .true.

contains

  subroutine advection_scalars(scheme)
    
    character(len=8), intent(in) :: scheme

    
    select case(trim(scheme))
    case('second')
        call scalarSecond
    case('third')
        call scalarThird
    case('fourth')
        call scalarFourth
    end select
        
  end subroutine advection_scalars
  
!-----------------------------------------------------------------------------
  subroutine scalarSecond
    use advf, only : mamaos, mamaos_x, mamaos_y, advtnd
    use stat, only : sflg, updtst
    use util, only : atob, get_avg3
    use grid, only : nzp, nxp, nyp, nxyzp, dxi, dyi, dzi_t, dzi_m, nscl, &
                     a_up, a_vp, a_wp, a_sp, a_st, newvar, dn0, &
                     a_scr1, a_scr2, nstep, dt
     
    real, dimension(:,:,:), allocatable    :: rhow
    real, dimension(:,:,:,:), allocatable :: flux
    integer :: ret
    integer :: n
    integer :: nghost(3)
    real    :: v1da(nzp)

    allocate(rhow(nzp, nxp, nyp), stat=ret)
    if(ret.ne.0) then
        print *, 'advect scalar second: auxiliary-vectors allocation failed'
        return
    endif
    nghost(:)=2
    call advection_make_rhow(nzp, nxp, nyp, a_wp, dn0, rhow)
     
    ! flux(direction, k, i, j)
    allocate(flux(3, nzp, nxp, nyp), stat=ret)
    if(ret.ne.0) then
        print *, 'advect scalar second: auxiliary-vectors allocation failed'
        return
    endif
     
     
    do n=4, nscl
        call newvar(n, istep=nstep)

        if (n.ge.6)  then
           call atob(nxyzp,a_sp,a_scr1)
           call mamaos_y(nzp,nxp,nyp,a_vp,a_sp,a_scr1,dyi,dt)
           call mamaos_x(nzp,nxp,nyp,a_up,a_sp,a_scr1,dxi,dt)
           call atob(nxyzp,a_wp,a_scr2)
           call mamaos(nzp,nxp,nyp,a_scr2,a_sp,a_scr1,dzi_t,dzi_m,dn0,dt,.false.)

           if (sflg) then
              call get_avg3(nzp,nxp,nyp,a_scr2,v1da)
              call updtst(nzp,'adv',n-3,v1da,1)
           end if

           call advtnd(nzp,nxp,nyp,a_sp,a_scr1,a_st,dt)
           cycle
        end if

        call scalarFluxSecond(nzp, nxp, nyp, a_up, a_vp, rhow, a_sp, flux)
                                   
        call scalarFluxDifferencesSecond(nzp, nxp, nyp, nghost, &
                   dxi, dyi, dzi_t, a_st, dn0, flux)
     
    end do
     
    deallocate(flux)
  
  end subroutine scalarSecond
!------------------------------------------------------------------------------
  subroutine scalarThird

! third order advection scheme. Fluxes are computed with 3rd order approximation, 
! whereas the differencing is done second order. See Wicker and Skamarock (2002) 
! Linda Schlemmer, 2012

    use advf, only : mamaos, mamaos_x, mamaos_y, advtnd
    use grid, only : nzp, nxp, nyp, nxyzp, dxi, dyi, dzi_t, dzi_m, nscl, &
                     a_up, a_vp, a_wp, a_sp, a_st, newvar, dn0, &
                     a_scr1, a_scr2, nstep, dt
    use stat, only : sflg, updtst
    use util, only : atob, get_avg3
     
    real, dimension(:,:,:), allocatable    :: rhow
    real, dimension(:,:,:), allocatable :: flux2z
    real, dimension(:,:,:), allocatable :: flux3x, flux3y
    integer :: ret
    integer :: n
    integer :: nghost(3),req(16)
    real    :: v1da(nzp)
     
    allocate(flux2z(nzp, nxp, nyp), stat=ret)
    flux2z(:,:,:)=0.0
    if(ret.ne.0) then
        print *, 'advect scalar third: auxiliary-vectors allocation failed'
        return
    endif
     
    allocate(flux3x(nzp, nxp, nyp), stat=ret)
    flux3x(:,:,:) = 0.0
    allocate(flux3y(nzp, nxp, nyp), stat=ret)
    flux3y(:,:,:) = 0.0
     
    nghost(:)=2

    allocate(rhow(nzp, nxp, nyp), stat=ret)
    rhow(:,:,:) = 0.0
   if(ret.ne.0) then
        print *, 'advect scalar third: auxiliary-vectors allocation failed'
        return
    endif
    nghost(:)=2
    call advection_make_rhow(nzp, nxp, nyp, a_wp, dn0, rhow)

    do n=4, nscl
         call newvar(n, istep=nstep)
         if (n.ge.8)  then
            call atob(nxyzp,a_sp,a_scr1)
            call mamaos_y(nzp,nxp,nyp,a_vp,a_sp,a_scr1,dyi,dt)
            call mamaos_x(nzp,nxp,nyp,a_up,a_sp,a_scr1,dxi,dt)
            call atob(nxyzp,a_wp,a_scr2)
            call mamaos(nzp,nxp,nyp,a_scr2,a_sp,a_scr1,dzi_t,dzi_m,dn0,dt,.false.)

            if (sflg) then
              call get_avg3(nzp,nxp,nyp,a_scr2,v1da)
              call updtst(nzp,'adv',n-3,v1da,1)
            end if

            call advtnd(nzp,nxp,nyp,a_sp,a_scr1,a_st,dt)
            cycle
         end if

         call scalarFluxThird(nzp, nxp, nyp, a_up, a_vp, rhow, a_sp, flux2z, flux3x, flux3y)

         call scalarFluxDifferencesSecond_2(nzp, nxp, nyp, nghost, &
                    dxi, dyi, dzi_t, a_st, dn0, flux2z, flux3x, flux3y)

    end do
     
    deallocate(flux2z)
    deallocate(flux3x)
    deallocate(flux3y)
  
  end subroutine scalarThird

  subroutine scalarFourth

    use mpi_interface, only : cyclics, cyclicc
    use advf, only : mamaos, mamaos_x, mamaos_y, advtnd
    use grid, only : nzp, nxp, nyp, nxyzp, dxi, dyi, dzi_t, dzi_m, nscl, &
                     a_up, a_vp, a_wp, a_sp, a_st, newvar, dn0, &
                     a_scr1, a_scr2, nstep, dt
    use stat, only : sflg, updtst
    use util, only : atob, get_avg3
     
    real, dimension(:,:,:), allocatable    :: rhow
    real, dimension(:,:,:,:), allocatable :: flux1!, flux3
    real, dimension(:,:,:), allocatable :: flux3a, flux3b, flux3c
    integer :: ret
    integer :: n
    integer :: nghost(3),req(16)
    real    :: v1da(nzp)
     
    ! flux(direction, k, i, j)
    allocate(flux1(3, nzp, nxp, nyp), stat=ret)
    flux1(:,:,:,:)=0.0
    if(ret.ne.0) then
        print *, 'advect scalar second: auxiliary-vectors allocation failed'
        return
    endif
     
!    allocate(flux3(3, nzp, nxp, nyp), stat=ret)
!    flux3(:,:,:,:) = 0.0
!    if(ret.ne.0) then
!        print *, 'advect scalar second: auxiliary-vectors allocation failed'
!        return
!    endif
    allocate(flux3a(nzp, nxp, nyp), stat=ret)
    flux3a(:,:,:) = 0.0
    allocate(flux3b(nzp, nxp, nyp), stat=ret)
    flux3b(:,:,:) = 0.0
    allocate(flux3c(nzp, nxp, nyp), stat=ret)
    flux3c(:,:,:) = 0.0
     
    nghost(:)=2

    allocate(rhow(nzp, nxp, nyp), stat=ret)
    rhow(:,:,:) = 0.0
   if(ret.ne.0) then
        print *, 'advect scalar second: auxiliary-vectors allocation failed'
        return
    endif
    nghost(:)=2
    call advection_make_rhow(nzp, nxp, nyp, a_wp, dn0, rhow)

    do n=4, nscl
         call newvar(n, istep=nstep)
         if (n.ge.6)  then
            call atob(nxyzp,a_sp,a_scr1)
            call mamaos_y(nzp,nxp,nyp,a_vp,a_sp,a_scr1,dyi,dt)
            call mamaos_x(nzp,nxp,nyp,a_up,a_sp,a_scr1,dxi,dt)
            call atob(nxyzp,a_wp,a_scr2)
            call mamaos(nzp,nxp,nyp,a_scr2,a_sp,a_scr1,dzi_t,dzi_m,dn0,dt,.false.)

            if (sflg) then
              call get_avg3(nzp,nxp,nyp,a_scr2,v1da)
              call updtst(nzp,'adv',n-3,v1da,1)
            end if

            call advtnd(nzp,nxp,nyp,a_sp,a_scr1,a_st,dt)
            cycle
         end if

         call scalarFluxFourth(nzp, nxp, nyp, a_up, a_vp, rhow, a_sp, flux1, flux3a, flux3b, flux3c)

!         flux3a(:,:,:) = flux3(1,:,:,:)
!         flux3b(:,:,:) = flux3(2,:,:,:)
!         flux3c(:,:,:) = flux3(3,:,:,:)

         call cyclics(nzp,nxp,nyp,flux3a,req)
         call cyclicc(nzp,nxp,nyp,flux3a,req)
         call cyclics(nzp,nxp,nyp,flux3b,req)
         call cyclicc(nzp,nxp,nyp,flux3b,req)
         call cyclics(nzp,nxp,nyp,flux3c,req)
         call cyclicc(nzp,nxp,nyp,flux3c,req)

         call scalarFluxDifferencesFourth(nzp, nxp, nyp, nghost, &
                    dxi, dyi, dzi_t, a_st, dn0, flux1, flux3a, flux3b, flux3c)

         call cyclics(nzp,nxp,nyp,a_st,req)
         call cyclicc(nzp,nxp,nyp,a_st,req)
     
    end do
     
    deallocate(flux1)
    deallocate(flux3a)
    deallocate(flux3b)
    deallocate(flux3c)
!    deallocate(flux3)
  
  end subroutine scalarFourth

  subroutine scalarFluxDifferencesSecond(nz, nx, ny, nghost, wx1, wy1, wz1, Zrhs, rho, flux)
     integer, intent(in) :: nz, nx, ny, nghost(3)
     real, intent(in) :: wx1, wy1
     real, dimension(nz), intent(in) ::rho, wz1
     real, dimension(3,nz,nx,ny), intent(in) :: flux
     real, dimension(nz,nx,ny), intent(inout) :: Zrhs
     
     integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
     
     ilo=nghost(1)+1
     jlo=nghost(2)+1
     klo=nghost(3)+1
               
     ihi=nx-nghost(1)
     jhi=ny-nghost(2)
     khi=nz-nghost(3)

     do j=jlo, jhi
         do i=ilo, ihi
             do k=klo, khi
                 Zrhs(k,i,j)=Zrhs(k,i,j) &
                             -(flux(1,k,i+1,j)-flux(1,k,i,j))*wx1      &
                             -(flux(2,k,i,j+1)-flux(2,k,i,j))*wy1      &
                             -(flux(3,k+1,i,j)-flux(3,k,i,j))*wz1(k)/rho(k) 
             end do
         end do
     end do

  end subroutine scalarFluxDifferencesSecond


  subroutine scalarFluxDifferencesSecond_2(nz, nx, ny, nghost, wx1, wy1, wz1, Zrhs, rho, flux2z, flux3x, flux3y)
     integer, intent(in) :: nz, nx, ny, nghost(3)
     real, intent(in) :: wx1, wy1
     real, dimension(nz), intent(in) ::rho, wz1
     real, dimension(nz,nx,ny), intent(in) :: flux2z, flux3x, flux3y
     real, dimension(nz,nx,ny), intent(inout) :: Zrhs
     
     integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
     
     ilo=nghost(1)+1
     jlo=nghost(2)+1
     klo=2
               
     ihi=nx-nghost(1)
     jhi=ny-nghost(2)
     khi=nz-1

     do j=jlo, jhi
         do i=ilo, ihi
             do k=klo, khi
                 Zrhs(k,i,j)=Zrhs(k,i,j) &
                             -(flux3x(k,i+1,j)-flux3x(k,i,j))*wx1      &
                             -(flux3y(k,i,j+1)-flux3y(k,i,j))*wy1      &
                             -(flux2z(k+1,i,j)-flux2z(k,i,j))*wz1(k)/rho(k)
             end do
         end do
     end do

  end subroutine scalarFluxDifferencesSecond_2


  subroutine scalarFluxDifferencesFourth(nz, nx, ny, nghost, dxi, dyi, dzi, Zrhs, rho, flux1, flux3a, flux3b, flux3c)
    integer, intent(in) :: nz, nx, ny, nghost(3)
    real, intent(in) :: dxi, dyi
    real, dimension(nz), intent(in) :: rho(nz), dzi(nz)
    real, dimension(3,nz,nx,ny), intent(in) :: flux1
    real, dimension(nz,nx,ny), intent(in) :: flux3a, flux3b, flux3c
    real, dimension(nz,nx,ny), intent(inout) :: Zrhs

    integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
    real :: wx1, wy1, wz1, wx3, wy3, wz3

    ilo=nghost(1)
    jlo=nghost(2)
    klo=nghost(3)

    ihi=nx-nghost(1)
    jhi=ny-nghost(2)
    khi=nz-nghost(3)

    wx1=+9.*dxi/8.
    wx3=-1.*dxi/(8.*3.)
    wy1=+9.*dyi/8.
    wy3=-1.*dyi/(8.*3.)

    do j=jlo, jhi
        do i=ilo, ihi
            do k=klo, khi
                wz1=+9.*dzi(k)/8.
                wz3=-1.*dzi(k)/(8.*3.)
                Zrhs(k,i,j)=Zrhs(k,i,j) &
                             -(flux1(1,k,i+1,j)-flux1(1,k,i,j))*wx1      &
                             -(flux1(2,k,i,j+1)-flux1(2,k,i,j))*wy1      &
                             -(flux1(3,k+1,i,j)-flux1(3,k,i,j))*wz1/rho(k)      &
                             -(flux3a(k,i+2,j)-flux3a(k,i-1,j))*wx3    &
                             -(flux3b(k,i,j+2)-flux3b(k,i,j-1))*wy3    &
                             -(flux3c(k+2,i,j)-flux3c(k-1,i,j))*wz3/rho(k)
             end do
         end do
     end do

  end subroutine scalarFluxDifferencesFourth


  subroutine scalarFluxSecond(nz, nx, ny, u, v, w, Z, flux)
     integer, intent(in) :: nz, nx, ny
     real, dimension(nz,nx,ny), intent(in) :: u, v, w, Z
     real, dimension(3,nz,nx,ny), intent(out) :: flux
     
     integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi

     ilo=2
     jlo=2
     klo=2
               
     ihi=nx-1
     jhi=ny-1
     khi=nz-1

     do j=jlo, jhi
         do i=ilo, ihi
             do k=klo, khi
                 flux(1,k,i,j)=u(k,i-1,j)*(Z(k,i-1,j)+Z(k,i,j))*0.5d0
                 flux(2,k,i,j)=v(k,i,j-1)*(Z(k,i,j-1)+Z(k,i,j))*0.5d0
                 flux(3,k,i,j)=w(k-1,i,j)*(Z(k-1,i,j)+Z(k,i,j))*0.5d0
             end do
         end do
     end do
     
  end subroutine scalarFluxSecond

  subroutine scalarFluxThird(nz, nx, ny, u, v, w, Z, flux2z, flux3x, flux3y)
    integer, intent(in) :: nz, nx, ny
    real, dimension(nz,nx,ny), intent(in) :: u, v, w, Z
    real, dimension(nz,nx,ny), intent(out) :: flux2z
    real, dimension(nz,nx,ny), intent(out) :: flux3x, flux3y
    real :: flux4x, flux4y

    integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi

    ilo=3
    jlo=3
    klo=1

    ihi=nx-1
    jhi=ny-1
    khi=nz

    do j=jlo,jhi
        do i=ilo,ihi
            do k=klo,khi

               flux4x = (7.0d0*(Z(k,i-1,j)+Z(k,i,j))-(Z(k,i-2,j)+Z(k,i+1,j)))*u(k,i-1,j)*0.0833333333d0
               flux4y = (7.0d0*(Z(k,i,j-1)+Z(k,i,j))-(Z(k,i,j-2)+Z(k,i,j+1)))*v(k,i,j-1)*0.0833333333d0

               flux3x(k,i,j) = flux4x - dabs(u(k,i-1,j))/12.0d0*(3.0d0*(Z(k,i,j)-Z(k,i-1,j))-(Z(k,i+1,j)-Z(k,i-2,j)))
               flux3y(k,i,j) = flux4y - dabs(v(k,i,j-1))/12.0d0*(3.0d0*(Z(k,i,j)-Z(k,i,j-1))-(Z(k,i,j+1)-Z(k,i,j-2)))

            end do

            flux2z(1,i,j)=0.0d0

            do k=klo+1,khi

               flux2z(k,i,j)=(Z(k-1,i,j)+Z(k,i,j))*w(k-1,i,j)*0.5d0

            end do
        end do
    end do

  end subroutine scalarFluxThird

  subroutine scalarFluxFourth(nz, nx, ny, u, v, w, Z, flux1, flux3a, flux3b, flux3c)
    integer, intent(in) :: nz, nx, ny
    real, dimension(nz,nx,ny), intent(in) :: u, v, w, Z
    real, dimension(3,nz,nx,ny), intent(out) :: flux1
    real, dimension(nz,nx,ny), intent(out) :: flux3a, flux3b, flux3c

    integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi

    ilo=3
    jlo=3
    klo=2

    ihi=nx-1
    jhi=ny-1
    khi=nz-1

    do j=jlo,jhi
        do i=ilo,ihi
            do k=klo,khi
               flux1(1,k,i,j)=(Z(k,i-1,j)+Z(k,i,j))*u(k,i-1,j)*0.5d0
               flux1(2,k,i,j)=(Z(k,i,j-1)+Z(k,i,j))*v(k,i,j-1)*0.5d0
               flux1(3,k,i,j)=(Z(k-1,i,j)+Z(k,i,j))*w(k-1,i,j)*0.5d0

               flux3a(k,i,j)=(Z(k,i-2,j)+Z(k,i+1,j))*u(k,i-1,j)*0.5d0
               flux3b(k,i,j)=(Z(k,i,j-2)+Z(k,i,j+1))*v(k,i,j-1)*0.5d0
            end do
            flux3c(1,i,j)=0.
            flux3c(2,i,j)=Z(3,i,j)*w(1,i,j)
            do k=klo+1,khi
               flux3c(k,i,j)=(Z(k-2,i,j)+Z(k+1,i,j))*w(k-1,i,j)*0.5d0
            end do
        end do
    end do

  end subroutine scalarFluxFourth


  subroutine advection_make_rhow(nz, nx, ny, w, rho, rhow)
    integer, intent(in) :: nz, nx, ny
    real, dimension(nz), intent(in) :: rho
    real, dimension(nz,nx,ny), intent(in) :: w
    real, dimension(nz,nx,ny), intent(out) :: rhow

    integer :: i, j, k

    do j=1,ny
        do i=1,nx
            do k=1,nz-1
                rhow(k,i,j)=0.5*(rho(k)+rho(k+1))*w(k,i,j)
            end do
        end do
    end do

  end subroutine advection_make_rhow

end module centered
