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
    double precision, dimension(:,:,:,:), allocatable :: flux
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
  subroutine scalarFourth
    use advf, only : mamaos, mamaos_x, mamaos_y, advtnd
    use grid, only : nzp, nxp, nyp, nxyzp, dxi, dyi, dzi_t, dzi_m, nscl, &
                     a_up, a_vp, a_wp, a_sp, a_st, newvar, dn0, &
                     a_scr1, a_scr2, nstep, dt
    use stat, only : sflg, updtst
    use util, only : atob, get_avg3
     
    real, dimension(:,:,:), allocatable    :: rhow
    double precision, dimension(:,:,:,:), allocatable :: flux1, flux3
    integer :: ret
    integer :: n
    integer :: nghost(3)
    real    :: v1da(nzp)
     
    ! flux(direction, k, i, j)
    allocate(flux1(3, nzp, nxp, nyp), stat=ret)
    if(ret.ne.0) then
        print *, 'advect scalar second: auxiliary-vectors allocation failed'
        return
    endif
     
    allocate(flux3(3, nzp, nxp, nyp), stat=ret)
    if(ret.ne.0) then
        print *, 'advect scalar second: auxiliary-vectors allocation failed'
        return
    endif
     
    nghost(:)=2

    allocate(rhow(nzp, nxp, nyp), stat=ret)
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

         call scalarFluxFourth(nzp, nxp, nyp, a_up, a_vp, rhow, a_sp, flux1, flux3)
                                   
         call scalarFluxDifferencesFourth(nzp, nxp, nyp, nghost, &
                    dxi, dyi, dzi_t, a_st, dn0, flux1, flux3)
     
    end do
     
    deallocate(flux1)
    deallocate(flux3)
  
  end subroutine scalarFourth

  subroutine scalarFluxDifferencesSecond(nz, nx, ny, nghost, wx1, wy1, wz1, Zrhs, rho, flux)
     integer, intent(in) :: nz, nx, ny, nghost(3)
     double precision, intent(in) :: wx1, wy1
     double precision, dimension(nz), intent(in) ::rho, wz1
     double precision, dimension(3,nz,nx,ny), intent(in) :: flux
     double precision, dimension(nz,nx,ny), intent(inout) :: Zrhs
     
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


  subroutine scalarFluxDifferencesFourth(nz, nx, ny, nghost, dxi, dyi, dzi, Zrhs, rho, flux1, flux3)
    integer, intent(in) :: nz, nx, ny, nghost(3)
    double precision, intent(in) :: dxi, dyi
    doubleprecision, dimension(nz), intent(in) :: rho(nz), dzi(nz)
    double precision, dimension(3,nz,nx,ny), intent(in) :: flux1, flux3
    double precision, dimension(nz,nx,ny), intent(inout) :: Zrhs

    integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi
    double precision :: wx1, wy1, wz1, wx3, wy3, wz3

    ilo=nghost(1)+1
    jlo=nghost(2)+1
    klo=nghost(3)+1

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
                             -(flux3(1,k,i+2,j)-flux3(1,k,i-1,j))*wx3    &
                             -(flux3(2,k,i,j+2)-flux3(2,k,i,j-1))*wy3    &
                             -(flux3(3,k+2,i,j)-flux3(3,k-1,i,j))*wz3/rho(k)
             end do
         end do
     end do

  end subroutine scalarFluxDifferencesFourth


  subroutine scalarFluxSecond(nz, nx, ny, u, v, w, Z, flux)
     integer, intent(in) :: nz, nx, ny
     double precision, dimension(nz,nx,ny), intent(in) :: u, v, w, Z
     double precision, dimension(3,nz,nx,ny), intent(out) :: flux
     
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

  subroutine scalarFluxFourth(nz, nx, ny, u, v, w, Z, flux1, flux3)
    integer, intent(in) :: nz, nx, ny
    double precision, dimension(nz,nx,ny), intent(in) :: u, v, w, Z
    double precision, dimension(3,nz,nx,ny), intent(out) :: flux1, flux3

    integer :: i, j, k, ilo, ihi, jlo, jhi, klo, khi

    ilo=3
    jlo=3
    klo=3

    ihi=nx-1
    jhi=ny-1
    khi=nz-1

    do j=jlo,jhi
        do i=ilo,ihi
            do k=klo,khi
               flux1(1,k,i,j)=(Z(k,i-1,j)+Z(k,i,j))*u(k,i-1,j)*0.5d0
               flux1(2,k,i,j)=(Z(k,i,j-1)+Z(k,i,j))*v(k,i,j-1)*0.5d0
               flux1(3,k,i,j)=(Z(k-1,i,j)+Z(k,i,j))*w(k-1,i,j)*0.5d0

               flux3(1,k,i,j)=(Z(k,i-2,j)+Z(k,i+1,j))*u(k,i-1,j)*0.5d0
               flux3(2,k,i,j)=(Z(k,i,j-2)+Z(k,i,j+1))*v(k,i,j-1)*0.5d0
               flux3(3,k,i,j)=(Z(k-2,i,j)+Z(k+1,i,j))*w(k-1,i,j)*0.5d0
            end do
        end do
    end do

  end subroutine scalarFluxFourth


  subroutine advection_make_rhow(nz, nx, ny, w, rho, rhow)
    integer, intent(in) :: nz, nx, ny
    double precision, dimension(nz), intent(in) :: rho
    double precision, dimension(nz,nx,ny), intent(in) :: w
    double precision, dimension(nz,nx,ny), intent(out) :: rhow

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
