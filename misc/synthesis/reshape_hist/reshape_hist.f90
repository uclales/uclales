Program reshape_hist

!-----------------------------------------------------------------------------!
!                                                                             !
! the program reads in restart files from UCLA-LES from nx1 x ny1 processors, !
! gathers the fields and redistributes them onto nx2 x ny2 fields and writes  !
! new restart files.                                                          !
!                                                                             !
! In order to save memory the programm reads only in a part of the restart    !
! files and puts that to the new files immediately.                           !
!                                                                             !
! Compilation: 1. compile the modules separately:                             !
!                 ifort -c -r8 read_hist_1.f90                                !
!              2. compile the main program and link the modules:              !
!                 ifort -r8 -o reshape_hist read_hist_1.o read_hist_srfc.o \  !
!                 read_hist_srfc_rad.o read_hist_2.o read_hist_3.o \          !
!                 write_hist_1.o write_hist_2.o write_hist_3.o \              !
!                 write_hist_srfc.o write_hist_srfc_rad.o reshape_hist.f90    !
!              Note: it is important to use -r8 since the restart files are   !
!                    double precision.                                        !
!                                                                             !
! Linda Schlemmer, July 2013                                                  !
!-----------------------------------------------------------------------------!

  use readhist_1
  use readhist_srfc
  use readhist_srfc_rad
  use readhist_2
  use readhist_3
  use writehist_1
  use writehist_srfc
  use writehist_srfc_rad
  use writehist_2
  use writehist_3

  implicit none

  integer           :: isfctyp = 0         ! surface flux parameterization type
  integer           :: iradtyp = 0         ! radiation model type
  integer           :: level   = 0         ! thermodynamic level
  integer           :: naddsc  = 0         ! number of additional scalars
  integer           :: ihtype
  integer           :: nz, nxt, nyt        ! grid size of gathered field
  integer           :: nx1, ny1, nxp1, nyp1! grid size of original field
  integer           :: nx2, ny2, nxp2, nyp2! grid size of new field
  real              :: th00, umean, vmean
  integer           :: nseed, n
 
  real              :: dt, psrf
  real              :: time
  real              :: timestamp
  character(len=40) :: hname
  integer           :: nscl = 4
  integer           :: iblank, nv2, nsmp
  logical           :: lwaterbudget, lcouvreux

  real, dimension(:), allocatable ::  &
       xt,  &
       xm,  &
       yt,  &
       ym,  &
       zt,  &
       zm,  &
       u0,  &
       v0,  &
       pi0, &
       pi1, &
       th0, &
       dn0, &
       rt0 

  real, dimension (:,:), allocatable :: &
       a_ustar, &
       a_tstar, &
       a_rstar, &
       cnd_acc, &  ! accumulated condensation  [kg/m2] (diagnostic for 2D output)
       cev_acc, &  ! accumulated evaporation of cloud water [kg/m2] (diagnostic for 2D output)
       rev_acc     ! accumulated evaporation of rainwater   [kg/m2] (diagnostic for 2D output)

  real, dimension (:,:,:), allocatable :: &
       prc_acc

  real, dimension (:,:), allocatable :: &
       svctr

  real, dimension (:,:,:),  allocatable :: & 
       a_pexnr, &
       a_tsoil, &
       a_phiw,  &
       a_flx

  real, dimension (:,:),    allocatable :: &
       a_tskin, &
       a_qskin, &
       a_Wl

  real, dimension (:,:,:), allocatable:: &
       a_xp
  integer, dimension(:), allocatable :: seed

  !--------------------------------------------------------------------------
  ! Read the parameters of the restart file
  !--------------------------------------------------------------------------


  print*,'Did you make subdirectory .out/ ?'
  read*,hname

  print*,'stem of the filename (e.g. rico)'
  read*,hname

  print*,'select type of restart file: 1 - .R'
  print*,'                             2 - .rst'
  print*,'                             3 - .time s'
  read*,ihtype

  select case(ihtype)
  case(1)
     hname = trim(hname)//'.R'
  case(2)
     hname = trim(hname)//'.rst'
  case(3)
     print*,'give time stamp'
     read*,timestamp
     iblank=index(hname,' ')
     write (hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(timestamp), 's'
  end select

  print*,'Number of processors nxp, nyp used for original run'
  read*,nxp1, nyp1

  print*,'nx, ny total points of field (without ghost cells)'
  read*,nxt, nyt
  nxt=nxt+4
  nyt=nyt+4
  print*,'Number of processors nxp, nyp for new run'
  read*,nxp2, nyp2

  if (mod(nxt-4,nxp2).ne.0) then
     print"(/'*** Wrong number of processors specified, ***'/'*** impossible to split the field, ***'/'*** chosse different nxp! ***'/)"
     stop
  end if
  if (mod(nyt-4,nyp2).ne.0) then
     print"(/'*** Wrong number of processors specified, ***'/'*** impossible to split the field, ***'/'*** chosse different nyp! ***'/)"
     stop
  end if

  print*,'nz'
  read*, nz
  print*,'level (0,1,2,3,4,5)'
  read*, level

  print*,'lwaterbudget (.true./.false.  default is .false.)'
  read*, lwaterbudget

  if(.not.lwaterbudget) lwaterbudget=.false.

  print*,'isfctyp (only important, if isfctyp=5)'
  read*,isfctyp

  print*,'lcouvreux (.true./.false.  default is .false.)'
  read*, lcouvreux

  if(.not.lcouvreux) lcouvreux=.false.

  if (level   > 0) nscl = nscl+1  ! qt only
  if (level   > 2) nscl = nscl+2  ! nr,qr
  if (level   > 3) nscl = nscl+4  ! ni,qi,qs,qg
  if (level   > 4) nscl = nscl+4  ! ns,ng,qh,nh (for Axel's two-moment scheme)

  if (lwaterbudget) then
     nscl = nscl+1 ! additional cloud water a_cld in the tracer array
  end if
  if (lcouvreux) then
     nscl = nscl+1 ! Additional radioactive scalar
  end if

  ! Determine size of subdomain of each processor

  nx1 = (nxt-4) / nxp1 +4
  ny1 = (nyt-4) / nxp1 +4
  nx2 = (nxt-4) / nxp2 +4
  ny2 = (nyt-4) / nxp2 +4

  ! Allocate the fields of the large gathered domain

  allocate (xt(nxt))
  allocate (yt(nyt))
  allocate (zt(nz))

  allocate (xm(nxt))
  allocate (ym(nyt))
  allocate (zm(nz))

  allocate (u0(nz),v0(nz),pi0(nz),pi1(nz),th0(nz),dn0(nz),rt0(nz))
  u0(:)  = 0.
  v0(:)  = 0.
  pi0(:) = 0.
  pi1(:) = 0.
  th0(:) = 0.
  dn0(:) = 0.
  rt0(:) = 0.
  
  allocate (a_pexnr(nz,nxt,nyt))
  a_pexnr(:,:,:) = 0.

  allocate (a_ustar(nxt,nyt),a_tstar(nxt,nyt),a_rstar(nxt,nyt))
  a_ustar(:,:) = 0.
  a_tstar(:,:) = 0.
  a_rstar(:,:) = 0.

  ! Read in the old restart files

  print*,'**************************************************'
  print*,'Reading in the first part of the old restart files'
  print*,'**************************************************'

  call read_hist_1(time, hname, nxp1, nyp1, nx1, ny1, nxt, nyt, nz, nscl,&
       umean, vmean, th00, level, isfctyp, lwaterbudget, iradtyp, xt, xm, yt,   &
       ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf, seed, nseed, dt, a_ustar,     &
       a_tstar, a_rstar, a_pexnr)

  ! Write out the new restart files

  print*,'***************************************************'
  print*,'Writing out the first part of the new restart files'
  print*,'***************************************************'

  call write_hist_1(time, hname, nxp2, nyp2, nx2, ny2, nxt, nyt, nz, nscl,&
       umean, vmean, th00, level, isfctyp, lwaterbudget, iradtyp, xt, xm, yt,   &
       ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf, seed, nseed, dt,      &
       a_ustar, a_tstar, a_rstar, a_pexnr)

  ! Deallocate the first part of the fields of the large gathered domain

  deallocate (xt)
  deallocate (yt)
  deallocate (zt)

  deallocate (xm)
  deallocate (ym)
  deallocate (zm)

  deallocate (u0,v0,pi0,pi1,th0,dn0,rt0)
  
  deallocate(seed)
  deallocate (a_ustar,a_tstar,a_rstar)
  deallocate (a_pexnr)

  ! check if the land surface model is used. If so, process the entries in the
  ! restart files

  if (isfctyp == 5) then

     allocate (a_tsoil(4,nxt,nyt))
     allocate (a_phiw (4,nxt,nyt))
     allocate (a_tskin (nxt,nyt))
     allocate (a_qskin (nxt,nyt))
     allocate (a_Wl    (nxt,nyt))

     a_tsoil(:,:,:)     = 0.
     a_phiw(:,:,:)      = 0.
     a_tskin(:,:)       = 0.
     a_qskin(:,:)       = 0.
     a_wl(:,:)          = 0.

     print*,'********************************************************'
     print*,'Reading in the surface fields from the old restart files'
     print*,'********************************************************'

     call read_hist_srfc(nxp1, nyp1, nx1, ny1, nxt, nyt,   &
          a_tsoil, a_phiw, a_tskin, a_qskin, a_wl)

     call write_hist_srfc(nxp2, nyp2, nx2, ny2, nxt, nyt,   &
          a_tsoil, a_phiw, a_tskin, a_qskin, a_wl)

     deallocate (a_tsoil)
     deallocate (a_phiw)
     deallocate (a_tskin)
     deallocate (a_qskin)
     deallocate (a_Wl)

     allocate (a_flx(100,nxt,nyt))

     a_flx(:,:,:) = 0.
     ! loop through the radiation fields separately, memory issues!
     do n=1,8
        call read_hist_srfc_rad(nxp1, nyp1, nx1, ny1, nxt, nyt, a_flx)
        call write_hist_srfc_rad(nxp2, nyp2, nx2, ny2, nxt, nyt, a_flx)
     end do
     deallocate (a_flx)

  end if

  ! in order to save memory the 3D fields are read and written
  ! one by one

  allocate (a_xp(nz,nxt,nyt))
  a_xp(:,:,:) = 0.

  print*,'***************************************************'
  print*,'Writing and reading of a_xp'
  print*,'***************************************************'

  do n=1,nscl
     print*,'n = ',n 
     call read_hist_2(nxp1, nyp1, nx1, ny1, nxt, nyt, nz, a_xp)
     call write_hist_2(nxp2, nyp2, nx2, ny2, nxt, nyt, nz, a_xp)
  end do

  deallocate (a_xp)

  ! process the remaining fields

  allocate (prc_acc(nxt,nyt,0))
  prc_acc(:,:,:) = 0.   ! accumulated precipitation for 2D output  [kg/m2]
  allocate (rev_acc(nxt,nyt))
  rev_acc(:,:) = 0.   ! accumulated evaporation of rain water    [kg/m2]
  if (lwaterbudget) then
     allocate (cnd_acc(nxt,nyt),cev_acc(nxt,nyt))
     cnd_acc(:,:) = 0.   ! accumulated condensation                 [kg/m2]
     cev_acc(:,:) = 0.   ! accumulated evaporation of cloud water   [kg/m2]
  end if

  print*,'**************************************************'
  print*,'Reading in the last part of the old restart files'
  print*,'**************************************************'

  if (.not.lwaterbudget) then
     call read_hist_3(hname, nxp1, nyp1, nx1, ny1, nxt, nyt, nz, level, &
          nv2, nsmp, svctr, prc_acc, rev_acc)
  else
     call read_hist_3(hname, nxp1, nyp1, nx1, ny1, nxt, nyt, nz, level, &
          nv2, nsmp, svctr, prc_acc, rev_acc, cnd_acc, cev_acc)
  end if

  ! Write out the new restart files

  print*,'***************************************************'
  print*,'Writing out the last part of the new restart files'
  print*,'***************************************************'

  if (.not.lwaterbudget) then
     call write_hist_3(hname, nxp2, nyp2, nx2, ny2, nxt, nyt, nz, level, &
          nv2, nsmp, svctr, prc_acc, rev_acc)
  else
     call write_hist_3(hname, nxp2, nyp2, nx2, ny2, nxt, nyt, nz, level, &
          nv2, nsmp, svctr, prc_acc, rev_acc, cnd_acc, cev_acc)
  end if

  deallocate (prc_acc)
  deallocate (rev_acc)
  if (lwaterbudget) then
     deallocate (cnd_acc,cev_acc)
  end if
  deallocate (svctr)


End Program reshape_hist
