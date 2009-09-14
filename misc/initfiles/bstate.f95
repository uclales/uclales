!
! This program generates a basic state for a case that depends on the value of
! icase.  The default is the DYCOMS-RF01 case, using a grid  which gives a 
! constant change in theta above the inversion.  icase = 0 gives the initial 
! state for the cumulus topped boundary layer with exponentially decreasing
! moisture profiles, icase=1 gives the RF02 initial conditions
!
Program bstate

  implicit none
  real, parameter :: cp = 1005.
  real, parameter :: Rd = 287.15
  real, parameter :: g  = 9.81
  real, parameter :: p0 = 100000.

  integer :: icase

  print *, "Enter case: CBL (1), RF02 (2), RF01 (default, or 0)"
  read  *, icase
  select case(icase)
  case(1)
     call cbl (290.0, 6.e-3, 0.6)
  case(2)
     call gcss8
  case default
     call gcss7
  end select

  stop
contains

  subroutine cbl(thv0,gamma,rh)

  implicit none
  real, intent (in)  :: gamma, thv0, rh

  real :: z, dz, prs, rv, th, t, thv, rv0
  integer :: k, iterate

  open (10,file='sound_in',status='unknown')
  dz = 5.0
  do k=1,150
     thv = thv0 + gamma * z
     prs = p0 * (1.-(g/(cp*gamma)) * log ( 1. + gamma * z/thv0))**(cp/Rd)
     rv  = 8.0e-3
     rv0 = 0.0e-3
     iterate = 0
     do while (iterate < 100 .and. abs(1.-rv0/rv) > 1.e-5)
        th  = (thv/(1.+0.608*rv))
        t   = th*(prs/p0)**(Rd/cp)
        rv0 = rv
        rv  = get_rv(prs,t,rh)
        iterate = iterate + 1
     end do
     if (k == 1) then
        write (10,'(5F10.3)') 1000., th, rv*1000., 0.0, 0.0
     else
        write (10,'(5F10.3)') z, th, rv*1000., 0.0, 0.0
     end if
     z = z + dz
     dz  = dz*1.025
  end do
  end subroutine cbl

  subroutine gcss7
    implicit none

    real, parameter :: dtheta = 0.15
    real    :: zz

    open (10,file='sound_in',status='unknown')
    write (10,'(5F10.3)') 1017.8, 289.00, 9.05, 7.0, -5.5
    write (10,'(5F10.3)') 839.95, 289.00, 9.05, 7.0, -5.5
    zz = 840.05
    do while (zz < 2500.)
       write (10,'(5F10.3)') zz, 297. + (zz - 840.)**(1./3.), 1.50, 7.0, -5.5
       zz = ((zz-840.)**(1./3.) + dtheta)**3. + 840.
    end do

    stop
  end subroutine gcss7

  subroutine gcss8
    implicit none

    real, parameter :: dtheta = 0.15
    real    :: zz, qq, uu, vv

    open (10,file='sound_in',status='unknown')
    write (10,'(5F10.3)') 1017.8, 288.70, 9.45, 6.50, -8.500
    write (10,'(5F10.3)') 794.95, 288.70, 9.45, 8.09, -4.525
    zz = 795.05
    do while (zz < 2500.)
       qq = 5.0-3.0*(1.-exp(-(zz-795.)/500.))
       uu = 6.50 + 2*zz/1000.
       vv =-8.50 + 5*zz/1000.
       write (10,'(5F10.3)') zz, 295. + (zz - 795.)**(1./3.), qq, uu, vv
       zz = ((zz-795.)**(1./3.) + dtheta)**3. + 840.
    end do

    stop
  end subroutine gcss8
  !
  ! ---------------------------------------------------------------------
  ! This function calculates the liquid saturation vapor mixing ratio as
  ! a function of temperature and pressure
  !
  real function get_rv(p,t,rh)

    implicit none
    real, intent (in) :: p, t, rh
    real, parameter :: c0=0.6105851e+03, c1=0.4440316e+02, c2=0.1430341e+01   &
         , c3= .2641412e-01, c4= .2995057e-03, c5= .2031998e-05               &
         , c6= .6936113e-08, c7= .2564861e-11, c8=-.3704404e-13 

    real ::  esl, x

    x=max(-80.,t-273.16)

    !     esl=612.2*exp(17.67*x/(t-29.65)) 
    esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    get_rv=.622*esl*rh/(p-esl*rh)

    return
  end function get_rv
  !

end Program bstate
