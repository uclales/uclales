!
! This program generates the time varying forcings for different cases. The defa! ult is the ASTEX case
!
Program lscale

  implicit none

  integer :: icase

  print *, "Enter case: ASTEX (0)"
  read  *, icase
  select case(icase)
  case(0)
     call astex 
  end select

  stop
contains

  subroutine astex
    implicit none

     real,dimension(4)        :: time_obs=(/0.,0.31,0.78,1.71/)
     real,dimension(3)        :: tgeo=(/10.,19.,36./)
     real,dimension(3)        :: ugeo_obs=(/-2.,-2.,-0.75/)
     real,dimension(3)        :: vgeo_obs=(/-10.,-10.,-4./)
     real,dimension(41)       :: div_var,sst_var
     real        :: hour_frac, del_sst,time_in,ugeo_grad,vgeo_grad
     real        :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
     integer     :: tcnt,n,tgeo_cnt,nln
     real        :: sst, div, u0,v0


div_var(36:41)=-1.
sst_var(36:41)=(/294.480,294.510,294.540,294.570,294.600,294.630/)

  open (1,file='lagr1_scalars_Chris',status='old',form='formatted')
       do nln=1,35
          read (1,*) tmp1,tmp2,tmp3,&
                             tmp4,tmp5,sst_var(nln),tmp6,div_var(nln)
          print *, sst_var(nln),div_var(nln)                 
       end do
       close (1)




    open (10,file='lscale_in',status='unknown')


    tcnt=2

    do n=0,40

    time_in=n*1.
    hour_frac=time_in/24.

   !sst
    sst=sst_var(n+1)

   !div 

    div=div_var(n+1)

    !geostrophic winds

    if (time_in .le. tgeo(2)) then 
     u0 = ugeo_obs(1)
     v0 = vgeo_obs(1)
    endif

    if (time_in .gt. tgeo(2)) then
     vgeo_grad = (vgeo_obs(3)-vgeo_obs(2))/(tgeo(3)-tgeo(2))
     v0 = vgeo_obs(2) + (time_in-tgeo(2))*vgeo_grad
     u0 = ugeo_obs(3) 
    endif

    if (time_in .ge. tgeo(3)) then 
     u0 = ugeo_obs(2)   
     v0 = vgeo_obs(3)
    endif

     
    write (10,'(5F10.3)') time_in, div, sst,u0,v0


    end do

    stop
  end subroutine astex


end Program lscale
