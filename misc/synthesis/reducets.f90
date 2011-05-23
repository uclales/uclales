        Program reduce

        use netcdf
        implicit none
! This program glues together the ts files from UCLA
! Compile with ifort -o reducets reducets.f90 -L/sw/sles9-x64/netcdf-3.6.2-intel/lib -lnetcdf
! -I/sw/sles9-x64/netcdf-3.6.2-intel/include on tornado
! Assume number of time steps nt smaller than 10000 and number of variables nv smaller than 35 ! otherwise need to be changes
        
        integer, parameter :: nv=50,nma=7,nmi=2,nsu=28
        integer :: nt
        character(100) stem,nm,nm2
        character(20) pref,name,dimname
        character(4) xsuf,ysuf
        character(1) dum1
        character(2) dum2
        integer nx,ny
        integer status,ncid,ndim,nvar,ncidw,timeid,varid
        integer flag(nv)
        character(len=10) maxnms(nma)
        character(len=10) minnms(nmi)
        character(len=10) sumnms(nsu)
        character(len=100) namout(nv),lname(nv),uname(nv)
        real, dimension(:,:), allocatable ::  var,varout
        integer, dimension(:,:), allocatable :: cnt
        integer i,j,k,kk,l, ntr


!0. Some preparations
!--------------------
!* maxnms contains the names where we take max, minnms where we take min, summns where we take !*sum; if nothing specified will take the mean; list can be extended if needed

!* Ask to get filenames
        call get_command_argument(1,stem)
        call get_command_argument(2,nx)
        call get_command_argument(3,ny)
!         print*,'Directory where files are'
!         read*,dirin
!         print*,'Files prefix'
!         read*,pref
!         print*,'Processors in x'
!         read*,nx
!         print*,'Processors in y'
!         read*,ny
! Open the first file to check the amoutn of timesteps; allocate the arrays
        nm=trim(stem)//".ts.00000000.nc"
        print*,trim(nm)
        status=nf90_open(nm,nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_inquire(ncid,unlimitedDimID=timeid)
        status=nf90_inquire_dimension(ncid,timeid,len=nt)
        status=nf90_close(ncid)
        print *, 'Number of timesteps ', nt
        allocate(var(nv,nt), varout(nv,nt), cnt(nv,nt))
        varout(:,:)=0.
        cnt(:,:)=0
        var(:,:)=0.

        maxnms(1)="cfl       "
        maxnms(2)="maxdiv    "
        maxnms(3)="wmax      "
        maxnms(4)="lmax     "
        maxnms(5)="bflxmx    "
        maxnms(6)="bflxrmx   "
        maxnms(7)="precip_m  "
        minnms(1)="bflxmn    "
        minnms(2)="bflxrmn   "
        sumnms(1)="wr_cs1    "
        sumnms(2)="wr_cs2    "
        sumnms(3)="wv_cs1    "
        sumnms(4)="wv_cs2    "
        sumnms(5)="wt_cs1    "
        sumnms(6)="wt_cs2    "
        sumnms(7)="rt_cs1    "
        sumnms(8)="rt_cs2    "
        sumnms(9)="rl_cs1    "
        sumnms(10)="rl_cs2    "
        sumnms(11)="tv_cs1    "
        sumnms(12)="tv_cs2    "
        sumnms(13)="tl_cs1    "
        sumnms(14)="tl_cs2    "
        sumnms(15)="w_cs1     "
        sumnms(16)="w_cs2     "
        sumnms(17)="cnt_cs1   "
        sumnms(18)="cnt_cs2   "
        sumnms(19)="nrain     "
        sumnms(20)="nrcnt     "
        sumnms(21)="hst_d1    "
        sumnms(22)="hst_d2    "
        sumnms(23)="hst_d3    "
        sumnms(24)="hst_d4    "
        sumnms(25)="dissum    "
        sumnms(26)="discnt    "
        sumnms(27)="dissum2   "
        sumnms(28)="discnt2   "



!1. Open each file and read variable
!----------------------------------

!*Loop over files and reconstruct filename
        do j=1,ny
         if (j-1.lt.10) then
          write(dum1,'(I1.1)') j-1
          ysuf="000"//dum1
         else
          write(dum2,'(I2.1)') j-1
          ysuf="00"//dum2
         end if
         do i=1,nx
          if (i-1.lt.10) then
           write(dum1,'(I1.1)') i-1
           xsuf="000"//dum1
          else
           write(dum2,'(I2.1)') i-1
           xsuf="00"//dum2
          end if
          nm=trim(stem)//".ts."//ysuf//xsuf//".nc"
          print*,trim(nm)
           
!*Open file
          status=nf90_open(nm,nf90_nowrite,ncid)
           if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!*Get Variable names, dims, attributes. Do this only once at it is the same for all files  
!*Also set the flag according to compute mean (0), max (2), min (3), sum (1)
          if ((i.eq.1).and.(j.eq.1)) then
           flag(:)=0
           status=nf90_inquire(ncid,ndim,nvar)
            if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           status=nf90_inquire_dimension(ncid,1,dimname,ntr)
            if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           print*,'There are ',ntr,' time steps and ',nvar, 'variables in the file'
           do k=1,nvar
            status=nf90_inquire_variable(ncid,k,name)
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            namout(k)=name
!*Assume the attribute names are longname and units
            status=nf90_get_att(ncid,k,'longname',lname(k))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_att(ncid,k,'units',uname(k))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            do kk=1,nma
              if (trim(name).eq.trim(maxnms(kk))) flag(k)=2
            end do
            do kk=1,nmi
              if (trim(name).eq.trim(minnms(kk))) flag(k)=3
            end do
            do kk=1,nsu
              if (trim(name).eq.trim(sumnms(kk))) flag(k)=1
            end do
           end do
          end if

!*Read variable in var
           do k=1,nvar
            status=nf90_get_var(ncid,k,var(k,1:ntr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           end do

!*Close file
          status=nf90_close(ncid)
           if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!2. Compute max, min, sum or average
!----------------------------------
          do k=1,nvar
           if (flag(k).le.1) then
            do l=1,ntr
             if (var(k,l).gt.-999.) then
              varout(k,l)=varout(k,l)+var(k,l)
              cnt(k,l)=cnt(k,l)+1
             end if
            end do
           else
            if (flag(k).eq.2) then
             do l=1,ntr
              if ((varout(k,l).eq.0).or.(var(k,l).gt.varout(k,l))) varout(k,l)=var(k,l)
             end do
            else
             do l=1,ntr
              if ((varout(k,l).eq.0).or.((var(k,l).gt.-999).and.(var(k,l).lt.varout(k,l)))) varout(k,l)=var(k,l)
             end do
            endif
           end if
          end do

         end do
        end do
!*End of looping over the files
  
       do k=1,nvar
         if (flag(k).eq.0) varout(k,:)=varout(k,:)/cnt(k,:)
      end do

!3 Write new variables in *.ts.nc
!---------------------------------

!* Enter define mode
      nm2=trim(stem)//".ts.nc"
      print*,nm2
      status=nf90_create(nm2, nf90_CLOBBER,ncidw)
       if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!* Define dimenstions
      status=nf90_def_dim(ncidw,'time',nf90_unlimited,timeid)
       if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!* Define and store variable
      do k=1,nvar
       status=nf90_def_var(ncidw,trim(namout(k)),NF90_float,(/timeid/),varid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_put_att(ncidw,varid,'longname',trim(lname(k)))
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_put_att(ncidw,varid,'_FillValue',-999.0)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_put_att(ncidw,varid,'units',trim(uname(k)))
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_enddef(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_put_var(ncidw,varid,varout(k,1:ntr))
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_redef(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end do
       status=nf90_enddef(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!* Close file
       status=nf90_close(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)

        end
