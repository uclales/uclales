begin

  nmb = 21600
  fname = "../bin/g9.dws"
  sffxps  = ".ps."
  sffxts  = ".ts."
  nfiles  = 16

; condencinsing all ts.nc files in one

  print ("cp "+fname+sffxts+"0000.nc "+fname+sffxts+"nc")
  system ("cp "+fname+sffxts+"0000.nc "+fname+sffxts+"nc")

  base_file = addfile(fname+sffxts+"nc","w")
  names = getfilevarnames (base_file)
 
  maxnms = (/"w2","precip_m","wmax","mxdiv","cfl"/)
  sumnms = (/"zi_bar","zi_var","zb_bar","zb_var","lwp_bar","lwp_var", \
           "cfrac","ndrop","tke","wstar","shf_bar","lhf_bar","precip", \
           "precip_v","ustar"/)
  varnms = (/"zi_var","zb_var","lwp_var","precip_v"/)
  avgnms = (/"zi_bar","zb_bar","lwp_bar","precip"/)

  do ij=1,nfiles-1
    if (ij .ge. 10) then
      test_file = addfile(fname+sffxts+"00"+ij+".nc","r")
    else 
      test_file = addfile(fname+sffxts+"000"+ij+".nc","r")
    end if
    do i = 0, dimsizes(names) -1
      notdone = True
      do j = 0, dimsizes(maxnms) -1 
        if (notdone .and. names(i).eq.maxnms(j)) then
          do k = 0, dimsizes(base_file->time)-1
            base_file->$names(i)$(k) = \
                (/max((/base_file->$names(i)$(k),test_file->$names(i)$(k)/))/)
          end do
          notdone = False
        end if
      end do
      if (notdone) then
        base_file->$names(i)$ = (/base_file->$names(i)$+test_file->$names(i)$/)
        notdone = False
      end if
    end do 
    delete(test_file)
  end do
  
  xx = 1.0 * nfiles
  do i = 0, dimsizes(names) -1
    notdone = True
    do j = 0, dimsizes(maxnms) -1 
      if (notdone .and. names(i).eq.maxnms(j)) then
        notdone = False
      end if
    end do
    if (notdone) then       
      base_file->$names(i)$ = (/base_file->$names(i)$/xx/)
    end if
  end do

; calculating variances and wstar
  do i = 0, dimsizes(names) -1
    notdone = True
    do j = 0, dimsizes(varnms) -1
      if (notdone .and. names(i).eq.varnms(j)) then
        base_file->$names(i)$ = (/base_file->$names(i)$ - \
            (base_file->$avgnms(j)$)*(base_file->$avgnms(j)$)/)
        notdone = False
      end if
    end do
    if (notdone .and. names(i).eq."wstar") then
      base_file->$names(i)$ = (/(base_file->$names(i)$)^(1./3.)/)
      notdone = False
    end if
  end do

  delete (base_file)
  delete (names)
  delete (maxnms)
  delete (sumnms)
  delete (varnms)
  delete (avgnms)

; condensing all ps.nc files in one

  print ("cp "+fname+sffxps+"0000.nc "+fname+sffxps+"nc")
  system ("cp "+fname+sffxps+"0000.nc "+fname+sffxps+"nc")

  base_file = addfile(fname+sffxps+"nc","w")
  names = getfilevarnames (base_file)
 
  sumnms = (/"dn0","u","v","thetal","rt","rl","rr","cfrac", \
         "ndrop_cl","ss","u_var","v_var","w_var","w_3","theta_2", \
         "rt_var","rl_var","rad_flx","precip","tot_tw","sfs_tw","tot_rw", \
         "sfs_rw","tot_uw","sfs_uw","tot_vw","sfs_vw","sgs_pr","sfs_tke", \
         "boy_prd","sfs_boy","shr_prd","sfs_shr","trans","diss"/)
  varnms = (/"u_var","v_var","theta_2","rt_var","rl_var"/)
  avgnms = (/"u","v","thetal","rt","rl"/)
  bsst = (/5.,-5.5,289.,0.,0./)

  do j=1,nfiles-1
    if (j .ge. 10) then
      test_file = addfile(fname+sffxps+"00"+j+".nc","r")
    else 
      test_file = addfile(fname+sffxps+"000"+j+".nc","r")
    end if
    do i = 0, dimsizes(names) -1 
      base_file->$names(i)$ = (/base_file->$names(i)$ + test_file->$names(i)$/)
    end do 
    delete(test_file)
  end do
  
  xx = 1.0 * nfiles
  do i = 0, dimsizes(names) -1
    notdone = True
    if (notdone) then
      base_file->$names(i)$ = (/base_file->$names(i)$/xx/)
    end if
  end do 

; calculating variances
  do i = 0, dimsizes(names) -1
    notdone = True
    do j = 0, dimsizes(varnms) -1
      if (notdone .and. names(i).eq.varnms(j)) then
        base_file->$names(i)$ = (/base_file->$names(i)$ - \
          (base_file->$avgnms(j)$-bsst(j))*(base_file->$avgnms(j)$-bsst(j))/)
        notdone = False
      end if
    end do
  end do

end
