load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
;
; this routine uses a function called barbell which I wrote and put in the file
; called myfuncs.load.  It makes an error bar like symbol to overlay on plots.
;

begin


  maxnms = (/"cfl","maxdiv","wmax","rlmax","bflxmx","bflxrmx","precip_m"/)
  minnms = (/"bflxmn","bflxrmn"/)
  sumnms = (/"wr_cs1","wr_cs2","wv_cs1","wv_cs2","wt_cs1","wt_cs2" \
	    ,"rt_cs1","rt_cs2","rl_cs1","rl_cs2","tv_cs1","tv_cs2" \
	    ,"tl_cs1","tl_cs2", "w_cs1", "w_cs2","cnt_cs1","cnt_cs2" /)

  stem = "rf01f"
  sffxps  = ".ps."
  sffxts  = ".ts."
  nxpart  = True
  nfiles  = 64

  if (nxpart) then
    nxfiles  = sqrt(nfiles)
    nyfiles  = nxfiles
  else 
    nxfiles = 1
    nyfiles = nfiles
  end if
  print ("cp "+stem+sffxts+"00000000.nc "+stem+sffxts+"nc")
  system("cp "+stem+sffxts+"00000000.nc "+stem+sffxts+"nc")
  print ("cp "+stem+sffxps+"00000000.nc "+stem+sffxps+"nc")
  system("cp "+stem+sffxps+"00000000.nc "+stem+sffxps+"nc")


  base_file = addfile(stem+sffxts+"nc","w")
  names = getfilevarnames (base_file)
 
  do i=0,nxfiles-1
     if (i .ge. 1000) then
       istem = stem + sffxts + i 
     else
       if (i .ge. 100) then
         istem = stem + sffxts + "0" + i
       else
         if(i .ge. 10) then
           istem = stem + sffxts + "00" + i
         else 
           istem = stem + sffxts + "000" + i
         end if
       end if
     end if
     if (.not.nxpart) then
       fname = stem + sffxts
     end if

     do j=max((/0,1-i/)),nyfiles-1
       if (j .ge. 1000) then
         fname = istem + j +".nc"
       else
         if (j .ge. 100) then
           fname = istem + "0"+j+".nc"
         else
           if(j .ge. 10) then
             fname = istem + "00"+j+".nc"
           else 
             fname = istem + "000"+j+".nc"
           end if
         end if
       end if

       print ("processing "+fname)
       test_file = addfile(fname,"r")
       do n = 0, dimsizes(names) -1
         notdone = True
         do k = 0, dimsizes(maxnms) -1 
           if (notdone .and. names(n).eq.maxnms(k)) then
             do k=0,dimsizes(base_file->time)-1
               base_file->$names(n)$(k) = \ 
               max( (/ base_file->$names(n)$(k),test_file->$names(n)$(k) /) )
             end do
             notdone = False
           end if
         end do
         do k = 0, dimsizes(minnms) -1 
           if (notdone .and. names(n).eq.minnms(k)) then
             do k=0,dimsizes(base_file->time)-1
               base_file->$names(n)$(k) = \ 
               min( (/ base_file->$names(n)$(k),test_file->$names(n)$(k) /) )
             end do
             notdone = False
           end if
         end do
         if (notdone) then
           base_file->$names(n)$ = base_file->$names(n)$+test_file->$names(n)$
         end if
       end do 
       delete(test_file)
     end do
   end do

   print ("final processing of ts files")
   xx = 1.0 * nfiles
   do n = 0, dimsizes(names) -1
     if (names(n).eq."time") then
       print ("xx "+xx+", "+base_file->$names(n)$)
     end if
     notdone = True
     do k = 0, dimsizes(maxnms) -1 
       if (notdone .and. names(n).eq.maxnms(k)) then
         notdone = False
       end if
     end do
     do k = 0, dimsizes(minnms) -1 
       if (notdone .and. names(n).eq.minnms(k)) then
         notdone = False
       end if
     end do
     if (notdone) then       
       base_file->$names(n)$ = base_file->$names(n)$/xx
     end if
   end do
   delete (names)
  
   base_file = addfile(stem+sffxps+"nc","w")
   names = getfilevarnames (base_file)
  
   do i=0,nxfiles-1
     if (i .ge. 1000) then
       istem = stem + sffxps + i
     else
       if (i .ge. 100) then
         istem = stem + sffxps + "0" + i
       else
         if(i .ge. 10) then
           istem = stem + sffxps + "00" + i
         else 
           istem = stem + sffxps + "000" + i
         end if
       end if
     end if

     if (.not.nxpart) then
       fname = stem + sffxts
     end if

     do j=max((/0,1-i/)),nyfiles-1
       if (j .ge. 1000) then
         fname = istem + j +".nc"
       else
         if (j .ge. 100) then
           fname = istem + "0"+j+".nc"
         else
           if(j .ge. 10) then
             fname = istem + "00"+j+".nc"
           else 
             fname = istem + "000"+j+".nc"
           end if
         end if
       end if

       test_file = addfile(fname,"r")
       print ("processing "+fname)

       do n = 0, dimsizes(names) -1 
         base_file->$names(n)$=base_file->$names(n)$ + test_file->$names(n)$
       end do 
       delete(test_file)
     end do
   end do
  
   xx = 1.0 * nfiles
   do n = 0, dimsizes(names) -1
     notdone = True
     do k = 0, dimsizes(sumnms) -1 
       if (notdone .and. names(n).eq.sumnms(k)) then
         notdone = False
       end if
     end do
     if (notdone) then
       base_file->$names(n)$ = base_file->$names(n)$/xx
     end if
   end do 
  
end


