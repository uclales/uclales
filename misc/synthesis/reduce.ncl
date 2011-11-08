load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
;
; this routine uses a function called barbell which I wrote and put in the file
; called myfuncs.load.  It makes an error bar like symbol to overlay on plots.
;

begin

  fname = "ccbl_256.25.06.15.100"
  sffxps  = ".ps."
  sffxts  = ".ts."
  nfiles  = 64

  maxnms = (/"cfl","maxdiv","wmax","rlmax","bflxmx","bflxrmx","precip_m"/)
  minnms = (/"bflxmn","bflxrmn"/)
  sumnms = (/"wr_cs1","wr_cs2","wv_cs1","wv_cs2","wt_cs1","wt_cs2" \
	    ,"rt_cs1","rt_cs2","rl_cs1","rl_cs2","tv_cs1","tv_cs2" \
	    ,"tl_cs1","tl_cs2", "w_cs1", "w_cs2","cnt_cs1","cnt_cs2" /)

  print ("cp "+fname+sffxts+"0000.nc "+fname+sffxts+"nc")
  system ("cp "+fname+sffxts+"0000.nc "+fname+sffxts+"nc")
  print ("cp "+fname+sffxps+"0000.nc "+fname+sffxps+"nc")
  system ("cp "+fname+sffxps+"0000.nc "+fname+sffxps+"nc")

  base_file = addfile(fname+sffxts+"nc","w")
  names = getfilevarnames (base_file)
 
  do j=1,nfiles-1
     print ("processing ts file "+j)
     if (j .ge. 10) then
        test_file = addfile(fname+sffxts+"00"+j+".nc","r")
     else 
        test_file = addfile(fname+sffxts+"000"+j+".nc","r")
     end if
     do i = 0, dimsizes(names) -1
	notdone = True
        do k = 0, dimsizes(maxnms) -1 
           if (notdone .and. names(i).eq.maxnms(k)) then
  	      do k=0,dimsizes(base_file->time)-1
                 base_file->$names(i)$(k) = \ 
	         max( (/ base_file->$names(i)$(k),test_file->$names(i)$(k) /) )
              end do
              notdone = False
           end if
        end do
        do k = 0, dimsizes(minnms) -1 
           if (notdone .and. names(i).eq.minnms(k)) then
  	      do k=0,dimsizes(base_file->time)-1
                 base_file->$names(i)$(k) = \ 
	         min( (/ base_file->$names(i)$(k),test_file->$names(i)$(k) /) )
              end do
              notdone = False
           end if
        end do
        if (notdone) then
           base_file->$names(i)$ = base_file->$names(i)$+test_file->$names(i)$
	end if
     end do 
     delete(test_file)
  end do

  print ("final processing of ts files")
  xx = 1.0 * nfiles
  do i = 0, dimsizes(names) -1
     notdone = True
     do k = 0, dimsizes(maxnms) -1 
        if (notdone .and. names(i).eq.maxnms(k)) then
           notdone = False
        end if
     end do
     do k = 0, dimsizes(minnms) -1 
        if (notdone .and. names(i).eq.minnms(k)) then
           notdone = False
        end if
     end do
     if (notdone) then       
 	base_file->$names(i)$ = base_file->$names(i)$/xx
     end if
  end do

  delete (names)

  base_file = addfile(fname+sffxps+"nc","w")
  names = getfilevarnames (base_file)
 
  do j=1,nfiles-1
     print ("processing ps file "+j)
     if (j .ge. 10) then
        test_file = addfile(fname+sffxps+"00"+j+".nc","r")
     else 
        test_file = addfile(fname+sffxps+"000"+j+".nc","r")
     end if
     do i = 0, dimsizes(names) -1 
           base_file->$names(i)$=base_file->$names(i)$ + test_file->$names(i)$
     end do 
     delete(test_file)
  end do

  xx = 1.0 * nfiles
  do i = 0, dimsizes(names) -1
     notdone = True
     do k = 0, dimsizes(sumnms) -1 
        if (notdone .and. names(i).eq.sumnms(k)) then
           notdone = False
        end if
     end do
     if (notdone) then
        base_file->$names(i)$ = base_file->$names(i)$/xx
     end if
  end do 

end




