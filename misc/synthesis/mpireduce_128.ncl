load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
;
; Routine for combining netcdf files for each of 'nfiles' subdomains (processors) into one file.
;
; Louise: I added nrain, nrcnt, hist_d1-4 to the variables in the time series netcdf file that need to be summed and not averaged.
; Malte: I added if statements in line 46-50 and 119-123 to allow nfiles bigger than 100 (l > 9). Loop indices must be put in manually (128proc: 8 x 16)

begin

  dirin = ""
  dirin2 = dirin+""
  fname = "rico_04_l2"

  sffxps  = ".ps."
  sffxts  = ".ts."
  nfiles  = 128

  maxnms = (/"cfl","maxdiv","wmax","rlmax","bflxmx","bflxrmx","precip_m"/)
  minnms = (/"bflxmn","bflxrmn"/)
  sumnms = (/"wr_cs1","wr_cs2","wv_cs1","wv_cs2","wt_cs1","wt_cs2" \
	    ,"rt_cs1","rt_cs2","rl_cs1","rl_cs2","tv_cs1","tv_cs2" \
	    ,"tl_cs1","tl_cs2", "w_cs1", "w_cs2","cnt_cs1","cnt_cs2","nrain","nrcnt","hst_d1","hst_d2","hst_d3","hst_d4","dissum","discnt","dissum2","discnt2"/)

 ; system("rm "+(/dirin+fname+sffxts+"nc"/))
 ; system("rm "+(/dirin+fname+sffxps+"nc"/))

  print ("cp "+dirin2+fname+sffxts+"00000000.nc "+dirin+fname+sffxts+"nc")
  system ("cp "+dirin2+fname+sffxts+"00000000.nc "+dirin+fname+sffxts+"nc")
  print ("cp "+dirin2+fname+sffxps+"00000000.nc "+dirin+fname+sffxps+"nc")
  system ("cp "+dirin2+fname+sffxps+"00000000.nc "+dirin+fname+sffxps+"nc")

  base_file = addfile(dirin+fname+sffxts+"nc","w")
  names = getfilevarnames (base_file)
 
 ; nf = sqrt(nfiles) 

  ;do j=0,nf-1
  ; do l = 0,nf-1 
  do j=0,7
   do l = 0,15 
     if .not.(j.eq.0.and.l.eq.0) then	
      print ("processing ts file "+j+" - "+l)
      
       if (l.lt.10) then
          test_file = addfile(dirin2+fname+sffxts+"000"+j+"000"+l+".nc","r")
       else
          test_file = addfile(dirin2+fname+sffxts+"000"+j+"00"+l+".nc","r")
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
     
     end if
   end do
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
     do k = 0, dimsizes(sumnms) -1 
        if (notdone .and. names(i).eq.sumnms(k)) then
           notdone = False
        end if
     end do
     if (notdone) then       
 	base_file->$names(i)$ = base_file->$names(i)$/xx
     end if
  end do

  delete (names)

  base_file = addfile(dirin+fname+sffxps+"nc","w")
  names = getfilevarnames (base_file)
 
  ;do j=0,nf-1
  ;do l = 0,nf-1 
  do j=0,7
  do l = 0,15
  
       if .not.(j.eq.0.and.l.eq.0) then	
       print ("processing ps file "+j+" - "+l)
      
       if (l.lt.10) then
          test_file = addfile(dirin2+fname+sffxps+"000"+j+"000"+l+".nc","r")
       else
          test_file = addfile(dirin2+fname+sffxps+"000"+j+"00"+l+".nc","r")
       end if

       do i = 0, dimsizes(names) -1 
           base_file->$names(i)$=base_file->$names(i)$ + test_file->$names(i)$
       end do 
       delete(test_file)
     end if
   end do
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




