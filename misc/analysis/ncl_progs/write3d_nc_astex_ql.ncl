;#########################################################################;
;                                                                         ;
;  File:        write3d_nc.ncl                                               ;
;                                                                         ;
;  Author:      Irina Sandu                                             ;
;                                                                         ;
;  Date:        November 09 2009                                              ;
;                                                                         ;
;                                                                         ;
;                                                                         ;
;  this  script reconstructs the 3d fileds from the files written on each 
; processor and writes them in a nc file
;                                                                         ;
;                                                                         ;
;#########################################################################;

;---------------------------------------------------------------------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "parallel_3d.ncl"


begin


  fname = "tr5m"
  stem = fname
  fi = addfile(stem+".00000000.nc","r")

  nprocs_t = fi@NPROCS
  nprocs=floattointeger(sqrt(nprocs_t))
  nxp = (dimsizes(fi->xt))*nprocs
  nyp = (dimsizes(fi->yt))*nprocs


  ;
  ; dimensions
  ;
  ntimes = dimsizes(fi->time)
print(ntimes)

 zt=fi->zt(:)
 zm=fi->zm(:)
 l_ini=fi->l(:,:,:,:)
printVarSummary(l_ini)

;initial time and nb of time steps required
req_time=82
in_time=0

time=fi->time(in_time:req_time+in_time-1)
print(time)

ql_3d=new((/req_time,nyp,nxp,dimsizes(zt)/), typeof(l_ini))
  ql_3d!0="time"
  ql_3d!3="zt"
 ql_3d&time=time
 ql_3d&zt=zt

  ;
  ; constants 
  do n=in_time,req_time+in_time-1
  

    ql = stitch_xy(stem,"l",n)
printVarSummary(ql)
     yt=ql&yt
     xt=ql&xt

    ql_3d(n-in_time,:,:,:)=ql
  ql_3d!1="yt"
  ql_3d!2="xt"
 ql_3d&yt=yt
 ql_3d&xt=xt
printVarSummary(ql_3d)

   end do


printVarSummary(ql_3d)
printMinMax(ql_3d,True)

;writes the cf in a netcdf file

; associated meta data: (a) coordinate variables  lat, lon       
; and (b) attributes

    filo = "tr5m_3d.nc"             ; Output file
       system("/bin/rm " + filo)    ; remove if exists
    fout  = addfile (filo, "c")  ; open output file

; explicitly declare file definition mode. Improve efficiency.

    setfileoption(fout,"DefineMode",True)
    
; create global attributes of the file
        fAtt               = True            ; assign file attributes
	fAtt@title         = "3d fields"  
	fAtt@source_file   = "trans1" 
	fAtt@Conventions   = "None"   
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; copy file attributes    

; predefine the coordinate variables and their dimensionality
        dimNames = (/"time","xt","xm","yt","zt","zm"/)  
        ;dimNames = (/"time","xt","xm","yt","ym","zt","zm"/)  
        dimSizes = (/dimsizes(time),nxp,nxp,nyp,dimsizes(zt),dimsizes(zm)/)  
        ;dimSizes = (/dimsizes(time),nxp,nxp,nyp,nyp,dimsizes(zt),dimsizes(zm)/)  
	dimUnlim = (/ True,False,False,False,False,False/)   
	;dimUnlim = (/ True,False,False,False,False,False,False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)


; predefine the the dimensionality of the variables to be written out
       filevardef(fout, "time"  ,typeof(time),getvardims(time))                       
       filevardef(fout, "xt"  ,typeof(xt),getvardims(xt))                          
       filevardef(fout, "yt"  ,typeof(yt),getvardims(yt))                          
       filevardef(fout, "zt"  ,typeof(zt),getvardims(zt))                          
       filevardef(fout, "zm"  ,typeof(zm),getvardims(zm))                          
       filevardef(fout, "ql"  ,typeof(ql_3d),getvardims(ql_3d))                          


       filevarattdef(fout,"time",time)                           ; copy T attributes
       filevarattdef(fout,"xt",xt)                           ; copy T attributes
       filevarattdef(fout,"yt",yt)                           ; copy T attributes
       filevarattdef(fout,"zt",zt)                           ; copy T attributes
       filevarattdef(fout,"zm",zm)                           ; copy T attributes
       filevarattdef(fout,"ql",ql_3d)                           ; copy T attributes


; explicitly exit file definition mode. **NOT REQUIRED**
        setfileoption(fout,"DefineMode",False)
	
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.

       fout->time    = (/time/)
       fout->xt    = (/xt/)
       fout->yt    = (/yt/)
       fout->zt    = (/zt/)
       fout->zm    = (/zm/)
      fout->ql    = (/ql_3d/)

end
