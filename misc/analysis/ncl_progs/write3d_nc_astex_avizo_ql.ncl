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


begin


  fname = "tr5m_3d.nc"

a=addfile(fname,"r")

time=a->time
time@units="seconds since 1992-06-14 00:00"
latitude_t=a->xt
longitude_t=a->yt
;longitude_m=a->ym
levels_m=a->zm(:)
levels_t=a->zt(:)

ql_3d=a->ql(:,:,:,:)
;qt_3d=a->qt
;qr_3d=a->qr
;w_3d=a->w
;thetal_3d=a->thetal

printVarSummary(ql_3d)
;printVarSummary(qt_3d)
;printVarSummary(qr_3d)
;printVarSummary(w_3d)
;printVarSummary(thetal_3d)


latitude_t@units="degrees_north"
latitude_t@longname="latitude"
longitude_t@units="degrees_east"
longitude_t@longname="longitude"
longitude_t@standard_name="longitude"
;longitude_m@units="degrees_east"
;longitude_m@longname="longitude"
;longitude_m@standard_name="longitude"

latitude_t!0="latitude_t"
latitude_t&latitude_t=latitude_t


longitude_t!0="longitude_t"
longitude_t&longitude_t=longitude_t

;longitude_m!0="longitude_m"
;longitude_m&longitude_m=longitude_m

levels_m@longname="levels"
levels_t@longname="levels"
levels_m!0="levels_m"
levels_m&levels_m=levels_m
levels_t!0="levels_t"
levels_t&levels_t=levels_t

printVarSummary(latitude_t)
printVarSummary(longitude_t)
printVarSummary(levels_m)
printVarSummary(levels_t)

ql_3d!3="levels_t"
ql_3d&levels_t=levels_t
ql_3d!2="latitude_t"
ql_3d&latitude_t=latitude_t
ql_3d!1="longitude_t"
ql_3d&longitude_t=longitude_t
printVarSummary(ql_3d)

ql_bis=ql_3d(time|:,levels_t|:,latitude_t|:,longitude_t|:)
printVarSummary(ql_bis)


;writes the cf in a netcdf file

; associated meta data: (a) coordinate variables  lat, lon       
; and (b) attributes

    filo = "tr5m_3d_avizo.nc"             ; Output file
       system("/bin/rm " + filo)    ; remove if exists
    fout  = addfile (filo, "c")  ; open output file

; explicitly declare file definition mode. Improve efficiency.

   ; setfileoption(fout,"DefineMode",True)
   
   ; make time and UNLIMITED dimension      ; recommended  for most applications
       filedimdef(fout,"time",-1,True) 

       fout->time    = time
       fout->latitude_t    = latitude_t
       fout->longitude_t    = longitude_t
     ;  fout->longitude_m    = longitude_m
       fout->levels_t    = levels_t
       fout->levels_m    = levels_m
       fout->ql    = ql_bis


; create global attributes of the file
        fAtt               = True            ; assign file attributes
	fAtt@title         = "3d fields"  
	fAtt@source_file   = "trans1" 
	fAtt@Conventions   = "CF-1.0"   
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; copy file attributes    


end
