;#########################################################################;
;                                                                         ;
;  File:        parallel.ncl                                              ;
;                                                                         ;
;  Author:      Bjorn Stevens
;                                                                         ;
;  Date:        June 17 2006                                              ;
;                                                                         ;
;                                                                         ;
;                                                                         ;
;  This file contains ncl functions for reconstructing fields distributed ;
;  across many files as a result of having been cacluated over many       ;
;  processors                                                             ;
;                                                                         ;
;  The function stitch_xy reconstructs a field over x & y given a file    ;
;      stem and a field name                                              ;
;                                                                         ;
;                                                                         ;
;#########################################################################;



;---------------------------------------------------------------------------;

; function stitch_xy fills a field over x & y

undef("stitch_xy")

function stitch_xy(file_stem:string,fld_name:string,itime:integer,k:integer)

begin

  ; data
  fi = addfile(file_stem+".0000.nc","r")

  nprocs = fi@NPROCS
  nxp = dimsizes(fi->xt)-4
  nyp = (dimsizes(fi->yt)-4)*nprocs

  print ("Stitching "+nyp+" by "+nxp+" x-y field from "+nprocs+" processors")
  
  fld = new ( (/nyp,nxp/), "float")
  fld@longname = fi->$fld_name$@longname
  fld@units = fi->$fld_name$@units
  do n=0,nprocs-1
    if (n.lt.10) then
      fname =file_stem+".000"+n+".nc" 
    else
      fname =file_stem+".00"+n+".nc" 
    end if
    fi = addfile(fname,"r")
    nyp = dimsizes(fi->yt)-4
    ja = n*nyp
    jb = ja + nyp - 1
    fld(ja:jb,:) = fi->$fld_name$(itime,2:nyp+1,2:nxp+1,k)
  end do

  return fld

end

;---------------------------------------------------------------------------;


