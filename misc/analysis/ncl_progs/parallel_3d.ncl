;#########################################################################;
;                                                                         ;
;  File:        parallel.ncl                                              ;
;                                                                         ;
;  Author:      Irina Sandu 
;                                                                         ;
;  Date:        November 2009                                              ;
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

function stitch_xy(file_stem:string,fld_name:string,itime:integer)

begin

  ; data
  fi = addfile(file_stem+".00000000.nc","r")

  nprocs_t = fi@NPROCS
  nprocs=floattointeger(sqrt(nprocs_t))
  nxp = (dimsizes(fi->xt))*nprocs
  nyp = (dimsizes(fi->yt))*nprocs
  nzp = (dimsizes(fi->zt))

  print ("Stitching "+nyp+" by "+nxp+" x-y field from "+nprocs_t+" processors")
  print ("for a file with "+nzp+" levels")
  
  fld = new ( (/nyp,nxp,nzp/), "float")
  fld@longname = fi->$fld_name$@longname
  fld@units = fi->$fld_name$@units
  do n=0,nprocs-1
  nn=0 
print("n"+n)
  if (n.lt.10) then
      fname =file_stem+".000"+n+"000"+nn+".nc" 
    else
      fname =file_stem+".00"+n+"000"+nn+".nc" 
    end if
    fi = addfile(fname,"r")
    nxp = dimsizes(fi->xt)
    ia = n*nxp
    ib = ia + nxp - 1
    nyp = dimsizes(fi->yt)
    ja = nn*nyp
    jb = ja + nyp - 1
print("ia ib "+ia+" "+ib)
print("ja jb "+ja+" "+jb)
    fld(ja:jb,ia:ib,:) = fi->$fld_name$(itime,:,:,:)
   do nn=1,nprocs-1
print("nn "+nn)
    if (nn.lt.10) then
     if (n.lt.10) then
      fname =file_stem+".000"+n+"000"+nn+".nc" 
     else
      fname =file_stem+".00"+n+"000"+nn+".nc" 
     end if
    else
     if (n.lt.10) then
      fname =file_stem+".000"+n+"00"+nn+".nc" 
     else
      fname =file_stem+".00"+n+"00"+nn+".nc" 
     end if
    end if
    fi = addfile(fname,"r")
    nyp = dimsizes(fi->yt)
    ja = nn*nyp
    jb = ja + nyp - 1
print("ja jb "+ja+" "+jb)
    fld(ja:jb,ia:ib,:) = fi->$fld_name$(itime,:,:,:)
  end do
  end do

  return fld

end

;---------------------------------------------------------------------------;


