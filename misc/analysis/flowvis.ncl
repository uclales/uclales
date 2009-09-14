;#########################################################################;
;                                                                         ;
;  File:        flowvis.ncl                                               ;
;                                                                         ;
;  Author:      Bjorn Stevens                                             ;
;                                                                         ;
;  Date:        July 28 2006                                              ;
;                                                                         ;
;                                                                         ;
;                                                                         ;
;  This file contains ncl functions for regridding and visualing flows    ;
;  it is in somewhat of a state of disrepair (usage) but useful things    ;
;  may be extractable from it, for instance, the newshape function which  ;
;  puts variable vertical grid data on a fixed grid.                      ;
;                                                                         ;
;                                                                         ;
;#########################################################################;

;---------------------------------------------------------------------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "./parallel.ncl"

undef ("new_shape")

function new_shape(x,z,jpt,npt,fi,vname)
local dims, nz, nx, dmy, fld_in, fld, jj
begin

  dims = dimsizes(fi->$vname$)
  
  nz = dimsizes(z)
  nx = dimsizes(x)
  fld = new ( (/nz,nx/), "float")
  
  dmy  = fi->$vname$
  jj = jpt + 2
  if (vname.eq."w")then
     fld_in = dmy(zm|:,xt|:,yt|jj,time|npt)
  else 
     fld_in = dmy(zt|:,xt|:,yt|jj,time|npt)
  end if
  fld=linint2(fi->xt,fi->zt,fld_in,False,x,z,0)

  fld!0 = "z"
  fld&z = z
  fld!1 = "x"
  fld&x = x
  return (fld)

end


begin

  pi = 2.*asin(1.)

  fname = "anim.25a"
  stem = "../runs/"+fname+"/"+fname
  fi = addfile(stem+".0000.nc","r")

  nprocs = fi@NPROCS

  zz = 30.
  kpt = closest_val(zz,fi->zt)

  wks = gsn_open_wks("ncgm",fname+round(zz,1))
  gsn_define_colormap (wks,"nrl_sirkes")
  
  cnres = True
  cnres@gsnFrame = False
  cnres@gsnDraw = False
  cnres@gsnSpreadColors     = True           ; use full range of colors
  cnres@gsnSpreadColorStart = 2
  cnres@gsnSpreadColorEnd   = -2
  cnres@lbLabelStride       = 4
  cnres@cnLinesOn = False
  cnres@cnLineLabelsOn = False
  cnres@cnRasterModeOn = True
  cnres@cnFillOn = True

  ;
  ; dimensions
  ;
  ntimes = dimsizes(fi->time)
  plots = new ( (/10/), "graphic")

  y = fi->xt(4:)
  x = fi->xt(4:)
  yn = dimsizes(y)-1
  xn = dimsizes(x)-1

  ;
  ; constants 
  do n=0,ntimes-1
    ;n = ntimes-3
    if (False) then
    fld_name = "rp"
    rp = stitch_xy(stem,fld_name,n,kpt)

    fld_name = "np"
    np = stitch_xy(stem,fld_name,n,kpt)
    np = mask(np,np.gt.1,True)

    Dm = ( 6. * rp/ (1000.*pi) * np )^(1./3.)     
    mu = 10*(1.+tanh(1200.*(Dm-1.5e-3)))
    Dp = (Dm^3/((mu+3.)*(mu+2.)*(mu+1.)))^(1./3.) 

    vbar = 9.65 - 9.65*exp(600.*25.e-6)*(1.+600*Dp)^(-(4.+mu))
    flux = vbar*rp*fi->dn0(kpt)*8640.  ;cm/day

    plots(0) = gsn_csm_contour(wks,flux,cnres)
    
    zz = 1500.
    kpt = closest_val(zz,fi->zt)
    rl = stitch_xy(stem,"rl",n,kpt)
    plots(1) = gsn_csm_contour(wks,rl,cnres)

    jindex = round(floor(1.*maxind(ndtooned(flux))/256.),3)
    end if
    jindex = 80
    iproc  = jindex/(256/nprocs)
    if (iproc.lt.100) then
      if (iproc.lt.10) then
        fi = addfile(stem+".000"+iproc+".nc","r")
      else
        fi = addfile(stem+".00"+iproc+".nc","r")
      end if
    else
      fi = addfile(stem+".0"+iproc+".nc","r")
    end if
    jpt = jindex - iproc*(256/nprocs)
  
    z  = fspan(0.,3500.,350)
    x  = fspan(fi->xt(2),fi->xt(257),256)

    if (True) then
    cnres@vpWidthF  = 0.6
    cnres@vpHeightF = cnres@vpWidthF * z(dimsizes(z)-1)/(x(dimsizes(x)-1)-x(0))
    cnres@cnLevelSelectionMode = "ManualLevels"

    fld = new_shape(x,z,jpt,n,fi,"rl")
    fld = mask(fld,fld.gt.1.e-5,True)*1000.
    print (" n "+n+", "+max(fld))

    cnres@cnMaxLevelValF =   2.5
    cnres@cnMinLevelValF =   0.0
    cnres@cnLevelSpacingF =  0.1
    plots(2) = gsn_csm_contour(wks,fld,cnres)
    delete(fld)

    fld = new_shape(x,z,jpt,n,fi,"rp")
    fld = mask(fld,fld.gt.1.e-9,True)*1000.

    cnres@cnMaxLevelValF =   2.5
    cnres@cnMinLevelValF =   0.0
    cnres@cnLevelSpacingF =  0.1
    plots(3) = gsn_csm_contour(wks,fld,cnres)
    delete(fld)

    fld = new_shape(x,z,jpt,n,fi,"w")
    cnres@cnLevelSelectionMode = "ManualLevels"
    cnres@cnMaxLevelValF =   6.
    cnres@cnMinLevelValF =  -6.
    cnres@cnLevelSpacingF =  0.25

    plots(4) = gsn_csm_contour(wks,fld,cnres)
    delete (fld)
  end if

  panel_res = True
  panel_res@gsnFrame = False
  panel_res@gsnMaximize = True
  panel_res@gsnPanelTop = 0.95
  panel_res@gsnPanelXWhiteSpacePercent = 4.
; gsn_panel(wks,plots(0:1),(/1,2/),panel_res)
; frame(wks)

  gsn_panel(wks,plots(2:4),(/3,1/),panel_res)
  frame(wks)
  end do

end

