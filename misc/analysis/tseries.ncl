load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   
load "~/Library/NCL/myfuncs.ncl"   
;

begin

  wks = gsn_open_wks("x11","tseries")
  
  stems = (/"dcbl_132x2","dcbl_256"/)
  tplts = new ( (/6/), "graphic")

  xyres = True
  xyres@gsnDraw = False
  xyres@gsnFrame = False
  xyres@tmXTOn = False
  xyres@tmXTBorderOn = False
  xyres@tmYROn = False
  xyres@tmYRBorderOn = False

  colors = (/"black","red","dodgerblue"/)
  t1 = 16200.
  t2 = 18000.

  do n=0,dimsizes(stems)-1
    xyres@xyLineColor = colors(n)
    tsfi = addfile(stems(n)+".ts.nc","r")
    psfi = addfile(stems(n)+".ps.nc","r")

    time = sqrt(tsfi->time)
    time@long_name = "sqrt(t) [s:S:1/2:N:]"
    wstar = avg(tsfi->wstar({t1:t2}))^(1./3.)
    wstar = avg(tsfi->zi1_bar({t1:t2})*100.*9.8/300./1200.)^(1./3.)
    print (wstar)
    xyres@trYMinF= 0.
    xyres@trYMaxF= 1750.

    if (n .eq. 0) then
      tplts(0) = gsn_csm_xy(wks,time,tsfi->zi1_bar,xyres)

      yy = dim_avg(psfi->w_var(zm|:,{time|t1:t2}))
      yy = mask(yy, yy.gt.0., True) / wstar^2
      tplts(1) = gsn_csm_xy(wks,yy,psfi->zm,xyres)
      yy = dim_avg(psfi->w_3(zm|:,{time|t1:t2}))/yy^1.5 /wstar^3
      tplts(2) = gsn_csm_xy(wks,yy,psfi->zm,xyres)
      yy = dim_avg(psfi->theta(zt|:,{time|t1:t2})) 
      tplts(3) = gsn_csm_xy(wks,yy,psfi->zt,xyres)
      yy = dim_avg(psfi->theta_2(zt|:,{time|t1:t2})) 
      tplts(4) = gsn_csm_xy(wks,yy,psfi->zt,xyres)
      yy = dim_avg(psfi->tot_tw(zm|:,{time|t1:t2})) 
      tplts(5) = gsn_csm_xy(wks,yy,psfi->zm,xyres)
    else
      ovrly = gsn_csm_xy(wks,time,tsfi->zi1_bar,xyres)
      overlay(tplts(0),ovrly)

      yy = dim_avg(psfi->w_var(zm|:,{time|t1:t2}))
      yy = mask(yy, yy.gt.0., True)  /wstar^2
      ovrly = gsn_csm_xy(wks,yy,psfi->zm,xyres)
      overlay(tplts(1),ovrly)
      yy = dim_avg(psfi->w_3(zm|:,{time|t1:t2}))/yy^1.5 /wstar^3
      ovrly = gsn_csm_xy(wks,yy,psfi->zm,xyres)
      overlay(tplts(2),ovrly)
      yy = dim_avg(psfi->theta(zt|:,{time|t1:t2})) 
      ovrly = gsn_csm_xy(wks,yy,psfi->zt,xyres)
      overlay(tplts(3),ovrly)
      yy = dim_avg(psfi->theta_2(zt|:,{time|t1:t2})) 
      ovrly = gsn_csm_xy(wks,yy,psfi->zt,xyres)
      overlay(tplts(4),ovrly)
      yy = dim_avg(psfi->tot_tw(zm|:,{time|t1:t2}))
      ovrly = gsn_csm_xy(wks,yy,psfi->zm,xyres)
      overlay(tplts(5),ovrly)
    end if
    delete (time)
    delete (yy)
  end do
   
  pres = True
  pres@gsnPanelRowSpec = True
  pres@gsnFrame = False
  pres@gsnMaximize = True
  shape = (/3,3/)
  gsn_panel(wks,tplts(0:5),shape,pres)
  frame(wks)

end






