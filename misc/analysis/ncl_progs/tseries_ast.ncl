;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"   

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;load "~/Library/NCL/myfuncs.ncl"   
;

begin

  wks = gsn_open_wks("ps","tseries_trDct")
 stems = (/"trDct"/)
; stems=(/"astx40"/)
 ; stems = (/"rf02","rf02kh","rf02nm","rf02R","rf02khR"/)
  tplts = new ( (/9/), "graphic")

  xyres = True
  xyres@gsnDraw = False
  xyres@gsnFrame = False
  xyres@tmXTOn = False
  xyres@tmXTBorderOn = False
  xyres@tmYROn = False
  xyres@tmYRBorderOn = False

  colors = (/"black","red","dodgerblue","green","red"/)
  dashes = (/0,0,0,0,1/)
 
  t1 = 3600.
  t2 = 7200.

  do n=0,dimsizes(stems)-1
    xyres@xyLineColor = colors(n)
    xyres@xyDashPattern = dashes(n)
    tsfi = addfile(stems(n)+".ts.nc","r")
    psfi = addfile(stems(n)+".ps.nc","r")

    time = tsfi->time/86400.
    time@long_name = "time (days)"
;    wstar = avg(tsfi->wstar({t1:t2}))^(1./3.)
;    wstar = avg(tsfi->zi1_bar({t1:t2})*100.*9.8/300./1200.)^(1./3.)
;    print (wstar)

;zcmn=tsfi->zcmn
;zbmn=tsfi->zbmn

;print(zcmn)
;print(zbmn)

;  if (any(isnan_ieee(zcmn))) then
;      value = -999.
;      replace_ieeenan (zcmn, value, 0)
;      zcmn@_FillValue = value
;  end if
;  if (any(isnan_ieee(zbmn))) then
;      value = -999.
;      replace_ieeenan (zbmn, value, 0)
;      zbmn@_FillValue = value
;  end if
; printMinMax(zcmn,True)
; printMinMax(zbmn,True)


    if (n .eq. 0) then
    xyres@trYMinF= 0.
    xyres@trYMaxF= 2500.
    xyres@tiYAxisString="inversion height (m)"
      tplts(0) = gsn_csm_xy(wks,time,tsfi->zi1_bar,xyres)

    xyres@tiYAxisString= "CT/CB (m)"
      tplts(1) = gsn_csm_xy(wks,time,tsfi->zc,xyres)
      ovrly = gsn_csm_xy(wks,time,tsfi->zb,xyres)
      overlay(tplts(1),ovrly)
;    xyres@xyLineColor = "yellow"
;      ovrly2 = gsn_csm_xy(wks,time,tsfi->zcmn,xyres)
;      overlay(tplts(1),ovrly2)
;      ovrly3 = gsn_csm_xy(wks,time,tsfi->zbmn,xyres)
;      overlay(tplts(1),ovrly3)
    xyres@xyLineColor = "black"
    xyres@trYMinF= 0.
    xyres@trYMaxF= 300.
    xyres@tiYAxisString= "lwp (g/m2)"
      tplts(2) = gsn_csm_xy(wks,time,tsfi->lwp_bar,xyres)
    xyres@trYMinF= 0.
    xyres@trYMaxF= 1.
   xyres@tiYAxisString="CF"
    tplts(3)= gsn_csm_xy(wks,time,tsfi->cfrac,xyres)
    xyres@trYMinF= 0.
    xyres@trYMaxF= 10.
    xyres@tiYAxisString= "rwp (g/m2)"
      tplts(4) = gsn_csm_xy(wks,time,tsfi->rwp_bar,xyres)
    xyres@trYMinF= 0.
    xyres@trYMaxF= 15.
    xyres@tiYAxisString="pp surface (W/m2)"
      tplts(5) = gsn_csm_xy(wks,time,(tsfi->prcp)*2.5e+06*1.223,xyres)
   printMinMax(tsfi->prcp,True)
    xyres@trYMinF= 294.
    xyres@trYMaxF= 306.
    xyres@tiYAxisString="T srf (K)"
      tplts(6) = gsn_csm_xy(wks,time,tsfi->tsrf,xyres)
    xyres@trYMinF= 0.
    xyres@trYMaxF= 20.
    xyres@tiYAxisString="SHF (W/m2)"
      tplts(7) = gsn_csm_xy(wks,time,tsfi->shf_bar,xyres)
    xyres@trYMinF= 50.
    xyres@trYMaxF= 200.
    xyres@tiYAxisString="LHF (W/m2)"
      tplts(8) = gsn_csm_xy(wks,time,tsfi->lhf_bar,xyres)

;print(tsfi->tsrf)
;printMinMax(tsfi->lhf_bar,True)
;printMinMax(tsfi->shf_bar,True)
    else
      ovrly = gsn_csm_xy(wks,time,tsfi->zi1_bar,xyres)
      overlay(tplts(0),ovrly)
      ovrly = gsn_csm_xy(wks,time,tsfi->zc,xyres)
      overlay(tplts(1),ovrly)
      ovrly1 = gsn_csm_xy(wks,time,tsfi->zb,xyres)
      overlay(tplts(1),ovrly1)
;if (n.eq. 1) then
;    xyres@xyLineColor = "yellow"
;      ovrly2 = gsn_csm_xy(wks,time,tsfi->zcmn,xyres)
;      overlay(tplts(1),ovrly2)
;      ovrly3 = gsn_csm_xy(wks,time,tsfi->zbmn,xyres)
;      overlay(tplts(1),ovrly3)
;    xyres@xyLineColor = "red"
;end if
      ovrly = gsn_csm_xy(wks,time,tsfi->lwp_bar,xyres)
      overlay(tplts(2),ovrly)
      ovrly = gsn_csm_xy(wks,time,tsfi->cfrac,xyres)
      overlay(tplts(3),ovrly)
      ovrly = gsn_csm_xy(wks,time,tsfi->rwp_bar,xyres)
      overlay(tplts(4),ovrly)
      ovrly = gsn_csm_xy(wks,time,(tsfi->prcp)*2.5e+06*1.223,xyres)
      overlay(tplts(5),ovrly)
      ovrly = gsn_csm_xy(wks,time,tsfi->tsrf,xyres)
      overlay(tplts(6),ovrly)
      ovrly = gsn_csm_xy(wks,time,tsfi->shf_bar,xyres)
      overlay(tplts(7),ovrly)
      ovrly = gsn_csm_xy(wks,time,tsfi->lhf_bar,xyres)
      overlay(tplts(8),ovrly)

;printMinMax(tsfi->prcp,True)
;printMinMax(tsfi->lhf_bar,True)
;printMinMax(tsfi->shf_bar,True)
    end if
    delete (time)
  end do
   
  pres = True
  pres@gsnPanelRowSpec = True
  pres@gsnFrame = False
  pres@gsnMaximize = True
  shape = (/3,3,3/)
  gsn_panel(wks,tplts(0:8),shape,pres)
  frame(wks)

end






