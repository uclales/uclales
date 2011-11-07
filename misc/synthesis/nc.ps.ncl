begin

;-----------
; input file
;-----------
; definition
  diri = "../bin/"
  fili = "g9.dws.ps.nc"
  fi = addfile(diri+fili,"r")

;------------------------
; variables in input file
;------------------------
; names
  vnm = (/"time","zt","zm","dn0","u","v","thetal","rt","rl","rr","cfrac", \
         "ndrop_cl","ss","u_var","v_var","w_var","w_3","theta_2",         \
         "rt_var","rl_var","rad_flx","precip","tot_tw","sfs_tw","tot_rw", \
         "sfs_rw","tot_uw","sfs_uw","tot_vw","sfs_vw","sgs_pr","sfs_tke", \
         "boy_prd","sfs_boy","shr_prd","trans","diss"/)

;------------
; output file
;------------
; definition
  filo = "rain-profiles.cdf"
  system("rm "+filo)
  fo = addfile(filo,"c")
; global attributes
  GlbAtt = True
  GlbAtt@PROCID = 0
  GlbAtt@NPROCS = 16
;  GlbAtt@NPTS = 16384
  GlbAtt@FillValue = -999
  GlbAtt@Source = "UCLA-LES (Verica Savic-Jovcic and Bjorn Stevens)"
  GlbAtt@history = "Created on "+systemfunc("date +%D")
  GlbAtt@title = "Drizzling RF02 case with cloud-droplet sedimentation" + \
  " - bulk microphysics based on Seifert and Beheng (2001)"
  fileattdef(fo,GlbAtt)

;--------------------------
; dimensions in output file
;--------------------------
; names
  vdn = (/"time","zt","zw"/)
  vds = (/dimsizes(fi->time),dimsizes(fi->zt)-1,dimsizes(fi->zm)/)
  vdu = (/True,False,False/)
  filedimdef(fo,vdn,vds,vdu)

;-------------------------
; variables in output file
;-------------------------
; names
  vn = (/"time","zt","zw","rho","u","v","thetal","qt","ql","qr","cfrac", \
  "ndrop_cld","ss","u_var","v_var","w_var","w_skw","thetal_var", \
  "qt_var","ql_var","rad_flx","precip","tot_tw","sgs_tw","tot_qw", \
  "sgs_qw","tot_uw","sgs_uw","tot_vw","sgs_vw","sgs_pr","sgs_tke", \
  "tot_boy","sgs_boy","tot_shr","transport","dissipation"/)
; longnames
  vln = (/"Time", \
  "Altitude of layer mid-points", \
  "Altitude of layer boundaries", \
  "Reference density", \
  "Zonal wind", \
  "Meridional wind", \
  "Liquid water potential temperature", \
  "Total water mixing ratio", \
  "Condensed water mixing ratio", \
  "Rain water mixing ratio",\
  "Fraction of grid cells with droplet concetration > 20 cm^-3", \
  "Mean cloud droplet concentration in grid cells with droplet concentration > 20 cm^-3", \
  "Mean supersaturation", \
  "Resolved variance of zonal wind", \
  "Resolved variance of meridional wind",  \
  "Resolved variance of vertical wind", \
  "Resolved <w'^3>",  \
  "Resolved variance of thetal", \
  "Resolved variance of qt", \
  "Resolved variance of ql", \
  "Net radiative flux", \
  "Precipitation flux (positive downward)", \
  "Total thetal flux, including subgrid-scale and precipitation", \
  "Subgrid-scale thetal flux", \
  "Total qt flux, including subgrid-scale and precipitation", \
  "Subgrid-scale qt flux", \
  "Total (sgs plus resolved) zonal momentum flux", \
  "Subgrid-scale zonal momentum flux", \
  "Total (sgs plus resolved) meridional momentum flux", \
  "Subgrid-scale meridional momentum flux", \
  "Subgrid-scale Prandtl number flux", \
  "Subgrid-scale TKE", \
  "Total (sgs plus resolved) buoyancy TKE production", \
  "Subgrid-scale buoyancy TKE production", \
  "Total (sgs plus resolved) shear TKE production", \
  "Resolved TKE transport (turbulent plus pressure)", \
  "TKE dissipation (explicit plus numerical)"/)
; units
  vu = (/"s","m","m","kg/m^3","m/s","m/s","K","g/kg","g/kg","g/kg","",  \
  "cm^-3","%","m^2/s^2","m^2/s^2","m^2/s^2","m^3/s^3","K^2",     \
  "g^2/kg^2","g^2/kg^2","W/m^2","W/m^2","W/m^2","W/m^2","W/m^2", \
  "W/m^2","kg/ms^2","kg/ms^2","kg/ms^2","kg/ms^2","","m^2/s^2",  \
  "m^2/s^3","m^2/s^3","m^2/s^3","m^2/s^3","m^2/s^3"/)
  
;-------------------------
; writing in the output file
;-------------------------
; dimensions

  do i = 0, 3
    dmy = (/fi->$vnm(i)$/)
    dmy@units = vu(i)
    dmy@long_name = vln(i)
    if (vnm(i).eq."time") then
      dmy!0 = "time"
      fo->time = dmy
    end if
    if (vnm(i).eq."zt") then
      dmy!0 = "zt"
      fo->zt = dmy(1:)
    end if    
    if (vnm(i).eq."zm") then
      dmy!0 = "zw"
      fo->zw = dmy
    end if
    if (vnm(i).eq."dn0") then
      dmyrho = new((/dimsizes(fi->time),dimsizes(fi->zt)-1/),float)
      dmyrho@units = dmy@units
      dmyrho@long_name = dmy@long_name
      delete(dmyrho@_FillValue)
      do j = 0, dimsizes(fi->time)-1
        dmyrho(j,:) = (/dmy(1:)/)
      end do
      dmyrho!0 = "time"
      dmyrho!1 = "zt"
      fo->rho = dmyrho
    end if
    delete(dmy)
  end do
;
; variables
;
  do i = 4, dimsizes(vn)-1
    dmy = fi->$vnm(i)$
    dmydim = getvardims(dmy)
    if (dmydim(1).eq."zt") then
      dmyZT = (/dmy(:,1:)/)
      dmyZT@units = vu(i)
      dmyZT@long_name = vln(i)
      dmyZT!0 = "time"
      dmyZT!1 = "zt"
      fo->$vn(i)$ = dmyZT
      delete(dmyZT)
    end if
    if (dmydim(1).eq."zm") then
      dmyZW = (/dmy/)
      dmyZW@units = vu(i)
      dmyZW@long_name = vln(i)
      dmyZW!0 = "time"
      dmyZW!1 = "zw"
      fo->$vn(i)$ = dmyZW
      delete(dmyZW)
    end if
    delete(dmy)
  end do
  
end
