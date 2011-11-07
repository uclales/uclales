begin

;-----------
; input file
;-----------
; definition
  diri = "../bin/"
  fili = "g9.dws.ts.nc"
  fi = addfile(diri+fili,"r")
;-------------------------
; variables in input file
;-------------------------
; names
  vnm = (/"time","zi_bar","zi_var","zb_bar","zb_var","lwp_bar","lwp_var", \
         "cfrac","ndrop","tke","wstar","w2","shf_bar","lhf_bar","precip", \
         "precip_v","precip_m","ustar","wmax","cfl","mxdiv"/)
;  print(" "+vnm)
  vnb = dimsizes(vnm)

;------------
; output file
;------------
; definition
  filo = "rain-scalars.cdf"
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
  GlbAtt@title = "Drizzling RF02 case with cloud-droplet sedimentation - bulk microphysics based on Seifert and Beheng (2001)"
  fileattdef(fo,GlbAtt)

;--------------------------
; dimensions in output file
;--------------------------
; names
  vdn = "time"
; sizes
  vds = 1000
; unlimited sizes
  vdu = True
; definition of dimensions in the output file
  filedimdef(fo,vdn,vds,vdu)

;-------------------------
; variables in output file
;-------------------------
; names
  vn = (/"time","zi","zi_var","zb","zb_var","lwp","lwp_var","cfrac", \
         "ndrop_cld","tke","wstar","w2_max","shf","vhf","precip", \
         "precip_var","precip_max","ustar","wmax","cfl","divmax"/)
; longnames
  vln = (/"Time", \
          "Mean height of qt=8g/kg surface", \
          "Variance of height of qt=8g/kg surface", \
          "Mean height of bottom of lowermost grid cell with droplet concentration > 20 cm^-3", \
          "Variance of height of bottom of lowermost grid cell with droplet concentration > 20 cm^-3", \
          "Mean liquid water path", \
          "Liquid water path variance", \
          "Fraction of columns with liquid water path > 20 g/m^2", \
          "Mean cloud droplet concentration in grid cells with droplet concentration > 20 cm^-3", \
          "Vertically integrated rho*TKE (sgs plus resolved)", \
          "Deardorf convective velocity scale", \
          "Maximum value of layer-averaged vertical velocity variance",  \
          "Mean upward surface sensible heat flux", \
          "Mean upward surface vapor latent heat flux",  \
          "Mean (downward) surface precipitation flux", \
          "Variance of (downward) surface precipitation flux", \
          "Maximum (downward) surface precipitation flux", \
          "Mean value of friction velocity", \
          "Domain wide maximum of vertical velocity", \
          "CFL stability criterion", \
          "Maximum divergence in Domain" /)
; units
  vu = (/"s","m","m^2","m","m^2","g/m^2","g^2/m^4","","cm^-3","kg/s^2", \
         "m/s","m^2/s^2","W/m^2","W/m^2","W/m^2","W^2/m^4","W/m^2", \
         "m/s","m/s","","s^-1" /)

;-------------------------
; write in the output file
;-------------------------
; dimensions
  fo->time = fi->time
; variables
  do i = 1, dimsizes(vn)-1
    dmy = (/fi->$vnm(i)$/)
    dmy@units = vu(i)
    dmy@long_name = vln(i)
    dmy!0 = "time"
    fo->$vn(i)$ = dmy
  end do

end
