#! /bin/csh

ncl plotfld.ncl 'fname=(/"dcbl.ts.nc"/)' \
'ffld=(/"zi1_bar","wmax","vtke"/)' 'foname="t1"' ncols=1 npages=1

ncl plotfld.ncl 'fname=(/"dcbl.ps.nc"/)' \
'ffld=(/"t","u","v","tot_tw","boy_prd","w_2","w_3","u_2","v_2"/)' 'foname="p1"' ncols=3 npages=1 tbeg=14400. tend=18000. 

 
