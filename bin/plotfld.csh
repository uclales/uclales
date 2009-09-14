#! /bin/csh

ncl plotfld.ncl 'fname=(/"rf01g.ts.nc","rf01k.ts.nc"/)' \
'ffld=(/"zi1_bar","zb","wmax","lmax","vtke","lwp_bar","lhf_bar","cfrac"/)' 'foname="t1"' ncols=2 npages=1

ncl plotfld.ncl 'fname=(/"rf01g.ps.nc","rf01k.ps.nc"/)' \
'ffld=(/"t","u","v","q","l","tot_tw","tot_qw","boy_prd","rflx","w_2","w_3","u_2","v_2"/)' 'foname="p1"' ncols=3 npages=1 tbeg=10800. tend=14400. ymax=1200. ymin=00.

