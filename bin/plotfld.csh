#! /bin/csh

#ncl plotfld.ncl 'fname=(/"rf01g.ts.nc","rf01k.ts.nc"/)' \
#'ffld=(/"zi1_bar","zb","wmax","lmax","vtke","lwp_bar","lhf_bar","cfrac"/)' 'foname="t1"' ncols=2 npages=1

ncl plotfld.ncl 'fname=(/"bubble.ps.nc"/)' \
'ffld=(/"t","l","rr","ice","snow","graupel"/)' 'foname="p1"' ncols=3 npages=1 tbeg=0. tend=600. ymax=18000. ymin=00.

