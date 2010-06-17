#! /bin/csh

#ncl plotfld.ncl 'fname=(/"rf01g.ts.nc","rf01k.ts.nc"/)' \
#'ffld=(/"zi1_bar","zb","wmax","lmax","vtke","lwp_bar","lhf_bar","cfrac"/)' 'foname="t1"' ncols=2 npages=1

#ncl plotfld.ncl 'fname=(/"bubble.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_dt.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_dt"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_s.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_s"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_g.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_g"' ncols=3 npages=1
ncl plotfld.ncl 'fname=(/"bubble_ra.ps.nc"/)' \
'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_ra"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_0000.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_0000"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_0060.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_0060"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_0600.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_0600"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_1800.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_1800"' ncols=3 npages=1
#ncl plotfld.ncl 'fname=(/"bubble_3600.ps.nc"/)' \
#'ffld=(/"ice","snow","graupel","l","rr"/)' 'foname="bubble_3600"' ncols=3 npages=1
