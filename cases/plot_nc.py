import numpy
from pylab import *
from netCDF4 import Dataset

close('all')

def find_nearest(array,value):
  idx = (abs(array-value)).argmin()
  return idx

ncfile = Dataset('test_uclales/dcbl64x64.particles.nc', 'r')
ncfile2 = Dataset('test_uclales/dcbl64x64.particlestat.nc', 'r')

if(False):
  # Open statistics at once..
  zts = ncfile2.variables["zt"][:]
  nps = ncfile2.variables["np"][:,:]
  ts  = ncfile2.variables["time"][:]
  up  = ncfile2.variables["u"][:,:]
  vp  = ncfile2.variables["v"][:,:]
  wp  = ncfile2.variables["w"][:,:]
  nts = size(ts)
  
  figure()
  subplot(121)
  for i in range(0,nts,10):
    plot(nps[i,:],zts,label=str(ts[i]))
  legend()

  subplot(122)
  for i in range(0,nts,10):
    plot(wp[i,:],zts,label=str(ts[i]))
  legend()

if(True):
  #####################
  for time in range(0,500,2):
  #####################
    # Open particle NetCDF file
    t      = ncfile.variables["time"][:]
    x      = ncfile.variables["x"][time,:]
    y      = ncfile.variables["y"][time,:]
    z      = ncfile.variables["z"][time,:]
    u      = ncfile.variables["u"][time,:]
    v      = ncfile.variables["v"][time,:]
    w      = ncfile.variables["w"][time,:]
    np     = size(x[:])
 
    print numpy.average(w)

    figure()
    scatter(x,z,marker='o',s=1)
    ylim(0,100)
    xlim(0,800)
    grid()
  
    name = 'figs/' + str(time).zfill(5) + '.png'
    savefig(name)
    close('all')


