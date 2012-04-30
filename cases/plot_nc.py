import numpy
from pylab import *
from netCDF4 import Dataset

close('all')

def find_nearest(array,value):
  idx = (abs(array-value)).argmin()
  return idx

ncfile = Dataset('test_uclales/dcbl.particles.nc', 'r')
ncfile2 = Dataset('test_uclales/dcbl.particlestat.nc', 'r')

if(False):
  # Open statistics at once..
  zts = ncfile2.variables["zt"][:]
  nps = ncfile2.variables["np"][:,:]
  ts  = ncfile2.variables["time"][:]
  nts = size(ts)
  
  figure()
  for i in range(0,nts):
    plot(nps[i,:],zts,label=str(ts[i]))
  legend()

if(True):
  #####################
  for time in range(0,50):
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
    plot(w,'o',ms=1)
 
    figure()
    scatter(x,z,marker='o',s=1)
    ylim(2000,2120)
    ylim(0,100)
    xlim(0,800)
  
    name = 'figs/' + str(time).zfill(5) + '.png'
    savefig(name)
    close('all')


