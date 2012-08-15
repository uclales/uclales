################################################################################
#                                                                              # 
#   Merge raw particle output from the Lagrangian particle module in UCLA-LES  #
#   Requires numpy, netCDF4, pylab, which are all provided on MPI systems      #
#   and certain DKRZ systems (at least Lizard) via module python/2.7-ve0       #
#    - Bart van Stratum, Aug 2012                                              #
#                                                                              #
#   Command-line input: script asks for the variables to be merged             #   
#   'all' puts all variables (x,y,z,u,v,w,t,tv,rt,rl) into one NetCDF file     #
#                                                                              #
#   TODO: fix long names                                                       #  
#                                                                              #
################################################################################

import numpy;
from netCDF4 import Dataset
from pylab import *
import time as time2
import sys

# USER INPUT:
pref   = 'arm'       # Path to particle output 
nfx    = 1           # Number of cores in x-direction
nfy    = 1           # Number of cores in y-direction
compr  = False       # Use zlib compression (saves ~20%, but script very slow)?

# No need to change below
# ----------------------------------------------
vars = raw_input('variables (e.g. x,y,u,t,.. or all): ').split(',')   # commandline input
if(vars[0] == 'all'):
  vars = (['x','y','z','u','v','w','t','tv','rt','rl'])
print 'Merging: ', vars

globalstart = time2.time()
# 1. Find total number of particles: 
np = 0
for i in range(nfx):
  for j in range(nfy):
    filein = pref + '.particles.' + str(i).zfill(4) + str(j).zfill(4) + '.nc'
    ncin   = Dataset(filein, 'r')
    np    += size(ncin.variables["particles"])

    if(i==0 and j==0):  # check if variables present
      for var in vars:
        try:
          varin = ncin.variables[var.strip()]
        except KeyError as e:
          sys.exit("One or more variables not found, stopping..")
print 'Found ', np, ' particles!'

# 2. Create new NetCDF
fileout = pref + '.particles.' + '_'.join(vars) + '.nc'
ncout   = Dataset(fileout, 'w', format='NETCDF4')
pdim    = ncout.createDimension('particles', np)
tdim    = ncout.createDimension('time'      ,None)
print 'output file = ',fileout

time            = ncout.createVariable('time', 'f4',('time',),zlib=compr)
time.units      = 's' 
time.longname   = 'Time' 
if 'x' in vars:
  x             = ncout.createVariable('x', 'f4',('time','particles',),zlib=compr)
  x.units       = ncin.variables["x"].units
  x.longname    = ncin.variables["x"].longname.rstrip(' \t\r\n\0') 
if 'y' in vars: 
  y             = ncout.createVariable('y', 'f4',('time','particles',),zlib=compr)
  y.units       = ncin.variables["y"].units
  y.longname    = ncin.variables["y"].longname 
if 'z' in vars: 
  z             = ncout.createVariable('z', 'f4',('time','particles',),zlib=compr)
  z.units       = ncin.variables["z"].units
  z.longname    = ncin.variables["z"].longname 
if 'u' in vars: 
  u             = ncout.createVariable('u', 'f4',('time','particles',),zlib=compr)
  u.units       = ncin.variables["u"].units
  u.longname    = ncin.variables["u"].longname 
if 'v' in vars: 
  v             = ncout.createVariable('v', 'f4',('time','particles',),zlib=compr)
  v.units       = ncin.variables["v"].units
  v.longname    = ncin.variables["v"].longname 
if 'w' in vars: 
  w             = ncout.createVariable('w', 'f4',('time','particles',),zlib=compr)
  w.units       = ncin.variables["w"].units
  w.longname    = ncin.variables["w"].longname 
if 't' in vars: 
  t             = ncout.createVariable('t', 'f4',('time','particles',),zlib=compr)
  t.units       = ncin.variables["t"].units
  t.longname    = ncin.variables["t"].longname 
if 'tv' in vars:
  tv            = ncout.createVariable('tv', 'f4',('time','particles',),zlib=compr)
  tv.units      = ncin.variables["tv"].units
  tv.longname   = ncin.variables["tv"].longname 
if 'rt' in vars:
  rt            = ncout.createVariable('rt', 'f4',('time','particles',),zlib=compr)
  rt.units      = ncin.variables["rt"].units
  rt.longname   = ncin.variables["rt"].longname 
if 'rl' in vars:
  rl            = ncout.createVariable('rl', 'f4',('time','particles',),zlib=compr)
  rl.units      = ncin.variables["rl"].units
  rl.longname   = ncin.variables["rl"].longname 
ncin.close()

# 3. Read in file-by-file and write to new NetCDF
loc = 0
for i in range(nfx):
  for j in range(nfy):
    start = time2.time()
    filein = pref + '.particles.' + str(i).zfill(4) + str(j).zfill(4) + '.nc'
    ncin   = Dataset(filein, 'r')
    npl    = size(ncin.variables["particles"])

    if(i==0 and j==0):
      time[:] = ncin.variables["time"][:]

    if 'x' in vars:
      x[:,loc:loc+npl]  = ncin.variables["x"][:,:] 
    if 'y' in vars:
      y[:,loc:loc+npl]  = ncin.variables["y"][:,:] 
    if 'z' in vars:
      z[:,loc:loc+npl]  = ncin.variables["z"][:,:] 
    if 'u' in vars:
      u[:,loc:loc+npl]  = ncin.variables["u"][:,:] 
    if 'v' in vars:
      v[:,loc:loc+npl]  = ncin.variables["v"][:,:] 
    if 'w' in vars:
      w[:,loc:loc+npl]  = ncin.variables["w"][:,:] 
    if 't' in vars:
      t[:,loc:loc+npl]  = ncin.variables["t"][:,:] 
    if 'tv' in vars:
      tv[:,loc:loc+npl] = ncin.variables["tv"][:,:]
    if 'rt' in vars:
      rt[:,loc:loc+npl] = ncin.variables["rt"][:,:]
    if 'rl' in vars:
      rl[:,loc:loc+npl] = ncin.variables["rl"][:,:]

    ncin.close()
    loc += npl
    print 'Processed: ',filein, ' in ', time2.time()-start, 's' 

ncout.close()
print 'Total user time: ', time2.time()-globalstart, 's' 

