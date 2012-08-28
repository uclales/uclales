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
#   TODO: add units, description, etc. to NetCDF output                        #  
#                                                                              #
################################################################################

import numpy;
from netCDF4 import Dataset
from pylab import *
import time
import sys

# USER INPUT:
pref   = 'arm'       # Path to particle output 
nfx    = 8           # Number of cores in x-direction
nfy    = 3           # Number of cores in y-direction

# No need to change below
# ----------------------------------------------
vars = raw_input('variables (e.g. x,y,u,t,.. or all): ').split(',')   # commandline input
if(vars[0] == 'all'):
  vars = (['x','y','z','u','v','w','t','tv','rt','rl'])
print 'Merging: ', vars

globalstart = time.time()
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

x  = ncout.createVariable('x', 'f4',('time','particles',)) if 'x'  in vars else 1
y  = ncout.createVariable('y', 'f4',('time','particles',)) if 'y'  in vars else 1
z  = ncout.createVariable('z', 'f4',('time','particles',)) if 'z'  in vars else 1
u  = ncout.createVariable('u', 'f4',('time','particles',)) if 'u'  in vars else 1
v  = ncout.createVariable('v', 'f4',('time','particles',)) if 'v'  in vars else 1
w  = ncout.createVariable('w', 'f4',('time','particles',)) if 'w'  in vars else 1
t  = ncout.createVariable('t', 'f4',('time','particles',)) if 't'  in vars else 1
tv = ncout.createVariable('tv','f4',('time','particles',)) if 'tv' in vars else 1
rt = ncout.createVariable('rt','f4',('time','particles',)) if 'rt' in vars else 1
rl = ncout.createVariable('rl','f4',('time','particles',)) if 'rl' in vars else 1

# 3. Read in file-by-file and write to new NetCDF
loc = 0
for i in range(nfx):
  for j in range(nfy):
    start = time.time()
    filein = pref + '.particles.' + str(i).zfill(4) + str(j).zfill(4) + '.nc'
    ncin   = Dataset(filein, 'r')
    npl    = size(ncin.variables["particles"])

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
    print 'Processed: ',filein, ' in ', time.time()-start, 's' 

ncout.close()
print 'Total user time: ', time.time()-globalstart, 's' 

