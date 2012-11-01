################################################################################
#                                                                              # 
#   Create density weighted uniform particle distribution for the Lagrangian   #
#   Particle Tracking Module in UCLA-LES. Since UCLA-LES uses the anelastic    #
#   equations, we need to start with a density weighted distribution.          # 
#   Requires only the math module,which is available by default in Python.     #
#   Call: "python make_particles.py", all options are hardcoded.               # 
#   Pylab is optional for plotting, numpy for easy math calculations           #
#   requires module python/2.7-ve0 (or similar) on MPI systems                 #
#    ~Bart vS, July 2012                                                       #
#                                                                              #
################################################################################
from math import floor
#import numpy
#from pylab import *

# USER INPUT:
tstart      = 3600.      # 'release' time of particles
nxy         = 64**2.     # Goal number of particles per level
nz          = 68         # Number of vertical levels
dz          = 25.        # Vertical grid spacing
xysize      = 1600.      # Domain size LES

ps          = 101540.     # From UCLALES Namelist
th00        = 299.8       #     "          "
p00         = 1.0e5       #     "          "

# No need to change below
# ----------------------------------------------
# CONSTANTS:
cp          = 1005.
R           = 287.04
g           = 9.8
p00i        = 1. / p00
rcp         = R/cp
cpr         = cp/R

# Calculate number of particles.
# Required before writing the output, first line should
# contain the total number of particles. 
npart       = 0
pi00        = cp * (ps * p00i)**rcp + g * (0.5 * dz) / th00
for k in range(nz):
  dz0       = -float(k+1)  * dz
  pi0       = pi00 + g * dz0 / th00
  dn0       = ((cp**(1.-cpr)) * p00) / (R * th00 * pi0**(1.-cpr))
  nxyl      = round((nxy * dn0)**0.5)
  npart    += nxyl**2. 

# Write to partstartpos
partfile    = open('partstartpos','w')
partfile.write(str(int(npart))+'\n')
z = dz / 2. 
for k in range(nz):
  dz0       = -float(k+1)  * dz
  pi0       = pi00 + g * dz0 / th00
  dn0       = ((cp**(1.-cpr)) * p00) / (R * th00 * pi0**(1.-cpr))
  nxyl      = round((nxy * dn0)**0.5)
  dxy       = xysize / float(nxyl)
  x         = dxy / 2.
  y         = dxy / 2.
  print 'k = %3i / %3i, z = %.2f, # = %.3f' % (k+1,nz,z,nxyl)
  #print '%3i  %.2f  %.3f' % (k+1,z,nxyl)
  for i in range(int(nxyl)):
    for j in range(int(nxyl)):
      #partfile.write('%i  %.3e  %.3e  %.3e \n'%(tstart,x,y,z))
      partfile.write('%i  %.5f  %.5f  %.5f \n'%(tstart,x,y,z))
      y += dxy
    y = dxy/2.
    x += dxy
  z += dz
  
partfile.close()
