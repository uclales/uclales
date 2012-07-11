import numpy
from pylab import *
# Create uniform initial particle distribution

  # ----------------------------------------------
  # NEW: number of particles on each level based
  #      on base state density
  # ----------------------------------------------
if(True):
  # INPUT:
  tstart      = 300.
  nxy         = 32**2.       # Goal # per level
  nz          = 118          # Number of vertical levels
  dz          = 25.          # Vertical grid spacing
  xysize      = 800.        # Domain size LES
  
  ps          = 101540.      # From UCLALES Namelist
  th00        = 299.8        #     "          "
  p00         = 1.0e5        #     "          "
  
  # ----------------------------------------------
  # CONSTANTS:
  cp          = 1005.
  R           = 287.04
  g           = 9.8
  p00i        = 1. / p00
  rcp         = R/cp
  cpr         = cp/R
  
  pi0         = numpy.zeros(nz)
  dn0         = numpy.zeros(nz)
  nxyl        = numpy.zeros(nz)
  z           = numpy.zeros(nz)
  
  npart       = 0
 
  # Calculate base state density & number of particles 
  pi00        = cp * (ps * p00i)**rcp + g * (0.5 * dz) / th00
  for k in range(nz):
    z[k]      = float(k+0.5) * dz
    dz0       = -float(k+1)  * dz
    pi0[k]    = pi00 + g * dz0 / th00
    dn0[k]    = ((cp**(1.-cpr)) * p00) / (R * th00 * pi0[k]**(1.-cpr))
    nxyl[k]   = floor((nxy * dn0[k])**0.5)
    npart    += nxyl[k]**2. 

  partfile = open('partstartpos','w')
  partfile.write(str(int(npart))+'\n')
  # Write z-dimension first, easier
  # for grouping particles by initial height
  # and id.
  for k in range(nz):
    dxy    = xysize / float(nxyl[k])
    x      = dxy / 2.
    y      = dxy / 2. 
    print 'k = %3i / %3i, # = %.4i^2' % (k+1,nz,nxyl[k])
    for i in range(int(nxyl[k])):
      for j in range(int(nxyl[k])):
        partfile.write('%i  %.3e  %.3e  %.3e \n'%(tstart,x,y,z[k]))
        y += dxy
      y = dxy/2.
      x += dxy
    x = dxy/2.
    y = dxy/2.
  partfile.close()


# ----------------------------------------------
# OLD: Uniform distribution
# ----------------------------------------------
if(False):
  # Input:
  # ---------------
  tstart = 100
  nx     = 200
  ny     = 200
  nz     = 200
  dx     = 25
  dy     = 25
  dz     = 25
  # ---------------
  
  # No need to change below.....
  npart  = nx * ny * nz
  x = (dx / 2.) 
  y = (dy / 2.) 
  z = (dz / 2.) 
  
  partfile = open('partstartpos','w')
  partfile.write(str(npart)+'\n')
  # Write z-dimension first, easier
  # for grouping particles by initial height
  # and id.
  for k in range(nz):
    print 'k = %3i / %3i' % (k+1,nz)
    for i in range(nx):
      for j in range(ny):
        partfile.write('%i  %.3e  %.3e  %.3e \n'%(tstart,x,y,z))
        y += dy
      y = dy/2.
      x += dx
    x = dx/2.
    y = dy/2.
    z += dz
  partfile.close()
