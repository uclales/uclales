# Create uniform initial particle distribution

# Input:
# ---------------
tstart = 0
nx     = 32
ny     = 32
nz     = 1
dx     = 25
dy     = 25
dz     = 20
# ---------------

# No need to change below.....
npart  = nx * ny * nz
x = (dx / 2.) #+ 0.01
y = (dy / 2.) #+ 0.01
z = (dz / 2.) #+ 0.01

partfile = open('partstartpos','w')
partfile.write(str(npart)+'\n')
for i in range(nx):
  print 'x = %3i / %3i' % (i+1,nx)
  for j in range(ny):
    for k in range(nz):
      partfile.write('%i  %.3e  %.3e  %.3e \n'%(tstart,x,y,z))
      z += dz
    z = dz / 2.  
    y += dy
  z = dz / 2.
  y = (dy / 2.) #+ 0.01
  x += dx 

partfile.close()
