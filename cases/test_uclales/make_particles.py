# Create uniform initial particle distribution

# Input:
# ---------------
tstart = 0
nx     = 64
ny     = 64
nz     = 40
dx     = 25
dy     = 25
dz     = 10
# ---------------

# No need to change below.....
npart  = nx * ny * nz
x = dx / 2.
y = dy / 2.
z = dz / 2.

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
  y = dy / 2.
  x += dx 

partfile.close()
