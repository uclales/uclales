# Create uniform initial particle distribution

# Input:
# ---------------
tstart = 100
nx     = 32
ny     = 32
nz     = 98
dx     = 25
dy     = 25
dz     = 25
# ---------------

# No need to change below.....
npart  = nx * ny * nz
x = (dx / 2.) #+ 0.01
y = (dy / 2.) #+ 0.01
z = (dz / 2.) #+ 0.01

np = 0

partfile = open('partstartpos','w')
partfile.write(str(npart)+'\n')
# Changed: write z-dimension first, easier
# for grouping particles by initial height
# and id.
for k in range(nz):
  print 'k = %3i / %3i' % (k+1,nz)
  for i in range(nx):
    for j in range(ny):
      partfile.write('%i  %.3e  %.3e  %.3e \n'%(tstart,x,y,z))
      np += 1
      y += dy
    y = dy/2.
    x += dx
  x = dx/2.
  y = dy/2.
  z += dz

print np

partfile.close()
