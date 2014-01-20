"""
Create sound_in and ls_flux_in for BOMEX in UCLA-LES
Bart van Stratum, Jan 2014
"""

import numpy as np
from pylab import *
import sys

close('all')

# --------------------------------------------------------
# From Siebesma et al, 2003:
zi  = ([0,520,1480,2000,3000])    # scalars input grid
qi  = ([17.,16.3,10.7,4.2,3.0])
ti  = ([298.7,298.7,302.4,308.2,311.85])

zi2 = ([0,700,3000])              # momentum input grid
ui  = ([-8.75,-8.75,-4.61])
vi  = ([0,0,0])

zi3 = ([0,1500,2100,3000])        # subsidence grid
wi  = ([0,-0.65,0,0])

zi4 = ([0,1500,3000])             # heating rate grid
Qri = ([-2.,-2.,0])

zi5 = ([0,300,500,3000])          # moisture adv grid
dqi = ([-1.2,-1.2,0,0])

wt  = 8.e-3    # prescribed in time
wq  = 5.2e-5   # prescribed in time
ps  = 1015.

# --------------------------------------------------------
# Output grid:
zo  = np.arange(0,3000.01,25)
# IMPORTANT, input UCLA in W/m2 so conversion needed:
dns = 1.136    
tvar = ([0,432000])
# --------------------------------------------------------

def interp(zor,znew,s):
  """
  zor  = original grid
  znew = grid to interpolate to
  s    = field to interpolate
  """

  snew = np.zeros_like(znew)
  kb   = 0
  for k in range(len(znew)):
    if(znew[k] > zor[kb+1]):
      kb += 1
    dsdz = (s[kb+1]-s[kb])/(zor[kb+1]-zor[kb])
    snew[k] = s[kb] + dsdz * (znew[k]-zor[kb])
  return snew 

# interpolate fields
qo  = interp(zi,zo,qi)    # total water mix ratio
to  = interp(zi,zo,ti)    # liquid water pot temp
uo  = interp(zi2,zo,ui)   # u
vo  = interp(zi2,zo,vi)   # v
wo  = interp(zi3,zo,wi)   # subsidence vel.
qro = interp(zi4,zo,Qri)  # heating rate
dqo = interp(zi5,zo,dqi)  # moist. advec.
ug  = -10. + 1.8e-3 * zo  # geostrophic wind
vg  = np.zeros_like(ug)   # geostrophic wind

# Write output fields

# ------------------------------------------------
# Input sounding (sound_in): 
# ------------------------------------------------
ss = open('sound_in','w')
for k in range(len(zo)):
  zp = ps if (k==0) else zo[k]   # for UCLA: lowest height == surface pressure!
  ss.write('{0:10.4F} {1:1.4E} {2:1.4E} {3:+1.4E} {4:+1.4E}\n'.format(zp,to[k],qo[k],uo[k],vo[k]))
ss.close()

# ------------------------------------------------
# Large-scale forcings (ls_flux_in):
# ------------------------------------------------
ls = open('ls_flux_in','w')
ls.write('# BOMEX case\n')
ls.write('# See Siebesma et al, 2003\n')

ls.write('#{0:>9s} {1:>10s} {2:>10s} {3:>10s} {4:>10s} {5:>10s}  \n'.format('time [s]','H [W/m2]','LE [W/m2]','Ts [K]','qs [g/kg]','Ps [pa]'))
for t in tvar:
  ls.write('{0:10.2F} {1:10.5E} {2:10.5E} {3:6.1F} {4:6.1F} {5:10.2F}\n'.format(t,wt*dns*1005.,wq*dns*2.45e6,-9999.,-9999.,ps*100.))

for t in tvar:
  ls.write('\n')
  ls.write('#{0:>9s} {1:>10s} {2:>10s} {3:>10s} {4:>10s} {5:>10s}  \n'.format('zf [m]','Ug [m/s]','Vg [m/s]','ws [m/s]','dq/dt','dT/dt'))
  ls.write('# {0:<10.2F} \n'.format(t))
  for k in range(len(zo)):
    ls.write('{0:+10.2F} {1:+10.5F} {2:+10.5F} {3:+10.5F} {4:+10.5E} {5:+10.5E}\n'.format(zo[k],ug[k],vg[k],wo[k]/100.,dqo[k]*1.e-8,qro[k]/86400.))
ls.close()


if(True):
  figure()
  subplot(231)
  plot(ti,zi,'o',label='input')
  plot(to,zo,'-',label='interpolated')
  legend(frameon=False)
  
  subplot(232)
  plot(qi,zi,'o',label='input')
  plot(qo,zo,'-',label='interpolated')
  legend(frameon=False)
  
  subplot(233)
  plot(ui,zi2,'o',label='input')
  plot(uo,zo,'-',label='interpolated')
  plot(ug,zo,'--',label='ugeo')
  plot(vg,zo,'--',label='ugeo')
  plot(vi,zi2,'o',label='input')
  plot(vo,zo,'-',label='interpolated')
  legend(frameon=False)
  
  subplot(234)
  plot(wi,zi3,'o',label='input')
  plot(wo,zo,'-',label='interpolated')
  legend(frameon=False)
  
  subplot(235)
  plot(Qri,zi4,'o',label='input')
  plot(qro,zo,'-',label='interpolated')
  legend(frameon=False)
  
  subplot(236)
  plot(dqi,zi5,'o',label='input')
  plot(dqo,zo,'-',label='interpolated')
  legend(frameon=False)
