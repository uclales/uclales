#
# automatically starts a poor man's parallelization for
# merging the particle output from UCLA-LES
# call by adding 'python auto_particlemerge.py' after 'poe ./uclales'
# no command-line input or change of the code is needed, automatically
# determines the case prefix and number of X,Y processes used.
# REQUIRES merge_particles in the same directory
#

import sys
import glob
import subprocess
import string
import os.path

# See if merge_particles exists, else quit.
if(os.path.isfile('merge_particles') != True):
  print 'merge_particles does not exists, STOPPING'
  sys.exit()

# Probe for prefix case
files = glob.glob('*.particles.0*')
pref = string.split(files[0],'.')[0]

# Get number of x&y processors
files.sort()
nxprocs = int(string.split(files[-1],'.')[2][:4])+1
nyprocs = int(string.split(files[-1],'.')[2][4:])+1

vars = (['x','y','z','u','v','w','t','tv','rt','rl'])

for var in vars:
  # Create jobscript
  jobscript = 'merge_' + var
  f = open(jobscript,'w')
  f.write('#!/client/bin/ksh \n')
  f.write('# @ shell = /client/bin/ksh \n')                      
  f.write('# @ job_type = serial \n')
  f.write('# @ job_name = mergeparticles_'+var+'.job \n')
  f.write('# @ output = $(job_name).o$(jobid) \n')
  f.write('# @ error = $(job_name).e$(jobid) \n')
  f.write('# @ queue \n')                                      
  f.write('./merge_particles ' + pref + ' ' + var + ' ' + str(nxprocs) + ' ' + str(nyprocs) + '\n') 
  f.close()

  # submit jobscript, wait for llsubmit to finish, delete jobscript
  p=subprocess.Popen("llsubmit "+jobscript,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  p.communicate()  # Wait for command to finish
  p=subprocess.Popen("rm "+jobscript,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
