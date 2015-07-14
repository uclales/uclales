#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset
import h5py as H
from glob import glob
import sys
import os
from time import sleep
from datetime import datetime

#try:
#    from filelock import FileLock
#    have_lock=True
#except Exception,e:
#    print 'Couldnt import FileLock:',e
have_lock=False
complvl=1

#--------------------------------------------------------------------------------------------------------------------------------
def write_nc(basename,varname, data, dims, attributes=None):
  fname= '{0:}'.format(basename)+'.merged.nc'
  try:
  
    if have_lock:
        start = datetime.now()
        lock = FileLock(fname+'.lock')
        lock.acquire() # wait until file is ready for opening
        print("Lock acquired. waited :: {}".format(datetime.now()-start) )

    if os.path.exists(fname):
      fmode='a'
    else:
      fmode='w'
    D=Dataset(fname,fmode)
  
    for dim in dims:
      if dim[0] not in D.dimensions: 
        D.createDimension(dim[0], len(dim[1]) )
        if dim[0] not in D.variables:
            D.createVariable(dim[0] ,'f4',(dim[0],) )
            D.variables[dim[0]][:] = dim[1]
            print 'write_netcdf: write dim:',dim,'shape',np.shape(dim[1])
  
    if varname not in D.variables:
        print 'write_netcdf: var: ',varname,' shape', np.shape(data)
        D.createVariable(varname, 'f4', [ d[0] for d in dims ] , zlib=True,least_significant_digit=6, complevel=complvl)
        D.variables[varname][:] = data
        if attributes!=None:
            print 'attributes',attributes
            [ D.variables[varname].setncattr(att[0],att[1]) for att in attributes ]
  
  except Exception,e:
    print 'error occured when we tried writing data to netcdf file',e
  finally:
    if 'D' in locals(): D.close()
    if have_lock:
        lock.release()

def exists_nc(basename, varname):
  try:
    fname= basename+'.merged.nc'
    if os.path.exists(fname):
      fmode='r'
    else:
      exists = False
      print 'Checking if {0:} exists in file {1:} :: File not found => {2:}'.format(varname,fname,exists)
      return exists 

    D=Dataset(fname,fmode)

    if varname in D.variables:
      exists=True
    else:
      exists=False

  except Exception,e:
    print 'error occured when we checked if variable exists',e
    exists=False
  finally:
    if 'D' in locals(): D.close()

  print 'Checking if {0:} exists in file {1:} :: {2:}'.format(varname,basename,exists)
  return exists
#--------------------------------------------------------------------------------------------------------------------------------

#-------------------------------Reduction functions------------------------------------------------------------------------------
reduc_functions={ 
        0:{ # coordinates
        # coordinate:
        'time'  : None,
        'zt'    : None,
        'zm'    : None,
        'xt'    : None,
        'xm'    : None,
        'yt'    : None,
        'ym'    : None,
        },

        1:{ # one dimensional variables
        # .ts. variables:
        'cfl'    : np.max ,
        'maxdiv' : np.max ,
        'zi1_bar': np.mean,
        'zi2_bar': np.mean,
        'zi3_bar': np.mean,
        'vtke'   : np.mean,
        'sfcbflx': np.mean,
        'wmax'   : np.max ,
        'tsrf'   : np.mean, 
        'ustar'  : np.mean, 
        'shf_bar': np.mean, 
        'lhf_bar' : np.mean,
        'zi_bar'  : np.mean,
        'lwp_bar' : np.mean,
        'lwp_var' : np.mean,
        'zc'      : np.max ,
        'zb'      : np.min ,
        'cfrac'   : np.mean,
        'lmax'    : np.max ,
        'albedo'  : np.mean,
        'rwp_bar' : np.mean,
        'prcp'    : np.mean,
        'pfrac'   : np.mean,
        'CCN'     : np.mean,
        'nrain'   : np.mean,
        'nrcnt'   : np.sum ,
        'zcmn'    : np.mean,
        'zbmn'    : np.mean,
        'tkeint'  : np.mean,
        'lflxut'  : np.mean,
        'lflxdt'  : np.mean,
        'sflxut'  : np.mean,
        'sflxdt'  : np.mean,
        'thl_int' : np.mean,
        'wvp_bar' : np.mean,
        'wvp_var' : np.mean,
        'iwp_bar' : np.mean,
        'iwp_var' : np.mean,
        'swp_bar' : np.mean,
        'swp_var' : np.mean,
        'gwp_bar' : np.mean,
        'gwp_var' : np.mean,
        'hwp_bar' : np.mean,
        'hwp_var' : np.mean,
        'Qnet'    : np.mean,
        'G0'      : np.mean,
        'tndskin' : np.mean,
        'ra'      : np.mean,
        'rsurf'   : np.mean,
        'rsveg'   : np.mean,
        'rssoil'  : np.mean,
        'tskinav' : np.mean,
        'qskinav' : np.mean,
        'obl'     : np.mean,
        'cliq'    : np.mean,
        'a_Wl'    : np.mean,
        'lflxutc' : np.mean,
        'sflxutc' : np.mean,
        'tsair'   : np.mean,
        'sflxds'  : np.mean,
        'sflxus'  : np.mean,
        'lflxds'  : np.mean,
        'lflxus'  : np.mean,
        'sflxdsc' : np.mean,
        'sflxusc' : np.mean,
        'lflxdsc' : np.mean,
        'lflxusc' : np.mean,

        # .ps. variables:
        'fsttm'   : np.mean,
        'lsttm'   : np.mean,
        'nsmp'    : np.mean,
        
        # 3d vars
        'dn0'   : np.mean,
        'u0'    : np.mean,
        'v0'    : np.mean,
        },

        2 : { # 2d variables
        # .ps. variables:
        'dn0'     : np.mean,
        'u0'      : np.mean,
        'v0'      : np.mean,
        'u'       : np.mean,
        'v'       : np.mean,
        't'       : np.mean,
        'p'       : np.mean,
        'u_2'     : np.mean,
        'v_2'     : np.mean,
        'w_2'     : np.mean,
        't_2'     : np.mean,
        'w_3'     : np.mean,
        't_3'     : np.mean,
        'tot_tw'  : np.mean,  # Total vertical flux of theta -- TODO should this be sum?
        'sfs_tw'  : np.mean,
        'tot_uw'  : np.mean,
        'sfs_uw'  : np.mean,
        'tot_vw'  : np.mean,
        'sfs_vw'  : np.mean,
        'tot_ww'  : np.mean,
        'sfs_ww'  : np.mean,
        'km'      : np.mean,
        'kh'      : np.mean,
        'lmbd'    : np.mean,
        'lmbde'   : np.mean,
        'sfs_tke' : np.mean,
        'sfs_boy' : np.mean,
        'sfs_shr' : np.mean,
        'boy_prd' : np.mean,
        'shr_prd' : np.mean,
        'trans'   : np.mean,
        'diss'    : np.mean,
        'dff_u'   : np.mean,
        'dff_v'   : np.mean,
        'dff_w'   : np.mean,
        'adv_u'   : np.mean,
        'adv_v'   : np.mean,
        'adv_w'   : np.mean,
        'prs_u'   : np.mean,
        'prs_v'   : np.mean,
        'prs_w'   : np.mean,
        'prd_uw'  : np.mean,
        'storage' : np.mean,
        'q'       : np.mean,
        'q_2'     : np.mean,
        'q_3'     : np.mean,
        'tot_qw'  : np.mean,
        'sfs_qw'  : np.mean,
        'rflx'    : np.mean,
        'rflx2'   : np.mean,
        'sflx'    : np.mean,
        'sflx2'   : np.mean,
        'l'       : np.mean,
        'l_2'     : np.mean,
        'l_3'     : np.mean,
        'tot_lw'  : np.mean,
        'sed_lw'  : np.mean,
        'cs1'     : np.mean,
        'cnt_cs1' : np.mean,
        'w_cs1'   : np.mean,
        'tl_cs1'  : np.mean,
        'tv_cs1'  : np.mean,
        'rt_cs1'  : np.mean,
        'rl_cs1'  : np.mean,
        'wt_cs1'  : np.mean,
        'wv_cs1'  : np.mean,
        'wr_cs1'  : np.mean,
        'cs2'     : np.mean,
        'cnt_cs2' : np.mean,
        'w_cs2'   : np.mean,
        'tl_cs2'  : np.mean,
        'tv_cs2'  : np.mean,
        'rt_cs2'  : np.mean,
        'rl_cs2'  : np.mean,
        'wt_cs2'  : np.mean,
        'wv_cs2'  : np.mean,
        'wr_cs2'  : np.mean,
        'Nc'      : np.mean,
        'Nr'      : np.mean,
        'rr'      : np.mean,
        'prc_r'   : np.mean,
        'evap'    : np.mean,
        'frc_prc' : np.mean,
        'prc_prc' : np.mean,
        'frc_ran' : np.mean,
        'hst_srf' : np.mean,
        'lflxu'   : np.mean,
        'lflxd'   : np.mean,
        'sflxu'   : np.mean,
        'sflxd'   : np.mean,
        'cdsed'   : np.mean,
        'i_nuc'   : np.mean,
        'ice'     : np.mean,
        'n_ice'   : np.mean,
        'snow'    : np.mean,
        'graupel' : np.mean,
        'rsup'    : np.mean,
        'prc_c'   : np.mean,
        'prc_i'   : np.mean,
        'prc_s'   : np.mean,
        'prc_g'   : np.mean,
        'prc_h'   : np.mean,
        'hail'    : np.mean,
        'qt_th'   : np.mean,
        's_1'     : np.mean,
        's_2'     : np.mean,
        's_3'     : np.mean,
        'RH'      : np.mean,
        'lwuca'   : np.mean,
        'lwdca'   : np.mean,
        'swuca'   : np.mean,
        'swdca'   : np.mean,
        },

        4: { # 4d variables
        # .3d. variables:
        'u'     : np.concatenate,
        'v'     : np.concatenate,
        'w'     : np.concatenate,
        't'     : np.concatenate,
        'p'     : np.concatenate,
        'q'     : np.concatenate,
        'l'     : np.concatenate,
        'r'     : np.concatenate,
        'n'     : np.concatenate,
        'a_rhl'  : np.concatenate,
        'a_rhs'  : np.concatenate,
        'rflx'  : np.concatenate,
        'lflxu' : np.concatenate,
        'lflxd' : np.concatenate,
        },
        
        }

#--------------------------------------------------------------------------------------------------------------------------------


maxtime=-1
#--------------------------------------------------------------------------------------------------------------------------------
def append_var(basename,varname,reduc_func=np.mean):
  global maxtime
  if exists_nc(basename, varname):
    print "variable already exists",varname
    return

  files=glob(basename+'.0*.nc')

  if 'coord_files' not in locals(): coord_files={}
  for f in sorted(files):
        try:
          ident,coord,ending = [ i[::-1] for i in f[::-1].split('.',2) ] [::-1] # split filename but use maximum 2 splits to retain file identifier if it contains dots.
        except:
          print "Couldnt split filename into ''ident,coord,ending'' :: ",f
          return -1

        x,y = ( int(coord[:4]), int(coord[4:]) )
        if (x,y) not in coord_files.keys(): coord_files[(x,y)] = { 'fname':f, }
#        print 'Addind to coords:',x,y

  nr_y = len(np.unique ([ k[1] for k in coord_files.keys() ]))
  nr_x = len(np.unique ([ k[0] for k in coord_files.keys() ]))


  for i in np.arange(nr_x):
    for j in np.arange(nr_y):

      try:
          del(td)
          del(yd)
          del(xd)
          del(zd)
      except Exception,e:
          pass

      D = Dataset(coord_files[(i,j)]['fname'] )

      attributes = dict([ [att,D.variables[varname].getncattr(att)] for att in ['longname','units'] ])
      ndim=len(D.variables[varname].dimensions[:])

      if varname in reduc_functions[0]: # this is coordinate variable....
          return # dont want to save this

      if varname in reduc_functions[ndim]:
          reduc_func = reduc_functions[ndim][varname]
      else:
          print 'Need to supply reduction function vor variable',varname,'({}d)'.format(ndim)
          sys.exit(-1)
#      print 'Using reduc function',reduc_func.__name__,' for variable ',varname,'({}d)'.format(ndim)

#      print 'reading data from file:',coord_files[(i,j)]['fname'],' coords ',i,j

      l4d = ndim==4
      l3d = ndim==3
      l2d = ndim==2
      l1d = ndim==1

      if l4d: td,yd,xd,zd = D.variables[varname].dimensions[:]
      if l3d: td,yd,xd    = D.variables[varname].dimensions[:]
      if l2d: td,zd       = D.variables[varname].dimensions[:]
      if l1d: td,         = D.variables[varname].dimensions[:]

      if maxtime==-1:
        if td!='time': continue
        maxtime = len(D.variables[ td ][:])
        print 'maxtime is',maxtime

      if varname not in coord_files[(i,j)].keys(): 
        coord_files[(i,j)][varname] = D.variables[varname][:maxtime]

#      print "Coords of variable:",varname,'::',D.variables[varname].dimensions[:],[l1d,l2d,l3d,l4d], np.shape( D.variables[varname] ), np.shape(coord_files[(i,j)][varname])

      coord_files[td] = D.variables[ td ][:maxtime]

      try:
        coord_files[zd] = D.variables[ zd ][:]
      except: # if there is no z axis, just use the indices
        if l4d: coord_files[zd] = np.arange( np.shape( D.variables[varname] )[3] )
        if l2d: coord_files[zd] = np.arange( np.shape( D.variables[varname] )[1] )

      if 'yd' in locals():
        if '{0:}.{1:}'.format(yd,j) not in coord_files.keys(): 
          coord_files['{0:}.{1:}'.format(yd,j)] = D.variables[ yd ][:] #save y-dimension

      if 'xd' in locals():
        if '{0:}.{1:}'.format(xd,i) not in coord_files.keys():
          coord_files['{0:}.{1:}'.format(xd,i)] = D.variables[ xd ][:] #save x-dimension

      D.close()

    # append individual arrays in y dimension
    if l3d or l4d: coord_files['concat_{0:}'.format(i)] = reduc_func  ( [ coord_files[(i,j)].pop(varname) for j in np.arange(nr_y) ], axis=1 )
    if l2d or l1d: coord_files['concat_{0:}'.format(i)] = reduc_func  ( [ coord_files[(i,j)].pop(varname) for j in np.arange(nr_y) ], axis=0 )

  # append individual arrays in x dimension
  if l3d or l4d: var = reduc_func ( [ coord_files.pop('concat_{0:}'.format(i)) for i in np.arange(nr_x) ], axis=2)
  if l2d or l1d: var = reduc_func ( [ coord_files.pop('concat_{0:}'.format(i)) for i in np.arange(nr_x) ], axis=0)

  # append coordinate arrays for x and y axis
  if 'yd' in locals(): 
#    import ipdb;ipdb.set_trace()
    coord_files[yd] = np.concatenate( [ coord_files.pop('{0:}.{1:}'.format(yd,j)) for j in np.arange(nr_y) ], axis=1 ) 
  if 'xd' in locals(): 
    coord_files[xd] = np.concatenate( [ coord_files.pop('{0:}.{1:}'.format(xd,i)) for i in np.arange(nr_x) ], axis=2 ) 

  if l4d:
    data = var.swapaxes(2,3).swapaxes(1,2) # from [time,y,x,z] to [time,z,y,x]
    dims = [ [td,coord_files[td]], [zd,coord_files[zd]], [yd,coord_files[yd]], [xd,coord_files[xd]] ]
  if l3d:
    data = var.swapaxes(1,2) # from [time,y,x] to [time,y,x]
    dims = [ [td,coord_files[td]], [yd,coord_files[yd]], [xd,coord_files[xd]] ]
  if l2d:
    data = var
    dims = [ [td,coord_files[td]], [zd,coord_files[zd]] ]
  if l1d:
    data = var
    dims = [ [td,coord_files[td]], ]


  write_nc(basename,varname, data, dims, attributes=attributes)
#--------------------------------------------------------------------------------------------------------------------------------

try:
  basename = str( sys.argv[1] )
except:
  sys.exit(-1)

files = sorted(glob(basename+'.0*.nc'))
print "Opening files:",files

# Check for all possible variables:
vars=[]
D = Dataset( files[0], 'r' )
for v in D.variables:
    print 'Found Variable:',v.__str__()
    # if len(D.variables[v].dimensions)>=idim: 
    vars.append( v.__str__() )
D.close()

try:
    do_vars=[]
    varnames = None
    varnames = str( sys.argv[2] ).split()
    if 'all' not in varnames:
        for varname in varnames:
            if varname in vars:
                do_vars.append(varname)
            else:
                raise Exception('variable not found')
        vars=do_vars

except Exception,e:
    print ''
    print ''
    print "Error occurred when we tried to find the correct variable (",varname,"): ",e
    print "Possible variables are: 'all' or one of the following"
    print ' '.join(vars)
    print ''
    print ''
    sys.exit(0)


for idim in [4,3,2,1]:
  for varname in vars:
    append_var(basename,varname)
