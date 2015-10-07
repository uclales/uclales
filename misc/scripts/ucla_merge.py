#!/usr/bin/env python
import numpy as np
from netCDF4 import Dataset
#import h5py as H
from glob import glob
import sys
import os
from time import sleep
from datetime import datetime

def nanaverage(a,weights):
    indices = ~np.isnan(a)
    return np.average(a[indices], weights=weights[indices])

#try:
#    from filelock import FileLock
#    have_lock=True
#except Exception,e:
#    print 'Couldnt import FileLock:',e
have_lock=False
complvl=3

#--------------------------------------------------------------------------------------------------------------------------------
def load_ncvar(basename,varname):
  fname= '{0:}'.format(basename)+'.merged.nc'
  try:
    D=Dataset(fname,'r')
    return D.variables[varname][:]
    
  except Exception,e:
    print 'Could not load variable {0:} from merged nc file because {1:}'.format(varname,e)
    raise(e)

  finally:
    D.close()
  
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
  
    try:
        for dim in dims:
          if dim[0] not in D.dimensions: 
            D.createDimension(dim[0], len(dim[1]) )
            if dim[0] not in D.variables:
                D.createVariable(dim[0] ,'f4',(dim[0],) )
                D.variables[dim[0]][:] = dim[1]
                print 'write_netcdf: write dim:',dim,'shape',np.shape(dim[1])
    except Exception,e:
        print 'Error happened writing dimensions',e
        raise(e)
  
    try:
        if varname not in D.variables:
            print 'write_netcdf: var: ',varname,' shape', np.shape(data),' nc dimensions:',[ d[0] for d in dims ]
            fillval = attributes.pop('_FillValue') if attributes!=None else None

            data[np.isnan(data)] = fillval

            D.createVariable(varname, 'f4', [ d[0] for d in dims ] , zlib=True,least_significant_digit=6, complevel=complvl,fill_value=fillval)
            D.variables[varname][:] = data

            try:
                if attributes!=None:
                    print 'attributes',attributes
                    [ D.variables[varname].setncattr(att,val) for att,val in attributes.iteritems() ]
            except Exception,e:
                print 'Error happened writing attributes',e
                raise(e)
    except Exception,e:
        print 'Error happened writing data',e
        raise(e)
  
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
        'cfl'     : np.nanmax ,
        'maxdiv'  : np.nanmax ,
        'zi1_bar' : np.nanmean,
        'zi2_bar' : np.nanmean,
        'zi3_bar' : np.nanmean,
        'vtke'    : np.nanmean,
        'sfcbflx' : np.nanmean,
        'wmax'    : np.nanmax ,
        'tsrf'    : np.nanmean, 
        'ustar'   : np.nanmean, 
        'shf_bar' : np.nanmean, 
        'lhf_bar' : np.nanmean,
        'zi_bar'  : np.nanmean,
        'lwp_bar' : np.nanmean,
        'lwp_var' : np.nanmean,
        'zc'      : np.nanmax ,
        'zb'      : np.nanmin ,
        'cfrac'   : np.nanmean,
        'lmax'    : np.nanmax ,
        'albedo'  : np.nanmean,
        'rwp_bar' : np.nanmean,
        'prcp'    : np.nanmean,
        'pfrac'   : np.nanmean,
        'CCN'     : np.nanmean,
        'nrain'   : np.nanmean,
        'nrcnt'   : np.nansum ,
        'zcmn'    : np.nanmean,
        'zbmn'    : np.nanmean,
        'tkeint'  : np.nanmean,
        'lflxut'  : np.nanmean,
        'lflxdt'  : np.nanmean,
        'sflxut'  : np.nanmean,
        'sflxdt'  : np.nanmean,
        'thl_int' : np.nanmean,
        'wvp_bar' : np.nanmean,
        'wvp_var' : np.nanmean,
        'iwp_bar' : np.nanmean,
        'iwp_var' : np.nanmean,
        'swp_bar' : np.nanmean,
        'swp_var' : np.nanmean,
        'gwp_bar' : np.nanmean,
        'gwp_var' : np.nanmean,
        'hwp_bar' : np.nanmean,
        'hwp_var' : np.nanmean,
        'Qnet'    : np.nanmean,
        'G0'      : np.nanmean,
        'tndskin' : np.nanmean,
        'ra'      : np.nanmean,
        'rsurf'   : np.nanmean,
        'rsveg'   : np.nanmean,
        'rssoil'  : np.nanmean,
        'tskinav' : np.nanmean,
        'qskinav' : np.nanmean,
        'obl'     : np.nanmean,
        'cliq'    : np.nanmean,
        'a_Wl'    : np.nanmean,
        'lflxutc' : np.nanmean,
        'sflxutc' : np.nanmean,
        'tsair'   : np.nanmean,
        'sflxds'  : np.nanmean,
        'sflxus'  : np.nanmean,
        'lflxds'  : np.nanmean,
        'lflxus'  : np.nanmean,
        'sflxdsc' : np.nanmean,
        'sflxusc' : np.nanmean,
        'lflxdsc' : np.nanmean,
        'lflxusc' : np.nanmean,

        # .ps. variables:
        'fsttm'   : np.nanmean,
        'lsttm'   : np.nanmean,
        'nsmp'    : np.nanmean,
        
        # 3d vars
        'dn0'   : np.nanmean,
        'u0'    : np.nanmean,
        'v0'    : np.nanmean,
        },

        2 : { # 2d variables
        # .ps. variables:
        'dn0'     : np.nanmean,
        'u0'      : np.nanmean,
        'v0'      : np.nanmean,
        'u'       : np.nanmean,
        'v'       : np.nanmean,
        't'       : np.nanmean,
        'p'       : np.nanmean,
        'u_2'     : np.nanmean,
        'v_2'     : np.nanmean,
        'w_2'     : np.nanmean,
        't_2'     : np.nanmean,
        'w_3'     : np.nanmean,
        't_3'     : np.nanmean,
        'tot_tw'  : np.nanmean,  # Total vertical flux of theta -- TODO should this be sum?
        'sfs_tw'  : np.nanmean,
        'tot_uw'  : np.nanmean,
        'sfs_uw'  : np.nanmean,
        'tot_vw'  : np.nanmean,
        'sfs_vw'  : np.nanmean,
        'tot_ww'  : np.nanmean,
        'sfs_ww'  : np.nanmean,
        'km'      : np.nanmean,
        'kh'      : np.nanmean,
        'lmbd'    : np.nanmean,
        'lmbde'   : np.nanmean,
        'sfs_tke' : np.nanmean,
        'sfs_boy' : np.nanmean,
        'sfs_shr' : np.nanmean,
        'boy_prd' : np.nanmean,
        'shr_prd' : np.nanmean,
        'trans'   : np.nanmean,
        'diss'    : np.nanmean,
        'dff_u'   : np.nanmean,
        'dff_v'   : np.nanmean,
        'dff_w'   : np.nanmean,
        'adv_u'   : np.nanmean,
        'adv_v'   : np.nanmean,
        'adv_w'   : np.nanmean,
        'prs_u'   : np.nanmean,
        'prs_v'   : np.nanmean,
        'prs_w'   : np.nanmean,
        'prd_uw'  : np.nanmean,
        'storage' : np.nanmean,
        'q'       : np.nanmean,
        'q_2'     : np.nanmean,
        'q_3'     : np.nanmean,
        'tot_qw'  : np.nanmean,
        'sfs_qw'  : np.nanmean,
        'rflx'    : np.nanmean,
        'rflx2'   : np.nanmean,
        'sflx'    : np.nanmean,
        'sflx2'   : np.nanmean,
        'l'       : np.nanmean,
        'l_2'     : np.nanmean,
        'l_3'     : np.nanmean,
        'tot_lw'  : np.nanmean,
        'sed_lw'  : np.nanmean,
        'cs1'     : np.nanmean,
        'cnt_cs1' : np.nansum,
        'w_cs1'   : np.nansum, 
        'tl_cs1'  : np.nansum, 
        'tv_cs1'  : np.nansum, 
        'rt_cs1'  : np.nansum, 
        'rl_cs1'  : np.nansum, 
        'wt_cs1'  : np.nansum, 
        'wv_cs1'  : np.nansum, 
        'wr_cs1'  : np.nansum, 
        'cs2'     : np.nanmean,
        'cnt_cs2' : np.nansum,
        'w_cs2'   : np.nansum, 
        'tl_cs2'  : np.nansum, 
        'tv_cs2'  : np.nansum, 
        'rt_cs2'  : np.nansum, 
        'rl_cs2'  : np.nansum, 
        'wt_cs2'  : np.nansum, 
        'wv_cs2'  : np.nansum, 
        'wr_cs2'  : np.nansum, 
        'Nc'      : np.nanmean,
        'Nr'      : np.nanmean,
        'rr'      : np.nanmean,
        'prc_r'   : np.nanmean,
        'evap'    : np.nanmean,
        'frc_prc' : np.nanmean,
        'prc_prc' : np.nanmean,
        'frc_ran' : np.nanmean,
        'hst_srf' : np.nanmean,
        'lflxu'   : np.nanmean,
        'lflxd'   : np.nanmean,
        'sflxu'   : np.nanmean,
        'sflxd'   : np.nanmean,
        'cdsed'   : np.nanmean,
        'i_nuc'   : np.nanmean,
        'ice'     : np.nanmean,
        'n_ice'   : np.nanmean,
        'snow'    : np.nanmean,
        'graupel' : np.nanmean,
        'rsup'    : np.nanmean,
        'prc_c'   : np.nanmean,
        'prc_i'   : np.nanmean,
        'prc_s'   : np.nanmean,
        'prc_g'   : np.nanmean,
        'prc_h'   : np.nanmean,
        'hail'    : np.nanmean,
        'qt_th'   : np.nanmean,
        's_1'     : np.nanmean,
        's_2'     : np.nanmean,
        's_3'     : np.nanmean,
        'RH'      : np.nanmean,
        'lwuca'   : np.nanmean,
        'lwdca'   : np.nanmean,
        'swuca'   : np.nanmean,
        'swdca'   : np.nanmean,


        },

        3: { # 3d variables
        # .2d. variables:
        'shf'     : np.concatenate,
        'lhf'     : np.concatenate,
        'ustars'  : np.concatenate,
        'a_tskin' : np.concatenate,
        'a_qskin' : np.concatenate,
        'a_Qnet'  : np.concatenate,
        'a_Qnet'  : np.concatenate,
        'a_G0'    : np.concatenate,
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
        'rice'  : np.concatenate,
        'nice'  : np.concatenate,
        'rsnow' : np.concatenate,
        'rgrp'  : np.concatenate,
        'nsnow' : np.concatenate,
        'ngrp'  : np.concatenate,
        'rhail' : np.concatenate,
        'nhail' : np.concatenate,
        'n'     : np.concatenate,
        'a_rhl' : np.concatenate,
        'a_rhs' : np.concatenate,
        'rflx'  : np.concatenate,
        'lflxu' : np.concatenate,
        'lflxd' : np.concatenate,

        'tsoil' : np.concatenate,
        'phiw'  : np.concatenate,
        },
        
        }

#--------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------
maxtime=-1
def append_var(basename,varname,reduc_func=np.mean):
  global maxtime
  if exists_nc(basename, varname):
    print "variable already exists",varname
    return
  else:
    print 'merging new variable',varname

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

      attributes = dict([ [att,D.variables[varname].getncattr(att)] for att in D.variables[varname].ncattrs() ])
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
        if td!='time':
          pass
        else:
          maxtime = len(D.variables[ td ][:])
          print 'maxtime is',maxtime

#      import ipdb;ipdb.set_trace()

      # set FillValues to NaN -- reduc functions will not care em then
      if type(D.variables[varname][:]) is np.ma.core.MaskedArray:
        valid = lambda x: np.where(x[:].mask==False,x[:].data,np.NaN)
      else:
        valid = lambda x: x

      # Save data into coord_files construct -- will later merge these
      data = valid( D.variables[varname][:maxtime] )

      longname = D.variables[varname].longname

      if varname not in coord_files[(i,j)].keys(): 
          coord_files[(i,j)][varname] = data

      # Save coordinates for variable:
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

      #print 'shape single variable',np.shape(coord_files[(i,j)][varname])

    #print 'shape single variable',np.shape( [ coord_files[(i,j)][varname] for j in np.arange(nr_y) ])

    # append individual arrays in y dimension
    if l3d or l4d: coord_files['concat_{0:}'.format(i)] = reduc_func  ( [ coord_files[(i,j)].pop(varname) for j in np.arange(nr_y) ], axis=1 )
    if l2d or l1d: coord_files['concat_{0:}'.format(i)] = reduc_func  ( [ coord_files[(i,j)].pop(varname) for j in np.arange(nr_y) ], axis=0 )

    #print 'shape yval', np.shape(coord_files['concat_{0:}'.format(i)])

  # append individual arrays in x dimension
  if l3d or l4d: var = reduc_func ( [ coord_files.pop('concat_{0:}'.format(i)) for i in np.arange(nr_x) ], axis=2)
  if l2d or l1d: var = reduc_func ( [ coord_files.pop('concat_{0:}'.format(i)) for i in np.arange(nr_x) ], axis=0)

  #print 'shape xval', np.shape(var)

  # append coordinate arrays for x and y axis
  if 'yd' in locals(): 
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

  # conditionally sampled data needs to be renormalized to global number of number of conditionally sampled data
  def renormalize_cs(data, longname, csvarname, basename):
    if ('over {0:}'.format(csvarname) in longname) or ('and {0:}'.format(csvarname) in longname):
      append_var(basename,'cnt_{0:}'.format(csvarname))
      Ncs = load_ncvar(basename,'cnt_{0:}'.format(csvarname)) # number of cs with dim: (time,zt)
      return data / Ncs
    else:
      return data


  data = renormalize_cs(data, longname, 'cs1', basename)
  data = renormalize_cs(data, longname, 'cs2', basename)

  write_nc(basename,varname, data, dims, attributes=attributes)
#--------------------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
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
        D = Dataset( files[0], 'r' )
        ndim = len(np.shape(D.variables[varname]))
        D.close()
        if ndim == idim:
            append_var(basename,varname)
