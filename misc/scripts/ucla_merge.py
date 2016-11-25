#!/usr/bin/env python
"""
Collects ucla-LES output data from multiple netCDF files
"""

import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from contextlib import closing

if "nanmean" not in dir(np):
    def nanmean(arr, *args, **kwargs):
        """ compute mean value ignoring nan entries """
        return np.nansum(arr, *args, **kwargs) / np.sum(np.isfinite(arr), *args, **kwargs)
    np.nanmean = nanmean

REDUCTION_FUNCTIONS = [
    (np.nanmax, ['cfl', 'lmax', 'maxdiv', 'wmax', 'zc']),
    (np.nanmean, [
        'CCN', 'G0', 'Nc', 'Nr', 'Qnet', 'RH', 'a_Wl', 'adv_u',
        'adv_v', 'adv_w', 'albedo', 'boy_prd', 'cdsed', 'cfrac',
        'cliq', 'cs1', 'cs2', 'dff_u', 'dff_v', 'dff_w', 'diss',
        'dn0', 'evap', 'frc_prc', 'frc_ran', 'fsttm', 'graupel',
        'gwp_bar', 'gwp_var', 'hail', 'hst_srf', 'hwp_bar', 'hwp_var',
        'i_nuc', 'ice', 'iwp_bar', 'iwp_var', 'kh', 'km', 'l', 'l_2',
        'l_3', 'lflxd', 'lflxds', 'lflxdsc', 'lflxdt', 'lflxu', 'lflxus',
        'lflxusc', 'lflxut', 'lflxutc', 'lhf_bar', 'lmbd', 'lmbde',
        'lsttm', 'lwdca', 'lwp_bar', 'lwp_var', 'lwuca', 'n_ice', 'nrain',
        'nsmp', 'obl', 'p', 'pfrac', 'prc_c', 'prc_g', 'prc_h', 'prc_i',
        'prc_prc', 'prc_r', 'prc_s', 'prcp', 'prd_uw', 'prs_u', 'prs_v',
        'prs_w', 'q', 'q_2', 'q_3', 'qskinav', 'qt_th', 'ra', 'rflx', 'rflx2',
        'rr', 'rssoil', 'rsup', 'rsurf', 'rsveg', 'rwp_bar', 's_1', 's_2',
        's_3', 'sed_lw', 'sfcbflx', 'sflx', 'sflx2', 'sflxd', 'sflxds',
        'sflxdsc', 'sflxdt', 'sflxu', 'sflxus', 'sflxusc', 'sflxut',
        'sflxutc', 'sfs_boy', 'sfs_qw', 'sfs_shr', 'sfs_tke', 'sfs_tw',
        'sfs_uw', 'sfs_vw', 'sfs_ww', 'shf_bar', 'shr_prd', 'snow', 'storage',
        'swdca', 'swp_bar', 'swp_var', 'swuca', 't', 't_2', 't_3', 'thl_int',
        'tkeint', 'tndskin', 'tot_lw', 'tot_qw', 'tot_tw', 'tot_uw', 'tot_vw',
        'tot_ww', 'trans', 'tsair', 'tskinav', 'tsrf', 'u', 'u0', 'u_2',
        'ustar', 'v', 'v0', 'v_2', 'vtke', 'w_2', 'w_3', 'wvp_bar', 'wvp_var',
        'zbmn', 'zcmn', 'zi1_bar', 'zi2_bar', 'zi3_bar', 'zi_bar',
        'time', 'xd', 'xm', 'xt', 'ym', 'yt', 'zm', 'zt']),
    (np.nanmin, ['zb']),
    (np.nansum, [
        'cnt_cs1', 'cnt_cs2', 'nrcnt', 'rl_cs1', 'rl_cs2', 'rt_cs1', 'rt_cs2',
        'tl_cs1', 'tl_cs2', 'tv_cs1', 'tv_cs2', 'w_cs1', 'w_cs2', 'wr_cs1',
        'wr_cs2', 'wt_cs1', 'wt_cs2', 'wv_cs1', 'wv_cs2'])
]

REDUCTION_FUNCTIONS_REV = {field: func
                           for func, fields in REDUCTION_FUNCTIONS
                           for field in fields}


def valid(array_like):
    """
    Returns an array or masked array as plain array
    with NaN on invalid fields.
    """
    if isinstance(array_like, np.ma.core.MaskedArray):
        return np.where(array_like.mask, np.NaN, array_like.data)
    else:
        return array_like


def calculate_transposition_rule(decomposition_dimensions, dimensions):
    """
    Calculates how a combined array of many subdomains with decomposition
    dimensions as first axes must be transposed such that a reshape will
    result in correct concatenation of the decomposed dimensions.

    Example:
    If an array with ``dimensions`` ('x', 'y', 'z') would have been decomposed
    into ('x', 'y'), the x axis of the subdomain must be swapped with the
    y axis of the decomposition. The correct transposition rule would be
    therefore: (0, 2, 1, 3, 4).
    """
    subdomain_positions = []
    for i, dim in enumerate(decomposition_dimensions):
        try:
            subdomain_positions.append((dimensions.index(dim), -1))
        except ValueError:
            subdomain_positions.append((-1, i))  # will not be concatenated
    dimension_positions = [(i, 0) for i in xrange(len(dimensions))]
    new_order = zip(
        *sorted(enumerate(subdomain_positions + dimension_positions),
                key=lambda x: x[1]))[0]
    print 'new_order:', new_order
    return new_order


def calculate_concatenated_shape(decomposition_dimensions,
                                 decomposition_shape,
                                 dimensions,
                                 shape):
    """
    Calculates the resulting shape of the concatenation of an array decomposed
    into subdomains along ``decomposition_dimensions`` into
    ``decomposition_shape`` subdomains.
    """
    print decomposition_dimensions, decomposition_shape
    remaining_size = reduce(lambda a, b: a * b,
                            [size
                             for name, size in zip(decomposition_dimensions,
                                                   decomposition_shape)
                             if name not in dimensions],
                            1)
    rescale_dict = dict(zip(decomposition_dimensions, decomposition_shape))
    new_shape = tuple(size * rescale_dict.get(name, 1)
                      for name, size in zip(dimensions, shape))
    print 'new size:', remaining_size, new_shape
    if remaining_size == 1:
        return new_shape
    else:
        return (remaining_size,) + new_shape


def get_nc_attrs(var):
    """
    Get all attributes of a netCDF4 variable
    """
    return {name: var.getncattr(name) for name in var.ncattrs()}


class NetCDFCollector(object):

    """
    Collects data from local subdomains
    """
    decomposition_dimensions = [
        ('xt', 'yt'),
        ('xm', 'ym')]
    subdomain_id_size = 4

    def __init__(self, basename, complevel=3, mintime=None, maxtime=None):
        self.basename = basename
        self.complevel = complevel
        self.files = sorted(glob(basename + '.0*.nc'))
        with closing(Dataset(self.files[0], 'r')) as dataset:
            self.variables = {name: {"shape": var.shape,
                                     "dimensions": var.dimensions,
                                     "dtype": var.dtype,
                                     "attributes": get_nc_attrs(var)}
                              for name, var in dataset.variables.items()}
        self.subdomain_shape = tuple(
            max(pos) + 1 for pos in zip(*map(self.file_coords, self.files)))
        ntimes = []
        for filename in self.files:
            with closing(Dataset(filename, 'r')) as dataset:
                try:
                    ntimes.append(len(dataset.variables['time']))
                except KeyError:
                    pass

        self.mintime = mintime
        if len(ntimes) > 0:
            self.maxtime = np.minimum(min(ntimes), maxtime) if maxtime!=None else min(ntimes)
        else:
            self.maxtime = None

        print 'time slice is from', self.mintime, 'to', self.maxtime

        subdomain_size = {dim: size
                          for dims in self.decomposition_dimensions
                          for dim, size in zip(dims, self.subdomain_shape)}
        self.dimensions = {
            dim: size * subdomain_size.get(dim, 1)
            for var in self.variables.values()
            for dim, size in zip(var["dimensions"], var["shape"])
        }
        openmode = 'a' if os.path.exists(self.merged_filename) else 'w'
        with closing(Dataset(self.merged_filename, openmode)) as dataset:
            for dim, size in self.dimensions.items():
                if dim not in dataset.dimensions:
                    dataset.createDimension(dim, size=size)

    @property
    def merged_filename(self):
        """
        The name of the merged netCDF file
        """
        return '%s.merged.nc' % self.basename

    def variable_exists(self, varname):
        """
        Checks if variable exists in output netCDF
        """
        with closing(Dataset(self.merged_filename)) as dataset:
            return varname in dataset.variables

    def file_coords(self, filename):
        """
        Calculates subdomain coordinates from filename
        """
        coord = list(reversed(filename.split('.')))[1]
        char_count = self.subdomain_id_size
        dimension_count = len(self.decomposition_dimensions)
        return tuple(map(int, [coord[i*char_count:(i+1)*char_count]
                               for i in xrange(dimension_count)]))

    def find_decomposition_dimenstions(self, dimensions):
        """
        Finds the name of the dimensions on which a variable
        defined on ``dimensions`` is decomposed.
        """
        possible_dimensions = map(set, zip(*self.decomposition_dimensions))
        dimensions = set(dimensions)

        dec_dimensions = []
        for i, dims in enumerate(possible_dimensions):
            dimset = dims & dimensions
            if len(dimset) > 1:
                raise ValueError('bad dimensions: %s' % str(dimensions))
            if len(dimset) == 0:
                dec_dimensions.append('__dummy_dimension_%d' % i)
            else:
                dec_dimensions.append(dimset.pop())
        return tuple(dec_dimensions)

    def load_variable(self, varname):
        """
        Loads a variable from all source files,
        can be indexed by subdomain position
        """
        dimensions = self.variables[varname]['dimensions'] # list of dim names
        slices = [ slice(None) for _ in dimensions ]       # create list of slices as [:,:,...]
        try:
            slices[dimensions.index("time")] = slice(self.mintime, self.maxtime, None) # restrict time slice to maxtime
        except ValueError: # dont have time axis, e.g. axis variables
            pass
        slices = tuple(slices)

        data = {}
        for filename in self.files:
            with closing(Dataset(filename, 'r')) as dataset:
                data[self.file_coords(filename)] = valid(
                    dataset.variables[varname][slices])
        return data

    def concatenate_data(self, varname, data):
        """
        Concatenates data as loaded by ``load_variable`` into one array.

        .. note::

            All subdomain axes which are not contained in the variable will be
            concatenated into an additional zeroth axis.
        """
        var_info = self.variables[varname]

        tmp_shape = list(var_info["shape"])
        for i,dimname in enumerate(var_info["dimensions"]):
            if dimname=='time':
                start = 0 if self.mintime==None else self.mintime
                tmp_shape[i] = self.maxtime - start
        tmp_shape = tuple(tmp_shape)

        intermediate_shape = self.subdomain_shape + tmp_shape
        print 'intermediate shape',intermediate_shape
        temp = np.zeros(intermediate_shape, dtype=var_info["dtype"])

        decomposition_dimensions = self.find_decomposition_dimenstions(
            var_info["dimensions"])

        for pos in data.keys():
            temp[pos] = data.pop(pos) # immediately remove data to save some much needed mem

        temp = temp.transpose(
            calculate_transposition_rule(decomposition_dimensions,
                                         var_info["dimensions"]))

        new_shape = calculate_concatenated_shape(decomposition_dimensions,
                                                  self.subdomain_shape,
                                                  var_info["dimensions"],
                                                  tmp_shape)
        temp = temp.reshape(new_shape)

        return temp

    def collect_variable(self,
                         varname,
                         reduction_function=None,
                         skip_if_existent=False):
        """
        Collects one variable ``varname`` from all netCDF-Files available.
        """
        if skip_if_existent and self.variable_exists(varname):
            print('Skipping {} bc. already exists'.format(varname))
            return
        print "collecting variable %s" % varname
        data = self.load_variable(varname)
        data = self.concatenate_data(varname, data)
        print "before reduction:", varname, data.shape
        need_reduction = len(data.shape) != len(
            self.variables[varname]["shape"])
        if need_reduction:
            if reduction_function is None:
                reduction_function = REDUCTION_FUNCTIONS_REV[varname]
            data = reduction_function(data, axis=0)
        print "after reduction:", varname, data.shape

        with closing(Dataset(self.merged_filename, 'a')) as dataset:
            attributes = self.variables[varname]["attributes"]
            fillvalue = attributes.get('_FillValue', -999)
            try:
                var = dataset.variables[varname]
            except KeyError:
                var = dataset.createVariable(
                    varname,
                    self.variables[varname]["dtype"],
                    self.variables[varname]["dimensions"],
                    fill_value=fillvalue,
                    complevel=self.complevel)
            for attr_name, attr_value in attributes.items():
                if attr_name == '_FillValue':
                    continue
                setattr(var, attr_name, attr_value)
            data[np.isnan(data)] = fillvalue

            # Get slice considering mintime and maxtime
            var_info = self.variables[varname]
            tmp_slices = [slice(None,None,None)]*len(var_info['dimensions'])
            for i,dimname in enumerate(var_info["dimensions"]):
                if dimname=='time': 
                    tmp_slices[i] = slice(self.mintime, self.maxtime, None)

            var[tmp_slices] = data


        print "done collecting variable %s" % varname

        # conditionally sampled data needs to be renormalized to global number
        # of number of conditionally sampled data
        # def renormalize_cs(data, longname, csvarname, basename):
        #     if ('over {0:}'.format(csvarname) in longname)
        #             or ('and {0:}'.format(csvarname) in longname):
        #         append_var(basename,'cnt_{0:}'.format(csvarname))
        #         # number of cs with dim: (time,zt)
        #         Ncs = load_ncvar(basename,'cnt_{0:}'.format(csvarname))
        #         return data / Ncs
        #     else:
        #         return data

        # data = renormalize_cs(data, longname, 'cs1', basename)
        # data = renormalize_cs(data, longname, 'cs2', basename)

    def collect(self, variables=None, skip_if_existent=False):
        """
        Collect all or selected ``variables``.
        """
        variable_names = set(self.variables)
        print 'Found Variables:', variable_names

        if variables is None:
            selected_variables = variable_names
        else:
            selected_variables = variables
            missing_variables = selected_variables - variable_names
            if len(missing_variables) > 0:
                raise ValueError('the following variables are missing: %s',
                                 str(missing_variables))

        map( lambda x: self.collect_variable(x, skip_if_existent=skip_if_existent),
            sorted(selected_variables,
                key=lambda x: len(self.variables[x]["shape"])))


def _main():
    import argparse

    parser = argparse.ArgumentParser(description='Merge uclales output files.')
    parser.add_argument('basename', type=str,
                        help='basename of the netCDF files')
    parser.add_argument('variables', type=str, nargs='+',
                        help='variables to extract (all: all of them)')
    parser.add_argument('-mintime', type=int, default=None,
                        help='minimum time index that is being loaded')
    parser.add_argument('-maxtime', type=int, default=None,
                        help='maximum time index that is being loaded')
    parser.add_argument('--skip_existent', action='store_true', default=False,
                        help='skip variable if it already exists in outfile')
    args = parser.parse_args()

    collector = NetCDFCollector(args.basename, mintime=args.mintime, maxtime=args.maxtime)
    selected_variables = set(args.variables)
    if 'all' in selected_variables:
        collector.collect(skip_if_existent=args.skip_existent)
    else:
        collector.collect(selected_variables, skip_if_existent=args.skip_existent)

if __name__ == '__main__':
    _main()
