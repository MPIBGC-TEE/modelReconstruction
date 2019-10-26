## Select only the necessary dimensions and variables from 
## '..Data/holger.new3sites.nc', take 'dz' from 'DZSOI.net', bring everything to
## 'JAMES.nc'

import numpy as np
import numpy.ma as ma

from netCDF4 import Dataset


if __name__ == '__main__':
    filename_all = '../Data/holger.new3sites.nc'
    filename_small = '../Data/JAMES.nc'
    nr_layers = 10
    len_lat = 3
    len_lon = 3

    ds_all = Dataset(filename_all, 'r')
    ds_small = Dataset(filename_small, 'w')

    ## copy global attributes
    attrs_all = ['title', 'comment', 'Conventions', 'history', 'source',
             'hostname', 'username', 'version', 'revision_id', 'case_title',
             'case_id', 'Surface_dataset', 'Initial_conditions_dataset',
             'PFT_physiological_constants_dataset',
             'NCO', 'nco_openmp_thread_number']
#    attrs_minimal = []
    for attr in attrs_all:
        setattr(ds_small, attr, ds_all.getncattr(attr))

    ## copy dimensions
    for dim_name in ['lndgrid', 'time']:
        dim = ds_all.dimensions[dim_name]
        ds_small.createDimension(dim.name, dim.size)

    ## create new dimensions
    ds_small.createDimension('levdcmp', nr_layers)
    
    ## time variable 
    times_all = ds_all.variables['time']
    times = ds_small.createVariable('time', 'f8', ('time',))
    times.long_name = times_all.long_name
    times.units = times_all.units
    times.calendar = times_all.calendar
    times[:] = times_all[:]

    ## depth variable 
    ds_depth = Dataset('../Data/DZSOI.nc', 'r')
    dz = ds_depth.variables['DZSOI']
    dz = dz[:nr_layers,0]
    ds_depth.close()
    dzs = ds_small.createVariable('dz', 'f4', ('levdcmp',), fill_value=1e+36)
    dzs[:] = dz
    dzs.long_name = 'soil layer thickness'
    dzs.unit = 'm'


    ## stocks and fluxes
    print('Converting stocks and fluxes to gridded structure...\n')

    filenames = ['../Data/variable_names_stocks.txt', 
                 '../Data/variable_names_fluxes.txt']
    for file_nr, filename in enumerate(filenames):
        with open(filename, 'r') as f:
            variable_name = f.readline()
            variable_name = variable_name.rstrip('\n')
            while variable_name:
                var_all = ds_all.variables[variable_name]
                print(variable_name)
                print(var_all.dimensions)
                var_small = ds_small.createVariable(
                    variable_name, 
                    var_all.dtype, 
                    ('lndgrid', 'levdcmp', 'time'),
                    fill_value = var_all._FillValue
                )
                var_small.long_name = var_all.long_name
                var_small.unit = var_all.units
                var_small.cell_methods = var_all.cell_methods
                print(var_small.shape, var_all.shape)
                var_small[...] = var_all[:,:nr_layers,:]

                variable_name = f.readline()
                variable_name = variable_name.rstrip('\n')

    ## finish and save
    print('\n...done.\n') 
    ds_all.close()
    ds_small.close()
