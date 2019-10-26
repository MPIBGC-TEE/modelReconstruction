## script to convert Qing's old .nc files with a 
## sequence of cells into a gridded structure (lat, lon)


import numpy as np
import numpy.ma as ma

from netCDF4 import Dataset


#def serial2grid_vr(var_s, nr_layers, lat_args, lon_args):
#    """From Qing's old list of sites to a grid.
#    """
#    data = np.zeros(
#        (var_s.shape[2], nr_layers, len(lat_args), len(lon_args)), 
#        dtype = var_s.dtype
#    )
#    mask = np.ones_like(data)
#
#    for site_nr in range(len(lat_args)):
#        data[:,:,site_nr,site_nr] = \
#            var_s[site_nr,:nr_layers,:].transpose()
#        mask[:,:,site_nr,site_nr] = 0
#    data = data[:,:,lat_args,:]
#    data = data[:,:,:,lon_args]
#    mask = mask[:,:,lat_args,:]
#    mask = mask[:,:,:,lon_args]
#    mdata = ma.masked_array(data, mask)
#
#    return mdata


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
#    ds_g.createDimension('level', nr_layers)
#    ds_g.createDimension('lat', len_lat)
#    ds_g.createDimension('lon', len_lon)
#    
    ## time variable with bounds
    times_all = ds_all.variables['time']
    times = ds_small.createVariable('time', 'f8', ('time',))
    times.long_name = times_all.long_name
    times.units = times_all.units
    times.calendar = times_all.calendar
#    times.bounds = 'time_bounds'
    times[:] = times_all[:]

#    time_bounds = ds_g.createVariable(
#        'time_bounds', 
#        'f8', 
#        ('time','hist_interval')
#    )
#    time_bounds.long_name = 'history time interval endpoints'
#    time_bounds[...] = ma.masked_array(
#        ds_s.variables['time_bounds'][...], 
#        mask=0
#    )

    ## depth variable with bounds
    ds_depth = Dataset('../Data/DZSOI.nc', 'r')
    dz = ds_depth.variables['DZSOI']
    dz = dz[:nr_layers,0]
    ds_depth.close()
#    depth = np.ndarray((len(dz)+1,), dz.dtype)
#    depth[0] = 0
#    depth[1:] = np.cumsum(dz)
#    depths = ds_g.createVariable('depth', 'f4', ('depth',), fill_value=1e+36)
#    depths[:] = depth
#    depths.long_name = 'soil layer boundaries'
#    depths.units = 'm'
#    depths.bounds = 'depth_bounds'
    dzs = ds_small.createVariable('dz', 'f4', ('levdcmp',), fill_value=1e+36)
    dzs[:] = dz
    dzs.long_name = 'soil layer thickness'
    dzs.unit = 'm'


#    depth_bounds = ds_g.createVariable(
#        'depth_bounds', 
#        'f4', 
#        ('level', 'hist_interval')
#    )
#    depth_bounds.long_name = 'depth level endpoints'
#    depth_bnds = np.ndarray((nr_layers,2))
#    depth_bnds[:,0] = np.roll(depth,1)
#    depth_bnds[0,0] = 0.0
#    depth_bnds[:,1] = depth
#    depth_bounds[...] = depth_bnds

#    ## latidude and longitude
#    lat_s = ds_s['lat']
#    lat = ds_g.createVariable(
#        'lat', 
#        lat_s.dtype, 
#        ('lat',),
#        fill_value = lat_s._FillValue
#    )
#    lat.long_name = lat_s.long_name
#    lat.units = lat_s.units
#    ## sort lat
#    lat_args = np.argsort(lat_s)
#    lat[:] = lat_s[lat_args]
#
#    lon_s = ds_s['lon']
#    lon = ds_g.createVariable(
#        'lon', 
#        lon_s.dtype, 
#        ('lon',),
#        fill_value = lon_s._FillValue
#    )
#    lon.long_name = lon_s.long_name
#    lon.units = lon_s.units
#    ## sort lon
#
#    ##### first move lon into -180 to +180 #####
#    ## ILAMP.Variable.Variable.__init__ tries to do so but is flawed
#    lon_s = lon_s[:].copy()
#    lon_s = (lon_s<=180)*lon_s + (lon_s>180)*(lon_s-360) + (lon_s<-180)*360
#    ##### done #####
#    lon_args = np.argsort(lon_s)
#    lon[:] = lon_s[lon_args]

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
#                var_small = var_all
                var_small.long_name = var_all.long_name
                var_small.unit = var_all.units
                var_small.cell_methods = var_all.cell_methods
                print(var_small.shape, var_all.shape)
                var_small[...] = var_all[:,:nr_layers,:]

#                var_g[:] = serial2grid_vr(
#                    var_s, 
#                    nr_layers, 
#                    lat_args, 
#                    lon_args
#                )
                
                variable_name = f.readline()
                variable_name = variable_name.rstrip('\n')

    ## finish and save
    print('\n...done.\n') 
    ds_all.close()
    ds_small.close()
