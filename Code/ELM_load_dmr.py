import numpy as np
import xarray as xr

from tqdm import tqdm

import discrete_model_run as DMR

################################################################################
#
# functions to load and prepare data, create discrete model run
#
################################################################################


def adapt_data(ds, cell_nr):
    ## pools
    pool_types = ['CWDC',
                  'LITR1C', 'LITR2C', 'LITR3C',
                  'SOIL1C', 'SOIL2C', 'SOIL3C']

    layers = [str(ly) for ly in range(10)]
    pools = []
    for pool_type in pool_types:
        pools += [pool_type + ly for ly in layers]

    nr_layers = len(layers)
    nr_pools = len(pools)

    # convert time unit of data from nanoseconds to days
    seconds_per_day = 86400
    data_times = np.array(ds.time, dtype='datetime64[D]')
    data_times_diff = np.array([1.0]*(len(data_times)-1))

    ## soil layer thickneses in m
    # have to get rid of coordinates, otherwise xarray won't multiply
    # because coords do not match --> make numpy objects out of it
    dzsoi_xr = ds.dz[:nr_layers]
    dzsoi = np.zeros_like(dzsoi_xr, dtype=np.float32)
    dzsoi[:] = dzsoi_xr[:]

    def convert_gCpm3_to_gCpm2(stock_data):
        # conversion from gC/m^3 to gC/m^2
        sd = np.zeros_like(stock_data, dtype=np.float32)
        sd[:] = stock_data[:]
        sdc = sd * dzsoi

        return sdc

    def convert_gCpm3ps_to_gCpm2pts(flux_data):
        # conversion from gC/m^3/s to gC/m^2/d
        fd = np.zeros_like(flux_data, dtype=np.float32)
        fd[:] = flux_data[:]
        fdc = fd * dzsoi * seconds_per_day

        # multiply with time difference
        flux_data_converted = (fdc[:,:].T * data_times_diff).T
        
        return flux_data_converted

    def convert_gCpm2ps_to_gCpm2pts(flux_data):
        # conversion from gC/m^3/s to gC/m^2/d
        fd = np.zeros_like(flux_data, dtype=np.float32)
        fd[:] = flux_data[:]
        fdc = fd * seconds_per_day

        # multiply with time difference
        flux_data_converted = (fdc[:,:].T * data_times_diff).T
        
        return flux_data_converted

    ## xs, stocks
    print('\nextracting C stocks data')
    xs = np.zeros((len(data_times), nr_pools))
    for pt_nr, pool_type in enumerate(tqdm(pool_types)):
        key = pool_type + '_vr'
        pool_type_vals = ds[key][cell_nr, :nr_layers, :].transpose()
        xs[:,(pt_nr*nr_layers):((pt_nr+1)*nr_layers)] = \
            convert_gCpm3_to_gCpm2(pool_type_vals)
    
#    key = 'LITR1C_vr'
#    pool_type_vals = ds[key][9,:nr_layers,cell_nr] * dzsoi
#    print('x0', pool_type_vals)
#    pool_type_vals = ds[key][10,:nr_layers,cell_nr] * dzsoi
#    print('x1', pool_type_vals)
#    input()

    ## Fs, internal fluxes
    def get_pool_nr(pool_name):
        return pools.index(pool_name)

    ## horizontal fluxes
    print('\nextracting decomposition fluxes')
    flux_names_dic = {('CWDC', 'LITR2C'): 'CWDC_TO_LITR2C_vr',
                      ('CWDC', 'LITR3C'): 'CWDC_TO_LITR3C_vr',
                      ('LITR1C', 'SOIL1C'): 'LITR1C_TO_SOIL1C_vr',
                      ('LITR2C', 'SOIL1C'): 'LITR2C_TO_SOIL1C_vr',
                      ('LITR3C', 'SOIL2C'): 'LITR3C_TO_SOIL2C_vr',
                      ('SOIL1C', 'SOIL2C'): 'SOIL1C_TO_SOIL2C_vr',
                      ('SOIL1C', 'SOIL3C'): 'SOIL1C_TO_SOIL3C_vr',
                      ('SOIL2C', 'SOIL1C'): 'SOIL2C_TO_SOIL1C_vr',
                      ('SOIL2C', 'SOIL3C'): 'SOIL2C_TO_SOIL3C_vr',
                      ('SOIL3C', 'SOIL1C'): 'SOIL3C_TO_SOIL1C_vr'}

    Fs = np.zeros((len(data_times)-1,nr_pools,nr_pools)) 
    for source_type, dest_type in tqdm(flux_names_dic.keys()):
        flux_name = flux_names_dic[(source_type, dest_type)]
        flux_data = ds[flux_name][cell_nr,:nr_layers,1:].transpose()
        flux_data_converted = convert_gCpm3ps_to_gCpm2pts(flux_data)

        for ly_nr, ly in enumerate(layers):
            j = get_pool_nr(source_type + ly)
            i = get_pool_nr(dest_type + ly)
            Fs[:,i,j] = flux_data_converted[:,ly_nr]

#    for ly in range(10):
#        print('F', ly, Fs[9,40+ly,10+ly])
#    input()
    
    # helper function for rs, us
    def collect_fluxes(flux_names_dic):
        fluxes = np.zeros((len(data_times)-1, nr_pools))

        for pt_nr, pool_type in enumerate(tqdm(pool_types)):
            if pool_type in flux_names_dic.keys():
                flux_names = flux_names_dic[pool_type]       
                flux_data = np.zeros((len(data_times)-1, nr_layers))
    
                for flux_name in flux_names:
                    data = ds[flux_name][cell_nr,:nr_layers,1:].transpose()
                    flux_data += data
    
                fluxes[:,(pt_nr*nr_layers):((pt_nr+1)*nr_layers)] \
                    = convert_gCpm3ps_to_gCpm2pts(flux_data)

        return fluxes

    ## rs, external outputs
    print('\nextracting external outputs')
    flux_names_dic = {'CWDC': ['CWD_HR_L2_vr', 'CWD_HR_L3_vr'],
                      'LITR1C': ['LITR1_HR_vr'],
                      'LITR2C': ['LITR2_HR_vr'],
                      'LITR3C': ['LITR3_HR_vr'],
                      'SOIL1C': ['SOIL1_HR_S2_vr', 'SOIL1_HR_S3_vr'],
                      'SOIL2C': ['SOIL2_HR_S1_vr', 'SOIL2_HR_S3_vr'],
                      'SOIL3C': ['SOIL3_HR_vr']}
    rs = collect_fluxes(flux_names_dic)
#    print('r', rs[9,10:20])
#    input()


    ## us, external inputs
    print('\nextracting external inputs')
    flux_names_dic = {'CWDC': ['fire_mortality_c_to_cwdc_col',
                               'gap_mortality_c_to_cwdc_col',
                               'harvest_c_to_cwdc_col'],
                      'LITR1C': ['gap_mortality_c_to_litr_met_c_col',
                                 'harvest_c_to_litr_met_c_col',
                                 'm_c_to_litr_met_fire_col',
                                 'phenology_c_to_litr_met_c_col'],
                      'LITR2C': ['gap_mortality_c_to_litr_cel_c_col',
                                 'harvest_c_to_litr_cel_c_col',
                                 'm_c_to_litr_cel_fire_col',
                                 'phenology_c_to_litr_cel_c_col'],
                      'LITR3C': ['gap_mortality_c_to_litr_lig_c_col',
                                 'harvest_c_to_litr_lig_c_col',
                                 'm_c_to_litr_lig_fire_col',
                                 'phenology_c_to_litr_lig_c_col']}
    us = collect_fluxes(flux_names_dic)
#    print('u', us[9,10:20])
#    input()


    ## vertical fluxes
    print('\nextracting vertical fluxes')

    for pt in tqdm(pool_types):
        pt_vert = pt[:-1] # fluxes have no C at the end: LITR1 instead of LITR1C
        down_flux_out_name = pt_vert + '_diffus_down'
        up_flux_in_name = pt_vert + '_diffus_up'

        down_flux_out_vr = convert_gCpm2ps_to_gCpm2pts(ds[down_flux_out_name][cell_nr,:nr_layers,1:].transpose())
        up_flux_in_vr = convert_gCpm2ps_to_gCpm2pts(ds[up_flux_in_name][cell_nr,:nr_layers,1:].transpose())


        ## Qing changed the code, fluxes might be negative, then in opposite direction
#        # losses to below
#        for ly_nr, ly in enumerate(layers[:-1]):
#            j = get_pool_nr(pt + ly)
#            i = j + 1
#            Fs[:,i,j] = down_flux_out_vr[:,ly_nr]
#
#        # leaching, fixme: get rid of it?
#        #j = get_pool_nr(pt + layers[-1])
#        #rs[:,j] += down_flux_out_vr[:,-1]
#
#        # gain from below
#        for ly_nr, ly in enumerate(layers[:-1]):
#            i = get_pool_nr(pt + ly)
#            j = i + 1
#            Fs[:,i,j] = up_flux_in_vr[:,ly_nr] 

        ## here's the try to correct for Qing's changes
    
        # losses to below
        for ly_nr, ly in enumerate(layers[:-1]):
            j = get_pool_nr(pt + ly)
            i = j + 1

            fluxes = down_flux_out_vr[:,ly_nr]
            for ti in range(fluxes.shape[0]):
                flux = fluxes[ti]
                if flux >=0 :
                    Fs[ti,i,j] = flux
                else:
                    Fs[ti,j,i] = -flux

        # gain from below
        for ly_nr, ly in enumerate(layers[:-1]):
            i = get_pool_nr(pt + ly)
            j = i + 1

            fluxes = up_flux_in_vr[:,ly_nr]
            for ti in range(fluxes.shape[0]):
                flux = fluxes[ti]
                if flux >=0:
                    Fs[ti,i,j] = flux
                else:
                    Fs[ti,j,i] = -flux


#    for ly in range(10):
#        print('v_down', ly, Fs[9,10+ly+1,10+ly])
#    input()

#    for ly in range(10):
#        print('v_up', ly, Fs[9,10+ly,10+ly+1])
#    input()

#    ## vertical fluxes, CWD does not have any            
#
#    # we need to compute it recursively as I did in my notes!
#    # option 1: use the TNDNCY fluxes (which are actually concentration 
#    #           changes, multiply by the layer thickness and obtain the fluxes
#    #           the top layer's flux can only go downwards, so it is obtained 
#    #           immediately, the layer below can now use this known additional 
#    #           input to compute its own (downward) flux, etc.
#
#    # option 2:  do not use TNDCY, but rather use the x values of
#    #            the next time step, find it below
#
#    # option 1 seems far more precise than option 2
#    for pool_type in tqdm(pool_types[1:]):
#        TNDNCY_flux_name = pool_type + '_TNDNCY_VERT_TRANS'
#        TNDNCY_flux_data = ds[TNDNCY_flux_name][1:,:nr_layers,cell_nr]
#        TNDNCY_flux_data_converted = convert_gCpm3ps_to_gCpm2pts(
#                                        TNDNCY_flux_data
#                                     )
#        flux_data_converted = np.zeros_like(TNDNCY_flux_data_converted)
#        flux_data_converted[:,0] = TNDNCY_flux_data_converted[:,0]
#        for ly_nr in range(1, nr_layers-1):
#            flux_data_converted[:,ly_nr] = \
#                TNDNCY_flux_data_converted[:,ly_nr] + \
#                flux_data_converted[:,ly_nr-1]
#            
#        for ly_nr, ly in enumerate(layers[:-1]):
#            for t, data in enumerate(flux_data_converted[:,ly_nr]):
#                #print(data)
#                #print(type(data))
#                if data < 0:
#                    # mass moves one layer down
#                    j = get_pool_nr(pool_type + ly)
#                    i = j + 1
#                    data = -data # fluxes are required to be nonnegative
#                else:
#                    # mass comes from one layer above
#                    i = get_pool_nr(pool_type + ly)
#                    j = i - 1
#
#                Fs[t,i,j] = data
#
#    # option 2 for vertical fluxes
#
#    new_Fs = Fs.copy()
#    for t in tqdm(range(len(Fs))):
#        for pt in pool_types[1:]:
#            old_data = 0
#            data = 0
#            for ly_nr, ly in enumerate(layers[:-1]):
#                pool = get_pool_nr(pt + ly)
#                #print('pt', pt, 'ly', ly, 'pool', pool)
#
#                inp_h = Fs[t,pool,:].sum(0) + us[t,pool]
#                #print('Fs_inp', Fs[t,pool,:]) 
#                out_h = Fs[t,:,pool].sum(0) + rs[t,pool]
#                #print('Fs_out', Fs[t,:,pool])
#                net_h = inp_h-out_h
#                #print('ion', inp_h, out_h, net_h)
#
#                data = xs[t+1,pool] - xs[t,pool] - net_h + old_data
#                #print('xd', xs[t+1,pool], xs[t,pool], old_data, data)
#                old_data = data
#                if data < 0:
#                    # mass moves one layer down
#                    j = pool
#                    i = j + 1
#                    data_abs = -data # fluxes are required to be nonnegative
#                else:
#                    # mass comes from one layer above
#                    i = pool
#                    j = i - 1
#                    data_abs = data
#
#                new_Fs[t,i,j] = data_abs
#    Fs = new_Fs.copy()


    ## xs_14C, 14C solutions
    print('\nextracting 14C solution data')
    xs_14C = np.zeros((len(data_times), nr_pools))
    for pt_nr, pool_type in enumerate(tqdm(pool_types)):
        key = 'C14_' + pool_type + '_vr'
        pool_type_vals = ds[key][cell_nr, :nr_layers, :].transpose()
        xs_14C[:,(pt_nr*nr_layers):((pt_nr+1)*nr_layers)] = \
            convert_gCpm3_to_gCpm2(pool_type_vals)
    

#    ## rs_14C, not vertically resolved, not by pool only total
    print('\nextracting external 14C outputs')
#    def convert_gCpm2ps_to_gCpm2pts(flux_data):
#        fdc = flux_data * seconds_per_day
#
#        # multiply with time difference
#        flux_data_converted = (fdc * data_times_diff).T
#        
#        return flux_data_converted
#
#    rs_14C_xr = convert_gCpm2ps_to_gCpm2pts(ds['C14_HR'][1:,cell_nr])
#    rs_14C = np.zeros_like(rs_14C_xr)
#    rs_14C[:] = rs_14C_xr[:]

    ## rs, external outputs
    flux_names_dic_14C = {'CWDC': ['CWD_C14HR_L2_vr', 'CWD_C14HR_L3_vr'],
                      'LITR1C': ['LITR1_C14HR_vr'],
                      'LITR2C': ['LITR2_C14HR_vr'],
                      'LITR3C': ['LITR3_C14HR_vr'],
                      'SOIL1C': ['SOIL1_C14HR_S2_vr', 'SOIL1_C14HR_S3_vr'],
                      'SOIL2C': ['SOIL2_C14HR_S1_vr', 'SOIL2_C14HR_S3_vr'],
                      'SOIL3C': ['SOIL3_C14HR_vr']}
#    rs = collect_fluxes(flux_names_dic)
    rs_14C = collect_fluxes(flux_names_dic_14C)


    ## us_14C, external inputs
    print('\nextracting external 14C inputs')
    flux_names_dic_14C = {'CWDC': ['C14_fire_mortality_c_to_cwdc_col',
                                   'C14_gap_mortality_c_to_cwdc_col',
                                   'C14_harvest_c_to_cwdc_col'],
                      'LITR1C': ['C14_gap_mortality_c_to_litr_met_c_col',
                                 'C14_harvest_c_to_litr_met_c_col',
                                 'C14_m_c_to_litr_met_fire_col',
                                 'C14_phenology_c_to_litr_met_c_col'],
                      'LITR2C': ['C14_gap_mortality_c_to_litr_cel_c_col',
                                 'C14_harvest_c_to_litr_cel_c_col',
                                 'C14_m_c_to_litr_cel_fire_col',
                                 'C14_phenology_c_to_litr_cel_c_col'],
                      'LITR3C': ['C14_gap_mortality_c_to_litr_lig_c_col',
                                 'C14_harvest_c_to_litr_lig_c_col',
                                 'C14_m_c_to_litr_lig_fire_col',
                                 'C14_phenology_c_to_litr_lig_c_col']}

    us_14C = collect_fluxes(flux_names_dic_14C)

    ## bring 14C stocks to the same order of magnitude as the 12C stocks
    xs_14C *= 1e12
#    print(xs_14C)
#    input()
    rs_14C *= 1e12
    us_14C *= 1e12

    return data_times, xs, Fs, rs, us, xs_14C, rs_14C, us_14C

def aggregate_days_to_time_step(days_per_ts, xs, Fs, rs, us, \
            xs_14C, rs_14C, us_14C):
    # combine days_per_ts days into one time step
    xs_agg = xs[np.arange(0, len(xs), days_per_ts)]
    Fs_agg = np.add.reduceat(Fs, np.arange(0, len(Fs), days_per_ts))
    Fs_agg = Fs_agg[0:(len(xs_agg)-1)]
    rs_agg = np.add.reduceat(rs, np.arange(0, len(rs), days_per_ts))
    rs_agg = rs_agg[0:(len(xs_agg)-1)]
    us_agg = np.add.reduceat(us, np.arange(0, len(us), days_per_ts))
    us_agg = us_agg[0:(len(xs_agg)-1)]

    xs_14C_agg = xs_14C[np.arange(0, len(xs_14C), days_per_ts)]
    rs_14C_agg = np.add.reduceat(rs_14C, np.arange(0, len(rs_14C), days_per_ts))
    rs_14C_agg = rs_14C_agg[0:(len(xs_14C_agg)-1)]
    us_14C_agg = np.add.reduceat(us_14C, np.arange(0, len(us_14C), days_per_ts))
    us_14C_agg = us_14C_agg[0:(len(xs_14C_agg)-1)]

    return xs_agg, Fs_agg, rs_agg, us_agg, xs_14C_agg, rs_14C_agg, us_14C_agg


def cut_off_delay(delay_in_ts, times_agg, xs_agg, Fs_agg, rs_agg, us_agg, \
            xs_14C_agg, rs_14C_agg, us_14C_agg):
    ## cut off the first delay time steps
    times_delayed = times_agg[delay_in_ts:]

    xs_agg_delayed = xs_agg[delay_in_ts:]
    Fs_agg_delayed = Fs_agg[delay_in_ts:]
    rs_agg_delayed = rs_agg[delay_in_ts:]
    us_agg_delayed = us_agg[delay_in_ts:]

    xs_14C_agg_delayed = xs_14C_agg[delay_in_ts:]
    rs_14C_agg_delayed = rs_14C_agg[delay_in_ts:]
    us_14C_agg_delayed = us_14C_agg[delay_in_ts:]

    return times_delayed, xs_agg_delayed, Fs_agg_delayed, rs_agg_delayed, \
           us_agg_delayed, xs_14C_agg_delayed, rs_14C_agg_delayed, \
            us_14C_agg_delayed

def compute_cache(ds, cell_nr, days_per_ts, delay_in_ts):
    data_times, xs, Fs, rs, us, xs_14C, rs_14C, us_14C = \
        adapt_data(ds, cell_nr)

    times = data_times[np.arange(0, len(data_times), days_per_ts)]

    xs, Fs, rs, us, xs_14C, rs_14C, us_14C = aggregate_days_to_time_step(
        days_per_ts, 
        xs, Fs, rs, us, 
        xs_14C, rs_14C, us_14C
    )

    times, xs, Fs, rs, us, xs_14C, rs_14C, us_14C = cut_off_delay(
        delay_in_ts, 
        times,
        xs, Fs, rs, us, 
        xs_14C, rs_14C, us_14C
    )

    cache = {'times': times,
             'xs': xs, 
             'Fs': Fs, 
             'rs': rs, 
             'us': us,
             'xs_14C': xs_14C, 
             'rs_14C': rs_14C,
             'us_14C': us_14C}

    return cache

def load_discrete_model_run(ds, days_per_ts, max_sim_length,
        delay_in_ts, cell_nr, pools):


    cache = compute_cache(ds, cell_nr, days_per_ts, delay_in_ts)

    ## unpack cache dictionary
    times = cache['times']
    xs = cache['xs']
    Fs = cache['Fs']
    rs = cache['rs']
    us = cache['us']
    xs_14C = cache['xs_14C']
    rs_14C = cache['rs_14C']
    us_14C = cache['us_14C']

    ## cut from the end to correct simulation length
    sim_length = max_sim_length
    times = times[:sim_length+1]
    xs = xs[:sim_length+1, pools]
    Fs = Fs[:sim_length, pools][:,:,pools]
    rs = rs[:sim_length, pools]
    us = us[:sim_length, pools]
    xs_14C = xs_14C[:sim_length+1, pools]
    rs_14C = rs_14C[:sim_length, pools]
    us_14C = us_14C[:sim_length, pools]

    ## finally initialize the model run
    start_values = xs[0]
    dmr = DMR.DiscreteModelRun.reconstruct_from_data(
        times,
        start_values,
        xs,
        Fs,
        rs,
        us
    )

    return dmr, xs, xs_14C, us_14C

