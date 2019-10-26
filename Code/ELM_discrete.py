## Repdrocude the ELM plots of JAMES manuscript
##
## Title:        Mathematical reconstruction of land carbon models from their
##               numerical output: computing soil radiocarbon from 12C dynamics
##
## Authors:      Holger Metzler, Qing Zhu, William Riley, Alison Hoyt,
##               Markus Müller, Carlos A. Sierra 
##
## Code:         Holger Metzler (10/2019)
##
## Requirements: matplotlib
##               numpy
##               pandas
##               tqdm
##               xarray
##               scipy
##               discrete_model_run.py (included)
##               ELM_load_dmr.py (included)              
##
## Input:        '../Data/JAMES.nc' (ELM model run data, maybe unzip first)
##               '../Data/C14Atm_NH.csv' (Atmospheric 14C values)
##
## Output:       To subfolder 'Output"
##
## Remark:       The codes seems to load the models unnecessarily often.
##               I had to do it this way because of memory restrictions on my 
##               machine.


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d

from ELM_load_dmr import load_discrete_model_run


## return all the pool numbers that belong to one type (e.g. CWD)
def get_pools_in_pool_type(pt, pools, layers):
    sp = pt*len(layers)
    fp = min(sp+len(layers), len(pools))
    pt_pools = pools[sp:fp]
    return pt_pools


## return all the pool numbers that belong to one layer 
def get_pools_in_layer(ly, pools, layers):
    ly_pools = np.arange(ly, len(pools), len(layers))
    return ly_pools


def load_dmr_and_dmr_14C(ds, cell_nr, cell_name, pools, parameter_values):
    print('\nPreparing %s model...' % cell_name)
    dmr, xs, xs_14C, us_14C = load_discrete_model_run(
                ds, 
                parameter_values['days_per_ts'], 
                parameter_values['max_sim_length'], 
                parameter_values['delay_in_ts'], 
                cell_nr, 
                pools
            )
    dmr.name = cell_name
    dmr.xs = xs
    dmr.xs_14C = xs_14C
    dmr.us_14C = us_14C

    ## create 14C DMR
    start_values_14C = dmr.xs_14C[0,:]
    dmr_14C = dmr.to_14C_only(
        start_values_14C, 
        dmr.us_14C, 
        decay_rate = parameter_values['decay_rate']
    )
    dmr_14C.name = cell_name + '_14C'

    return dmr, dmr_14C


################################################################################


def plot_depth_profile(ds, cell_names, pools, parameter_values, filename):
    plot_times = np.array([0,365*70,365*110]) // days_per_ts
    fig = plt.figure(figsize=(7*len(cell_names), 10)) # (width, height)

    fontsize = 40
    for ax_nr, cell_nr in enumerate([2,0,1]):
        cell_name = cell_names[cell_nr]
        dmr, dmr_14C = load_dmr_and_dmr_14C(
            ds, 
            cell_nr, 
            cell_name, 
            pools, 
            parameter_values
        )
        soln = dmr.solve() 
        soln_14C = dmr_14C.solve()
        xs = dmr.xs
        xs_14C = dmr.xs_14C
    
        ## plot Delta 14C depth profile
        dzsoi = ds['dz']
        depths = np.cumsum(dzsoi)

        ax = fig.add_subplot(1, len(cell_names), ax_nr + 1)
        ax.text(0.9,0.1,
                cell_name, 
        #        pad = 10,
                fontdict = {'fontsize': 20, 'fontweight': 'bold'},
                transform = ax.transAxes,
                horizontalalignment = 'right'
        )
        ax.set_xlim([-1100,200])
        ax.set_ylim([-4.0,0.0])

        for nr, time in enumerate(plot_times):
            year = int(1900 + time // (365/days_per_ts))

            Delta_14C = np.zeros_like(np.arange(nr_layers))
            Delta_14C_ELM = np.zeros_like(np.arange(nr_layers))
            for ly in range(nr_layers):
                ly_pools = np.arange(nr_pool_types)*nr_layers + ly

                ratio = soln_14C[time, ly_pools].sum(0)\
                        /soln[time, ly_pools].sum(0)
                Delta_14C[ly] = (ratio-1)*1000
        
                ratio_ELM = xs_14C[time, ly_pools].sum(0)\
                        /xs[time, ly_pools].sum(0)
                Delta_14C_ELM[ly] = (ratio_ELM-1)*1000
        
            c = ['orange', 'green', 'blue'][nr]
            if (ax_nr == 0):
                ax.plot(
                    Delta_14C, 
                    -depths, 
                    'o', 
                    color = c, 
                    label = 'DTA, ' + str(year), 
                    markersize=5
                )
                ax.plot(
                    Delta_14C_ELM, 
                    -depths, 
                    '-', 
                    color = c, 
                    label = 'ELMv1-ECA, ' + str(year), 
                    markersize=5
                )
            else:
                ax.plot(
                    Delta_14C, 
                    -depths, 
                    'o', 
                    color = c, 
                    label = '_', 
                    markersize=5
                )
                ax.plot(
                    Delta_14C_ELM, 
                    -depths, 
                    '-', 
                    color = c, 
                    label = '_', 
                    markersize=5
                )
            
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(r'$\Delta{}^{14}$C (‰)', fontsize=40, labelpad=20)

        if ax_nr == 0: ax.set_ylabel('depth (m)', fontsize=40)
        
        #for item in [ax.title, ax.xaxis.label, ax.yaxis.label]:
        #    item.set_fontsize(14)
    
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(18)
        
        if ax_nr == 0:
            ax.legend(fontsize=18, loc=2)
    
    fig.savefig(filename, tight_layout=True)
    plt.close(fig)


def print_errors(ds, cell_names, pools, parameter_values, filename):
    aggregate_pool_types = [np.array([0]), 
                            np.array([1,2,3]),
                            np.array([4,5,6])]
    aggregate_pool_type_names = ['CWD', 'LITTER', 'SOIL']


    with open(filename, 'w') as f:
        for cell_nr, cell_name in enumerate(cell_names):
            dmr, dmr_14C = load_dmr_and_dmr_14C(
                ds, 
                cell_nr, 
                cell_name, 
                pools, 
                parameter_values
            )
            soln = dmr.solve() 
            soln_14C = dmr_14C.solve()
            xs = dmr.xs
            xs_14C = dmr.xs_14C
    
            f.write(dmr.name + '\n' + '-'*len(dmr.name) + '\n')
    
            order = [2,0,1]
            for pt_nr in order:
                agg_pt = aggregate_pool_types[pt_nr]
                agg_pools_l = []
                for pt in agg_pt:
                    agg_pools_l.append(
                        np.arange(pt*nr_layers, (pt+1)*nr_layers, 1)
                    )
    
                agg_pools = np.array(agg_pools_l).flatten()
                pt_agg_name = aggregate_pool_type_names[pt_nr]
    
                C12_agg = soln[:,agg_pools].sum(1)
                C12_agg_ELM = xs[:,agg_pools].sum(1)
                diff = C12_agg-C12_agg_ELM
                rel_err = (diff/C12_agg_ELM) * 100
                rel_err_abs_C12 = np.abs(diff/C12_agg_ELM) * 100

                line = 'C12 mean absolute relative error (%s): %f' \
                        % (pt_agg_name, np.mean(rel_err_abs_C12))
                f.write(line + '\n')
    
                C14_agg = soln_14C[:,agg_pools].sum(1)
                C14_agg_ELM = xs_14C[:, agg_pools].sum(1)
                diff_14C = C14_agg-C14_agg_ELM
                rel_err_14C = (diff_14C/C14_agg_ELM) * 100
                rel_err_abs_14C = np.abs(diff_14C/C14_agg_ELM) * 100

                line = 'C14 mean absolute relative error (%s): %f' \
                        % (pt_agg_name, np.mean(rel_err_abs_14C))
                f.write(line + '\n\n')
        

def plot_14C_rel_err(ds, cell_nr, cell_name, pools, parameter_values, filename):
    aggregate_pool_types = [np.array([0]), 
                            np.array([1,2,3]),
                            np.array([4,5,6])]
    aggregate_pool_type_names = ['CWD', 'LITTER', 'SOIL']

    fontsize = 20
    fig_14C_rel_err = plt.figure(figsize=(12,7)) # (width, height)
    
    ## plot 14C through time
    dmr, dmr_14C = load_dmr_and_dmr_14C(
        ds, 
        cell_nr, 
        cell_name, 
        pools, 
        parameter_values
    )
    soln = dmr.solve() 
    soln_14C = dmr_14C.solve()
    xs_14C = dmr.xs_14C

    ax_14C_rel_err = fig_14C_rel_err.add_subplot(1, 1, 1)
    ax_14C_rel_err.set_title(dmr.name, 
                 pad = -30,
                 fontdict = {'fontsize': 30, 'fontweight': 'bold'}
                )

    order = [2,0,1]
    colors = ['orange', 'green', 'blue']
    for nr in order:
        agg_pt = aggregate_pool_types[nr]
        agg_pools_l = []
        for pt in agg_pt:
            agg_pools_l.append(np.arange(pt*nr_layers, (pt+1)*nr_layers, 1)) 

        agg_pools = np.array(agg_pools_l).flatten()
        pt_agg_name = aggregate_pool_type_names[nr]

        C14_agg = soln_14C[:,agg_pools].sum(1)
        C14_agg_ELM = xs_14C[:, agg_pools].sum(1)
        diff_14C = C14_agg-C14_agg_ELM
        rel_err_14C = (diff_14C/C14_agg_ELM) * 100
            
        ax_14C_rel_err.plot(
            dmr.times, 
            rel_err_14C, 
            '-', 
            label = pt_agg_name, 
            c=colors[nr]
        )

        ax_14C_rel_err.set_xlim([dmr.times[0], dmr.times[-1]])
        ax_14C_rel_err.set_xticks(np.array([1920,1940,1960,1980, 2000])*365+150)
        ax_14C_rel_err.set_ylim([-0.25, 0.15])
        ax_14C_rel_err.set_ylabel(r'${}^{14}$C relative error ($\%$)')
        ax_14C_rel_err.legend(fontsize=14, loc=2)
#        ax_14C_rel_err.set_xlabel('time (yr)')

        for item in [
            ax_14C_rel_err.title, 
            ax_14C_rel_err.xaxis.label, 
            ax_14C_rel_err.yaxis.label
        ]:
            item.set_fontsize(fontsize)
    
        for item in (
            ax_14C_rel_err.get_xticklabels() + ax_14C_rel_err.get_yticklabels()
        ):
            item.set_fontsize(fontsize)

    fig_14C_rel_err.savefig(filename, tight_layout=True)
    plt.close(fig_14C_rel_err)


def plot_Delta_14C_through_time(ds, cell_names, pools,\
                                parameter_values, filename):
    fig_pt = plt.figure(figsize=(12,7*len(cell_names))) # (width, height)

    fontsize = 20
    for ax_nr, cell_nr in enumerate([2,0,1]):
        cell_name = cell_names[cell_nr]
        dmr, dmr_14C = load_dmr_and_dmr_14C(
            ds, 
            cell_nr, 
            cell_name, 
            pools, 
            parameter_values
        )
        soln = dmr.solve() 
        soln_14C = dmr_14C.solve()
        xs = dmr.xs
        xs_14C = dmr.xs_14C

        ## Delta 14C per pool type over time
        ax_pt = fig_pt.add_subplot(len(cell_names), 1, ax_nr + 1)
        ax_pt.set_title(cell_names[cell_nr], 
                     pad = -30,
                     fontdict = {'fontsize': 40, 'fontweight': 'bold'}
                    )
        ax_pt.axhline(c='black', lw=0.5)

        ## plot atmospheric Delta 14C
        atm_delta_14C = np.loadtxt(
            '../Data/C14Atm_NH.csv', 
            skiprows=1, 
            delimiter=','
        )
        Fa_func = interp1d(
            atm_delta_14C[:,0], 
            atm_delta_14C[:,1], 
            fill_value='extrapolate')
        days = pd.DatetimeIndex(dmr.times)
        days_as_float = days.year+(days.dayofyear-1)/365

        ax_pt.plot(
            dmr.times, 
            [Fa_func(days_as_float[i]) for i in range(len(dmr.times))], 
            label='Atmosphere')

        ## plot pools
        for pt in pool_types:
            pt_pools = get_pools_in_pool_type(pt, pools, layers)
            
            # delta 14C, ELM
            ratio = xs_14C[:,pt_pools].sum(1)/xs[:,pt_pools].sum(1)
            delta_14C_ELM = (ratio-1)*1000
            line, = ax_pt.plot(dmr.times, delta_14C_ELM, label=pool_names[pt])

            # delta 14C, my computation
            color = line.get_color()
            ratio = soln_14C[:,pt_pools].sum(1)/soln[:,pt_pools].sum(1)
            delta_14C = (ratio-1)*1000
            ax_pt.plot(dmr.times, 
                delta_14C, 
                label='_',
                dashes=[1,4],
                lw=5,
                c=color
                )

        ax_pt.set_xlim([dmr.times[0], dmr.times[-1]])
        ax_pt.set_xticks(np.array([1920,1940,1960,1980, 2000])*365+150)
        ax_pt.set_ylabel(r'$\Delta^{14}$C$\,$(‰)')

        ax_pt.set_ylim([-800,1200])

        if ax_nr == 0:
            ax_pt.legend(fontsize=14, loc=2)

        if ax_nr == 2:
            ax_pt.set_xlabel('time (yr)')

        for item in [ax_pt.title, ax_pt.xaxis.label, ax_pt.yaxis.label]:
            item.set_fontsize(fontsize)
    
        for item in (ax_pt.get_xticklabels() + ax_pt.get_yticklabels()):
            item.set_fontsize(fontsize)

    fig_pt.savefig(filename, tight_layout=True)
    plt.close(fig_pt)


################################################################################
#
# Main program
#
################################################################################


if __name__ == '__main__':
    ## we have 70 pools:
    ##  00-09: 10 layers of pool type CWD
    ##  ...
    ##  60-69: 10 layers of pool type SOIL3C    
    pools = np.arange(70)
    pool_types = np.arange(7)
    pool_names = ['CWD', 
                  'LITTER 1', 'LITTER 2', 'LITTER 3',
                  'SOIL 1', 'SOIL 2', 'SOIL 3']

    nr_pool_types = len(pool_types)
    layers = range(len(pools)//len(pool_types))
    nr_layers = len(layers)


    ############################################################################
    #
    # Options to load and prepare the model output correctly
    #
    ############################################################################

    ## netCDF file to be used
    filename = '../Data/JAMES.nc' # 2019-00-22, 3 new sites

    ## number of days per discrete time step
    days_per_ts = 10 # up to 15 is possible for this data set

    ## number of time steps to start later 
    delay_in_ts = 0

    ## the number of time steps the simulation should last
    max_sim_length = (365*110)//days_per_ts - delay_in_ts
    
    ## daily 14C decay rate    
    decay_rate = np.log(2)/5568.0/365.0


    ##### end of global options #####


    ## load dataset, prepare parameter structure
    ds = xr.open_dataset(filename)
    print('Dataset %s loaded\n' % filename)

    cell_names = ['Temperate forest', 'Tropical forest', 'Boreal forest']

    parameter_values = {
        'days_per_ts': days_per_ts,
        'delay_in_ts': delay_in_ts,
        'max_sim_length': max_sim_length,
        'decay_rate': decay_rate
    }
    

    ############################################################################
    #
    # Produce data and figures for JAMES paper
    #
    ############################################################################


    print('\nComputing error table...\n')
    print_errors(
        ds, 
        cell_names, 
        pools, 
        parameter_values,
        'Output/relative_errors_%02d.txt' % days_per_ts
    )

    print('\nPlotting relative 14C errors...\n')
    cell_nr = 0
    plot_14C_rel_err(
        ds,
        cell_nr,
        cell_names[cell_nr],
        pools,
        parameter_values,
        'Output/C14_through_time_rel_err_%02d.pdf' % days_per_ts
    )

    print('\nPlotting Delta 14C depth profiles...\n')
    plot_depth_profile(
        ds, 
        cell_names, 
        pools, 
        parameter_values,
        'Output/depth_profile_%02d.pdf' % days_per_ts
    )

    print('\nPlotting Delta 14C through time...\n')
    plot_Delta_14C_through_time(
        ds, 
        cell_names, 
        pools,
        parameter_values, 
        'Output/Delta_14C_through_time_per_pools_%02d.pdf' % days_per_ts
    )
    
