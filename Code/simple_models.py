## Reproduce the simple model plots of JAMES manuscript
##
## Title:        Mathematical reconstruction of land carbon models from their
##               numerical output: computing soil radiocarbon from 12C dynamics
##
## Authors:      Holger Metzler, Qing Zhu, William Riley, Alison Hoyt,
##               Markus Müller, Carlos A. Sierra 
##
## Code:         Holger Metzler (10/2019)
##
## Requirements: CompartmentalSystems 
##                  (https://github.com/MPIBGC-TEE/CompartmentalSystems)
##               matplotlib
##               numpy
##               sympy
##               tqdm
##               xarray
##               scipy
##               model_run.py (included)
##
## Input:       
##
## Output:       To subfolder 'Output"


import matplotlib.pyplot as plt
import numpy as np

from sympy import Function, Matrix, sin, symbols

from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel
from CompartmentalSystems.smooth_model_run import SmoothModelRun

from model_run import ModelRun


################################################################################


def plot_approximations_one_axis(smr, nr_data_points_list, filename):
    fig = plt.figure(figsize=(30,15))
    ax = fig.add_subplot(1,1,1)
    fontsize = 40

    soln = smr.solve()
    ax.plot(smr.times, soln, c='black', label='original', lw=3)
    for nr, nr_data_points in enumerate(nr_data_points_list):
        print('Number of data points:', nr_data_points)
        srm = smr.model
        times = smr.times

        ## create fake discrete data
        data_times = np.round(np.linspace(times[0], times[-1], nr_data_points),5)
        delta_t = (times[-1]-times[0]) / (nr_data_points-1)
        print('Delta t:', delta_t)
        xs, Fs, rs, data_us = smr._fake_discretized_output(data_times)
        print('fake data created')
    
        time_symbol = symbols('t')
        mr_pwc = ModelRun.reconstruct_from_data(
            time_symbol,
            data_times, 
            start_values, 
            times, 
            xs, 
            Fs,
            rs, 
            data_us)
        print('mr_pwc created')
        soln_pwc = mr_pwc.solve()

        bias = np.abs(soln.sum(1)-soln_pwc.sum(1))/soln.sum(1) * 100
        print('mean annual bias, %2.2f per cent:' % delta_t, np.mean(bias)) 

        colors = ['red', 'blue', 'olive', 'black']
        ax.plot(smr.times, soln_pwc.sum(1), 
                c=colors[nr], 
                dashes=[3,4],
                label='CTA, $\Delta t=%2d$ yr' % delta_t, lw=3)

        print()        

    ax.set_xlim([times[0], times[-1]])
    ax.set_ylim([max(0, ax.get_ylim()[0]), ax.get_ylim()[1]])

    for item in [ax.title, ax.xaxis.label, ax.yaxis.label]:
        item.set_fontsize(fontsize)

    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)


    ax.legend(fontsize=fontsize)
    ax.set_xlabel('time (yr)')
    ax.set_ylabel('C')
    
    fig.savefig(filename, tight_layout=True)
    plt.close(fig)


def plot_approximations_n_axes(smr, nr_data_points_list, figsize, filename):
    nr_pools = smr.nr_pools
    nr_axes = len(nr_data_points_list)
    fig = plt.figure(figsize=figsize)
    fontsize = 14

    soln = smr.solve()
    for nr, nr_data_points in enumerate(nr_data_points_list):
        print('Number of data points:', nr_data_points)
        srm = smr.model
        times = smr.times

        ## create fake discrete data
        data_times = np.round(np.linspace(times[0], times[-1], nr_data_points),5)
        xs, Fs, rs, data_us = smr._fake_discretized_output(data_times)
#        print('xs', xs)
#        print('Fs', Fs)
#        print('rs', rs)
#        print('data_us', data_us)

        print('fake data created')
    
        time_symbol = symbols('t')
        mr_pwc = ModelRun.reconstruct_from_data(
            time_symbol,
            data_times, 
            start_values, 
            times, 
            xs, 
            Fs,
            rs, 
            data_us)
        print('mr_pwc created')

        soln_pwc = mr_pwc.solve()

        ## plot solution

        ax = fig.add_subplot(nr_axes, 1, nr+1)
        ax.set_title(r'$n=%d$ data points' % nr_data_points, 
                     pad = 20,
                    )

        if nr_pools == 1:
            colors = ['red', 'blue', 'olive', 'black']
            ax.plot(smr.times, soln.sum(1), c='black', label='original', lw=3)
            ax.plot(
                smr.times, 
                soln_pwc.sum(1), 
                dashes=[3,4], 
                c='red', 
                label='approximation', 
            lw=3)

            ax.set_xlim([times[0], times[-1]])
            ax.set_ylim([max(0, ax.get_ylim()[0]), ax.get_ylim()[1]])
            ax.vlines(data_times, ax.get_ylim()[0], ax.get_ylim()[1], alpha=1.0, lw=0.5)

            if nr == 0:
                ax.legend(fontsize=12)
            if nr+1 == nr_axes:
                ax.set_xlabel('time (yr)')
            ax.set_ylabel('C')
        else:
            for pool in range(nr_pools):
                sv = srm.state_vector[pool].name
                colors = ['blue', 'olive', 'red']
                ax.plot(smr.times, soln[:,pool], c=colors[pool], label='$'+sv+'$', lw=3)
                if pool == nr_pools-1:
                    l = 'approximation'
                else:
                    l = '_'
                ax.plot(smr.times, soln_pwc[:,pool], dashes=[3,4], c='red', label=l, lw=3)

            ax.set_xlim([times[0], times[-1]])
            ax.set_ylim([max(0, ax.get_ylim()[0]), ax.get_ylim()[1]])
            ax.vlines(data_times, ax.get_ylim()[0], ax.get_ylim()[1], alpha=1.0, lw=0.5)

            if nr == 0:
                ax.legend(fontsize=12)
            if nr+1 == nr_axes:
                ax.set_xlabel('time (yr)')
            ax.set_ylabel('carbon content')

        for item in [ax.title, ax.xaxis.label, ax.yaxis.label]:
            item.set_fontsize(fontsize)

        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)

    plt.subplots_adjust(hspace=0.6)
    fig.savefig(filename, tight_layout=True)
    plt.close(fig)



################################################################################


if __name__ == '__main__':
    print('One-pool model, oscillating, Delta t = 30,15,10\n')
    ## create ReservoirModel
    C = symbols('C')
    state_vector = Matrix(1, 1, [C]) 
    t = symbols('t')
    lamda = symbols('lamda')
    B = Matrix([[-lamda]])
    u = Matrix(1, 1, [3+sin(t/20)])

    srm = SmoothReservoirModel.from_B_u(state_vector, t, B, u)

    ## create ModelRun
    lamda_val = 0.04
    par_set = {'lamda': lamda_val}
    ss = (-B**(-1)*u).subs(par_set)
    start_values = np.array([40])
    times = np.linspace(1919, 2009, 901)
    
    smr = SmoothModelRun(srm, par_set, start_values, times)

    ## use three different numbers of data points
    nr_data_points_list = [4,7,10]
    plot_approximations_one_axis(
        smr, 
        nr_data_points_list,
        'Output/interpol_pwc_%1d.pdf' % smr.nr_pools
    )


    ############################################################################
    print('One-pool model, autonomous, n = 3,6\n')

    ## create ReservoirModel
    C = symbols('C')
    state_vector = Matrix(1, 1, [C]) 
    t = symbols('t')
    lamda = symbols('lamda')
    B = Matrix([[-lamda]])
    u = Matrix(1,1, [3])

    srm = SmoothReservoirModel.from_B_u(state_vector, t, B, u)

    ## create ModelRun
    lamda_val = 0.04
    par_set = {'lamda': lamda_val}
    ss = (-B**(-1)*u).subs(par_set)
    start_values = np.array([40])
    times = np.linspace(1909, 2009, 1001)
    
    smr = SmoothModelRun(srm, par_set, start_values, times)
    
    ## use two different numbers of data points
    nr_data_points_list = [3, 6]
    plot_approximations_n_axes(
        smr, 
        nr_data_points_list,
        (14,3*len(nr_data_points_list)),
        'Output/interpol_pwc_%1d_auton.pdf' % smr.nr_pools
    )


    ############################################################################
    print('Two-pool model, nonautonomous, n = 3,6,11,21\n')

    ## create ReservoirModel
    C_1, C_2 = symbols('C_1 C_2')
    state_vector = Matrix(2, 1, [C_1, C_2]) 
    t = symbols('t')
    lamda_1 = Function('lamda_1')(t)
    lamda_2 = Function('lamda_2')(t)
    B = Matrix([[-lamda_1,      0.5*lamda_2],
                [ lamda_1,         -lamda_2]])
    u = Matrix(2, 1, [0.3,0.5*0.5*sin(1/100*t)])

    srm = SmoothReservoirModel.from_B_u(state_vector, t, B, u)

    ## create ModelRun
    def rate_1(t):
        return 0.1+1/20*np.sin(1/10*t)

    def rate_2(t):
        return 0.05

    par_set = {}
    func_set = {lamda_1: rate_1, lamda_2: rate_2}
    start_values = np.array([40, 0])
    times = np.linspace(1909, 2009, 1001)
    
    smr = SmoothModelRun(srm, par_set, start_values, times, func_set)
    
    ## use four different numbers of data points
    nr_data_points_list = [3, 6, 11, 21]
    plot_approximations_n_axes(
        smr, 
        nr_data_points_list,
        (14,4*len(nr_data_points_list)),
        'Output/interpol_pwc_%1d.pdf' % smr.nr_pools
    )


