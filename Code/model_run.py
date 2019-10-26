import numpy as np

from numpy.linalg import pinv
from scipy.linalg import expm, inv
from scipy.optimize import root
from sympy import Function, Matrix, symbols, zeros
from tqdm import tqdm

from CompartmentalSystems.smooth_model_run import SmoothModelRun
from CompartmentalSystems.smooth_reservoir_model import SmoothReservoirModel


################################################################################


class Error(Exception):
    """Generic error occurring in this module."""
    pass


################################################################################


class ModelRun(object):
    def __init__(self, model, parameter_set, 
                        start_values, times, disc_pts, func_set=None):

        ## add discrete points to times
        times = np.append(np.asarray(times), disc_pts)
        times = np.sort(np.unique(times))
        #times = list(sorted(set(times)))

        ## add start and end time to discrete points
        discp_pts = np.asarray(disc_pts)
        disc_pts = np.append(disc_pts, [times[0], times[-1]])
        disc_pts = np.sort(np.unique(disc_pts))

        smrs = []
        k_start_values = np.asarray(start_values)
        print('creating pwc ModelRun')
        for k in tqdm(range(len(disc_pts)-1)):
            k_t0 = disc_pts[k]
            k_t1 = disc_pts[k+1]
            k_times = times[(times >= k_t0) & (times <= k_t1)]
            k_smr = SmoothModelRun(
                model, 
                parameter_set, 
                k_start_values, 
                k_times, 
                func_set
            )
            smrs.append(k_smr)

            k_soln = k_smr.solve()
            k_start_values = k_soln[-1]

        self.smrs = smrs

        self.model = model
        #self.nr_pools = model.nr_pools
        self.parameter_set = parameter_set
        self.start_values = start_values
        self.times = times
        self.disc_pts = disc_pts
        if func_set is None:
            func_set = dict()
        self.func_set = func_set

    @classmethod
    def reconstruct_from_data(cls, time_symbol, data_times, start_values, times, xs, 
                                 Fs, rs, data_us):
        nr_pools = len(start_values)
        strlen=len(str(nr_pools))
        pool_str = lambda i: ("{:0"+str(strlen)+"d}").format(i)
        par_set = {}
        func_set = dict()
    
        us = cls.reconstruct_us(data_times, data_us)
        print('creating us')
        u_funcs = cls.u_pwc(data_times, us)
        print('creating Bs')
        Bs = cls.reconstruct_Bs(data_times, xs, Fs, rs, us)
        B_funcs = cls.B_pwc(data_times, Bs)
    
        srm_generic = cls.create_srm_generic(time_symbol, Bs, us)
        time_symbol = srm_generic.time_symbol
    
        def func_maker_u(pool):
            def func(s):
                return u_funcs[pool](s)
            return func
    
        v = srm_generic.external_inputs
        for i in range(nr_pools):
            if v[i].is_Function:
                func_set['u_'+pool_str(i)+'('+time_symbol.name+')'] = func_maker_u(i)
    
    
        def func_maker_B(p_to, p_from):
            def func(s):
                return B_funcs[(p_to, p_from)](s)
            return func
    
        M = srm_generic.compartmental_matrix
        for j in range(nr_pools):
            for i in range(nr_pools):
                if M[i,j].is_Function:
                    func_set['b_'+pool_str(i)+pool_str(j)+'('+time_symbol.name+')'] = \
                        func_maker_B(i,j)
    
        #### why not here a smr directly? ####
        mr_pwc = cls(
            srm_generic, 
            par_set, 
            start_values, 
            times,
            data_times,
            func_set)

        mr_pwc.us = us
        mr_pwc.Bs = Bs

        return mr_pwc

    @classmethod
    def reconstruct_B_surrogate(cls, dt, x0, F, r, u, B0):
        ## reconstruct a B that meets F and r possibly well,
        ## B0 is some initial guess for the optimization
        nr_pools = len(x0)
    
        ## integrate x by trapezoidal rule
        def x_trapz(tr_times, B):
            def x(tau):
                M = expm(tau*B)
                #x = M @ x0 + inv(B) @ (-np.identity(nr_pools)+M) @ u
                x = M @ x0 + pinv(B) @ (-np.identity(nr_pools)+M) @ u
                return x
    
            xs = np.array(list(map(x, tr_times)))
            return np.trapz(xs, tr_times, axis=0)
    
        ## convert parameter vector to compartmental matrix
        def pars_to_matrix(pars):
            B = np.zeros((nr_pools**2))
            B[int_indices.tolist()] = pars[:len(int_indices)].copy()
            B = B.reshape((nr_pools,nr_pools))
            d = B.sum(0)
            d[(ext_indices-nr_pools**2).tolist()] += pars[len(int_indices):]
            B[np.diag_indices(nr_pools)] = -d
    
            return B
    
        ## function to minimize difference vector of
        ## internal fluxes F and outfluxes r
        def g_tr(pars):
            B = pars_to_matrix(pars)
            tr_times = np.linspace(0, dt, 11)
            int_x = x_trapz(tr_times, B)
    
            res1_int = F.reshape((nr_pools**2,))[int_indices.tolist()]
            res2_int = pars[:len(int_indices)] * int_x[(int_indices % nr_pools).tolist()]
            res_int = np.abs(res1_int-res2_int)
    
            res1_ext = r[(ext_indices-nr_pools**2).tolist()] 
            res2_ext = pars[len(int_indices):] * int_x[(ext_indices-nr_pools**2).tolist()]
            res_ext = np.abs(res1_ext-res2_ext)
    
            res = np.append(res_int, res_ext)
#            print(pars, res1_int, res2_int, res1_ext, res2_ext, res)
    
            return res
    
        ## wrapper around g_tr that keeps parameter within boundaries
        def constrainedFunction(x, f, lower, upper, minIncr=0.001):
            x = np.asarray(x)
            lower = np.asarray(lower)
            upper = np.asarray(upper)
    
            xBorder = np.where(x<lower, lower, x)
            xBorder = np.where(x>upper, upper, xBorder)
            fBorder = f(xBorder)
            distFromBorder = (np.sum(np.where(x<lower, lower-x, 0.))
                          +np.sum(np.where(x>upper, x-upper, 0.)))
            return (fBorder + (fBorder
                           +np.where(fBorder>0, minIncr, -minIncr))*distFromBorder) 
    
        lbounds = [0]*(nr_pools**2 + nr_pools)
        for i in range(nr_pools):
            lbounds[i*nr_pools+i] = -10
    
        ubounds = [10]*(nr_pools**2 + nr_pools)
        for i in range(nr_pools):
            ubounds[i*nr_pools+i] = 0
    
        par_indices = []
        for i in range(nr_pools):
            for j in range(nr_pools):
                if (F[i,j] > 0):
                    par_indices.append(i*nr_pools+j)
        for i in range(nr_pools):
            if r[i] > 0:
                par_indices.append(nr_pools**2+i)
    
        par_indices = np.array(par_indices)
        int_indices = par_indices[par_indices<nr_pools**2]
        ext_indices = par_indices[par_indices>=nr_pools**2]
    
        lbounds = np.array(lbounds)[par_indices.tolist()]
        ubounds = np.array(ubounds)[par_indices.tolist()]
        A0 = np.append(B0.reshape((nr_pools**2,)), -B0.sum(0))
        pars0 = A0[par_indices.tolist()]
    
        y = root(
            constrainedFunction, 
            x0=pars0,
            args=(g_tr, lbounds, ubounds)
        )
    
        B = pars_to_matrix(y.x)
    
        return B

    @classmethod
    def reconstruct_Bs(cls, data_times, xs, Fs, rs, us):
        nr_pools = len(xs[0])
    
        def guess_B0(dt, x_approx, F, r):
            nr_pools = len(x_approx)
            B = np.identity(nr_pools)
        
            # construct off-diagonals
            for j in range(nr_pools):
                if x_approx[j] != 0:
                    B[:,j] = F[:,j] / x_approx[j] / dt
                else:
                    B[:,j] = 0
        
            # construct diagonals
            for j in range(nr_pools):
                if x_approx[j] != 0:
                    B[j,j] = - (sum(B[:,j]) - B[j,j] + r[j] / x_approx[j] / dt)
                else:
                    B[j,j] = -1
        
            return B
    
        x = xs[0]
        Bs = np.zeros((len(data_times)-1, nr_pools, nr_pools))
        for k in tqdm(range(len(data_times)-1)):
            dt = data_times[k+1] - data_times[k]
            #dt = dt.item().days ##### to be removed or made safe
            B0 = guess_B0(dt, (xs[k]+xs[k+1])/2, Fs[k], rs[k])
            B = cls.reconstruct_B_surrogate(dt, x, Fs[k], rs[k], us[k], B0)
            M = expm(dt*B)
            x = M @ x + inv(B) @ (-np.identity(nr_pools)+M) @ us[k]
    
            Bs[k,:,:] = B
    
        return Bs

    @classmethod
    def B_pwc(cls, data_times, Bs):
        def func_maker(i,j):
            def func(t):
                index = np.where(data_times<=t)[0][-1]
                index = min(index, Bs.shape[0]-1)
                return Bs[index,i,j]
            return func    

        nr_pools = Bs[0].shape[0]
        B_funcs = dict()
        for j in range(nr_pools):
            for i in range(nr_pools):
                B_funcs[(i,j)] = func_maker(i,j)

        return B_funcs

    @classmethod
    def reconstruct_u(cls, dt, data_u):
        return data_u / dt

    @classmethod
    def reconstruct_us(cls, data_times, data_us):
        us = np.zeros_like(data_us)
        for k in range(len(data_times)-1):
            dt = data_times[k+1] - data_times[k]
            #dt = dt.item().days ##### to be removed or made safe
            us[k,:] = cls.reconstruct_u(dt, data_us[k])

        return us

    @classmethod
    def u_pwc(cls, data_times, us):
        def func_maker(i):
            def func(t):
                index = np.where(data_times<=t)[0][-1]
                index = min(index, us.shape[0]-1)
                return us[index,i]
            return func    
    
        nr_pools = us[0].shape[0]
        u_funcs = []
        for i in range(nr_pools):
            u_funcs.append(func_maker(i))
    
        return u_funcs

    @classmethod
    def create_srm_generic(cls, time_symbol, Bs, us):
        nr_pools = Bs.shape[1]
        strlen = len(str(nr_pools))
        pool_str = lambda i: ("{:0"+str(strlen)+"d}").format(i)
    
        state_vector_generic = Matrix(nr_pools, 1, 
            [symbols('x_'+pool_str(i)) for i in range(nr_pools)])
    
        def is_constant_Bs(i,j):
            c = Bs[0,i,j]
            diff = Bs[:,i,j]-c
    
            if len(diff[diff == 0]) == len(Bs):
                res = True
            else:
                res = False
    
            return res
    
        B_generic = zeros(nr_pools, nr_pools)
        for j in range(nr_pools):
            for i in range(nr_pools):
                if not is_constant_Bs(i,j):
                    B_generic[i,j] = Function('b_'+pool_str(i)+pool_str(j))(time_symbol)
                else:
                    B_generic[i,j] = Bs[0,i,j]
    
        def is_constant_us(i):
            c = us[0,i]
            diff = us[:,i]-c
    
            if len(diff[diff == 0]) == len(us):
                res = True
            else:
                res = False
    
            return res
    
        u_generic = zeros(nr_pools, 1)
        for i in range(nr_pools):
            if not is_constant_us(i):
                u_generic[i] = Function('u_'+pool_str(i))(time_symbol)
            else:
                u_generic[i] = us[0,i]
    
        srm_generic = SmoothReservoirModel.from_B_u(
            state_vector_generic,
            time_symbol,
            B_generic,
            u_generic)
    
        return srm_generic

    def solve(self):
        soln = self.smrs[0].solve()
        for k in range(1, len(self.disc_pts)-1):
            soln = soln[:-1]
            k_soln = self.smrs[k].solve()
            soln = np.concatenate((soln, self.smrs[k].solve()), axis=0)

        return soln










