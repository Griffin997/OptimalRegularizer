#################################
# Importing Libraries:
#################################
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import golden
from scipy.optimize import fminbound
import addcopyfighandler
import time
import h5py
from tqdm import trange

import multiprocessing as mp
from multiprocessing import Pool, freeze_support
from multiprocessing import set_start_method
import concurrent.futures
from itertools import product


###Example fminbound call:
#resource for fminbound: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html
# FMB_biX_std = fminbound(lambda TI: calc_std(TI, TE_array, S_biX_4p, 2), 200, 600, xtol = 10**-1, full_output = True, disp = 3)
# print("Estimated Nullpoint via T21 std of biX results:" + convert_fullOutput(FMB_biX_std))
# print("     Error = {:.2f} :: Relative Error = {:.2f}%".format(FMB_biX_std[0] - TI1star,np.abs(TI1star-FMB_biX_std[0])/TI1star*100))


#################################
# Setting Parameters:
#################################

#### Grid Parameters
n_elements = 128
d_start = 4
d_increment = 4
TE_array = np.linspace(d_start, d_increment*n_elements, n_elements)
assert(np.diff(TE_array)[0]==d_increment)

#### Model Parameters
c1_set = np.array([0.1]) #np.array([0.1, 0.3, 0.5, 0.7, 0.9])
c2_set = 1 - c1_set
T21_set = np.array([10, 20]) #np.array([10, 20, 30, 40, 50, 60])
T22_set = np.array([70, 85]) #np.array([70, 85, 100, 120, 150, 200])

total_size = c1_set.shape[0]*c2_set.shape[0]*T21_set.shape[0]*T22_set.shape[0]

#### Process Parameters
SNR_fixed  = 100    #SNR value used to generate noise realizations
num_iters  = 100     #Number of noise realizations for a single parameter sample
rand_seed  = 10     #Seed for rng if we want rng to be fixed
num_starts = 10     #Number of random starts used in NLLS
low_power = -7      #Low power of lambda used
high_power = 3      #High power of lambda used
ob_weight = 100     #The magnitude to scale the T2 parameters by


#################################
# Signal Related Functions:
#################################

def G(t, con_1, con_2, tau_1, tau_2): 
    function = con_1*np.exp(-t/tau_1) + con_2*np.exp(-t/tau_2)
    return function

def G_tilde(lam, SA = 1, weight = ob_weight):
    #SA defines the signal amplitude, defaults to 1 for simulated data
    def Gt_lam(t, con1, con2, tau1, tau2):
        return np.append(G(t, con1, con2, tau1, tau2), [lam*con1/SA, lam*con2/SA, lam*tau1/weight, lam*tau2/weight])
    return Gt_lam

def add_noise(data,SNR = SNR_fixed):
    sigma = 1/SNR #np.max(np.abs(data))/SNR
    noise = np.random.normal(0,sigma,data.shape)
    noised_data = data + noise
    return noised_data


#################################
# Process Functions:
#################################

def get_func_bounds():

    lower_bound = (0,0,1,1)
    upper_bound = (1,1,300,1000)

    return lower_bound, upper_bound

def check_param_order(popt):
    #Reshaping of array to ensure that the parameter pairs all end up in the appropriate place - ensures that T22 > T21
    if (popt[-1] < popt[-2]): #We want by convention to make sure that T21 is <= T22
        for pi in range(np.size(popt)//2):
            p_hold = popt[2*pi]
            popt[2*pi] = popt[2*pi+1]
            popt[2*pi+1] = p_hold
    return popt

def estimate_parameters(noisey_sig, lam, init_p, tdata = TE_array):
    noisey_tilde = np.append(noisey_sig, [0,0,0,0])
    lb, ub = get_func_bounds()

    popt, _ = curve_fit(G_tilde(lam), tdata, noisey_tilde, bounds = (lb, ub), p0 = init_p, max_nfev = 1500)
    popt = check_param_order(popt)
        
    return popt

def estimate_parameters_multistart( noisey_sig, lam, tdata = TE_array, starts = num_starts):
    noisey_tilde = np.append(noisey_sig, [0,0,0,0])
    lb, ub = get_func_bounds(function)

    recreated_popt = []
    recreated_curves_RSS = np.zeros((starts,1))   
    
    for start in range(starts):
        init_p = tuple(np.add(np.subtract(ub,lb)*np.random.uniform(0,1,np.size(lb)),lb))
        popt, _ = curve_fit(G_tilde(lam), tdata, noisey_tilde, bounds = (0, ub), p0 = init_p, max_nfev = 1500)
        popt = check_param_order(popt)
        recreated_popt.append(popt)
        recreated_curves_RSS[start] = np.sum((noisey_sig - G(tdata, *recreated_popt[start]))**2)

    popt = recreated_popt[np.argmin(recreated_curves_RSS),:]
        
    return popt


def minLambda_objFunc(lamb, noisey_data, p_true, agg_arr = np.array([1, 1,1/ob_weight,1/ob_weight])):

    p_est = np.array(estimate_parameters(noisey_data, 10**lamb, p_true))
    obj_value = np.linalg.norm(agg_arr*(p_est-p_true))

    return obj_value

def generate_noiseySig_lambda_pairs_fminbound(params, tdat = TE_array, iterations = num_iters, lb = low_power, ub = high_power):

    print(f"Building {iterations} noisey signals for params:{params}")
    signal = G(tdat, *params)

    noised_array = []
    lamb_opt_array = []
    for _ in range(iterations):
        noised_sig = add_noise(signal)
        noised_array.append(noised_sig)

        lamb_pow = fminbound(lambda lamb: minLambda_objFunc(lamb, noised_sig, params), lb, ub, xtol = 10**-3)
        lamb_opt_array.append(10**lamb_pow)

    return (noised_array, lamb_opt_array)


#################################
# Code Body:
#################################

if __name__ == '__main__':
    freeze_support()
    mp.set_start_method('spawn', force=  True)

    data_array = []
    lamb_array = []
    results = []
    counter = 0
    start = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i_c1 in c1_set:
            for i_c2 in c2_set:
                for i_T21 in T21_set:
                    for i_T22 in T22_set:
                        counter += 1

                        if counter%10 == 0:
                            print(f"{counter} of {total_size} param samples complete")

                        sample_params = (i_c1, i_c2, i_T21, i_T22)

                        results.append(executor.submit(generate_noiseySig_lambda_pairs_fminbound, sample_params))

        for f in concurrent.futures.as_completed(results):
            data_array.append(f.result()[0])
            lamb_array.append(f.result()[1])

    data_full = np.array(np.concatenate(data_array))
    lamb_full = np.array(np.concatenate(lamb_array))

    stop = time.perf_counter()

    print(f"Total time = {stop - start} second(s)")
    print(f"Avg time per iteration = {(stop - start)/total_size} second(s)")

    with h5py.File('C://co//NIA//OptimalRegularizer//DataSets//first_try_16Sep22.hdf5','a') as f: 
        #Save a data set
        dset = f.create_dataset('signals',data=data_full)
        dset = f.create_dataset('lambdas',data=lamb_full)
