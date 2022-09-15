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
c1_set = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
c2_set = 1 - c1_set
T21_set = np.array([10, 20, 30, 40, 50, 60])
T22_set = np.array([70, 85, 100, 120, 150, 200])

#### Process Parameters
SNR_fixed  = 100    #SNR value used to generate noise realizations
num_iters  = 10     #Number of noise realizations for a single parameter sample
rand_seed  = 10     #Seed for rng if we want rng to be fixed
num_starts = 10     #Number of random starts used in NLLS
low_power = -7      #Low power of lambda used
high_power = 3      #High power of lambda used


#################################
# Signal Related Functions:
#################################

def S_biX_4p(TE, d1, d2, T21, T22):
    exp1 = d1*np.exp(-TE/T21)
    exp2 = d2*np.exp(-TE/T22)
    return exp1 + exp2

def add_noise(data,SNR = SNR_fixed):
    sigma = 1/SNR #np.max(np.abs(data))/SNR
    noise = np.random.normal(0,sigma,data.shape)
    noised_data = data + noise
    return noised_data


#################################
# Process Functions:
#################################

def generate_noiseySig_lambda_pairs(params, tdat = TE_array, iterations = num_iters):

    signal = S_biX_4p(tdat, *params)

    noised_array = []
    for iter in range(iterations):
        noised_sig = add_noise(signal)
        noised_array.append(noised_sig)

        lamb_opt = fminbound(minLambda_objFunc, )





def minLambda_objFunc(lamb, tdata, function, starts = num_starts):

    np.linalg.norm(agg_arr*(est_high-p_true))


    data = generate_noised_sigs(TI, tdata)

    iteration_RSS_values = np.zeros((data.shape[0],1))
    
    for i in range(data.shape[0]):
        noisey_signal = data[i,:]

        recreated_curves = np.zeros((starts,np.size(tdata)))
        recreated_curves_RSS = np.zeros((starts,1))
        recreated_popt = []                     #general to accomodate for different parameter sizes
        for start in range(starts):
            recreated_popt.append(estimate_NLLS(noisey_signal, tdata, function))
            recreated_curves[start,:] = function(tdata, *recreated_popt[start])
            recreated_curves_RSS[start] = np.sum((noisey_signal - recreated_curves[start,:])**2)

        iteration_RSS_values[i] = np.min(recreated_curves_RSS)

    return np.mean(iteration_RSS_values)