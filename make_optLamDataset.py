#################################
# Importing Libraries:
#################################
from scipy.optimize import curve_fit
from scipy.optimize import golden
from scipy.optimize import fminbound
import addcopyfighandler
import statistics
import math
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
TE_array = np.arange(8, 512, 8) #ms units


#### Model parameters
c1 = 0.5
c2 = 0.5
T21 = 45
T22 = 200
T11 = 600
T12 = 1200

true_params = np.array([T11, T12, c1, c2, T21, T22])

#### Nullpoint Values
TI1star = np.log(2)*T11
TI2star = np.log(2)*T12

#### Process parameters
SNR_fixed = 100
num_iters = 10
rand_seed = 10
num_starts = 10




def calc_lambda_error(TI, tdata, function, starts = num_starts):

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