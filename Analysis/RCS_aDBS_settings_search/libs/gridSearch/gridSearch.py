"""
Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/05/2022
*Note, as the script is updated, the author may not be accurate anymore*

This script calculates the FFT from raw RCS time domain data.
The script takes in two inputs, the first is the path to the raw RCS time domain data, 
and the second is the path to the RCS settings used during the same recording.
The settings from the second argument will determine how the FFT is calculated.

Dependencies:   rcssim      [=] RCS simulation code written by former postdoc Tanner Dixon from Simon Little's lab
                numpy       [=] Common library and can be obtained using pip
                pandas      [=] Common library and can be obtained using pip

Inputs:         -data       [=] Flag to notify the program the next input will be the path to the data file to convert
                -settings   [=] Flag to notify the program the next input will be the path to the settings file of the recording
                -sim_type   [=] Optional input to determine what to calculate. Option 1 (default): "fft" creates and saves a table of fft values for each window. 
                                Option 2: "pb" appends power band values to the raw RCS data table and settings that were read in.
                -save_path  [=] Where to save the tables that are generated
"""

# Custom module, not sure if needed
# from libs.optimizingFunctions import thresholdAccuracy as ta # function to evaluate new settings
from libs.bayesOpt.bayesOpt import createParamSpace as cps

# import libraries and modules
import numpy as np
import pandas as pd

# import timer for code timing
import time



"""
Module functions
"""
def runGS(func_2_optimize, RCS_data: pd.DataFrame, event_timings:pd.DataFrame,gs_options:dict):
    """
    This function runs the Bayesian optimization algorithm. The data follows this flow:
        (1) Create a GP with the tested parameters optimizing for,
        (2) Minimize the acquisition function to get the next sample point. (A sign change is added such that if maximizing,
            the minimization still works),
        (3) Evaluate the next sample point on the function being optimized,
        (4) Update the list of sampled parameters and function output,
        (5) Repeat until number of iterations is complete

    Inputs:     n_itr               [=] The number of iterations the Bayesian optimization should do (Optional)
                gp_params           [=] Parameters for the Gaussian process regressor object (Optional)
                initial_samples     [=] Initial evaluations of the function to optimize to cut down on repeats. 
                                        nxd where n is the number of samples evaluated and d is the dimensions to optimize for. 
                initial_outcomes    [=] The accuracy value of the optimizing function evaluated at the initial samples above. 

    Outputs:
    """
    """
    This function runs the Bayesian optimization algorithm. The data follows this flow:
        (1) Create a GP with the tested parameters optimizing for,
        (2) Minimize the acquisition function to get the next sample point. (A sign change is added such that if maximizing,
            the minimization still works),
        (3) Evaluate the next sample point on the function being optimized,
        (4) Update the list of sampled parameters and function output,
        (5) Repeat until number of iterations is complete

    Inputs:     n_itr               [=] The number of iterations the Bayesian optimization should do (Optional)
                gp_params           [=] Parameters for the Gaussian process regressor object (Optional)
                initial_samples     [=] Initial evaluations of the function to optimize to cut down on repeats. 
                                        nxd where n is the number of samples evaluated and d is the dimensions to optimize for. 
                initial_outcomes    [=] The accuracy value of the optimizing function evaluated at the initial samples above. 

    Outputs:
    """
    # Go through variable inputs
    for arg, value in gs_options.items():
        match arg:
            case "param_combos":
                param_combos = value
            case "space_type":
                space_type = value
            case "event_strings":
                event_strings = value   # FUTURE CAPABILITY, RIGHT NOW ONLY DOUBLE SUPPORT TIME

    # Set default values for variables needed to run Bayesian optimization
    # Grab number of parameters of the search space
    if "param_combos" not in locals():
        param_combos = cps(RCS_data,type=space_type)
    
    # Allocate memory to hold sampled parameters and their outcomes.
    Y = np.zeros((param_combos.shape[0],1))

    # Run through all parameter combos
    RCS_time = RCS_data.time
    for i in range(param_combos.shape[0]):
        # Timer to see how long it takes to run each iteration
        if "key" in gs_options.keys():
            channel = gs_options["key"]
            print(f"Running key{channel} parameter combo {i+1}/{param_combos.shape[0]}")
        else:
            print(f"Running parameter combo {i}/{param_combos.shape[0]}...", end = " ")
            param_tic = time.perf_counter()

        # Evaluate the next sample point
        # Can use a different type of rounding because this is not a numpy array...
        pb_data = RCS_data.iloc[:,param_combos[i,0].astype(int)].to_numpy()
        threshold_val = param_combos[i,1]
        Y[i] = func_2_optimize(RCS_time,pb_data,event_timings,threshold_val)

        # write out how long this iteration took
        if "key" not in gs_options.keys():
            param_toc = time.perf_counter()
            print(f"took {param_toc-param_tic} seconds")

    # Create output
    result_table = pd.DataFrame(np.hstack((param_combos,Y)),columns = ["Frequency Band Ind","Threshold","Accuracy"])

    # Return sampled values and outputs. 
    # Possibly add in the EI and GP output for each iteration as well.
    return result_table