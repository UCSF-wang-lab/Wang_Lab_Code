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

# import libraries and modules
import numpy as np
import pandas as pd

# import timer for code timing
import time



"""
Module functions
"""
def runGS(func_2_optimize, search_bounds: np.array, search_bounds_resolution: np.array, RCS_data: pd.DataFrame, event_timings:pd.DataFrame, event_strings: list = None):
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
    # Grab number of parameters of the search space
    n_param = search_bounds.shape[0]
    
    # Allocate memory for a list
    param_space = [None] * n_param

    # Calculate the parameter space as two separate np arrays and save to the list
    for i in range(0,n_param):
        param_space[i] = np.arange(search_bounds[i,0],search_bounds[i,1],search_bounds_resolution[i]).tolist()

    # Create a [?,n_param] array of evaluation points using the parameter space to be evaluated by the acquisition function
    if n_param == 2:
        param_combos = np.array([[p1,p2] for p1 in param_space[0] for p2 in param_space[1]])
    elif n_param == 3:
        param_combos = np.array([[p1,p2,p3] for p1 in param_space[0] for p2 in param_space[1] for p3 in param_space[2]])

    # Go through all parameter combos and calculate the accuracy
    combo_accuracy = np.zeros((param_combos.shape[0],1))
    RCS_time = RCS_data.time
    for i in range(param_combos.shape[0]):
        print(f"Calculating accuracy of parameter combo {i+1}/{param_combos.shape[0]}", end = " ")
        combo_tic = time.perf_counter()

        pb_data = RCS_data.iloc[:,param_combos[i,0].astype(int)].to_numpy()
        threshold_val = param_combos[i,1]
        combo_accuracy[i] = func_2_optimize(RCS_time,pb_data,event_timings,threshold_val)

        combo_toc = time.perf_counter()
        print(f"took {combo_toc-combo_tic} seconds")

    # Create output
    result_table = pd.DataFrame(np.hstack((param_combos,combo_accuracy)),columns = ["Frequency Band Ind","Threshold","Accuracy"])

    # Return sampled values and outputs. 
    # Possibly add in the EI and GP output for each iteration as well.
    return result_table, param_combos, combo_accuracy