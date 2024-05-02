"""
Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/14/2023
*Note, as the script is updated, the author may not be accurate anymore*

Bayesian Optimization iteration code to determine the best weight and threshold value to classify
fast aDBS states


Dependencies:   numpy       [=] Common library and can be obtained using pip
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
def createParamSpace(parameter_ranges: list,parameter_resolutions: list):
    # Get number of parameters
    n_param = len(parameter_ranges)
    param_space = [None]*n_param

    # Create spaces
    for i in range(0,n_param):
        param_space[i] = np.arange(parameter_ranges[i][0],parameter_ranges[i][1],parameter_resolutions[i])

    # Create a [?,n_param] array of evaluation points using the parameter space to be evaluated by the acquisition function
    if n_param == 1:
        param_combos = np.array(param_space).reshape((-1,1))
    elif n_param == 2:
        param_combos = np.array([[p1,p2] for p1 in param_space[0] for p2 in param_space[1]])
    elif n_param == 3:
        param_combos = np.array([[p1,p2,p3] for p1 in param_space[0] for p2 in param_space[1] for p3 in param_space[2]])

    return param_combos


def runBF(cost_function, RCS_data: pd.DataFrame, bf_options:dict):
    """
    This function runs the Bayesian optimization algorithm. The data follows this flow:
        (1) Create a Gaussian Process (GP) with the tested parameters optimizing for,
        (2) Minimize the acquisition function to get the next sample point. (A sign change is added such that if maximizing,
            the minimization still works),
        (3) Evaluate the next sample point on the function being optimized,
        (4) Update the list of sampled parameters and function output,
        (5) Repeat until number of iterations is complete

    Inputs:     cost_function   [=] Function to determine how well the sample point performs
                RCS_data        [=] Table of LFP neural data and associated class. Rows are observations. 1 column is mandatory, "Class". 
                                    At least one column name must have the form keyX_YY_YYYYYYYY, where X can be a number from 0-3 and Y can take any
                                    numeric value. The first two YY can be a single Y if the frequency band is under 10.
                bo_options      [=] Dictionary that describes the search space and GP hyperparameters

    Outputs:
    """
    # Create parameter space to run Bayesian optimization
    param_space = createParamSpace(bf_options["parameter_space_range"],bf_options["parameter_space_resolution"])

    # Go through all possible parameter combinations
    thresholds = np.zeros(param_space.shape[0])     # Hold all threshold values used
    accuracies = np.zeros(param_space.shape[0])     # Hold accuracies of all parameter combos
    gp_accuracies = np.zeros(param_space.shape[0])  # Hold gait phase specific accuracies of all parameter combos
    for i in range(param_space.shape[0]):
        # details = bf_options["hemisphere"] + ": "
        # for x in bf_options["freq_bands"]:
        #     details += x + " "
        # print(f"{details} : {param_space[i]} ({i+1}/{param_space.shape[0]})")

        # itr_tic = time.perf_counter()
        threshold = np.round(np.quantile(np.sum(RCS_data[:,1:],axis = 1),param_space[i][0]),2)
        
        power = np.zeros(RCS_data.shape[0])
        for j in range(1,param_space.shape[1]):
            power += (RCS_data[:,j]*param_space[i][j])
        
        loss,accuracy,gp_accuracy = cost_function(threshold,power,RCS_data[:,0])
        accuracies[i] = accuracy
        thresholds[i] = threshold
        gp_accuracies[i] = gp_accuracy

        # itr_toc = time.perf_counter()
        # print(f"took {itr_toc-itr_tic} seconds")

    # Create output
    column_names = []
    column_names.append("Threshold")
    for i in range(param_space.shape[1]-1):
        column_names.append(f"Weight{i+1}")
    column_names.append("Accuracy")
    column_names.append("GaitPhaseAccuracy")

    result_table = pd.DataFrame(np.hstack((thresholds.reshape(-1,1),param_space[:,1:],accuracies.reshape(-1,1),gp_accuracies.reshape(-1,1))),columns = column_names)
    # hemisphere_vec = [bf_options["hemisphere"]]*len(thresholds)
    result_table.insert(0,"Hemisphere",bf_options["hemisphere"],True)
    for i in range(len(bf_options["freq_bands"])-1,-1,-1):
        # freq_band_vec = [bf_options["freq_bands"][i]]*len(thresholds)
        result_table.insert(1,f"FreqBand{i+1}",bf_options["freq_bands"][i])

    # Filter out values where Accuracy and GaitPhaseAccuracy is not at least 0.5
    filt_ind = np.where((result_table.Accuracy>=0.5) & (result_table.GaitPhaseAccuracy>=0.5))
    result_table = result_table.iloc[filt_ind]


    # Return sampled values and outputs. 
    # Possibly add in the EI and GP output for each iteration as well.
    return result_table