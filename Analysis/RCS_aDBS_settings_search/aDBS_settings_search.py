"""
Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/05/2022
*Note, as the script is updated, the author may not be accurate anymore*

This script calculates the FFT from raw RCS time domain data.
The script takes in two inputs, the first is the path to the raw RCS time domain data, 
and the second is the path to the RCS settings used during the same recording.
The settings from the second argument will determine how the FFT is calculated.

Dependencies:   rcssim          [=] RCS simulation code written by former postdoc Tanner Dixon from Simon Little's lab
                numpy           [=] Common library and can be obtained using pip
                pandas          [=] Common library and can be obtained using pip

Inputs:         -data           [=] Flag to notify the program the next input will be the path to the data file to convert
                -gait_events    [=] Flag to notify the program the next input will be the path to the gait events used to determine the time period to stimulate
                -search_bounds  [=] Flag to notify the program the next input will be a string representation of a search space 2D matrix. For example,
                                    [(1,3),(5,7)]. In this example, there are optimizing for 2 variables (2 rows) with the first variable has a minimum value of 1 
                                    and a maximum value of 3, and the second variable has a minimum value of 5 and a maximum value of 7.
                -n_power_bands  [=] Optional input to control how many power bands it considers when trying to find the optimal band/threshold combination. 
                                    Currently, this option is not available and passing anything in will do nothing to change anything.
                -save_path      [=] Where to save the tables that are generated
"""

# import libraries and modules
import sys
import os
from ast import literal_eval
import csv
import time
# import datetime
import numpy as np
import pandas as pd
from libs.gridSearch import gridSearch as gs # grid search class
from libs.bayesOpt import bayesOpt as bo # Bayesopt class
from libs.optimizingFunctions import thresholdAccuracy as ta # function to evaluate new settings

"""
Helper functions
"""
# NEED TO REDO TO FIT IN THIS CONTEXT
def loadData(file_path: str,type: str = "data") -> pd.DataFrame:
    if type == "data":
        raw_input = pd.read_csv(file_path)
        return raw_input
    elif type == "gait_events":
        gait_events = pd.read_csv(file_path)
        return gait_events
    else:
        sys.exit('Data type flag not supported.')



def BOItrDF(X_samples,Y_samples,search_bounds,search_bounds_resolution,n_initial_samples: int = 0,gp_params: dict = None,itr_spacing: list = None, itr_reps: list = None):
    # Get number of samples
    n_samples = X_samples.shape[0]
    
    # Get number of parameters
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

    # Get iteration numbers to calculate the Gaussian process and acquisition function
    n_itr_indexes = np.zeros(1000,dtype=int)   # Allocate memory
    count = 0
    idx_value = 0
    if itr_spacing is not None:
        for i in range(len(itr_spacing)):
            if count == 0:
                idx_value = n_initial_samples + itr_spacing[i] - 1
            else:
                idx_value = n_itr_indexes[count-1] + itr_spacing[i]
            spacing = np.arange(idx_value,idx_value+itr_spacing[i]*itr_reps[i],itr_spacing[i])
            n_itr_indexes[count:count + itr_reps[i]] = spacing
            count += itr_reps[i]
    else:
        n_itr_indexes = np.arange(n_initial_samples-1,n_samples,1) # subtract one so you get the initial random sampling of the space
        count = n_samples-n_initial_samples + 1 # add one so you get the correct number of iterations + the initial random sampling

    n_itr_indexes = n_itr_indexes[0:count]

    # Allocate memory for the Gaussian process and acquisition function iterations to calculate
    # GP_itr_values = np.zeros((count,param_combos.shape[0]))
    # AF_itr_values = np.zeros((count,param_combos.shape[0]))

    # Allocate memory for a matrix output of all the iterations to plot
    GP_mat = np.zeros((count,param_combos.shape[0]+1))
    AF_mat = np.zeros((count,param_combos.shape[0]+1))

    # Create GP model
    if gp_params is not None:
        gp_model = bo.createGP(gp_params = gp_params)
    else:
        gp_model = bo.createGP()

    # Loop through the iterations.
    for i in range(0,count):
        # Timer to see how long it takes to recreate each Gaussian process and acquisition function surface
        print(f"Recreating GP and AF of iteration {n_itr_indexes[0] + 1}...", end = " ")
        recreate_tic = time.perf_counter()

        # Fit the model to current sampled points
        curr_Y = Y_samples[0:n_itr_indexes[i]]
        gp_model.fit(X = X_samples[0:n_itr_indexes[i],:], y = curr_Y)

        # Predict the values of the parameter space
        mu,std = gp_model.predict(param_combos, return_std=True)

        # Calculate acquisition funtion values
        af_vals = bo.expectedImprovement(param_combos,curr_Y,gp_model,n_param,min_fun=False)

        # Save outputs to a matrix
        GP_mat[i,0] = n_itr_indexes[i]
        GP_mat[i,1:param_combos.shape[0]+1] = mu
        AF_mat[i,0] = n_itr_indexes[i]
        AF_mat[i,1:param_combos.shape[0]+1] = af_vals

        # write out how long this recreation took
        recreate_toc = time.perf_counter()
        print(f"took {recreate_toc-recreate_tic} seconds")

    # # Add in iteration number to the beginning
    # GP_mat = np.hstack((n_itr_indexes.reshape(len(n_itr_indexes),1),GP_itr_values))
    # AF_mat = np.hstack((n_itr_indexes.reshape(len(n_itr_indexes),1),AF_itr_values))

    return param_combos, GP_mat, AF_mat



def save_multiple_csv(data,n_entries,base_path: str = os.getcwd(),header_names: np.array = None):
    file_index = 0
    data_index = 0
    end_of_data = False
    while not end_of_data:
        file_name = "%s_%d.csv" % (base_path,file_index)

        # Save current chunk
        np.savetxt(file_name,data[data_index:data_index+n_entries,:],delimiter = ",")

        if data_index+n_entries >= data.shape[0]:
            end_of_data = True

        # Update the data index number for the next file
        data_index += n_entries

        # Update file number counter
        file_index+=1




# main function to load in data and run the Bayesian optimization function
def main(n_power_bands: int = 1,save_path: str = os.getcwd()):
    # timing measurement
    tic = time.perf_counter()

    # grab command line arguments
    args = sys.argv[1:]

    # load in data and grab RCS amplitude gains
    for i in range(0,len(args),2):
        if args[i] == "-data":
            rcs_power_data = loadData(args[i+1],type="data")
            data_file_name = os.path.basename(args[i+1])        # grabs the current filename to be used as a filename.
        elif args[i] == "-gait_events":
            gait_events = loadData(args[i+1],type="gait_events")
        elif args[i] == "-n_power_bands":
            n_power_bands = args[i+1]       # Set number of power bands to consider for a control signal. Leave at 1 for now.
        elif args[i] == "-search_bounds":
            search_bounds_str = args[i+1]
            search_bounds = np.array(literal_eval(search_bounds_str))
        elif args[i] == "-search_bounds_resolution":
            search_bounds_resolution_str = args[i+1]
            search_bounds_resolution = np.array(literal_eval(search_bounds_resolution_str))
        elif args[i] == "-save_path":
            save_path = args[i+1]           # Set the save path if given, otherwise will use the current working directory of the code

    # generate save name
    # curr_time = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_results.csv"
    save_name_GP_itr = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_GP_surf"
    save_name_AF_itr = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_AF_surf"
    save_name_param_combos = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_param_combos.csv"


    # Run Bayesian optimization
    # event_strings = []    # TODO add in options to specify gait phase
    gp_options = {'nu':0.5,'length_scale':np.array([1e-2,2.5]),'length_scale_bounds':"fixed"}
    # n_itrs = 55269  # 10% of the total parameter space
    n_itrs = 100
    # out_DF,raw_X,raw_Y = bo.runBO(ta.calcThresholdAccuracyDST,search_bounds,search_bounds_resolution,rcs_power_data,gait_events,n_itr=n_itrs,gp_params=gp_options,n_random_initial_samples=6141) # n_random_initial_sample is 1% of total parameter space
    out_DF,raw_X,raw_Y = bo.runBO(ta.calcThresholdAccuracyDST,search_bounds,search_bounds_resolution,rcs_power_data,gait_events,n_itr=n_itrs,gp_params=gp_options,n_random_initial_samples=100) # n_random_initial_sample is 1% of total parameter space
    out_DF.to_csv(save_name_results,index = False)

    # Generate dataframes of the Gaussian process and acquisition function so I can plot in MATLAB
    param_combos,GP_itr,AF_itr = BOItrDF(raw_X,raw_Y,search_bounds,search_bounds_resolution,n_initial_samples=out_DF.shape[0]-n_itrs,gp_params = gp_options,itr_spacing = [10000,1000,100,10,1], itr_reps = [4,7,18,28,48])
    # param_combos,GP_itr,AF_itr = BOItrDF(raw_X,raw_Y,search_bounds,search_bounds_resolution,n_initial_samples=out_DF.shape[0]-n_itrs,gp_params = gp_options)

    # Save the iteration outputs and parameter combos that were evaluated
    # np.savetxt(save_name_param_combos,param_combos,delimiter = ",")
    save_multiple_csv(GP_itr,n_entries = 50,base_path = save_name_GP_itr)
    save_multiple_csv(AF_itr,n_entries = 50,base_path = save_name_AF_itr)

    # Run grid search
    # out_DF= gs.runGS(ta.calcThresholdAccuracyDST,search_bounds,search_bounds_resolution,rcs_power_data,gait_events)
    # out_DF.to_csv(save_name_results,index = False)

    # End run time
    toc = time.perf_counter()

    print(f"Code ran in {toc-tic:0.4f} seconds")


"""
Entry point into the script
"""
if __name__=="__main__":
    main()