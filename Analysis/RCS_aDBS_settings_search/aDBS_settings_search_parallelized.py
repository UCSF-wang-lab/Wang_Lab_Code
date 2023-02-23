"""
Original Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/05/2022
*Note, as the script is updated, the author may not be accurate anymore*

This script calculates the FFT from raw RCS time domain data.
The script takes in two inputs, the first is the path to the raw RCS time domain data, 
and the second is the path to the RCS settings used during the same recording.
The settings from the second argument will determine how the FFT is calculated.

Dependencies:   rcssim          [=] RCS simulation code written by former postdoc Tanner Dixon from Simon Little's lab
                numpy           [=] Common library and can be obtained using pip
                pandas          [=] Common library and can be obtained using pip
                itertools       [=] Used to create combinations of frequency bands

Inputs:         -data                       [=] Flag to notify the program the next input will be the path to the data file to convert
                -gait_events                [=] Flag to notify the program the next input will be the path to the gait events used to determine the time period to stimulate
                -search_bounds              [=] Flag to notify the program the next input will be a string representation of a search space 2D matrix. For example,
                                                [(1,3),(5,7)]. In this example, there are optimizing for 2 variables (2 rows) with the first variable has a minimum value of 1 
                                                and a maximum value of 3, and the second variable has a minimum value of 5 and a maximum value of 7.
                -serach_bounds_resolution   [=]
                -n_power_bands              [=] Optional input to control how many power bands it considers when trying to find the optimal band/threshold combination. 
                                                Currently, this option is not available and passing anything in will do nothing to change anything.
                -save_path                  [=] Where to save the tables that are generated
"""

# import libraries and modules
import sys
import os
from ast import literal_eval
import time # code timing / debugging purposes
import numpy as np
import pandas as pd
import itertools
from libs.gridSearch import gridSearch as gs # grid search class
from libs.bayesOpt import bayesOpt as bo # Bayesopt class
from libs.optimizingFunctions import thresholdAccuracy as ta # function to evaluate new settings
import multiprocessing as mp

"""
Helper functions
"""

# loadData
#   Loads in data to use for aDBS settings search. The type of data can either be power band data, estimated using the 
#   "rcssim" package, or gait events. Currently, both types of data are read in the same way, but are kept separate
#   to provide flexibility in the future.
#
#   Inputs:
#       file_path   [=] The path to the data file
#       type        [=] The type of data being passed in, either power band or gait events.
#
#   Returns:
#       The data within the file as a panda's DataFrame object
def loadData(file_path: str,type: str = "data") -> pd.DataFrame:
    if type == "data":
        raw_input = pd.read_csv(file_path)
        return raw_input
    elif type == "gait_events":
        gait_events = pd.read_csv(file_path)
        return gait_events
    else:
        sys.exit('Data type flag not supported.')


def createPowerBands(rcs_power_data,frequency_limits):

    # Get column names, which contain the frequency value
    col_names = rcs_power_data.columns

    # Indices where the frequency value is within the range
    key0_inds = np.zeros(512,dtype=int)
    key0_count = 0
    key1_inds = np.zeros(512,dtype=int)
    key1_count = 0
    key2_inds = np.zeros(512,dtype=int)
    key2_count = 0
    key3_inds = np.zeros(512,dtype=int)
    key3_count = 0

    # Loop through column names and see if they should be include in the power band creations
    count = 0
    for i in col_names:
        if "key" in i:
            freq_val_str = i[5:len(i)]
            freq_val_str = freq_val_str.replace("_",".")
            freq_val = literal_eval(freq_val_str)
            if (freq_val >= frequency_limits[0]) & (freq_val <= frequency_limits[1]):
                match i[0:4]:
                    case "key0":
                        key0_inds[key0_count] = count
                        key0_count += 1
                    case "key1":
                        key1_inds[key1_count] = count
                        key1_count += 1
                    case "key2":
                        key2_inds[key2_count] = count
                        key2_count += 1
                    case "key3":
                        key3_inds[key3_count] = count
                        key3_count += 1
        count += 1

    # Remove extra indices that do not matter
    key0_inds = key0_inds[0:key0_count]
    key1_inds = key1_inds[0:key1_count]
    key2_inds = key2_inds[0:key2_count]
    key3_inds = key3_inds[0:key3_count]

    # Generate all combinations of indexes from each key
    key0_combos = list(itertools.combinations_with_replacement(key0_inds,2))
    key1_combos = list(itertools.combinations_with_replacement(key1_inds,2))
    key2_combos = list(itertools.combinations_with_replacement(key2_inds,2))
    key3_combos = list(itertools.combinations_with_replacement(key3_inds,2))

    # Allocate memory for power band column names
    key0_pb_col_names = [None]*len(key0_combos)
    key1_pb_col_names = [None]*len(key1_combos)
    key2_pb_col_names = [None]*len(key2_combos)
    key3_pb_col_names = [None]*len(key3_combos)

    # Allocate memory for power band matrices
    key0_pb_mat = np.zeros((rcs_power_data.shape[0],len(key0_combos)))
    key1_pb_mat = np.zeros((rcs_power_data.shape[0],len(key1_combos)))
    key2_pb_mat = np.zeros((rcs_power_data.shape[0],len(key2_combos)))
    key3_pb_mat = np.zeros((rcs_power_data.shape[0],len(key3_combos)))

    # Create power band matrices
    print(f"Creating key0 power bands ...", end = " ")
    tic = time.perf_counter()
    for i in range(len(key0_combos)):
        str = "key0_" + rcs_power_data.columns[key0_combos[i][0]][5:len(rcs_power_data.columns[key0_combos[i][0]])] + "__" + rcs_power_data.columns[key0_combos[i][1]][5:len(rcs_power_data.columns[key0_combos[i][1]])]
        key0_pb_col_names[i] = str
        for j in range(rcs_power_data.shape[0]):
            key0_pb_mat[j,i] = sum(rcs_power_data.values[j,key0_combos[i][0]:key0_combos[i][1]+1])
    key0_pb_col_names.insert(0,"time")
    key0_pb_dataframe = pd.DataFrame(np.hstack((rcs_power_data.time.to_numpy().reshape(-1,1),key0_pb_mat)),columns = key0_pb_col_names)
    toc = time.perf_counter()
    print(f"Code ran in {toc-tic:0.4f} seconds")

    print(f"Creating key1 power bands ...", end = " ")
    tic = time.perf_counter()
    for i in range(len(key1_combos)):
        str = "key1_" + rcs_power_data.columns[key1_combos[i][0]][5:len(rcs_power_data.columns[key1_combos[i][0]])] + "__" + rcs_power_data.columns[key1_combos[i][1]][5:len(rcs_power_data.columns[key1_combos[i][1]])]
        key1_pb_col_names[i] = str
        for j in range(rcs_power_data.shape[0]):
            key1_pb_mat[j,i] = sum(rcs_power_data.values[j,key1_combos[i][0]:key1_combos[i][1]+1])
    key1_pb_col_names.insert(0,"time")
    key1_pb_dataframe = pd.DataFrame(np.hstack((rcs_power_data.time.to_numpy().reshape(-1,1),key1_pb_mat)),columns = key1_pb_col_names)
    toc = time.perf_counter()
    print(f"Code ran in {toc-tic:0.4f} seconds")

    print(f"Creating key2 power bands ...", end = " ")
    tic = time.perf_counter()
    for i in range(len(key2_combos)):
        str = "key2_" + rcs_power_data.columns[key2_combos[i][0]][5:len(rcs_power_data.columns[key2_combos[i][0]])] + "__" + rcs_power_data.columns[key2_combos[i][1]][5:len(rcs_power_data.columns[key2_combos[i][1]])]
        key2_pb_col_names[i] = str
        for j in range(rcs_power_data.shape[0]):
            key2_pb_mat[j,i] = sum(rcs_power_data.values[j,key2_combos[i][0]:key2_combos[i][1]+1])
    key2_pb_col_names.insert(0,"time")
    key2_pb_dataframe = pd.DataFrame(np.hstack((rcs_power_data.time.to_numpy().reshape(-1,1),key2_pb_mat)),columns = key2_pb_col_names)
    toc = time.perf_counter()
    print(f"Code ran in {toc-tic:0.4f} seconds")

    print(f"Creating key3 power bands ...", end = " ")
    tic = time.perf_counter()
    for i in range(len(key3_combos)):
        str = "key3_" + rcs_power_data.columns[key3_combos[i][0]][5:len(rcs_power_data.columns[key3_combos[i][0]])] + "__" + rcs_power_data.columns[key3_combos[i][1]][5:len(rcs_power_data.columns[key3_combos[i][1]])]
        key3_pb_col_names[i] = str
        for j in range(rcs_power_data.shape[0]):
            key3_pb_mat[j,i] = sum(rcs_power_data.values[j,key3_combos[i][0]:key3_combos[i][1]+1])
    key3_pb_col_names.insert(0,"time")
    key3_pb_dataframe = pd.DataFrame(np.hstack((rcs_power_data.time.to_numpy().reshape(-1,1),key3_pb_mat)),columns = key3_pb_col_names)
    toc = time.perf_counter()
    print(f"Code ran in {toc-tic:0.4f} seconds")

    return key0_pb_dataframe,key1_pb_dataframe,key2_pb_dataframe,key3_pb_dataframe


# BOItrDF
#   Calculates the Gaussian process and Acquisition function at a certain itration of Bayesian optimization.
#
#   Inputs:
#       X_samples                   [=] Parameter pairs that have been evaluated. Shape is an (n x m) matrix, 
#                                       where "n" is the number of all evaluations, and "m" is the number of 
#                                       parameters being optimized over.
#       Y_samples                   [=] Evaluations of the function being optimized over at the X_sample points. 
#                                       Shape is a (n x 1) matrix.
#       search_bounds               [=] The min and max values of the parameters to optimize over. Shape is an
#                                       (p x 2) matrix, where p is the number of parameters to optimize over.
#                                       The first parameter is the index numbers in the raw data table to look
#                                       through.  
#       search_bounds_resolution    [=] The spacing between the min and max value of the parameters above.
#       n_initial_samples           [=] How many initial samples there were before starting Bayesian optimization.
#       gp_params                   [=] Parameters of the Gaussian process used during the Bayesian optimization.
#       itr_spacing                 [=] The spacing of the iterations. For example: [1000,100,10,1].
#       itr_reps                    [=] How many times to repeat the itration spacing above. Continuing the example,
#                                       itr_reps = [8,5,4,10] and the spacing above will result in the function
#                                       generating 
#   Outputs:
#
def BOItrDF(X_samples,Y_samples,param_combos,n_initial_samples: int = 0,gp_params: dict = None,itr_spacing: list = None, itr_reps: list = None):
    # Get number of samples
    n_samples = X_samples.shape[0]
    
    # Get number of parameters
    n_param = X_samples.shape[1]

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
        print(f"Recreating GP and AF of iteration {n_itr_indexes[i] + 1}...", end = " ")
        recreate_tic = time.perf_counter()

        # Fit the model to current sampled points
        curr_Y = Y_samples[0:n_itr_indexes[i]]
        gp_model.fit(X = X_samples[0:n_itr_indexes[i],:], y = curr_Y)

        # Get current best Y value
        curr_best_Y = max(curr_Y)

        # Predict the values of the parameter space
        mu,std = gp_model.predict(param_combos, return_std=True)

        # Calculate acquisition funtion values
        af_vals = bo.expectedImprovement(param_combos,curr_best_Y,gp_model,n_param,min_fun=False)

        # Save outputs to a matrix
        GP_mat[i,0] = n_itr_indexes[i]
        GP_mat[i,1:param_combos.shape[0]+1] = mu
        AF_mat[i,0] = n_itr_indexes[i]
        AF_mat[i,1:param_combos.shape[0]+1] = af_vals

        # write out how long this recreation took
        recreate_toc = time.perf_counter()
        print(f"took {recreate_toc-recreate_tic} seconds")

    return GP_mat, AF_mat


# save_multiple_csv
#   Saves the output of the Gaussian process and acquisition function to csv files. If there are more than
#   "n_entries" entries, this function will break up the data into csv files containing "n_entries" per file.
#
#   Inputs:
#       data            [=] Matrix of data from either Gaussian process or acquisition function recreation.
#       n_entries       [=] The number of recreations to save per csv file.
#       base_path       [=] The generic save name and folder location to use. This function will append 
#                           numbers to the end of this name to indicate which file came first.
#   Returns:
#       Nothing.
#
def save_multiple_csv(data,n_entries,base_path: str = os.getcwd()):
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



def runSearch(search_arguments):
    if "bayesOpt" in search_arguments["method_function"].__module__:
        bo_options = {dict_key: search_arguments[dict_key] for dict_key in ("n_itr","gp_params","min_fun","param_combos","key")}
        out_DF,raw_X,raw_Y = search_arguments["method_function"](search_arguments["optimizing_function"],search_arguments["data"],search_arguments["gait_events"],bo_options)

        # Save data to save progress
        out_DF.to_csv(search_arguments["save_name_results"],index = False)

        # Generate dataframes of the Gaussian process and acquisition function so I can plot in MATLAB
        GP_itr,AF_itr = BOItrDF(raw_X,raw_Y,n_initial_samples=out_DF.shape[0]-search_arguments["n_itr"],gp_params = search_arguments["gp_options"],itr_spacing = [100,50,10,1], itr_reps = [8,2,7,30])

        # Save the iteration outputs and parameter combos that were evaluated
        np.savetxt(search_arguments["save_name_param_combos"],search_arguments["param_combos"],delimiter = ",")
        save_multiple_csv(GP_itr,n_entries = 50,base_path = search_arguments["save_name_GP_itr"])
        save_multiple_csv(AF_itr,n_entries = 50,base_path = search_arguments["save_name_AF_itr"])

    elif "gridSearch" in search_arguments["method_function"].__module__:
        gs_options = {dict_key: search_arguments[dict_key] for dict_key in ("param_combos","key")}
        out_DF = search_arguments["method_function"](search_arguments["optimizing_function"],search_arguments["data"],search_arguments["gait_events"],gs_options)
        out_DF.to_csv(search_arguments["save_name_results"],index = False)



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
        elif args[i] == "-search_type":
            search_type = args[i+1]
        elif args[i] == "-n_power_bands":
            n_power_bands = args[i+1]       # Set number of power bands to consider for a control signal. Leave at 1 for now. TODO, add in two power band search
        elif args[i] == "-search_bounds":
            search_bounds_str = args[i+1]       # First row is always the frequency limits to search through
            search_bounds = np.array(literal_eval(search_bounds_str))
        elif args[i] == "-search_bounds_resolution":
            search_bounds_resolution_str = args[i+1]
            search_bounds_resolution = np.array(literal_eval(search_bounds_resolution_str))
        elif args[i] == "-save_path":
            save_path = args[i+1]           # Set the save path if given, otherwise will use the current working directory of the code


    # Create power bands within the first search bound limits
    key0_pb_DF,key1_pb_DF,key2_pb_DF,key3_pb_DF = createPowerBands(rcs_power_data,search_bounds[0,:])

    # generate base save name
    save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_results.csv"
    save_name_GP_itr = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_GP_surf"
    save_name_AF_itr = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_AF_surf"
    save_name_param_combos = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_param_combos.csv"

    # Setup parallel processor pool
    processing_pool = mp.Pool(4)

    if search_type == "bayesian_optimization":
        # Set Bayesian optimization parameters
        # event_strings = []    # TODO add in options to specify gait phase
        # gp_options = {'nu':0.5,'length_scale':np.array([1e-2,2.5]),'length_scale_bounds':"fixed"}
        # gp_options = {'nu':1.5,'length_scale':np.array([1,2.5]),'length_scale_bounds':"fixed"}
        n_itr = 1000
        gp_options = {"nu":1.5,"length_scale":np.array([1,2.5]),"length_scale_bounds":"fixed"}

        bo_option_list = [None] * 4

        for i in range(4):
            param_combos = bo.createParamSpace(eval("key" + str(i) + "_pb_DF"),type = "log")
            bo_options = {"method_function":bo.runBO,"optimizing_function":ta.calcThresholdAccuracyDST,"gait_events":gait_events,"n_itr":n_itr,"gp_params":gp_options,"min_fun":False}
            bo_options["param_combos"] = param_combos
            bo_options["data"] = eval("key" + str(i) + "_pb_DF")
            bo_options["n_itr"] = n_itr
            bo_options["key"] = i
            bo_options["save_name_results"] = save_name_results.replace("full_spec","Bayes_Opt_key"+str(i))
            bo_options["save_name_GP_itr"] = save_name_GP_itr.replace("full_spec","Bayes_Opt_key"+str(i))
            bo_options["save_name_AF_itr"] = save_name_AF_itr.replace("full_spec","Bayes_Opt_key"+str(i))
            bo_options["save_name_param_combos"] = save_name_param_combos.replace("full_spec","Bayes_Opt_key"+str(i))
            bo_option_list[i] = bo_options

        runSearch(bo.runBO,bo_option_list[0])
        # with processing_pool as pool:
        #     pool.map(runSearch,bo_option_list)

        # End run time
        toc = time.perf_counter()
        print(f"Code ran in {toc-tic:0.4f} seconds")
        
    elif search_type == "grid_search":
        grid_search_option_list = [None] * 4

        for i in range(4):
            param_combos = bo.createParamSpace(eval("key" + str(i) + "_pb_DF"),type = "log")
            grid_search_options = {"method_function":gs.runGS,"optimizing_function":ta.calcThresholdAccuracyDST,"gait_events":gait_events}
            grid_search_options["param_combos"] = param_combos
            grid_search_options["data"] = eval("key" + str(i) + "_pb_DF")
            grid_search_options["key"] = i
            grid_search_options["save_name_results"] = save_name_results.replace("full_spec","Grid_Search_key"+str(i))
            grid_search_option_list[i] = grid_search_options
            
        # runSearch(gs.runGS,grid_search_option_list[0])
        with processing_pool as pool:
            pool.map(runSearch,grid_search_option_list)


    # Close the processing pool
    # processing_pool.Close()




"""
Entry point into the script
"""
if __name__=="__main__":
    main()