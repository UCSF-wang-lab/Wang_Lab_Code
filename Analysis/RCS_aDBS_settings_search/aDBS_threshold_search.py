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
import csv
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


def createPowerBands(rcs_power_data,keys,frequency_range,band_type):
    # If number of elements in keys is larger than the number
    # of rows in the frequency range, repeat the elements to match
    # the number of elements in keys.
    if len(keys) > frequency_range.shape[0]:
        frequency_range = np.repeat(frequency_range,len(keys),axis=0)

    # Allocate memory
    key_pb_dataframe_dict = dict.fromkeys(keys)

    # Get column names, which contain the frequency value
    col_names = rcs_power_data.columns

    # Loop through column names and see if they should be include in the power band creations
    if band_type == "all_combo":
        for k in range(len(keys)):
            key_str = keys[k]
            key_inds = np.zeros(512,dtype=int)
            key_count = 0
            count = 0
            print(f"Creating {key_str} power bands ...", end = " ")
            for i in col_names:
                if key_str in i:
                    freq_val_str = i[5:len(i)]
                    freq_val_str = freq_val_str.replace("_",".")
                    freq_val = literal_eval(freq_val_str)
                    if (freq_val >= frequency_range[k][0]) & (freq_val <= frequency_range[k][1]):
                            key_inds[key_count] = count
                            key_count += 1
                count += 1

            # Remove extra indices that do not matter
            key_inds = key_inds[0:key_count]

            # Generate all combinations of indexes from each key
            key_combos = list(itertools.combinations_with_replacement(key_inds,2))

            # Allocate memory for power band column names
            key_pb_col_names = [None]*len(key_combos)

            # Allocate memory for power band matrices
            key_pb_mat = np.zeros((rcs_power_data.shape[0],len(key_combos)))

            # Create power band matrix
            tic = time.perf_counter()
            for i in range(len(key_combos)):
                str = key_str + "_" + rcs_power_data.columns[key_combos[i][0]][5:len(rcs_power_data.columns[key_combos[i][0]])] + "__" + rcs_power_data.columns[key_combos[i][1]][5:len(rcs_power_data.columns[key_combos[i][1]])]
                key_pb_col_names[i] = str
                for j in range(rcs_power_data.shape[0]):
                    key_pb_mat[j,i] = sum(rcs_power_data.values[j,key_combos[i][0]:key_combos[i][1]+1])
            key_pb_col_names.insert(0,"time")
            key_pb_dataframe = pd.DataFrame(np.hstack((rcs_power_data.time.to_numpy().reshape(-1,1),key_pb_mat)),columns = key_pb_col_names)
            toc = time.perf_counter()
            print(f"Code ran in {toc-tic:0.4f} seconds")

            key_pb_dataframe_dict[key_str] = key_pb_dataframe

    elif band_type == "specific_freq":
        # Allocate memory for power band matrices and column names
        key_pb_mat = np.zeros((rcs_power_data.shape[0],1))
        key_pb_col_names = [None]*2
        key_pb_col_names[0] = "time"

        # Create power band matrices
        for k in range(len(keys)):
            key_str = keys[k]
            key_inds = np.zeros(512,dtype=int)
            key_count = 0
            count = 0
            print(f"Extracting {keys[k]}: {frequency_range[k,0]} - {frequency_range[k,1]} power band ...", end = " ")
            tic = time.perf_counter()
            for i in col_names:
                if key_str in i:
                    freq_val_str = i[5:len(i)]
                    freq_val_str = freq_val_str.replace("_",".")
                    freq_val = literal_eval(freq_val_str)
                    if (freq_val >= frequency_range[k][0]) & (freq_val <= frequency_range[k][1]):
                            key_inds[key_count] = count
                            key_count += 1
                count += 1

            # Remove extra indices that do not matter
            key_inds = key_inds[0:key_count]

            # Extract the frequency band
            key_pb_col_names[1] = key_str + "_" + rcs_power_data.columns[key_inds[0]][5:len(rcs_power_data.columns[key_inds[0]])] + "__" + rcs_power_data.columns[key_inds[-1]][5:len(rcs_power_data.columns[key_inds[-1]])]
            key_pb_mat = np.sum(rcs_power_data.values[:,key_inds[0]:key_inds[-1]],axis=1).reshape((rcs_power_data.shape[0],1))
            key_pb_dataframe = pd.DataFrame(np.hstack((rcs_power_data.time.to_numpy().reshape(-1,1),key_pb_mat)),columns = key_pb_col_names)
            toc = time.perf_counter()
            print(f"Code ran in {toc-tic:0.4f} seconds")

            key_pb_dataframe_dict[key_str] = key_pb_dataframe

    return key_pb_dataframe_dict



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
        bo_options = {dict_key: search_arguments[dict_key] for dict_key in ("n_itr","gp_params","min_fun","param_combos","key","search_type")}
        out_DF,X,Y_weighted,Y_dst,Y_raw = search_arguments["method_function"](search_arguments["optimizing_function"],search_arguments["data"],search_arguments["gait_events"],bo_options)

        # Save data to save progress
        out_DF.to_csv(search_arguments["save_name_results"],index = False)

        # Generate dataframes of the Gaussian process and acquisition function so I can plot in MATLAB
        # GP_itr,AF_itr = bo.BOItrDF(X,Y_weighted,n_initial_samples=out_DF.shape[0]-search_arguments["n_itr"],gp_params = search_arguments["gp_params"],itr_spacing = [100,50,10,1], itr_reps = [8,2,7,30])
        GP_itr,AF_itr = bo.BOItrDF(X,Y_weighted,search_arguments["param_combos"],n_initial_samples=out_DF.shape[0]-search_arguments["n_itr"],gp_params = search_arguments["gp_params"])

        # Save the iteration outputs and parameter combos that were evaluated
        np.savetxt(search_arguments["save_name_param_combos"],search_arguments["param_combos"],delimiter = ",")
        save_multiple_csv(GP_itr,n_entries = 100,base_path = search_arguments["save_name_GP_itr"])
        save_multiple_csv(AF_itr,n_entries = 100,base_path = search_arguments["save_name_AF_itr"])

    elif "gridSearch" in search_arguments["method_function"].__module__:
        gs_options = {dict_key: search_arguments[dict_key] for dict_key in ("param_combos","key","search_type")}
        out_DF = search_arguments["method_function"](search_arguments["optimizing_function"],search_arguments["data"],search_arguments["gait_events"],gs_options)
        out_DF.to_csv(search_arguments["save_name_results"],index = False)

        column_names = list(search_arguments["data"].columns)
        myfile = open(search_arguments["save_name_freq_bands"],"w+",newline='')
        with myfile:
            w = csv.writer(myfile)
            for name in column_names:
                w.writerow([name])



# main function to start aDBS settings search
def main(n_power_bands: int = 1,parallized_flag: bool = False,save_path: str = os.getcwd()):

    # grab command line arguments
    args = sys.argv[1:]

    # load in data and grab RCS amplitude gains
    for i in range(0,len(args),2):
        if args[i] == "-data":
            rcs_power_data = loadData(args[i+1],type="data")
            data_file_name = os.path.basename(args[i+1])        # grabs the current filename to be used as a filename.
        elif args[i] == "-gait_events":
            gait_events = loadData(args[i+1],type="gait_events")
        elif args[i] == "-search_method":     # Options: bayesian_optimization|grid_search
            search_method = args[i+1]
        elif args[i] == "-search_type":       # Options: all_combo|specific_freq
            search_type = args[i+1]
        elif args[i] == "-n_power_bands":
            n_power_bands = args[i+1]       # Set number of power bands to consider for a control signal. Leave at 1 for now. TODO, add in two power band search
        elif args[i] == "-frequency_range":
            frequency_range_str = args[i+1]       # Depends on search_type value
            frequency_range = np.array(literal_eval(frequency_range_str))
        elif args[i] == "-parameter_space_spacing_type":    # Depends on search_bounds does not have to be set, will assume some default values
            spacing_type_str = args[i+1]
        elif args[i] == "-parameter_space_spacing_n":
            spacing_str = args[i+1]
        elif args[i] == "-keys": # Options key0|key1|key2|key3|all. To include more than one key, used a "," to separate keys
            if args[i+1] == "all":
                search_keys_str = ["key0","key1","key2","key3"]
            else:
                search_keys_str = args[i+1]
        elif args[i] == "-parallized":  # boolean string
            parallized_flag = args[i+1]
        elif args[i] == "-save_path":
            save_path = args[i+1]           # Set the save path if given, otherwise will use the current working directory of the code

    # Set spacing and search key values
    spacing = np.array(literal_eval(spacing_str))
    spacing_type = spacing_type_str.split(",")
    search_keys = search_keys_str.split(",")

    # check if parallized_flag exist, if not, set to false
    if "parallized_flag" not in locals():
        parallized_flag = "False"

    # Create power bands
    if (frequency_range.shape[0] == len(search_keys)) | (frequency_range.shape[0] == 1):
        pb_DataFrame_dict = createPowerBands(rcs_power_data,search_keys,frequency_range,search_type)
    else:
        sys.exit("Number of frequency ranges given is more than the number of keys to search.")

    # generate base save name
    if search_method == "bayesian_optimization":
        save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_results.csv"
        save_name_GP_itr = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_GP_surf"
        save_name_AF_itr = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_AF_surf"
        save_name_param_combos = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_param_combos.csv"
    elif search_method == "grid_search":
        save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_results.csv"

    # timing measurement for debugging
    tic = time.perf_counter()

    if search_method == "bayesian_optimization":
        # Set Bayesian optimization parameters
        # event_strings = []    # TODO add in options to specify gait phase
        # gp_options = {'nu':0.5,'length_scale':np.array([1e-2,2.5]),'length_scale_bounds':"fixed"}
        # gp_options = {'nu':1.5,'length_scale':np.array([1,2.5]),'length_scale_bounds':"fixed"}

        if search_type == "all_combo":
            gp_options = {"nu":0.5,"length_scale":np.array([1,2.5]),"length_scale_bounds":"fixed"}
        elif search_type == "specific_freq":
            gp_options = {"nu":0.5,"length_scale":np.array([2.5]),"length_scale_bounds":"fixed"}

        if parallized_flag == "True":
            bo_option_list = [None] * len(search_keys)

        for i in range(len(search_keys)):
            # Reset spacing type and spacing values
            spacing_type = spacing_type_str.split(",")
            spacing = np.array(literal_eval(spacing_str))

            # Set search parameters
            bo_options = {"method_function":bo.runBO,"optimizing_function":ta.calcThresholdAccuracyDST,"gait_events":gait_events,"gp_params":gp_options,"min_fun":False}
            bo_options["data"] = pb_DataFrame_dict[search_keys[i]]
            bo_options["key"] = search_keys[i]
            bo_options["search_type"] = search_type
            bo_options["save_name_results"] = save_name_results.replace("full_spec","Bayes_Opt_"+search_keys[i])
            bo_options["save_name_GP_itr"] = save_name_GP_itr.replace("full_spec","Bayes_Opt_"+search_keys[i])
            bo_options["save_name_AF_itr"] = save_name_AF_itr.replace("full_spec","Bayes_Opt_"+search_keys[i])
            bo_options["save_name_param_combos"] = save_name_param_combos.replace("full_spec","Bayes_Opt_"+search_keys[i])

            # Set parameter space options
            data_quantiles = np.quantile(pb_DataFrame_dict[search_keys[i]].values[:,1:pb_DataFrame_dict[search_keys[i]].shape[1]],[0.05,0.95])
            if search_type == "all_combo":
                spacing_type.insert(0,"linear")
                spacing = np.insert(spacing,0,pb_DataFrame_dict[search_keys[i]].shape[1]-1,axis=0)
                parameter_range = np.zeros((2,2))
                parameter_range[0,:] = [1,pb_DataFrame_dict[search_keys[i]].shape[1]-1]
                if (data_quantiles[0] == 0) & (spacing_type[1] == "log"):
                    data_quantiles[0] = 1e-3
                    parameter_range[1,:] = data_quantiles
                else:
                    parameter_range[1,:] = data_quantiles

                param_combos = bo.createParamSpace(parameter_range,spacing_type,n_spacings=spacing)
                bo_options["param_combos"] = param_combos
            elif search_type == "specific_freq":
                parameter_range = data_quantiles.reshape((-1,2))
                param_combos = bo.createParamSpace(parameter_range,spacing_type,n_spacings=spacing)
                bo_options["param_combos"] = param_combos

            # Set the number of iterations to be 50% of the total number of parameter combos
            bo_options["n_itr"] = round(param_combos.shape[0]*0.5)

            if parallized_flag == "True":
                bo_option_list[i] = bo_options
            else:
                runSearch(bo_options)

        if parallized_flag == "True":
            processing_pool = mp.Pool(len(search_keys))
            with processing_pool as pool:
                pool.map(runSearch,bo_option_list)

        # End run time
        toc = time.perf_counter()
        print(f"Code ran in {toc-tic:0.4f} seconds")
        
    elif search_method == "grid_search":
        if parallized_flag == "True":
            grid_search_option_list = [None] * len(search_keys)

        for i in range(len(search_keys)):
            # Set search options
            grid_search_options = {"method_function":gs.runGS,"optimizing_function":ta.calcThresholdAccuracyDST,"gait_events":gait_events}
            grid_search_options["data"] = pb_DataFrame_dict[search_keys[i]]
            grid_search_options["key"] = search_keys[i]
            grid_search_options["save_name_results"] = save_name_results.replace("full_spec","Grid_Search_"+search_keys[i])
            grid_search_options["save_name_freq_bands"] = grid_search_options["save_name_results"].replace("results","freq_bands")
            grid_search_options["search_type"] = search_type

            # Set parameter space options
            data_quantiles = np.quantile(pb_DataFrame_dict[search_keys[i]].values[:,1:pb_DataFrame_dict[search_keys[i]].shape[1]],[0.05,0.95])
            if search_type == "all_combo":
                spacing_type.insert(0,"linear")
                spacing = np.insert(spacing,0,pb_DataFrame_dict[search_keys[i]].shape[1]-1,axis=0)
                parameter_range = np.zeros((2,2))
                parameter_range[0,:] = [1,pb_DataFrame_dict[search_keys[i]].shape[1]-1]
                if (data_quantiles[0] == 0) & (spacing_type[1] == "log"):
                    data_quantiles[0] = 1e-3
                    parameter_range[1,:] = data_quantiles
                else:
                    parameter_range[1,:] = data_quantiles

                param_combos = bo.createParamSpace(parameter_range,spacing_type,n_spacings=spacing)
                grid_search_options["param_combos"] = param_combos
            elif search_type == "specific_freq":
                parameter_range = data_quantiles.reshape((-1,2))
                param_combos = bo.createParamSpace(parameter_range,spacing_type,n_spacings=spacing)
                grid_search_options["param_combos"] = param_combos

            if parallized_flag == "True":
                grid_search_option_list[i] = grid_search_options
            else:
                runSearch(grid_search_options)            
            
        if parallized_flag == "True":
            processing_pool = mp.Pool(len(search_keys))
            with processing_pool as pool:
                pool.map(runSearch,grid_search_option_list)

        # End run time
        toc = time.perf_counter()
        print(f"Code ran in {toc-tic:0.4f} seconds")

    # Close the processing pool
    if parallized_flag == "True":
        processing_pool.Close()



"""
Entry point into the script
"""
if __name__=="__main__":
    main()