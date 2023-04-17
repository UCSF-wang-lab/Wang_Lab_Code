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
from libs.optimizingFunctions import thresholdAccuracy as ta # function to evaluate new settings
from libs.powerSim import powerSimRCS as psr # RCS power sim file
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
    elif type == "settings":
        col_names = pd.read_csv(file_path,nrows = 0).columns
        types_dict = {"fs_td":int,"FFT_size":int,"FFT_interval":int,"FFT_bitshift":int,"rise_time":int,"fall_time":int}
        types_dict.update({col: str for col in col_names if col not in types_dict})
        raw_input = pd.read_csv(file_path,dtype = types_dict)
        formatted_input = raw_input.copy()
        formatted_input.amp_gains[0] = np.asarray(formatted_input.amp_gains[0].split('|'),dtype = int)
        formatted_input.power_bands[0] = np.asarray(formatted_input.power_bands[0].split('|'),dtype = float)
        formatted_input.power_bands[0] = formatted_input.power_bands[0].reshape(int(len(formatted_input.power_bands[0])/2),2)
        return formatted_input
    elif type == "gait_events":
        gait_events = pd.read_csv(file_path)
        return gait_events
    else:
        sys.exit('Data type flag not supported.')



def sortGaitEvents(gait_events,start_event):
    if gait_events.columns[0] == start_event:
        return gait_events
    else:
        match start_event:
            case "LHS":
                gait_event_order = ["LHS","RTO","RHS","LTO"]
            case "RTO":
                gait_event_order = ["RTO","RHS","LTO","LHS"]
            case "RHS":
                gait_event_order = ["RHS","LTO","LHS","RTO"]
            case "LTO":
                gait_event_order = ["LTO","LHS","RTO","RHS"]

        shift_ind = list(gait_events.columns).index(start_event)
        sorted_gait_events = np.empty((gait_events.shape[0]+1,gait_events.shape[1]))
        sorted_gait_events.fill(np.nan)

        # Create initial table
        for i in range(len(gait_event_order)-shift_ind):        
            sorted_gait_events[1:,i] = gait_events[gait_event_order[i]]

        for j in range(len(gait_event_order)-shift_ind,len(gait_event_order)):
            sorted_gait_events[0:-1,j] = gait_events[gait_event_order[j]]

        # Adjust table so that events are time reasonable
        end_of_table = False
        count = 0
        while not end_of_table:
            for k in range(1,sorted_gait_events.shape[1]):
                event_diff = sorted_gait_events[count,k] - sorted_gait_events[count,k-1]
                if event_diff > 1.5:
                    X = np.array(sorted_gait_events[0:count+1,:])
                    Y = np.array(sorted_gait_events[count+1:sorted_gait_events.shape[0],:])
                    Z = np.array(X[X.shape[0]-1,:])

                    X[X.shape[0]-1,range(k,X.shape[1])] = np.nan
                    Z[0:k] = np.nan

                    sorted_gait_events = np.vstack((X,Z,Y))

            if count + 1 >= sorted_gait_events.shape[0]:
                end_of_table = True
            
            count+=1
            
        
        # Remove rows with all nan
        remove_inds = np.where(np.isnan(sorted_gait_events).all(axis=1))
        sorted_gait_events = np.delete(sorted_gait_events,remove_inds,axis=0)
        sorted_gait_event_table = pd.DataFrame(sorted_gait_events,columns = gait_event_order)
        
        return sorted_gait_event_table



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
    gs_options = {dict_key: search_arguments[dict_key] for dict_key in ("param_combos","key","search_type","gait_phase","search_method","FFT_size","FFT_interval","FFT_bitshift","FFT_overlap")}
    out_DF = search_arguments["method_function"](search_arguments["optimizing_function"],search_arguments["data"],search_arguments["gait_events"],gs_options)
    out_DF["FFT_overlap"] = search_arguments["FFT_overlap"]
    out_DF["FFT_size"] = search_arguments["FFT_size"]
    out_DF["FFT_interval"] = search_arguments["FFT_interval"]
    out_DF["FFT_bitshift"] = search_arguments["FFT_bitshift"]

    if search_arguments["search_method"] == "grid_search_modified":
        out_DF.rename(columns = {'Frequency Band Ind':'Frequency Band'},inplace = True)
        out_DF["Frequency Band"] = list(search_arguments["data"].columns)[1]
    else:
        column_names = list(search_arguments["data"].columns)
        myfile = open(search_arguments["save_name_freq_bands"],"w+",newline='')
        with myfile:
            w = csv.writer(myfile)
            for name in column_names:
                w.writerow([name])

    out_DF.to_csv(search_arguments["save_name_results"],mode = "a", header = not os.path.exists(search_arguments["save_name_results"]),index = False)

        



# main function to start aDBS settings search
def main(n_power_bands: int = 1,parallized_flag: bool = False,save_path: str = os.getcwd()):

    # grab command line arguments
    args = sys.argv[1:]

    # load in data and grab RCS amplitude gains
    for i in range(0,len(args),2):
        if args[i] == "-data":
            rcs_power_data = loadData(args[i+1],type="data")
            data_file_name = os.path.basename(args[i+1])        # grabs the current filename to be used as a filename.
        elif args[i] == "-settings":
            rcs_base_settings = loadData(args[i+1],type="settings")
        elif args[i] == "-gait_events_table":
            gait_events = loadData(args[i+1],type="gait_events")
        elif args[i] == "-gait_events":
            gait_events_str = args[i+1]
        elif args[i] == "-gait_phase":
            gait_phase = args[i+1]
        elif args[i] == "-search_method":     # Options: bayesian_optimization|grid_search
            search_method = args[i+1]
        elif args[i] == "-search_type":       # Options: all_combo|specific_freq
            search_type = args[i+1]
        elif args[i] == "-frequency_range":
            frequency_range_str = args[i+1]       # Depends on search_type value
            frequency_range = np.array(literal_eval(frequency_range_str))
        elif args[i] == "-keys": # Options key0|key1|key2|key3|all. To include more than one key, used a "," to separate keys
            if args[i+1] == "all":
                search_keys_str = ["key0","key1","key2","key3"]
            else:
                search_keys_str = args[i+1]
        elif args[i] == "-parallized":  # boolean string
            parallized_flag = args[i+1]
        elif args[i] == "-save_path":
            save_path = args[i+1]           # Set the save path if given, otherwise will use the current working directory of the code

    # Set search key values
    search_keys = search_keys_str.split(",")

    # check if parallized_flag exist, if not, set to false
    if "parallized_flag" not in locals():
        parallized_flag = "False"

    # generate base save name
    save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_" + search_method + "_results.csv"

    # Trim the gait event table to the relavant events
    gait_events_sorted = sortGaitEvents(gait_events,gait_events_str.split(",")[0])
    gait_events_trimmed = gait_events_sorted[gait_events_str.split(",")]

    # lambda function to calculate fft interval for simulation
    fft_point_size = [64,256,1024]
    overlap_percent = [0.5,0.6,0.7,0.8,0.9]
    bit_shift = [0,1,2,3,4,5,6,7]
    power_sim_param_combos = list(itertools.product(overlap_percent,fft_point_size,bit_shift))
    fft_int_calc = lambda x,y: ((1-x)*y*1000)/500

    # timing measurement for debugging
    tic = time.perf_counter()

    # Loop through each key and create the data to be searched over
    for combo in power_sim_param_combos:
        fft_int_val = fft_int_calc(combo[0],combo[1])
        # Check to make sure fft_inverval value is greater than or equal to the minimum value.
        if fft_int_val >= 50:
            rcs_settings = rcs_base_settings.copy()
            rcs_settings.FFT_size = combo[1]
            rcs_settings.FFT_interval = fft_int_val
            rcs_settings.FFT_bitshift = combo[2]
            
            # Calculate RCS power
            rcs_power_sim = psr.calcRCSFFT(rcs_power_data,rcs_settings)

            # Create power bands
            if (frequency_range.shape[0] == len(search_keys)) | (frequency_range.shape[0] == 1):
                pb_DataFrame_dict = createPowerBands(rcs_power_sim,search_keys,frequency_range,search_type)
            else:
                sys.exit("Number of frequency ranges given is more than the number of keys to search.")

            if parallized_flag == "True":
                n_total_freq_bands = 0
                for i in range(len(search_keys)):
                    n_total_freq_bands += pb_DataFrame_dict[search_keys[i]].shape[1]-1
                grid_search_option_list = [None] * n_total_freq_bands
        
            count = 0
            for i in range(len(search_keys)):
                # Loop through each frequency band and create a search parameter for it.
                for j in range(1,pb_DataFrame_dict[search_keys[i]].shape[1]):

                    # Set search options
                    grid_search_options = {"method_function":gs.runGS,"gait_events":gait_events_trimmed}
                    grid_search_options["data"] = pb_DataFrame_dict[search_keys[i]].iloc[:,[0,j]]
                    grid_search_options["key"] = search_keys[i]
                    grid_search_options["save_name_results"] = save_name_results.replace("results",search_keys[i]+"_results")
                    grid_search_options["save_name_freq_bands"] = grid_search_options["save_name_results"].replace("results","freq_bands")
                    grid_search_options["search_type"] = search_type
                    grid_search_options["search_method"] = search_method
                    grid_search_options["FFT_overlap"] = combo[0]
                    grid_search_options["FFT_size"] = combo[1]
                    grid_search_options["FFT_interval"] = fft_int_val
                    grid_search_options["FFT_bitshift"] = combo[2]

                    match gait_phase:
                        case "DST":
                            grid_search_options["optimizing_function"] = ta.calcThresholdAccuracyDST
                            grid_search_options["gait_phase"] = 'DST'
                        case "Swing":
                            grid_search_options["optimizing_function"] = ta.calcThresholdAccuracySwingPhase
                            grid_search_options["gait_phase"] = 'Swing'
                        case "Stance":
                            grid_search_options["optimizing_function"] = ta.calcThresholdAccuracyStancePhase
                            grid_search_options["gait_phase"] = 'Stance'

                    # Set parameter space options
                    data_quantiles = np.quantile(grid_search_options["data"].iloc[:,1],np.arange(0.10,0.95,0.025))
                    param_combos = np.hstack((np.ones((data_quantiles.size,1)),data_quantiles.reshape((-1,1))))
                    grid_search_options["param_combos"] = param_combos

                    if parallized_flag == "True":
                        grid_search_option_list[count] = grid_search_options
                        count += 1
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
        processing_pool.close()



"""
Entry point into the script
"""
if __name__=="__main__":
    main()