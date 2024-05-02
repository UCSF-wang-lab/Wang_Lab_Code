import sys
import os
import numpy as np
import pandas as pd
import itertools
from ast import literal_eval



# createPowerBands
#   Generate power band sets for all keys under considerations. The power bands the index of the correct column, represented as an integer
#
#   Inputs:     col_names       [=] All the column names of a pandas dataframe.
#               included_keys   [=] All the keys under consideration as a list.
#               n_features      [=] How many power bands to piece together
#
#   Outputs:    powerBands      [=] List of all power band sets given the included keys and number of features
#               class_col_ind   [=] Column index of the class column. Useful when generating smaller dataframes for optimization
def createPowerBands(col_names,included_keys,n_features):
    powerBands = []
    singlePowerBands = []
    for key in included_keys:
        includeColumns = []
        for i,j in enumerate(list(col_names)):
            if key in j:
                includeColumns.append(i)
            elif (j == "Class") & (j not in locals()):
                class_col_ind = i
        singlePowerBands += (list(itertools.combinations(includeColumns,2)))
    
    exclude_inds = []
    for row,band in enumerate(singlePowerBands):
        if band[1]-band[0] >= 64:
            exclude_inds.append(row)

    singlePowerBandsFilt = []
    for row,band in enumerate(singlePowerBands):
        if row not in exclude_inds:
            singlePowerBandsFilt.append(band)

    powerBands = list(itertools.combinations(singlePowerBandsFilt,n_features))

    return powerBands, class_col_ind



def main():
    # grab command line arguments
    args = sys.argv[1:]

    # Parses optional inputs after the first three inputs
    for i in range(0,len(args),2):
        if args[i] == "-data":
            rcs_power_data = pd.read_csv(args[i+1])         # read in csv table
            data_file_name = os.path.basename(args[i+1])    # grabs the current filename to be used as a filename for saving.
        elif args[i] == "-n_features":     # Can be any value 1-4
            n_features = int(args[i+1])
        elif args[i] == "-keys_left": # Options key0|key1|key2|key3|all. To include more than one key, used a "," to separate keys. This will be used to filter out columns that are not in the key list.
            if args[i+1] == "all":
                search_keys_str_left = "key0,key1,key2,key3"
            else:
                search_keys_str_left = args[i+1]
        elif args[i] == "-keys_right": # Options key0|key1|key2|key3|all. To include more than one key, used a "," to separate keys. This will be used to filter out columns that are not in the key list.
            if args[i+1] == "all":
                search_keys_str_right = "key0,key1,key2,key3"
            else:
                search_keys_str_right = args[i+1]   
        elif args[i] == "-save_path":
            save_path = args[i+1]       # Set the save path if given, otherwise will use the current working directory of the code

    # Default values
    if "n_features" not in locals():    # set number of features to 1 by default
        n_features = 1

    if "search_keys_str_left" not in locals():
        search_keys_str_left = "key0"
        search_keys_left = search_keys_str_left.split(",") # Set search key values
    else:
        search_keys_left = search_keys_str_left.split(",") # Set search key values

    if "search_keys_str_right" not in locals():
        search_keys_str_right = "key0"
        search_keys_right = search_keys_str_right.split(",") # Set search key values
    else:
        search_keys_right = search_keys_str_right.split(",") # Set search key values
    
    if "save_path" not in locals():
        save_path = os.getcwd()
        save_name = save_path + "/" + data_file_name[0:len(data_file_name)-4] + f"_{n_features}features" + "_aDBS_search_list.csv"   # generate base save name
    else:
        save_name = save_path + "/" + data_file_name[0:len(data_file_name)-4] + f"_{n_features}features" + "_aDBS_search_list.csv"   # generate base save name

    # Create combos of FFT settings
    FFT_combos = np.array([[p1,p2,p3] for p1 in np.unique(rcs_power_data.Bitshift) for p2 in np.unique(rcs_power_data.FFTInterval) for p3 in np.unique(rcs_power_data.OverlapPercent)])
        
    # # Determine which columns to consider when creating parameter space
    # # and then create list of sets to optimize
    # hemisphere_vec = []
    # fft_int_vec = []
    # fft_overlap_vec = []
    # fft_bitshift_vec = []
    # param_set_str_vec = []
    # if "left" in set(rcs_power_data.Side):
    #     power_band_combos_left,class_col_ind = createPowerBands(rcs_power_data.columns,search_keys_left,n_features)
    #     for fft_combo in FFT_combos:
    #         include_ind = np.where((rcs_power_data.Side == "left") & (rcs_power_data.Bitshift == fft_combo[0]) & (rcs_power_data.FFTInterval == fft_combo[1]) & (rcs_power_data.OverlapPercent == fft_combo[2]))
    #         if include_ind[0].size>0:
    #             for param_set in power_band_combos_left:
    #                 hemisphere_vec.append("left")
    #                 fft_int_vec.append(fft_combo[1])
    #                 fft_overlap_vec.append(fft_combo[2])
    #                 fft_bitshift_vec.append(fft_combo[0])
    #                 param_set_str_vec.append(str([class_col_ind] + list(param_set)))

    # if "right" in set(rcs_power_data.Side):
    #     power_band_combos_right,class_col_ind = createPowerBands(rcs_power_data.columns,search_keys_right,n_features)
    #     for fft_combo in FFT_combos:
    #         include_ind = np.where((rcs_power_data.Side == "right") & (rcs_power_data.Bitshift == fft_combo[0]) & (rcs_power_data.FFTInterval == fft_combo[1]) & (rcs_power_data.OverlapPercent == fft_combo[2]))
    #         if include_ind[0].size>0:
    #             for param_set in power_band_combos_right:
    #                 hemisphere_vec.append("right")
    #                 fft_int_vec.append(fft_combo[1])
    #                 fft_overlap_vec.append(fft_combo[2])
    #                 fft_bitshift_vec.append(fft_combo[0])
    #                 param_set_str_vec.append(str([class_col_ind] + list(param_set)))

    # aDBS_search_DF = pd.DataFrame(np.column_stack([hemisphere_vec,fft_bitshift_vec,fft_int_vec,fft_overlap_vec,param_set_str_vec]),columns = ['Hemisphere','FFTBitshift','FFTInt','FFTOverlap','FreqBandCols'])
    # aDBS_search_DF.to_csv(save_name,index=False)

    # Determine which columns to consider when creating parameter space
    # and then write combo to file

    search_space_file = open(save_name,'w')
    search_space_file.write("Hemisphere,FFTBitshift,FFTInt,FFTOverlap,FreqBandCols\n")
    if "left" in set(rcs_power_data.Side):
        power_band_combos_left,class_col_ind = createPowerBands(rcs_power_data.columns,search_keys_left,n_features)
        for fft_combo in FFT_combos:
            include_ind = np.where((rcs_power_data.Side == "left") & (rcs_power_data.Bitshift == fft_combo[0]) & (rcs_power_data.FFTInterval == fft_combo[1]) & (rcs_power_data.OverlapPercent == fft_combo[2]))
            if include_ind[0].size>0:
                for param_set in power_band_combos_left:
                    line_part1 = f"left,{fft_combo[0]},{fft_combo[1]},{fft_combo[2]},"
                    line_part2 = "\"" + str([class_col_ind] + list(param_set)) + "\"\n"
                    line = line_part1 + line_part2
                    search_space_file.write(line)
                    # hemisphere_vec.append("left")
                    # fft_int_vec.append(fft_combo[1])
                    # fft_overlap_vec.append(fft_combo[2])
                    # fft_bitshift_vec.append(fft_combo[0])
                    # param_set_str_vec.append(str([class_col_ind] + list(param_set)))

    if "right" in set(rcs_power_data.Side):
        power_band_combos_right,class_col_ind = createPowerBands(rcs_power_data.columns,search_keys_right,n_features)
        for fft_combo in FFT_combos:
            include_ind = np.where((rcs_power_data.Side == "right") & (rcs_power_data.Bitshift == fft_combo[0]) & (rcs_power_data.FFTInterval == fft_combo[1]) & (rcs_power_data.OverlapPercent == fft_combo[2]))
            if include_ind[0].size>0:
                for param_set in power_band_combos_right:
                    line_part1 = f"right,{fft_combo[0]},{fft_combo[1]},{fft_combo[2]},"
                    line_part2 = "\"" + str([class_col_ind] + list(param_set)) + "\"\n"
                    line = line_part1 + line_part2
                    search_space_file.write(line)
                    # hemisphere_vec.append("right")
                    # fft_int_vec.append(fft_combo[1])
                    # fft_overlap_vec.append(fft_combo[2])
                    # fft_bitshift_vec.append(fft_combo[0])
                    # param_set_str_vec.append(str([class_col_ind] + list(param_set)))
    
    search_space_file.close()

"""
Entry point into the script
"""
if __name__=="__main__":
    main()