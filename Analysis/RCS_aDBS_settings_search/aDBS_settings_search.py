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
import sys
import os
import numpy as np
import pandas as pd
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
    else:
        sys.exit('In proper data load type passed in.')



# main function to load in data and do things
def main(n_power_bands: int = 1,save_path: str = os.getcwd()):
    # grab command line arguments
    args = sys.argv[1:]

    # load in data and grab RCS amplitude gains
    for i in range(0,len(args),2):
        if args[i] == "-data":
            rcs_power_data = loadData(args[i+1],type="data")
            data_file_name = os.path.basename(args[i+1])
        elif args[i] == "-n_power_bands":
            n_power_bands = args[i+1]
        elif args[i] == "-save_path":
            save_path = args[i+1]

    # generate save name
    save_name = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_BO_surf.csv"

    # Check what type of simulation to run and save output
    out_DF = bo.runBO()
    out_DF.to_csv(save_name,index = False)



"""
Entry point into the script
"""
if __name__=="__main__":
    main()