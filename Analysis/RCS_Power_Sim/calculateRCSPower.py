"""
Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   10/27/2022
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
import glob
import numpy as np
import pandas as pd
import itertools
from libs.rcssim import rcs_sim as rcs # Custom module. Grab latest version from github

"""
Helper functions
"""
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
        if not pd.isnull(formatted_input.power_bands[0]):
            formatted_input.power_bands[0] = np.asarray(formatted_input.power_bands[0].split('|'),dtype = float)
            formatted_input.power_bands[0] = formatted_input.power_bands[0].reshape(int(len(formatted_input.power_bands[0])/2),2)
        return formatted_input
    elif type == "gait_events":
        gait_events = pd.read_csv(file_path)
        return gait_events
    else:
        sys.exit('Data type flag not supported.')


def calcRCSFFT(rcs_data: pd.DataFrame,rcs_settings: pd.DataFrame,pb: bool=False, fs = None, nfft = None, fft_int = None, fft_bitshift = None) -> pd.DataFrame:

    if fs is None:
        fs = rcs_settings.fs_td[0]

    if nfft is None:
        nfft = rcs_settings.FFT_size[0]

    if fft_int is None:
        fft_int = rcs_settings.FFT_interval[0]

    if fft_bitshift is None:
        fft_bitshift = rcs_settings.FFT_bitshift[0]


    # grab hann window values used by the RCS device
    hann_win = rcs.create_hann_window(nfft, percent=100)

    if pb:
        rcs_data_sim = rcs_data[['time','key0','key1','key2','key3']].copy()
        rcs_data_sim['pb1_key0'] = np.nan*np.ones(np.shape(rcs_data['key0']))
        rcs_data_sim['pb2_key0'] = np.nan*np.ones(np.shape(rcs_data['key0']))
        rcs_data_sim['pb1_key1'] = np.nan*np.ones(np.shape(rcs_data['key1']))
        rcs_data_sim['pb2_key1'] = np.nan*np.ones(np.shape(rcs_data['key1']))
        rcs_data_sim['pb1_key2'] = np.nan*np.ones(np.shape(rcs_data['key2']))
        rcs_data_sim['pb2_key2'] = np.nan*np.ones(np.shape(rcs_data['key2']))
        rcs_data_sim['pb1_key3'] = np.nan*np.ones(np.shape(rcs_data['key3']))
        rcs_data_sim['pb2_key3'] = np.nan*np.ones(np.shape(rcs_data['key3']))

        # compute power bands for all time domain channels
        time_domain_channel_names = [col for col in rcs_data.columns if 'key' in col]
        count = 0
        for channel in time_domain_channel_names:
            if np.sum(np.isnan(rcs_data[channel])) != rcs_data.shape[0]:
                data_rcs_unit = rcs.transform_mv_to_rcs(rcs_data_sim[channel].values,rcs_settings.amp_gains[0][count])
                data_rcs_fft, time_pb = rcs.td_to_fft(data_rcs_unit,
                                                    rcs_data_sim['time'].values,
                                                    fs,    # sample rate of the neural data
                                                    nfft, # fft size set on the RCS
                                                    fft_int, # update interval set on the RCS
                                                    hann_win,
                                                    output_in_mv = False)
                
                data_rcs_pb = rcs.fft_to_pb(data_rcs_fft,
                                            fs,
                                            nfft,
                                            fft_bitshift,
                                            band_edges_hz = rcs_settings.power_bands[0][count*2,],
                                            input_is_mv = False)
                pb_sample_mask = np.isin(rcs_data_sim.time,time_pb)
                rcs_data_sim.loc[pb_sample_mask,'pb1_'+channel] = data_rcs_pb

                data_rcs_pb = rcs.fft_to_pb(data_rcs_fft,
                                            fs,
                                            nfft,
                                            fft_bitshift,
                                            band_edges_hz = rcs_settings.power_bands[0][count*2+1,],
                                            input_is_mv = False)
                pb_sample_mask = np.isin(rcs_data_sim.time,time_pb)
                rcs_data_sim.loc[pb_sample_mask,'pb2_'+channel] = data_rcs_pb
            count = count+1
    else:
        # compute power bands for all time domain channels
        rcs_data_sim = pd.DataFrame()
        time_domain_channel_names = [col for col in rcs_data.columns if 'key' in col]
        count = 0
        for channel in time_domain_channel_names:
            if np.sum(np.isnan(rcs_data[channel])) != rcs_data.shape[0]:
                data_rcs_unit = rcs.transform_mv_to_rcs(rcs_data[channel].values,rcs_settings.amp_gains[0][count])
                data_rcs_fft, time_pb = rcs.td_to_fft(data_rcs_unit,
                                                    rcs_data['time'].values,
                                                    fs,    # sample rate of the neural data
                                                    nfft, # fft size set on the RCS
                                                    fft_int, # update interval set on the RCS
                                                    hann_win,
                                                    output_in_mv = False)
                
                data_rcs_power = rcs.fft_to_pb(data_rcs_fft,
                                            fs,
                                            nfft,
                                            fft_bitshift,
                                            input_is_mv = False)
                
                key_freq_names = generateKeyFreqNames(fs,nfft,channel)
                curr_df = pd.DataFrame(data_rcs_power,columns = key_freq_names)
                rcs_data_sim = pd.concat([rcs_data_sim,curr_df],axis=1)
            count = count+1
            
        rcs_data_sim.insert(loc=0,column = "time",value = time_pb)  # add in time to dataframe

    return rcs_data_sim

    
def generateKeyFreqNames(fs: int, nfft: int, key: str):
    # generate center frequencies
    center_freqs = np.arange(nfft/2) * fs/nfft

    # convert to string
    center_freqs = center_freqs.astype(str)
    center_freqs = np.char.replace(center_freqs,".","_")
    center_freqs = [key + "_" + freq for freq in center_freqs]
    return center_freqs



def runSim(data: pd.DataFrame = None, settings: pd.DataFrame = None,options: dict = None,save_path: str = None, file_name: str = None):
    # Calculate fft interval if not defined
    fft_int_calc = lambda x,y: ((1-x)*y*1000)/500 

    # Loop through parameter combos
    for combo in options["param_combos"]:
        # Get fft interval for a specific overlap percentage
        fft_int = int(np.round(fft_int_calc(combo[0],combo[1])))
        
        if fft_int > 150 and fft_int < 500:
            # generate save name
            sim_info = "_fs_%d_nfft_%d_fftbitshift_%d_overlap_%d_fftint_%d" % (options["fs"],combo[1],combo[2],combo[0]*100,fft_int)
            if options["sim_type"] == "fft":
                save_name = save_path + "/" + file_name[0:len(file_name)-4] + sim_info + "_full_spec.csv"
            else:
                save_name = save_path + "/" + file_name[0:len(file_name)-4] + sim_info + "_pb.csv"


            # Check what type of simulation to run and save output
            if options["sim_type"] == "fft":
                out_DF = calcRCSFFT(data,settings,False,options["fs"],combo[1],fft_int,combo[2])
                out_DF.to_csv(save_name,index = False)
                print("Wrote out file: %s" % save_name)
            elif options["sim_type"] == "pb":
                out_DF = calcRCSFFT(data,settings,True,options["fs"],combo[1],fft_int,combo[2])
                out_DF.to_csv(save_name,index = False)
                print("Wrote out file: %s" % save_name)




# main function to load in data and do things
def main(sim_type: str = "fft",save_path: str = os.getcwd()):
    # grab command line arguments
    args = sys.argv[1:]

    # load in data and grab RCS amplitude gains
    for i in range(0,len(args),2):
        if args[i] == "-data":
            rcs_data = loadData(args[i+1],type="data")
            data_file_name = os.path.basename(args[i+1])
        elif args[i] == "-settings":
            rcs_settings = loadData(args[i+1],type="settings")
        elif args[i] == "-data_folder":
            data_folder = args[i+1]
        elif args[i] == "-fs":
            fs = int(args[i+1])
        elif args[i] == "-nfft":
            nfft = int(args[i+1])                # Can be 64, 256, or 1024
        elif args[i] == "-fft_int":
            fft_int = int(args[i+1])
        elif args[i] == "-fft_bitshift":
            fft_bitshift = int(args[i+1])
        elif args[i] == "-overlap":
            overlap_percent = float(args[i+1])
        elif args[i] == "-sim_type":
            sim_type = args[i+1]
        elif args[i] == "-save_path":
            save_path = args[i+1]

    # # debugging prints
    # print(args[0])          # should be RCS data path
    # print(args[1])          # should be RCS settings path
    # print('sucessfullly loaded data')   

    # default values
    if 'fs' not in locals():
        fs = 500

    if 'nfft' not in locals():
        nfft = (64,256,1024)

    if 'fft_bitshift' not in locals():
        fft_bitshift = (0,1,2,3,4,5,6,7)

    if 'overlap_percent' not in locals():
        overlap_percent = (0.5,0.6,0.7,0.8,0.9)

    if 'save_path' not in locals():
        save_path = os.getcwd()

    # Create a list of all possible combos of nfft, bitshift, and overlap percentage
    power_sim_param_combos = list(itertools.product(overlap_percent,nfft,fft_bitshift))

    # Create dictionary to pass into the sim function
    option_dict = {
        "sim_type": sim_type,
        "param_combos": power_sim_param_combos,
        "fs" : fs
    }

    # if a folder was given instead of a single file loop through all files in the folder
    if ('rcs_data' in locals()) and ('rcs_settings' in locals()) and ('data_folder' not in locals()):
        runSim(data = rcs_data, settings = rcs_settings,options = option_dict,save_path = save_path, file_name = data_file_name)
    elif ('rcs_data' not in locals()) and ('rcs_settings' not in locals()) and ('data_folder' in locals()):
        left_file_list = glob.glob(data_folder + "/*LEFT.csv")
        right_file_list = glob.glob(data_folder + "/*RIGHT.csv")
        file_list = left_file_list+right_file_list
        for file in file_list:
            rcs_data = loadData(file,type="data")
            rcs_settings = loadData(file[0:len(file)-4]+"_SETTINGS.csv",type = "settings")
            runSim(data = rcs_data, settings=rcs_settings, options = option_dict,save_path = save_path,file_name=os.path.basename(file))
    else:
        raise ValueError("Function input error. Either have the rcs data and setting file as input with the data folder not defined, or have the data folder defined but not rcs data and setting file.")

   

""" 
Entry point into the script
"""
if __name__=="__main__":
    main()