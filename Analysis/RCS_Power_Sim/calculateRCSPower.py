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
import numpy as np
import pandas as pd
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
        formatted_input.power_bands[0] = np.asarray(formatted_input.power_bands[0].split('|'),dtype = float)
        formatted_input.power_bands[0] = formatted_input.power_bands[0].reshape(int(len(formatted_input.power_bands[0])/2),2)
        return formatted_input
    else:
        sys.exit('In proper data load type passed in.')


def calcRCSFFT(rcs_data: pd.DataFrame,rcs_settings: pd.DataFrame,pb: bool=False) -> pd.DataFrame:
    # grab hann window values used by the RCS device
    hann_win = rcs.create_hann_window(rcs_settings.FFT_size[0], percent=100)

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
            data_rcs_unit = rcs.transform_mv_to_rcs(rcs_data_sim[channel].values,rcs_settings.amp_gains[0][count])
            data_rcs_fft, time_pb = rcs.td_to_fft(data_rcs_unit,
                                                rcs_data_sim['time'].values,
                                                rcs_settings.fs_td[0],    # sample rate of the neural data
                                                rcs_settings.FFT_size[0], # fft size set on the RCS
                                                rcs_settings.FFT_interval[0], # update interval set on the RCS
                                                hann_win,
                                                output_in_mv = False)
            
            data_rcs_pb = rcs.fft_to_pb(data_rcs_fft,
                                        rcs_settings.fs_td[0],
                                        rcs_settings.FFT_size[0],
                                        rcs_settings.FFT_bitshift[0],
                                        band_edges_hz = rcs_settings.power_bands[0][count*2,],
                                        input_is_mv = False)
            pb_sample_mask = np.isin(rcs_data_sim.time,time_pb)
            rcs_data_sim.loc[pb_sample_mask,'pb1_'+channel] = data_rcs_pb

            data_rcs_pb = rcs.fft_to_pb(data_rcs_fft,
                                        rcs_settings.fs_td[0],
                                        rcs_settings.FFT_size[0],
                                        rcs_settings.FFT_bitshift[0],
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
            data_rcs_unit = rcs.transform_mv_to_rcs(rcs_data[channel].values,rcs_settings.amp_gains[0][count])
            data_rcs_fft, time_pb = rcs.td_to_fft(data_rcs_unit,
                                                rcs_data['time'].values,
                                                rcs_settings.fs_td[0],    # sample rate of the neural data
                                                rcs_settings.FFT_size[0], # fft size set on the RCS
                                                rcs_settings.FFT_interval[0], # update interval set on the RCS
                                                hann_win,
                                                output_in_mv = False)
            
            data_rcs_power = rcs.fft_to_pb(data_rcs_fft,
                                        rcs_settings.fs_td[0],
                                        rcs_settings.FFT_size[0],
                                        rcs_settings.FFT_bitshift[0],
                                        input_is_mv = False)
            
            key_freq_names = generateKeyFreqNames(rcs_settings.fs_td[0],rcs_settings.FFT_size[0],channel)
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
        elif args[i] == "-sim_type":
            sim_type = args[i+1]
        elif args[i] == "-save_path":
            save_path = args[i+1]

    # debugging prints
    print(args[0])          # should be RCS data path
    print(args[1])          # should be RCS settings path
    print('sucessfullly loaded data')   

    # generate save name
    if sim_type == "fft":
        save_name = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_full_spec.csv"
    else:
        save_name = save_path + "/" + data_file_name[0:len(data_file_name)-4] + "_pb.csv"

    # Check what type of simulation to run and save output
    if sim_type == "fft":
        out_DF = calcRCSFFT(rcs_data,rcs_settings)
        out_DF.to_csv(save_name,index = False)
    elif sim_type == "pb":
        out_DF = calcRCSFFT(rcs_data,rcs_settings,pb = True)
        out_DF.to_csv(save_name,index = False)



""" 
Entry point into the script
"""
if __name__=="__main__":
    main()