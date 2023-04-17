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