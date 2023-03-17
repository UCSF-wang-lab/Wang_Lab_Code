"""
Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/08/2022
*Note, as the script is updated, the author may not be accurate anymore*

This module defines a couple of different ways to calculate the accuracy of a given threshold and other hyperparameter values.

Dependencies:   numpy       [=] Common library and can be obtained using pip
                pandas      [=] Common library and can be obtained using pip
"""

import numpy as np
import pandas as pd

def calcThresholdAccuracySwingPhase(RCS_time,RCS_data,event_timings,threshold_val,flip_vals: bool = False):
    """
    Inputs:     parameter_values[=] The parameter values to use when calculating accuracy. [1,n] size, where n is the number of parameters
                RCS_data        [=] Neural power data that is used to check the parameter values.
                event_timings   [=] The timings of the gait event to determine if the threshold is accurate or not.
                RCS_fs          [=] RCS neural data sampling rate. Optional. Default is set to 500 Hz.

    Outputs:    accuracy        [=] The accuracy of detecting a toe off or heel strike from the neural data. Ranges from 0 to 1
    """

    if flip_vals:
        RCS_data = -RCS_data
        threshold_val = -threshold_val

    # Get number of data points and gait events
    n_data_points = RCS_data.shape[0]
    n_events = event_timings.shape[0]

    # initialize arrays
    detector_state = np.zeros((n_data_points,1))
    correct_state = np.zeros((n_data_points,1))
    correct_swing = np.zeros((n_events,1))

    for i in range(len(RCS_data)):
        # Determine the detector state
        if RCS_data[i] >= threshold_val:
            detector_state[i] = 1
        else:
            detector_state[i] = 0

        # Check to see if the current time point is within the double support
        # period. Then, compare if the state is correct.
        A = RCS_time[i] >= event_timings.iloc[:,0]
        B = RCS_time[i] <= event_timings.iloc[:,1]
        C = np.logical_and(A,B)

        if sum(C) >= 1 and detector_state[i] == 1:
            correct_state[i] = 1;   # detector should be in stim state since it is within the double support period
            swing_ind = np.where(C==True)[0][0]
            correct_swing[swing_ind] = 1
        elif sum(C) >= 1 and detector_state[i] == 0:
            correct_state[i] = 0;   # detector should be in stim state, but it is not
            swing_ind = np.where(C==True)[0][0]
            correct_swing[swing_ind] = 0
        elif sum(C) == 0 and detector_state[i] == 1:
            correct_state[i] = 0;   # detector in stim state, but it shouldn't be because it is not during double support period
        elif sum(C) == 0 and detector_state[i] == 0:
            correct_state[i] = 1;   # detector is not stim state, this is correct
    
    # Determine how many nan's there are. These should not be counted in the accuracy calculations
    D = np.where(event_timings.iloc[:,0].isnull().values==True)[0]
    E = np.where(event_timings.iloc[:,1].isnull().values==True)[0]
    n_nans = len(list(set().union(D,E)))

    # Determine how many double support times would not have been detected anyway
    event_not_detected = np.zeros((event_timings.shape[0],1))
    for j in range(event_timings.shape[0]):
        if not (np.isnan(event_timings.iloc[j,0])) | (np.isnan(event_timings.iloc[j,1])):
            F = RCS_time >= event_timings.iloc[j,0]
            G = RCS_time <= event_timings.iloc[j,1]
            if sum(F&G) == 0:
                event_not_detected[j] = 1

    full_thresh_accuracy = sum(correct_state)/len(correct_state)
    swing_thresh_accuracy = sum(correct_swing)/(len(correct_swing)-n_nans-sum(event_not_detected))

    # weighted accuracy. puts more emphasis that the state should be correct during the double support time period
    # than during other time periods
    overall_thresh_accuracy = 0.5*swing_thresh_accuracy + 0.5*full_thresh_accuracy    

    return overall_thresh_accuracy, full_thresh_accuracy, swing_thresh_accuracy



def calcThresholdAccuracyStancePhase(RCS_time,RCS_data,event_timings,threshold_val,flip_vals: bool = False):
    """
    Inputs:     parameter_values[=] The parameter values to use when calculating accuracy. [1,n] size, where n is the number of parameters
                RCS_data        [=] Neural power data that is used to check the parameter values.
                event_timings   [=] The timings of the gait event to determine if the threshold is accurate or not.
                RCS_fs          [=] RCS neural data sampling rate. Optional. Default is set to 500 Hz.

    Outputs:    accuracy        [=] The accuracy of detecting a toe off or heel strike from the neural data. Ranges from 0 to 1
    """
    return 0



def calcThresholdAccuracyDST(RCS_time,RCS_data,event_timings,threshold_val,flip_vals: bool = False):
    """
    This optmizing function calculates the accuracy of a single threshold to change states during double support time.

    Inputs:     parameter_values[=] The parameter values to use when calculating accuracy. [1,n] size, where n is the number of parameters
                RCS_data        [=] Neural power data that is used to check the parameter values.
                event_timings   [=] The timings of RHS to LTO (left hemisphere data) or LHS to RTO (right hemisphere data)
                                    to determine if the threshold is accurate or not. nx2 
                RCS_fs          [=] RCS neural data sampling rate. Optional. Default is set to 500 Hz.

    Outputs:    accuracy        [=] The accuracy of detecting a toe off or heel strike from the neural data. Ranges from 0 to 1
    """

    if flip_vals:
        RCS_data = -RCS_data
        threshold_val = -threshold_val

    # Get number of data points and gait events
    n_data_points = RCS_data.shape[0]
    n_events = event_timings.shape[0]

    # initialize arrays
    detector_state = np.zeros((n_data_points,1))
    correct_state = np.zeros((n_data_points,1))
    correct_dst = np.zeros((n_events,1))

    for i in range(len(RCS_data)):
        # Determine the detector state
        if RCS_data[i] >= threshold_val:
            detector_state[i] = 1
        else:
            detector_state[i] = 0

        # Check to see if the current time point is within the double support
        # period. Then, compare if the state is correct.
        # temp3 = gait_events[["RHS","LTO"]].to_numpy()
        A = RCS_time[i] >= event_timings.iloc[:,0]
        B = RCS_time[i] <= event_timings.iloc[:,1]
        C = np.logical_and(A,B)


        if sum(C) >= 1 and detector_state[i] == 1:
            correct_state[i] = 1;   # detector should be in stim state since it is within the double support period
            dst_ind = np.where(C==True)[0][0]
            correct_dst[dst_ind] = 1
        elif sum(C) >= 1 and detector_state[i] == 0:
            correct_state[i] = 0;   # detector should be in stim state, but it is not
            dst_ind = np.where(C==True)[0][0]
            correct_dst[dst_ind] = 0
        elif sum(C) == 0 and detector_state[i] == 1:
            correct_state[i] = 0;   # detector in stim state, but it shouldn't be because it is not during double support period
        elif sum(C) == 0 and detector_state[i] == 0:
            correct_state[i] = 1;   # detector is not stim state, this is correct


    # Determine how many nan's there are. These should not be counted in the accuracy calculations
    D = np.where(event_timings.iloc[:,0].isnull().values==True)[0]
    E = np.where(event_timings.iloc[:,1].isnull().values==True)[0]
    n_nans = len(list(set().union(D,E)))

    # Determine how many double support times would not have been detected anyway
    event_not_detected = np.zeros((event_timings.shape[0],1))
    for j in range(event_timings.shape[0]):
        if not (np.isnan(event_timings.iloc[j,0])) | (np.isnan(event_timings.iloc[j,1])):
            F = RCS_time >= event_timings.iloc[j,0]
            G = RCS_time <= event_timings.iloc[j,1]
            if sum(F&G) == 0:
                event_not_detected[j] = 1

    full_thresh_accuracy = sum(correct_state)/len(correct_state)
    dst_thresh_accuracy = sum(correct_dst)/(len(correct_dst)-n_nans-sum(event_not_detected))

    # weighted accuracy. puts more emphasis that the state should be correct during the double support time period
    # than during other time periods
    overall_thresh_accuracy = 0.5*dst_thresh_accuracy + 0.5*full_thresh_accuracy    

    return overall_thresh_accuracy, full_thresh_accuracy, dst_thresh_accuracy