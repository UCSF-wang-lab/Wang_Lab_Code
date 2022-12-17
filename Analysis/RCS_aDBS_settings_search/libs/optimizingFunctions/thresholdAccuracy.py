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

def calcThresholdAccuracyTO2HS(parameter_values,RCS_data,event_timings,RCS_fs: int = 500):
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
                                    to determine if the threshold is accurate or not.
                RCS_fs          [=] RCS neural data sampling rate. Optional. Default is set to 500 Hz.

    Outputs:    accuracy        [=] The accuracy of detecting a toe off or heel strike from the neural data. Ranges from 0 to 1
    """

    threshold_val = 5

    if flip_vals:
        RCS_data = -RCS_data
        threshold_val = -threshold_val

    # Get number of data points and gait events
    n_data_points = RCS_data.shape[0]
    n_events = event_timings.shape[0]

    # initialize arrays
    detector_state = np.empty((n_data_points,1))
    correct_state = np.empty((n_data_points,1))
    correct_dst = np.empty((n_events,1))

    # Count variable to fill in array that tracks if the state is correct specifically during double support time
    count = 1

    for i in range(len(RCS_data)):
        # Determine the detector state
        if RCS_data[i] >= threshold_val:
            detector_state[i] = 1
        else:
            detector_state[i] = 0

        # Check to see if the current time point is within the double support
        # period. Then, compare if the state is correct.
        # temp3 = gait_events[["RHS","LTO"]].to_numpy()
        A = RCS_time[i] >= event_timings[:,0]
        B = RCS_time[i] <= event_timings[:,1]
        C = np.logical_and(A,B)

        if sum(C) >= 1 and detector_state[i] == 1:
            correct_state[i] = 1;   # detector should be in stim state since it is within the double support period
            correct_dst[count] = 1
            count = count + 1       # increment count variable
        elif sum(C) >= 1 and detector_state(i) == 0:
            correct_state[i] = 0;   # detector should be in stim state, but it is not
            correct_dst[count] = 0
            count = count + 1       # increment count
        elif sum(C) == 0 and detector_state[i] == 1:
            correct_state[i] = 0;   # detector in stim state, but it shouldn't be because it is not during double support period
        elif sum(C) == 0 and detector_state(i) == 0:
            correct_state[i] = 1;   # detector is not stim state, this is correct

    full_thresh_accuracy = sum(correct_state)/len(correct_state)
    dst_thresh_accuracy = sum(correct_dst)/len(correct_dst)

    # weighted accuracy. puts more emphasis that the state should be correct during the double support time period
    # than during other time periods
    overall_thresh_accuracy = 0.7*dst_thresh_accuracy + 0.3*full_thresh_accuracy    

    return overall_thresh_accuracy