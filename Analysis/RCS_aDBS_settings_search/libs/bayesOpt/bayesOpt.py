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
import numpy as np
import pandas as pd
import sklearn.gaussian_process as gp

# import a specific function to use
from scipy.stats import norm 
from scipy.optimize import minimize



"""
Module functions
"""
def expectedImprovement(x_samples,y_evaluated,gp_model,min_fun: bool = False):
    """
    Expected Improvement Acquisition Function
        Calculates the expected improvement value given the mean and standard deviation of a
        Gaussian process evaluated at all possible samples of the search space. A sign change is added
        because we will use the minimizer function to determine the minimum of the acquisition function.
        When minimizing, this means a negative needs to be added such that we can detect the smaller values
        using the minimizer. When maximizing, no sign change is needed.

    Inputs
        x_samples       [=] Samples to calculate the expected improvement on. Can be [1,n_parameters] or [1...,n_parameters]
        y_evaluated     [=] Values of all evaluated Current evaluated point values
        gp_model        [=] Gaussian process regressor object
        minFun          [=] Boolean flag to denote if the function should be calculated such that you
                            obtain the correct min or max of the acquisition function.
    
    Returns the entire expected improvement calculation. Does not automatically calculates what point to 
        sample next.
    """

    # Check to see if minimizing or maximizing the function. If minimizing the function,
    # flip the values so that you're "maximizing" the function. Also grabs the current best value.
    if min_fun:
        curr_best = np.min(y_evaluated)
        flip_val = -1
    else:
        curr_best = np.max(y_evaluated)
        flip_val = 1

    # Grab GP mean and std of x_sample points
    x_sample_mean,x_sample_std = gp_model.predict(x_samples,return_std = True)


    # Error checking to see if dividing by 0. 
    with np.errstate(divide='ignore'):
        A = flip_val * (x_sample_mean - curr_best)
        B = A/x_sample_std

        cdf_val = norm.cdf(B)
        pdf_val = norm.pdf(B)
        ei_val = A * cdf_val + x_sample_std * pdf_val
        ei_val[x_sample_std == 0.0] == 0.0              # Set areas where the standard deviation is 0 to 0. 

    # flip the value again so in correct sign
    return -1 * ei_val


def getNextSamplePoint(acq_fun,gp_model,y_evaluated,search_bounds,min_fun: bool = False, search_restarts: int = 100):
    """
    This function determines the next point to sample by minimizing the acquisition function.
        The acquisition function does take into account a sign flip such that minimizing the
        acquisition function still produces the maximum value.

    Inputs:

    Outputs:
    """
    next_point = None
    next_point_val = None
    n_params = search_bounds.shape[0]   # determine how many parameters we are search through. Each row is a different bound for a parameter.

    # pick random starting points to minimize the acquisition function
    initial_samples = np.random.uniform(search_bounds[:,0],search_bounds[:,1],size = (search_restarts,n_params))

    # Loop through initial sample points to determine the minimal value of the acquisition function
    for initial_sample in initial_samples:
        curr_results = minimize(fun = acq_fun, x0 = initial_sample,bounds = search_bounds, method='L-BFGS-B', args=(gp_model,y_evaluated,min_fun,n_params))

        # Check to see if the current results ended up with a lower value than the current best
        if next_point is None:
            next_point = curr_results.x
            next_point_val = curr_results.fun
        elif next_point_val > curr_results.fun:
            next_point = curr_results.x
            next_point_val = curr_results.fun

    return next_point


def createGP(gp_params: dict = None):
    """
    This function generates a Gaussian process regressor object. The object
        is created using passed in parameters, such as the length constant and
        the covariance function, but also will set default settings if no parameters
        are passed in.

    Inputs:
        gp_params   [=] Dictionary object. Optional.
    
    Returns a Gaussian process regressor object with default or set parameters
    """
    if gp_params is not None:
        gp_model = gp.GaussianProcessRegressor(**gp_params)
    else:
        # Set Matern covariance function.
        # Use default values. It will optimize the length scale, 
        # but smoothness parameter nu is fixed at 3/2 which is good for our purposes.
        cov_fun = gp.kernels.Matern()
        gp_model= gp.GaussianProcessRegressor(kernel = cov_fun,n_restarts_optimizer=10,normalize_y=True)
    
    return gp_model


def runBO(func_2_optimize: function, search_bounds: np.array, RCS_data: pd.DataFrame, event_timings:pd.DataFrame, n_itr: int = 15, gp_params: dict = None, initial_samples: list = None, initial_outcomes: list = None):
    """
    This function runs the Bayesian optimization algorithm. The data follows this flow:
        (1) Create a GP with the tested parameters optimizing for,
        (2) Minimize the acquisition function to get the next sample point. (A sign change is added such that if maximizing,
            the minimization still works),
        (3) Evaluate the next sample point on the function being optimized,
        (4) Update the list of sampled parameters and function output,
        (5) Repeat until number of iterations is complete

    Inputs:     n_itr               [=] The number of iterations the Bayesian optimization should do (Optional)
                gp_params           [=] Parameters for the Gaussian process regressor object (Optional)
                initial_samples     [=]
                initial_outcomes    [=]

    Outputs:
    """
    # Grab number of parameters of the search space
    n_param = search_bounds.shape[0]

    # Create Gaussian process regression object
    gp_model = createGP()

    # Fit a Gaussian process using the initial samples and outcomes
    if initial_samples is None:
        # If there were no initial samples, randomly select 3 parameters in the search space to evaluate
        initial_samples = np.random.uniform(search_bounds[:,0],search_bounds[:,1],(3,n_param))
        initial_samples[:,1] = initial_samples[:,1].round()

        # Unlikely, but check to make sure none of the random samples are repeats
        # TODO

        # Evaluate the randomly selected parameters
        initial_outcomes = np.empty((3,1))
        RCS_time = RCS_data.time
        for i in range(initial_samples.shape[0]):
            pb_data = RCS_data.iloc[:,initial_samples[i,1]].to_numpy()
            threshold_val = initial_samples[i,0]
            initial_outcomes[i] = func_2_optimize(RCS_time,pb_data,event_timings,threshold_val)

    # Convert list to numpy arrays
    X = np.array(initial_samples)
    Y = np.array(initial_outcomes)

    for itr in range(n_itr):
        # Fit Gaussian process regressor with current sampled points and outcomes
        gp_model.fit(X,Y)

        # Determine next sample point
        next_sample_point = getNextSamplePoint(expectedImprovement,gp_model,Y,search_bounds,True)

        # Evaluate the next sample point
        pb_data = RCS_data.iloc[:,round(next_sample_point[0,1])].to_numpy()
        threshold_val = next_sample_point[0,0]
        new_outcome = func_2_optimize(RCS_time,pb_data,event_timings,threshold_val)

        # Update the list of evaluated sample points and outcomes
        X = np.append(X,next_sample_point)
        Y = np.append(Y,new_outcome)

    # Return sampled values and outputs. 
    # Possibly add in the EI and GP output for each iteration as well.
    return X,Y