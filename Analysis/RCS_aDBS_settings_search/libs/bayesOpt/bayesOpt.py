"""
Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/14/2023
*Note, as the script is updated, the author may not be accurate anymore*

Bayesian Optimization iteration code to determine the best weight and threshold value to classify
fast aDBS states


Dependencies:   numpy       [=] Common library and can be obtained using pip
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

# import timer for code timing
import time



"""
Module functions
"""
def expectedImprovement(x_samples,curr_best,gp_model,n_params,min_fun: bool = False):
    """
    Expected Improvement Acquisition Function
        Calculates the expected improvement value given the mean and standard deviation of a
        Gaussian process evaluated at all possible samples of the search space. A sign change is added
        because we will use the minimizer function to determine the minimum of the acquisition function.
        When minimizing, this means a negative needs to be added such that we can detect the smaller values
        using the minimizer. When maximizing, no sign change is needed.

    Inputs
        x_samples       [=] Samples to calculate the expected improvement on. Can be [1,n_parameters] or [1...,n_parameters]
        curr_best       [=] "Best" value of the current sampled points. Best is relative if you want the min or max value. 
        gp_model        [=] Gaussian process regressor object
        n_params        [=] Number of parameters(?)
        minFun          [=] Boolean flag to denote if the function should be calculated such that you
                            obtain the correct min or max of the acquisition function.
    
    Returns the entire expected improvement calculation. Does not automatically calculates what point to 
        sample next.
    """

    # Check to see if minimizing or maximizing the function. If minimizing the function,
    # flip the values so that you're "maximizing" the function. Also grabs the current best value.
    if min_fun:
        flip_val = -1
    else:
        flip_val = 1

    # Grab GP mean and std of x_sample points
    x_sample_mean,x_sample_std = gp_model.predict(x_samples.reshape(-1,n_params),return_std = True)


    # Error checking to see if dividing by 0. 
    with np.errstate(divide='ignore'):
        A = flip_val * (x_sample_mean - curr_best)
        B = A/x_sample_std

        cdf_val = norm.cdf(B)
        pdf_val = norm.pdf(B)
        ei_val = A * cdf_val + x_sample_std * pdf_val
        ei_val[x_sample_std == 0.0] == 0.0              # Set areas where the standard deviation is 0 to 0. 

    # flip the value again so in correct sign
    return ei_val



def getNextSamplePoint(acq_fun,gp_model,curr_best,param_combos,min_fun: bool = False):
    """
    This function determines the next point to sample by minimizing the acquisition function.
        The acquisition function does take into account a sign flip such that minimizing the
        acquisition function still produces the maximum value.

    Inputs:     acq_fun         [=] The acquisition function to determine the next parameter to sample. 
                gp_model        [=] A gaussian process model of the points evaluated already.
                curr_best       [=] Current best parameter combo
                param_combos    [=] 
                min_fun         [=] Boolean flag to determine if you are minimizing the funtion (True) or if you're trying to maximize the function (False).

    Outputs:    next_point      [=] The next sample point to evaluate the optimizing function.
    """
    next_point = None
    next_point_val = None
    n_param = param_combos.shape[1]   # determine how many parameters we are search through. Each row is a different bound for a parameter.
        
    # Calculate the expected improvement of the parameter space
    ei_vals = acq_fun(param_combos,curr_best,gp_model,n_param)

    # Locate the index of the maximum expected improvement
    next_point = param_combos[ei_vals.argmax(axis=0)]

    return next_point



# BOItrDF
#   Calculates the Gaussian process and Acquisition function at a certain itration of Bayesian optimization.
#
#   Inputs:
#       X_samples                   [=] Parameter pairs that have been evaluated. Shape is an (n x m) matrix, 
#                                       where "n" is the number of all evaluations, and "m" is the number of 
#                                       parameters being optimized over.
#       Y_samples                   [=] Evaluations of the function being optimized over at the X_sample points. 
#                                       Shape is a (n x 1) matrix.
#       search_bounds               [=] The min and max values of the parameters to optimize over. Shape is an
#                                       (p x 2) matrix, where p is the number of parameters to optimize over.
#                                       The first parameter is the index numbers in the raw data table to look
#                                       through.  
#       search_bounds_resolution    [=] The spacing between the min and max value of the parameters above.
#       n_initial_samples           [=] How many initial samples there were before starting Bayesian optimization.
#       gp_params                   [=] Parameters of the Gaussian process used during the Bayesian optimization.
#       itr_spacing                 [=] The spacing of the iterations. For example: [1000,100,10,1].
#       itr_reps                    [=] How many times to repeat the itration spacing above. Continuing the example,
#                                       itr_reps = [8,5,4,10] and the spacing above will result in the function
#                                       generating 
#   Outputs:
#
def BOItrDF(X_samples,Y_samples,param_combos,n_initial_samples: int = 0,gp_params: dict = None,itr_spacing: list = None, itr_reps: list = None):
    # Get number of samples
    n_samples = X_samples.shape[0]
    
    # Get number of parameters
    n_param = X_samples.shape[1]

    # Get iteration numbers to calculate the Gaussian process and acquisition function
    n_itr_indexes = np.zeros(1000,dtype=int)   # Allocate memory
    count = 0
    idx_value = 0
    if itr_spacing is not None:
        for i in range(len(itr_spacing)):
            if count == 0:
                idx_value = n_initial_samples + itr_spacing[i] - 1
            else:
                idx_value = n_itr_indexes[count-1] + itr_spacing[i]
            spacing = np.arange(idx_value,idx_value+itr_spacing[i]*itr_reps[i],itr_spacing[i])
            n_itr_indexes[count:count + itr_reps[i]] = spacing
            count += itr_reps[i]
    else:
        n_itr_indexes = np.arange(n_initial_samples,n_samples+1,1) # Have to add one as the way we are indexing lower only takes the range from [x1,x2) where x2 is not included.
        count = n_samples-n_initial_samples + 1 # add one so you get the correct number of iterations + the initial random sampling

    n_itr_indexes = n_itr_indexes[0:count]

    # Allocate memory for a matrix output of all the iterations to plot
    GP_mat = np.zeros((count,param_combos.shape[0]+1))
    AF_mat = np.zeros((count,param_combos.shape[0]+1))

    # Create GP model
    if gp_params is not None:
        gp_model = createGP(gp_params = gp_params)
    else:
        gp_model = createGP()

    # Loop through the iterations.
    for i in range(0,count):
        # Timer to see how long it takes to recreate each Gaussian process and acquisition function surface
        print(f"Recreating GP and AF of iteration {n_itr_indexes[i]}...", end = " ")
        recreate_tic = time.perf_counter()

        # Fit the model to current sampled points
        curr_Y = Y_samples[0:n_itr_indexes[i]]
        gp_model.fit(X = X_samples[0:n_itr_indexes[i],:], y = curr_Y)

        # Get current best Y value
        curr_best_Y = max(curr_Y)

        # Predict the values of the parameter space
        mu,std = gp_model.predict(param_combos, return_std=True)

        # Calculate acquisition funtion values
        af_vals = expectedImprovement(param_combos,curr_best_Y,gp_model,n_param,min_fun=False)

        # Save outputs to a matrix
        GP_mat[i,0] = n_itr_indexes[i]
        GP_mat[i,1:param_combos.shape[0]+1] = mu
        AF_mat[i,0] = n_itr_indexes[i]
        AF_mat[i,1:param_combos.shape[0]+1] = af_vals

        # write out how long this recreation took
        recreate_toc = time.perf_counter()
        print(f"took {recreate_toc-recreate_tic} seconds")

    return GP_mat, AF_mat



def createGP(gp_params: dict = None, norm_data: bool = True):
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
        cov_fun = gp.kernels.Matern(**gp_params)
    else:
        # Set Matern covariance function.
        # Use default values. It will optimize the length scale, 
        # but smoothness parameter nu is fixed at 3/2 which is good for our purposes.
        cov_fun = gp.kernels.Matern()
    
    gp_model = gp.GaussianProcessRegressor(kernel = cov_fun,n_restarts_optimizer=0,normalize_y=norm_data,copy_X_train=False)
    return gp_model



def createParamSpace(parameter_ranges: list,parameter_resolutions: list):
    # Get number of parameters
    n_param = len(parameter_ranges)
    param_space = [None]*n_param

    # Create spaces
    for i in range(0,n_param):
        param_space[i] = np.arange(parameter_ranges[i][0],parameter_ranges[i][1],parameter_resolutions[i])

    # Create a [?,n_param] array of evaluation points using the parameter space to be evaluated by the acquisition function
    if n_param == 1:
        param_combos = np.array(param_space).reshape((-1,1))
    elif n_param == 2:
        param_combos = np.array([[p1,p2] for p1 in param_space[0] for p2 in param_space[1]])
    elif n_param == 3:
        param_combos = np.array([[p1,p2,p3] for p1 in param_space[0] for p2 in param_space[1] for p3 in param_space[2]])

    return param_combos


def runBO(cost_function, RCS_data: pd.DataFrame, bo_options:dict):
    """
    This function runs the Bayesian optimization algorithm. The data follows this flow:
        (1) Create a Gaussian Process (GP) with the tested parameters optimizing for,
        (2) Minimize the acquisition function to get the next sample point. (A sign change is added such that if maximizing,
            the minimization still works),
        (3) Evaluate the next sample point on the function being optimized,
        (4) Update the list of sampled parameters and function output,
        (5) Repeat until number of iterations is complete

    Inputs:     cost_function   [=] Function to determine how well the sample point performs
                RCS_data        [=] Table of LFP neural data and associated class. Rows are observations. 1 column is mandatory, "Class". 
                                    At least one column name must have the form keyX_YY_YYYYYYYY, where X can be a number from 0-3 and Y can take any
                                    numeric value. The first two YY can be a single Y if the frequency band is under 10.
                bo_options      [=] Dictionary that describes the search space and GP hyperparameters

    Outputs:
    """
    # Go through variable inputs
    for arg, value in bo_options.items():
        match arg:
            case "parameter_space_range":
                parameter_space_range = value
            case "parameter_space_resolution":
                parameter_space_resolution = value
            case "n_itr":
                n_itr = value
            case "gp_params":
                gp_params = value
            case "initial_samples":
                initial_samples = value
            case "initial_outcomes":
                initial_outcomes = value

    # Create parameter space to run Bayesian optimization
    param_space = createParamSpace(parameter_space_range,parameter_space_resolution)

    # Generate and evaluate initial samples if none have been passed in
    if 'initial_samples' not in locals():
        np.random.seed(123) # For reproducability

        # Generate initial samples ~15% of the number of possible parameter combinations
        # n_initial_samples = int(np.floor(param_space.shape[0]*0.15))
        n_initial_samples = int(np.floor(param_space.shape[0]*0.05))
        initial_sample_ind = np.floor(param_space.shape[0]*np.random.random_sample((n_initial_samples)))
        initial_samples = param_space[initial_sample_ind.astype(int),]

    initial_outcomes = np.zeros(n_initial_samples)
    for i in range(0,n_initial_samples):
        details = bo_options["hemisphere"] + ": "
        for x in bo_options["freq_bands"]:
            details += x + " "
        print(f"Running {details} initial random sample {i+1}/{n_initial_samples}")
        
        # Calculate threshold and power band values
        itr_tic = time.perf_counter()
        threshold = np.round(np.quantile(np.sum(RCS_data[:,1:],axis = 1),param_space[int(initial_sample_ind[i])][0]),2)
        power = np.zeros(RCS_data.shape[0])
        for j in range(1,param_space.shape[1]):
            power += (RCS_data[:,j]*param_space[int(initial_sample_ind[i])][j])
        
        loss = cost_function(threshold,power,RCS_data[:,0])
        initial_outcomes[i] = loss

        itr_toc = time.perf_counter()
        print(f"took {itr_toc-itr_tic} seconds")
    
    # Allocate memory to hold sampled parameters and their outcomes. Populates with the all previously evaluated parametesr
    X = np.zeros((initial_samples.shape[0]+n_itr,initial_samples.shape[1]))
    X[0:initial_samples.shape[0],:] = initial_samples

    Y = np.zeros((initial_samples.shape[0]+n_itr,1))    
    Y[0:initial_samples.shape[0]] = initial_outcomes.reshape(-1,1)

    # Determine current best outcome
    curr_best = np.min(initial_outcomes)

    # Create Gaussian process regression object
    if gp_params is None:
        gp_model = createGP()
    else:
        gp_model = createGP(gp_params)

    # Iterate through Bayesian optimization n_itr times
    for itr in range(n_itr):
        # Output iteration details
        print(f"{details} - Running iteration {itr+1}/{n_itr}")
        itr_tic = time.perf_counter()

        # Fit Gaussian process regressor with current sampled points and outcomes
        gp_model.fit(X[0:initial_samples.shape[0]+itr,:],Y[0:initial_samples.shape[0]+itr])

        # Determine next sample point
        # next_sample_point = getNextSamplePoint_OLD(expectedImprovement_OLD,gp_model,Y,search_bounds)
        next_sample_point = getNextSamplePoint(expectedImprovement,gp_model,curr_best,param_space)

        # Evaluate the next sample point
        threshold = np.round(np.quantile(np.sum(RCS_data[:,1:],axis = 1),param_space[int(initial_sample_ind[i])][0]),2)
        power = np.zeros(RCS_data.shape[0])
        for j in range(1,param_space.shape[1]):
            power += (RCS_data[:,j]*param_space[int(initial_sample_ind[i])][j])
        
        new_outcome = cost_function(threshold,power,RCS_data[:,0])

        # Update the list of evaluated sample points and outcomes
        if "search_type" in bo_options.keys():
            if bo_options["search_type"] == "all_combo":
                next_sample_point[0] = round(next_sample_point[0])
                next_sample_point = next_sample_point.reshape(1,2)
        else:
            next_sample_point[0] = round(next_sample_point[0])
            next_sample_point = next_sample_point.reshape(1,2)

        X[initial_samples.shape[0]+itr,:] = next_sample_point
        Y[initial_samples.shape[0]+itr] = new_outcome

        # If new parameter combo results in a lower loss, update current best outcome
        if new_outcome < curr_best:
            curr_best = new_outcome

        # write out how long this iteration took
        itr_toc = time.perf_counter()
        print(f"took {itr_toc-itr_tic} seconds")


    # Create output
    column_names = [None]*(n_param+3)
    for i in range(n_param):
        column_names[i] = "x%d" % (i+1)

    match bo_options["gait_phase"]:
        case "DST":
            column_names[n_param:n_param+3] = ["Weighted Accuracy"]
        case "Swing":
            column_names[n_param:n_param+3] = ["Weighted Accuracy"]
        case "Stance":    
            column_names[n_param:n_param+3] = ["Weighted Accuracy"]
        case _:
            column_names[n_param:n_param+3] = ["Weighted Accuracy"]

    result_table = pd.DataFrame(np.hstack((X,Y)),columns = column_names)

    # Return sampled values and outputs. 
    # Possibly add in the EI and GP output for each iteration as well.
    return result_table, X, Y