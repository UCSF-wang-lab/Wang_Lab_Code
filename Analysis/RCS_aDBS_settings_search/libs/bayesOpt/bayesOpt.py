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

# import timer for code timing
import time



"""
Module functions
"""
def expectedImprovement_OLD(x_samples,y_evaluated,gp_model,n_params,min_fun: bool = False):
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
    return -1 * ei_val


def getNextSamplePoint_OLD(acq_fun,gp_model,y_evaluated,search_bounds,min_fun: bool = False, search_restarts: int = 100):
    """
    This function determines the next point to sample by minimizing the acquisition function.
        The acquisition function does take into account a sign flip such that minimizing the
        acquisition function still produces the maximum value.

    Inputs:     acq_fun         [=] The acquisition function to determine the next parameter to sample. 
                gp_model        [=] A gaussian process model of the points evaluated already.
                y_evaluated     [=] The acquisition function value at the evaluated sample points. 
                search_bounds   [=] The search range to consider the next sample points. Should be the same throughout the optimization.
                min_fun         [=] Boolean flag to determine if you are minimizing the funtion (True) or if you're trying to maximize the function (False).
                search_restarts [=] How many random samples points to evaluate in the acquisition function before deciding on a point.

    Outputs:    next_point      [=] The next sample point to evaluate the optimizing function.
    """
    next_point = None
    next_point_val = None
    n_params = search_bounds.shape[0]   # determine how many parameters we are search through. Each row is a different bound for a parameter.

    # pick random starting points to minimize the acquisition function
    initial_samples = np.random.uniform(search_bounds[:,0],search_bounds[:,1],size = (search_restarts,n_params))

    # Loop through initial sample points to determine the minimal value of the acquisition function
    for initial_sample in initial_samples:
        initial_sample = initial_sample.reshape(1,-1)   # This is needed because this is how the gp code references a single sample
        curr_results = minimize(fun = acq_fun, x0 = initial_sample,bounds = search_bounds, method='L-BFGS-B', args=(y_evaluated,gp_model,n_params,min_fun))

        # Check to see if the current results ended up with a lower value than the current best
        if next_point is None:
            next_point = curr_results.x
            next_point_val = curr_results.fun
        elif next_point_val > curr_results.fun:
            next_point = curr_results.x
            next_point_val = curr_results.fun

    return next_point



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
                search_bounds   [=] The search range to consider the next sample points. Should be the same throughout the optimization.
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
        n_itr_indexes = np.arange(n_initial_samples-1,n_samples,1) # subtract one so you get the initial random sampling of the space
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
        print(f"Recreating GP and AF of iteration {n_itr_indexes[i] + 1}...", end = " ")
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



def createParamSpace(parameter_range: np.array = None,spacing: str = "linear",n_spacings: np.array = None):
    # Get number of parameters
    n_param = parameter_range.shape[0]
    param_space = [None] * n_param  # Allocate memory for an empty list

    # Check the n_spacings variable. If single value, will apply it to all parameters, 
    # otherwise it should have the same length as the number of parameters.
    if n_spacings is None or (isinstance(n_spacings,int)) or (len(n_spacings)<n_param & len(n_spacings == 1)):
        n_spacings = np.repeat(n_spacings,n_param)
    
    # Check the spacing variable. If single value, will apply it to all parameters, 
    # otherwise it should have the same length as the number of parameters.
    if spacing is None or ((len(spacing)<n_param) & (len(spacing) == 1)):
        spacing = np.repeat(spacing,n_param)

    for i in range(n_param):
        if spacing[i] == "linear":
            param_space[i] = np.linspace(parameter_range[i,0],parameter_range[i,1],n_spacings[i])
        elif spacing[i] == "log":
            log_transform_limit = np.log10(parameter_range[i,:]).round(2)
            param_space[i] = np.logspace(log_transform_limit[0],log_transform_limit[1],num=n_spacings[i]).round(2)

    # if search_bounds is None:
    #     search_bounds = np.zeros((2,2))
        
    #     # set first limit to all columns. First column is time so ignore.
    #     search_bounds[0,0] = 1
    #     search_bounds[0,1] = RCS_data.shape[1]

    #     # set second limit to the 5th and 95th quartile. search bounds should cover ~90% of the data
    #     data_quantiles = np.quantile(RCS_data.values[:,1:RCS_data.shape[1]],[0.05,0.95])
    #     if data_quantiles[0] == 0:
    #         data_quantiles[0] = 1e-3
    #         search_bounds[1,:] = data_quantiles
    #     else:
    #         search_bounds[1,:] = data_quantiles

    #     search_bounds_resolution = np.zeros(2)
    #     search_bounds_resolution[0] = 1
    #     search_bounds_resolution[1] = 5

    # Create a [?,n_param] array of evaluation points using the parameter space to be evaluated by the acquisition function
    if n_param == 1:
        param_combos = np.array(param_space).reshape((-1,1))
    elif n_param == 2:
        param_combos = np.array([[p1,p2] for p1 in param_space[0] for p2 in param_space[1]])
    elif n_param == 3:
        param_combos = np.array([[p1,p2,p3] for p1 in param_space[0] for p2 in param_space[1] for p3 in param_space[2]])

    return param_combos


def runBO(func_2_optimize, RCS_data: pd.DataFrame, event_timings:pd.DataFrame,bo_options:dict):
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
                initial_samples     [=] Initial evaluations of the function to optimize to cut down on repeats. 
                                        nxd where n is the number of samples evaluated and d is the dimensions to optimize for. 
                initial_outcomes    [=] The accuracy value of the optimizing function evaluated at the initial samples above. 

    Outputs:
    """
    # Go through variable inputs
    for arg, value in bo_options.items():
        match arg:
            case "search_bounds":
                search_bounds = value
            case "search_bounds_resolution":
                search_bounds_resolution = value
            case "param_combos":
                param_combos = value
            case "event_strings":
                event_strings = value   # FUTURE CAPABILITY, RIGHT NOW ONLY DOUBLE SUPPORT TIME
            case "n_itr":
                n_itr = value
            case "gp_params":
                gp_params = value
            case "initial_samples":
                initial_samples = value
            case "initial_outcomes":
                initial_outcomes = value
            case "min_fun":
                min_fun = value

    # Set default values for variables needed to run Bayesian optimization
    # Grab number of parameters of the search space
    if param_combos is not None:
        n_param = param_combos.shape[1]
    else:
        n_param = search_bounds.shape[0]

    # Generate and evaluate initial samples if none have been passed in
    if 'initial_samples' not in locals():
        np.random.seed(123) # For reproducability

        if "search_type" in bo_options.keys():
            if bo_options["search_type"] == "all_combo":
                initial_samples_X = np.arange(1,RCS_data.shape[1])
                initial_samples_X = initial_samples_X.reshape(RCS_data.shape[1]-1,1)
            elif bo_options["search_type"] == "specific_freq":
                initial_samples_X = np.repeat(1,round(param_combos.shape[0]*0.1))
                initial_samples_X = initial_samples_X.reshape(-1,1)
        else:
            initial_samples_X = np.arange(1,RCS_data.shape[1])
            initial_samples_X = initial_samples_X.reshape(RCS_data.shape[1]-1,1)

        if param_combos is None:
            initial_samples_Y = np.random.uniform(search_bounds[1,0],search_bounds[1,1],(initial_samples_X.shape[0],1))
            initial_samples_Y = initial_samples_Y.round()
        else:
            if "search_type" in bo_options.keys():
                if bo_options["search_type"] == "all_combo":
                    threshold_unique = np.unique(param_combos[:,1])
                    initial_samples_Y = np.random.choice(threshold_unique,initial_samples_X.shape[0])
                    initial_samples_Y = initial_samples_Y.reshape(RCS_data.shape[1]-1,1)
                elif bo_options["search_type"] == "specific_freq":
                    initial_samples_Y = np.random.choice(param_combos.reshape((param_combos.shape[0],)),initial_samples_X.shape[0])
                    initial_samples_Y = initial_samples_Y.reshape(initial_samples_X.shape[0],1)
            else:
                threshold_unique = np.unique(param_combos[:,1])
                initial_samples_Y = np.random.choice(threshold_unique,initial_samples_X.shape[0])
                initial_samples_Y = initial_samples_Y.reshape(RCS_data.shape[1]-1,1)

        if "search_type" in bo_options.keys():
            if bo_options["search_type"] == "all_combo":
                initial_samples = np.hstack((initial_samples_X,initial_samples_Y))
            elif bo_options["search_type"] == "specific_freq":
                initial_samples = initial_samples_Y

        else:
            initial_samples = np.hstack((initial_samples_X,initial_samples_Y))

    initial_outcomes = np.zeros((initial_samples.shape[0],1))
    initial_dst_outcomes = np.zeros((initial_samples.shape[0],1))
    initial_full_outcomes = np.zeros((initial_samples.shape[0],1))
    RCS_time = RCS_data.time
    for i in range(initial_samples.shape[0]):
        if "key" in bo_options.keys():
            channel = bo_options["key"]
            print(f"Running {channel} initial random sample {i+1}/{initial_samples.shape[0]}")
        else:
            print(f"Running initial random sample {i+1}/{initial_samples.shape[0]}...", end = " ")
            itr_tic = time.perf_counter()

        if "search_type" in bo_options.keys():
            if bo_options["search_type"] == "all_combo":
                pb_data = RCS_data.iloc[:,param_combos[i,0].astype(int)].to_numpy()
                threshold_val = param_combos[i,1]
            elif bo_options["search_type"] == "specific_freq":
                pb_data = RCS_data.iloc[:,1].to_numpy()
                threshold_val = initial_samples[i]
        else:
            pb_data = RCS_data.iloc[:,initial_samples[i,0].astype(int)].to_numpy()
            threshold_val = initial_samples[i,1]
        
        weighted_accuracy,dst_accuracy,full_accuracy = func_2_optimize(RCS_time,pb_data,event_timings,threshold_val)

        initial_outcomes[i] = weighted_accuracy
        initial_dst_outcomes[i] = dst_accuracy
        initial_full_outcomes[i] = full_accuracy

        if "key" not in bo_options.keys():
            itr_toc = time.perf_counter()
            print(f"took {itr_toc-itr_tic} seconds")
    
    # Allocate memory to hold sampled parameters and their outcomes. Populates with the all previously evaluated parametesr
    X = np.zeros((initial_samples.shape[0]+n_itr,n_param))
    Y = np.zeros((initial_samples.shape[0]+n_itr,1))
    Y_dst = np.zeros((initial_samples.shape[0]+n_itr,1))
    Y_full = np.zeros((initial_samples.shape[0]+n_itr,1))
    
    X[0:initial_samples.shape[0],:] = initial_samples
    Y[0:initial_samples.shape[0]] = initial_outcomes
    Y_dst[0:initial_samples.shape[0]] = initial_dst_outcomes
    Y_full[0:initial_samples.shape[0]] = initial_full_outcomes

    # Determine current best outcome
    if min_fun:
        curr_best = np.min(initial_outcomes)
    else:
        curr_best = np.max(initial_outcomes)

    # Create Gaussian process regression object
    if gp_params is None:
        gp_model = createGP()
    else:
        gp_model = createGP(gp_params)

    # Calculate all possible parameter combinations in a linear space. 
    if param_combos is None:
        param_combos = createParamSpace(RCS_data,type="log")

    for itr in range(n_itr):
        # Timer to see how long it takes to run each iteration
        if "key" in bo_options.keys():
            channel = bo_options["key"]
            print(f"{channel} - Running iteration {itr+1}/{n_itr}")
        else:
            print(f"Running iteration {itr}...", end = " ")
            itr_tic = time.perf_counter()

        # Fit Gaussian process regressor with current sampled points and outcomes
        gp_model.fit(X[0:initial_samples.shape[0]+itr,:],Y[0:initial_samples.shape[0]+itr])

        # Determine next sample point
        # next_sample_point = getNextSamplePoint_OLD(expectedImprovement_OLD,gp_model,Y,search_bounds)
        next_sample_point = getNextSamplePoint(expectedImprovement,gp_model,curr_best,param_combos,min_fun)

        # Evaluate the next sample point
        # Can use a different type of rounding because this is not a numpy array...
        if "search_type" in bo_options.keys():
            if bo_options["search_type"] == "all_combo":
                pb_data = RCS_data.iloc[:,round(next_sample_point[0])].to_numpy()
                threshold_val = next_sample_point[1]
            elif bo_options["search_type"] == "specific_freq":
                pb_data = RCS_data.iloc[:,1].to_numpy()
                threshold_val = next_sample_point[0]
        else:
            pb_data = RCS_data.iloc[:,round(next_sample_point[0])].to_numpy()
            threshold_val = next_sample_point[1]
        
        new_outcome,new_full,new_dst = func_2_optimize(RCS_time,pb_data,event_timings,threshold_val)

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
        Y_dst[initial_samples.shape[0]+itr] = new_dst
        Y_full[initial_samples.shape[0]+itr] = new_full
        # X = np.append(X,next_sample_point,axis=0)
        # Y = np.append(Y,new_outcome.reshape(1,1),axis=0)

        if min_fun:
            if new_outcome < curr_best:
                curr_best = new_outcome
        else:
            if new_outcome > curr_best:
                curr_best = new_outcome

        # write out how long this iteration took
        if "key" not in bo_options.keys():
            itr_toc = time.perf_counter()
            print(f"took {itr_toc-itr_tic} seconds")


    # Create output
    column_names = [None]*(n_param+3)
    for i in range(n_param):
        column_names[i] = "x%d" % (i+1)

    match bo_options["gait_phase"]:
        case "DST":
            column_names[n_param:n_param+3] = ["Weighted Accuracy","DST Accuracy","Full Accuracy"]
        case "Swing":
            column_names[n_param:n_param+3] = ["Weighted Accuracy","Swing Accuracy","Full Accuracy"]
        case "Stance":    
            column_names[n_param:n_param+3] = ["Weighted Accuracy","Stance Accuracy","Full Accuracy"]
        case _:
            column_names[n_param:n_param+3] = ["Weighted Accuracy","Gait Phase Accuracy","Full Accuracy"]

    result_table = pd.DataFrame(np.hstack((X,Y,Y_dst,Y_full)),columns = column_names)

    # Return sampled values and outputs. 
    # Possibly add in the EI and GP output for each iteration as well.
    return result_table, X, Y, Y_dst, Y_full