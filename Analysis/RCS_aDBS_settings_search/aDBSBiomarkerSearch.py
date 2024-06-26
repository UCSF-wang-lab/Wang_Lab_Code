"""
Original Author: Kenneth Louie, PhD (Doris Wang Lab)
Date:   12/13/2022
*Note, as the script is updated, the author may not be accurate anymore*

aDBSBiomarkerSearch (aDBS = adaptive deep brain stimulation):
    Evaluates all frequency band combos for the optimal weighting and threshold value to classify n states.
    Up to 4 frequency bands can be combined to find the optimal classification. Optimization is performed using
    Bayesian optimization.  

Dependencies:   numpy           [=] Common library and can be obtained using pip
                pandas          [=] Common library and can be obtained using pip
                itertools       [=] Used to create combinations of frequency bands

Inputs:         -data                       [=] Flag to notify the program the next input will be the path to the data file to convert
                -n_power_bands              [=] Optional input to control how many power bands it considers when trying to find the optimal band/threshold combination. 
                -BO_nItr                    [=] Maximum number of iterations to run the Bayesian optimization optimizer
                -parallelize                [=] Boolean flag to indicate whether to run the search parallized. 
                -save_path                  [=] Where to save the tables that are generated
"""

# import libraries and modules
import sys
import os
from ast import literal_eval
import time # code timing / debugging purposes
import numpy as np
import pandas as pd
from libs.costFunctions import thresholdAccuracy as ta # function to evaluate new settings
# from libs.bayesOpt import bayesOpt as bo # Bayesopt class
# from libs.bruteForce import bruteForce as bf # Brute force class

"""
Helper functions
"""
def writeResults(write_queue,save_name,pipe_child):
    # Status is a multiprocessing pipe object. Indicate this process has started.
    if pipe_child is not None:
        pipe_child.send('started')

    while True:
        try:
            # Grab next thing in queue
            obj = write_queue.get_nowait()

            # If it is a string and it say "kill", it means all jobs are done and to break this loop
            if type(obj) is str:
                if obj == "kill":
                    # print('Writer kill message received',flush=True)
                    pipe_child.send('completed')
                elif obj == "close":
                    break
            elif type(obj) is pd.DataFrame:
                obj.to_csv(save_name,mode = "a", header = not os.path.exists(save_name),index = False)
            else:
                pass
        except Exception:
            pass
    return "Writer closed."
        


def runSearch(search_arguments,file_name = None,write_queue = None,pipe_child = None,save_name = None, rerun = False, data = None):
    # Status is a multiprocessing pipe object. Indicate this process has started.
    if pipe_child is not None:
        pipe_child.send('started')

    # Load in data
    if data is None and file_name is not None:
        data = pd.read_csv(file_name)

    # Filter data to included rows
    include_ind = np.where((data.Side == search_arguments["hemisphere"]) & (data.Bitshift == search_arguments["bitshift"]) & (data.FFTInterval == search_arguments["fft_int"]) & (data.OverlapPercent == search_arguments["fft_overlap"]))
    data = data.loc[include_ind]

    # Initialize numpy array
    powerBandSum = np.zeros((data.shape[0],len(search_arguments["param_set"])))
    powerBandSum[:,0] = data.iloc[:,search_arguments["param_set"][0]]

    # Sum up values for each power band
    for x in range(1,len(search_arguments["param_set"])):
        powerBand = data.iloc[:,search_arguments["param_set"][x][0]:search_arguments["param_set"][x][1]+1].sum(axis = 1, numeric_only = True)
        powerBandSum[:,x] = powerBand

    # Create string of frequency bands for labeling
    freq_bands = [None]*(len(search_arguments["param_set"])-1)
    for i in range(1,len(search_arguments["param_set"])):
        freq_bands[i-1] = data.columns[search_arguments["param_set"][i][0]]+"--"+data.columns[search_arguments["param_set"][i][1]]

    # Pull out filter variables
    fft_size = np.unique(data.NFFT)[0]

    # Check which search method is being used
    if "bayesOpt" in search_arguments["search_method"].__module__:
        # Create dictionary of Bayesian optimization parameters
        # gp_options = {"nu":0.5,"length_scale":np.array([1,2.5]),"length_scale_bounds":"fixed"}
        gp_options = {"nu":0.5,"length_scale":np.concatenate((np.array([0.05]),np.repeat(0.5,len(search_arguments["param_set"])-1)))}
        options = {dict_key: search_arguments[dict_key] for dict_key in ("parameter_space_range","parameter_space_resolution","hemisphere")}
        options["gp_params"] = gp_options
        options["n_itr"] = 500
        options["freq_bands"] = freq_bands
    elif "bruteForce" in search_arguments["search_method"].__module__:
        options = {dict_key: search_arguments[dict_key] for dict_key in ("parameter_space_range","parameter_space_resolution","hemisphere")}
        options["freq_bands"] = freq_bands

    out_DF = search_arguments["search_method"](search_arguments["costFunction"],powerBandSum,options)
    
    if out_DF.empty is False:
        threshold_column_ind = np.where(out_DF.columns == "Threshold")[0][0]
        out_DF.insert(int(threshold_column_ind),"FFT_overlap",search_arguments["fft_overlap"])
        out_DF.insert(int(threshold_column_ind),"FFT_size",fft_size)
        out_DF.insert(int(threshold_column_ind),"FFT_Int",search_arguments["fft_int"])
        out_DF.insert(int(threshold_column_ind),"Bitshift",search_arguments["bitshift"])

        # Add data frame to write queue if exist
        if write_queue is not None:
            write_queue.put(out_DF)
        else:
            out_DF.to_csv(save_name,mode = "a", header = not os.path.exists(save_name),index = False)

    # Create detail string   
    details = str(search_arguments["worker_num"]) + " " + search_arguments["hemisphere"] + f" - fft size = {fft_size}, " + "fft interval = " + str(search_arguments["fft_int"]) + ", fft overlap = " + str(search_arguments["fft_overlap"]) + ", bitshift = " + str(search_arguments["bitshift"]) + ": " # for debugging
    # details = search_arguments["hemisphere"] + f" - fft size = {fft_size}, " + "fft interval = " + str(search_arguments["fft_int"]) + ", fft overlap = " + str(search_arguments["fft_overlap"]) + ", bitshift = " + str(search_arguments["bitshift"]) + ": "
    for x in freq_bands:
        details += x + " "
    details = details + "completed."

    # print message of settings that finished running 
    print(details,flush=True)

    # If this function is called using a multiprocessing pool, close the pipe to indicate process is complete.
    if pipe_child is not None:
        pipe_child.send('completed')

    # return message
    return details



"""
Entry point into the script
"""
if __name__=="__main__":
    # grab command line arguments
    args = sys.argv[1:]

    ##### Parses optional inputs after the first three inputs #####
    for i in range(0,len(args),2):
        if args[i] == "-data":
            file_name = args[i+1]
            data_file_name = os.path.basename(file_name)    # grabs the current filename to be used as a filename for saving.
        elif args[i] == "-search_list":
            search_list = pd.read_csv(args[i+1])
        elif args[i] == "-n_features":     # Can be any value 1-4
            n_features = int(args[i+1])
        elif args[i] == "-parameter_space_resolution":
            # number of elements needs to be n_features+1. First element is the threshold resolution as the number of linearly spaced values, 
            # the next elements are the resolution of the weight vector
            parameter_space_resolution = literal_eval(args[i+1])
        elif args[i] == "-parameter_space_range":
            # string should indicate a list of pairs. First element is the range for the threshold represented as a percentage. For example, if you have 2 features and the first input is (-0.80,0.80),
            # the range would end up being -0.80*sum(max(features)),0.8*sum(max(features))
            parameter_space_range = literal_eval(args[i+1])   
        elif args[i] == "-search_method":
            if args[i+1] == "brute-force":
                from libs.bruteForce import bruteForce as bf # Brute force class
                search_method = bf.runBF
            elif args[i+1] == "bayes-opt":
                from libs.bayesOpt import bayesOpt as bo # Bayesopt class
                search_method = bo.runBO
        elif args[i] == "-parallel":   # indicate if the search should be parallized. 
            if args[i+1].lower() == "true":
                parallel = True
            else:
                parallel = False
        elif args[i] == "-chunk_size":
            chunk_size = int(args[i+1])
        elif args[i] == "-starting_chunk_number":
            chunk_start = int(args[i+1])
        elif args[i] == "-save_path":
            save_path = args[i+1]       # Set the save path if given, otherwise will use the current working directory of the code

    ##### Default values #####
    if "data" not in locals():
        file_name = "C:/Users/Wang Lab/Documents/gait_RCS_02_AggregateSpecData_fs500_nfft_1024.csv"
        data_file_name = os.path.basename(file_name)

    if "search_list" not in locals():
        search_list = pd.read_csv("C:/Users/Wang Lab/Documents/gait_RCS_02_AggregateSpecData_fs500_nfft_1024_1features_aDBS_search_list.csv")

    if "parallel" not in locals():     # sets parallized flag to false if not passed in by the user
        parallel = "False"
    
    if "n_features" not in locals():    # set number of features to 1 by default
        n_features = 1

    if "parameter_space_range" not in locals():
        parameter_space_range = [(0.15,0.85)]
        for i in range(0,n_features):
            parameter_space_range.append((1.0,1.01))

    if "parameter_space_resolution" not in locals():
        parameter_space_resolution = [0.05]
        for i in range(0,n_features):
            parameter_space_resolution.append(0.05)

    if "search_method" not in locals():
        from libs.bruteForce import bruteForce as bf # Brute force class
        search_method = bf.runBF

    if "chunk_size" not in locals():
        chunk_size = 500

    if "chunk_start" not in locals():
        chunk_start = 1

    if "save_path" not in locals():
        # save_path = os.getcwd()
        save_path = "C:/Users/Wang Lab/Documents"
        save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + f"_{n_features}features" + "_search_results_v2.csv"   # generate base save name
    else:
        save_name_results = save_path + "/" + data_file_name[0:len(data_file_name)-4] + f"_{n_features}features" + "_search_results_v2.csv"   # generate base save name

    ##### Run through each line in search space file #####
    progress_file_name = "C:/Users/Wang Lab/Documents/gait_RCS_02_nfft_1024_progress.txt"
    if not os.path.exists(progress_file_name):
        progress_file = open(progress_file_name,"w")
        progress_file.write("Starting aDBS biomarker search.\n")
    else:
        progress_file = open(progress_file_name,"a")
        progress_file.write(f"Restarting aDBS biomarker search at chunk: {chunk_start}.\n")
    progress_file.close()
    
    search_list_inds = [*range(0,search_list.shape[0],chunk_size)]
    if search_list_inds[len(search_list_inds)-1] != search_list.shape[0]:
        search_list_inds.append(search_list.shape[0])

    data = pd.read_csv(file_name)

    for i in range(0,len(search_list_inds)-1):
        if i+1 >= chunk_start:
            chunk_tic = time.perf_counter()
            print(f"Starting chunk {i+1}/{len(search_list_inds)-1}: rows {search_list_inds[i]+1} to {search_list_inds[i+1]}")
            if parallel == True:
                # Import multiprocessing library
                if "mp" not in sys.modules:
                    import multiprocessing as mp

                print("Parallel processing starting.",flush=True)

                # Create a manager to queue up writing output to a file and to manage data
                if "manager" not in locals():
                    manager = mp.Manager()
                write_queue = manager.Queue() 
                
                # Create processing pool
                # pool_size = 9           # Use this when running on local iMac (uses ~75% of available theads)
                pool_size = 12          # Use this when running on storage windows computer (uses 75% of available threads)
                # pool_size = 40          # Use this when running on compute server only (50% of available threads)     
                processing_pool = mp.Pool(processes = pool_size, maxtasksperchild=int(np.ceil(chunk_size/pool_size))+5)        

                # Create empty list to hold parallel processing jobs
                workers = []

                # Start writing queue first
                pipe_parent,pipe_child = mp.Pipe()
                processing_pool.apply_async(writeResults,(write_queue,save_name_results,pipe_child))
                worker = {'type': "writer",'pipe': pipe_parent, 'started': False, 'completed': False, 'error_count':0}
                workers.append(worker)
                print("Process starting: result writer.",flush=True)

            # Open search list file and start workers for each row
            print("Starting jobs.",flush=True)
            count = 1
            for j in range(search_list_inds[i],search_list_inds[i+1]):
                search_dict = {"costFunction":ta.calcThresholdAccuracySwingPhase,"search_method": search_method,
                                        "parameter_space_range":parameter_space_range,"parameter_space_resolution":parameter_space_resolution,
                                        "hemisphere":search_list.Hemisphere[j],
                                        "bitshift":search_list.FFTBitshift[j],
                                        "fft_int":search_list.FFTInt[j],
                                        "fft_overlap":search_list.FFTOverlap[j],
                                        "param_set":literal_eval(search_list.FreqBandCols[j]),
                                        "worker_num":count}
                if parallel == True:
                    pipe_parent,pipe_child = mp.Pipe()
                    processing_pool.apply_async(runSearch,(search_dict,file_name,write_queue,pipe_child))
                    worker = {'type': "worker",'pipe': pipe_parent, 'started': False, 'completed': False,'error_count': 0,'rerun': False, 'worker_num': count}
                    workers.append(worker)
                else:
                    runSearch(search_dict,file_name=file_name,save_name=save_name_results)
                
                count+=1

            if parallel:
                # Parse through workers until they have completed
                while True:
                    # Check status of each worker
                    for worker_num,worker in enumerate(workers):
                        if worker.get('type') == "worker":
                            # If the worker has started, continue to check the other workers
                            if worker.get('started') and worker.get('completed'):
                                continue
                            # If the worker has not started or completed, poll to see if it has started or completed during the check
                            pipe_parent = worker.get('pipe')
                            try: 
                                if pipe_parent.poll(0.1):
                                    message = pipe_parent.recv()
                                    if message == "started":
                                        worker['started'] = True
                                    elif message == "completed":
                                        worker['completed'] = True
                                        worker['pipe'].close()
                            except:
                                if (worker['started'] and worker['error_count']<=10) or (not worker['started']):
                                    # Could not poll worker. Worker may have ended prematurely or not started at all.
                                    worker['error_count'] +=1

                    # Check to see if all unfinished jobs has errored out at least 10 times. If true, rerun each job until finished
                    for worker_num,worker in enumerate(workers):
                        if (worker['type'] == "worker") and (not worker['completed']) and (worker['error_count']>10):
                            worker['rerun'] = True
                            worker['completed'] = True

                    # Check how many jobs have started
                    jobs_in_queue = len(workers)-1
                    for worker in workers:
                        if worker.get('completed'):
                            jobs_in_queue-=1
                    
                    # All jobs completed. Send kill signal to writing process
                    if not jobs_in_queue:
                        write_queue.put("kill")
                        break

                # Check to make sure write process child send a completed message
                write_pipe_parent = workers[0].get('pipe')
                while not workers[0]['completed']:
                    try:
                        if write_pipe_parent.poll(0.1):
                            message = write_pipe_parent.recv()
                            if message == "completed":
                                workers[0]['completed'] = True
                                write_queue.put("close")
                                write_queue.close()
                                write_pipe_parent.close()
                                break
                    except:
                        if workers[0]['error_count']>100:
                            # Assume writer completed job and closed by garbage collector before kill signal sent
                            workers[0]['completed'] = True
                            break
                        else:
                            workers[0]['error_count'] += 1

                if workers[0]['completed']:
                    print("All parallel jobs completed.",flush=True)
                    processing_pool.close()
                    print("Processing pool closed.",flush=True)
                    processing_pool.terminate()
                    print("All parallel jobs terminated.",flush=True)
                    processing_pool.join()

            # Check and run workers manually that didn't finish using parallel processing
            print("Checking for workers that need to be rerun.",flush = True)
            for worker in workers:
                if worker["type"] == "worker" and worker["rerun"]:
                    search_dict = {"costFunction":ta.calcThresholdAccuracySwingPhase,"search_method": search_method,
                                "parameter_space_range":parameter_space_range,"parameter_space_resolution":parameter_space_resolution,
                                "hemisphere":search_list.Hemisphere[search_list_inds[i]+(worker['worker_num']-1)],
                                "bitshift":search_list.FFTBitshift[search_list_inds[i]+(worker['worker_num']-1)],
                                "fft_int":search_list.FFTInt[search_list_inds[i]+(worker['worker_num']-1)],
                                "fft_overlap":search_list.FFTOverlap[search_list_inds[i]+(worker['worker_num']-1)],
                                "param_set":literal_eval(search_list.FreqBandCols[search_list_inds[i]+(worker['worker_num']-1)]),
                                "worker_num":worker['worker_num']}
                    print(f"Rerunning worker {worker['worker_num']} manually.",flush=True)
                    runSearch(search_arguments=search_dict,data=data,save_name=save_name_results, rerun = True)
                    print(f"Worker {worker['worker_num']} completed manual run.",flush=True)
            print("Worker reruns complete.",flush=True)
        
        print(f"Chunk {i+1} finished.")
        chunk_toc = time.perf_counter()
        print(f"Chunk ran in {chunk_toc-chunk_tic:0.4f} seconds")    
        progress_file = open("C:/Users/Wang Lab/Documents/gait_RCS_02_nfft_1024_progress.txt","a")
        progress_file.write(f"Chunk {i+1} finished.\n")
        progress_file.close()

        # Pause calculations to let garbage collection happen
        time.sleep(5)