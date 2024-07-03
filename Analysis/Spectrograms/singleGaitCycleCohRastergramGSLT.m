function singleGaitCycleCohRastergramGSLT(condition,unilat,varargin)
% INPUTS: 
%       Required
%               condition[=]    String describing the condition analyzed, 
%                               e.g. OFFmed_ONdbs. It creates appropriate subfolders
%                               to save the plots in.
%
%               unilat[=]       0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral
% 
%       Optional
%               dataPath[=]     The path to the folder containing the
%                               coherence data of the 1st condition.
%
%                   keys[=]     Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%               cohPairs[=]     Pairs to compute coherence between. The
%                               list of pairs will be used across all
%                               files. Default is to do all possible
%                               unilateral pairs.
%
%               savePlot[=]     Boolean option to save the resulting merged plots.
%                               Default is false.
%
%
for i = 1:2:(nargin-2)
    switch varargin{i}
        case 'keys'
            keys = varargin{i+1};
        case 'cohPairs'
            cohPairs = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('cohPairs','var') || isempty(cohPairs)
    cohPairs = nchoosek(keys,2);
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    correct_path = false;
    while ~correct_path
        parent_dir = uigetdir();
        if ismac
            split_parent_dir = strsplit(parent_dir,'/');
        elseif ispc
            split_parent_dir = strsplit(parent_dir,'\');
        else
            error('Platform not supported');
        end
            
        if strcmp(split_parent_dir{end},'Data')
            correct_path = true;
        else
            warning('Please select the folder called "Data"');
        end
    end
    
    plot_save_path = fullfile(parent_dir,'Figures','GSLT','Coherence_figs_adaptive','Rastergrams',condition);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load the data
% Fetch & load coherence data for the first condition
if ~exist('dataPath','var') || isempty(dataPath)
    [~,dataPath] = uigetfile('normCohUD.mat','Select the data file of the 1st condition.');
end
data_cond = load(fullfile(dataPath, 'normCohUD.mat'));
fc_cond = load(fullfile(dataPath, 'allfc.mat'));
toeOffs_cond = load(fullfile(dataPath, 'allToeOffsUD.mat'));

%%  Ensure that the number of toe-offs is the same as the number of step epochs
if unilat==0 || unilat==-1
    num_TO_left = cell2mat(cellfun(@(x) length(x),toeOffs_cond.alltoeOffsUD.Left,'UniformOutput',false));
    assert(all(num_TO_left==num_TO_left(1)));

    num_StepEpochs_left = cell2mat(cellfun(@(x) length(x),data_cond.normCohUD.Left,'UniformOutput',false));
    assert(all(num_StepEpochs_left==num_StepEpochs_left(1)));

    assert(num_TO_left(1)==num_StepEpochs_left(1));
end
if unilat==0 || unilat==1
    num_TO_right = cell2mat(cellfun(@(x) length(x),toeOffs_cond.alltoeOffsUD.Right,'UniformOutput',false));
    assert(all(num_TO_right==num_TO_right(1)));

    num_StepEpochs_right = cell2mat(cellfun(@(x) length(x),data_cond.normCohUD.Right,'UniformOutput',false));
    assert(all(num_StepEpochs_right==num_StepEpochs_right(1)));

    assert(num_TO_right(1)==num_StepEpochs_right(1));
end

%% Define the frequency bands
freq_vec = squeeze(fc_cond.allfc.Left{1}(:,:,1));
freq_inds{1} = find(freq_vec>=2 & freq_vec<4); % delta
freq_inds{2} = find(freq_vec>=4 & freq_vec<8); % theta
freq_inds{3} = find(freq_vec>=8 & freq_vec<=12); % alpha
freq_inds{4} = find(freq_vec>=13 & freq_vec<30); % beta
freq_inds{5} = find(freq_vec>=30 & freq_vec<70); % gamma

% Left
if unilat==0 || unilat==-1
    for i=1:length(data_cond.normCohUD.Left)
        data_cohPair = data_cond.normCohUD.Left{i};
        toeOffs_cohPair = toeOffs_cond.alltoeOffsUD.Left{i};
        epoch_lengths = cell2mat(cellfun(@(x) size(x,2),data_cohPair,'UniformOutput',false));

        %% Sort the step segments
        [epoch_lengths_sorted,sortIdx] = sort(epoch_lengths,'descend');
        max_epoch_len = epoch_lengths_sorted(1);

        data_cohPair_sorted = data_cohPair(sortIdx);
        toeOffs_sorted = toeOffs_cohPair(sortIdx);

        %% Populate the array below with the sorted coherence step epochs
        data_cohPair_sorted_mat_left = nan(length(epoch_lengths_sorted),max_epoch_len,length(freq_inds));
        right_sr = 500;
        tvec_left = (0:max_epoch_len-1).*(1/right_sr);

        %% Calculation and Plotting
        for f = 1:length(freq_inds)
            for epoch=1:length(epoch_lengths_sorted)
                data_cohPair_sorted_epoch =  data_cohPair_sorted{epoch};
                for t=1:size(data_cohPair_sorted_epoch,2)
                    data_cohPair_sorted_mat_left(epoch,t,f) = sum(data_cohPair_sorted_epoch(freq_inds{f},t),1);
                end
            end

            % fig_vec(end+1) = figure();
            figure;
            pcolor(tvec_left,size(data_cohPair_sorted_mat_left,1):-1:1,data_cohPair_sorted_mat_left(:,:,f));
            % title({['Left ',chan_names{i},' ',freq_bands_names{j}]});
            shading interp
            
            hold on;
            a = scatter(toeOffs_sorted,length(epoch_lengths_sorted):-1:1,10,'MarkerEdgeColor','r','MarkerFaceColor','r','DisplayName','RTO');

        end
        
    end
end

% Right
if unilat==0 || unilat==1
    for i=1:length(data_cond.normCohUD.Right)
        data_cohPair = data_cond.normCohUD.Right{i};
        toeOffs_cohPair = toeOffs_cond.alltoeOffs.Right{i};
        epoch_lengths = cell2mat(cellfun(@(x) size(x,2),data_cohPair,'UniformOutput',false));
        
        %% Sort the step segments
        [epoch_lengths_sorted,sortIdx] = sort(epoch_lengths,'descend');
        max_epoch_len = epoch_lengths_sorted(1);

        data_cohPair_sorted = data_cohPair(sortIdx);
        toeOffs_sorted = toeOffs_cohPair(sortIdx);

        %% Plotting
        
    end
end

end