function calcGSLTperformance(processed_cirris_fname,sequence_vec,varargin)
% Calculates and plots GSLT performance across blocks
%
% INPUTS:    
%       Required:
%             processed_cirris_fname[=]     processed Cirris csv file name
%
%                       sequence_vec[=]     boolean vector, with length equal to the
%                                           number of task blocks (random vs sequenced). 
%                                           It contains 0 for random blocks and 1 for 
%                                           sequenced blocks.  
%
%       Optional:
%                           savePlot[=]     Boolean option to save the resulting plots.
%                                           Default is false.
%
%                           saveData[=]     Boolean option to save the resulting data.
%                                           Default is false.
%
%                           savePlot[=]     String describing the current
%                                           task condition.
%
% OUTPUTS:   NONE
%
% Example call: calcGSLTperformance([],[0,0,1,1,1,0,1],'savePlot',1,'condStr','ONmed_OFFdbs');
%
% Author:   Eleni Patelaki
% Date:     11/3/2023

%% Parse optional arguments
for i = 1:2:nargin-2
    switch varargin{i}
        case 'savePlot'
            savePlot = varargin{i+1};
        case 'saveData'
            saveData = varargin{i+1};
        case 'condStr'
            condStr = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('condStr','var') || isempty(condStr)
    condStr = '';
end

% Determine parent directory
correct_path = false;
while ~correct_path
    parent_dir = uigetdir([],'Select the "Data" folder');
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

% If savePlot=1, define the plot_save_path
if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save the plot by default.
else
    % plot_save_path = fullfile(parent_dir,'Analysis Data','GSLT','Coherence_data_adaptive','Behavioral_performance');
    plot_save_path = fullfile(parent_dir,'Figures','GSLT','Behavioral_performance');
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

% If saveData=1, define the plot_data_path
if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;   % Does not save the data by default.
else
    data_save_path = fullfile(parent_dir,'Analysis Data','GSLT','Behavioral_performance');
    if ~exist(data_save_path, 'dir')
       mkdir(data_save_path);
    end
end


% If the name of the proccesed Cirris csv hasn't been provided, get it through the UI option.
if ~exist('processed_cirris_fname','var')||isempty(processed_cirris_fname)
    [proc_cirris_fname,proc_cirris_path] = uigetfile('g*_processed_targets_aligned.csv');
    procCirrisData = readtable(fullfile(proc_cirris_path,proc_cirris_fname));
end

% If the name of the sequence csv hasn't been provided, get it through the UI option.
if ~exist('sequence_vec','var')||isempty(sequence_vec)
   sequence_vec = [0,0,1,1,1,0,1];
end

%% Calculate
% First, check if we have any sequenced block data in the first place.
% If not, then calculate the hit rate for sequenced vs random.
hitcount_random = zeros(1,length(sequence_vec));

if procCirrisData.TargetNumber(end)<243
    warning('Patient did not make it to the sequenced blocks. No sequenced data available');  
    hitcount_seq = nan(1,length(sequence_vec));
else
    hitcount_seq = zeros(1,length(sequence_vec));
end

for i = 1:height(procCirrisData)
    block_num  = floor(procCirrisData.TargetNumber(i)/120)+1;
    if sequence_vec(block_num)
        hitcount_seq(block_num) = hitcount_seq(block_num) + 1;
    else
        hitcount_random(block_num) = hitcount_random(block_num) + 1;
    end         
end

% Calculate the target counts
target_counts = zeros(1,length(sequence_vec));
target_counts(1:floor((procCirrisData.TargetNumber(end)+1)/120)) = 120;
if (floor((procCirrisData.TargetNumber(end)+1)/120))<length(sequence_vec)
    target_count_interrupted_block = mod(procCirrisData.TargetNumber(end),120)+1;
    if target_count_interrupted_block>12 % 10% of the target count of a full block
        target_counts(floor((procCirrisData.TargetNumber(end)+1)/120)+1) = target_count_interrupted_block;
    end
end

% Calculate the hit rates
hitrate_random = hitcount_random./target_counts;
hitrate_seq = hitcount_seq./target_counts;

%% Test for differences in the hit rates between R and S blocks

 % Observed data
n1 = sum(hitcount_seq);
N1 = sum(target_counts(sequence_vec==1));
n2 = sum(hitcount_random);
N2 = sum(target_counts(sequence_vec==0));

% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2);
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
[h,p,stats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2);


%% Test for the learning effect, by testing differences in the hit rates between block 5 (S) and 6 (R), and block 6 and 7 (S)
if ~isnan(hitcount_seq(5)) && hitcount_seq(5)~=0 && ~isnan(hitcount_random(6)) && hitcount_random(6)~=0
    n1 = hitcount_seq(7)+hitcount_seq(5);
    N1 = target_counts(7)+target_counts(5);
    n2 = hitcount_random(6);
    N2 = target_counts(6);

    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    % Chi-square test, by hand
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    [h,p,stats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2);
end

%% Visualize the results
figure;
if procCirrisData.TargetNumber(end)<243
    h_rand = bar(hitrate_random);
    h_rand.FaceColor = 'flat';
else
    hitrate_all = hitrate_random + hitrate_seq;
    h_all = bar(hitrate_all);
    h_all.FaceColor = 'flat';
    h_all.CData(sequence_vec==1,:) = repmat([0.7 0 0],length(find(sequence_vec==1)),1);
end

string_sequence_vec = num2cell(sequence_vec);
string_sequence_vec(sequence_vec==1)={'S'};
string_sequence_vec(sequence_vec==0)={'R'};

set(gca,'XtickLabels',string_sequence_vec);
set(gca,'FontSize',15);

% Save plots
if savePlot
    saveas(gcf,fullfile(plot_save_path,strcat('hit_rate_',condStr,'.fig')));
    saveas(gcf,fullfile(plot_save_path,strcat('hit_rate_',condStr,'.tiff')));
end

% Save data
if saveData
    save(fullfile(data_save_path,strcat('hitrate_all_',condStr,'.mat')),'hitrate_all');
end

end