function calcGSLTperformance(processed_cirris_fname,sequence_vec)
% Calculates and plots GSLT performance across blocks
%
% INPUTS:    
%       Required:
%           processed_cirris_fname  [=]     processed Cirris csv file name
%
%           sequence_vec            [=]     boolean vector, with length equal to the
%                                           number of task blocks (random vs learning). 
%                                           It contains 0 for random blocks and 1 for 
%                                           learning blocks.                 
% OUTPUTS:   NONE
%
% Example call: calcGSLTperformance([],[0,0,1,1,1,0,1]);
%
% Author:   Eleni Patelaki
% Date:     11/3/2023

% If the name of the proccesed Cirris csv hasn't been provided, get it through the UI option.
if ~exist('processed_cirris_fname','var')||isempty(processed_cirris_fname)
    [proc_cirris_fname,proc_cirris_path] = uigetfile('g*_processed_targets_aligned.csv');
    procCirrisData = readtable(fullfile(proc_cirris_path,proc_cirris_fname));
end

% If the name of the sequence csv hasn't been provided, get it through the UI option.
if ~exist('sequence_vec','var')||isempty(sequence_vec)
   sequence_vec = [0,0,1,1,1,0,1];
end

% First, check if we have any learning data in the first place.
% If not, then calculate the hit rate for learning vs random.
hitcount_random = zeros(1,length(sequence_vec));

if procCirrisData.TargetNumber(end)<243
    warning('Patient did not make it to the learning blocks. No learning data available');  
    hitcount_learn = nan(1,length(sequence_vec));
else
    hitcount_learn = zeros(1,length(sequence_vec));
end

for i = 1:height(procCirrisData)
    block_num  = floor(procCirrisData.TargetNumber(i)/120)+1;
    if sequence_vec(block_num)
        hitcount_learn(block_num) = hitcount_learn(block_num) + 1;
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
hitrate_learn = hitcount_learn./target_counts;

% Visualize the results
figure;
if procCirrisData.TargetNumber(end)<243
    h_rand = bar(hitrate_random);
    h_rand.FaceColor = 'flat';
else
    hitrate_all = hitrate_random + hitrate_learn;
    h_all = bar(hitrate_all);
    h_all.FaceColor = 'flat';
    h_all.CData(sequence_vec==1,:) = repmat([0.7 0 0],length(find(sequence_vec==1)),1);
end

string_sequence_vec = num2cell(sequence_vec);
string_sequence_vec(sequence_vec==1)={'S'};
string_sequence_vec(sequence_vec==0)={'R'};

set(gca,'XtickLabels',string_sequence_vec);
set(gca,'FontSize',15);
end