function calcLearningResponses(processed_cirris_fname,sequence_vec,varargin)
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

% hitcount_random = zeros(1,length(sequence_vec));
% 
% if procCirrisData.TargetNumber(end)<243
%     warning('Patient did not make it to the sequenced blocks. No sequenced data available');  
%     hitcount_seq = nan(1,length(sequence_vec));
% else
%     hitcount_seq = zeros(1,length(sequence_vec));
% end
Responses = zeros(1,120);
for i = 1:height(procCirrisData)
    block_num  = floor(procCirrisData.TargetNumber(i)/120)+1;

    if block_num == 5
        Responses(mod(procCirrisData.TargetNumber(i),120)+1) = 1;
    end         
end

%% Run the learning probability estimation code
MaxResponse = 1;
BackgroundProb = .5;

runanalysisv2(Responses, MaxResponse, BackgroundProb);

load resultsindividual.mat;

fontsize1 = 20;
linewidth1 = 1;

t=1:size(p,2)-1; 

figure(1);  clf;

%% Visualize the results
subplot(111);  
h = plot(t, pmid(2:end),'b-'); set(h, 'LineWidth',2);
hold on;
plot(t, p05(2:end),'b', t, p95(2:end), 'b', 'LineWidth', linewidth1);
if(MaxResponse == 1)
  hold on; [y, x] = find(Responses > 0);
  h = plot(x,y,'o'); set(h, 'MarkerFaceColor',[.9 .9 .9]);
  set(h, 'MarkerEdgeColor', 'k' ,'MarkerSize', 4);
  hold on; [y, x] = find(Responses == 0);
  h = plot(x,zeros(size(x)),'o'); set(h, 'MarkerFaceColor', [.9 .9 .9]);
  set(h, 'MarkerEdgeColor', 'k','MarkerSize', 4);
  axis([1 t(end)  0 1.0]);
else
  hold on; plot(t, Responses./MaxResponse,'ko');
  axis([1 t(end)  0 1]);
end
plot([1 t(end)], [BackgroundProb  BackgroundProb ], 'b', 'LineWidth', linewidth1);
title(['Learning trial = ' num2str(cback)],'fontsize',fontsize1);
xlabel('Trial Number','fontsize',fontsize1)
ylabel('Probability of a correct response','fontsize',fontsize1)
set(gca,'tickdir','out','fontsize',fontsize1), box off



% % Save plots
% if savePlot
%     saveas(gcf,fullfile(plot_save_path,strcat('hit_rate_',condStr,'.fig')));
%     saveas(gcf,fullfile(plot_save_path,strcat('hit_rate_',condStr,'.tiff')));
% end
% 
% % Save data
% if saveData
%     save(fullfile(data_save_path,strcat('hitrate_all_',condStr,'.mat')),'hitrate_all');
% end
% 
% end