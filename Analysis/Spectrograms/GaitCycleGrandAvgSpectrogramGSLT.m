function GaitCycleGrandAvgSpectrogramGSLT(fileList,varargin)
% Calculates the grand average gait cycle average spectrogram from all the
% files passed into this function. Multiple optional inputs allows
% customization of the computation including resolution of the gait cycle,
% and normalization strategy.
%
% INPUTS:  Required
%                   fileList[=] list of strings in a cell array. Each
%                               string is the file path to the .mat file of
%                               the aligned data with gait events.
%          Optional
%               nPercentBins[=] The number of bins to break up the gait
%                               cycle in. Default 100, meaning each bin is
%                               1% of the gait cycle.
%
%             stepStartEvent[=] Which gait event to be the start of the
%                               step. Can be left heel strike (LHS) or 
%                               right heel strike (RHS). Default is LHS.
%
%               stepEndEvent[=] Which gait event to be the end of the
%                               step. For example, if the start event is
%                               LHS then end event is RHS,
%                               respectively. If the start event is RHS,
%                               then end event is Default is LHS,
%                               respecitvely. Due to the structure of the
%                               task, the step can only start and end at
%                               heel strike events.
%
%                     normBy[=] What to normalize the gait cycle average
%                               too. Can be 3 options: none, baseline,
%                               average_during_walking. Default is none.
%
%                   normType[=] How to calculate the normalization. Can be
%                               3 options: none, percent_change, or zscore.
%                               Default is none.
%
%               geRangeTable[=] Limits the number of gait events to
%                               consider for the gait cycle average. This
%                               is because the gait event markings may
%                               include when the patient is walking to a
%                               chair at the end of the trial. Default is
%                               to consider all gait events.
%               
%                       keys[=] Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%                   swapKeys[=] 2x<number of files> cell array that
%                               indicates if keys need to be swapped.
%                               This is so the code can be backward
%                               compatible with Phil's RCS patient's
%                               recording settings.
%
%                   savePlot[=] Boolean option to save the resulting plots.
%                               Default is false.
%
%                   saveData[=] Boolean option to save the spectrum data.
%                               Default is false.
%
%                multpartNum[=] In case the dataset has been split into
%                               multiple parts, give the part number analyzed
%                               
%
%   Example call:
%           RCS03 = '/Users/klouie/Documents/Backup Data/RCS03_OG_after_task_OFF_STIM_w_Gait_Events_Julia.mat';
%           RCS14 = '/Users/klouie/Documents/Backup Data/RCS14_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
%           gRCS01 = '/Users/klouie/Documents/Backup Data/gait_RCS_01_OG_ON_Meds_Trial1_w_Gait_Events_Ken.mat';
%           gRCS02 = '/Users/klouie/Documents/Backup Data/gait_RCS_02_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
%           gRCS03 = '/Users/klouie/Documents/Backup Data/gait_RCS_03_OG_OFF_Stim_ON_Meds_Trial1_w_Gait_Events_Julia.mat';
%           fileList = {RCS03,RCS14,gRCS01,gRCS02,gRCS03};
%           keys = {'key0','key1','key2'};
%           swapKeys = {{'key2','key3'},{'key2','key3'},{},{},{}};
%           GaitCycleGrandAvgSpectrogramGSLT(fileList,'gcStartEvent','LHS','normBy','average_during_walking','normType','zscore','keys',keys,'swapKeys',swapKeys)
%
% Date:     03/11/2024
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Contributor: Eleni Patelaki (eleni.patelaki@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Parse optional arguments
for i = 1:2:nargin-1
    switch varargin{i}
        case 'nPercentBins'
            nPercentBins = varargin{i+1};
        case 'stepStartEvent'
            stepStartEvent = varargin{i+1};
        case 'stepEndEvent'
            stepEndEvent = varargin{i+1};
        case 'normBy'
            normBy = varargin{i+1};
        case 'normType'
            normType = varargin{i+1};
        case 'geRangeTable'
            geRangeTable = varargin{i+1};
        case 'keys'
            keys = varargin{i+1};
        case 'swapKeys'
            swapKeys = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
        case 'saveData'
            saveData = varargin{i+1};
        case 'condition'
            condition = varargin{i+1};
        case 'multpartNum'
            multpartNum = varargin{i+1};
    end
end

%% Determine parent directory
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

%% Set default values if not passed in by user
if ~exist('nPercentBins','var') || isempty(nPercentBins)
    nPercentBins = 100;
end

if ~exist('stepStartEvent','var') || isempty(stepStartEvent)
    stepStartEvent = 'RHS';
end

if ~exist('stepEndEvent','var') || isempty(stepEndEvent)
    stepEndEvent = 'LHS';
end

if ~exist('normBy','var') || isempty(normBy)
    normBy = 'none';
end

if ~exist('normType','var') || isempty(normType)
    normType = 'none';
end

if ~exist('geRangeTable','var') || isempty(geRangeTable)
    geRangeTable = [];
end

if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('swapKeys','var') || isempty('swapKeys')
    swapKeys = [];
end

if ~exist('condition','var')
    condition = '';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    plot_save_path = fullfile(parent_dir,'Figures','GSLT_spectrum_figs_adaptive',condition);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;   % Does not save data by default.
else
    data_save_path = fullfile(parent_dir,'Analysis Data','aDBS','GSLT_spectrum_data_adaptive',condition);
    if ~exist(data_save_path, 'dir')
       mkdir(data_save_path);
    end
end

%% Set up variables for each recording key
% To hold gait cycle data
step.Left = cell(1,length(fileList)*length(keys));
step.Right = cell(1,length(fileList)*length(keys));

% To hold normalization data
normalization.Left = cell(1,length(fileList)*length(keys));
normalization.Right = cell(1,length(fileList)*length(keys));

% To hold output frequency of wavelet ouptut (fc stands for frequency
% cell)
fc.Left = cell(1,length(fileList)*length(keys));
fc.Right = cell(1,length(fileList)*length(keys));

% Hold normalized data
normalized.Left = cell(1,length(fileList)*length(keys));
normalized.Right = cell(1,length(fileList)*length(keys));

% Hold grand averages
grandAverage.Left = cell(1,length(fileList)*length(keys));
grandAverage.Right = cell(1,length(fileList)*length(keys));

% Hold toe off values (in % GC)
% The goal is to identify the start and end point of the pre-swing phase
toeOff.Left = cell(1,length(fileList)*length(keys));
toeOff.Right = cell(1,length(fileList)*length(keys));

%% Go through all files and extract data
for i = 1:length(fileList)
    load(fileList{i});
    gaitEventsSorted = sortGaitEventsGSLT(aligned_data.gait_events,stepStartEvent);
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
    % Check if neural data keys need to be swapped
    if ~isempty(swapKeys) && ~isempty(swapKeys{i})
        signalAnalysisData = swapData(signalAnalysisData,swapKeys{i}{1},swapKeys{i}{2});
    end
    
    if isempty(geRangeTable)
        geRange = [1,height(gaitEventsSorted)];
    elseif isinf(geRangeTable(i,2)) && geRangeTable(i,1) == 1
        geRange = [1,height(gaitEventsSorted)];
    else
        start_ind = find(gaitEventsSorted.(stepStartEvent) > geRangeTable(i,1),1,'first');
        end_ind = find(gaitEventsSorted.(stepEndEvent) < geRangeTable(i,2),1,'last');
        geRange = [start_ind,end_ind];
    end
    
    if isfield(signalAnalysisData,'Left')
        for j = 1:length(keys)
            x = abs(signalAnalysisData.Left.Values{j});
            y = signalAnalysisData.Left.Freq_Values{j};
            t = signalAnalysisData.Left.Time{j};
            
            % Truncate portions of matrices x,y that correspond to
            % frequencues lower than 0.1Hz and higher than 140Hz
            x = x(y>=0.1&y<=140,:);
            y = y(y>=0.1&y<=140,:);
            
            walking_start_ind = find(t >= min(gaitEventsSorted{geRange(1),:})-1,1,'first');
            walking_end_ind = find(t <= max(gaitEventsSorted{geRange(2),:}),1,'last');

            step_mat_left = zeros(length(y),nPercentBins,1);
            count = 1;

            for k = geRange(1):geRange(2)
                if ~isnan(gaitEventsSorted.(stepStartEvent)(k)) && ~isnan(gaitEventsSorted.(stepEndEvent)(k)) && ...
                   ((gaitEventsSorted.(stepEndEvent)(k)-gaitEventsSorted.(stepStartEvent)(k)) < 2) && ...
                   ((gaitEventsSorted.(stepEndEvent)(k)-gaitEventsSorted.(stepStartEvent)(k)) > 0.3) && ...
                   gaitEventsSorted.('AdaptiveRight')(k)==0

                    start_ind = find(abs(t-gaitEventsSorted.(stepStartEvent)(k))<0.001);
                    end_ind = find(abs(t-gaitEventsSorted.(stepEndEvent)(k))<0.001);

                    if (start_ind>end_ind)
                        warning('Something wrong with the indices');
                    end
                    
                    if ~isempty(start_ind) && ~isempty(end_ind) && (start_ind<end_ind)
                        data_snip = x(:,start_ind:end_ind);

                        if sum(isinf(data_snip),'all') == 0
                            percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                            for m = 1:length(percent_inds)-1
                                if m == 1
                                    step_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                                else
                                    step_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                                end
                            end
                            count = count + 1;
                        end

                        % Store the toe-offs
                        if strcmp(stepStartEvent,'LHS')
                            toeOff_ind = find(abs(t-gaitEventsSorted.('RTO')(k))<0.001);
                            if ~isempty(toeOff_ind) && toeOff_ind>start_ind && toeOff_ind<end_ind
                                toeOff.Left{(i-1)*length(keys)+j}(end+1)  = (toeOff_ind-start_ind)/(end_ind-start_ind)*100;
                            end
                        elseif strcmp(stepStartEvent,'RHS')
                            toeOff_ind = find(abs(t-gaitEventsSorted.('LTO')(k))<0.001);
                            if ~isempty(toeOff_ind) && toeOff_ind>start_ind && toeOff_ind<end_ind
                                toeOff.Left{(i-1)*length(keys)+j}(end+1)  = (toeOff_ind-start_ind)/(end_ind-start_ind)*100;
                            end
                        end
                    end
                end
            end
            
            step.Left{(i-1)*length(keys)+j} = step_mat_left;
            if strcmp(normBy,'average_during_walking')
                normalization.Left{(i-1)*length(keys)+j} = [mean(x(:,walking_start_ind:walking_end_ind),2,'omitnan'),...
                        std(x(:,walking_start_ind:walking_end_ind),0,2,'omitnan'),...
                        median(x(:,walking_start_ind:walking_end_ind),2,'omitnan')];
            elseif strcmp(normBy,'baseline')
                normalization.Left{(i-1)*length(keys)+j} = [mean(x(:,start_ind:end_ind),2,'omitnan'),...
                        std(x(:,start_ind:end_ind),0,2,'omitnan'),...
                        median(x(:,start_ind:end_ind),2,'omitnan')];
            end
            fc.Left{(i-1)*length(keys)+j} = y;
            
            clear x y t
        end
    end
    
    if isfield(signalAnalysisData,'Right')
        for j = 1:length(keys)
            x = abs(signalAnalysisData.Right.Values{j});
            y = signalAnalysisData.Right.Freq_Values{j};
            t = signalAnalysisData.Left.Time{j};
            
            % Truncate portions of matrices x,y that correspond to
            % frequencues lower than 0.1Hz and higher than 140Hz
            x = x(y>=0.1&y<=140,:);
            y = y(y>=0.1&y<=140,:);

            walking_start_ind = find(t >= min(gaitEventsSorted{geRange(1),:})-1,1,'first');
            walking_end_ind = find(t <= max(gaitEventsSorted{geRange(2),:}),1,'last');

            step_mat_right = zeros(length(y),nPercentBins,1);
            count = 1;

            for k = geRange(1):geRange(2)
                if ~isnan(gaitEventsSorted.(stepStartEvent)(k)) && ~isnan(gaitEventsSorted.(stepEndEvent)(k)) && ...
                   ((gaitEventsSorted.(stepEndEvent)(k)-gaitEventsSorted.(stepStartEvent)(k)) < 2) && ...
                   ((gaitEventsSorted.(stepEndEvent)(k)-gaitEventsSorted.(stepStartEvent)(k)) > 0.3) && ...
                   gaitEventsSorted.('AdaptiveRight')(k)==0
               
                    start_ind = find(abs(t-gaitEventsSorted.(stepStartEvent)(k))<0.001);
                    end_ind = find(abs(t-gaitEventsSorted.(stepEndEvent)(k))<0.001);
                    
                    if (start_ind>end_ind)
                        warning('Something wrong with the indices');
                    end 
                    
                    if ~isempty(start_ind) && ~isempty(end_ind) && (start_ind<end_ind)
                        data_snip = x(:,start_ind:end_ind);

                        if sum(isinf(data_snip),'all') == 0
                            percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                            for m = 1:length(percent_inds)-1
                                if m == 1
                                    step_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                                else
                                    step_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                                end
                            end
                            count = count + 1;
                        end

                         % Store the toe-offs
                        if strcmp(stepStartEvent,'LHS')
                            toeOff_ind = find(abs(t-gaitEventsSorted.('RTO')(k))<0.001);
                            if ~isempty(toeOff_ind) && toeOff_ind>start_ind && toeOff_ind<end_ind
                                toeOff.Right{(i-1)*length(keys)+j}(end+1)  = (toeOff_ind-start_ind)/(end_ind-start_ind)*100;
                            end
                        elseif strcmp(stepStartEvent,'RHS')
                            toeOff_ind = find(abs(t-gaitEventsSorted.('LTO')(k))<0.001);
                            if ~isempty(toeOff_ind) && toeOff_ind>start_ind && toeOff_ind<end_ind
                                toeOff.Right{(i-1)*length(keys)+j}(end+1)  = (toeOff_ind-start_ind)/(end_ind-start_ind)*100;
                            end
                        end
                    end
                end
            end

            step.Right{(i-1)*length(keys)+j} = step_mat_right;
            if strcmp(normBy,'average_during_walking')
                normalization.Right{(i-1)*length(keys)+j} = [mean(x(:,walking_start_ind:walking_end_ind),2,'omitnan'),...
                        std(x(:,walking_start_ind:walking_end_ind),0,2,'omitnan'),...
                        median(x(:,walking_start_ind:walking_end_ind),2,'omitnan')];
            elseif strcmp(normBy,'baseline')
                normalization.Right{(i-1)*length(keys)+j} = [mean(x(:,start_ind:end_ind),2,'omitnan'),...
                        std(x(:,start_ind:end_ind),0,2,'omitnan'),...
                        median(x(:,start_ind:end_ind),2,'omitnan')];
            end
            fc.Right{(i-1)*length(keys)+j} = y;
            
            clear x y t
        end
    end
end

% Normalize and average all together
count = 1;
for i = 1:length(keys)
    for j = i:length(keys):length(fileList)*length(keys)
        if ~isempty(step.Left{j})
            normalization_mat = nan(size(step.Left{j}));

            for k = 1:size(step.Left{j},3)
                if strcmp(normType,'zscore')
                    normalization_mat(:,:,k) = normalizeData(step.Left{j}(:,:,k),'zscore',normalization.Left{j});
                elseif strcmp(normType,'rm_baseline') 
                    normalization_mat(:,:,k) = normalizeData(step.Left{j}(:,:,k),'rm_baseline',normalization.Left{j});
                elseif strcmp(normType,'percent_change')
                    % TODO
                end
            end
            normalized.Left{count} = cat(3,normalized.Left{count},normalization_mat);
        end
    end
    count = count + 1;
end
    
count = 1;
for i = 1:length(keys)
    for j = i:length(keys):length(fileList)*length(keys)
        if ~isempty(step.Right{j})
            normalization_mat = nan(size(step.Right{j}));

            for k = 1:size(step.Right{j},3)
                if strcmp(normType,'zscore')
                    normalization_mat(:,:,k) = normalizeData(step.Right{j}(:,:,k),'zscore',normalization.Right{j});
                elseif strcmp(normType,'rm_baseline') 
                    normalization_mat(:,:,k) = normalizeData(step.Right{j}(:,:,k),'rm_baseline',normalization.Right{j});
                elseif strcmp(normType,'percent_change') 
                    % TODO
                end
            end
            normalized.Right{count} = cat(3,normalized.Right{count},normalization_mat);
        end
    end
    count = count + 1;
end

% Save spectrum, frequency, toe-off & normalization data
if saveData
    if ~exist('multpartNum','var') || isempty(multpartNum)
        save(fullfile(data_save_path,'step_spec.mat'),'step','-v7.3');
        save(fullfile(data_save_path,'fc_spec.mat'),'fc','-v7.3');
        save(fullfile(data_save_path,'toeoff'),'toeOff','-v7.3');
        save(fullfile(data_save_path,'normalization.mat'),'normalization','-v7.3');
    else
        save(fullfile(data_save_path,strcat('step_spec_part',num2str(multpartNum),'.mat')),'step','-v7.3');
        save(fullfile(data_save_path,strcat('fc_spec_part',num2str(multpartNum),'.mat')),'fc','-v7.3');
        save(fullfile(data_save_path,strcat('toeoff_part',num2str(multpartNum),'.mat')),'toeOff','-v7.3');
        save(fullfile(data_save_path,strcat('normalization_part',num2str(multpartNum),'.mat')),'normalization','-v7.3');
    end
end

% calculate grand averages
for i = 1:length(normalized.Left)
    grandAverage.Left{i} = mean(normalized.Left{i},3);
end

for i = 1:length(normalized.Right)
    grandAverage.Right{i} = mean(normalized.Right{i},3);
end

%% Plotting
for i = 1:length(grandAverage.Left)

    freq_vec = fc.Left{1};

    figure;
    ax = pcolor(1:100,log2(freq_vec),grandAverage.Left{i});
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    xline(mean(toeOff.Left{i}),'LineStyle','--','LineWidth',1.5);
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title({'Left';sprintf('%s Grand Average Spectrum',keys{i})});
    
    % Save plots
    if savePlot    
        if ~exist('multpartNum','var') || isempty(multpartNum)
            saveas(gcf,fullfile(plot_save_path,strcat('Left_',keys{i}),'.fig'));
            saveas(gcf,fullfile(plot_save_path,strcat('Left_',keys{i}),'.tiff'));
        else
            if ~exist(fullfile(plot_save_path,strcat('part_',num2str(multpartNum))), 'dir')
               mkdir(fullfile(plot_save_path,strcat('part_',num2str(multpartNum))));
            end
            saveas(gcf,fullfile(plot_save_path,strcat('part_',num2str(multpartNum)),strcat('Left_',keys{i},'.fig')));
            saveas(gcf,fullfile(plot_save_path,strcat('part_',num2str(multpartNum)),strcat('Left_',keys{i},'.tiff'))); 
        end
    end
end

for i = 1:length(grandAverage.Right)

    freq_vec = fc.Right{1};

    figure;
    ax = pcolor(1:100,log2(freq_vec),grandAverage.Right{i});
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    xline(mean(toeOff.Right{i}),'LineStyle','--','LineWidth',1.5);
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title({'Right';sprintf('%s Grand Average Spectrum',keys{i})});
    
    % Save plots
    if savePlot
        if ~exist('multpartNum','var') || isempty(multpartNum)
            saveas(gcf,fullfile(plot_save_path,strcat('Right_',keys{i}),'.fig'));
            saveas(gcf,fullfile(plot_save_path,strcat('Right_',keys{i}),'.tiff'));
        else
            if ~exist(fullfile(plot_save_path,strcat('part_',num2str(multpartNum))), 'dir')
               mkdir(fullfile(plot_save_path,strcat('part_',num2str(multpartNum))));
            end
            saveas(gcf,fullfile(plot_save_path,strcat('part_',num2str(multpartNum)),strcat('Right_',keys{i},'.fig')));
            saveas(gcf,fullfile(plot_save_path,strcat('part_',num2str(multpartNum)),strcat('Right_',keys{i},'.tiff'))); 
        end
    end
end
end

function newData = swapData(oldData,swapKey1,swapKey2)
newData = oldData;

leftDataExist = sum(contains(fields(oldData),'left_LFP'));
rightDataExist = sum(contains(fields(oldData),'right_LFP'));

if leftDataExist
    temp = oldData.left_LFP_table.(swapKey1);
    newData.left_LFP_table.(swapKey1) = newData.left_LFP_table.(swapKey2);
    newData.left_LFP_table.(swapKey2) = temp;
end

if rightDataExist
    temp = oldData.right_LFP_table.(swapKey1);
    newData.right_LFP_table.(swapKey1) = newData.right_LFP_table.(swapKey2);
    newData.right_LFP_table.(swapKey2) = temp;
end
end