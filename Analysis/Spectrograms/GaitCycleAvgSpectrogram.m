function GaitCycleAvgSpectrogram(aligned_data,signalAnalysisData,varargin)
%% GaitCycleAvgSpectrogram
% Calculates the average gait cycle spectrogram for the aligned data passed
% in. The average gait cycle can be normalized so that it is comparable
% across patients. 
%
% INPUTS:  Required
%               aligned_data        [=] Struct containing all the aligned
%                                       data from the trial of interests. 
%
%               signalAnalysisData  [=] Processed aligned data above.
%                                       Can be either processed using the
%                                       short time Fourier transform
%                                       (calcRCS_STFT function) or
%                                       continuous wavelet transformation
%                                       (calcRCS_CWT function).
%
%          Optional
%               cycleStartEvent     [=] Gait event that will be set as the
%                                       beginning of the gait cycle. Can be
%                                       left heel strike ('LHS'), right
%                                       heel strike ('RHS'), left toe off
%                                       ('LTO'), and right toe off ('RTO').
%                                       Default is left heel strike. 
%
%               nPercentbins        [=] Number of bins to break up the gait
%                                       cycle into. Default is 100 bins.
%
%               normalizeBy         [=] How to normalize the average gait
%                                       cycle data. Can normalize using
%                                       baseline data (baseline), before
%                                       they are walking, or the average
%                                       during the entire walking period
%                                       (average_during walking). You do
%                                       not need to normalize the gait
%                                       cycle ('none'). Default is none.
%
%               normalizationType   [=] How to compute the normalization.
%                                       Can be the percent change
%                                       ('percent_change'), or the z-score
%                                       ('zscore'). You do not need to
%                                       perform a normalization ('none').
%                                       Default is none.
%
%               subjectID           [=] String variable of the name of the
%                                       subject being analyzed.
%
%               savePlot            [=] Boolean option to save the 
%                                       resulting plot. Default is false.
%
%   Example call:
%           load(<filename>)
%           B = calcRCS_CWT(aligned_data,[],1,0.9,[]);
%           GaitCycleAvgSpectrogram(aligned_data,B,'normalizeBy','average_during_walking','normalizationType','zscore');
%
% Date:     05/25/2022
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Project:  MJFF aDBS Gait
% cycle_start_event,n_percent_bins,baseline_data,subjectID,save_flag)

for i = 1:2:nargin-2
    switch varargin{i}
        case 'cycleStartEvent'
            cycleStartEvent = varargin{i+1};
        case 'nPercentBins'
            nPercentBins = varargin{i+1};
        case 'normalizeBy'
            normalizeBy = varargin{i+1};   % Can be none|baseline|average_during_walking
        case 'normalizationType'
            normalizationType = varargin{i+1}; % Can be none|percent_change|zscore
        case 'baselineTime'
            baselineTime = varargin{i+1};  % Only valid if "normalization_by" is set to "baseline"
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

if ~exist('cycleStartEvent','var')
    cycleStartEvent = 'LHS';
end

if ~exist('nPercentBins','var')
    nPercentBins = 100;
end

if ~exist('baselineTime','var')
    baselineTime = [];
end

if ~exist('normalizeBy','var')
    normalizeBy = 'none';
end

if ~exist('normalizationType','var')
    normalizationType = 'none';
end

if strcmp(normalizeBy,'none') && ~strcmp(normalizationType,'none')
    warning('No input for normalize_by. Defauling to average during walking.');
    normalizeBy = 'average_during_walking';
end

if ~strcmp(normalizeBy,'none') && strcmp(normalizationType,'none')
    warning('No input for normalization_type. Defauling to zscore.');
    normalizationType = 'zscore';
end

if strcmp(normalizeBy,'baseline') && ~exist('baselineTime','var')
    error('Normalization by baseline, but no baseline data was passed in.');
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

if ~exist('savePlot','var')
    savePlot = 0;
end

%% Average gait cycle value
% re-sort gait events with the starting event
gait_events_sorted = sortGaitEvents(aligned_data.gait_events,cycleStartEvent);

% Go through all valid gait cycles
% Left
if isfield(signalAnalysisData,'Left')
    if isfield(signalAnalysisData.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    gait_cycle_avg.Left = cell(1,length(signalAnalysisData.Left.Chan_Names));
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{i}),nPercentBins,1);
        count = 1;
        for j = 1:height(gait_events_sorted)-1
            if ~isnan(gait_events_sorted.(cycleStartEvent)(j)) && ~isnan(gait_events_sorted.(cycleStartEvent)(j+1)) && (diff(gait_events_sorted.(cycleStartEvent)(j:j+1)) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Left.Time{i}-gait_events_sorted.(cycleStartEvent)(j)));
                [~,end_ind] = min(abs(signalAnalysisData.Left.Time{i}-gait_events_sorted.(cycleStartEvent)(j+1)));
                
                if isfield(signalAnalysisData.Left,'PSD')
                    data_snip = 20*log10(abs(signalAnalysisData.Left.PSD{i}(:,start_ind:end_ind)));
                else
                    data_snip = abs(signalAnalysisData.Left.Values{i}(:,start_ind:end_ind));
                end
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for k = 1:length(percent_inds)-1
                        if k == 1
                            gait_cycle_mat_left(:,k,count) = mean(data_snip(:,percent_inds(k):percent_inds(k+1)),2);
                        else
                            gait_cycle_mat_left(:,k,count) = mean(data_snip(:,percent_inds(k)+1:percent_inds(k+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycle_avg.Left{i} = mean(gait_cycle_mat_left,3);
    end
end

% Right
if isfield(signalAnalysisData,'Right')
    if isfield(signalAnalysisData.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    gait_cycle_avg.Right = cell(1,length(signalAnalysisData.Right.Chan_Names));
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{i}),nPercentBins,1);
        count = 1;
        for j = 1:height(gait_events_sorted)-1
            if ~isnan(gait_events_sorted.(cycleStartEvent)(j)) && ~isnan(gait_events_sorted.(cycleStartEvent)(j+1)) && (diff(gait_events_sorted.(cycleStartEvent)(j:j+1)) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Right.Time{i}-gait_events_sorted.(cycleStartEvent)(j)));
                [~,end_ind] = min(abs(signalAnalysisData.Right.Time{i}-gait_events_sorted.(cycleStartEvent)(j+1)));
                if isfield(signalAnalysisData.Right,'PSD')
                    data_snip = 20*log10(abs(signalAnalysisData.Right.PSD{i}(:,start_ind:end_ind)));
                else
                    data_snip = abs(signalAnalysisData.Right.Values{i}(:,start_ind:end_ind));
                end
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for k = 1:length(percent_inds)-1
                        if k == 1
                            gait_cycle_mat_right(:,k,count) = mean(data_snip(:,percent_inds(k):percent_inds(k+1)),2);
                        else
                            gait_cycle_mat_right(:,k,count) = mean(data_snip(:,percent_inds(k)+1:percent_inds(k+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycle_avg.Right{i} = mean(gait_cycle_mat_right,3);
    end
end

%% Normalize if set
if ~strcmp(normalizeBy,'none')
    normalization = [];
    if strcmp(normalizeBy,'average_during_walking')
        if isfield(signalAnalysisData,'Left')
            normalization.Left = cell(1,length(signalAnalysisData.Left.Chan_Names));
            walking_start_ind = find(signalAnalysisData.Left.Time{i} >= min(gait_events_sorted{1,:})-1,1,'first');
            walking_end_ind = find(signalAnalysisData.Left.Time{i} <= max(gait_events_sorted{end,:}),1,'last');
            for i = 1:length(signalAnalysisData.Left.Chan_Names)
                normalization.Left{i} = signalAnalysisData.Left.Values{i}(:,walking_start_ind:walking_end_ind);
            end
        end
        if isfield(signalAnalysisData,'Right')
            normalization.Right = cell(1,length(signalAnalysisData.Right.Chan_Names));
            walking_start_ind = find(signalAnalysisData.Right.Time{1} >= min(gait_events_sorted{1,:})-1,1,'first');
            walking_end_ind = find(signalAnalysisData.Right.Time{1} <= max(gait_events_sorted{end,:}),1,'last');
            for i = 1:length(signalAnalysisData.Right.Chan_Names)
                normalization.Right{i} = signalAnalysisData.Right.Values{i}(:,walking_start_ind:walking_end_ind);
            end
        end
    elseif strcmp(normalizeBy,'baseline')
        if isfield(signalAnalysisData,'Left')
            normalization.Left = cell(1,length(signalAnalysisData.Left.Chan_Names));
            baseline_start_ind = find(signalAnalysisData.Left.Time{i} >= baselineTime(1),1,'first');
            baseline_end_ind = find(signalAnalysisData.Left.Time{i} <= baselineTime(2),1,'last');
            for i = 1:length(signalAnalysisData.Left.Chan_Names)
                normalization.Left{i} = signalAnalysisData.Left.Values{i}(:,baseline_start_ind:baseline_end_ind);
            end
        end
        if isfield(signalAnalysisData,'Right')
            normalization.Right = cell(1,length(signalAnalysisData.Right.Chan_Names));
            baseline_start_ind = find(signalAnalysisData.Right.Time{i} >= baselineTime(1),1,'first');
            baseline_end_ind = find(signalAnalysisData.Right.Time{i} <= baselineTime(2),1,'last');
            for i = 1:length(signalAnalysisData.Right.Chan_Names)
                normalization.Right{i} = signalAnalysisData.Right.Values{i}(:,baseline_start_ind:baseline_end_ind);
            end
        end
    end
end

if strcmp(normalizationType,'percent_change')
    if isfield(gait_cycle_avg,'Left')
        for i = 1:length(gait_cycle_avg.Left)
            mu = mean(abs(normalization.Left{i}),2);
            gait_cycle_avg.Left{i} = (gait_cycle_avg.Left{i}-mu)./mu;
        end
    end
    
    if isfield(gait_cycle_avg,'Right')
        for i = 1:length(gait_cycle_avg.Right)
            mu = mean(abs(normalization.Right{i}),2);
            gait_cycle_avg.Right{i} = (gait_cycle_avg.Right{i}-mu)./mu;
        end
    end
elseif strcmp(normalizationType,'zscore')
    if isfield(gait_cycle_avg,'Left')
        for i = 1:length(gait_cycle_avg.Left)
            mu = mean(abs(normalization.Left{i}),2);
            sigma = std(abs(normalization.Left{i}),0,2);
            gait_cycle_avg.Left{i} = (gait_cycle_avg.Left{i}-mu)./sigma;
        end
    end
    
    if isfield(gait_cycle_avg,'Right')
        for i = 1:length(gait_cycle_avg.Right)
            mu = mean(abs(normalization.Right{i}),2);
            sigma = std(abs(normalization.Right{i}),0,2);
            gait_cycle_avg.Right{i} = (gait_cycle_avg.Right{i}-mu)./sigma;
        end
    end
end

%% Plot

fig_vec = [];
if isfield(gait_cycle_avg,'Left')
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signalAnalysisData.Left,'PSD')
            if sum(isnan(gait_cycle_avg.Left{i}),'all') == 0
                ax = pcolor(1:100,signalAnalysisData.Left.Freq_Values{i},gait_cycle_avg.Left{i});
                ylim([2.5,50]);
                ax.Parent.XTick = [1,10:10:100];
                ax.Parent.XTickLabel = {'0','10','20','30','40','50','60','70','80','90','100'};
            end
        else
            if sum(isnan(gait_cycle_avg.Left{i}),'all') == 0
                ax = pcolor(1:100,log2(signalAnalysisData.Left.Freq_Values{i}),gait_cycle_avg.Left{i});
                ticks = logspace(log10(2.5),log10(50),10);
                ax.Parent.YTick = log2(ticks);
                ax.Parent.YTickLabel = ticks;
                ylim([log2(2.5),log2(50)]);
                ax.Parent.XTick = [1,10:10:100];
                ax.Parent.XTickLabel = {'0','10','20','30','40','50','60','70','80','90','100'};
                if strcmp(normalizeBy,'baseline') && strcmp(normalizationType,'zscore')
                    caxis([-2,2]);
                end
            end
        end
        shading interp;
        title({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i}});
        xlabel('Gait Cycle %');
        ylabel('Frequency (Hz)');
    end
end

if isfield(gait_cycle_avg,'Right')
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signalAnalysisData.Right,'PSD')
            if sum(isnan(gait_cycle_avg.Right{i}),'all') == 0
                ax = pcolor(1:100,signalAnalysisData.Right.Freq_Values{i},gait_cycle_avg.Right{i});
                ylim([2.5,50]);
                ax.Parent.XTick = [1,10:10:100];
                ax.Parent.XTickLabel = {'0','10','20','30','40','50','60','70','80','90','100'};
            end
        else
            if sum(isnan(gait_cycle_avg.Right{i}),'all') == 0
                ax = pcolor(1:100,log2(signalAnalysisData.Right.Freq_Values{i}),gait_cycle_avg.Right{i});
                ticks = logspace(log10(2.5),log10(50),10);
                ax.Parent.YTick = log2(ticks);
                ax.Parent.YTickLabel = ticks;
                ylim([log2(2.5),log2(50)]);
                ax.Parent.XTick = [1,10:10:100];
                ax.Parent.XTickLabel = {'0','10','20','30','40','50','60','70','80','90','100'};
                if strcmp(normalizeBy,'baseline') && strcmp(normalizationType,'zscore')
                    caxis([-2,2]);
                end
            end
        end
        shading interp;
        title({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i}});
        xlabel('Gait Cycle %');
        ylabel('Frequency (Hz)');
    end
end

%% Save plots
if savePlot
    save_dir = uigetdir();
    
    figure_format(6,6,10);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec'))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM']))
    end
    
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED']))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED']))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT'))
            mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT'))
            mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{n}));
            end
        end
    end
    
    for i = 1:length(fig_vec)
        curr_axes = gca(fig_vec(i));
        save_name = [];
        for j = 1:length(curr_axes.Title.String)
            if isempty(save_name)
                save_name = curr_axes.Title.String{j};
            else
                save_name = [save_name,' ', curr_axes.Title.String{j}];
            end
        end
        
        if isfield(aligned_data,'trial_num') && ~isempty(aligned_data.trial_num)
            save_name = [save_name,' ',sprintf('Trial%i',aligned_data.trial_num)];
        end
        
        if ~strcmp(normalizationType,'none')
            save_name = [save_name,' ',normalizationType];
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{1},strrep(save_name,' ','_')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end
end