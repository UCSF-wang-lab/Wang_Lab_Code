function GaitEventAvgCanonicalBandPower(aligned_data,signalAnalysisData,varargin)
%% GaitEventAvgCanonicalBandPower
% Calculates the average power in the canonical frequency bands for
% different gait events. The canonical frequency bands are delta (0.5-3
% Hz), theta (4-7 Hz), alpha (8-12 Hz), beta (13-30 Hz), low gamma (30-50
% Hz). Beta is also broken into high (21-30 Hz) and low (13-20 Hz) sub
% bands. The average power can also have error bars associated with the
% plot and can be either the standard deviation or standard error.
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
%               prePostTime         [=] How far back and forward to look
%                                       after a gait event. Default is +-1
%                                       second.
%
%               eventsToCompare     [=] Cell array of gait events to
%                                       compare on the same plot. Default
%                                       is to show all on the same plot. 
%
%               spreadType          [=] How the error of the average should
%                                       be calculated. Can be either the
%                                       standard deviation (STD; default)
%                                       or standard error (SE).
%
%               subjectID           [=] String variable of the name of the
%                                       subject being analyzed.
%
%               savePlot            [=] Boolean option to save the 
%                                       resulting plot. Default is false.
%
%   Example call:
%           load(<filename>)
%           A = calcRCS_STFT(aligned_data,[],1,0.9,[]);
%           GaitEventAvgCanonicalBandPower(aligned_data,A,'prePostTime',1,'eventsToCompare',{{'LHS','RTO','RHS','LTO'}},'spreadType','SE','subjectID','gRCS02','savePlot',0)
%
% Date:     05/25/2022
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Option variables
for i = 1:2:nargin-2
    switch varargin{i}
        case 'prePostTime'
            prePostTime = varargin{i+1};
        case 'eventsToCompare'
            eventsToCompare = varargin{i+1};
        case 'spreadType'
            spreadType = varargin{i+1};
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('prePostTime','var') || isempty(prePostTime)
    prePostTime = 1;
end

if ~exist('eventsToCompare','var') || isempty(eventsToCompare)
    eventsToCompare{1} = {'LTO','RTO'};
    eventsToCompare{2} = {'LHS','RHS'};
end

if ~exist('spreadType','var') || isempty(spreadType)
    spreadType = 'STD';
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var') || isempty(savePlot)
    savePlot = 0;
end

%% Extract data
if isfield(signalAnalysisData,'Left')
    left_sr = uniquetol(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate,1);
    time_res_left = uniquetol(diff(signalAnalysisData.Left.Time{1}),1);
    if isfield(signalAnalysisData.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    average_power.Left = {};
    std_power.Left = {};
    [~,band_names] = getCanonicalFreqBandInd(signalAnalysisData.Left.Freq_Values{1});
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        band_inds = getCanonicalFreqBandInd(signalAnalysisData.Left.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Left,aligned_data.gait_events.Properties.VariableNames{j})
                average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                std_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                std_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
            end
            vals = [];
            vals_phase = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signalAnalysisData.Left.Time{i}-(event_time-prePostTime)));
                    [~,min_ind_post] = min(abs(signalAnalysisData.Left.Time{i}-(event_time+prePostTime)));
                    [~,event_ind] = min(abs(signalAnalysisData.Left.Time{i}-event_time));
                    if isfield(signalAnalysisData.Left,'PSD')
                        temp = 10*log10(abs(signalAnalysisData.Left.PSD{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signalAnalysisData.Left.PSD{i}(:,event_ind));
                    else
                        temp = abs(signalAnalysisData.Left.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signalAnalysisData.Left.Values{i}(:,event_ind));
                    end
                    
                    if size(temp,2) == round((prePostTime*2)/time_res_left)+1
%                     if size(temp,2) == left_sr*2+1
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals,1));
            average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            
            if strcmp(spreadType,'SE')
                std_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals,0,1))./sqrt(size(vals,1));
                std_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./sqrt(size(vals_phase,1));
            else
                std_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals,0,1));
                std_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
            end
        end
    end
end

if isfield(signalAnalysisData,'Right')
    right_sr = uniquetol(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate,1);
    time_res_right = uniquetol(diff(signalAnalysisData.Right.Time{1}),1);
    if isfield(signalAnalysisData.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    average_power.Right = {};
    std_power.Right = {};
    [~,band_names] = getCanonicalFreqBandInd(signalAnalysisData.Right.Freq_Values{1});
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        band_inds = getCanonicalFreqBandInd(signalAnalysisData.Right.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Right,aligned_data.gait_events.Properties.VariableNames{j})
                average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                std_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                std_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
            end
            vals = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signalAnalysisData.Right.Time{i}-(event_time-prePostTime)));
                    [~,min_ind_post] = min(abs(signalAnalysisData.Right.Time{i}-(event_time+prePostTime)));
                    [~,event_ind] = min(abs(signalAnalysisData.Right.Time{i}-event_time));
                    
                    if isfield(signalAnalysisData.Right,'PSD')
                        temp = 10*log10(abs(signalAnalysisData.Right.PSD{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signalAnalysisData.Right.PSD{i}(:,event_ind));
                    else
                        temp = abs(signalAnalysisData.Right.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signalAnalysisData.Right.Values{i}(:,event_ind));
                    end
                    
                    if size(temp,2) == round((prePostTime*2)/time_res_right)+1
%                     if size(temp,2) == right_sr*2+1
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals,1));
            average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            
            if strcmp(spreadType,'SE')
                std_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals,0,1))./size(vals,1);
                std_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./size(vals_phase,1);
            else
                std_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals,0,1));
                std_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
            end
        end
    end
end

%% Plot
fig_vec = [];
colors = CBMap('GaitEvents',4);

if isfield(average_power,'Left')
%     time_vec = -1:1/unique(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate):1;
    time_vec = -1:time_res_left:1;
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        for j = 1:length(eventsToCompare)
            fig_vec(end+1) = figure;
            for k = 1:length(band_names)
                ax_hand = subplot(2,4,k);
                for m = 1:length(eventsToCompare{j})
                    mean_std_plot(time_vec,average_power.Left.(eventsToCompare{j}{m}){i}(:,k),std_power.Left.(eventsToCompare{j}{m}){i}(:,k),ax_hand,colors.(eventsToCompare{j}{m}),[]);
                end
                hold(ax_hand,'on');
                xline(0,'--k');
                hold(ax_hand,'off');
                xlabel('Time (s)');
                
                if isfield(signalAnalysisData.Left,'PSD')
                    ylabel('db/Hz');
                else
                    ylabel('Magnitude');
                end
                
                title(band_names{k});
            end
            sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
        end
    end
end

if isfield(average_power,'Right')
%     time_vec = -1:1/unique(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate):1;
    time_vec = -1:time_res_right:1;
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        for j = 1:length(eventsToCompare)
            fig_vec(end+1) = figure;
            for k = 1:length(band_names)
                ax_hand = subplot(2,4,k);
                for m = 1:length(eventsToCompare{j})
                    mean_std_plot(time_vec,average_power.Right.(eventsToCompare{j}{m}){i}(:,k),std_power.Right.(eventsToCompare{j}{m}){i}(:,k),ax_hand,colors.(eventsToCompare{j}{m}),[]);
                end
                hold(ax_hand,'on');
                xline(0,'--k');
                hold(ax_hand,'off');
                xlabel('Time (s)');
                
                if isfield(signalAnalysisData.Right,'PSD')
                    ylabel('db/Hz');
                else
                    ylabel('Magnitude');
                end
                title(band_names{k});
            end
            sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
        end
    end
end

%% Save plots
if savePlot
    save_dir = uigetdir();
    
    figure_format(12,8,12);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgEventPower'))
        mkdir(fullfile(save_dir,'AvgEventPower'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED']))
        mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED']))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT'))
            mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT'))
            mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{n}));
            end
        end
    end
    
    for i = 1:length(fig_vec)
        curr_axes = gca(fig_vec(i));
        save_name = [];
        for j = 1:length(curr_axes.Parent.Children(1).String)
            if isempty(save_name)
                save_name = curr_axes.Parent.Children(1).String{j};
            else
                save_name = [save_name,' ', curr_axes.Parent.Children(1).String{j}];
            end
        end
        
        if isfield(aligned_data,'trial_num')
            save_name = [save_name,' ',sprintf('Trial%i',aligned_data.trial_num)];
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end
end

function plot_title = createPlotTitle(s)
plot_title = s{1};

for i = 2:length(s)
    plot_title = plot_title + " vs. " + s{i};
end
end