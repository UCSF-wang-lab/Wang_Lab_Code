function varargout = GaitEventAvgCanonicalBandPowerPhase(aligned_data,signalAnalysisData,varargin)
%% GaitEventAvgCanonicalBandPowerPhase
% Calculates the average power and/or phase in the canonical frequency
% bands for different gait events. The canonical frequency bands are delta
% (0.5-3 Hz), theta (4-7 Hz), alpha (8-12 Hz), beta (13-30 Hz), low gamma
% (30-50 Hz). Beta is also broken into high (21-30 Hz) and low (13-20 Hz)
% sub bands. The average phase can also have error bars associated with the
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
%               plotType            [=] How to represent the average power
%                                       and/or phase. Can either be line 
%                                       graph ("power" or "phase"),
%                                       boxplot ("power_boxplot" or
%                                       "phase_boxplot"), boxplot with
%                                       individal points shown
%                                       ("power_boxplot_scatter" or
%                                       "phase_boxplot_scatter"). Phase
%                                       have an extra option of showing the
%                                       plots as a polar plot
%                                       ("phase_polar"). One
%                                       other option is to plot all types
%                                       ("all"). Default is all.
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
%           GaitEventAvgCanonicalBandPowerPhase(aligned_data,A,'prePostTime',1,'eventsToCompare',{{'LHS','RTO','RHS','LTO'}},'spreadType','SE','subjectID','gRCS02','savePlot',0)
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
        case 'plotType'
            plotType = varargin{i+1};
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
    eventsToCompare{1} = {'LHS','RTO','RHS','LTO'};
end

if ~exist('spreadType','var') || isempty(spreadType)
    spreadType = 'STD';
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('plotType','var') || isempty(plotType)
    plotType = 'all';  % all, power, power_boxplot, power_boxplot_scatter, phase, phase_boxplot, phase_boxplot_scatter, phase_polar
end

if ~exist('savePlot','var') || isempty(savePlot)
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
    power_at_event.Left = {};
    average_power.Left = {};
    average_power_time.Left = {};
    phase_at_event.Left = {};
    average_phase.Left = {};
    average_phase_time.Left = {};
    range_power.Left = {};
    range_power_time.Left = {};
    range_phase.Left = {};
    range_phase_time.Left = {};
    [~,band_names] = getFreqBandInd(signalAnalysisData.Left.Freq_Values{1});
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        band_inds = getFreqBandInd(signalAnalysisData.Left.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Left,aligned_data.gait_events.Properties.VariableNames{j})
                power_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                average_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                phase_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                average_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
            end
            vals_power = [];
            vals_power_time = [];
            vals_phase = [];
            vals_phase_time = [];
            count = 1;
            count_time = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signalAnalysisData.Left.Time{i}-(event_time-prePostTime)));
                    [~,min_ind_post] = min(abs(signalAnalysisData.Left.Time{i}-(event_time+prePostTime)));
                    [~,event_ind] = min(abs(signalAnalysisData.Left.Time{i}-event_time));
                    
                    if isfield(signalAnalysisData.Left,'PSD')
                        temp = 20*log10(abs(signalAnalysisData.Left.Values{i}(:,event_ind)));
                        temp1 = 20*log10(abs(signalAnalysisData.Left.Values{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signalAnalysisData.Left.Values{i}(:,event_ind));
                        temp3 = angle(signalAnalysisData.Left.Values{i}(:,min_ind_pre:min_ind_post));
                    else
                        temp = abs(signalAnalysisData.Left.Values{i}(:,event_ind));
                        temp1 = abs(signalAnalysisData.Left.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signalAnalysisData.Left.Values{i}(:,event_ind));
                        temp3 = angle(signalAnalysisData.Left.Values{i}(:,min_ind_pre:min_ind_post));
                    end
                    
                    if size(temp1,2) == round((prePostTime*2)/time_res_left)+1
                        if sum(isinf(temp1),'all') == 0
                            for m = 1:length(band_names)
                                vals_power_time(count_time,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count_time,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                            end
                            count_time = count_time + 1;
                        else
                            fprintf("dropped packet");
                        end
                        
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            power_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_power);
            average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power,1));
            average_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power_time,1));
            phase_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_phase);
            average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            average_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase_time,1));
            
            if strcmp(spreadType,'SE')
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1))./sqrt(size(vals_power,1));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1))./sqrt(size(vals_power_time,1));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./sqrt(size(vals_phase,1));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1))./sqrt(size(vals_phase_time,1));
            else
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1));
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
    power_at_event.Right = {};
    average_power.Right = {};
    average_power_time.Right = {};
    phase_at_event.Right = {};
    average_phase.Right = {};
    average_phase_time.Right = {};
    range_power.Right = {};
    range_power_time.Right = {};
    range_phase.Right = {};
    range_phase_time.Right = {};
    [~,band_names] = getFreqBandInd(signalAnalysisData.Right.Freq_Values{1});
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        band_inds = getFreqBandInd(signalAnalysisData.Right.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Right,aligned_data.gait_events.Properties.VariableNames{j})
                power_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                average_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                phase_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                average_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
            end
            vals_power = [];
            vals_power_time = [];
            vals_phase = [];
            vals_phase_time = [];
            count = 1;
            count_time = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signalAnalysisData.Right.Time{i}-(event_time-prePostTime)));
                    [~,min_ind_post] = min(abs(signalAnalysisData.Right.Time{i}-(event_time+prePostTime)));
                    [~,event_ind] = min(abs(signalAnalysisData.Right.Time{i}-event_time));
                    
                    if isfield(signalAnalysisData.Right,'PSD')
                        temp = 20*log10(abs(signalAnalysisData.Right.Values{i}(:,event_ind)));
                        temp1 = 20*log10(abs(signalAnalysisData.Right.Values{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signalAnalysisData.Right.Values{i}(:,event_ind));
                        temp3 = angle(signalAnalysisData.Right.Values{i}(:,min_ind_pre:min_ind_post));
                    else
                        temp = 20*log10(abs(signalAnalysisData.Right.Values{i}(:,event_ind)));
                        temp1 = abs(signalAnalysisData.Right.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signalAnalysisData.Right.Values{i}(:,event_ind));
                        temp3 = angle(signalAnalysisData.Right.Values{i}(:,min_ind_pre:min_ind_post));
                    end
                    
                    if size(temp1,2) == round((prePostTime*2)/time_res_right)+1
                        if sum(isinf(temp1),'all') == 0
                            for m = 1:length(band_names)
                                vals_power_time(count_time,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count_time,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                            end
                            count_time = count_time + 1;
                        else
                            fprintf("dropped packet");
                        end
                        
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            power_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_power);
            average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power,1));
            average_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power_time,1));
            phase_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_phase);
            average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            average_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase_time,1));
            
            if strcmp(spreadType,'SE')
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1))./sqrt(size(vals_power,1));
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1))./sqrt(size(vals_power_time,1));
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./sqrt(size(vals_phase,1));
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1))./sqrt(size(vals_phase_time,1));
            else
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1));
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1));
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1));
            end
        end
    end
end

%% Plot
fig_vec = [];
colors = CBMap('GaitEvents',4);

if strcmp(plotType,'all') || strcmp(plotType,'power')
    if isfield(average_power_time,'Left')
        time_vec = -1:time_res_left:1;
        for i = 1:length(signalAnalysisData.Left.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(eventsToCompare{j})
                        mean_std_plot(time_vec,average_power_time.Left.(eventsToCompare{j}{m}){i}(:,k),range_power_time.Left.(eventsToCompare{j}{m}){i}(:,k),ax_hand,colors.(eventsToCompare{j}{m}),[]);
                        %                     mean_std_plot(time_vec,average_phase.Left.(event_compare{j}{m}){i}(:,k),std_phase.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signalAnalysisData.Left,'PSD')
                        ylabel('mV^2'); % mV^2
                    else
                        ylabel('Magnitude');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
    
    if isfield(average_power_time,'Right')
        time_vec = -1:time_res_right:1;
        for i = 1:length(signalAnalysisData.Right.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(eventsToCompare{j})
                        mean_std_plot(time_vec,average_power_time.Right.(eventsToCompare{j}{m}){i}(:,k),range_power_time.Right.(eventsToCompare{j}{m}){i}(:,k),ax_hand,colors.(eventsToCompare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signalAnalysisData.Right,'PSD')
                        ylabel('mV^2'); % mV^2
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
end

if strcmp(plotType,'all') || strcmp(plotType,'phase')
    if isfield(average_phase_time,'Left')
        time_vec = -1:time_res_left:1;
        for i = 1:length(signalAnalysisData.Left.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(eventsToCompare{j})
                        mean_std_plot(time_vec,average_phase_time.Left.(eventsToCompare{j}{m}){i}(:,k),range_phase_time.Left.(eventsToCompare{j}{m}){i}(:,k),ax_hand,colors.(eventsToCompare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signalAnalysisData.Left,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
    
    if isfield(average_phase_time,'Right')
        time_vec = -1:time_res_right:1;
        for i = 1:length(signalAnalysisData.Right.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(eventsToCompare{j})
                        mean_std_plot(time_vec,average_phase_time.Right.(eventsToCompare{j}{m}){i}(:,k),range_phase_time.Right.(eventsToCompare{j}{m}){i}(:,k),ax_hand,colors.(eventsToCompare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signalAnalysisData.Right,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
end

if strcmp(plotType,'all') || contains(plotType,'power_boxplot')
    if isfield(power_at_event,'Left')
        for i = 1:length(signalAnalysisData.Left.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(eventsToCompare{j})
                        box_vals = [box_vals;power_at_event.Left.(eventsToCompare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({eventsToCompare{j}{m}},size(power_at_event.Left.(eventsToCompare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plotType,'scatter')
                        for n = 1:length(eventsToCompare{j})
                            inds = strcmp(box_group,eventsToCompare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(eventsToCompare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    
                    if isfield(signalAnalysisData.Left,'PSD')
                        ylabel('mV^2');
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
    
    if isfield(power_at_event,'Right')
        for i = 1:length(signalAnalysisData.Right.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(eventsToCompare{j})
                        box_vals = [box_vals;power_at_event.Right.(eventsToCompare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({eventsToCompare{j}{m}},size(power_at_event.Right.(eventsToCompare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plotType,'scatter')
                        for n = 1:length(eventsToCompare{j})
                            inds = strcmp(box_group,eventsToCompare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(eventsToCompare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    if isfield(signalAnalysisData.Right,'PSD')
                        ylabel('mV^2');
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
end

if strcmp(plotType,'all') || contains(plotType,'phase_boxplot')
    if isfield(phase_at_event,'Left')
        for i = 1:length(signalAnalysisData.Left.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(eventsToCompare{j})
                        box_vals = [box_vals;phase_at_event.Left.(eventsToCompare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({eventsToCompare{j}{m}},size(phase_at_event.Left.(eventsToCompare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plotType,'scatter')
                        for n = 1:length(eventsToCompare{j})
                            inds = strcmp(box_group,eventsToCompare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(eventsToCompare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    ylabel('Phase');
                    ylim([-pi,pi]);
                    yticks([-pi:pi/4:pi]);
                    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
                    set(gca,'TickLabelInterpreter','tex')
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
    
    if isfield(phase_at_event,'Right')
        for i = 1:length(signalAnalysisData.Right.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(eventsToCompare{j})
                        box_vals = [box_vals;phase_at_event.Right.(eventsToCompare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({eventsToCompare{j}{m}},size(phase_at_event.Right.(eventsToCompare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plotType,'scatter')
                        for n = 1:length(eventsToCompare{j})
                            inds = strcmp(box_group,eventsToCompare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(eventsToCompare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    ylabel('Phase');
                    ylim([-pi,pi]);
                    yticks([-pi:pi/4:pi]);
                    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
                    set(gca,'TickLabelInterpreter','tex')
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
end

if strcmp(plotType,'all') || strcmp(plotType,'phase_polar')
    if isfield(phase_at_event,'Left')
        for i = 1:length(signalAnalysisData.Left.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    polar_handle = polaraxes('Units',ax_hand.Units,'Position',ax_hand.Position);
                    delete(ax_hand);
                    hold(polar_handle,'on');
                    for m = 1:length(eventsToCompare{j})
                        polarhistogram(polar_handle,wrapTo2Pi(phase_at_event.Left.(eventsToCompare{j}{m}){i}(:,k)),10,'FaceColor',colors.(eventsToCompare{j}{m}),'FaceAlpha',0.3);
                    end
                    hold(polar_handle,'off');
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
    
    if isfield(phase_at_event,'Right')
        for i = 1:length(signalAnalysisData.Right.Chan_Names)
            for j = 1:length(eventsToCompare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    polar_handle = polaraxes('Units',ax_hand.Units,'Position',ax_hand.Position);
                    delete(ax_hand);
                    hold(polar_handle,'on');
                    for m = 1:length(eventsToCompare{j})
                        polarhistogram(polar_handle,phase_at_event.Right.(eventsToCompare{j}{m}){i}(:,k),10,'FaceColor',colors.(eventsToCompare{j}{m}),'FaceAlpha',0.3);
                    end
                    hold(polar_handle,'off');
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i};createPlotTitle(eventsToCompare{j})});
            end
        end
    end
end

%% Save plots (NEED TO CHANGE; commented out for now)
if savePlot
    save_dir = uigetdir();
    
    figure_format(12,8,12);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgEventPower'))
        mkdir(fullfile(save_dir,'AvgEventPower'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM']))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'FT'))
            mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'CWT'))
            mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}));
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
        
        if contains(plotType,'boxplot') || contains(plotType,'phase_polar')
            save_name = [save_name,' ', plotType];
        else
            save_name = [save_name,' ', plotType, ' ', spreadType];
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end

varargout = {power_at_event, average_power, average_power_time,...
    phase_at_event, average_phase, average_phase_time,...
    range_power, range_power_time,...
    range_phase, range_phase_time};
end

function plot_title = createPlotTitle(s)
plot_title = s{1};

for i = 2:length(s)
    plot_title = plot_title + " vs. " + s{i};
end
end