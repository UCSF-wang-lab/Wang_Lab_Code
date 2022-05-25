function GaitEventPSD(aligned_data,signalAnalysisData,varargin)
%% GaitEventPSD
% Calculates the average PSD for each gait event in each recorded area. The
% resulting plots show the average PSD value in a thicker line, and also
% shows the PSD values for each gait event in a lighter shade.
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
%               subjectID           [=] String variable of the name of the
%                                       subject being analyzed.
%
%               savePlot            [=] Boolean option to save the 
%                                       resulting plot. Default is false.
%
%   Example call:
%           load(<filename>)
%           A = calcRCS_STFT(aligned_data,[],1,0.9,[]);
%           GaitEventsPSD(aligned_data,A);
%
% Date:     05/25/2022
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Option variables
for i = 1:2:nargin-2
    switch varargin{i}
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;
end

%% PSD at each gait event
% Left
if isfield(signalAnalysisData,'Left')
    PSD_gait_events.Left = {};
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(PSD_gait_events.Left,aligned_data.gait_events.Properties.VariableNames{j})
                PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Left.Chan_Names));
            end
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind] = min(abs(signalAnalysisData.Left.Time{i}-event_time));
                    power_values = signalAnalysisData.Left.PSD{i}(:,min_ind);
                    if sum(power_values == 0) == 0
                        if isempty(PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}))
                            PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = power_values;
                        else
                            PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i}(:,end+1) = power_values;
                        end
                    end
                end
            end
        end
    end
end

% Right
if isfield(signalAnalysisData,'Right')
    PSD_gait_events.Right = {};
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(PSD_gait_events.Right,aligned_data.gait_events.Properties.VariableNames{j})
                PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signalAnalysisData.Right.Chan_Names));
            end
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind] = min(abs(signalAnalysisData.Right.Time{i}-event_time));
                    power_values = signalAnalysisData.Right.PSD{i}(:,min_ind);
                    if sum(power_values == 0) == 0
                        if isempty(PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}))
                            PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = power_values;
                        else
                            PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i}(:,end+1) = power_values;
                        end
                    end
                end
            end
        end
    end
end

%% Plot
fig_vec = [];
if isfield(PSD_gait_events,'Left')
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        fig_vec(end+1) = figure;
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            curr_event = aligned_data.gait_events.Properties.VariableNames{j};
            freq_mat = repmat(signalAnalysisData.Left.Freq_Values{i},1,size(PSD_gait_events.Left.(curr_event){i},2));
            avg_power = mean(10*log10(abs(PSD_gait_events.Left.(curr_event){i})),2);
            
            switch curr_event
                case 'LTO'
                    subplot(2,2,1);
                case 'LHS'
                    subplot(2,2,2);
                case 'RTO'
                    subplot(2,2,3);
                case 'RHS'
                    subplot(2,2,4);
            end
            
            plot(freq_mat,10*log10(abs(PSD_gait_events.Left.(curr_event){i})),'Color',[0,0,0,0.2]);
            hold on;
            plot(freq_mat(:,1),avg_power,'-k','linewidth',1.5);
            hold off;
            title(aligned_data.gait_events.Properties.VariableNames{j});
            xlim([2.5,40])
            xlabel('Frequency (Hz)');
            ylabel('db/Hz');
        end
        sgtitle({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i}});
    end
end

% Right
if isfield(PSD_gait_events,'Right')
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        fig_vec(end+1) = figure;
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            curr_event = aligned_data.gait_events.Properties.VariableNames{j};
            freq_mat = repmat(signalAnalysisData.Right.Freq_Values{i},1,size(PSD_gait_events.Right.(curr_event){i},2));
            avg_power = mean(10*log10(abs(PSD_gait_events.Right.(curr_event){i})),2);
            
            switch curr_event
                case 'LTO'
                    subplot(2,2,1);
                case 'LHS'
                    subplot(2,2,2);
                case 'RTO'
                    subplot(2,2,3);
                case 'RHS'
                    subplot(2,2,4);
            end
            
            plot(freq_mat,10*log10(abs(PSD_gait_events.Right.(curr_event){i})),'Color',[0,0,0,0.2]);
            hold on;
            plot(freq_mat(:,1),avg_power,'-k','linewidth',1.5);
            hold off;
            title(aligned_data.gait_events.Properties.VariableNames{j});
            xlim([2.5,40])
            xlabel('Frequency (Hz)');
            ylabel('db/Hz');
        end
        sgtitle({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i}});
    end
end

if savePlot
    save_dir = uigetdir();
    
    figure_format(6,6,10);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'PSD'))
        mkdir(fullfile(save_dir,'PSD'));
    end
    
    if ~isfolder(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM']))
    end
    
    if ~isfolder(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT'))
        mkdir(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT'))
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if ~isfolder(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}))
            mkdir(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}));
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
        
        savefig(fig_vec(i),fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{1},strrep(save_name,' ','_')));
        for k = 2:length(folders_to_check)
            print(fig_vec(i),[fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
        end
        
    end
end
end