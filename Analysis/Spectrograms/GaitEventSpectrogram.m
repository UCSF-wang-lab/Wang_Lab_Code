function GaitEventSpectrogram(aligned_data,signalAnalysisData,varargin)
%% GaitEventSpectrogram
% Calculates the spectrogram of the entire aligned data passed in. Then,
% the spectrogram will be marked for gait events, denoted by a vertical
% line at the time of the gait event. The aligned data must have a
% "gait_events" field.
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
%               eventsToMark        [=] Cell array of gait events to mark
%                                       on the spectrogram. Default is to
%                                       plot all gait events.
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
%           GaitEventSpectrogram(aligned_data,A);
%
% Date:     05/25/2022
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Option variables
for i = 1:2:nargin-2
    switch varargin{i}
        case 'eventsToMark'
            eventsToMark = varargin{i+1};
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('eventsToMark','var') || isempty(eventsToMark)
    eventsToMark = {'LTO','LHS','RTO','RHS'};
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;
end

%% Plotting
colors = CBMap('GaitEvents',4);
fig_vec = [];

% Left
if isfield(signalAnalysisData,'Left')
    if isfield(signalAnalysisData.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    for i = 1:length(signalAnalysisData.Left.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signalAnalysisData.Left,'PSD')
            pcolor(signalAnalysisData.Left.Time{i},signalAnalysisData.Left.Freq_Values{i},20*log10(abs(signalAnalysisData.Left.PSD{i})));
            shading flat;
            ylim([2.5,50]);
            A = caxis;
            caxis(A.*0.80);
        else
            ax = pcolor(signalAnalysisData.Left.Time{i},log2(signalAnalysisData.Left.Freq_Values{i}),abs(signalAnalysisData.Left.Values{i}));
            ax.EdgeAlpha = 0;
            ticks = logspace(log10(2.5),log10(50),10);
            ax.Parent.YTick = log2(ticks);
            ax.Parent.YTickLabel = ticks;
            ylim([log2(2.5),log2(50)]);
        end
        
        hold on;
        axes_vec = [];
        for j = 1:length(eventsToMark)
            for k = 1:height(aligned_data.gait_events)
                temp = aligned_data.gait_events.(eventsToMark{j})(k);
                if ~isnan(temp)
                    if length(axes_vec) < j
                        axes_vec(j) = xline(temp,'Color',colors.(eventsToMark{j}),'DisplayName',eventsToMark{j});
                    else
                        xline(temp,'Color',colors.(eventsToMark{j}));
                    end
                end
            end
        end
        hold off;
        legend(axes_vec);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({[subjectID,' Left'];signalAnalysisData.Left.Chan_Names{i}});
    end
end

% Right
if isfield(signalAnalysisData,'Right')
    if isfield(signalAnalysisData.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    for i = 1:length(signalAnalysisData.Right.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signalAnalysisData.Right,'PSD')
            pcolor(signalAnalysisData.Right.Time{i},signalAnalysisData.Right.Freq_Values{i},20*log10(abs(signalAnalysisData.Right.PSD{i})));
            shading flat;
            ylim([2.5,50]);
            A = caxis;
            caxis(A.*0.80);
        else
            ax = pcolor(signalAnalysisData.Right.Time{i},log2(signalAnalysisData.Right.Freq_Values{i}),abs(signalAnalysisData.Right.Values{i}));
            ax.EdgeAlpha = 0;
            ticks = logspace(log10(2.5),log10(50),10);
            ax.Parent.YTick = log2(ticks);
            ax.Parent.YTickLabel = ticks;
            ylim([log2(2.5),log2(50)]);
        end
        
        hold on;
        axes_vec = [];
        for j = 1:length(eventsToMark)
            for k = 1:height(aligned_data.gait_events)
                temp = aligned_data.gait_events.(eventsToMark{j})(k);
                if ~isnan(temp)
                    if length(axes_vec) < j
                        axes_vec(j) = xline(temp,'Color',colors.(eventsToMark{j}),'DisplayName',eventsToMark{j});
                    else
                        xline(temp,'Color',colors.(eventsToMark{j}));
                    end
                end
            end
        end
        hold off;
        legend(axes_vec);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({[subjectID,' Right'];signalAnalysisData.Right.Chan_Names{i}});
    end
end


%% Save plots
if savePlot
    save_dir = uigetdir();
    
%     figure_format(6,6,10);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'Spectrogram'))
        mkdir(fullfile(save_dir,'Spectrogram'));
    end
    
    if ~isfolder(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM']))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'FT'))
            mkdir(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'CWT'))
            mkdir(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}));
            end
        end
    end
    
    for i = 1:length(fig_vec)
        curr_axes = gca(fig_vec(i));
        legend off;
        save_name = [];
        for j = 1:length(curr_axes.Title.String)
            if isempty(save_name)
                save_name = curr_axes.Title.String{j};
            else
                save_name = [save_name,' ', curr_axes.Title.String{j}];
            end
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{1},strrep(save_name,' ','_')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'Spectrogram',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        end
    end
end
end