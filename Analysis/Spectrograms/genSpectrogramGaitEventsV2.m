function genSpectrogramGaitEventsV2(aligned_data,signal_analysis_data,events_to_mark,subjectID,save_flag)
if ~exist('events_to_mark','var')
    events_to_mark = {'LTO','LHS','RTO','RHS'};
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var')
    save_flag = 0;
end

%% Plotting
colors = CBMap('GaitEvents',4);
fig_vec = [];

% Left
if isfield(signal_analysis_data,'Left')
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signal_analysis_data.Left,'PSD')
            pcolor(signal_analysis_data.Left.Time{i},signal_analysis_data.Left.Freq_Values{i},10*log10(abs(signal_analysis_data.Left.PSD{i})));
            shading flat;
            ylim([2.5,50]);
            A = caxis;
            caxis(A.*0.80);
        else
            ax = pcolor(signal_analysis_data.Left.Time{i},log2(signal_analysis_data.Left.Freq_Values{i}),abs(signal_analysis_data.Left.Values{i}));
            ax.EdgeAlpha = 0;
            ticks = logspace(log10(2.5),log10(50),10);
            ax.Parent.YTick = log2(ticks);
            ax.Parent.YTickLabel = ticks;
            ylim([log2(2.5),log2(50)]);
        end
        
        hold on;
        axes_vec = [];
        for j = 1:length(events_to_mark)
            for k = 1:height(aligned_data.gait_events)
                temp = aligned_data.gait_events.(events_to_mark{j})(k);
                if ~isnan(temp)
                    if length(axes_vec) < j
                        axes_vec(j) = xline(temp,'Color',colors.(events_to_mark{j}),'DisplayName',events_to_mark{j});
                    else
                        xline(temp,'Color',colors.(events_to_mark{j}));
                    end
                end
            end
        end
        hold off;
        legend(axes_vec);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i}});
    end
end

% Right
if isfield(signal_analysis_data,'Right')
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signal_analysis_data.Right,'PSD')
            pcolor(signal_analysis_data.Right.Time{i},signal_analysis_data.Right.Freq_Values{i},10*log10(abs(signal_analysis_data.Right.PSD{i})));
            shading flat;
            ylim([2.5,50]);
            A = caxis;
            caxis(A.*0.80);
        else
            ax = pcolor(signal_analysis_data.Right.Time{i},log2(signal_analysis_data.Right.Freq_Values{i}),abs(signal_analysis_data.Right.Values{i}));
            ax.EdgeAlpha = 0;
            ticks = logspace(log10(2.5),log10(50),10);
            ax.Parent.YTick = log2(ticks);
            ax.Parent.YTickLabel = ticks;
            ylim([log2(2.5),log2(50)]);
        end
        
        hold on;
        axes_vec = [];
        for j = 1:length(events_to_mark)
            for k = 1:height(aligned_data.gait_events)
                temp = aligned_data.gait_events.(events_to_mark{j})(k);
                if ~isnan(temp)
                    if length(axes_vec) < j
                        axes_vec(j) = xline(temp,'Color',colors.(events_to_mark{j}),'DisplayName',events_to_mark{j});
                    else
                        xline(temp,'Color',colors.(events_to_mark{j}));
                    end
                end
            end
        end
        hold off;
        legend(axes_vec);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i}});
    end
end


%% Save plots
if save_flag
    save_dir = uigetdir();
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'Spectrogram'))
        mkdir(fullfile(save_dir,'Spectrogram'));
    end
    
    if isfield(signal_analysis_data,'PSD')
        if ~isfolder(fullfile(save_dir,'Spectrogram','FT'))
            mkdir(fullfile(save_dir,'Spectrogram','FT'))
        end
    else
        if ~isfolder(fullfile(save_dir,'Spectrogram','CWT'))
            mkdir(fullfile(save_dir,'Spectrogram','CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    for n = 1:length(folders_to_check)
        if isfield(signal_analysis_data,'PSD')
            if ~isfolder(fullfile(save_dir,'Spectrogram','FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'Spectrogram','FT',folders_to_check{n}));
            end
        else
            if ~isfolder(fullfile(save_dir,'Spectrogram','CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'Spectrogram','CWT',folders_to_check{n}));
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
        
        if isfield(signal_analysis_data,'PSD')
            savefig(fig_vec(i),fullfile(save_dir,'Spectrogram','FT',folders_to_check{1},save_name));
        else
            savefig(fig_vec(i),fullfile(save_dir,'Spectrogram','CWT',folders_to_check{1},save_name));
        end
    end
end
end