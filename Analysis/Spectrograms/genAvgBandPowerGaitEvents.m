function genAvgBandPowerGaitEvents(aligned_data,signal_analysis_data,pre_post_time,event_compare,subjectID,save_flag)

if ~exist('pre_post_time','var') || isempty(pre_post_time)
    pre_post_time = 1;
end

if ~exist('event_compare','var') || isempty(event_compare)
    event_compare{1} = {'LTO','RTO'};
    event_compare{2} = {'LHS','RHS'};
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var') || isempty(save_flag)
    save_flag = 0;
end

%% Extract data
if isfield(signal_analysis_data,'Left')
    average_power.Left = {};
    std_power.Left = {};
    [~,band_names] = getFreqBandInd(signal_analysis_data.Left.Freq_Values{1});
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        band_inds = getFreqBandInd(signal_analysis_data.Left.Freq_Values{i});
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Left,aligned_data.gait_events.Properties.VariableNames{j})
                average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                std_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
            end
            vals = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signal_analysis_data.Left.Time{i}-(event_time-pre_post_time)));
                    [~,min_ind_post] = min(abs(signal_analysis_data.Left.Time{i}-(event_time+pre_post_time)));
                    
                    temp = signal_analysis_data.Left.PSD{i}(:,min_ind_pre:min_ind_post);
                    if sum(temp==0,'all') == 0
                        for m = 1:length(band_names) 
                            vals(count,:,m) = mean(10*log10(abs(temp(band_inds(m,1):band_inds(m,2),:))));
                        end
                        count = count + 1;
                    end
                end
            end
            average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals,1));
            std_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals,0,1));
        end
    end
end

if isfield(signal_analysis_data,'Right')
    average_power.Right = {};
    std_power.Right = {};
    [~,band_names] = getFreqBandInd(signal_analysis_data.Right.Freq_Values{1});
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        band_inds = getFreqBandInd(signal_analysis_data.Right.Freq_Values{i});
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Right,aligned_data.gait_events.Properties.VariableNames{j})
                average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                std_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
            end
            vals = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signal_analysis_data.Right.Time{i}-(event_time-pre_post_time)));
                    [~,min_ind_post] = min(abs(signal_analysis_data.Right.Time{i}-(event_time+pre_post_time)));
                    
                    temp = signal_analysis_data.Right.PSD{i}(:,min_ind_pre:min_ind_post);
                    if sum(temp==0,'all') == 0
                        for m = 1:length(band_names) 
                            vals(count,:,m) = mean(10*log10(abs(temp(band_inds(m,1):band_inds(m,2),:))));
                        end
                        count = count + 1;
                    end
                end
            end
            average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals,1));
            std_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals,0,1));
        end
    end
end

%% Plot
fig_vec = [];
colors = CBMap('GaitEvents',4);

time_vec = -1:1/unique(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate):1;
for i = 1:length(signal_analysis_data.Left.Chan_Names)
    for j = 1:length(event_compare)
        fig_vec(end+1) = figure;
        for k = 1:length(band_names)
            ax_hand = subplot(2,4,k);
            for m = 1:length(event_compare{j})
                mean_std_plot(time_vec,average_power.Left.(event_compare{j}{m}){i}(:,k),std_power.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
            end
            hold(ax_hand,'on');
            xline(0,'--k');
            hold(ax_hand,'off');
            xlabel('Time (s)');
            ylabel('db/Hz');
            title(band_names{k});
        end
        sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
    end
end

time_vec = -1:1/unique(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate):1;
for i = 1:length(signal_analysis_data.Right.Chan_Names)
    for j = 1:length(event_compare)
        fig_vec(end+1) = figure;
        for k = 1:length(band_names)
            ax_hand = subplot(2,4,k);
            for m = 1:length(event_compare{j})
                mean_std_plot(time_vec,average_power.Right.(event_compare{j}{m}){i}(:,k),std_power.Right.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
            end
            hold(ax_hand,'on');
            xline(0,'--k');
            hold(ax_hand,'off');
            xlabel('Time (s)');
            ylabel('db/Hz');
            title(band_names{k});
        end
        sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
    end
end

%% Save plots
if save_flag
    save_dir = uigetdir();
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgPower'))
        mkdir(fullfile(save_dir,'AvgPower'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgPower','FT'))
        mkdir(fullfile(save_dir,'AvgPower','FT'))
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    for n = 1:length(folders_to_check)
        if ~isfolder(fullfile(save_dir,'AvgPower','FT',folders_to_check{n}))
            mkdir(fullfile(save_dir,'AvgPower','FT',folders_to_check{n}));
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
        
        savefig(fig_vec(i),fullfile(save_dir,'AvgPower','FT',folders_to_check{1},save_name));
    end
end
end

function plot_title = createPlotTitle(s)
plot_title = s{1};

for i = 2:length(s)
    plot_title = plot_title + " vs. " + s{i};
end
end