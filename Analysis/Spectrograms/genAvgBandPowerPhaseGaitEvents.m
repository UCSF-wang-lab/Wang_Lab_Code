function genAvgBandPowerPhaseGaitEvents(aligned_data,signal_analysis_data,pre_post_time,event_compare,spread_type,subjectID,plot_type,save_flag)

if ~exist('pre_post_time','var') || isempty(pre_post_time)
    pre_post_time = 1;
end

if ~exist('event_compare','var') || isempty(event_compare)
    event_compare{1} = {'LTO','RTO'};
    event_compare{2} = {'LHS','RHS'};
end

if ~exist('spread_type','var') || isempty(spread_type)
    spread_type = 'STD';
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('plot_type','var') || isempty(plot_type)
    plot_type = 'all';  % all, power, power_boxplot, phase, phase_boxplot
end

if ~exist('save_flag','var') || isempty(save_flag)
    save_flag = 0;
end

%% Extract data
if isfield(signal_analysis_data,'Left')
    left_sr = uniquetol(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate,1);
    time_res_left = uniquetol(diff(signal_analysis_data.Left.Time{1}),1);
    if isfield(signal_analysis_data.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    average_power.Left = {};
    average_power_time.Left = {};
    average_phase.Left = {};
    average_phase_time.Left = {};
    range_power.Left = {};
    range_power_time.Left = {};
    range_phase.Left = {};
    range_phase_time.Left = {};
    [~,band_names] = getFreqBandInd(signal_analysis_data.Left.Freq_Values{1});
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        band_inds = getFreqBandInd(signal_analysis_data.Left.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Left,aligned_data.gait_events.Properties.VariableNames{j})
                average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
            end
            vals_power = [];
            vals_power_time = [];
            vals_phase = [];
            vals_phase_time = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signal_analysis_data.Left.Time{i}-(event_time-pre_post_time)));
                    [~,min_ind_post] = min(abs(signal_analysis_data.Left.Time{i}-(event_time+pre_post_time)));
                    [~,event_ind] = min(abs(signal_analysis_data.Left.Time{i}-event_time));
                    if isfield(signal_analysis_data.Left,'PSD')
                        temp = 10*log10(abs(signal_analysis_data.Left.Values{i}(:,event_ind)));
                        temp1 = 10*log10(abs(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post));
                    else
                        temp = abs(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp1 = abs(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post));
                    end
                    
                    if size(temp1,2) == round((pre_post_time*2)/time_res_left)+1
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_power_time(count,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power,1));
            average_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power_time,1));
            average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            average_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase_time,1));
            
            if strcmp(spread_type,'SE')
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1))./size(vals_power,1);
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1))./size(vals_power_time,1);
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./size(vals_phase,1);
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1))./size(vals_phase_time,1);
            else
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1));
            end
        end
    end
end

if isfield(signal_analysis_data,'Right')
    right_sr = uniquetol(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate,1);
    time_res_right = uniquetol(diff(signal_analysis_data.Right.Time{1}),1);
    if isfield(signal_analysis_data.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    average_power.Right = {};
    average_power_time.Right = {};
    average_phase.Right = {};
    average_phase_time.Right = {};
    range_power.Right = {};
    range_power_time.Right = {};
    range_phase.Right = {};
    range_phase_time.Right = {};
    [~,band_names] = getFreqBandInd(signal_analysis_data.Right.Freq_Values{1});
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        band_inds = getFreqBandInd(signal_analysis_data.Right.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Right,aligned_data.gait_events.Properties.VariableNames{j})
                average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
            end
            vals_power = [];
            vals_power_time = [];
            vals_phase = [];
            vals_phase_time = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signal_analysis_data.Right.Time{i}-(event_time-pre_post_time)));
                    [~,min_ind_post] = min(abs(signal_analysis_data.Right.Time{i}-(event_time+pre_post_time)));
                    [~,event_ind] = min(abs(signal_analysis_data.Right.Time{i}-event_time));
                    
                    if isfield(signal_analysis_data.Right,'PSD')
                        temp = 10*log10(abs(signal_analysis_data.Right.Values{i}(:,event_ind)));
                        temp1 = 10*log10(abs(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signal_analysis_data.Right.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post));
                    else
                        temp = 10*log10(abs(signal_analysis_data.Right.Values{i}(:,event_ind)));
                        temp1 = abs(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signal_analysis_data.Right.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post));
                    end
                    
                    if size(temp1,2) == round((pre_post_time*2)/time_res_right)+1
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_power_time(count,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power,1));
            average_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power_time,1));
            average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            average_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase_time,1));
            
            if strcmp(spread_type,'SE')
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1))./size(vals_power,1);
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1))./size(vals_power_time,1);
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./size(vals_phase,1);
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1))./size(vals_phase_time,1);
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

if strcmp(plot_type,'all') || strcmp(plot_type,'power')
    if isfield(average_power_time,'Left')
        time_vec = -1:time_res_left:1;
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_power_time.Left.(event_compare{j}{m}){i}(:,k),range_power_time.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                        %                     mean_std_plot(time_vec,average_phase.Left.(event_compare{j}{m}){i}(:,k),std_phase.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Left,'PSD')
                        ylabel('db');
                    else
                        ylabel('Magnitude');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(average_power_time,'Right')
        time_vec = -1:time_res_right:1;
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_power_time.Right.(event_compare{j}{m}){i}(:,k),range_power_time.Right.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Right,'PSD')
                        ylabel('db');
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

if strcmp(plot_type,'all') || strcmp(plot_type,'phase')
    if isfield(average_phase_time,'Left')
        time_vec = -1:time_res_left:1;
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_phase_time.Left.(event_compare{j}{m}){i}(:,k),range_phase_time.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Left,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(average_phase_time,'Right')
        time_vec = -1:time_res_right:1;
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_phase_time.Right.(event_compare{j}{m}){i}(:,k),range_phase_time.Right.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Right,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

if strcmp(plot_type,'all') || strcmp(plot_type,'power_boxplot')
    if isfield(average_phase_time,'Left')
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(event_compare{j})
%                         box_vals = [box_vals;average_power(
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Left,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    Z = average_phase.Left.LHS{1};
    Y = Z(:);
    X = transpose(repelem(band_names,size(Z,1)));
    title('Left LHS');
    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
    set(gca,'TickLabelInterpreter','tex')
end

if strcmp(plot_type,'all') || strcmp(plot_type,'phase_boxplot')
end

%% Save plots
if save_flag
    save_dir = uigetdir();
    
    figure_format(12,8,12);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgEventPower'))
        mkdir(fullfile(save_dir,'AvgEventPower'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition))
        mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT'))
            mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT'))
            mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{n}));
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
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
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