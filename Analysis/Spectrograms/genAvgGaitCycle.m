function genAvgGaitCycle(aligned_data,signal_analysis_data,cycle_start_event,n_percent_bins,subjectID,save_flag)
if ~exist('cycle_start_event','var') || isempty(cycle_start_event)
    cycle_start_event = 'LHS';
end

if ~exist('n_percent_bins','var') || isempty(n_percent_bins)
    n_percent_bins = 100;
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var') || isempty(save_flag)
    save_flag = 0;
end 

%% Average gait cycle value
% re-sort gait events with the starting event
gait_events_sorted = sort_gait_events(aligned_data.gait_events,cycle_start_event);

% Go through all valid gait cycles
% Left
if isfield(signal_analysis_data,'Left')
    if isfield(signal_analysis_data.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    gait_cycle_avg.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        gait_cycle_mat_left = zeros(length(signal_analysis_data.Left.Freq_Values{i}),n_percent_bins,1);
        count = 1;
        for j = 1:height(gait_events_sorted)-1
            if ~isnan(gait_events_sorted.(cycle_start_event)(j)) && ~isnan(gait_events_sorted.(cycle_start_event)(j+1)) && (diff(gait_events_sorted.(cycle_start_event)(j:j+1)) < 2)
                [~,start_ind] = min(abs(signal_analysis_data.Left.Time{i}-gait_events_sorted.(cycle_start_event)(j)));
                [~,end_ind] = min(abs(signal_analysis_data.Left.Time{i}-gait_events_sorted.(cycle_start_event)(j+1)));
                
                if isfield(signal_analysis_data.Left,'PSD')
                    data_snip = 10*log10(abs(signal_analysis_data.Left.PSD{i}(:,start_ind:end_ind)));
                else
                    data_snip = abs(signal_analysis_data.Left.Values{i}(:,start_ind:end_ind));
                end
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
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
if isfield(signal_analysis_data,'Right')
    if isfield(signal_analysis_data.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    gait_cycle_avg.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        gait_cycle_mat_right = zeros(length(signal_analysis_data.Right.Freq_Values{i}),n_percent_bins,1);
        count = 1;
        for j = 1:height(gait_events_sorted)-1
            if ~isnan(gait_events_sorted.(cycle_start_event)(j)) && ~isnan(gait_events_sorted.(cycle_start_event)(j+1)) && (diff(gait_events_sorted.(cycle_start_event)(j:j+1)) < 2)
                [~,start_ind] = min(abs(signal_analysis_data.Right.Time{i}-gait_events_sorted.(cycle_start_event)(j)));
                [~,end_ind] = min(abs(signal_analysis_data.Right.Time{i}-gait_events_sorted.(cycle_start_event)(j+1)));
                if isfield(signal_analysis_data.Right,'PSD')
                    data_snip = 10*log10(abs(signal_analysis_data.Right.PSD{i}(:,start_ind:end_ind)));
                else
                    data_snip = abs(signal_analysis_data.Right.Values{i}(:,start_ind:end_ind));
                end
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
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

%% Plot
fig_vec = [];
for i = 1:length(signal_analysis_data.Left.Chan_Names)
    fig_vec(end+1) = figure;
    if isfield(signal_analysis_data.Left,'PSD')
        pcolor(1:100,signal_analysis_data.Left.Freq_Values{i},gait_cycle_avg.Left{i});
        ylim([2.5,50]);
    else
        ax = pcolor(1:100,log2(signal_analysis_data.Left.Freq_Values{i}),gait_cycle_avg.Left{i});
        ticks = logspace(log10(2.5),log10(50),10);
        ax.Parent.YTick = log2(ticks);
        ax.Parent.YTickLabel = ticks;
        ylim([log2(2.5),log2(50)]);
    end
    shading interp;
    title({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i}});
    xlabel('Gait Cycle %');
    ylabel('Frequency (Hz)');
end

for i = 1:length(signal_analysis_data.Right.Chan_Names)
    fig_vec(end+1) = figure;
    if isfield(signal_analysis_data.Right,'PSD')
        pcolor(1:100,signal_analysis_data.Right.Freq_Values{i},gait_cycle_avg.Right{i});
        ylim([2.5,50]);
    else
        ax = pcolor(1:100,log2(signal_analysis_data.Right.Freq_Values{i}),gait_cycle_avg.Right{i});
        ticks = logspace(log10(2.5),log10(50),10);
        ax.Parent.YTick = log2(ticks);
        ax.Parent.YTickLabel = ticks;
        ylim([log2(2.5),log2(50)]);
    end
    shading interp;
    title({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i}});
    xlabel('Gait Cycle %');
    ylabel('Frequency (Hz)');
end

%% Save plots
if save_flag
    save_dir = uigetdir();
    
    figure_format(6,6,10);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec'))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'FT'))
            mkdir(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'CWT'))
            mkdir(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'CWT',folders_to_check{n}));
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
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'FT',folders_to_check{1},strrep(save_name,' ','_')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        end
        
        %         for k = 2:length(folders_to_check)
        %             if strcmp(analysis_type,'FT')
        %                 print(fig_vec(i),[fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'FT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
        %             elseif strcmp(analysis_type,'CWT')
        %                 print(fig_vec(i),[fullfile(save_dir,'AvgGaitCycleSpec',aligned_data.stim_condition,'CWT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
        %             end
        %         end
    end
end
end

function sorted_gait_events = sort_gait_events(gait_events,cycle_start_event)
if strcmp(gait_events.Properties.VariableNames{1},cycle_start_event)
    sorted_gait_events = gait_events;
else
    gait_event_order = [];
    switch cycle_start_event
        case 'LHS'
            gait_event_order = {'LHS','RTO','RHS','LTO'};
        case 'LTO'
            gait_event_order = {'LTO','LHS','RTO','RHS'};
        case 'RHS'
            gait_event_order = {'RHS','LTO','LHS','RTO'};
        case 'RTO'
            gait_event_order = {'RTO','RHS','LTO','LHS'};
    end
    
    shift_ind = find(cellfun(@(x) strcmp(x,cycle_start_event),gait_events.Properties.VariableNames))-1;
    
    sorted_gait_events = nan(1+height(gait_events),4);
    for i = 1:length(gait_event_order)-shift_ind
        sorted_gait_events(2:end,i) = gait_events.(gait_event_order{i});
    end
    
    for j = length(gait_event_order)-(shift_ind-1):length(gait_event_order)
        sorted_gait_events(1:end-1,j) = gait_events.(gait_event_order{j});
    end
    
    sorted_gait_events = array2table(sorted_gait_events,'VariableNames',gait_event_order);
end
end