function genAvgGaitCycleV2(aligned_data,signal_analysis_data,varargin)
% cycle_start_event,n_percent_bins,baseline_data,subjectID,save_flag)

for i = 1:2:nargin-2
    switch varargin{i}
        case 'cycle_start_event'
            cycle_start_event = varargin{i+1};
        case 'n_percent_bins'
            n_percent_bins = varargin{i+1};
        case 'normalize_by'
            normalize_by = varargin{i+1};   % Can be none|baseline|average_during_walking
        case 'normalization_type'
            normalization_type = varargin{i+1}; % Can be none|percent_change|zscore
        case 'baseline_time'
            baseline_time = varargin{i+1};  % Only valid if "normalization_by" is set to "baseline"
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'save_flag'
            save_flag = varargin{i+1};
    end
end

if ~exist('cycle_start_event','var')
    cycle_start_event = 'LHS';
end

if ~exist('n_percent_bins','var')
    n_percent_bins = 100;
end

if ~exist('baseline_time','var')
    baseline_time = [];
end

if ~exist('normalize_by','var')
    normalize_by = 'none';
end

if ~exist('normalization_type','var')
    normalization_type = 'none';
end

if strcmp(normalize_by,'none') && ~strcmp(normalization_type,'none')
    warning('No input for normalize_by. Defauling to average during walking.');
    normalize_by = 'average_during_walking';
end

if ~strcmp(normalize_by,'none') && strcmp(normalization_type,'none')
    warning('No input for normalization_type. Defauling to zscore.');
    normalization_type = 'zscore';
end

if strcmp(normalize_by,'baseline') && ~exist('baseline_time','var')
    error('Normalization by baseline, but no baseline data was passed in.');
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var')
    save_flag = 0;
end

%% Average gait cycle value
% re-sort gait events with the starting event
gait_events_sorted = sortGaitEvents(aligned_data.gait_events,cycle_start_event);

% Go through all valid gait cycles
% Left
if isfield(signal_analysis_data,'Left')
    if isfield(signal_analysis_data.Left,'PSD')
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
                    data_snip = 20*log10(abs(signal_analysis_data.Left.PSD{i}(:,start_ind:end_ind)));
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
                    data_snip = 20*log10(abs(signal_analysis_data.Right.PSD{i}(:,start_ind:end_ind)));
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

%% Normalize if set
if ~strcmp(normalize_by,'none')
    normalization = [];
    if strcmp(normalize_by,'average_during_walking')
        if isfield(signal_analysis_data,'Left')
            normalization.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
            walking_start_ind = find(signal_analysis_data.Left.Time{i} >= min(gait_events_sorted{1,:})-1,1,'first');
            walking_end_ind = find(signal_analysis_data.Left.Time{i} <= max(gait_events_sorted{end,:}),1,'last');
            for i = 1:length(signal_analysis_data.Left.Chan_Names)
                normalization.Left{i} = signal_analysis_data.Left.Values{i}(:,walking_start_ind:walking_end_ind);
            end
        end
        if isfield(signal_analysis_data,'Right')
            normalization.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
            walking_start_ind = find(signal_analysis_data.Right.Time{1} >= min(gait_events_sorted{1,:})-1,1,'first');
            walking_end_ind = find(signal_analysis_data.Right.Time{1} <= max(gait_events_sorted{end,:}),1,'last');
            for i = 1:length(signal_analysis_data.Right.Chan_Names)
                normalization.Right{i} = signal_analysis_data.Right.Values{i}(:,walking_start_ind:walking_end_ind);
            end
        end
    elseif strcmp(normalize_by,'baseline')
        if isfield(signal_analysis_data,'Left')
            normalization.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
            baseline_start_ind = find(signal_analysis_data.Left.Time{i} >= baseline_time(1),1,'first');
            baseline_end_ind = find(signal_analysis_data.Left.Time{i} <= baseline_time(2),1,'last');
            for i = 1:length(signal_analysis_data.Left.Chan_Names)
                normalization.Left{i} = signal_analysis_data.Left.Values{i}(:,baseline_start_ind:baseline_end_ind);
            end
        end
        if isfield(signal_analysis_data,'Right')
            normalization.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
            baseline_start_ind = find(signal_analysis_data.Right.Time{i} >= baseline_time(1),1,'first');
            baseline_end_ind = find(signal_analysis_data.Right.Time{i} <= baseline_time(2),1,'last');
            for i = 1:length(signal_analysis_data.Right.Chan_Names)
                normalization.Right{i} = signal_analysis_data.Right.Values{i}(:,baseline_start_ind:baseline_end_ind);
            end
        end
    end
end

if strcmp(normalization_type,'percent_change')
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
elseif strcmp(normalization_type,'zscore')
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
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signal_analysis_data.Left,'PSD')
            if sum(isnan(gait_cycle_avg.Left{i}),'all') == 0
                pcolor(1:100,signal_analysis_data.Left.Freq_Values{i},gait_cycle_avg.Left{i});
                ylim([2.5,50]);
            end
        else
            if sum(isnan(gait_cycle_avg.Left{i}),'all') == 0
                ax = pcolor(1:100,log2(signal_analysis_data.Left.Freq_Values{i}),gait_cycle_avg.Left{i});
                ticks = logspace(log10(2.5),log10(50),10);
                ax.Parent.YTick = log2(ticks);
                ax.Parent.YTickLabel = ticks;
                ylim([log2(2.5),log2(50)]);
                if strcmp(normalize_by,'baseline') && strcmp(normalization_type,'zscore')
                    caxis([-2,2]);
                end
            end
        end
        shading interp;
        title({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i}});
        xlabel('Gait Cycle %');
        ylabel('Frequency (Hz)');
    end
end

if isfield(gait_cycle_avg,'Right')
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signal_analysis_data.Right,'PSD')
            if sum(isnan(gait_cycle_avg.Right{i}),'all') == 0
                pcolor(1:100,signal_analysis_data.Right.Freq_Values{i},gait_cycle_avg.Right{i});
                ylim([2.5,50]);
            end
        else
            if sum(isnan(gait_cycle_avg.Right{i}),'all') == 0
                ax = pcolor(1:100,log2(signal_analysis_data.Right.Freq_Values{i}),gait_cycle_avg.Right{i});
                ticks = logspace(log10(2.5),log10(50),10);
                ax.Parent.YTick = log2(ticks);
                ax.Parent.YTickLabel = ticks;
                ylim([log2(2.5),log2(50)]);
                if strcmp(normalize_by,'baseline') && strcmp(normalization_type,'zscore')
                    caxis([-2,2]);
                end
            end
        end
        shading interp;
        title({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i}});
        xlabel('Gait Cycle %');
        ylabel('Frequency (Hz)');
    end
end

%% Save plots
if save_flag
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
        
        if ~strcmp(normalization_type,'none')
            save_name = [save_name,' ',normalization_type];
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