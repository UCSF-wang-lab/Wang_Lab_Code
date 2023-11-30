function singleGaitCycleRastergram(file_names,varargin)

%% go through optional inputs
for i = 1:2:nargin-1
    switch varargin{i}
        case 'analysisType'
            analysis_type = varargin{i+1};
        case 'normalizationType'
            normalizationType = varargin{i+1};    
        case 'turn_threshold'
            turn_threshold = varargin{i+1};
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

if ~exist('analysis_type','var')
    analysis_type = 'CWT';
end

if ~exist('normalizationType','var')
    normalizationType = 'zscore';
end

if ~exist('turn_threshold','var')
    turn_threshold = 30;
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

if ~exist('savePlot','var')
    savePlot = false;
end

% colors for gait events
marker_colors = CBMap('GaitEvents');

% set frequency bands and channel names
freq_bands = [4,8;...       % theta
    8,13;...      % alpha
    13,30;...     % beta
    13,20;...     % low beta
    20,30;...     % high beta
    30,50];       % low gamma
freq_bands_names = {'Theta','Alpha','Beta','Low Beta','High Beta','Low Gamma','Custom'};
chan_names = {'+2-0','+3-1','+9-8','+11-10'};

% varibles to hold all the data from all files passed in
left_single_gait_cycle_cell_all = [];
right_single_gait_cycle_cell_all = [];

left_RTO_relative_time = [];
left_RHS_relative_time = [];
left_LTO_relative_time = [];
right_RTO_relative_time = [];
right_RHS_relative_time = [];
right_LTO_relative_time = [];

%% Loop through files passed in
for f = 1:length(file_names)
    % load in data and calculate cwt
    load(file_names{f});

    % Calculate spectral data
    if strcmp(analysis_type,'FT')
        lfp_spec = calcRCS_STFT(aligned_data,[],1,0.9,[]);
    elseif strcmp(analysis_type,'CWT')
        lfp_spec = calcRCS_CWT(aligned_data);
    end


    % detect and remove turns from inclusion
    gait_events_turns_removed = removeGaitCyclesTurns(aligned_data.Xsens,aligned_data.gait_events,turn_threshold);

    % sort gait events and extract valid gait cycles
    gait_events_ordered = sortGaitEvents(gait_events_turns_removed,'LHS');


    % extract power data
    if isfield(lfp_spec,'Left')
        % determine valid gait cycles
        gc_start_search = find(gait_events_ordered.LHS>=lfp_spec.Left.Time{1}(1),1,'first');
        gc_end_search = find(gait_events_ordered.LHS<=lfp_spec.Left.Time{1}(end),1,'last');
        gait_events_ordered_trim = gait_events_ordered(gc_start_search:gc_end_search,:);

        gc_test = ~isnan(gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LHS(1:end-1));
        gc_test2 = (gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LHS(1:end-1)) < 1.8;
        valid_gc = gc_test&gc_test2;
        left_single_gait_cycle_cell = cell(size(freq_bands,1),sum(valid_gc),4);   % freq_band x gait cycle x recording area

        count = 1;
        for j = 1:length(valid_gc)
            if valid_gc(j)
                for k = 1:length(lfp_spec.Left.Chan_Names)
                    time_inds = getTimeInds(lfp_spec.Left.Time{k},gait_events_ordered_trim.LHS(j),gait_events_ordered_trim.LHS(j+1));
                    for m = 1:size(freq_bands,1)
                        freq_inds = getFreqInds(lfp_spec.Left.Freq_Values{k},freq_bands(m,1),freq_bands(m,2));
                        gc_data = lfp_spec.Left.Values{k}([freq_inds(1):freq_inds(2)],[time_inds(1):time_inds(2)]);
                        left_single_gait_cycle_cell{m,count,k} = mean(abs(gc_data));
                    end
                end
                count = count + 1;
            end
        end

        % normalize data
        left_single_gait_cycle_cell_norm = normalizeSingleGaitCycle(left_single_gait_cycle_cell,lfp_spec.Left,gait_events_ordered_trim,freq_bands,normalizationType);
        left_single_gait_cycle_cell_all = [left_single_gait_cycle_cell_all,left_single_gait_cycle_cell_norm];

        % calc event time markers from LHS
        left_RTO_relative_time = [left_RTO_relative_time;gait_events_ordered_trim.('RTO')(valid_gc)-gait_events_ordered_trim.LHS(valid_gc)];
        left_RHS_relative_time = [left_RHS_relative_time;gait_events_ordered_trim.('RHS')(valid_gc)-gait_events_ordered_trim.LHS(valid_gc)];
        left_LTO_relative_time = [left_LTO_relative_time;gait_events_ordered_trim.('LTO')(valid_gc)-gait_events_ordered_trim.LHS(valid_gc)];
    end

    if isfield(lfp_spec,'Right')
        % determine valid gait cycles
        gc_start_search = find(gait_events_ordered.LHS>=lfp_spec.Right.Time{1}(1),1,'first');
        gc_end_search = find(gait_events_ordered.LHS<=lfp_spec.Right.Time{1}(end),1,'last');
        gait_events_ordered_trim = gait_events_ordered(gc_start_search:gc_end_search,:);

        gc_test = ~isnan(gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LHS(1:end-1));
        gc_test2 = (gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LHS(1:end-1)) < 1.8;
        valid_gc = gc_test&gc_test2;

        right_single_gait_cycle_cell = cell(size(freq_bands,1),sum(valid_gc),4);   % freq_band x gait cycle x recording area
        
        count = 1;
        for j = 1:length(valid_gc)
            if valid_gc(j)
                for k = 1:length(lfp_spec.Right.Chan_Names)
                    time_inds = getTimeInds(lfp_spec.Right.Time{k},gait_events_ordered_trim.LHS(j),gait_events_ordered_trim.LHS(j+1));
                    for m = 1:size(freq_bands,1)
                        freq_inds = getFreqInds(lfp_spec.Right.Freq_Values{k},freq_bands(m,1),freq_bands(m,2));
                        gc_data = lfp_spec.Right.Values{k}([freq_inds(1):freq_inds(2)],[time_inds(1):time_inds(2)]);
                        right_single_gait_cycle_cell{m,count,k} = mean(abs(gc_data));
                    end
                end
                count = count + 1;
            end
        end
        % normalize data
        right_single_gait_cycle_cell_norm = normalizeSingleGaitCycle(right_single_gait_cycle_cell,lfp_spec.Right,gait_events_ordered_trim,freq_bands,normalizationType);
        right_single_gait_cycle_cell_all = [right_single_gait_cycle_cell_all,right_single_gait_cycle_cell_norm];

        % calc event time markers from LHS
        right_RTO_relative_time = [right_RTO_relative_time;gait_events_ordered_trim.('RTO')(valid_gc)-gait_events_ordered_trim.LHS(valid_gc)];
        right_RHS_relative_time = [right_RHS_relative_time;gait_events_ordered_trim.('RHS')(valid_gc)-gait_events_ordered_trim.LHS(valid_gc)];
        right_LTO_relative_time = [right_LTO_relative_time;gait_events_ordered_trim.('LTO')(valid_gc)-gait_events_ordered_trim.LHS(valid_gc)];
    end
end

%% Determine which gait cycle is the longest and setup visualization matrix
if isfield(lfp_spec,'Left')
    temp = left_single_gait_cycle_cell_all(1,:,1);
    gc_lengths = cellfun(@(x)size(x,2),temp);
    max_gc_length = max(gc_lengths);
    [gc_sort_lengths,gc_sort_inds] = sort(gc_lengths,'descend');
    left_single_gait_cycle_mat = nan(length(gc_lengths),max_gc_length,length(freq_bands_names),4);       % gait cycle x time x freq band x channel

    for i = 1:length(lfp_spec.Left.Chan_Names)
        for j = 1:length(freq_bands_names)
            for k = 1:length(gc_sort_inds)
                left_single_gait_cycle_mat(k,1:gc_sort_lengths(k),j,i) = left_single_gait_cycle_cell_all{j,gc_sort_inds(k),i};
            end
        end
    end

    left_RTO_relative_time = left_RTO_relative_time(gc_sort_inds);
    left_RHS_relative_time = left_RHS_relative_time(gc_sort_inds);
    left_LTO_relative_time = left_LTO_relative_time(gc_sort_inds);

    left_single_gait_cycle_struct = removeArtifactGC(left_single_gait_cycle_mat,left_RTO_relative_time,left_RHS_relative_time,left_LTO_relative_time);
    left_sr = unique(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate);
    left_single_gait_cycle_time_vec = (0:max_gc_length-1).*(1/left_sr);
end

if isfield(lfp_spec,'Right')
    temp = right_single_gait_cycle_cell_all(1,:,1);
    gc_lengths = cellfun(@(x)size(x,2),temp);
    max_gc_length = max(gc_lengths);
    [gc_sort_lengths,gc_sort_inds] = sort(gc_lengths,'descend');
    right_single_gait_cycle_mat = nan(length(gc_lengths),max_gc_length,length(freq_bands_names),4);       % gait cycle x time x freq band x channel

    for i = 1:length(lfp_spec.Right.Chan_Names)
        for j = 1:length(freq_bands_names)
            for k = 1:length(gc_sort_inds)
                right_single_gait_cycle_mat(k,1:gc_sort_lengths(k),j,i) = right_single_gait_cycle_cell_all{j,gc_sort_inds(k),i};
            end
        end
    end

    right_RTO_relative_time = right_RTO_relative_time(gc_sort_inds);
    right_RHS_relative_time = right_RHS_relative_time(gc_sort_inds);
    right_LTO_relative_time = right_LTO_relative_time(gc_sort_inds);

    right_single_gait_cycle_struct = removeArtifactGC(right_single_gait_cycle_mat,right_RTO_relative_time,right_RHS_relative_time,right_LTO_relative_time);
    right_sr = unique(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate);
    right_single_gait_cycle_time_vec = (0:max_gc_length-1).*(1/right_sr);
end

%% Plot the data
fig_vec = [];
if exist('left_single_gait_cycle_struct','var')
    for i = 1:size(left_single_gait_cycle_struct.data,2)        % channel
        for j = 1:size(left_single_gait_cycle_struct.data,1)    % freq band
            fig_vec(end+1) = figure();
            pcolor(left_single_gait_cycle_time_vec,size(left_single_gait_cycle_struct.data{j,i},1):-1:1,left_single_gait_cycle_struct.data{j,i});
            title({subjectID;['Left ',chan_names{i},' ',freq_bands_names{j}]});
            shading flat
            hold on;
            a(1) = scatter(left_single_gait_cycle_struct.RTO_timings{j,i},length(left_single_gait_cycle_struct.RTO_timings{j,i}):-1:1,20,'MarkerEdgeColor',marker_colors.RTO,'MarkerFaceColor',marker_colors.RTO,'DisplayName','RTO');
            a(2) = scatter(left_single_gait_cycle_struct.RHS_timings{j,i},length(left_single_gait_cycle_struct.RHS_timings{j,i}):-1:1,20,'MarkerEdgeColor',marker_colors.RHS,'MarkerFaceColor',marker_colors.RHS,'DisplayName','RHS');
            a(3) = scatter(left_single_gait_cycle_struct.LTO_timings{j,i},length(left_single_gait_cycle_struct.LTO_timings{j,i}):-1:1,20,'MarkerEdgeColor',marker_colors.LTO,'MarkerFaceColor',marker_colors.LTO,'DisplayName','LTO');

            if mod(size(left_single_gait_cycle_struct.data{j,i},1)-1,5)~=0
                yticks([1:5:size(left_single_gait_cycle_struct.data{j,i},1),size(left_single_gait_cycle_struct.data{j,i},1)]);
            else
                yticks([1:5:size(left_single_gait_cycle_struct.data{j,i},1)])
            end

            ylabel('Gait Cycle');
            xlabel('Time (sec)');
            colorbar;
            caxis([-3,3])
        end
    end
end

if exist('right_single_gait_cycle_struct','var')
    for i = 1:size(right_single_gait_cycle_struct.data,2)       % channel
        for j = 1:size(right_single_gait_cycle_struct.data,1)   % freq band
            fig_vec(end+1) = figure();
            pcolor(right_single_gait_cycle_time_vec,size(right_single_gait_cycle_struct.data{j,i},1):-1:1,right_single_gait_cycle_struct.data{j,i});
            title({subjectID;['Right ',chan_names{i},' ',freq_bands_names{j}]});
            shading flat
            hold on;
            a(1) = scatter(right_single_gait_cycle_struct.RTO_timings{j,i},length(right_single_gait_cycle_struct.RTO_timings{j,i}):-1:1,20,'MarkerEdgeColor',marker_colors.RTO,'MarkerFaceColor',marker_colors.RTO,'DisplayName','RTO');
            a(2) = scatter(right_single_gait_cycle_struct.RHS_timings{j,i},length(right_single_gait_cycle_struct.RHS_timings{j,i}):-1:1,20,'MarkerEdgeColor',marker_colors.RHS,'MarkerFaceColor',marker_colors.RHS,'DisplayName','RHS');
            a(3) = scatter(right_single_gait_cycle_struct.LTO_timings{j,i},length(right_single_gait_cycle_struct.LTO_timings{j,i}):-1:1,20,'MarkerEdgeColor',marker_colors.LTO,'MarkerFaceColor',marker_colors.LTO,'DisplayName','LTO');

            if mod(size(right_single_gait_cycle_struct.data{j,i},1)-1,5)~=0
                yticks([1:5:size(right_single_gait_cycle_struct.data{j,i},1),size(right_single_gait_cycle_struct.data{j,i},1)]);
            else
                yticks([1:5:size(right_single_gait_cycle_struct.data{j,i},1)])
            end

            ylabel('Gait Cycle');
            xlabel('Time (sec)');
            colorbar;
            caxis([-3,3])
        end
    end
end

%% Save plots
if savePlot
    save_dir = uigetdir();

    figure_format(6,6,10);

    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram'))
        mkdir(fullfile(save_dir,'singleGaitCycleRastergram'));
    end

    if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM']))
    end

    if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED']))
        mkdir(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED']))
    end

    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT'))
            mkdir(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT'))
            mkdir(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT'))
        end
    end

    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{n}));
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

        if length(file_names) == 1
            if isfield(aligned_data,'trial_num') && ~isempty(aligned_data.trial_num)
                save_name = [save_name,' ',sprintf('Trial%i',aligned_data.trial_num)];
            end
        end

        if ~strcmp(normalizationType,'none')
            save_name = [save_name,' ',normalizationType];
        end

        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{1},strrep(save_name,' ','_')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        end

        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'FT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'singleGaitCycleRastergram',[aligned_data.stim_condition,'_STIM'],[aligned_data.med_condition,'_MED'],'CWT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end

end

function time_indexes = getTimeInds(time_vec,start_time,end_time)
[~,start_ind] = min(abs(time_vec-start_time));
[~,end_ind] = min(abs(time_vec-end_time));
time_indexes = sort([start_ind,end_ind]);
end

function freq_inds = getFreqInds(freq_vec,start_freq,end_freq)
[~,first_freq_ind] = min(abs(freq_vec-start_freq));
[~,second_freq_ind] = min(abs(freq_vec-end_freq));
freq_inds = sort([first_freq_ind,second_freq_ind]);
end

function norm_data = normalizeSingleGaitCycle(gait_cycle_cell,full_spec,gait_events,freq_bands,norm_type)

% check to make sure gait cycle cell has the same number of channels
if size(gait_cycle_cell,3) ~= size(full_spec.Values,2)
    error('Mismatched number of channels');
end

% check to make sure gait cycle cell has the same number of frequency bands
if size(gait_cycle_cell,1) ~= size(freq_bands,1)
    error('Mismatched number of frequency bands analyzed');
end

norm_data = cell(size(gait_cycle_cell));


% Get normalization info
start_mvnt = min(gait_events{1,:});
end_mvnt = max(gait_events{end,:});
time_inds = getTimeInds(full_spec.Time{1},start_mvnt,end_mvnt);

normalization_mean = nan(size(freq_bands,1),size(gait_cycle_cell,3));
normalization_sd = nan(size(freq_bands,1),size(gait_cycle_cell,3));
normalization_var = nan(size(freq_bands,1),size(gait_cycle_cell,3));

for i = 1:size(freq_bands,1)
    for j = 1:size(gait_cycle_cell,3)
        freq_inds = getFreqInds(full_spec.Freq_Values{j},freq_bands(i,1),freq_bands(i,2));

        % Go through movement data and remove artifacts
        data_snip = abs(full_spec.Values{j}(freq_inds,time_inds(1):time_inds(2)));
        temp_mean = mean(data_snip,'all');
        temp_sd = std(mean(data_snip));
        remove_inds = or(mean(data_snip) > (temp_mean + 2*temp_sd),mean(data_snip) < (temp_mean - 2*temp_sd));
        data_clean = mean(data_snip);
        data_clean(remove_inds) = [];

        normalization_mean(i,j) = mean(data_clean);
        normalization_sd(i,j) = std(data_clean);
        normalization_var(i,j) = var(data_clean);
    end
end

% Normalize the data
for i = 1:size(gait_cycle_cell,1)
    for j = 1:size(gait_cycle_cell,3)
        for k = 1:size(gait_cycle_cell,2)
            if strcmp(norm_type,'zscore')
                norm_data{i,k,j} = (gait_cycle_cell{i,k,j}-normalization_mean(i,j))/normalization_sd(i,j);
            end
        end
    end
end
end

function data_artifact_removed_struct = removeArtifactGC(data,RTO_timings,RHS_timings,LTO_timings)
data_artifact_removed = cell(size(data,3),size(data,4));
RTO_timings_artifact_removed = cell(size(data,3),size(data,4));
RHS_timings_artifact_removed = cell(size(data,3),size(data,4));
LTO_timings_artifact_removed = cell(size(data,3),size(data,4));

remove_mat = zeros(size(data,1),size(data,3),size(data,4));
for i = 1:size(data,4)          % channel
    for j = 1:size(data,3)      % frequency band
        remove_ind = zeros(size(data,1),1);
        for k = 1:size(data,1)  % gait cycle
            if any(data(k,:,j,i)>8)
                remove_ind(k,1) = 1;
            end
        end
        data_artifact_removed{j,i} = data(~remove_ind,:,j,i);
        RTO_timings_artifact_removed{j,i} = RTO_timings(~remove_ind);
        RHS_timings_artifact_removed{j,i} = RHS_timings(~remove_ind);
        LTO_timings_artifact_removed{j,i} = LTO_timings(~remove_ind);
    end
end

data_artifact_removed_struct.data = data_artifact_removed;
data_artifact_removed_struct.RTO_timings = RTO_timings_artifact_removed;
data_artifact_removed_struct.RHS_timings = RHS_timings_artifact_removed;
data_artifact_removed_struct.LTO_timings = LTO_timings_artifact_removed;

end