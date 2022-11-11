function singleGaitCycleRastergram(file_name,varargin)

% go through optional inputs
for i = 1:2:nargin-2
    switch varargin{i}
        case 'gait_cycle_ignore'
            sorted_gc_ignore = varargin{i+1};
    end
end

if ~exist('sorted_gc_ignore','var')
    sorted_gc_ignore = [];
end


% load in data and calculate cwt
load(file_name);
lfp_cwt = calcRCS_CWT(aligned_data);

% set frequency bands and channel names
freq_bands = [1,4;...       % delta
              4,8;...       % theta
              8,13;...      % alpha
              13,30;...     % beta
              30,50];       % low gamma
freq_bands_names = {'Delta','Theta','Alpha','Beta','Low Gamma'};
chan_names = {'+2-0','+3-1','+9-8','+11-10'};

% sort gait events and extract valid gait cycles
gait_events_sorted = sortGaitEvents(aligned_data.gait_events,'LHS');

gc_test = ~isnan(gait_events_sorted.LHS(2:end)-gait_events_sorted.LHS(1:end-1));
gc_test2 = (gait_events_sorted.LHS(2:end)-gait_events_sorted.LHS(1:end-1)) < 1.8;
valid_gc = and(gc_test,gc_test2);
left_single_gait_cycle_cell = cell(size(freq_bands,1),sum(valid_gc),4);   % freq_band x gait cycle x recording area
right_single_gait_cycle_cell = cell(size(freq_bands,1),sum(valid_gc),4);   % freq_band x gait cycle x recording area

% extract power data
if isfield(lfp_cwt,'Left')
    count = 1;
    for j = 1:length(valid_gc)
        if valid_gc(j)
            for k = 1:length(lfp_cwt.Left.Chan_Names)
                time_inds = getTimeInds(lfp_cwt.Left.Time{k},gait_events_sorted.LHS(j),gait_events_sorted.LHS(j+1));
                for m = 1:size(freq_bands,1)
                    freq_inds = getFreqInds(lfp_cwt.Left.Freq_Values{k},freq_bands(m,1),freq_bands(m,2));
                    gc_data = lfp_cwt.Left.Values{k}([freq_inds(1):freq_inds(2)],[time_inds(1):time_inds(2)]);
                    left_single_gait_cycle_cell{m,count,k} = mean(abs(gc_data));
                end
            end
            count = count + 1;
        end
    end

%     temp = {left_single_gait_cycle_cell{1,:,1}};
%     gc_lengths = cellfun(@(x)size(x,2),temp);
%     max_gc_length = max(gc_lengths);
%     [gc_sort_lengths,gc_sort_inds] = sort(gc_lengths,'descend');

    % normalize data
    left_single_gait_cycle_cell_norm = normalizeSingleGaitCycle(left_single_gait_cycle_cell,lfp_cwt.Left,gait_events_sorted,freq_bands,'zscore');
end

if isfield(lfp_cwt,'Right')
    count = 1;
    for j = 1:length(valid_gc)
        if valid_gc(j)
            for k = 1:length(lfp_cwt.Right.Chan_Names)
                time_inds = getTimeInds(lfp_cwt.Right.Time{k},gait_events_sorted.LHS(j),gait_events_sorted.LHS(j+1));
                for m = 1:size(freq_bands,1)
                    freq_inds = getFreqInds(lfp_cwt.Right.Freq_Values{k},freq_bands(m,1),freq_bands(m,2));
                    gc_data = lfp_cwt.Right.Values{k}([freq_inds(1):freq_inds(2)],[time_inds(1):time_inds(2)]);
                    right_single_gait_cycle_cell{m,count,k} = mean(abs(gc_data));
                end
            end
            count = count + 1;
        end
    end
    % normalize data
    right_single_gait_cycle_cell_norm = normalizeSingleGaitCycle(right_single_gait_cycle_cell,lfp_cwt.Right,gait_events_sorted,freq_bands,'zscore');
end

% Determine which gait cycle is the longest and setup visualization matrix
if isfield(lfp_cwt,'Left')
    temp = {left_single_gait_cycle_cell{1,:,1}};
    gc_lengths = cellfun(@(x)size(x,2),temp);
    max_gc_length = max(gc_lengths);
    [gc_sort_lengths,gc_sort_inds] = sort(gc_lengths,'descend');
    left_single_gait_cycle_mat = nan(sum(valid_gc),max_gc_length,length(freq_bands_names),4);       % gait cycle x time x freq band x channel

    for i = 1:length(lfp_cwt.Left.Chan_Names)
        for j = 1:length(freq_bands_names)
            for k = 1:length(gc_sort_inds)
                left_single_gait_cycle_mat(k,1:gc_sort_lengths(k),j,i) = left_single_gait_cycle_cell_norm{j,gc_sort_inds(k),i};
            end
        end
    end

    left_sr = unique(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate);
    left_single_gait_cycle_time_vec = (0:max_gc_length-1).*(1/left_sr);
end


if isfield(lfp_cwt,'Right')
    temp = {right_single_gait_cycle_cell{1,:,1}};
    gc_lengths = cellfun(@(x)size(x,2),temp);
    max_gc_length = max(gc_lengths);
    [gc_sort_lengths,gc_sort_inds] = sort(gc_lengths,'descend');
    right_single_gait_cycle_mat = nan(sum(valid_gc),max_gc_length,length(freq_bands_names),4);      % gait cycle x time x freq band x channel

    for i = 1:length(lfp_cwt.Right.Chan_Names)
        for j = 1:length(freq_bands_names)
            for k = 1:length(gc_sort_inds)
                right_single_gait_cycle_mat(k,1:gc_sort_lengths(k),j,i) = right_single_gait_cycle_cell_norm{j,gc_sort_inds(k),i};
            end
        end
    end

    right_sr = unique(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate);
    right_single_gait_cycle_time_vec = (0:max_gc_length-1).*(1/right_sr);
end

% calc event time markers from LHS
RTO_relative_time = gait_events_sorted.('RTO')(valid_gc)-gait_events_sorted.LHS(valid_gc);
RHS_relative_time = gait_events_sorted.('RHS')(valid_gc)-gait_events_sorted.LHS(valid_gc);
LTO_relative_time = gait_events_sorted.('LTO')(valid_gc)-gait_events_sorted.LHS(valid_gc);

% colors for gait events
marker_colors = CBMap('GaitEvents');

% Check to see if there are any gait cycles to ignore and remove them
if ~isempty(sorted_gc_ignore)
    if exist('left_single_gait_cycle_mat','var')
        for i = 1:length(sorted_gc_ignore)
            left_single_gait_cycle_mat(sorted_gc_ignore(i),:,:,:) = [];
        end
    end

    if exist('right_single_gait_cycle_mat','var')
        for i = 1:length(sorted_gc_ignore)
            right_single_gait_cycle_mat(sorted_gc_ignore(i),:,:,:) = [];
        end
    end
end

% Plot the data
if exist('left_single_gait_cycle_mat','var')
    for i = 1:size(left_single_gait_cycle_mat,4)        % channel
        for j = 1:size(left_single_gait_cycle_mat,3)    % freq band
            figure;
            pcolor(left_single_gait_cycle_time_vec,sum(valid_gc):-1:1,left_single_gait_cycle_mat(:,:,j,i));
            title(['Left ',chan_names{i},' ',freq_bands_names{j}]);
            shading interp
            hold on;
            scatter(RTO_relative_time(gc_sort_inds),sum(valid_gc):-1:1,'MarkerEdgeColor',marker_colors.RTO,'MarkerFaceColor',marker_colors.RTO);
            scatter(RHS_relative_time(gc_sort_inds),sum(valid_gc):-1:1,'MarkerEdgeColor',marker_colors.RHS,'MarkerFaceColor',marker_colors.RHS);
            scatter(LTO_relative_time(gc_sort_inds),sum(valid_gc):-1:1,'MarkerEdgeColor',marker_colors.LTO,'MarkerFaceColor',marker_colors.LTO);
            yticks([1:5:sum(valid_gc),sum(valid_gc)]);
            ylabel('Gait Cycle');
            xlabel('Time (sec)');
        end
    end
end

if exist('right_single_gait_cycle_mat','var')
    for i = 1:size(right_single_gait_cycle_mat,4)       % channel
        for j = 1:size(right_single_gait_cycle_mat,3)   % freq band
            figure;
            pcolor(right_single_gait_cycle_time_vec,sum(valid_gc):-1:1,right_single_gait_cycle_mat(:,:,j,i));
            title(['Right ',chan_names{i},' ',freq_bands_names{j}]);
            shading interp
            hold on;
            scatter(RTO_relative_time(gc_sort_inds),sum(valid_gc):-1:1,'MarkerEdgeColor',marker_colors.RTO,'MarkerFaceColor',marker_colors.RTO);
            scatter(RHS_relative_time(gc_sort_inds),sum(valid_gc):-1:1,'MarkerEdgeColor',marker_colors.RHS,'MarkerFaceColor',marker_colors.RHS);
            scatter(LTO_relative_time(gc_sort_inds),sum(valid_gc):-1:1,'MarkerEdgeColor',marker_colors.LTO,'MarkerFaceColor',marker_colors.LTO);
            yticks([1:5:sum(valid_gc),sum(valid_gc)]);
            ylabel('Gait Cycle');
            xlabel('Time (sec)');
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
        normalization_sd(i,j) = std(mean(data_clean));
        normalization_var(i,j) = var(mean(data_clean));
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