function [data_power_table,data_phase_table] = extractDataFeaturesForML(files,subject_ID,stim_target,visit_number,varargin)

% % GPi
% RCS03 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS03/Gait/2020-10-05/Data/Aligned Data/RCS03_OG_after_task_OFF_STIM_w_Gait_Events_Julia.mat'}; % Using "after task" but the task had a lot of issues. Unlikely they learned.
% RCS14 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS14/2021-04-12/Data/Aligned Data/RCS14_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat'};
% gRCS01 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_01/v3_2021-12-20_pre-programming/Data/Aligned Data/gait_RCS_01_OG_ON_Meds_Trial1_w_Gait_Events_Ken.mat'};
% gRCS02 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_02/v4_2022-01-14_pre-programming_second/Data/Aligned Data/gait_RCS_02_OG_OFF_STIM_OFF_MEDS_w_Gait_Events_Julia.mat';...
%     '/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_02/v4_2022-01-14_pre-programming_second/Data/Aligned Data/gait_RCS_02_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat'};
% gpi_patients = [RCS03;RCS14;gRCS01;gRCS02];
% 
% % STN
% RCS07 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/2019-12-17/Data/Aligned Data/RCS07_OG_before_task_OFF_STIM_w_Gait_Events.mat'};
% RCS12 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS12/10 Visit/Data/Aligned Data/RCS12_OG_OFF_STIM_w_Gait_Events_Julia.mat'};
% RCS15 = {'/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS15/2021-07-29/Data/Aligned Data/RCS15_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat'};
% stn_patients = [RCS07;RCS12;RCS15];
% 
% patient_files = [stn_patients;gpi_patients];
% patient_names = [{'RCS07'};{'RCS12'};{'RCS15'};{'RCS03'};{'RCS14'};repelem({'gRCS01'},length(gRCS01),1);repelem({'gRCS02'},length(gRCS02),1)];
% dbs_target = [repelem({'STN'},length(stn_patients),1);repelem({'GPi'},length(gpi_patients),1)];
% visit_number = [1;1;1;1;1;...       % Phil's patients
%     3;...                           % gRCS01
%     3;3];                           % gRCS02
% extractDataFeaturesForML(patient_files,patient_names,dbs_target,visit_number,'STFT_overlap',0.9);
% extractDataFeaturesForML(patient_files,patient_names,dbs_target,visit_number,'power_type','peri','adj_time',0.05,'STFT_overlap',0.9);
% extractDataFeaturesForML(patient_files,patient_names,dbs_target,visit_number,'power_type','peri','adj_time',0.1,'STFT_overlap',0.9);
% extractDataFeaturesForML(patient_files,patient_names,dbs_target,visit_number,'power_type','pre','STFT_overlap',0.9);

%% Initiation
if ~exist('files','var') || isempty(files)
    error('No files to parse');
end

if ~exist('subject_ID','var') || isempty(subject_ID)
    error('Files must have subject identifiers');
end

if ~exist('stim_target','var') || isempty(stim_target)
    error('Files must have stimulation target indentifiers');
end

if ~exist('visit_number','var') || isempty(visit_number)
    error('Need to know visit numbers to distinguish data.');
end

for i = 1:2:nargin-4
    switch varargin{i}
        case 'previous_table'
            old_table = varargin{i+1};
        case 'freq_search_limits'
            freq_limit = varargin{i+1};
        case 'STFT_overlap'
            STFT_overlap = varargin{i+1};
        case 'export_data'
            export_data = varargin{i+1};
        case 'remove_nan'
            remove_nan = varargin{i+1};
        case 'power_type'
            power_type = varargin{i+1};
        case 'adj_time'
            adj_time = varargin{i+1};  % In seconds
    end
end

if ~exist('freq_limit','var')
    freq_limit = [2.5,50];
end

if ~exist('STFT_overlap','var')
    STFT_overlap = 0.9;
end

if ~exist('export_data','var')
    export_data = true;
end

if ~exist('remove_nan','var')
    remove_nan = true;
end

if ~exist('power_type','var')
    power_type = 'instantaneous';   % instantaneous|peri|pre
end

if (strcmp(power_type,'peri') || strcmp(power_type,'pre')) && ~exist('adj_time','var')
    warning(sprintf('For %s power around gait event, adj_time is not set. Will use patient inter-heel-strike interval',power_type));
    calc_adj_time = true;
    adj_time = [];
else
    calc_adj_time = false;
end

if ~strcmp(power_type,'instantaneous')
    if length(adj_time) == 1
        adj_time = repelem(adj_time,length(files));
    elseif ~isempty(adj_time) && length(adj_time) ~= length(files)
        error('Number of adjustment time does not match number of files. To use one value for all files, input only one value.');
    end
end

%% output variables
data_power_table = [];
data_phase_table = [];

%% inbetween loop variables
data_power_mat = [];
data_phase_mat = [];
ids = {};
targets = {};
side = {};
contact = {};
event = {};
gait_cycle = [];
med_state = {};
stim_state = {};
visit = [];

if strcmp(power_type,'peri') || strcmp(power_type,'pre')
    time_range = [];
end

freq_vals = [];

%% Extract data
printStatus();

for i = 1:length(files)
    % Load data
    printStatus(sprintf('Loading file: %s',files{i}));
    load(files{i});
    signal_analysis_data = calcRCS_STFT(aligned_data,[],1,STFT_overlap,[]); % set STFT_overlap to 0.996 for old version
    
    % gait event search range
    gait_event_range = getGaitEventRange(subject_ID{i},aligned_data,files{i});
    
    
    % If adj_time is not set and power type is peri or pre event, use
    % interstep or 25% (15%?) of gait cycle
    if (strcmp(power_type,'peri') || strcmp(power_type,'pre')) && calc_adj_time
%         adj_time(i) = calcAdjTime(aligned_data.gait_events,'inter-heel-strike');
        adj_time(i) = calcAdjTime(aligned_data.gait_events,'percent_gait_cycle');
    end
    
    % Check if frequency bins are different from previous table
    if exist('old_table','var')
        if ~isSameBinsAsPrevious(signal_analysis_data,old_table)
            response = questdlg({'New frequency features in the old table do not match.';'Create new table with new fequency features?';
                '(will not overwrite old table)'}, ...
                '','Yes','No','No');  
            if strcmp(response,'No')
                error('Aborting parsing');
            end
        end
    end    
    
    % Check if Left field exist
    if isfield(signal_analysis_data,'Left')
        printStatus('Extracting data from left contacts.');
        [freq_bin_inds,freq_vals_new] = genFreqBinPairs(signal_analysis_data.Left.Freq_Values{1},freq_limit);
        if ~isempty(freq_vals)
            if sum(freq_vals==freq_vals_new) ~= length(freq_vals)
                error('Frequency vals are different.');
            end
        else
            freq_vals = freq_vals_new;
        end
        
        n_channels = length(signal_analysis_data.Left.Chan_Names);
        n_event_names = length(aligned_data.gait_events.Properties.VariableNames);
        n_events = height(aligned_data.gait_events);
        n_freq_bins = size(freq_bin_inds,1);
        temp_power_data = nan(n_channels*n_event_names*n_events,n_freq_bins);
        temp_phase_data = nan(n_channels*n_event_names*n_events,n_freq_bins);
        for j = 1:n_channels
            printStatus(sprintf('Starting extraction Left %s',signal_analysis_data.Left.Chan_Names{j}));
            for k = 1:n_event_names
                for m = 1:n_events
                    if aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m) >= gait_event_range(1) && ...
                            aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m) <= gait_event_range(2)
                        temp_power = nan(1,n_freq_bins);
                        temp_phase = nan(1,n_freq_bins);
                        event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m);
                        if ~isnan(event_time)
                            if strcmp(power_type,'instantaneous')
                                [~,event_ind] = min(abs(signal_analysis_data.Left.Time{j}-event_time));
                            elseif strcmp(power_type,'peri')
                                [~,start_ind] = min(abs(signal_analysis_data.Left.Time{j}-(event_time-adj_time(i))));
                                [~,end_ind] = min(abs(signal_analysis_data.Left.Time{j}-(event_time+adj_time(i))));
                                event_ind = start_ind:end_ind;
                            elseif strcmp(power_type,'pre')
                                [~,start_ind] = min(abs(signal_analysis_data.Left.Time{j}-(event_time-adj_time(i))));
                                [~,end_ind] = min(abs(signal_analysis_data.Left.Time{j}-(event_time)));
                                event_ind = start_ind:end_ind;
                            end
                            
                            if isfield(signal_analysis_data.Left,'PSD')
                                power_snip = 10*log10(abs(signal_analysis_data.Left.Values{j}(:,event_ind)));
                                phase_snip = angle(signal_analysis_data.Left.Values{j}(:,event_ind));
                            else
                                power_snip = abs(signal_analysis_data.Left.Values{j}(:,event_ind));
                                phase_snip = angle(signal_analysis_data.Left.Values{j}(:,event_ind));
                            end
                            
                            if sum(isinf(power_snip),'all') == 0
                                parfor n = 1:n_freq_bins
                                    mean_power = mean(power_snip(freq_bin_inds(n,1):freq_bin_inds(n,2),:),'all');
                                    mean_phase = mean(phase_snip(freq_bin_inds(n,1):freq_bin_inds(n,2),:),'all');
                                    
                                    if ~isinf(mean_power)
                                        temp_power(1,n) = mean_power;
                                        temp_phase(1,n) = mean_phase;
                                    end
                                end
                            end
                            
                            temp_power_data((j-1)*(n_event_names*n_events)+(k-1)*n_events+m,:) = temp_power;
                            temp_phase_data((j-1)*(n_event_names*n_events)+(k-1)*n_events+m,:) = temp_phase;
                        end
                    end
                end
            end
            printStatus(sprintf('Completed extraction of Left %s',signal_analysis_data.Left.Chan_Names{j}));
        end
        
%         temp_id = repelem(subject_ID(i),size(temp_power_data,1),1);
%         temp_stim_targets = repelem(stim_target(i),size(temp_power_data,1),1);
%         temp_side = repelem({'L'},size(temp_power_data,1),1);
%         temp_contact = transpose(repelem(signal_analysis_data.Left.Chan_Names,n_events*n_event_names));
%         temp_event = repmat(transpose(repelem(aligned_data.gait_events.Properties.VariableNames,n_events)),n_channels,1);
        
        data_power_mat = [data_power_mat;temp_power_data];
        data_phase_mat = [data_phase_mat;temp_phase_data];
        ids = [ids;repelem(subject_ID(i),size(temp_power_data,1),1)];
        targets = [targets;repelem(stim_target(i),size(temp_power_data,1),1)];
        side = [side;repelem({'L'},size(temp_power_data,1),1)];
        contact = [contact;transpose(repelem(signal_analysis_data.Left.Chan_Names,n_events*n_event_names))];
        event = [event;repmat(transpose(repelem(aligned_data.gait_events.Properties.VariableNames,n_events)),n_channels,1)];
        gait_cycle = [gait_cycle;repmat(transpose(1:n_events),length(aligned_data.gait_events.Properties.VariableNames)*n_channels,1)];
        visit = [visit;repelem(visit_number(i),size(temp_power_data,1),1)];
        
        if isfield(aligned_data,'stim_condition')
            stim_state = [stim_state;repmat({aligned_data.stim_condition},size(temp_power_data,1),1)];
        else % assume off stim
            stim_state = [stim_state;repmat({'OFF'},size(temp_power_data,1),1)];
        end
        
        if isfield(aligned_data,'med_condition')
            med_state = [med_state;repmat({aligned_data.med_condition},size(temp_power_data,1),1)];
        else % assume on med
            med_state = [med_state;repmat({'ON'},size(temp_power_data,1),1)];
        end
        
        if strcmp(power_type,'peri') || strcmp(power_type,'pre')
            time_range = [time_range;repelem(adj_time(i),size(temp_power_data,1),1)];
        end
        
        printStatus('');
    end
    
    
    % Check if Right field exist
    if isfield(signal_analysis_data,'Right')
        printStatus('Extracting data from right contacts.');
        [freq_bin_inds,freq_vals_new] = genFreqBinPairs(signal_analysis_data.Right.Freq_Values{1},freq_limit);
        if ~isempty(freq_vals)
            if sum(freq_vals==freq_vals_new) ~= length(freq_vals)
                error('Frequency vals are different.');
            end
        else
            freq_vals = freq_vals_new;
        end
        
        n_channels = length(signal_analysis_data.Right.Chan_Names);
        n_event_names = length(aligned_data.gait_events.Properties.VariableNames);
        n_events = height(aligned_data.gait_events);
        n_freq_bins = size(freq_bin_inds,1);
        temp_power_data = nan(n_channels*n_event_names*n_events,n_freq_bins);
        temp_phase_data = nan(n_channels*n_event_names*n_events,n_freq_bins);
        for j = 1:n_channels
            printStatus(sprintf('Starting extraction Right %s',signal_analysis_data.Right.Chan_Names{j}));
            for k = 1:n_event_names
                for m = 1:n_events
                    temp_power = nan(1,n_freq_bins);
                    temp_phase = nan(1,n_freq_bins);
                    event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m);
                    if ~isnan(event_time)
                        if strcmp(power_type,'instantaneous')
                            [~,event_ind] = min(abs(signal_analysis_data.Right.Time{j}-event_time));
                        elseif strcmp(power_type,'peri')
                            [~,start_ind] = min(abs(signal_analysis_data.Right.Time{j}-(event_time-adj_time(i))));
                            [~,end_ind] = min(abs(signal_analysis_data.Right.Time{j}-(event_time+adj_time(i))));
                            event_ind = start_ind:end_ind;
                        elseif strcmp(power_type,'pre')
                            [~,start_ind] = min(abs(signal_analysis_data.Right.Time{j}-(event_time-adj_time(i))));
                            [~,end_ind] = min(abs(signal_analysis_data.Right.Time{j}-(event_time)));
                            event_ind = start_ind:end_ind;
                        end
                        
                        if isfield(signal_analysis_data.Right,'PSD')
                            power_snip = 10*log10(abs(signal_analysis_data.Right.Values{j}(:,event_ind)));
                            phase_snip = angle(signal_analysis_data.Right.Values{j}(:,event_ind));
                        else
                            power_snip = abs(signal_analysis_data.Right.Values{j}(:,event_ind));
                            phase_snip = angle(signal_analysis_data.Right.Values{j}(:,event_ind));
                        end
                        
                        if sum(isinf(power_snip),'all') == 0
                            parfor n = 1:n_freq_bins
                                mean_power = mean(power_snip(freq_bin_inds(n,1):freq_bin_inds(n,2),:),'all');
                                mean_phase = mean(phase_snip(freq_bin_inds(n,1):freq_bin_inds(n,2),:),'all');
                                
                                if ~isinf(mean_power)
                                    temp_power(1,n) = mean_power;
                                    temp_phase(1,n) = mean_phase;
                                end
                            end
                        end
                        temp_power_data((j-1)*(n_event_names*n_events)+(k-1)*n_events+m,:) = temp_power;
                        temp_phase_data((j-1)*(n_event_names*n_events)+(k-1)*n_events+m,:) = temp_phase;
                    end
                end
            end
            printStatus(sprintf('Completed extraction of Right %s',signal_analysis_data.Right.Chan_Names{j}));
        end
        
%         temp_id = repelem(subject_ID(i),size(temp_power_data,1),1);
%         temp_stim_targets = repelem(stim_target(i),size(temp_power_data,1),1);
%         temp_side = repelem({'L'},size(temp_power_data,1),1);
%         temp_contact = transpose(repelem(signal_analysis_data.Left.Chan_Names,n_events*n_event_names));
%         temp_event = repmat(transpose(repelem(aligned_data.gait_events.Properties.VariableNames,n_events)),n_channels,1);
        
        data_power_mat = [data_power_mat;temp_power_data];
        data_phase_mat = [data_phase_mat;temp_phase_data];
        ids = [ids;repelem(subject_ID(i),size(temp_power_data,1),1)];
        targets = [targets;repelem(stim_target(i),size(temp_power_data,1),1)];
        side = [side;repelem({'R'},size(temp_power_data,1),1)];
        contact = [contact;transpose(repelem(signal_analysis_data.Right.Chan_Names,n_events*n_event_names))];
        event = [event;repmat(transpose(repelem(aligned_data.gait_events.Properties.VariableNames,n_events)),n_channels,1)];
        gait_cycle = [gait_cycle;repmat(transpose(1:n_events),length(aligned_data.gait_events.Properties.VariableNames)*n_channels,1)];
        visit = [visit;repelem(visit_number(i),size(temp_power_data,1),1)];
        
        if isfield(aligned_data,'stim_condition')
            stim_state = [stim_state;repmat({aligned_data.stim_condition},size(temp_power_data,1),1)];
        else % assume off stim
            stim_state = [stim_state;repmat({'OFF'},size(temp_power_data,1),1)];
        end
        
        if isfield(aligned_data,'med_condition')
            med_state = [med_state;repmat({aligned_data.med_condition},size(temp_power_data,1),1)];
        else % assume on med
            med_state = [med_state;repmat({'ON'},size(temp_power_data,1),1)];
        end
        
        if strcmp(power_type,'peri') || strcmp(power_type,'pre')
            time_range = [time_range;repelem(adj_time(i),size(temp_power_data,1),1)];
        end
        
        printStatus('');
    end
end

%% Clean up data
if remove_nan
    printStatus('Removing nan from extracted data.');
    nan_inds = ismissing(data_power_mat);
    remove_inds = find(nan_inds(:,1));
    data_power_mat(remove_inds,:) = [];
    data_phase_mat(remove_inds,:) = [];
    ids(remove_inds) = [];
    targets(remove_inds) = [];
    side(remove_inds) = [];
    contact(remove_inds) = [];
    event(remove_inds) = [];
    gait_cycle(remove_inds) = [];
    visit(remove_inds) = [];
    med_state(remove_inds) = [];
    stim_state(remove_inds) = [];
    
    if strcmp(power_type,'peri') || strcmp(power_type,'pre')
        time_range(remove_inds) = [];
    end

    printStatus('');
end

%% Create table
printStatus('Creating power data table');
col_names = createTableVariableNames(freq_bin_inds,freq_vals);
data_power_table = array2table(data_power_mat,'VariableNames',col_names);

if strcmp(power_type,'peri') || strcmp(power_type,'pre')
    data_power_table = addvars(data_power_table,ids,visit,targets,med_state,stim_state,side,contact,event,gait_cycle,time_range,'Before',1,'NewVariableNames',{'Subject_ID','Visit','Stim_Target','Med_State','Stim_State','Side','Contact','Gait_Event','Gait_Cycle','Time_Range'});
else
    data_power_table = addvars(data_power_table,ids,visit,targets,med_state,stim_state,side,contact,event,gait_cycle,'Before',1,'NewVariableNames',{'Subject_ID','Visit','Stim_Target','Med_State','Stim_State','Side','Contact','Gait_Event','Gait_Cycle'});
end

printStatus('Creating phase data table');
data_phase_table = array2table(data_phase_mat,'VariableNames',col_names);

if strcmp(power_type,'peri') || strcmp(power_type,'pre')
    data_phase_table = addvars(data_phase_table,ids,visit,targets,med_state,stim_state,side,contact,event,gait_cycle,time_range,'Before',1,'NewVariableNames',{'Subject_ID','Visit','Stim_Target','Med_State','Stim_State','Side','Contact','Gait_Event','Gait_Cycle','Time_Range'});
else
    data_phase_table = addvars(data_phase_table,ids,visit,targets,med_state,stim_state,side,contact,event,gait_cycle,'Before',1,'NewVariableNames',{'Subject_ID','Visit','Stim_Target','Med_State','Stim_State','Side','Contact','Gait_Event','Gait_Cycle'});
end

%% Save data
if export_data
    [save_file,save_path] = uiputfile({'*.csv'});
    if strcmp(power_type,'instantaneous')
        save_name = fullfile(save_path,[datestr(datetime('today','Format','yyyy-MM-dd'),'yyyy-mm-dd'),'_',strrep(save_file,'.csv',['_',power_type,'_power.csv'])]);
    elseif strcmp(power_type,'peri') || strcmp(power_type,'pre')
        if length(unique(adj_time)) == 1
            save_name = fullfile(save_path,[datestr(datetime('today','Format','yyyy-MM-dd'),'yyyy-mm-dd'),'_',strrep(save_file,'.csv',['_',power_type,'_',strrep(num2str(adj_time(1)),'.','__'),'_power.csv'])]);
        else
            save_name = fullfile(save_path,[datestr(datetime('today','Format','yyyy-MM-dd'),'yyyy-mm-dd'),'_',strrep(save_file,'.csv',['_',power_type,'_power.csv'])]);
        end
    end
    
    writetable(data_power_table,save_name);
%     save(strrep(save_name,'.csv','.mat'),'data_power_table'); % Files are too big and can't save as a .mat file
    printStatus(sprintf('Power table saved to: %s',save_path));
    
    save_name = strrep(save_name,'power','phase');
    writetable(data_phase_table,save_name);
%     save(strrep(save_name,'.csv','.mat'),'data_phase_table'); % Files are too big and can't save as a .mat file
    printStatus(sprintf('Phase table saved to: %s',save_path));
end

printStatus('Extraction complete.');
printStatus();
end

function printStatus(output_string)
if ~exist('output_string','var')
    fprintf('%s\n',repmat('%',1,50));
else
    fprintf('%s\n',output_string);
end
end

function val = isSameBinsAsPrevious(signal_analysis_data,old_table)
val = true;
end

function [freq_bin_inds,freqs] = genFreqBinPairs(freq_vals,freq_lim)
if freq_vals(2) < freq_vals(1)
    freq_vals_sorted = flipud(freq_vals);
else
    freq_vals_sorted = freq_vals;
end

ind_start = find(freq_vals_sorted >= freq_lim(1),1,'first');
ind_end = find(freq_vals_sorted <= freq_lim(2),1,'last');
n = ind_end - ind_start + 1;
skips = 1:n-1;
tot_pairs = (n^2+n)/2;
freqs = freq_vals_sorted(ind_start:ind_end);
freq_bin_inds = nan(tot_pairs,2);

% Crazy indexing. Look up triangle numbers to recalculate these equations
f = @(x) n*(x+1)-(x^2+x)/2;
freq_bin_inds(1:n,:) = repmat([ind_start:ind_end]',1,2);
for i = 1:length(skips)
    freq_bin_inds(f(i-1)+1:f(i),1) = transpose(ind_start:ind_end-skips(i));
    freq_bin_inds(f(i-1)+1:f(i),2) = freq_bin_inds(f(i-1)+1:f(i),1) + skips(i);
end

if freq_vals(2) < freq_vals(1)
    freq_bin_inds = abs(freq_bin_inds-length(freq_vals));
    freq_bin_inds = fliplr(freq_bin_inds);
end
end

function table_variable_names = createTableVariableNames(pair_inds,freq_vals)
table_variable_names = cell(1,size(pair_inds,1));

if min(pair_inds(:,1)) > 1 && length(min(pair_inds(:,1)):max(pair_inds(:,1))) == length(freq_vals)
    pair_inds = pair_inds - min(pair_inds(:,1))+1;
end

for i = 1:size(pair_inds,1)
    var_name_start = strrep(num2str(freq_vals(pair_inds(i,1))),'.','_');
    var_name_end = strrep(num2str(freq_vals(pair_inds(i,2))),'.','_');
    
    table_variable_names{i} = ['f',var_name_start,'__',var_name_end];
end
end

function adjustment_time = calcAdjTime(gait_event_table,type)
sorted_gait_event_table = sort_gait_events(gait_event_table,'LHS');

switch type
    case 'inter-heel-strike'
        data_snip = [sorted_gait_event_table.LHS,sorted_gait_event_table.RHS];
        
        found_match = 0;
        count = 1;
        while ~found_match
            if ~isnan(data_snip(count,1)) && ~isnan(data_snip(count,2))
                found_match = 1;
            end
            count = count + 1;
        end
        adjustment_time = mean(diff(data_snip,[],2),'omitnan');
        
    case 'percent_gait_cycle'
        gait_cycle_percentages = calcGaitCyclePercentage(sorted_gait_event_table);
        gait_cycle_times = diff(gait_event_table.LHS);
        inds = find(gait_cycle_times<2);
        
        adjustment_percent = min([1-gait_cycle_percentages.LTO,...
            gait_cycle_percentages.RHS-gait_cycle_percentages.RTO,...
            gait_cycle_percentages.RTO-gait_cycle_percentages.LHS,...
            gait_cycle_percentages.LTO-gait_cycle_percentages.RHS]);
        adjustment_time = mean(gait_cycle_times(inds))*adjustment_percent;
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

function gait_cycle_percentages = calcGaitCyclePercentage(gait_event_table)
gait_cycle_percentages.LHS = [];
gait_cycle_percentages.RTO = [];
gait_cycle_percentages.RHS = [];
gait_cycle_percentages.LTO = [];

gait_cycle_mat = [];
count = 1;
for i = 1:height(gait_event_table)-1
    if ~isnan(gait_event_table.LHS(i)) && ~isnan(gait_event_table.LHS(i+1))
        gait_cycle_time = gait_event_table.LHS(i+1)-gait_event_table.LHS(i);
        gait_cycle_percentages.LHS(count) = (gait_event_table.LHS(i)-gait_event_table.LHS(i))/gait_cycle_time;
        gait_cycle_percentages.RTO(count) = (gait_event_table.RTO(i)-gait_event_table.LHS(i))/gait_cycle_time;
        gait_cycle_percentages.RHS(count) = (gait_event_table.RHS(i)-gait_event_table.LHS(i))/gait_cycle_time;
        gait_cycle_percentages.LTO(count) = (gait_event_table.LTO(i)-gait_event_table.LHS(i))/gait_cycle_time;
        count = count + 1;
    end
end

gait_cycle_percentages.LHS = mean(gait_cycle_percentages.LHS,'omitnan');
gait_cycle_percentages.RTO = mean(gait_cycle_percentages.RTO,'omitnan');
gait_cycle_percentages.RHS = mean(gait_cycle_percentages.RHS,'omitnan');
gait_cycle_percentages.LTO = mean(gait_cycle_percentages.LTO,'omitnan');

end

function gait_event_range = getGaitEventRange(subjectID,aligned_data,file)

gait_event_range = [];
if strcmp(subjectID,'RCS07')
    gait_event_range(1) = 0;
    gait_event_range(2) = max(aligned_data.gait_events{end,:});
elseif strcmp(subjectID,'RCS12')
    gait_event_range(1) = 51;
    gait_event_range(2) = 234;
elseif strcmp(subjectID,'RCS15')
    gait_event_range(1) = 0;
    gait_event_range(2) = max(aligned_data.gait_events{end,:});
elseif strcmp(subjectID,'RCS03')
    gait_event_range(1) = 0;
    gait_event_range(2) = max(aligned_data.gait_events{end,:});
elseif strcmp(subjectID,'RCS14')
    gait_event_range(1) = 0;
    gait_event_range(2) = max(aligned_data.gait_events{end,:});
elseif strcmp(subjectID,'gRCS01')
    gait_event_range(1) = 0;
    gait_event_range(2) = max(aligned_data.gait_events{end,:});
elseif strcmp(subjectID,'gRCS02')
    gait_event_range(1) = 0;
    gait_event_range(2) = max(aligned_data.gait_events{end,:});
end

end