function data_table = extractGaitCyclePercentCoherencePhase(files,subject_ID,stim_target,varargin)
%%
%   This function extracts the coherence between two contact pairs from the
%   same hemisphere for each frequency and gait cycle percentage during 
%   walking and in one of two conditions, rest prior to walking and the 
%   average coherence at a given frequency for each gait cycle. Stores the
%   data as a CSV file. Note, this CSV file will be GB in size.

%% EXAMPLE FUNCTION RUN
% % GPi
% RCS03 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS03/Gait/2020-10-05/Data/Aligned Data/RCS03_OG_after_task_OFF_STIM_w_Gait_Events_Julia.mat'; % Using "after task" but the task had a lot of issues. Unlikely they learned.
% RCS14 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS14/2021-04-12/Data/Aligned Data/RCS14_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
% gpi_patients = {RCS03,RCS14};
% 
% % STN
% RCS07 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/2019-12-17/Data/Aligned Data/RCS07_OG_before_task_OFF_STIM_w_Gait_Events.mat';
% RCS12 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS12/10 Visit/Data/Aligned Data/RCS12_OG_OFF_STIM_w_Gait_Events_Julia.mat';
% RCS15 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS15/2021-07-29/Data/Aligned Data/RCS15_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
% stn_patients = {RCS07,RCS12,RCS15};
% 
% files = [stn_patients,gpi_patients];
% subject_ID = {'RCS07','RCS12','RCS15','RCS03','RCS14'};
% 
% extractGaitCyclePercentCoherencePhase(files,subject_ID,{'STN';'STN';'STN';'GPi';'GPi'});

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

for i = 1:2:nargin-3
    switch varargin{i}
        case 'export_data'
            export_data = varargin{i+1};
        case 'remove_nan'
            remove_nan = varargin{i+1};
        case 'n_percent_bins'
            n_percent_bins = varargin{i+1};
        case 'cycle_start_event'
            cycle_start_event = varargin{i+1};
        case 'elements_to_compare'
            elements_to_compare = varargin{i+1};
        case 'pairs'
            pairs = varargin{i+1};
    end
end

if ~exist('export_data','var')
    export_data = true;
end

if ~exist('remove_nan','var')
    remove_nan = true;
end

if ~exist('n_percent_bins','var')
    n_percent_bins = 100;
end

if ~exist('cycle_start_event','var')
    cycle_start_event = 'LHS';
end

if ~exist('elements_to_compare','var')
    elements_to_compare = [];
end

if ~exist('pairs','var')
    pairs = [];
end

%% Data extraction
state = {};
coherence = [];
phase = [];
frequency = [];
gait_cycle_percent = [];
gait_cycle_num = [];
contact = {};
dbs_target = {};
subject = {};
for i = 1:length(files)
    % Load data
    printStatus(sprintf('Loading file: %s',files{i}));
    load(files{i});
    
    % Sort gait events
    gait_events_sorted = sortGaitEvents(aligned_data.gait_events,cycle_start_event);
    
    % Gait event search range
    gait_event_range = [];
    if contains(files{i},'RCS07')
        gait_event_range(1) = 1;
        gait_event_range(2) = height(gait_events_sorted);
    elseif contains(files{i},'RCS12')
        gait_event_range(1) = find(gait_events_sorted.LHS > 51,1,'first');
        gait_event_range(2) = find(gait_events_sorted.LHS < 234,1,'last');
    elseif contains(files{i},'RCS15')
        gait_event_range(1) = 1;
        gait_event_range(2) = height(gait_events_sorted);
    elseif contains(files{i},'RCS03')
        gait_event_range(1) = 1;
        gait_event_range(2) = height(gait_events_sorted);
    elseif contains(files{i},'RCS14')
        gait_event_range(1) = 1;
        gait_event_range(2) = height(gait_events_sorted);
    end
    
    % Determine what pairs to run
    if isfield(aligned_data,'left_LFP_table') && isempty(elements_to_compare) && isempty(pairs)
        elements_to_compare = [elements_to_compare,1:4];
    end
    
    if isfield(aligned_data,'right_LFP_table') && isempty(elements_to_compare) && isempty(pairs)
        elements_to_compare = [elements_to_compare,5:8];
    end
    
    if isempty(pairs)
        pairs = nchoosek(elements_to_compare,2);
    end
    contact_names = createContactVariableName(pairs);
    
    for j = 1:size(pairs,1)
        if pairsValid(pairs(j,:),aligned_data)
            % Determine the indexing
            switch pairs(j,1)
                case {1,5}
                    contact1 = 'key0';
                case {2,6}
                    contact1 = 'key1';
                case {3,7}
                    contact1 = 'key2';
                case {4,8}
                    contact1 = 'key3';
            end
            
            switch pairs(j,2)
                case {1,5}
                    contact2 = 'key0';
                case {2,6}
                    contact2 = 'key1';
                case {3,7}
                    contact2 = 'key2';
                case {4,8}
                    contact2 = 'key3';
            end
            
            if pairs(j,1) <= 4
                side1 = 'left_LFP_table';
            else
                side1 = 'right_LFP_table';
            end
            
            if pairs(j,2) <= 4
                side2 = 'left_LFP_table';
            else
                side2 = 'right_LFP_table';
            end
            
            % Determine index limits
            n_val1 = height(aligned_data.(side1));
            n_val2 = height(aligned_data.(side2));
            [n_val,min_ind] = min([n_val1,n_val2]);
            
            if min_ind == 1
                if contains(side1,'left')
                    time_vec = aligned_data.left_taxis;
                else
                    time_vec = aligned_data.right_taxis;
                end
            else
                if contains(side2,'left')
                    time_vec = aligned_data.left_taxis;
                else
                    time_vec = aligned_data.right_taxis;
                end
            end
            
            % Get sample rate
            sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
            
            % Calculate wavelet coherence and extract gait cycle values
            [x,wcs,y] = wcoherence(aligned_data.(side1).(contact1)(1:n_val),aligned_data.(side2).(contact2)(1:n_val),sr,'VoicesPerOctave',10);
            
            % Extract coherence
            printStatus(sprintf('Extracting data from pair: %s %s to %s %s',side1,contact1,side2,contact2));
            gait_cycle_mat = zeros(length(y),n_percent_bins,1);
            gait_cycle_phase_mat = zeros(length(y),n_percent_bins,1);
            coherence_temp = [];
            phase_temp = [];
            average_coherence_gait_cycle_mat = zeros(length(y),1);
            average_coherence_temp = [];
            
            count = 1;
            for k = gait_event_range(1):gait_event_range(2)-1
                if ~isnan(gait_events_sorted.(cycle_start_event)(k)) && ~isnan(gait_events_sorted.(cycle_start_event)(k+1)) && (diff(gait_events_sorted.(cycle_start_event)(k:k+1)) < 2)
                    [~,start_ind] = min(abs(time_vec-gait_events_sorted.(cycle_start_event)(k)));
                    [~,end_ind] = min(abs(time_vec-gait_events_sorted.(cycle_start_event)(k+1)));
                    data_snip = x(:,start_ind:end_ind);
                    data_snip_phase = wcs(:,start_ind:end_ind);
                    
                    if sum(isinf(data_snip),'all') == 0
                        average_coherence_gait_cycle_mat(:,1,count) = mean(data_snip,2);
                        percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                        for m = 1:length(percent_inds)-1
                            if m == 1
                                gait_cycle_mat(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                                gait_cycle_phase_mat(:,m,count) = mean(angle(data_snip_phase(:,percent_inds(m):percent_inds(m+1))),2);
                            else
                                gait_cycle_mat(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                                gait_cycle_phase_mat(:,m,count) = mean(angle(data_snip_phase(:,percent_inds(m)+1:percent_inds(m+1))),2);
                            end
                        end
                        count = count + 1;
                    end
                end
            end
            
            % Re-formulate coherence matrices into a large 1 column vector
            % then populate other column variables with the correct number of
            % values
            n_gait_cycles = size(gait_cycle_mat,3);
            for n = 1:n_gait_cycles
                for o = 1:n_percent_bins
                    coherence_temp = [coherence_temp;gait_cycle_mat(:,o,n)];
                    phase_temp = [phase_temp;gait_cycle_phase_mat(:,o,n)];
                end
            end
            coherence = [coherence;coherence_temp];
            phase = [phase;phase_temp];
            frequency = [frequency;repmat(repmat(round(y,4),n_percent_bins,1),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(y)*n_percent_bins)'];
            gait_cycle_percent = [gait_cycle_percent;repmat(repelem(1:100,length(y))',n_gait_cycles,1)];
            state = [state;repmat({'Walking'},length(y)*n_percent_bins*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(y)*n_percent_bins*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(y)*n_percent_bins*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(y)*n_percent_bins*n_gait_cycles,1)];
            
            
            n_gait_cycles = size(average_coherence_gait_cycle_mat,3);
            for n = 1:n_gait_cycles
                average_coherence_temp = [average_coherence_temp;average_coherence_gait_cycle_mat(:,1,n)];
            end
            coherence = [coherence;average_coherence_temp];
            phase = [phase;nan(size(average_coherence_temp))];
            frequency = [frequency;repmat(round(y,4),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(y))'];
            gait_cycle_percent = [gait_cycle_percent;nan(length(y)*n_gait_cycles,1)];
            state = [state;repmat({'Average'},length(y)*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(y)*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(y)*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(y)*n_gait_cycles,1)];
        end
    end
    printStatus(' ');
end
printStatus('Data extraction complete.');
printStatus(' ');
printStatus();

%% Clean up data
if remove_nan
    printStatus('Removing nan from extracted data.');
    nan_inds = ismissing(coherence);
    remove_inds = find(nan_inds(:,1));
    coherence(remove_inds,:) = [];
    phase(remove_inds,:) = [];
    subject(remove_inds) = [];
    dbs_target(remove_inds) = [];
    contact(remove_inds) = [];
    state(remove_inds) = [];
    gait_cycle_num(remove_inds) = [];
    gait_cycle_percent(remove_inds) = [];
    frequency(remove_inds) = [];
    
    % repeat for phase
    nan_inds = ismissing(phase);
    remove_inds = find(nan_inds(:,1));
    coherence(remove_inds,:) = [];
    phase(remove_inds,:) = [];
    subject(remove_inds) = [];
    dbs_target(remove_inds) = [];
    contact(remove_inds) = [];
    state(remove_inds) = [];
    gait_cycle_num(remove_inds) = [];
    gait_cycle_percent(remove_inds) = [];
    frequency(remove_inds) = [];

    printStatus('Data cleaning complete.');
    printStatus(' ');
    printStatus();
end

%% Create table
printStatus('Creating data table...');
data_table = table(subject,dbs_target,contact,state,gait_cycle_num,gait_cycle_percent,frequency,coherence,phase,'VariableNames',{'Subject','Target','Contact','State','GC','GCP','Frequency','Power','Phase'});
printStatus('Table creation complete.');
printStatus(' ');
printStatus();

%% Save table
if export_data
    printStatus('Saving data table');
    [save_file,save_path] = uiputfile({'*.csv'});
    save_name = fullfile(save_path,save_file);
    
    writetable(data_table,save_name);  
    save(strrep(save_name,'.csv','.mat'),'data_table');
    printStatus(sprintf('Gait Cycle Coherence and Phase table saved to: %s',save_path));
end

printStatus('Save successful.');
printStatus(' ');
printStatus();
end

%% Helper functions
function printStatus(output_string)
if ~exist('output_string','var')
    fprintf('%s\n',repmat('%',1,50));
else
    fprintf('%s\n',output_string);
end
end


function contact_output_name = createContactVariableName(contact_pairs)
% create cell of contact names for coherence signal pairs
contact_output_name = cell(1,size(contact_pairs,1));

for i = 1:size(contact_pairs,1)
    % Determine the indexing
    switch contact_pairs(i,1)
        case {1,5}
            contact1 = 'key0';
        case {2,6}
            contact1 = 'key1';
        case {3,7}
            contact1 = 'key2';
        case {4,8}
            contact1 = 'key3';
    end
    
    switch contact_pairs(i,2)
        case {1,5}
            contact2 = 'key0';
        case {2,6}
            contact2 = 'key1';
        case {3,7}
            contact2 = 'key2';
        case {4,8}
            contact2 = 'key3';
    end
    
    if contact_pairs(i,1) <= 4
        side1 = 'L';
    else
        side1 = 'R';
    end
    
    if contact_pairs(i,2) <= 4
        side2 = 'L';
    else
        side2 = 'R';
    end
    
    contact_output_name{i} = [side1,contact1,'_',side2,contact2];
end
end


function valid = pairsValid(pairs,reference)
valid = false;

if (pairs(1) <= 4 && isfield(reference,'left_LFP_table')) && (pairs(2) <= 4 && isfield(reference,'left_LFP_table'))
    valid = true;
elseif (pairs(1) <= 4 && isfield(reference,'left_LFP_table')) && (pairs(2) >= 5 && isfield(reference,'right_LFP_table'))
    valid = true;
elseif (pairs(1) >= 5 && isfield(reference,'right_LFP_table')) && (pairs(2) >= 5 && isfield(reference,'right_LFP_table'))
    valid = true;
elseif (pairs(1) >= 5 && isfield(reference,'right_LFP_table')) && (pairs(2) <= 4 && isfield(reference,'left_LFP_table'))
    valid = true;
end
end