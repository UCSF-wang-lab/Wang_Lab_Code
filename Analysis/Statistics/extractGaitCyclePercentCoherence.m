function data_table = extractGaitCyclePercentCoherence(files,subject_ID,stim_target,varargin)
%% extractGaitCyclePercentCoherence
%   This function extracts the coherence between two contact pairs from the
%   same hemisphere for each frequency and gait cycle percentage during 
%   walking and in one of two conditions, rest prior to walking and the 
%   average coherence at a given frequency for each gait cycle. Stores the
%   data as a CSV file. Note, this CSV file will be GB in size.

%% EXAMPLE FUNCTION RUN
% RCS03 = '/Users/klouie/Documents/Backup Data/RCS03_OG_after_task_OFF_STIM_w_Gait_Events_Julia.mat';
% RCS07 = '/Users/klouie/Documents/Backup Data/RCS07_OG_before_task_OFF_STIM_w_Gait_Events.mat';
% RCS12 = '/Users/klouie/Documents/Backup Data/RCS12_OG_OFF_STIM_w_Gait_Events_Julia.mat';
% RCS14 = '/Users/klouie/Documents/Backup Data/RCS14_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
% RCS15 = '/Users/klouie/Documents/Backup Data/RCS15_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
% gRCS01 = '/Users/klouie/Documents/Backup Data/gait_RCS_01_OG_ON_Meds_Trial1_w_Gait_Events_Ken.mat';
% gRCS02 = '/Users/klouie/Documents/Backup Data/gait_RCS_02_OG_OFF_STIM_ON_MEDS_w_Gait_Events_Julia.mat';
% gRCS03 = '/Users/klouie/Documents/Backup Data/gait_RCS_03_OG_OFF_Stim_ON_Meds_Trial1_w_Gait_Events_Julia.mat';
% allFiles = {RCS07;RCS12;RCS15;RCS03;RCS14;gRCS01;gRCS02;gRCS03}
% subject_ID = {'RCS07','RCS12','RCS15','RCS03','RCS14','gRCS01','gRCS02','gRCS03'}
% stim_target = {'STN';'STN';'STN';'GPi';'GPi';'GPi';'GPi';'GPi'}
% geRangeTable = [1,inf;51,234;1,inf;1,inf;1,inf;1,inf;1,inf;1,inf]
% extractGaitCyclePercentCoherence(allFiles,subject_ID,stim_target,'geRangeTable',geRangeTable);

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
        case 'gcStartEvent'
            gcStartEvent = varargin{i+1};
        case 'geRangeTable'
            geRangeTable = varargin{i+1};
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

if ~exist('gcStartEvent','var')
    gcStartEvent = 'LHS';
end

if ~exist('geRangeTable','var') || isempty(geRangeTable)
    geRangeTable = [];
end

%% Data extraction
state = {};
coherence = [];
frequency = [];
gait_cycle_percent = [];
gait_cycle_num = [];
contact = {};
side = {};
dbs_target = {};
subject = {};
for i = 1:length(files)
    % Load data
    printStatus(sprintf('Loading file: %s',files{i}));
    load(files{i});
    
    % Sort gait events
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,gcStartEvent);
    
    % Gait event search range
    if isempty(geRangeTable)
        geRange = [1,height(gaitEventsSorted)];
    elseif isinf(geRangeTable(i,2)) && geRangeTable(i,1) == 1
        geRange = [1,height(gaitEventsSorted)];
    else
        start_ind = find(gaitEventsSorted.(gcStartEvent) > geRangeTable(i,1),1,'first');
        end_ind = find(gaitEventsSorted.(gcStartEvent) < geRangeTable(i,2),1,'last');
        geRange = [start_ind,end_ind];
    end
    
    % Determine what pairs to run
    elements_to_compare = [];
    if isfield(aligned_data,'left_LFP_table')
        elements_to_compare = [elements_to_compare,1:4];
    end
    
    if isfield(aligned_data,'right_LFP_table')
        elements_to_compare = [elements_to_compare,5:8];
    end
    
    pairs = nchoosek(elements_to_compare,2);
    contact_names = createContactVariableName(pairs);
    
    for j = 1:size(pairs,1)
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
        [x,~,y] = wcoherence(aligned_data.(side1).(contact1)(1:n_val),aligned_data.(side2).(contact2)(1:n_val),sr,'VoicesPerOctave',10);
        
        % Extract coherence
        printStatus(sprintf('Extracting data from pair: %s %s to %s %s',side1,contact1,side2,contact2));
        gait_cycle_mat = zeros(length(y),n_percent_bins,1);
        coherence_temp = [];
        average_coherence_gait_cycle_mat = zeros(length(y),1);
        average_coherence_temp = [];
        
        count = 1;
        for k = geRange(1):geRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)(k:k+1)) < 2)
                [~,start_ind] = min(abs(time_vec-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(time_vec-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = x(:,start_ind:end_ind);
                
                if sum(isinf(data_snip),'all') == 0
                    average_coherence_gait_cycle_mat(:,1,count) = mean(data_snip,2);
                    percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
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
            end
        end
        coherence = [coherence;coherence_temp];
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
        frequency = [frequency;repmat(round(y,4),n_gait_cycles,1)];
        gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(y))'];
        gait_cycle_percent = [gait_cycle_percent;nan(length(y)*n_gait_cycles,1)];
        state = [state;repmat({'Average'},length(y)*n_gait_cycles,1)];
        contact = [contact;repmat(contact_names(j),length(y)*n_gait_cycles,1)];
        dbs_target = [dbs_target;repmat(stim_target(i),length(y)*n_gait_cycles,1)];
        subject = [subject;repmat(subject_ID(i),length(y)*n_gait_cycles,1)];
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
data_table = table(subject,dbs_target,contact,state,gait_cycle_num,gait_cycle_percent,frequency,coherence,'VariableNames',{'Subject','Target','Contact','State','GC','GCP','Frequency','Power'});
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
    printStatus(sprintf('Gait Cycle Coherence table saved to: %s',save_path));
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