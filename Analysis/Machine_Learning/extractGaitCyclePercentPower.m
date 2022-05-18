function data_table = extractGaitCyclePercentPower(files,subject_ID,stim_target,varargin)
%%
%   This function extracts the power in each frequency and gait cycle
%   percentage during walking and in one of two conditions, rest prior to
%   walking and the average power at a given frequency for each gait cycle.
%

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
% extractGaitCyclePercentPowerV2(files,subject_ID,{'STN';'STN';'STN';'GPi';'GPi'});

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

%% Data extraction
state = {};
power = [];
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
    
    % Sort gait events and calculate CWT
    gait_events_sorted = sortGaitEvents(aligned_data.gait_events,'LHS');
    signal_analysis_data = calcRCS_CWT(aligned_data);
    
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
    
    if isfield(signal_analysis_data,'Left')
        walking_start_ind = find(signal_analysis_data.Left.Time{1} >= min(gait_events_sorted{1,:})-1,1,'first');
        walking_end_ind = find(signal_analysis_data.Left.Time{1} <= max(gait_events_sorted{end,:}),1,'last');
        contact_names = createContactVariableName(signal_analysis_data.Left.Chan_Names);
        
        printStatus('Extracting data from left contacts.');
        for j = 1:length(contact_names)
            printStatus(sprintf('Extracting data from left %s.',contact_names{j}));
            
            % Walking power
            gait_cycle_mat_left = zeros(length(signal_analysis_data.Left.Freq_Values{j}),n_percent_bins,1);
            power_temp = [];
            average_gait_cycle_power_mat = zeros(length(signal_analysis_data.Left.Freq_Values{j}),1);
            average_power_temp = [];
            count = 1;
            for k = gait_event_range(1):gait_event_range(2)-1
                if ~isnan(gait_events_sorted.(cycle_start_event)(k)) && ~isnan(gait_events_sorted.(cycle_start_event)(k+1)) && (diff(gait_events_sorted.(cycle_start_event)([k,k+1])) < 2)
                    [~,start_ind] = min(abs(signal_analysis_data.Left.Time{j}-gait_events_sorted.(cycle_start_event)(k)));
                    [~,end_ind] = min(abs(signal_analysis_data.Left.Time{j}-gait_events_sorted.(cycle_start_event)(k+1)));
                    data_snip = abs(signal_analysis_data.Left.Values{j}(:,start_ind:end_ind));
                    
                    if sum(isinf(data_snip),'all') == 0
                        average_gait_cycle_power_mat(:,1,count) = mean(data_snip,2);
                        percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                        for m = 1:length(percent_inds)-1
                            if m == 1
                                gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                            else
                                gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                            end
                        end
                        count = count + 1;
                    end
                end
            end
            
            % Store gait cycle percentages
            n_gait_cycles = size(gait_cycle_mat_left,3);
            for x = 1:n_gait_cycles
                for y = 1:n_percent_bins
                    power_temp = [power_temp;gait_cycle_mat_left(:,y,x)];
                end
            end
            power = [power;power_temp];
            frequency = [frequency;repmat(repmat(round(signal_analysis_data.Left.Freq_Values{j},4),n_percent_bins,1),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signal_analysis_data.Left.Freq_Values{j})*n_percent_bins)'];
            gait_cycle_percent = [gait_cycle_percent;repmat(repelem(1:100,length(signal_analysis_data.Left.Freq_Values{j}))',n_gait_cycles,1)];
            state = [state;repmat({'Walking'},length(signal_analysis_data.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signal_analysis_data.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            side = [side;repmat({'L'},length(signal_analysis_data.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signal_analysis_data.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signal_analysis_data.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            
            % Store average power at each frequency per gait cycle
            n_gait_cycles = size(average_gait_cycle_power_mat,3);
            for x = 1:n_gait_cycles
                average_power_temp = [average_power_temp;average_gait_cycle_power_mat(:,1,x)];
            end
            power = [power;average_power_temp];
            frequency = [frequency;repmat(round(signal_analysis_data.Left.Freq_Values{j},4),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signal_analysis_data.Left.Freq_Values{j}))'];
            gait_cycle_percent = [gait_cycle_percent;nan(length(signal_analysis_data.Left.Freq_Values{j})*n_gait_cycles,1)];
            state = [state;repmat({'Average'},length(signal_analysis_data.Left.Freq_Values{j})*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signal_analysis_data.Left.Freq_Values{j})*n_gait_cycles,1)];
            side = [side;repmat({'L'},length(signal_analysis_data.Left.Freq_Values{j})*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signal_analysis_data.Left.Freq_Values{j})*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signal_analysis_data.Left.Freq_Values{j})*n_gait_cycles,1)];
        end
        printStatus(' ');
    end
    
    if isfield(signal_analysis_data,'Right')
        walking_start_ind = find(signal_analysis_data.Right.Time{1} >= min(gait_events_sorted{1,:})-1,1,'first');
        walking_end_ind = find(signal_analysis_data.Right.Time{1} <= max(gait_events_sorted{end,:}),1,'last');
        contact_names = createContactVariableName(signal_analysis_data.Right.Chan_Names);
        
        printStatus('Extracting data from right contacts.');
        for j = 1:length(contact_names)
            printStatus(sprintf('Extracting data from right %s.',contact_names{j}));
            
            % Walking power
            gait_cycle_mat_right = zeros(length(signal_analysis_data.Right.Freq_Values{j}),n_percent_bins,1);
            power_temp = [];
            average_gait_cycle_power_mat = zeros(length(signal_analysis_data.Right.Freq_Values{j}),1);
            average_power_temp = [];
            count = 1;
            for k = gait_event_range(1):gait_event_range(2)-1
                if ~isnan(gait_events_sorted.(cycle_start_event)(k)) && ~isnan(gait_events_sorted.(cycle_start_event)(k+1)) && (diff(gait_events_sorted.(cycle_start_event)([k,k+1])) < 2)
                    [~,start_ind] = min(abs(signal_analysis_data.Right.Time{j}-gait_events_sorted.(cycle_start_event)(k)));
                    [~,end_ind] = min(abs(signal_analysis_data.Right.Time{j}-gait_events_sorted.(cycle_start_event)(k+1)));
                    data_snip = abs(signal_analysis_data.Right.Values{j}(:,start_ind:end_ind));
                    
                    if sum(isinf(data_snip),'all') == 0
                        average_gait_cycle_power_mat(:,1,count) = mean(data_snip,2);
                        percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                        for m = 1:length(percent_inds)-1
                            if m == 1
                                gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                            else
                                gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                            end
                        end
                        count = count + 1;
                    end
                end
            end
            
            % Store gait cycle percentages
            n_gait_cycles = size(gait_cycle_mat_right,3);
            for x = 1:n_gait_cycles
                for y = 1:n_percent_bins
                    power_temp = [power_temp;gait_cycle_mat_right(:,y,x)];
                end
            end
            power = [power;power_temp];
            frequency = [frequency;repmat(repmat(round(signal_analysis_data.Right.Freq_Values{j},4),n_percent_bins,1),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signal_analysis_data.Right.Freq_Values{j})*n_percent_bins)'];
            gait_cycle_percent = [gait_cycle_percent;repmat(repelem(1:100,length(signal_analysis_data.Right.Freq_Values{j}))',n_gait_cycles,1)];
            state = [state;repmat({'Walking'},length(signal_analysis_data.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signal_analysis_data.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            side = [side;repmat({'R'},length(signal_analysis_data.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signal_analysis_data.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signal_analysis_data.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            
            % Store average power at each frequency per gait cycle
            n_gait_cycles = size(average_gait_cycle_power_mat,3);
            for x = 1:n_gait_cycles
                average_power_temp = [average_power_temp;average_gait_cycle_power_mat(:,1,x)];
            end
            power = [power;average_power_temp];
            frequency = [frequency;repmat(round(signal_analysis_data.Right.Freq_Values{j},4),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signal_analysis_data.Right.Freq_Values{j}))'];
            gait_cycle_percent = [gait_cycle_percent;nan(length(signal_analysis_data.Right.Freq_Values{j})*n_gait_cycles,1)];
            state = [state;repmat({'Average'},length(signal_analysis_data.Right.Freq_Values{j})*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signal_analysis_data.Right.Freq_Values{j})*n_gait_cycles,1)];
            side = [side;repmat({'R'},length(signal_analysis_data.Right.Freq_Values{j})*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signal_analysis_data.Right.Freq_Values{j})*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signal_analysis_data.Right.Freq_Values{j})*n_gait_cycles,1)];
        end
    end
    
    % Uncomment line below to debug
    % [length(power),length(frequency),length(gait_cycle_percent),length(state),length(contact),length(side),length(dbs_target),length(subject)]
    printStatus(' ');
end
printStatus('Data extraction complete.');
printStatus(' ');
printStatus();

%% Clean up data
if remove_nan
    printStatus('Removing nan from extracted data.');
    nan_inds = ismissing(power);
    remove_inds = find(nan_inds(:,1));
    power(remove_inds,:) = [];
    subject(remove_inds) = [];
    dbs_target(remove_inds) = [];
    side(remove_inds) = [];
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
data_table = table(subject,dbs_target,side,contact,state,gait_cycle_num,gait_cycle_percent,frequency,power,'VariableNames',{'Subject','Target','Side','Contact','State','GC','GCP','Frequency','Power'});
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
    printStatus(sprintf('Gait Cycle Power table saved to: %s',save_path));
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

function contact_output_name = createContactVariableName(contact_name)
contact_output_name = cell(1,length(contact_name));

for i = 1:length(contact_name)
    plus_ind = strfind(contact_name{i},'+');
    first_space_ind = strfind(contact_name{i},' ');
    contact_output_name{i} = contact_name{i}(plus_ind:first_space_ind-1);
end
end