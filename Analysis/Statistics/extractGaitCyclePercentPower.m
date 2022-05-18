function data_table = extractGaitCyclePercentPower(files,subject_ID,stim_target,varargin)
%%
%   This function extracts the power in each frequency and gait cycle
%   percentage during walking and in one of two conditions, rest prior to
%   walking and the average power at a given frequency for each gait cycle.
%

% extractGaitCyclePercentPower(allFiles,subject_ID,stim_target,'geRangeTable',geRangeTable);

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

if ~exist('cycle_start_event','var')
    cycle_start_event = 'LHS';
end

if ~exist('geRangeTable','var') || isempty(geRangeTable)
    geRangeTable = [];
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
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,'LHS');
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
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
    
    if isfield(signalAnalysisData,'Left')
        walking_start_ind = find(signalAnalysisData.Left.Time{1} >= min(gaitEventsSorted{1,:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Left.Time{1} <= max(gaitEventsSorted{end,:}),1,'last');
        contact_names = createContactVariableName(signalAnalysisData.Left.Chan_Names);
        
        printStatus('Extracting data from left contacts.');
        for j = 1:length(contact_names)
            printStatus(sprintf('Extracting data from left %s.',contact_names{j}));
            
            % Walking power
            gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{j}),n_percent_bins,1);
            power_temp = [];
            average_gait_cycle_power_mat = zeros(length(signalAnalysisData.Left.Freq_Values{j}),1);
            average_power_temp = [];
            count = 1;
            for k = geRange(1):geRange(2)-1
                if ~isnan(gaitEventsSorted.(cycle_start_event)(k)) && ~isnan(gaitEventsSorted.(cycle_start_event)(k+1)) && (diff(gaitEventsSorted.(cycle_start_event)([k,k+1])) < 2)
                    [~,start_ind] = min(abs(signalAnalysisData.Left.Time{j}-gaitEventsSorted.(cycle_start_event)(k)));
                    [~,end_ind] = min(abs(signalAnalysisData.Left.Time{j}-gaitEventsSorted.(cycle_start_event)(k+1)));
                    data_snip = abs(signalAnalysisData.Left.Values{j}(:,start_ind:end_ind));
                    
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
            frequency = [frequency;repmat(repmat(round(signalAnalysisData.Left.Freq_Values{j},4),n_percent_bins,1),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signalAnalysisData.Left.Freq_Values{j})*n_percent_bins)'];
            gait_cycle_percent = [gait_cycle_percent;repmat(repelem(1:100,length(signalAnalysisData.Left.Freq_Values{j}))',n_gait_cycles,1)];
            state = [state;repmat({'Walking'},length(signalAnalysisData.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signalAnalysisData.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            side = [side;repmat({'L'},length(signalAnalysisData.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signalAnalysisData.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signalAnalysisData.Left.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            
            % Store average power at each frequency per gait cycle
            n_gait_cycles = size(average_gait_cycle_power_mat,3);
            for x = 1:n_gait_cycles
                average_power_temp = [average_power_temp;average_gait_cycle_power_mat(:,1,x)];
            end
            power = [power;average_power_temp];
            frequency = [frequency;repmat(round(signalAnalysisData.Left.Freq_Values{j},4),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signalAnalysisData.Left.Freq_Values{j}))'];
            gait_cycle_percent = [gait_cycle_percent;nan(length(signalAnalysisData.Left.Freq_Values{j})*n_gait_cycles,1)];
            state = [state;repmat({'Average'},length(signalAnalysisData.Left.Freq_Values{j})*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signalAnalysisData.Left.Freq_Values{j})*n_gait_cycles,1)];
            side = [side;repmat({'L'},length(signalAnalysisData.Left.Freq_Values{j})*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signalAnalysisData.Left.Freq_Values{j})*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signalAnalysisData.Left.Freq_Values{j})*n_gait_cycles,1)];
        end
        printStatus(' ');
    end
    
    if isfield(signalAnalysisData,'Right')
        walking_start_ind = find(signalAnalysisData.Right.Time{1} >= min(gaitEventsSorted{1,:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Right.Time{1} <= max(gaitEventsSorted{end,:}),1,'last');
        contact_names = createContactVariableName(signalAnalysisData.Right.Chan_Names);
        
        printStatus('Extracting data from right contacts.');
        for j = 1:length(contact_names)
            printStatus(sprintf('Extracting data from right %s.',contact_names{j}));
            
            % Walking power
            gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{j}),n_percent_bins,1);
            power_temp = [];
            average_gait_cycle_power_mat = zeros(length(signalAnalysisData.Right.Freq_Values{j}),1);
            average_power_temp = [];
            count = 1;
            for k = geRange(1):geRange(2)-1
                if ~isnan(gaitEventsSorted.(cycle_start_event)(k)) && ~isnan(gaitEventsSorted.(cycle_start_event)(k+1)) && (diff(gaitEventsSorted.(cycle_start_event)([k,k+1])) < 2)
                    [~,start_ind] = min(abs(signalAnalysisData.Right.Time{j}-gaitEventsSorted.(cycle_start_event)(k)));
                    [~,end_ind] = min(abs(signalAnalysisData.Right.Time{j}-gaitEventsSorted.(cycle_start_event)(k+1)));
                    data_snip = abs(signalAnalysisData.Right.Values{j}(:,start_ind:end_ind));
                    
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
            frequency = [frequency;repmat(repmat(round(signalAnalysisData.Right.Freq_Values{j},4),n_percent_bins,1),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signalAnalysisData.Right.Freq_Values{j})*n_percent_bins)'];
            gait_cycle_percent = [gait_cycle_percent;repmat(repelem(1:100,length(signalAnalysisData.Right.Freq_Values{j}))',n_gait_cycles,1)];
            state = [state;repmat({'Walking'},length(signalAnalysisData.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signalAnalysisData.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            side = [side;repmat({'R'},length(signalAnalysisData.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signalAnalysisData.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signalAnalysisData.Right.Freq_Values{j})*n_percent_bins*n_gait_cycles,1)];
            
            % Store average power at each frequency per gait cycle
            n_gait_cycles = size(average_gait_cycle_power_mat,3);
            for x = 1:n_gait_cycles
                average_power_temp = [average_power_temp;average_gait_cycle_power_mat(:,1,x)];
            end
            power = [power;average_power_temp];
            frequency = [frequency;repmat(round(signalAnalysisData.Right.Freq_Values{j},4),n_gait_cycles,1)];
            gait_cycle_num = [gait_cycle_num;repelem(1:n_gait_cycles,length(signalAnalysisData.Right.Freq_Values{j}))'];
            gait_cycle_percent = [gait_cycle_percent;nan(length(signalAnalysisData.Right.Freq_Values{j})*n_gait_cycles,1)];
            state = [state;repmat({'Average'},length(signalAnalysisData.Right.Freq_Values{j})*n_gait_cycles,1)];
            contact = [contact;repmat(contact_names(j),length(signalAnalysisData.Right.Freq_Values{j})*n_gait_cycles,1)];
            side = [side;repmat({'R'},length(signalAnalysisData.Right.Freq_Values{j})*n_gait_cycles,1)];
            dbs_target = [dbs_target;repmat(stim_target(i),length(signalAnalysisData.Right.Freq_Values{j})*n_gait_cycles,1)];
            subject = [subject;repmat(subject_ID(i),length(signalAnalysisData.Right.Freq_Values{j})*n_gait_cycles,1)];
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
    
    % Need to figure out how to save in the expanded memory format
%     writetable(data_table,save_name);
%     save(strrep(save_name,'.csv','.mat'),'data_table');
%     printStatus(sprintf('Gait Cycle Power table saved to: %s',save_path));
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