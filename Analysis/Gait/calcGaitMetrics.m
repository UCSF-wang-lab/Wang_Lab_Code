function [gait_metrics_table,cadence,gait_speed] = calcGaitMetrics(filename,varargin)
% Calculates spatial and temporal gait metrics from Delsys or Xsens data
%
% INPUTS:   
%           filesname           [=] Matrix of gait events. Must be 4
%
% OPTIONAL INPUTS:
%
%           level_type          [=] "none"/"single". Changes the way gait
%                                   events are automatically detected and
%                                   how gait metrics are calculated.
%                                   Default is "single".
%           visit_name          [=] String to denote visit type. Examples:
%                                   preop, preprogramming, 
%                                   preprogrammingUnilat,
%                                   preprogrammingBilat, dbsOpt,
%                                   dbsOptUnilat,dbsOptBilat,aDBSXXX. For
%                                   aDBS session, the XXX denotes visit
%                                   number.
%
%           trial_num           [=] Trial number if this task is repeated.
%           
%           med_state           [=] String input to denote medication state
%                                   of data. Can be on/off/low/na. Default
%                                   is na.
%
%           stim_state          [=] String input to denote stimulation
%                                   state. Can be clinical/aDBS_XXX/off/na.
%                                   If the option is aDBS_XXX, the XXX can
%                                   be changed to denote more information
%                                   about the aDBS settings. Default is na.
%
% OUTPUTS:  gait_metrics_table  [=] Table of gait events by gait cycle.
%           cadence             [=] Cadence of patient during walking.
%           gait_speed          [=] Gait speed of patient during walking.
%
% EXAMPLE:
%           calcGaitMetrics('xsens_filename',filename2,'level_type','single','sample_tol',5,'preview',true);
%
% Author:   Kenneth Louie
% Date:     06/27/2023

% Example: 
% out_table = calcGaitMetrics('xsens_filename',xsens_filename,'level_type','single','flip_data',true,'baseline_nSamples',4500,'sample_tol',3,'preview',true);
% filename = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/Pre-Processed Data/Xsens /Xsens Raw/MVNX/RCS07_Gait_12_17_19_walking_back_and_forth_after_task.csv'
% filename2 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/Pre-Processed Data/Xsens /Xsens Raw/MVNX/RCS07_Gait_12_17_19_walking_back_and_forth_before_task.csv'
%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nargin/2
    switch varargin{i*2-1}
        case 'level_type'
            level_type = varargin{i*2};
        case 'visit_name'
            visit_name = varargin{i*2};
        case 'med_state'
            med_state = varargin{i*2};
        case 'stim_state'
            stim_state = varargin{i*2};
        case 'trial_num'
            trial_num = varargin{i*2};
        case 'subjectID'
            subjectID = varargin{i*2};
        case 'single_direction'
            single_direction = varargin{i*2};
        case 'direction_turn_time_range'
            direction_turn_time_range = varargin{i*2};
        case 'save_data'
            save_data = varargin{i*2};
    end
    
end

%% Load data
if ~isempty(filename)
    load(filename);
    xsens_data = aligned_data.Xsens;
    delsys_data = aligned_data.Delsys;
    gait_events = aligned_data.gait_events;
    
    if isfield(aligned_data,'med_condition')
        med_state = aligned_data.med_condition;
    else
        med_state = 'na';
    end
    
    if isfield(aligned_data,'stim_condition')
        stim_state = aligned_data.stim_condition;
    else
        stim_state = 'na';
    end

    if isfield(aligned_data,'trial_num')
        trial_num = aligned_data.trial_num;
    end
end

if (~exist('xsens_data','var') || isempty(xsens_data)) && exist('xsens_filename','var')
    try
        xsens_data = readtable(xsens_filename);
    catch
        error('Xsens file does not exist, data is empty, or wrong filename.');
    end
end

if (~exist('delsys_data','var') || isempty(delsys_data)) && exist('delsys_filename','var')
    try
        temp_data = load(delsys_filename);
        delsys_data = temp_data.out_struct;
    catch
        error('Delsys file does not exist, data is empty, or wrong filename.');
    end
end

%% Set warnings if xsens dataset is missing
if ~exist('xsens_data','var') && exist('delsys_data','var')
    warning(sprintf('Gait Metric Calculation. \n No Xsens data detected. \n Delsys data detected \n Only temporal gait metrics will be calculated.'));
end

%% Default values
if ~exist('level_type','var') || isempty(level_type)
    level_type = 'single';
end

if ~exist('single_direction','var') || isempty(single_direction)
    single_direction = false;
end

if ~exist('direction_turn_time_range','var') || isempty(direction_turn_time_range)
    direction_turn_time_range = [0,aligned_data.Xsens.Time(end)];
end

if ~exist('visit_name','var') || isempty(visit_name)
    visit_name = 'na';
end

if ~exist('med_state','var') || isempty(med_state)
    med_state = 'na';
end

if ~exist('stim_state','var') || isempty(stim_state)
    stim_state = 'na';
end

if ~exist('trial_num','var') || isempty(trial_num)
    trial_num = 1;
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'na';
end

if ~exist('save_data','var') || isempty(save_data)
    save_data = false;
end

%% Calculate gait matrics and build output table
% out_table = array2table(nan(n_gait_events*2,11));
% out_table.Properties.VariableNames = {'Gait_Cycle','Side',...
%     'Step_Time','Step_Length','Step_Width',...
%     'Stride_Time','Stride_Length',...
%     'Swing_Time','Stance_Time',...
%     'Single_Limb_Support','Double_Limb_Support'};
% out_table.Gait_Cycle = repelem(1:n_gait_events,2)';
% out_table.Side = repmat(['R';'L'],n_gait_events,1);

if exist('xsens_data','var') && ~isempty(xsens_data)
    [gait_metrics_table,cadence,gait_speed] = getMetrics(gait_events,xsens_data,level_type,single_direction,direction_turn_time_range);
elseif ~exist('xsens_data','var') && exist('delsys_data','var') && ~isempty(delsys_data)
    [gait_metrics_table,cadence,gait_speed] = getMetrics(gait_events,delsys_data,level_type,single_direction,direction_turn_time_range);
end

SubjectID = repmat({subjectID},height(gait_metrics_table),1);
VisitName = repmat({visit_name},height(gait_metrics_table),1);
MedState = repmat({med_state},height(gait_metrics_table),1);
StimState = repmat({stim_state},height(gait_metrics_table),1);
TrialNum = repmat(trial_num,height(gait_metrics_table),1);

gait_metrics_table = addvars(gait_metrics_table,SubjectID,VisitName,MedState,StimState,TrialNum,'Before','GaitCycle');

if save_data
    % gait metrics table
    [folder,file] = fileparts(filename);
    folder = strrep(folder,"Aligned Data","Analysis Data/Gait_Metrics");
    if ~exist(folder,'dir')
        mkdir(folder);
    end

    idx = strfind(file,"w");
    file = [file(1:idx-1),'gait_metrics.csv'];

    writetable(gait_metrics_table,fullfile(folder,file));
    
    % cadence and speed
    other_table = table(subjectID,{visit_name},{med_state},{stim_state},trial_num,cadence,gait_speed,'VariableNames',{'SubjectID','VisitName','MedState','StimState','TrialNum','Cadence (steps/min)','Speed (m/s)'});
    file = strrep(file,'gait_metrics','gait_summary');
    writetable(other_table,fullfile(folder,file));

end

end

function [gait_metrics_table,cadence,gait_speed] = getMetrics(gait_events,foot_data,level_type,single_direction,direction_turn_time_range)
% [LHS,RTO,RHS,LTO]
% step length, width, time
% stride length, time
% swing time
% stance time
% double support time

if istable(foot_data)
    % Order gait events
    gait_events_turns_removed = removeGaitCyclesTurns(foot_data,gait_events);
    gait_events_ordered = sortGaitEvents(gait_events_turns_removed,'LHS');
    
    % Determine which gait cycles to consider
    gc_start_search = find(gait_events_ordered.LHS>=foot_data.Time(1),1,'first');
    gc_end_search = find(gait_events_ordered.LHS<=foot_data.Time(end),1,'last');
    gait_events_ordered_trim = gait_events_ordered(gc_start_search:gc_end_search,:);
elseif isstruct(foot_data)
    % Order gait events
    gait_events_ordered = sortGaitEvents(gait_events,'LHS');
end

if ~isempty(gait_events_ordered_trim)
    n_gait_events = size(gait_events_ordered_trim,1);
    step_length = nan(n_gait_events*2,1);
    step_time = nan(n_gait_events*2,1);
    stride_length = nan(n_gait_events*2,1);
    stride_time = nan(n_gait_events*2,1);
    swing_time = nan(n_gait_events*2,1);
    stance_time = nan(n_gait_events*2,1);
    dst = nan(n_gait_events*2,1);

    if strcmp(level_type,'none')
    end

    if strcmp(level_type,'single')
        % [LHS,RTO,RHS,LTO]
        % Step metrics using left leg as the line of progression and vector projections
        for i = 1:height(gait_events_ordered_trim)-1
            if (~isnan(gait_events_ordered_trim.RHS(i)) && ~isnan(gait_events_ordered_trim.LHS(i+1))) && ((gait_events_ordered_trim.LHS(i+1)-gait_events_ordered_trim.RHS(i))<1)
                [~,left_foot_pos1_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.LHS(i)));
                [~,right_foot_pos1_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.RHS(i)));
                [~,left_foot_pos2_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.LHS(i+1)));

                left_foot_pos1 = [foot_data.LeftFoot_PosX(left_foot_pos1_ind),foot_data.LeftFoot_PosY(left_foot_pos1_ind)];
                right_foot_pos1 = [foot_data.RightFoot_PosX(right_foot_pos1_ind),foot_data.RightFoot_PosY(right_foot_pos1_ind)];
                left_foot_pos2 = [foot_data.LeftFoot_PosX(left_foot_pos2_ind),foot_data.LeftFoot_PosY(left_foot_pos2_ind)];

                % Right step
                left_foot_seg = left_foot_pos2-left_foot_pos1;
                left_right_foot_seg = right_foot_pos1-left_foot_pos1;
                projected_point = left_foot_pos1 + dot(left_right_foot_seg,left_foot_seg)/dot(left_foot_seg,left_foot_seg)*left_foot_seg;
                step_length(i*2) = sqrt(sum((projected_point-left_foot_pos1).^2));

                % Left Stride length
                stride_length((i*2)-1) = sqrt(sum((left_foot_pos2-left_foot_pos1).^2));

                % Left Step length
                step_length((i*2)-1) = stride_length((i*2)-1)-step_length(i*2);
            end
        end

        % Step metrics using pelvis as the line of progression


        % Stride length
        % Right
        for i = 1:height(gait_events_ordered_trim)-1
            if (~isnan(gait_events_ordered_trim.RHS(i)) && ~isnan(gait_events_ordered_trim.RHS(i+1))) && ((gait_events_ordered_trim.RHS(i+1)-gait_events_ordered_trim.RHS(i))<1.5)
                [~,RHS_pre_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.RHS(i)));
                [~,RHS_post_ind] = min(abs(foot_data.Time-gait_events_ordered_trim.RHS(i+1)));
                x_axis_diff = (foot_data.RightFoot_PosX(RHS_post_ind)-foot_data.RightFoot_PosX(RHS_pre_ind))^2;
                y_axis_diff = (foot_data.RightFoot_PosY(RHS_post_ind)-foot_data.RightFoot_PosY(RHS_pre_ind))^2;
                z_axis_diff = (foot_data.RightFoot_PosZ(RHS_post_ind)-foot_data.RightFoot_PosZ(RHS_pre_ind))^2;
                stride_length(i*2) = sqrt(x_axis_diff+y_axis_diff+z_axis_diff);
            end
        end

        % Temporal metrics
        % Left
        step_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.RHS(1:end-1);     % Left step time = time from right heel strike to left heel-strike
        stride_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LHS(1:end-1);   % Left stride time = time from left heel strike to next left heel-strike
        swing_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LTO(1:end-1);    % Left swing time = Left toe-off to left heel-strike
        stance_time(1:2:n_gait_events*2) = gait_events_ordered_trim.LTO-gait_events_ordered_trim.LHS;                       % Left stance time = left heel-strike to left toe-off
        
        % Right
        step_time(2:2:(n_gait_events-1)*2) = gait_events_ordered_trim.RHS(2:end)-gait_events_ordered_trim.LHS(1:end-1);     % Left step time = time from right heel strike to left heel-strike
        stride_time(2:2:(n_gait_events-1)*2) = gait_events_ordered_trim.RHS(2:end)-gait_events_ordered_trim.RHS(1:end-1);
        swing_time(2:2:n_gait_events*2) = gait_events_ordered_trim.RHS-gait_events_ordered_trim.RTO;
        stance_time(2:2:(n_gait_events-1)*2) = gait_events_ordered_trim.RTO(2:end)-gait_events_ordered_trim.RHS(1:end-1);

        % Double support time
        dst(1:2:n_gait_events*2) = gait_events_ordered_trim.RTO-gait_events_ordered_trim.LHS;
        dst(2:2:n_gait_events*2) = gait_events_ordered_trim.LTO-gait_events_ordered_trim.RHS;
    end

    % Filter values that are outside what is normal
    step_length(step_length>0.70*mean(stride_length,'omitnan') | step_length < 0) = nan;
    step_time(step_time>1 | step_time < 0) = nan;

    stride_length(stride_length>2 | stride_length < 0) = nan;
    stride_time(stride_time>2 | stride_time < 0) = nan;

    stance_time(stance_time>1 | stance_time < 0) = nan;
    swing_time(swing_time>1 | swing_time < 0) = nan;
    dst(dst>0.5 | dst < 0) = nan;

    side = repmat({'L';'R'},n_gait_events,1);
    gait_cycle = repelem(1:n_gait_events,2)';

    % Create table
    gait_metrics_table = table(gait_cycle,side,step_length,step_time,stride_length,stride_time,swing_time,stance_time,dst,'VariableNames',{'GaitCycle','Side','StepLength','StepTime','StrideLength','StrideTime','SwingTime','StanceTime','DoubleSupportTime'});

    % Find when the patient switches direction
    if ~single_direction
        A = gait_events_ordered_trim{:,:}';
        B = rmmissing(A(:));
        ind = find(diff(B)>10);
        direction_change_time = B(ind);
        direction_change_time = direction_change_time((direction_change_time>=direction_turn_time_range(1)) & (direction_change_time<=direction_turn_time_range(2)));
    else
        direction_change_time = [];
    end

    % Calculate cadence (can only calculate with xsens right now)
    % NEED TO UPDATE CALCULATION TO BE MORE ACCURATE
    if istable(foot_data)
        if ~isempty(direction_change_time)
            n_steps = sum(gait_events_ordered_trim.LHS<=direction_change_time)+sum(gait_events_ordered_trim.RHS<=direction_change_time);
            walk_duration = max([max(gait_events_ordered_trim.LHS(gait_events_ordered_trim.LHS<=direction_change_time)),...
                max(gait_events_ordered_trim.RHS(gait_events_ordered_trim.RHS<=direction_change_time))]) - ...
                min([min(gait_events_ordered_trim.LHS(gait_events_ordered_trim.LHS<=direction_change_time)),...
                min(gait_events_ordered_trim.RHS(gait_events_ordered_trim.RHS<=direction_change_time))]);
            cadence(1) = (n_steps/walk_duration)*60;

            n_steps = sum(gait_events_ordered_trim.LHS>=direction_change_time)+sum(gait_events_ordered_trim.RHS>=direction_change_time);
            walk_duration = max([max(gait_events_ordered_trim.LHS(gait_events_ordered_trim.LHS>=direction_change_time)),...
                max(gait_events_ordered_trim.RHS(gait_events_ordered_trim.RHS>=direction_change_time))]) - ...
                min([min(gait_events_ordered_trim.LHS(gait_events_ordered_trim.LHS>=direction_change_time)),...
                min(gait_events_ordered_trim.RHS(gait_events_ordered_trim.RHS>=direction_change_time))]);
            cadence(2) = (n_steps/walk_duration)*60;
        else
            n_steps = sum(~isnan(gait_events_ordered_trim.LHS))+sum(~isnan(gait_events_ordered_trim.RHS));
            walk_duration = max([max(gait_events_ordered_trim.LHS),max(gait_events_ordered_trim.RHS)]) - min([min(gait_events_ordered_trim.LHS),min(gait_events_ordered_trim.RHS)]);
            cadence = (n_steps/walk_duration)*60;
        end
        cadence = mean(cadence);
    end


    % Calculate gait speed (can only calculate with xsens right now)
    if istable(foot_data)
        if ~isempty(direction_change_time)
            walking_distance = 0;

            walk_start_ind = find(foot_data.Time>=min(gait_events_ordered_trim{1,:}),1,'first');
            walk_end_ind = find(foot_data.Time<=direction_change_time,1,'last');
            current_posX = foot_data.Pelvis_PosX(walk_start_ind);
            current_posY = foot_data.Pelvis_PosY(walk_start_ind);
            for i = walk_start_ind+1:walk_end_ind
                deltaX = foot_data.Pelvis_PosX(i)-current_posX;
                deltaY = foot_data.Pelvis_PosY(i)-current_posY;
                deltaD = sqrt(deltaX^2 + deltaY^2);
                walking_distance = walking_distance + deltaD;

                current_posX = foot_data.Pelvis_PosX(i);
                current_posY = foot_data.Pelvis_PosY(i);
            end
            gait_speed(1) = walking_distance/(foot_data.Time(walk_end_ind)-foot_data.Time(walk_start_ind));

            walk_start_ind = find(foot_data.Time>=direction_change_time,1,'first');
            walk_end_ind = find(foot_data.Time<=max(gait_events_ordered_trim{end,:}),1,'last');
            walking_distance = 0;
            current_posX = foot_data.Pelvis_PosX(walk_start_ind);
            current_posY = foot_data.Pelvis_PosY(walk_start_ind);
            for i = walk_start_ind+1:walk_end_ind
                deltaX = foot_data.Pelvis_PosX(i)-current_posX;
                deltaY = foot_data.Pelvis_PosY(i)-current_posY;
                deltaD = sqrt(deltaX^2 + deltaY^2);
                walking_distance = walking_distance + deltaD;

                current_posX = foot_data.Pelvis_PosX(i);
                current_posY = foot_data.Pelvis_PosY(i);
            end
            gait_speed(2) = walking_distance/(foot_data.Time(walk_end_ind)-foot_data.Time(walk_start_ind));
        else
            walk_start_ind = find(foot_data.Time>=min(gait_events_ordered_trim{1,:}),1,'first');
            walk_end_ind = find(foot_data.Time<=max(gait_events_ordered_trim{end,:}),1,'last');
            walking_distance = 0;
            current_posX = foot_data.Pelvis_PosX(walk_start_ind);
            current_posY = foot_data.Pelvis_PosY(walk_start_ind);
            for i = walk_start_ind+1:walk_end_ind
                deltaX = foot_data.Pelvis_PosX(i)-current_posX;
                deltaY = foot_data.Pelvis_PosY(i)-current_posY;
                deltaD = sqrt(deltaX^2 + deltaY^2);
                walking_distance = walking_distance + deltaD;

                current_posX = foot_data.Pelvis_PosX(i);
                current_posY = foot_data.Pelvis_PosY(i);
            end
            gait_speed = walking_distance/(foot_data.Time(walk_end_ind)-foot_data.Time(walk_start_ind));
        end
        gait_speed = mean(gait_speed);
    end
else
    error("No Gait Events within this file.");
end

if mean(stride_length,'omitnan') < 0.75
    warning(sprintf("Stride length are abnormally low: %4.2fm\n Consider changing walking axis.",mean(stride_length,'omitnan')));
    figure;
    
    subplot(2,1,1); 
    plot(foot_data.Time,foot_data.LeftFoot_PosX); 
    hold on;
    plot(foot_data.Time,foot_data.RightFoot_PosX);
    title('X axis')

    subplot(2,1,2);
    plot(foot_data.Time,foot_data.LeftFoot_PosY); 
    hold on;
    plot(foot_data.Time,foot_data.RightFoot_PosY);
    title('Y axis');
end

end

function gait_events_turns_removed = removeGaitCyclesTurns(xsens_data,gait_events)
% Check to see if last row of xsens data is nan
if isnan(xsens_data{end,2})
    xsens_data(end,:) = [];
end

% Filter and convert to degrees
[b,a]=butter(4,(1.5/30));
pelvis_data = abs(filtfilt(b,a,xsens_data.Pelvis_angVelZ)*57.3);

% Find all threshold crossings in both directions
filter_threshold = 30;
pos_ind = find((pelvis_data(2:end) > filter_threshold) & (pelvis_data(1:end-1) < filter_threshold));
neg_ind = find((pelvis_data(2:end) < filter_threshold) & (pelvis_data(1:end-1) > filter_threshold));

% Filter through the positive and negative threshold crossings to remove
% erroneous crossings
remove_ind = [];
count_consideration = min([length(pelvis_data)-pos_ind(end),100]);
for i = 1:length(pos_ind)
    base_value = pelvis_data(pos_ind(i));
    for j = 1:count_consideration
        if pelvis_data(pos_ind(i)+j)<filter_threshold
            remove_ind(end+1) = i;
            break;
        end
    end
end
pos_ind(remove_ind) = [];

remove_ind = [];
count_consideration = min([length(pelvis_data)-neg_ind(end),60]);
for i = 1:length(neg_ind)
    base_value = pelvis_data(neg_ind(i));
    for j = 1:count_consideration
        if pelvis_data(neg_ind(i)+j)>filter_threshold
            remove_ind(end+1) = i;
            break;
        end
    end
end
neg_ind(remove_ind) = [];

min_length = min([length(pos_ind),length(neg_ind)]);
turn_times = [xsens_data.Time(pos_ind(1:min_length)),xsens_data.Time(neg_ind(1:min_length))];
remove_ind = [];
for i = 1:height(gait_events)
    for j = 1:size(turn_times,1)
        if any(gait_events{i,:}>=turn_times(j,1) & gait_events{i,:} <= turn_times(j,2))
            remove_ind(end+1) = i;
        end
    end
end

gait_events_turns_removed = gait_events;
gait_events_turns_removed(remove_ind,:) = [];
end