function [out_table,cadence,gait_speed] = calcGaitMetrics(varargin)
% Author:   Kenneth Louie
% Date:     12/16/2020
%
% Calculates spatial and temporal gait metrics from Delsys or Xsens data
%
% INPUTS:   gait_events         [=] Matrix of gait events. 4 columns in the
%                                   following order: LHS,RTO,RHS,LTO.
%           gait_event_type     [=] "samples"/"time". What type of data the
%                                   gait event matrix contains.
%           gait_event_source   [=] "xsens"/"delsys". Source of gait event
%                                   matrix. 
%           xsens_data          [=] Xsens table of data from csv file. 
%           xsens_filename      [=] Xsens csv file to read in. If
%                                   gait_events is empty, automatic gait
%                                   event detection occurs.
%           level_type          [=] "none"/"single". Changes the way gait
%                                   events are automatically detected and
%                                   how gait metrics are calculated.
%           flip_data           [=] true/false. Used to flip the 
%                                   coordinates to accurately calculate 
%                                   gait metrics. Should be set to true if 
%                                   on the xsens recording the mannequin 
%                                   body is walking into the screen from 
%                                   the origin. If walking to the left or 
%                                   right, should be set to false.
%           baseline_nSamples   [=] How many samples to ignore in the
%                                   beginning of the recording. Used to
%                                   remove false gait events if using
%                                   automatic detection.
%           sample_tol          [=] How many samples in forward/backward
%                                   direction to search for Z peak, which
%                                   helps indicate gait events. Only used
%                                   if using automatic gait event
%                                   detection.
%           preview             [=] Show a preview of automatic detected
%                                   gait events.
%
% OUTPUTS:  out_table           [=] Table of gait events by gait cycle.
%           cadence             [=] Cadence of patient during walking.
%           gait_speed          [=] Gait speed of patient during walking.
%
<<<<<<< HEAD
% EXAMPLE:
%           calcGaitMetrics('xsens_filename',filename2,'level_type','single','sample_tol',5,'preview',true);
%
% Author:   Kenneth Louie
% Date:     12/16/2020

% Example: 
% out_table = calcGaitMetrics('xsens_filename',xsens_filename,'level_type','single','flip_data',true,'baseline_nSamples',4500,'sample_tol',3,'preview',true);
% filename = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/Pre-Processed Data/Xsens /Xsens Raw/MVNX/RCS07_Gait_12_17_19_walking_back_and_forth_after_task.csv'
% filename2 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/Pre-Processed Data/Xsens /Xsens Raw/MVNX/RCS07_Gait_12_17_19_walking_back_and_forth_before_task.csv'
%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nargin/2
    switch varargin{i*2-1}
        case 'gait_events'
            gait_events = varargin{i*2};
        case 'gait_event_type'
            gait_event_type = varargin{i*2};
        case 'gait_event_source'
            gait_event_source = varargin{i*2};
        case 'xsens_data'
            xsens_data = varargin{i*2};
        case 'xsens_filename'
            xsens_filename = varargin{i*2};
        case 'level_type'
            level_type = varargin{i*2};
        case 'flip_data'
            % Used to flip the coordinate.
            flip_data = varargin{i*2};
        case 'baseline_nSamples'
            baseline_nSamples = varargin{i*2};
        case 'sample_tol'
            sample_tol = varargin{i*2};
        case 'preview'
            preview = varargin{i*2};
    end
    
end

% Load data
if (~exist('xsens_data','var') || isempty(xsens_data)) && exist('xsens_filename','var')
    try
        xsens_data = readtable(xsens_filename);
    catch
        error('Xsens file does not exist or data is empty. Cannot calculate gait metrics.');
    end
else
    error('No Xsens data found.');
end

% Default values
if ~exist('gait_event_type','var') || isempty(gait_event_type)
    gait_event_type = 'sample';
end

if ~exist('gait_event_source','var') || isempty(gait_event_source)
    gait_event_source = 'Xsens';
end

if ~exist('level_type','var') || isempty(level_type)
    level_type = 'single';
end

if ~exist('flip_data','var') || isempty(flip_data)
    flip_data = false;
end

if ~exist('baseline_nSamples','var') || isempty(baseline_nSamples)
    baseline_nSamples = 1000;
end

if ~exist('sample_tol','var') || isempty(sample_tol)
    sample_tol = 0;
end

if ~exist('preview','var') || isempty(preview)
    preview = false;
end

% Detect gait events from xsens position data if gait events input is empty
% or doesn't exist
foot_data = [xsens_data.LeftFoot_PosX,xsens_data.LeftFoot_PosY,xsens_data.LeftFoot_PosZ...
    xsens_data.RightFoot_PosX,xsens_data.RightFoot_PosY,xsens_data.RightFoot_PosZ];
if (~exist('gait_events','var') || isempty(gait_events))
    if flip_data
        temp = [foot_data(:,2),foot_data(:,5)];
        foot_data(:,2) = foot_data(:,1);
        foot_data(:,5) = foot_data(:,4);
        foot_data(:,1) = temp(:,1);
        foot_data(:,4) = temp(:,2);
    end
    
    if strcmp(level_type,'single')
        gait_events = getGaitEvents(foot_data,'single',sample_tol,baseline_nSamples,preview);
    elseif strcmp(level_type,'none')
        gait_events = getGaitEvents(foot_data,'none',sample_tol,baseline_nSamples,preview);
    end
end

% 0 indicates there was not a detected gait event there using the detector
% above. Set it to nan.
fill_ind = gait_events == 0;
gait_events(fill_ind) = nan;
n_gait_events = size(gait_events,1);

% Calculate gait matrics and build output table
out_table = array2table(nan(n_gait_events*2,11));
out_table.Properties.VariableNames = {'Gait_Cycle','Side',...
    'Step_Time','Step_Length','Step_Width',...
    'Stride_Time','Stride_Length',...
    'Swing_Time','Stance_Time',...
    'Single_Limb_Support','Double_Limb_Support'};
out_table.Gait_Cycle = repelem(1:n_gait_events,2)';
out_table.Side = repmat(['R';'L'],n_gait_events,1);

[step_time,step_length,step_width,swing_time,stance_time,sls,dls] = getStepMetrics(gait_events,gait_event_type,gait_event_source,foot_data,level_type);
out_table.Step_Time = step_time;
out_table.Step_Length = step_length;
out_table.Step_Width = step_width;
out_table.Swing_Time = swing_time;
out_table.Stance_Time = stance_time;
out_table.Single_Limb_Support = sls;
out_table.Double_Limb_Support = dls;


% Calculate stride length and time from step metrics (allows to calculate
% stride metrics regardless of level)
out_table.Stride_Time(2:2:end) = out_table.Step_Time(1:2:n_gait_events*2) + out_table.Step_Time(2:2:n_gait_events*2);   % Left
out_table.Stride_Time(3:2:end) = out_table.Step_Time(2:2:n_gait_events*2-1) + out_table.Step_Time(3:2:n_gait_events*2); % Right

out_table.Stride_Length(2:2:end) = out_table.Step_Length(1:2:n_gait_events*2) + out_table.Step_Length(2:2:n_gait_events*2);   % Left
out_table.Stride_Length(3:2:end) = out_table.Step_Length(2:2:n_gait_events*2-1) + out_table.Step_Length(3:2:n_gait_events*2); % Right

cadence = [];
gait_speed = [];

end

function gait_events = getGaitEvents(foot_data,level_type,sample_tol,baseline_nSamples,preview)
% Determine peaks and valleys
left_foot_prevX = diff(foot_data(1:end-1,1));
left_foot_prevY = diff(foot_data(1:end-1,2));
left_foot_prevZ = diff(foot_data(1:end-1,3));
left_foot_postX = diff(foot_data(2:end,1));
left_foot_postY = diff(foot_data(2:end,2));
left_foot_postZ = diff(foot_data(2:end,3));

right_foot_prevX = diff(foot_data(1:end-1,4));
right_foot_prevY = diff(foot_data(1:end-1,5));
right_foot_prevZ = diff(foot_data(1:end-1,6));
right_foot_postX = diff(foot_data(2:end,4));
right_foot_postY = diff(foot_data(2:end,5));
right_foot_postZ = diff(foot_data(2:end,6));

% Single level
if strcmp(level_type,'single')
    % Determine heel strikes by looking for flat slopes
    n_points = 8;
    for i = 1:size(foot_data,1)-n_points
        left_coeff = polyfit(1:1+n_points-1,foot_data(i:i+n_points-1,2),1);
        right_coeff = polyfit(1:1+n_points-1,foot_data(i:i+n_points-1,5),1);
        
        left_slopes(i,1) = left_coeff(1);
        right_slopes(i,1) = right_coeff(1);
    end
    
    left_baseline_Y_avg = mean(foot_data(1:1000,2));
    left_min_Y_dist = 0.15; % m
    
    right_baseline_Y_avg = mean(foot_data(1:1000,5));
    right_min_Y_dist = 0.15; % m
    
    LHS_Y = abs(left_slopes) < 1e-3 & foot_data(1:end-n_points,2) > left_baseline_Y_avg + left_min_Y_dist;
    RHS_Y = abs(right_slopes) < 1e-3 & foot_data(1:end-n_points,5) > right_baseline_Y_avg + right_min_Y_dist;
    
    LHS_Z = (left_foot_prevZ > 0 & left_foot_postZ < 0);
    LHS_Z = [0;LHS_Z(1:end-1)]; % Shift to get correct index
    
    RHS_Z = (right_foot_prevZ > 0 & right_foot_postZ < 0);
    RHS_Z = [0;RHS_Z(1:end-1)]; % Shift to get correct index
    
    LHS = zeros(size(LHS_Z,1),2+sample_tol*2);
    RHS = zeros(size(RHS_Z,1),2+sample_tol*2);
    
    LHS(:,1) = [LHS_Y;zeros(length(LHS_Z)-length(LHS_Y),1)];
    LHS(:,2) = and(LHS(:,1),LHS_Z);
    
    RHS(:,1) = [RHS_Y;zeros(length(RHS_Z)-length(RHS_Y),1)];
    RHS(:,2) = and(RHS(:,1),RHS_Z);
    
    for i = 1:sample_tol
        shifted_backward = [LHS_Z(i+1:end);zeros(i,1)];
        shifted_forward = [zeros(i,1);LHS_Z(1:end-i)];
        
        LHS(:,2+(i*2-1)) = and(LHS(:,1),shifted_backward);
        LHS(:,2+i*2) = and(LHS(:,1),shifted_forward);
        
        shifted_backward = [RHS_Z(i+1:end);zeros(i,1)];
        shifted_forward = [zeros(i,1);RHS_Z(1:end-i)];
        
        RHS(:,2+(i*2-1)) = and(RHS(:,1),shifted_backward);
        RHS(:,2+i*2) = and(RHS(:,1),shifted_forward);
    end
    
    %     mean_left_Y = abs(mean(foot_data(:,2)));
    %     std_left_Y = std(foot_data(:,2));
    
    % Remove detected heel strikes that are not above the threshold or
    % occur too soon after a previous detected heel strike
    LHS_detected = find(sum(LHS,2) >= 2);
    LHS_valid = LHS_detected;
    
    LHS_avg_stride_time_tol = ceil(0.4*60);  % ~400 ms tolerance
    LHS_remove_inds = diff(LHS_valid) < LHS_avg_stride_time_tol;
    LHS_valid([false;LHS_remove_inds]) = [];
    
    LHS_avg_stride_time = mean(diff(LHS_valid));
    LHS_std_stride_time = std(diff(LHS_valid));
    
    LHS_remove_inds2 = diff(LHS_valid) < LHS_avg_stride_time - LHS_std_stride_time;
    LHS_valid([false;LHS_remove_inds2]) = [];
    LHS_valid(LHS_valid<baseline_nSamples) = [];
    
    % Remove detected heel strikes that are not above the threshold or
    % occur too soon after a previous detected heel strike
    RHS_detected = find(sum(RHS,2) >= 2);
    RHS_valid = RHS_detected;
    
    RHS_avg_stride_time_tol = ceil(0.4*60);
    RHS_remove_inds = diff(RHS_valid) < RHS_avg_stride_time_tol;
    RHS_valid([false;RHS_remove_inds]) = [];
    
    RHS_avg_stride_time = mean(diff(RHS_valid));
    RHS_std_stride_time = std(diff(RHS_valid));
    
    RHS_remove_inds2 = diff(RHS_valid) < RHS_avg_stride_time - RHS_std_stride_time;
    RHS_valid([false;RHS_remove_inds2]) = [];
    RHS_valid(RHS_valid<baseline_nSamples) = [];
    
end

% No level (none)
if strcmp(level_type,'none')
    % LHS detection
    LHS_forward_X = (left_foot_prevX > 0 & left_foot_postX > 0);
    LHS_forward_Y = (left_foot_prevY > 0 & left_foot_postY < 0);
    LHS_forward = and(LHS_forward_X,LHS_forward_Y);
    
    LHS_backward_X = (left_foot_prevX < 0 & left_foot_postX < 0);
    LHS_backward_Y = (left_foot_prevY < 0 & left_foot_postY > 0);
    LHS_backward = and(LHS_backward_X,LHS_backward_Y);
    
    LHS_Z = (left_foot_prevZ > 0 & left_foot_postZ < 0);
    
    LHS = zeros(size(LHS_Z,1),2+sample_tol*2);
    LHS(:,1) = or(LHS_forward,LHS_backward);
    LHS(:,2) = and(LHS(:,1),LHS_Z);
    
    for i = 1:sample_tol
        shifted_backward = [LHS_Z(i+1:end);zeros(i,1)];
        shifted_forward = [zeros(i,1);LHS_Z(1:end-i)];
        
        LHS(:,2+(i*2-1)) = and(LHS(:,1),shifted_backward);
        LHS(:,2+i*2) = and(LHS(:,1),shifted_forward);
    end
    
    % Remove detected heel strikes that are not above the threshold or
    % occur too soon after a previous detected heel strike
    mean_left_Y = abs(mean(foot_data(:,2)));
    std_left_Y = std(foot_data(:,2));
    
    LHS_detected = find(sum(LHS,2) >= 2);
    LHS_valid_inds = abs(foot_data(LHS_detected+1,2)) > mean_left_Y + std_left_Y;
    LHS_valid = LHS_detected(LHS_valid_inds)+1;
    
    LHS_avg_stride_time = mean(diff(LHS_valid));
    LHS_std_stride_time = std(diff(LHS_valid));
    LHS_remove_inds = diff(LHS_valid) < LHS_avg_stride_time-LHS_std_stride_time;
    LHS_valid([false;LHS_remove_inds]) = [];
    
    % RHS detection
    RHS_forward_X = (right_foot_prevX < 0 & right_foot_postX < 0);
    RHS_forward_Y = (right_foot_prevY > 0 & right_foot_postY < 0);
    RHS_forward = and(RHS_forward_X,RHS_forward_Y);
    
    RHS_backward_X = (right_foot_prevX > 0 & right_foot_postX > 0);
    RHS_backward_Y = (right_foot_prevY < 0 & right_foot_postY > 0);
    RHS_backward = and(RHS_backward_X,RHS_backward_Y);
    
    RHS_Z = (right_foot_prevZ > 0 & right_foot_postZ < 0);
    
    RHS = zeros(size(RHS_Z,1),2+sample_tol*2);
    RHS(:,1) = or(RHS_forward,RHS_backward);
    RHS(:,2) = and(RHS(:,1),RHS_Z);
    
    for i = 1:sample_tol
        shifted_backward = [RHS_Z(i+1:end);zeros(i,1)];
        shifted_forward = [zeros(i,1);RHS_Z(1:end-i)];
        
        RHS(:,2+(i*2-1)) = and(RHS(:,1),shifted_backward);
        RHS(:,2+i*2) = and(RHS(:,1),shifted_forward);
    end
    
    % Remove detected heel strikes that are not above the threshold or
    % occur too soon after a previous detected heel strike
    mean_right_Y = abs(mean(foot_data(:,5)));
    std_right_Y = std(foot_data(:,5));
    
    RHS_detected = find(sum(RHS,2) >= 2);
    RHS_valid_inds = abs(foot_data(RHS_detected+1,5)) > mean_right_Y + std_right_Y;
    RHS_valid = RHS_detected(RHS_valid_inds)+1;
    
    RHS_avg_stride_time = mean(diff(RHS_valid));
    RHS_std_stride_time = std(diff(RHS_valid));
    RHS_remove_inds = diff(RHS_valid) < RHS_avg_stride_time-RHS_std_stride_time;
    RHS_valid([false;RHS_remove_inds]) = [];
end

% Preview gait event detected
if preview
    figure;
    subplot(1,2,1); hold on;
    plot(foot_data(:,1),foot_data(:,2));
    scatter(foot_data(LHS_valid,1),foot_data(LHS_valid,2),'^k','filled');
    xlabel('X position'); ylabel ('Y position');
    title('Left foot');
    
    subplot(1,2,2); hold on;
    plot(foot_data(:,4),foot_data(:,5));
    scatter(foot_data(RHS_valid,4),foot_data(RHS_valid,5),'^k','filled');
    xlabel('X position'); ylabel ('Y position');
    title('Right foot');
    
    sgtitle('Cycle plots');
    
    figure;
    h(1) = subplot(2,1,1); hold on;
    plot(foot_data(:,1))
    plot(foot_data(:,2))
    plot(foot_data(:,3));
    scatter(LHS_valid,foot_data(LHS_valid,2),'^k','filled')
    title('Left foot');
    
    h(2) = subplot(2,1,2); hold on;
    plot(foot_data(:,4))
    plot(foot_data(:,5))
    plot(foot_data(:,6));
    scatter(RHS_valid,foot_data(RHS_valid,5),'^k','filled')
    title('Right foot');
    
    sgtitle('Sample plots');
    
    linkaxes(h,'x');
end

% Sort detected heel strikes into gait cycles
n_LHS = length(LHS_valid);
n_RHS = length(RHS_valid);
gait_events = [];

left_count = 1;
right_count = 1;
gait_cycle_count = 1;
gait_cycle_tol = ceil(mean([LHS_avg_stride_time,RHS_avg_stride_time])*0.75);
loop_max = 1000;
loop_count = 1;
while (left_count <= n_LHS || right_count <= n_RHS ) && loop_count < loop_max
    if left_count <= n_LHS
        cur_left = LHS_valid(left_count);
    else
        cur_left = [];
    end
    
    if right_count <= n_RHS
        cur_right = RHS_valid(right_count);
    else
        cur_right = [];
    end
    
    if ~isempty(cur_left) && ~isempty(cur_right)
        if cur_right < cur_left
            gait_events(gait_cycle_count,2) = cur_right;
            right_count = right_count + 1;
            gait_cycle_count = gait_cycle_count + 1;
        elseif cur_left < cur_right && abs(cur_right-cur_left) < gait_cycle_tol
            gait_events(gait_cycle_count,:) = [cur_left,cur_right];
            left_count = left_count + 1;
            right_count = right_count + 1;
            gait_cycle_count = gait_cycle_count + 1;
        else
            gait_events(gait_cycle_count,1) = cur_left;
            left_count = left_count + 1;
            gait_cycle_count = gait_cycle_count + 1;
        end
    elseif ~isempty(cur_left) && isempty(cur_right)
        gait_events(gait_cycle_count,1) = cur_left;
        left_count = left_count + 1;
        gait_cycle_count = gait_cycle_count + 1;
    else
        gait_events(gait_cycle_count,2) = cur_right;
        right_count = right_count + 1;
        gait_cycle_count = gait_cycle_count + 1;
    end
    
    loop_count = loop_count + 1;
end

gait_events = [gait_events(:,1),nan(size(gait_events,1),1),gait_events(:,2),nan(size(gait_events,1),1)];
end

function [step_time,step_length,step_width,swing_time,stance_time,sls,dls] = getStepMetrics(gait_events,gait_event_type,gait_event_source,foot_data,level_type)
% [LHS,RTO,RHS,LTO]

% Right step length first since LHS start gait cycle convention
n_gait_events = size(gait_events,1);
step_time = nan(n_gait_events*2,1);
step_length = nan(n_gait_events*2,1);
step_width = nan(n_gait_events*2,1);
swing_time = nan(n_gait_events*2,1);
stance_time = nan(n_gait_events*2,1);
sls = nan(n_gait_events*2,1);
dls = nan(n_gait_events*2,1);

if strcmp(gait_event_type,'time')
    % Convert to samples
    gait_events_spatial = round(gait_events.*60);
    time_correction = 1;
elseif strcmp(gait_event_type,'sample')
    if strcmpi(gait_event_source,'Delsys')
        error('Gait events from Delsys should be passed in as time values.');
    elseif strcmpi(gait_event_source,'Xsens')
        gait_events_spatial = gait_events;
    end
    time_correction = 60;
end

if strcmp(level_type,'none')
    for i = 1:n_gait_events
        % i*2 = left metrics
        % i*2-1 = right metrics
        
        % out of bounds
        if i+1 > n_gait_events
            step_time(i*2) = nan;
            swing_time(i*2) = nan;
            sls(i*2) = nan;
            stance_time(i*2-1) = nan;
        else
            step_time(i*2) = (gait_events(i+1,1)-gait_events(i,3))/time_correction;
            swing_time(i*2) = (gait_events(i+1,1)-gait_events(i,4))/time_correction;
            sls(i*2) = (gait_events(i+1,1)-gait_events(i,4))/time_correction;
            stance_time(i*2-1) = (gait_events(i+1,2)-gait_events(i,3))/time_correction;
            
            if isnan(gait_events(i+1,1))
                step_length(i*2) = nan;
                step_width(i*2) = nan;
            else
                step_length(i*2) = abs(foot_data(gait_events_spatial(i+1,1),2)-foot_data(gait_events_spatial(i+1,1),5));
                step_width(i*2) = abs(foot_data(gait_events_spatial(i+1,1),1)-foot_data(gait_events_spatial(i+1,1),4));
            end
        end
        
        % Can't index with nan
        if isnan(gait_events(i,3))
            step_length(i*2-1) = nan;
            step_width(i*2-1) = nan;
        else
            step_length(i*2-1) = abs(foot_data(gait_events_spatial(i,3),5)-foot_data(gait_events_spatial(i,3),2));  % Right
            step_width(i*2-1) = abs(foot_data(gait_events_spatial(i,3),4)-foot_data(gait_events_spatial(i,3),1));   % Right
        end
        
        % Only looks at current gait cycle and doesn't matter if nan
        stance_time(i*2) = (gait_events(i,4)-gait_events(i,1))/time_correction;
        dls(i*2) = (gait_events(i,4)-gait_events(i,3))/time_correction;
        step_time(i*2-1) = (gait_events(i,3)-gait_events(i,1))/time_correction;
        swing_time(i*2-1) = (gait_events(i,3)-gait_events(i,2))/time_correction;
        sls(i*2-1) = (gait_events(i,3)-gait_events(i,2))/time_correction;
        dls(i*2-1) = (gait_events(i,2)-gait_events(i,1))/time_correction;
    end
end

if strcmp(level_type,'single')
    for i = 1:n_gait_events
        if i+1 > n_gait_events || isnan(gait_events(i+1,1)) || isnan(gait_events(i,3))
            step_time(i*2) = nan;
            step_length(i*2) = nan;
            step_width(i*2)= nan;
        else
            step_time(i*2) = (gait_events(i+1,1)-gait_events(i,3))/time_correction;                 % Left
            step_length(i*2) = abs(foot_data(gait_events_spatial(i+1,1),2)-foot_data(gait_events_spatial(i,3),5));  % Left
            step_width(i*2) = abs(foot_data(gait_events_spatial(i+1,1),1)-foot_data(gait_events_spatial(i,3),4));   % Left
        end
        
        
        if isnan(gait_events(i,3)) || isnan(gait_events(i,1))
            step_time(i*2-1) = nan;
            step_length(i*2-1) = nan;
            step_width(i*2-1) = nan;
        else
            step_time(i*2-1) = (gait_events(i,3)-gait_events(i,1))/time_correction;                 % Right
            step_length(i*2-1) = abs(foot_data(gait_events_spatial(i,3),5)-foot_data(gait_events_spatial(i,1),2));  % Right
            step_width(i*2-1) = abs(foot_data(gait_events_spatial(i,3),4)-foot_data(gait_events_spatial(i,1),1));   % Right
        end
        
        if i+1 < n_gait_events
            swing_time(i*2) = (gait_events(i+1,1)-gait_events(i,4))/time_correction;     % Left
            sls(i*2) = (gait_events(i+1,1)-gait_events(i,4))/time_correction;            % Left
            stance_time(i*2-1) = (gait_events(i+1,2)-gait_events(i,3))/time_correction;  % Right
        end
        
        stance_time(i*2) = (gait_events(i,4)-gait_events(i,1))/time_correction;      % Left
        dls(i*2) = (gait_events(i,4)-gait_events(i,3))/time_correction;              % Left
        
        swing_time(i*2-1) = (gait_events(i,3)-gait_events(i,2))/time_correction;     % Right
        sls(i*2-1) = (gait_events(i,3)-gait_events(i,2))/time_correction;            % Right
        dls(i*2-1) = (gait_events(i,2)-gait_events(i,1))/time_correction;            % Right
    end
end
end