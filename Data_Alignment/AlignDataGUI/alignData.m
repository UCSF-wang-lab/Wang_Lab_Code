function alignData(src,event)
main = src.Parent.Parent;
aligned_data = [];

if main.UserData.pre_alignment_time < 0
    pre_align_time = 0;
    addEvent('Cannot have a negative pre alignment time. Setting pre alignment time to 0.');
else
    pre_align_time = main.UserData.pre_alignment_time;
end

if main.UserData.post_alignment_time < 0
    post_align_time = inf;
    addEvent('Cannot have a negative post alignment time. Setting post alignment time to end of data file.');
else
    post_align_time = main.UserData.post_alignment_time;
end

% Add event
addEvent('Aligning data using marked points...');

%% Trim RC+S data
% Left LFP
if ~isempty(main.UserData.LFP_data) && isfield(main.UserData.LFP_data,'Left')
    addEvent('Left LFP...');
    
    if ~isnan(main.UserData.alignment_times(1))         % Left LFP mark
        align_time = main.UserData.alignment_times(1);
    elseif ~isnan(main.UserData.alignment_times(3))     % Left Accel mark
        temp1 = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        temp2 = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        alignment_adjust = seconds(temp2-temp1);
        
        align_time = main.UserData.alignment_times(3) + alignment_adjust;
    else
        addEvent('No alignment point for Left RC+S data.');
    end
    
    lfp_time = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    lfp_time = seconds(lfp_time-lfp_time(1));
    [left_taxis,left_LFP_table] = trimData(lfp_time,main.UserData.LFP_data.Left.timeDomainDataTable,[align_time-pre_align_time,align_time+post_align_time]);
        
    aligned_data.left_taxis = left_taxis;
    aligned_data.left_LFP_table = left_LFP_table;
    
    addEvent('Complete!',1);
end

% Left Accel
if ~isempty(main.UserData.Accel_data) && isfield(main.UserData.Accel_data,'Left')
    addEvent('Left Acceleration...');
    
    if ~isnan(main.UserData.alignment_times(3))         % Left Accel mark
        align_time = main.UserData.alignment_times(3);
    elseif ~isnan(main.UserData.alignment_times(1))     % Left LFP mark
        temp1 = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        temp2 = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        alignment_adjust = seconds(temp2-temp1);
        
        align_time = main.UserData.alignment_times(1) + alignment_adjust;
    else
        addEvent('No alignment point for Left RC+S data.');
    end
    
    accel_time = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    accel_time = seconds(accel_time-accel_time(1));
    [left_accel_taxis,left_Accel_table] = trimData(accel_time,main.UserData.Accel_data.Left.accelDataTable,[align_time-pre_align_time,align_time+post_align_time]);
        
    aligned_data.left_accel_taxis = left_accel_taxis;
    aligned_data.left_Accel_table = left_Accel_table;
    
    addEvent('Complete!',1);
end

% Right LFP
if ~isempty(main.UserData.LFP_data) && isfield(main.UserData.LFP_data,'Right')
    addEvent('Right LFP...');
    
    if ~isnan(main.UserData.alignment_times(2))         % Right LFP mark
        align_time = main.UserData.alignment_times(2);
    elseif ~isnan(main.UserData.alignment_times(4))     % Right Accel mark
        temp1 = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        temp2 = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        alignment_adjust = seconds(temp2-temp1);
        
        align_time = main.UserData.alignment_times(4) + alignment_adjust;
    else
        addEvent('No alignment point for Right RC+S data.');
    end
    
    lfp_time = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    lfp_time = seconds(lfp_time-lfp_time(1));
    [right_taxis,right_LFP_table] = trimData(lfp_time,main.UserData.LFP_data.Right.timeDomainDataTable,[align_time-pre_align_time,align_time+post_align_time]);
        
    aligned_data.right_taxis = right_taxis;
    aligned_data.right_LFP_table = right_LFP_table;
    
    addEvent('Complete!',1);
end

% Right Accel
if ~isempty(main.UserData.Accel_data) && isfield(main.UserData.Accel_data,'Right')
    addEvent('Right Acceleration...');
    
    if ~isnan(main.UserData.alignment_times(4))         % Right Accel mark
        align_time = main.UserData.alignment_times(4);
    elseif ~isnan(main.UserData.alignment_times(2))     % Right LFP mark
        temp1 = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        temp2 = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        alignment_adjust = seconds(temp2-temp1);
        
        align_time = main.UserData.alignment_times(2) + alignment_adjust;
    else
        addEvent('No alignment point for Right RC+S data.');
    end
    
    accel_time = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
    accel_time = seconds(accel_time-accel_time(1));
    [right_accel_taxis,right_Accel_table] = trimData(accel_time,main.UserData.Accel_data.Right.accelDataTable,[align_time-pre_align_time,align_time+post_align_time]);
        
    aligned_data.right_accel_taxis = right_accel_taxis;
    aligned_data.right_Accel_table = right_Accel_table;
    
    addEvent('Complete!',1);
end

%% Align adaptive data to LFP and Acceleration
addEvent('Checking for adaptive data...');
if isfield(main.UserData,'LogTable')
    % Left adaptive data
    if isfield(main.UserData.LogTable,'Left')
        if height(main.UserData.LogTable.Left.LogTable.adaptive) > 0
            addEvent('Left adaptive data...');
            % Time adjustment to alignement source
            if ~isnan(main.UserData.alignment_times(1))         % Left LFP mark
                temp1 = datetime(main.UserData.LogTable.Left.LogTable.adaptive.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                alignment_adjust = seconds(temp2-temp1);
                
                align_time = main.UserData.alignment_times(1) + alignment_adjust;
                
            elseif ~isnan(main.UserData.alignment_times(3))     % Left Accel mark
                temp1 = datetime(main.UserData.LogTable.Left.LogTable.adaptive.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                alignment_adjust = seconds(temp2-temp1);
                
                align_time = main.UserData.alignment_times(3) + alignment_adjust;
            else
                addEvent('No alignment point for Left RC+S data.');
            end
            
            adaptive_time = datetime(main.UserData.LogTable.Left.LogTable.adaptive.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            adaptive_time = seconds(adaptive_time-adaptive_time(1));
            [left_adaptive_taxis,left_adaptive_table] = trimData(adaptive_time,main.UserData.LogTable.Left.LogTable.adaptive,[align_time-pre_align_time,align_time+post_align_time]);
            
            aligned_data.left_adaptive_taxis = left_adaptive_taxis;
            aligned_data.left_adaptive_table = left_adaptive_table;
            
            addEvent('Complete!',1);
        end
        
    % Right adpative data
    elseif isfield(main.UserData.LogTable,'Right')
        if height(main.UserData.LogTable.Right.LogTable.adaptive) > 0
            addEvent('Right adaptive data...');
            % Time adjustment ot alignment source
            if ~isnan(main.UserData.alignment_times(2))         % Right LFP mark
                temp1 = datetime(main.UserData.LogTable.Right.LogTable.adaptive.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                alignment_adjust = seconds(temp2-temp1);
                
                align_time = main.UserData.alignment_times(2) + alignment_adjust;
                
            elseif ~isnan(main.UserData.alignment_times(4))     % Right Accel mark
                temp1 = datetime(main.UserData.LogTable.Right.LogTable.adaptive.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                alignment_adjust = seconds(temp2-temp1);
                
                align_time = main.UserData.alignment_times(4) + alignment_adjust;
            else
                addEvent('No alignment point for Right RC+S data.');
            end
            
            adaptive_time = datetime(main.UserData.LogTable.Right.LogTable.adaptive.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
            adaptive_time = seconds(adaptive_time-adaptive_time(1));
            [right_adaptive_taxis,right_adaptive_table] = trimData(adaptive_time,main.UserData.LogTable.Right.LogTable.adaptive,[align_time-pre_align_time,align_time+post_align_time]);
            
            aligned_data.right_adaptive_taxis = right_adaptive_taxis;
            aligned_data.right_adaptive_table = right_adaptive_table;
            
            addEvent('Complete!',1);
        end
    end
end

%% Trim Other data
% Delsys (if not marked, assumed same time as Xsens or Force plate)
if ~isempty(main.UserData.Delsys_data)
    addEvent('Delsys...');
    if ~isnan(main.UserData.alignment_times(5))
        align_time = main.UserData.alignment_times(5);
    elseif ~isnan(main.UserData.alignment_times(6))
        align_time = main.UserData.alignment_times(6);
    elseif ~isnan(main.UserData.ailgnment_times(7))
        align_time = main.UserData.alignment_times(7);
    else
        addEvent('Cannot align. Must have Delsys, Xsens, or Force plate alignment point marked.',1);
        return;
    end
    
    names = fieldnames(main.UserData.Delsys_data.out_struct.Time);
    for i = 1:length(names)
        time_vec = main.UserData.Delsys_data.out_struct.Time.(names{i});
        data = main.UserData.Delsys_data.out_struct.Data.(names{i});
        [new_time,new_data] = trimData(time_vec,data,[align_time-pre_align_time,align_time+post_align_time]);
        Delsys.Time.(names{i}) = new_time;
        Delsys.Data.(names{i}) = new_data;
    end
    
    Delsys.srates = main.UserData.Delsys_data.out_struct.srates;
    Delsys.Chan_names = main.UserData.Delsys_data.out_struct.Chan_names;
    addEvent('Complete!',1);
    
    aligned_data.Delsys = Delsys;
end

% Xsens (if not marked, assumed same time as Force plate or Delsys)
if ~isempty(main.UserData.Xsens_data)
    addEvent('Xsens...');
    if ~isnan(main.UserData.alignment_times(6))
        align_time = main.UserData.alignment_times(6);
    elseif ~isnan(main.UserData.alignment_times(5))
        align_time = main.UserData.alignment_times(5);
    elseif ~isnan(main.UserData.alignment_times(7))
        align_time = main.UserData.alignment_times(7);
    else
        addEvent('Cannot align. Must have Delsys, Xsens, or Force plate alignment point marked.',1);
        return;
    end
    
    [~,Xsens] = trimData([],main.UserData.Xsens_data,[align_time-pre_align_time,align_time+post_align_time]);
    
    addEvent('Complete!',1);
    
    aligned_data.Xsens = Xsens;
end

% Force plate (if not marked, assumed same time as Xsens or Delsys)
if ~isempty(main.UserData.FP_data)
    addEvent('Force plate...');
    if ~isnan(main.UserData.alignment_times(7))
        align_time = main.UserData.alignment_times(7);
    elseif ~isnan(main.UserData.alignment_times(5))
        align_time = main.UserData.alignment_times(5);
    elseif ~isnan(main.UserData.alignment_times(6))
        align_time = main.UserData.alignment_times(6);
    else
        addEvent('Cannot align. Must have Delsys, Xsens, or Force plate alignment point marked.',1);
        return;
    end
    
    time_vec = (0:height(main.UserData.FP_data)-1)/1000;
    [FP_time,FP] = trimData(time_vec,main.UserData.FP_data,[align_time-pre_align_time,align_time+post_align_time]);
    FP.Time = FP_time';
    
    addEvent('Complete!',1);
    
    aligned_data.FP = FP;
end

% Rover
if ~isempty(main.UserData.Rover_data)
    addEvent('Left Rover...');
    if ~isnan(main.UserData.alignment_times(8))
        align_time = main.UserData.alignment_times(8);
    else
        addEvent('No alignment point for Left Rover data.',1);
        return;
    end

    [~,Rover] = trimData([],main.UserData.Rover_data.Left,[align_time,pre_align_time,post_align_time]);
    addEvent('Complete!',1);
    aligned_data.Rover.Left = Rover;

    addEvent('Right Rover...');
    if ~isnan(main.UserData.alignment_times(9))
        align_time = main.UserData.alignment_times(9);
    else
        addEvent('No alignment point for Right Rover data.',1);
        return;
    end
    [~,Rover] = trimData([],main.UserData.Rover_data.Right,[align_time,pre_align_time,post_align_time]);
    addEvent('Complete!',1);
    aligned_data.Rover.Right = Rover;
end

%% Add alignment to UserData struct
addEvent('Alignment complete!');
main.UserData.aligned_data = aligned_data;
end

%% HELPER FUNCTIONS
function varargout = trimData(time_vec,data,alignment_time)

if ~isempty(time_vec) && ~isempty(data)
    if alignment_time(1)<0
        start_alignment_ind = 1;
    else
        [~,start_alignment_ind] = min(abs(time_vec-alignment_time(1)));
    end
    
    if isinf(alignment_time(2))
        end_alignment_ind = length(time_vec);
    else
        [~,end_alignment_ind] = min(abs(time_vec-alignment_time(2)));
    end
    
    aligned_time_vec = time_vec(start_alignment_ind:end_alignment_ind)-alignment_time(1);
    aligned_data = data(start_alignment_ind:end_alignment_ind,:);
elseif isempty(time_vec) && istable(data) && sum(ismember(data.Properties.VariableNames,'Time')) == 1
    if alignment_time(1)<0
        start_alignment_ind = 1;
    else
        [~,start_alignment_ind] = min(abs(data.Time-alignment_time(1)));
    end
    
    if isinf(alignment_time(2))
        if ~isempty(time_vec) && ~isempty(data)
            end_alignment_ind = length(time_vec);
        else
            end_alignment_ind = length(data.Time);
        end
    else
        [~,end_alignment_ind] = min(abs(data.Time-alignment_time(2)));
    end
    
    aligned_time_vec = [];
    aligned_data = data(start_alignment_ind:end_alignment_ind,:);
    aligned_data.Time = aligned_data.Time-alignment_time(1);
elseif isempty(time_vec) && istable(data) && sum(ismember(data.Properties.VariableNames,'DateTime')) == 1   % Special Rover data case.
    % alignment_time variable has three entries, the alignment point as a
    % data index and the pre and post time in seconds.
    align_time = data.DateTime(alignment_time(1));
    
    % Calculate pre and post alignment time
    pre_align_time = align_time-seconds(alignment_time(2));
    post_align_time = align_time+seconds(alignment_time(3));

    % Calculate pre and post align ind
    [~,pre_align_ind] = min(abs(data.DateTime-pre_align_time));
    if ~isinf(post_align_time)
        [~,post_align_ind] = min(abs(data.DateTime-post_align_time));
    else
        post_align_ind = height(data);
    end

    aligned_time_vec = [];
    aligned_data = data(pre_align_ind:post_align_ind,:);
end

varargout{1} = aligned_time_vec;
varargout{2} = aligned_data;
end