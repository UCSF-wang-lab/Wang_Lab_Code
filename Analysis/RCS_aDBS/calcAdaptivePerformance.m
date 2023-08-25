function performanceOutput = calcAdaptivePerformance(data,phase,flipped_state)
% phase input can be contralateral swing (cSwing), ipsilateral swing
% (iSwing), contralateral stance (cStance), or ipsilateral stance (iStance)
%
% flipped_state is a 1x2 vector. The first element determines if the the
% state is flipped for the left side, and the second element determines if
% the state is flipped for the right side.

switch phase
    case 'cSwing'
        gait_events_left = sortGaitEvents(data.gait_events,'RTO');
        gait_events_right = sortGaitEvents(data.gait_events,'LTO');
    case 'iSwing'
        gait_events_left = sortGaitEvents(data.gait_events,'LTO');
        gait_events_right = sortGaitEvents(data.gait_events,'RTO');
    case 'cStance'
        gait_events_left = sortGaitEvents(data.gait_events,'RHS');
        gait_events_right = sortGaitEvents(data.gait_events,'LHS');
    case 'iStance'
        gait_events_left = sortGaitEvents(data.gait_events,'LHS');
        gait_events_right = sortGaitEvents(data.gait_events,'RHS');
    case 'dst'
        gait_events_left = sortGaitEvents(data.gait_events,'LHS');
        gait_events_right = sortGaitEvents(data.gait_events,'RHS');
end

if isfield(data,'left_adaptive_taxis')
    left_correct_state = zeros(height(data.left_adaptive_table),1);
    for i = 1:height(data.left_adaptive_table)
        switch phase
            case 'cSwing'
                if sum((data.left_adaptive_taxis(i) >= gait_events_left.RTO) & (data.left_adaptive_taxis(i) <= gait_events_left.RHS)) > 0
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        left_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        left_correct_state(i) = 1;
                    end
                end
            case 'iSwing'
                if sum((data.left_adaptive_taxis(i) >= gait_events_left.LTO) & (data.left_adaptive_taxis(i) <= gait_events_left.LHS)) > 0
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        left_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        left_correct_state(i) = 1;
                    end
                end
            case 'cStance'
                if sum((data.left_adaptive_taxis(i) >= gait_events_left.RHS) & (data.left_adaptive_taxis(i) <= gait_events_left.RTO)) > 0
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        left_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        left_correct_state(i) = 1;
                    end
                end
            case 'iStance'
                if sum((data.left_adaptive_taxis(i) >= gait_events_left.LHS) & (data.left_adaptive_taxis(i) <= gait_events_left.LTO)) > 0
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        left_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        left_correct_state(i) = 1;
                    end
                end
            case 'dst'
                if sum((data.left_adaptive_taxis(i) >= gait_events_left.LHS) & (data.left_adaptive_taxis(i) <= gait_events_left.RTO)) > 0
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        left_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.left_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        left_correct_state(i) = 1;
                    end
                end
        end
    end
    
    if flipped_state(1) == 1
        left_correct_state = ~left_correct_state;
    end
    
    performanceOutput.leftAccuracy = sum(left_correct_state)/length(left_correct_state);
    performanceOutput.leftStates = left_correct_state;
end

if isfield(data,'right_adaptive_taxis')
    right_correct_state = zeros(height(data.right_adaptive_table),1);
    for i = 1:height(data.right_adaptive_table)
        switch phase
            case 'cSwing'
                if sum((data.right_adaptive_taxis(i) >= gait_events_right.LTO) & (data.right_adaptive_taxis(i) <= gait_events_right.LHS)) > 0
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        right_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        right_correct_state(i) = 1;
                    end
                end
            case 'iSwing'
                if sum((data.right_adaptive_taxis(i) >= gait_events_right.RTO) & (data.right_adaptive_taxis(i) <= gait_events_right.RHS)) > 0
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        right_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        right_correct_state(i) = 1;
                    end
                end
            case 'cStance'
                if sum((data.right_adaptive_taxis(i) >= gait_events_right.LHS) & (data.right_adaptive_taxis(i) <= gait_events_right.LTO)) > 0
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        right_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        right_correct_state(i) = 1;
                    end
                end
            case 'iStance'
                if sum((data.right_adaptive_taxis(i) >= gait_events_right.RHS) & (data.right_adaptive_taxis(i) <= gait_events_right.RTO)) > 0
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        right_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        right_correct_state(i) = 1;
                    end
                end
            case 'dst'
                if sum((data.right_adaptive_taxis(i) >= gait_events_right.RHS) & (data.right_adaptive_taxis(i) <= gait_events_right.LTO)) > 0
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 1')
                        right_correct_state(i) = 1;
                    end
                else
                    if strcmp(data.right_adaptive_table.CurrentAdaptiveState(i),'State 0')
                        right_correct_state(i) = 1;
                    end
                end
        end
    end
    
    if flipped_state(2) == 1
        right_correct_state = ~right_correct_state;
    end
    
    performanceOutput.rightAccuracy = sum(right_correct_state)/length(right_correct_state);
    performanceOutput.rightStates = right_correct_state;
end

% Calculate gait phase accuracy
% Left
if isfield(data,'left_adaptive_taxis')
    start_gait_cycle = find(data.left_adaptive_taxis(1)<gait_events_left{:,1},1,'first');
    end_gait_cycle = find(data.left_adaptive_taxis(end)<gait_events_left{:,1},1,'first');
    
    if isempty(start_gait_cycle)
        start_gait_cycle = 1;
    end
    
    if isempty(end_gait_cycle)
        end_gait_cycle = height(gait_events_left);
    end
    
    left_correct_gait_phase = nan(length(start_gait_cycle:end_gait_cycle),1);
    for i = start_gait_cycle:end_gait_cycle
        switch phase
            case 'cSwing'
                gait_phase_start_ind = find(data.left_adaptive_taxis>gait_events_left.RTO(i),1,'first');
                gait_phase_end_ind = find(data.left_adaptive_taxis<gait_events_left.RHS(i),1,'last');
                
                adaptive_states = data.left_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    left_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    left_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
                
            case 'iSwing'
                gait_phase_start_ind = find(data.left_adaptive_taxis>gait_events_left.LTO(i),1,'first');
                gait_phase_end_ind = find(data.left_adaptive_taxis<gait_events_left.LHS(i),1,'last');
                
                adaptive_states = data.left_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    left_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    left_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
            case 'cStance'
                gait_phase_start_ind = find(data.left_adaptive_taxis>gait_events_left.RHS(i),1,'first');
                gait_phase_end_ind = find(data.left_adaptive_taxis<gait_events_left.RTO(i),1,'last');
                
                adaptive_states = data.left_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    left_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    left_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
            case 'iStance'
                gait_phase_start_ind = find(data.left_adaptive_taxis>gait_events_left.LHS(i),1,'first');
                gait_phase_end_ind = find(data.left_adaptive_taxis<gait_events_left.LTO(i),1,'last');
                
                adaptive_states = data.left_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    left_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    left_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
            case 'dst'
                gait_phase_start_ind = find(data.left_adaptive_taxis>gait_events_left.LHS(i),1,'first');
                gait_phase_end_ind = find(data.left_adaptive_taxis<gait_events_left.RTO(i),1,'last');
                
                adaptive_states = data.left_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    left_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    left_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
        end
    end
    performanceOutput.leftPhaseAccuracy = sum(left_correct_gait_phase)/length(left_correct_gait_phase);
    performanceOutput.leftPhaseStates = left_correct_gait_phase;
end

% Right
if isfield(data,'right_adaptive_taxis')
    start_gait_cycle = find(data.right_adaptive_taxis(1)<gait_events_right{:,1},1,'first');
    end_gait_cycle = find(data.right_adaptive_taxis(end)<gait_events_right{:,1},1,'first');
    
    if isempty(start_gait_cycle)
        start_gait_cycle = 1;
    end
    
    if isempty(end_gait_cycle)
        end_gait_cycle = height(gait_events_right);
    end
    
    right_correct_gait_phase = nan(length(start_gait_cycle:end_gait_cycle),1);
    for i = start_gait_cycle:end_gait_cycle
        switch phase
            case 'cSwing'
                gait_phase_start_ind = find(data.right_adaptive_taxis>gait_events_right.LTO(i),1,'first');
                gait_phase_end_ind = find(data.right_adaptive_taxis<gait_events_right.LHS(i),1,'last');
                
                adaptive_states = data.right_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    right_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    right_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
                
            case 'iSwing'
                gait_phase_start_ind = find(data.right_adaptive_taxis>gait_events_right.RTO(i),1,'first');
                gait_phase_end_ind = find(data.right_adaptive_taxis<gait_events_right.RHS(i),1,'last');
                
                adaptive_states = data.right_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    right_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    right_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
            case 'cStance'
                gait_phase_start_ind = find(data.right_adaptive_taxis>gait_events_right.LHS(i),1,'first');
                gait_phase_end_ind = find(data.right_adaptive_taxis<gait_events_right.LTO(i),1,'last');
                
                adaptive_states = data.right_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    right_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    right_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
            case 'iStance'
                gait_phase_start_ind = find(data.right_adaptive_taxis>gait_events_right.RHS(i),1,'first');
                gait_phase_end_ind = find(data.right_adaptive_taxis<gait_events_right.RTO(i),1,'last');
                
                adaptive_states = data.right_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    right_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    right_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
            case 'dst'
                gait_phase_start_ind = find(data.right_adaptive_taxis>gait_events_right.RHS(i),1,'first');
                gait_phase_end_ind = find(data.right_adaptive_taxis<gait_events_right.LTO(i),1,'last');
                
                adaptive_states = data.right_adaptive_table.CurrentAdaptiveState(gait_phase_start_ind:gait_phase_end_ind);
                
                if sum(strcmp(adaptive_states,'State 1')) > 0
                    right_correct_gait_phase(i-start_gait_cycle+1) = 1;
                else
                    right_correct_gait_phase(i-start_gait_cycle+1) = 0;
                end
        end
    end
    performanceOutput.rightPhaseAccuracy = sum(right_correct_gait_phase)/length(right_correct_gait_phase);
    performanceOutput.rightPhaseStates = right_correct_gait_phase;
end

end