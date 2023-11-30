function performanceTable = aggregateAdaptivePerformance(base_folder,settings_num,gait_phase,flip_state)

base_search_string = [base_folder,'/*OG_*_Setting_%i*'];

accuracy = [];
side = [];
setting = [];
accuracy_type = [];

for j = 1:length(settings_num)
    search_string = sprintf(base_search_string,settings_num(j));
    file_list = dir(search_string);
    
    allStates_left = [];
    allStates_right = [];
    allGaitPhaseStates_left = [];
    allGaitPhaseStates_right = [];
    for i = 1:length(file_list)
        if contains(file_list(i).name,'Gait_Events')
            load(fullfile(file_list(i).folder,file_list(i).name));
            performanceOutput = calcAdaptivePerformance(aligned_data,gait_phase,flip_state);
            
            allStates_left = [allStates_left;performanceOutput.leftStates];
            allStates_right = [allStates_right;performanceOutput.rightStates];
            allGaitPhaseStates_left = [allGaitPhaseStates_left;performanceOutput.leftPhaseStates];
            allGaitPhaseStates_right = [allGaitPhaseStates_right;performanceOutput.rightPhaseStates];
        end
    end
    
    accuracy_left = sum(allStates_left)/length(allStates_left);
    accuracy_right = sum(allStates_right)/length(allStates_right);
    accuracyGaitPhase_left = sum(allGaitPhaseStates_left)/length(allGaitPhaseStates_left);
    accuracyGaitPhase_right = sum(allGaitPhaseStates_right)/length(allGaitPhaseStates_right);
    
    accuracy = [accuracy;accuracy_left;accuracyGaitPhase_left;accuracy_right;accuracyGaitPhase_right];
    side = [side;'L';'L';'R';'R'];
    setting = [setting;settings_num(j);settings_num(j);settings_num(j);settings_num(j)];
    accuracy_type = [accuracy_type;"Overall";"Phase";"Overall";"Phase"];
end

performanceTable = table(setting,side,accuracy_type,accuracy,'VariableNames',{'Setting','Side','AccuracyType','Accuracy'});

end