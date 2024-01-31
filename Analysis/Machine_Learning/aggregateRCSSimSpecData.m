function aggregateTables = aggregateRCSSimSpecData(patientFolder)
patientFolder = '/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_01';
folderList = dir([patientFolder,'/*v*']);
subjectID = patientFolder(find(patientFolder=='/',1,'last')+1:end);
aggregateTables = [];

% Go through all of the patient's folders
for i = 1:length(folderList)
    % Check to make sure the Rover folder is not included
    if ~strcmp(folderList(i).name,'Rover')
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('Parsing folder: %s\n',folderList(i).name);
        
        % Grab all the spectrogram and mat files
        spec_files = dir(fullfile(folderList(i).folder,folderList(i).name,'Data','Analysis Data','aDBS','*bitshift*'));
        mat_files = dir(fullfile(folderList(i).folder,folderList(i).name,'Data','Aligned Data','*w_Gait_Events*'));
        
        if ~isempty(spec_files)
            % Trim number of mat files to the most recent ones
            mat_files = trimMatFiles(mat_files);
            
            % Loop through mat files and associated spec files and fill in a
            % table
            for j = 1:length(mat_files)
                % Load in mat file
                load(fullfile(folderList(1).folder,folderList(i).name,'Data/Aligned Data',mat_files{j}));
                
                % Check which spec files are associated with this file
                specPatientState = sprintf('%s_stim_%s_med',aligned_data.stim_condition,aligned_data.med_condition);
                specFileMatch = cellfun(@(x)contains(x,specPatientState,'IgnoreCase',true),{spec_files.name});
                specFileNames = {spec_files(specFileMatch).name};
                
                % Used to keep track of spec params. Reduces additional
                % computations
                prev_fs = 0;
                prev_nfft = 0;
                prev_fft_int = 0;
                prev_hemisphere = 0;
                
                % Create table of all spec files that are related to each
                % other
                for k = 1:length(specFileNames)
                    fprintf('Reading file: %s\n',specFileNames{k});
                    
                    % Load in spectral data
                    specData = readtable(fullfile(folderList(1).folder,folderList(i).name,'Data/Analysis Data/aDBS',specFileNames{k}));
                    
                    % Grab current spec parameters and hemisphere
                    specParam = strsplit(specFileNames{k},'_');
                    fs = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'fs'),specParam))+1)));
                    nfft = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'nfft'),specParam))+1)));
                    fft_int = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'fftint'),specParam))+1)));
                    overlap = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'overlap'),specParam))+1)));
                    bitshift = round(str2double(specParam(find(cellfun(@(x) strcmp(x,'fftbitshift'),specParam))+1)));
                    
                    if contains(specFileNames{k},'left','IgnoreCase',true)
                        hemisphere = "left";
                    else
                        hemisphere = "right";
                    end
                    
                    % Check to see if this calculation was done previously
                    % but for a different bitshift or hemisphere. If it
                    % wasn't, calculate it for the new fft settings.
%                     sameSpecParams = (fs ~= prev_fs && nfft ~= prev_nfft && fft_int ~= prev_fft_int && strcmp(hemisphere,prev_hemisphere));
                    sameSpecParams = 0;
                    if ~sameSpecParams
                        state = zeros(height(specData),1);
                        
                        % Update values
                        prev_fs = fs;
                        prev_nfft = nfft;
                        prev_fft_int = fft_int;
                        prev_hemisphere = hemisphere;
                        
                        % New parameters, calculate the correct states base
                        % on gait events.
                        sortedGaitEvents = sortGaitEvents(aligned_data.gait_events,'RTO');
                        for m = 1:height(specData)
                            % Specific for swing phase. Future will change
                            % to make this generic.
                            if strcmp(hemisphere,'left')
                                toeOff = (specData.time(m)>=sortedGaitEvents.RTO) & (specData.time(m)<=sortedGaitEvents.RHS);
                                if sum(toeOff) == 1
                                    state(m) = 1;
                                end
                            else
                                toeOff = (specData.time(m)>=sortedGaitEvents.LTO) & (specData.time(m)<=sortedGaitEvents.LHS);
                                if sum(toeOff) == 1
                                    state(m) = 1;
                                end
                            end
                        end
                    end
                    
                    % Create table with states and values. Always
                    % recalculate as verification on the RCS recording
                    % contacts are the same. 
                    % includeColumns = cellfun(@isempty,cellfun(@(x)regexp(x,'time|_[0-3]_|_[1-2][0-9][0-9]_'),specData.Properties.VariableNames,'UniformOutput',false));
                    includeColumns = ~cellfun(@isempty,cellfun(@(x)regexp(x,'_[4-9]_|_[0-5][0-9]_'),specData.Properties.VariableNames,'UniformOutput',false));
                    currTable = specData(:,includeColumns);
                    
                    % Other columns
                    subjectIDVec = repmat({subjectID},height(currTable),1);
                    
                    switch subjectID
                        case 'gait_RCS_01'
                            if contains(folderList(i).name,'v4')
                                visitNameVec = repmat({'dbsOptBilateral'},height(currTable),1);
                            end
                        case 'gait_RCS_02'
                            if contains(folderList(i).name,'v4')
                                visitNameVec = repmat({'preprogrammingUnilateral'},height(currTable),1);
                            elseif contains(folderList(i).name,'v7')
                                visitNameVec = repmat({'preprogrammingBilateral'},height(currTable),1);
                            elseif contains(folderList(i).name,'v5')
                                visitNameVec = repmat({'dbsOptUnilateral'},height(currTable),1);
                            elseif contains(folderList(i).name,'v8')
                                visitNameVec = repmat({'dbsOptBilateral'},height(currTable),1);
                            end
                        case 'gait_RCS_04'
                            if contains(folderList(i).name,'v3')
                                visitNameVec = repmat({'preprogrammingBilateral'},height(currTable),1);
                            elseif contains(folderList(i).name,'v5')
                                visitNameVec = repmat({'dbsOptBilateral'},height(currTable),1);
                            end
                        case 'gait_RCS_05'
                            if contains(folderList(i).name,'v3')
                                visitNameVec = repmat({'preprogrammingUnilateral'},height(currTable),1);
                            elseif contains(folderList(i).name,'v4')
                                visitNameVec = repmat({'dbsOptUnilateral'},height(currTable),1);
                            end
                        otherwise
                    end
                    
                    medStateVec = repmat({aligned_data.med_condition},height(currTable),1);
                    stimStateVec = repmat({aligned_data.stim_condition},height(currTable),1);
                    hemisphereVec = repmat({hemisphere},height(currTable),1);
                    fsVec = repmat(fs,height(currTable),1);
                    nfftVec = repmat(nfft,height(currTable),1);
                    overlapVec = repmat(overlap,height(currTable),1);
                    fftIntVec = repmat(fft_int,height(currTable),1);
                    bitshiftVec = repmat(bitshift,height(currTable),1);
                    extraTable = table(subjectIDVec,visitNameVec,medStateVec,stimStateVec,hemisphereVec,fsVec,nfftVec,overlapVec,fftIntVec,bitshiftVec,state,...
                        'VariableNames',{'SubjectID','VisitName','MedState','StimState','Side','FS','NFFT','OverlapPercent','FFTInterval','Bitshift','Class'});
                    outTable = [extraTable,currTable];
                    
                    aggregateTableName = sprintf('%s_AggregateSpecData_fs%d_nfft%d',subjectID,fs,nfft);
                    
                    if isempty(aggregateTables) || sum(cellfun(@(x) strcmp(x,aggregateTableName),fields(aggregateTables))) == 0
                        aggregateTables.(aggregateTableName) = outTable;
                    elseif sum(cellfun(@(x) strcmp(x,aggregateTableName),fields(aggregateTables))) == 1
                        aggregateTables.(aggregateTableName) = [aggregateTables.(aggregateTableName);outTable];
                    end
                end
            end
        end
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    end
end
end

function trimmedFiles = trimMatFiles(mat_files)
trimmedFiles = [];

off_stim_off_med_files = matchFiles(mat_files,'off_stim_off_meds');
off_stim_on_med_files = matchFiles(mat_files,'off_stim_on_meds');
on_stim_off_med_files = matchFiles(mat_files,'on_stim_off_meds');
on_stim_on_med_files = matchFiles(mat_files,'on_stim_on_meds');

trimmedFiles = [off_stim_off_med_files;off_stim_on_med_files;on_stim_off_med_files';on_stim_on_med_files'];
end

function matchedFiles = matchFiles(file_names,name2match)
matchedFiles = [];

fileMatch = cellfun(@(x)contains(x,name2match,'IgnoreCase',true),{file_names.name});
fileMatchInd = find(fileMatch);

% Check if there are multiple trials for a single med and stim state
if ~isempty(fileMatchInd)
    fileNames = {file_names(fileMatchInd).name};
    multipleTrials = sum(cellfun(@(x) contains(x,'Trial2','IgnoreCase',true),fileNames));
    
    if multipleTrials>0
        for i = 1:2
            % Find all files for Trial#i
            trialMatch = cellfun(@(x) contains(x,sprintf('Trial%d',i),'IgnoreCase',true),fileNames);
            trialFileInd = find(trialMatch);
%             trialFileNames = fileNames(trialFileInd);
            
            % Determine most recent file
            fileDates = cellfun(@datetime,{file_names(fileMatchInd(trialFileInd)).date});
            [~,mostRecentFileInd] = max(fileDates);
            matchedFiles{end+1} = file_names(trialFileInd(mostRecentFileInd)).name;
        end
    else
        % Determine most recent file
        fileDates = cellfun(@datetime,{file_names(fileMatchInd).date});
        [~,mostRecentFileInd] = max(fileDates);
        matchedFiles = {file_names(fileMatchInd(mostRecentFileInd)).name};
    end
end
end