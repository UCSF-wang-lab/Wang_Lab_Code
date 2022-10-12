function outTable = processGIBinaryData(filePath,savePath)

fid = fopen(filePath);
data = fread(fid,[3,inf],'double');

trialType = {};
interTrialState = {};

for i = 1:size(data,2)
    trialType{i,1} = getTrialType(data(1,i));
    interTrialState{i,1} = getTrialState(data(2,i));
end

% outTable = array2table(data','VariableNames',{'GI Trial Type','Inter-trial State','Computer Time'});
outTable = table(trialType,interTrialState,data(3,:)','VariableNames',{'Trial Type','Inter-trial State','Event Time'});

[~,filename] = fileparts(filePath);
writetable(outTable,fullfile(savePath,[filename,'.csv']));
save(fullfile(savePath,[filename,'.mat']),"outTable");
end

function trialType = getTrialType(val)
    switch val
        case 1
            trialType = 'Cued';
        case 2
            trialType = 'Self-Initiated';
        case 3
            trialType = 'Prepare-and-Go';
    end
end

function trialState = getTrialState(val)
    switch val
        case 0
            trialState = 'Trial Start';
        case 1
            trialState = 'Stop';
        case 2
            trialState = 'Prepare';
        case 3
            trialState = 'Go';
        case 4
            trialState = 'Trial End';
    end
end