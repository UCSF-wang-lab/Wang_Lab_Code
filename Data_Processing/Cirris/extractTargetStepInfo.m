function extractTargetStepInfo(filename,savepath)
% Extracts only the target information from the massive csv file output of
% the Cirris app
%
% Input:    filename    [=] The csv output of Cirris
%           savepath    [=] Save name and path of the extract target
%                           information
% Output:   NONE
%
% Author:   Kenneth Louie
% Date:     10/12/20

fid = fopen(filename);
targetSection = 0;
TargetNum = [];
Side = {};
StepModifier = [];
Success = [];
TaskTimer = [];
while(~feof(fid))
    line = fgetl(fid);
    if contains(line,'Start task')
        targetSection = 1;
%         fprintf(line + "\n");
    end
    
    if targetSection
        splits = strsplit(line,',');
        if (contains(line,'Left') || contains(line,'Right')) && ...
                ~isnan(str2double(splits{1}(2:end)))
%             fprintf(line + "\n");
            TargetNum(end+1) = str2double(splits{1}(2:end));
            Side{end+1} = splits{2};
            StepModifier(end+1) = str2double(splits{3});
            Success(end+1) = str2double(splits{4}(1:end-1));
            TaskTimer(end+1) = str2double(splits{5});
        end
    end
end

outTable = table(TargetNum(:),Side(:),StepModifier(:),Success(:),TaskTimer(:),...
    'VariableNames',{'Target Number','Side','Step Modifier','Success','Task Timer'});

if ~exist('savepath','var')
    [A,B,~] = fileparts(filename);
    savepath = fullfile(A,B)+"_all_targets.csv";
end


try
    writetable(outTable,savepath);
catch
    error('Unable to write table to file');
end


end