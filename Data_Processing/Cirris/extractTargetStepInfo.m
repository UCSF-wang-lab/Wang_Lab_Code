function extractTargetStepInfo(filename,savepath,align_timediff)
% Extracts only the target information from the massive csv file output of
% the Cirris app
%
% Input:    filename    [=] The csv output of Cirris
%           savepath    [=] Save name and path of the extract target
%                           information
%       align_timediff  [=] Time difference to account for, between the
%                           Cirris and the Xsens stream (Cirris minus Xsens diff)
% Output:       NONE
% Example call: extractTargetStepInfo([],[],407.377)
%
% Author:      Kenneth Louie
% Contributor: Eleni Patelaki
% Date:     10/12/20

if isempty(filename)
    [filename,filepath] = uigetfile('*.csv');
else
    if ispc
        fileparts = strsplit(filename,'\');
    elseif ismac
        fileparts = strsplit(filename,'/');
    else
        error('Platform not currently supported.');
    end
    filepath  = fullfile(fileparts{1:end-1});
    filename =  fullfile(fileparts{end});
end

fid = fopen(fullfile(filepath,filename));

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
            TaskTimer(end+1) = str2double(splits{5})-align_timediff;
        end
    end
end

outTable = table(TargetNum(:),Side(:),StepModifier(:),Success(:),TaskTimer(:),...
    'VariableNames',{'Target Number','Side','Step Modifier','Success','Task Timer'});

if ~exist('savepath','var')||isempty(savepath)
    % [~,fname,~] = fileparts(filename);
    fname = extractBefore(filename,'.');
    savepath = fullfile(filepath,strcat(fname,'_all_targets_aligned.csv'));
end


try
    writetable(outTable,savepath);
catch
    error('Unable to write table to file');
end


end