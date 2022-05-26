function outTable = parseRawRCSFolder(folderPath)
%% parseRawRCSFolder
% Parses through all session folders of RCS recordings. This function will
% create a table of all session folders with their approximate start time
% and any "extra_comments" in the EventLog.json file.
%
% INPUT:
%       folderPath  [=] The folder of raw RCS json files to parse through.
%
% OUTPUT:
%       outTable    [=] Table of the session files, approximate start time
%                       of the session, and any extra comments by the user.
%
% Example call:
%       f = '/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_02/v5_2022-05-04_dbsOpt/Data/Raw Data/RC+S/Left INS';
%       A = parseRawRCSFolder(f);
%
% Author: Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Date:   05/26/2022

% Check if a folder path was passed in. If not, will prompt the user for a
% directory to look through.
if ~exist('folderPath','var') || isempty(folderPath)
    folderPath = uigetdir();
end

% Grab list of session folders in the directory
folderList = dir(folderPath);
sessionFolders = {folderList(contains({folderList.name},'Session')).name};

% Extract RCS serial number
temp = dir(fullfile(folderPath,sessionFolders{1}));
deviceName = temp(contains({temp.name},'Device')).name;


% Go through all found session folders and extract comments. 
extraCommentsArray = {};
sessionArray = {};
sessionTime = {};
for i = 1:length(sessionFolders)
    fid = fopen(fullfile(folderPath,sessionFolders{i},deviceName,'EventLog.json'));
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);

    % Add entry for each session
    extraCommentsArray = [extraCommentsArray;{'New Session'}];
    sessionArray = [sessionArray;repelem({sessionFolders{i}},1,1)];
    sessionTime = [sessionTime;repelem(convertSessionToDateTime(sessionFolders{i}),1,1)];

    % Go through all extra comments
    for j = 1:length(val)
        if strcmp(val(j).Event.EventType,'extra_comments')
            comments = strsplit(val(j).Event.EventSubType,newline);
            extraCommentsArray = [extraCommentsArray;comments'];
            sessionArray = [sessionArray;repelem({sessionFolders{i}},length(comments),1)];
            sessionTime = [sessionTime;repelem(convertSessionToDateTime(sessionFolders{i}),length(comments),1)];
        end


    end
end
deviceArray = repelem({deviceName},length(extraCommentsArray),1);

% Create the output table
outTable = table(deviceArray,sessionArray,sessionTime,extraCommentsArray,'VariableNames',{'Device Name','Session','Session Approx. Start Time','Comments'});

end

%% HELPER FUNCTIONS
function convertedDateTime = convertSessionToDateTime(sessionString)
% convertSessionToDateTime
% Helper function to convert the session time in UTC to local time (i.e.
% pacific time).

dt = datetime(str2num(sessionString(8:end)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
dtSplit = strsplit(string(dt),' ');

if str2num(dtSplit{2})-7 < 0
    timeCorrection = abs(str2num(dtSplit{2}(1:2))-7);
    dtSplit{2}(1:2) = num2str(24-timeCorrection);
else
    dtSplit{2}(1:2) = num2str(str2num(dtSplit{2}(1:2))-7);
end

convertedDateTime = dtSplit(1) + " " + dtSplit(2);
end