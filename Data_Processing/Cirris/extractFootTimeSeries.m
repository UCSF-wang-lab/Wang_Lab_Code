function extractFootTimeSeries()
% Extracts only the (position and speed) time series information from the 
% massive csv file output of the Cirris app
%
% INPUT:   filename    [=] The csv output of Cirris
%
% Author:   Eleni Patelaki
% Date:     10/25/23

[filename,filepath] = uigetfile('*.csv');
fid = fopen(fullfile(filepath,filename));
timeframeSection = 0;
rfootflag = 0;
lfootflag = 0;
frameIdx = [];
frameTime = [];
frameIdx_L = [];
posX_L = [];
posY_L = [];
posZ_L = [];
speedX_L = [];
speedY_L = [];
speedZ_L = [];
frameIdx_R = [];
posX_R = [];
posY_R = [];
posZ_R = [];
speedX_R = [];
speedY_R = [];
speedZ_R = [];

while(~feof(fid))
    line = fgetl(fid);
    if contains(line,'Frame Index') %&& isempty(frameIdx)
        timeframeSection = 1;
    end
    if contains(line,'PDGait_Left_foot') %&& isempty(frameIdx_L)
        lfootflag = 1;
    end
    if contains(line,'PDGait_Right_foot') %&& isempty(frameIdx_R)
        rfootflag = 1;
    end
    if timeframeSection && ~lfootflag && ~rfootflag
        splits = strsplit(line,',');
        if ~isnan(str2double(splits{1}))
            if mod(splits{1},1) == 0 
                frameIdx(end+1) = str2double(splits{1});
                frameTime(end+1) = str2double(splits{2});
            end
        end
    end
    
    if timeframeSection && lfootflag
        splits = strsplit(line,',');
        if ~isnan(str2double(splits{1}))
            if (mod(splits{1},1)==0) & (numel(frameIdx_L)<numel(frameIdx))
                frameIdx_L(end+1) = str2double(splits{1});
                posX_L(end+1) = str2double(splits{2});
                posY_L(end+1) = str2double(splits{3});
                posZ_L(end+1) = str2double(splits{4});
                speedX_L(end+1) = str2double(splits{5});
                speedY_L(end+1) = str2double(splits{6});
                speedZ_L(end+1) = str2double(splits{7});
            end
        end
    end
    
    if timeframeSection && rfootflag
        splits = strsplit(line,',');
        if ~isnan(str2double(splits{1}))
            if (mod(splits{1},1)==0) & (numel(frameIdx_R)<numel(frameIdx))
                frameIdx_R(end+1) = str2double(splits{1});
                posX_R(end+1) = str2double(splits{2});
                posY_R(end+1) = str2double(splits{3});
                posZ_R(end+1) = str2double(splits{4});
                speedX_R(end+1) = str2double(splits{5});
                speedY_R(end+1) = str2double(splits{6});
                speedZ_R(end+1) = str2double(splits{7});
            end
        end
    end
end


assert(isequal(frameIdx,frameIdx_R)&&isequal(frameIdx,frameIdx_L))

timeL=frameTime;
timeR=frameTime;

lfootPos = table(timeL',posX_L',posY_L',posZ_L',speedX_L',speedY_L',speedZ_L','VariableNames',{'Time','PosX','PosY','PosZ','SpeedX','SpeedY','SpeedZ'});
rfootPos = table(timeR',posX_R',posY_R',posZ_R',speedX_R',speedY_R',speedZ_R','VariableNames',{'Time','PosX','PosY','PosZ','SpeedX','SpeedY','SpeedZ'});

if ~exist('savepath','var')
    savepath = uigetdir();
    [~,base_fname,~] = fileparts(filename);
    csv_rfoot = fullfile(savepath,strcat(base_fname,'_RFoot.csv'));
    csv_lfoot = fullfile(savepath,strcat(base_fname,'_LFoot.csv'));
end


try
    writetable(rfootPos,csv_rfoot);
    writetable(lfootPos,csv_lfoot);
catch
    error('Unable to write table to file');
end


end