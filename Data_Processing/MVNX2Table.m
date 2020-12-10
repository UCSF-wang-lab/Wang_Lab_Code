function xsens_table = MVNX2Table(filename,savename)
%% Processes the converted mvnx files from Xsens to .mat and .csv
% The .mat file follows the same structure as previous .mat files.
% The .csv file outputs only the segment and joint data with no meta-data
% contained such as sampling rate.
%
% Dependencies:     load_mvnx - This file is from Xsens
%
% Inputs:   filename    [=] The file to be loaded in (.mvnx)
%           savename    [=] The name to save the .mat and .csv file. This
%                           must include the correct directory to save the
%                           data. If not passed in or an empty vector is
%                           passed, automatically save in the same folder
%                           as the file.
% Outputs:  xsens_table [=] Extracted xsens values
% 
% Author:   Kenneth Louie
% Date:     10/12/20
%% to mat
% load data
tree = load_mvnx(filename);

% read some basic data from the file
mvnxVersion = tree;
fileComments = tree.subject.comment;

%read some basic properties of the subject;
frameRate = tree.subject.frameRate;
suitLabel = tree.subject.label;
originalFilename = tree.subject.originalFilename;
recDate = tree.subject.recDate;
segmentCount = tree.subject.segmentCount;

%retrieve sensor labels
%creates a struct with sensor data
if isfield(tree.subject,'sensors') && isstruct(tree.subject.sensors)
    sensorData = tree.subject.sensors.sensor;
end

%retrieve segment labels
%creates a struct with segment definitions
if isfield(tree.subject,'segments') && isstruct(tree.subject.segments)
    segmentData = tree.subject.segments.segment;
end

%retrieve the data frames from the subject
nSamples = length(tree.subject.frames.frame);

%pre allocate some memory for the position of Segment1
p_Segment1 = zeros(nSamples,3);
%read the data from the structure e.g. segment 1
if isfield(tree.subject.frames.frame(1),'position')
    for i=[1:nSamples]
        p_Segment1(i,:)=tree.subject.frames.frame(i).position(1:3);
    end
    
    %Plot (for debugging)
    figure('name','Position of first segment')
    plot(p_Segment1)
    figure('name','Position of first segment in 3D')
    plot3(p_Segment1(:,1),p_Segment1(:,2),p_Segment1(:,3));
end

[default_path,default_filename,default_ext] = fileparts(filename);
if exist('savename','var') && ~isempty(savename)
    save_name = fullfile(savename,[default_filename,'.mat']);
else
    save_name = fullfile(default_path,[default_filename,'.mat']);
end
save(save_name,'fileComments','frameRate','mvnxVersion','nSamples',...
    'originalFilename','p_Segment1','recDate','segmentCount',...
    'segmentData','sensorData','suitLabel','tree');


%% to CSV
% Only export the raw data, no extra meta-data
nFrames = length(unique([tree.subject.frames.frame.index]));
nSegments = tree.subject.segmentCount;
nJoints = length({tree.subject.joints.joint.label});
seg_pos = nan(nFrames,nSegments*3);          % XYZ sequence
seg_vel = nan(nFrames,nSegments*3);         % XYZ sequence
seg_accel = nan(nFrames,nSegments*3);       % XYZ sequence
seg_aVel = nan(nFrames,nSegments*3);        % XYZ sequence (angular vel)
seg_aAccel = nan(nFrames,nSegments*3);      % XYZ sequence (angular accel)
joint_angles = nan(nFrames,nJoints*3);      % ZXY sequence
time_vec = zeros(nFrames,1);                % Global time since same rate

currInd = 1;
for i = 1:length(tree.subject.frames.frame)
    if ~isempty(tree.subject.frames.frame(i).index)
        % Segment position matrix
        seg_pos(currInd,:) = tree.subject.frames.frame(i).position;
        
        % Segment veclocity matrix
        seg_vel(currInd,:) = tree.subject.frames.frame(i).velocity;
        
        % Segment acceleration matrix
        seg_accel(currInd,:) = tree.subject.frames.frame(i).acceleration;
        
        % Segment angular velocity matrix
        seg_aVel(currInd,:) = tree.subject.frames.frame(i).angularVelocity;
        
        % Segment angular acceleration matrix
        seg_aAccel(currInd,:) = tree.subject.frames.frame(i).angularAcceleration;
        
        % Joint angle matrix
        joint_angles(currInd,:) = tree.subject.frames.frame(i).jointAngle;
        
        currInd = currInd + 1;
    end
end

time_vec = 0:nFrames-1;
time_vec = time_vec'.*(1/tree.subject.frameRate);



% Create header cell array
header = cell(1,1+nSegments*3*5+nJoints*3);
header{1} = 'Time';
currInd = 2;
for i = 1:length({tree.subject.segments.segment.label})
    
    header{currInd} = [tree.subject.segments.segment(i).label,'_PosX'];
    header{currInd+nSegments*3*1} = [tree.subject.segments.segment(i).label,'_VelX'];
    header{currInd+nSegments*3*2} = [tree.subject.segments.segment(i).label,'_AccX'];
    header{currInd+nSegments*3*3} = [tree.subject.segments.segment(i).label,'_angVelX'];
    header{currInd+nSegments*3*4} = [tree.subject.segments.segment(i).label,'_angAccX'];
    
    
    header{currInd+1} = [tree.subject.segments.segment(i).label,'_PosY'];
    header{currInd+1+nSegments*3*1} = [tree.subject.segments.segment(i).label,'_VelY'];
    header{currInd+1+nSegments*3*2} = [tree.subject.segments.segment(i).label,'_AccY'];
    header{currInd+1+nSegments*3*3} = [tree.subject.segments.segment(i).label,'_angVelY'];
    header{currInd+1+nSegments*3*4} = [tree.subject.segments.segment(i).label,'_angAccY'];
    
    header{currInd+2} = [tree.subject.segments.segment(i).label,'_PosZ'];
    header{currInd+2+nSegments*3*1} = [tree.subject.segments.segment(i).label,'_VelZ'];
    header{currInd+2+nSegments*3*2} = [tree.subject.segments.segment(i).label,'_AccZ'];
    header{currInd+2+nSegments*3*3} = [tree.subject.segments.segment(i).label,'_angVelZ'];
    header{currInd+2+nSegments*3*4} = [tree.subject.segments.segment(i).label,'_angAccZ'];
    
    % Update count
    currInd = currInd + 3;
end

for j = 1:length({tree.subject.joints.joint.label})
    header{(j*3-2)+1+nSegments*3*5} = [tree.subject.joints.joint(j).label,'_Z'];
    header{(j*3-1)+1+nSegments*3*5} = [tree.subject.joints.joint(j).label,'_X'];
    header{(j*3)+1+nSegments*3*5} = [tree.subject.joints.joint(j).label,'_Y'];
end

xsens_table = array2table([time_vec,seg_pos,seg_vel,seg_accel,seg_aVel,seg_aAccel,joint_angles],'VariableNames',header);
writetable(xsens_table,[save_name(1:end-4),'.csv']);

end