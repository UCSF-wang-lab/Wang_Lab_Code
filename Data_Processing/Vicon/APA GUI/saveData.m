function saveData(src,obj)
%
%
% Author:   Kenneth Louie
% Date:     02/10/2021

if ~isfield(src.Parent.UserData,'APAMetrics')
    calcAPAMetrics();
end

cueTime = src.Parent.UserData.Marked_events(1);
mvmtOnsetTime = src.Parent.UserData.Marked_events(2);
stanceToeOffTime = src.Parent.UserData.Marked_events(3);
Vicon_time = src.Parent.UserData.Vicon_time;
CoP_startInd = src.Parent.UserData.CoP_startInd;
CoP_endInd = src.Parent.UserData.CoP_endInd;
L1_ind = src.Parent.UserData.Marked_events(4);
L2_ind = src.Parent.UserData.Marked_events(5);

S1_x_displacement = src.Parent.UserData.APAMetrics(1,1);
S1_y_displacement = src.Parent.UserData.APAMetrics(2,1);

S2_x_displacement = src.Parent.UserData.APAMetrics(1,2);
S2_y_displacement = src.Parent.UserData.APAMetrics(2,2);

S3_x_displacement = src.Parent.UserData.APAMetrics(1,3);
S3_y_displacement = src.Parent.UserData.APAMetrics(2,3);

src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;'Saving data...'];
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
pause(0.5);

% Save mat file
filename = [src.Parent.UserData.save_path,src.Parent.UserData.data_file(1:end-4),'_CoP_Markings'];
save([filename,'.mat'],'cueTime','mvmtOnsetTime','stanceToeOffTime',...
    'Vicon_time','CoP_startInd','CoP_endInd',...
    'L1_ind','L2_ind',...
    'S1_x_displacement','S1_y_displacement',...
    'S2_x_displacement','S2_y_displacement',...
    'S3_x_displacement','S3_y_displacement');

output_string = sprintf('Saved .mat file: %s',filename);
src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;output_string];
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
pause(0.5);

% Save txt file
fid = fopen([filename,'.txt'],'w+');
fprintf(fid,'Alignment file used: %s\n',src.Parent.UserData.data_file);
fprintf(fid,'Cue time: %3.6f\n',cueTime);
fprintf(fid,'Movement onset time: %3.6f\n',mvmtOnsetTime);
fprintf(fid,'Stance leg toe-off time: %3.6f\n',stanceToeOffTime);
fprintf(fid,'CoP Vicon Start Ind: %i\n',CoP_startInd);
fprintf(fid,'CoP Vicon End Ind: %i\n',CoP_endInd);
fprintf(fid,'Landmark #1 index: %i\n',L1_ind);
fprintf(fid,'Landmark #2 index: %i\n',L2_ind);
fprintf(fid,'S1 X Displacement: %3.6f\n',S1_x_displacement);
fprintf(fid,'S1 Y Displacement: %3.6f\n',S1_y_displacement);
fprintf(fid,'S2 X Displacement: %3.6f\n',S2_x_displacement);
fprintf(fid,'S2 Y Displacement: %3.6f\n',S2_y_displacement);
fprintf(fid,'S3 X Displacement: %3.6f\n',S3_x_displacement);
fprintf(fid,'S3 Y Displacement: %3.6f\n',S3_y_displacement);
fclose(fid);

output_string = sprintf('Saved .txt file: %s\nSaving finished.',filename);
src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;output_string];
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
end
