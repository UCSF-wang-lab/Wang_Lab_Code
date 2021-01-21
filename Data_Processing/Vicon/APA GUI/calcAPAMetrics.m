function calcAPAMetrics(src,obj)
%
%
% Author:   Kenneth Louie
% Date:     12/29/20

CoP_startInd = src.Parent.UserData.CoP_startInd;
CoP_endInd = src.Parent.UserData.CoP_endInd;

L1_ind = src.Parent.UserData.Marked_events(4);
L2_ind = src.Parent.UserData.Marked_events(5);

xpos_start = src.Parent.UserData.Vicon_data.Cx_mm(CoP_startInd);
ypos_start = src.Parent.UserData.Vicon_data.Cy_mm(CoP_endInd);

xpos_L1 = src.Parent.UserData.Vicon_data.Cx_mm(CoP_startInd+L1_ind);
ypos_L1 = src.Parent.UserData.Vicon_data.Cy_mm(CoP_startInd+L1_ind);

xpos_L2 = src.Parent.UserData.Vicon_data.Cx_mm(CoP_startInd+L2_ind);
ypos_L2 = src.Parent.UserData.Vicon_data.Cy_mm(CoP_startInd+L2_ind);

xpos_end = src.Parent.UserData.Vicon_data.Cx_mm(CoP_endInd);
ypos_end = src.Parent.UserData.Vicon_data.Cy_mm(CoP_endInd);

S1_x_displacement = abs(xpos_start-xpos_L1);
S1_y_displacement = abs(ypos_start-ypos_L1);

S2_x_displacement = abs(xpos_L1-xpos_L2);
S2_y_displacement = abs(ypos_L1-ypos_L2);

S3_x_displacement = abs(xpos_L2-xpos_end);
S3_y_displacement = abs(ypos_L2-ypos_end);

src.Parent.UserData.APAMetrics = [S1_x_displacement,S2_x_displacement,S3_x_displacement;...
    S1_y_displacement,S2_y_displacement,S3_y_displacement];

src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;...
    'APA metrics calculated.'];
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
end