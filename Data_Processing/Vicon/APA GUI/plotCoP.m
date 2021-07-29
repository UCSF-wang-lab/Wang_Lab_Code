function plotCoP(src,obj)
% Author:   Kenneth Louie
% Date:     12/29/20

if isnan(src.Parent.UserData.Marked_events(2)) || isnan(src.Parent.UserData.Marked_events(3))
    src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;...
        'Not all makers set.'];
    src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    return;
end

vicon_time = (0:height(src.Parent.UserData.Vicon_data)-1)./1000;

% Find approximate index
[~,startInd] = min(abs(vicon_time-src.Parent.UserData.Marked_events(2)));
[~,endInd] = min(abs(vicon_time-src.Parent.UserData.Marked_events(3)));
src.Parent.UserData.CoP_startInd = startInd;
src.Parent.UserData.CoP_endInd = endInd;

if sum(contains(src.Parent.UserData.Vicon_data.Properties.VariableNames,'Fx_N')) == 1
        single_force_plate = true;
end

if single_force_plate
    cop_ml = src.Parent.UserData.Vicon_data.Cx_mm(startInd:endInd);
    cop_ap = src.Parent.UserData.Vicon_data.Cy_mm(startInd:endInd);
else
    cop_init_y = ((src.Parent.UserData.Vicon_data.Cx_mm_FP2(1)-src.Parent.UserData.Vicon_data.Cx_mm_FP1(1))/2) + src.Parent.UserData.Vicon_data.Cx_mm_FP1(1);
    cop_ap = ((src.Parent.UserData.Vicon_data.Cx_mm_FP2-src.Parent.UserData.Vicon_data.Cx_mm_FP1)/2) + src.Parent.UserData.Vicon_data.Cx_mm_FP1;
    cop_ap = cop_ap-cop_init;
    cop_ap_trim_ind = find(src.Parent.UserData.Vicon_data.Cx_mm_FP1>464,1,'first');
    
    cop_init_x = mean([src.Parent.UserData.Vicon_data.Cy_mm_FP1(1),src.Parent.UserData.Vicon_data.Cy_mm_FP2(1)]);
    cop_ml = mean([src.Parent.UserData.Vicon_data.Cy_mm_FP1,src.Parent.UserData.Vicon_data.Cy_mm_FP2],2);
    cop_ml = cop_ml-cop_init_x;
end

% Plot
cla(src.Parent.UserData.plot_axes(2));
hold(src.Parent.UserData.plot_axes(2),'on');
plot(src.Parent.UserData.plot_axes(2),cop_ml,cop_ap,'-k','linewidth',1.5);
scatter(src.Parent.UserData.plot_axes(2),cop_ml(1),cop_ap(1),60,'ok','filled');

src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;...
    'Plotted CoP data'];
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);

xlabel(src.Parent.UserData.plot_axes(2),'Mediolateral (mm)','FontWeight','bold');
ylabel(src.Parent.UserData.plot_axes(2),'Anteroposterior (mm)','FontWeight','bold');
end