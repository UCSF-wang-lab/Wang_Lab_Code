function updateForcePlot(src,event)

main_window = findobj('Tag','main_window');

% Clear old plot
yyaxis(main_window.UserData.plot_axes(1),'left');
cla(main_window.UserData.plot_axes(1));


if contains(src.Tag,'FP1')
    FP1 = ['FP Right - ',src.String{src.Value}];
    % Get baseline values
    if strcmp(src.String{src.Value},'X')
        baseline_FP1 = mean(main_window.UserData.Vicon_data.Fx_N_FP1(1:10000));
        FP1_data = main_window.UserData.Vicon_data.Fx_N_FP1;
    elseif strcmp(src.String{src.Value},'Y')
        baseline_FP1 = mean(main_window.UserData.Vicon_data.Fy_N_FP1(1:10000));
        FP1_data = main_window.UserData.Vicon_data.Fy_N_FP1;
    else
        baseline_FP1 = mean(main_window.UserData.Vicon_data.Fz_N_FP1(1:10000));
        FP1_data = main_window.UserData.Vicon_data.Fz_N_FP1;
    end
    
    FP_other = findobj('Tag','Plot_Line_FP2');
    FP2 = ['FP Left - ',FP_other.String{FP_other.Value}];
    if strcmp(FP_other.String{FP_other.Value},'X')
        baseline_FP2 = mean(main_window.UserData.Vicon_data.Fx_N_FP2(1:10000));
        FP2_data = main_window.UserData.Vicon_data.Fx_N_FP2;
    elseif strcmp(FP_other.String{FP_other.Value},'Y')
        baseline_FP2 = mean(main_window.UserData.Vicon_data.Fy_N_FP2(1:10000));
        FP2_data = main_window.UserData.Vicon_data.Fy_N_FP2;
    else
        baseline_FP2 = mean(main_window.UserData.Vicon_data.Fz_N_FP2(1:10000));
        FP2_data = main_window.UserData.Vicon_data.Fz_N_FP2;
    end
else
    FP2 = ['FP Left - ',src.String{src.Value}];
    % Get baseline values
    if strcmp(src.String{src.Value},'X')
        baseline_FP2 = mean(main_window.UserData.Vicon_data.Fx_N_FP2(1:10000));
        FP2_data = main_window.UserData.Vicon_data.Fx_N_FP2;
    elseif strcmp(src.String{src.Value},'Y')
        baseline_FP2 = mean(main_window.UserData.Vicon_data.Fy_N_FP2(1:10000));
        FP2_data = main_window.UserData.Vicon_data.Fy_N_FP2;
    else
        baseline_FP2 = mean(main_window.UserData.Vicon_data.Fz_N_FP2(1:10000));
        FP2_data = main_window.UserData.Vicon_data.Fz_N_FP2;
    end
    
    FP_other = findobj('Tag','Plot_Line_FP1');
    FP1 = ['FP Right - ',FP_other.String{FP_other.Value}];
    if strcmp(FP_other.String{FP_other.Value},'X')
        baseline_FP1 = mean(main_window.UserData.Vicon_data.Fx_N_FP1(1:10000));
        FP1_data = main_window.UserData.Vicon_data.Fx_N_FP1;
    elseif strcmp(FP_other.String{FP_other.Value},'Y')
        baseline_FP1 = mean(main_window.UserData.Vicon_data.Fy_N_FP1(1:10000));
        FP1_data = main_window.UserData.Vicon_data.Fy_N_FP1;
    else
        baseline_FP1 = mean(main_window.UserData.Vicon_data.Fz_N_FP1(1:10000));
        FP1_data = main_window.UserData.Vicon_data.Fz_N_FP1;
    end
end

% Plot the data
leg_elem = [];
vicon_time = main_window.UserData.Vicon_time;
colors = colororder;
hold(main_window.UserData.plot_axes(1),'on');

leg_elem(1) = plot(main_window.UserData.plot_axes(1),...
    vicon_time,...
    FP1_data-baseline_FP1,'Color',colors(1,:),...
    'LineStyle','-','DisplayName',FP1);

leg_elem(2) = plot(main_window.UserData.plot_axes(1),...
    vicon_time,...
    FP2_data-baseline_FP2,'Color',colors(2,:),...
    'LineStyle','-','DisplayName',FP2);

if ~isempty(main_window.UserData.Delsys_data)
    TTL_ind = find(contains(main_window.UserData.Delsys_data.Chan_names,'TTL'));
    chan_name = main_window.UserData.Delsys_data.Chan_names{TTL_ind};
    
    yyaxis(main_window.UserData.plot_axes(1),'right');
    cla(main_window.UserData.plot_axes(1));
    
    leg_elem(3) = plot(main_window.UserData.plot_axes(1),...
        main_window.UserData.Delsys_data.Time.(chan_name),...
        main_window.UserData.Delsys_data.Data.(chan_name),...
        'Color',colors(4,:),'DisplayName','Delsys TTL Pulse');
end

xlabel(main_window.UserData.plot_axes(1),'Time (s)','FontWeight','bold');
axes(main_window.UserData.plot_axes(1));

% Plot old markers
if sum(isnan(main_window.UserData.Marked_events)) ~= length(main_window.UserData.Marked_events)
    for i = 1:length(main_window.UserData.Marked_events)
        if ~isnan(main_window.UserData.Marked_events(i))
            if i == 1
                h = xline(main_window.UserData.plot_axes(1),main_window.UserData.Marked_events(1),'--k','DisplayName','Cue');
                main_window.UserData.plot_markers(1) = h;
                leg_elem(end+1) = h;
            elseif i == 2
                h = xline(main_window.UserData.plot_axes(1),main_window.UserData.Marked_events(2),'--g','DisplayName','Movement onset');
                main_window.UserData.plot_markers(2) = h;
                leg_elem(end+1) = h;
            else
                h = xline(main_window.UserData.plot_axes(1),main_window.UserData.Marked_events(3),'--r','DisplayName','Stance leg toe-off');
                main_window.UserData.plot_markers(3) = h;
                leg_elem(end+1) = h;
            end
        end
    end
end

legend(leg_elem);
end