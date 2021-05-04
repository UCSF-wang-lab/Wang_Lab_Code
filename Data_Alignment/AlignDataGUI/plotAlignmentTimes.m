function plotAlignmentTimes()
main = findobj('Tag','AlignData_main_window');
for i = 1:length(main.UserData.plot_axes)
    plot_title = main.UserData.plot_axes(i).Title.String;
    if contains(plot_title,'Left')  && contains(plot_title,'LFP')
        if ~isnan(main.UserData.alignment_times(1))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(1),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Right')  && contains(plot_title,'LFP')
        if ~isnan(main.UserData.alignment_times(2))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(2),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Left')  && contains(plot_title,'Accel')
        if ~isnan(main.UserData.alignment_times(3))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(3),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Right')  && contains(plot_title,'Accel')
        if ~isnan(main.UserData.alignment_times(4))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(4),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Delsys')
        if ~isnan(main.UserData.alignment_times(5))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(5),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Xsens')
        if ~isnan(main.UserData.alignment_times(6))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(6),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Force')
        if ~isnan(main.UserData.alignment_times(7))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(7),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Teensey')
        if ~isnan(main.UserData.alignment_times(8))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.alignment_times(8),'--r','Tag','alignment_marker');
        end
    end
end
end

function deletePreviousMarks(axes)
for i = 1:length(axes.Children)
    if strcmp(axes.Children(i).Tag,'alignment_marker')
        delete(axes.Children(i));
        return;
    end
end
end