function plotAlignmentTimes()
main = findobj('Tag','AlignData_main_window');
rcs_time_menu = findobj('Tag','RCSTimeUnit');
for i = 1:length(main.UserData.plot_axes)
    plot_title = main.UserData.plot_axes(i).Title.String;
    if contains(plot_title,'Left')  && contains(plot_title,'LFP')
        if ~isnan(main.UserData.alignment_times(1))
            deletePreviousMarks(main.UserData.plot_axes(i));
            if (floor(main.UserData.alignment_times(1))~=main.UserData.alignment_times(1)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                xline(main.UserData.plot_axes(i),main.UserData.alignment_times(1),'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(1))~=main.UserData.alignment_times(1)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                temp1 = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = seconds(temp1-temp1(1));
                [~,ind] = min(abs(temp2-main.UserData.alignment_times(1)));
                mark_time = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(ind),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(1))==main.UserData.alignment_times(1)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                time_mark = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(main.UserData.alignment_times(1)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),time_mark,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(1))==main.UserData.alignment_times(1)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                temp1 = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(main.UserData.alignment_times(1)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                mark_time = seconds(temp1-temp2);
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            end
        end
    elseif contains(plot_title,'Right')  && contains(plot_title,'LFP')
        if ~isnan(main.UserData.alignment_times(2))
            deletePreviousMarks(main.UserData.plot_axes(i));
            if (floor(main.UserData.alignment_times(2))~=main.UserData.alignment_times(2)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                xline(main.UserData.plot_axes(i),main.UserData.alignment_times(2),'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(2))~=main.UserData.alignment_times(2)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                temp1 = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = seconds(temp1-temp1(1));
                [~,ind] = min(abs(temp2-main.UserData.alignment_times(2)));
                mark_time = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(ind),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(2))==main.UserData.alignment_times(2)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                time_mark = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(main.UserData.alignment_times(2)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),time_mark,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(2))==main.UserData.alignment_times(2)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                temp1 = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(main.UserData.alignment_times(2)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                mark_time = seconds(temp1-temp2);
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            end
        end
    elseif contains(plot_title,'Left')  && contains(plot_title,'Accel') && ~contains(plot_title,'Linear')
        if ~isnan(main.UserData.alignment_times(3))
            deletePreviousMarks(main.UserData.plot_axes(i));
            if (floor(main.UserData.alignment_times(3))~=main.UserData.alignment_times(3)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                xline(main.UserData.plot_axes(i),main.UserData.alignment_times(3),'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(3))~=main.UserData.alignment_times(3)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                temp1 = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = seconds(temp1-temp1(1));
                [~,ind] = min(abs(temp2-main.UserData.alignment_times(3)));
                mark_time = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(ind),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(3))==main.UserData.alignment_times(3)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                time_mark = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(main.UserData.alignment_times(3)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),time_mark,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(3))==main.UserData.alignment_times(3)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                temp1 = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(main.UserData.alignment_times(3)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.Accel_data.Left.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                mark_time = seconds(temp1-temp2);
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            end
        end
    elseif contains(plot_title,'Right')  && contains(plot_title,'Accel') && ~contains(plot_title,'Linear')
        if ~isnan(main.UserData.alignment_times(4))
            deletePreviousMarks(main.UserData.plot_axes(i));
            if (floor(main.UserData.alignment_times(4))~=main.UserData.alignment_times(4)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                xline(main.UserData.plot_axes(i),main.UserData.alignment_times(4),'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(4))~=main.UserData.alignment_times(4)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                temp1 = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = seconds(temp1-temp1(1));
                [~,ind] = min(abs(temp2-main.UserData.alignment_times(4)));
                mark_time = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(ind),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(4))==main.UserData.alignment_times(4)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'Date Time')
                time_mark = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(main.UserData.alignment_times(4)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                xline(main.UserData.plot_axes(i),time_mark,'--r','Tag','alignment_marker');
            elseif (floor(main.UserData.alignment_times(4))==main.UserData.alignment_times(4)) && strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')
                temp1 = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(main.UserData.alignment_times(4)),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                temp2 = datetime(main.UserData.Accel_data.Right.accelDataTable.DerivedTime(1),'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
                mark_time = seconds(temp1-temp2);
                xline(main.UserData.plot_axes(i),mark_time,'--r','Tag','alignment_marker');
            end
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
    elseif contains(plot_title,'Left')  && contains(plot_title,'Rover')
        if ~isnan(main.UserData.alignment_times(8))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.Rover_data.Left.DateTime(main.UserData.alignment_times(8)),'--r','Tag','alignment_marker');
        end
    elseif contains(plot_title,'Right')  && contains(plot_title,'Rover')
        if ~isnan(main.UserData.alignment_times(9))
            deletePreviousMarks(main.UserData.plot_axes(i));
            xline(main.UserData.plot_axes(i),main.UserData.Rover_data.Right.DateTime(main.UserData.alignment_times(9)),'--r','Tag','alignment_marker');
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