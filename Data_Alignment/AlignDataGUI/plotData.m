function plotData(src,event)

% Determine plot window
parent_obj_title = src.Parent.Title;
if contains(parent_obj_title,'1')
    window = 1;
elseif contains(parent_obj_title,'2')
    window = 2;
elseif contains(parent_obj_title,'3')
    window = 3;
elseif contains(parent_obj_title,'4')
    window = 4;
end

% Determine window data source
window_data_source = src.Parent.Parent.Parent.UserData.plot_options.sources(window,1).String{src.Parent.Parent.Parent.UserData.plot_options.sources(window,1).Value};
if contains(window_data_source,'LFP')
    if contains(window_data_source,'Left')
        A = src.Parent.Parent.Parent.UserData.LFP_data.Left.timeDomainDataTable.DerivedTime;
        B = datetime(A,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        time = seconds(B-B(1));
        data = src.Parent.Parent.Parent.UserData.LFP_data.Left.timeDomainDataTable.(src.String{src.Value});
    elseif contains(window_data_source,'Right')
        A = src.Parent.Parent.Parent.UserData.LFP_data.Right.timeDomainDataTable.DerivedTime;
        B = datetime(A,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        time = seconds(B-B(1));
        data = src.Parent.Parent.Parent.UserData.LFP_data.Right.timeDomainDataTable.(src.String{src.Value});
    end
elseif contains(window_data_source,'Accel')
    if contains(window_data_source,'Left')
        A = src.Parent.Parent.Parent.UserData.Accel_data.Left.accelDataTable.DerivedTime;
        B = datetime(A,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        time = seconds(B-B(1));
        data = src.Parent.Parent.Parent.UserData.Accel_data.Left.accelDataTable.(src.String{src.Value});
    elseif contains(window_data_source,'Right')
        A = src.Parent.Parent.Parent.UserData.Accel_data.Right.accelDataTable.DerivedTime;
        B = datetime(A,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        time = seconds(B-B(1));
        data = src.Parent.Parent.Parent.UserData.Accel_data.Right.accelDataTable.(src.String{src.Value});
    end
elseif contains(window_data_source,'Delsys')
    time = src.Parent.Parent.Parent.UserData.Delsys_data.out_struct.Time.(src.String{src.Value});
    data = src.Parent.Parent.Parent.UserData.Delsys_data.out_struct.Data.(src.String{src.Value});
elseif contains(window_data_source,'Xsens')
    time = src.Parent.Parent.Parent.UserData.Xsens_data.Time;
    data = src.Parent.Parent.Parent.UserData.Xsens_data.(src.String{src.Value});
elseif contains(window_data_source,'Force')
    time = (0:height(src.Parent.Parent.Parent.UserData.FP_data)-1)/1000;
    data = src.Parent.Parent.Parent.UserData.FP_data.(src.String{src.Value});
elseif contains(window_data_source,'Teensey')
    time = src.Parent.Parent.Parent.UserData.Teensey_data.Time;
    data = src.Parent.Parent.Parent.UserData.Teensey_data.(src.String{src.Value});
end

% Plot data
cla(src.Parent.Parent.Parent.UserData.plot_axes(window));
plot(src.Parent.Parent.Parent.UserData.plot_axes(window),time,data,'-k');
xlabel(src.Parent.Parent.Parent.UserData.plot_axes(window),'Time (s)');
title(src.Parent.Parent.Parent.UserData.plot_axes(window),[window_data_source,': ',src.String{src.Value}],'Interpreter','none');

% Plot annotations 
plotAlignmentTimes;

% TODO: Add y label based on the data 
end