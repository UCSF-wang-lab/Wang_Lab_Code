function parseFiles(subjectID,varargin)
%% TODO
% Parse that looks at recording folders and looks in the events to see what
% the recording contains

for i = 1:2:nargin-3
    switch varargin{i}
        case 'export_data'
            export_data = varargin{i+1};
        case 'remove_nan'
            remove_nan = varargin{i+1};
        case 'n_percent_bins'
            n_percent_bins = varargin{i+1};
        case 'baseline_time'
            baseline_time = varargin{i+1};
        case 'cycle_start_event'
            cycle_start_event = varargin{i+1};
    end
end

datetime(A,'ConvertFrom','epochtime','TicksPerSecond',1e3,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
end

