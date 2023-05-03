function FormatRoverCSV(filename)

if ~exist('filename','var') || isempty(filename)
    [file,path] = uigetfile('*.csv');
    filename = fullfile(path,file);
end

try
    fid = fopen(filename);
catch
    error('Full file path should be passed in');
end

% Estimate number of values and get extra info
sw_ver = [];
header_names = {};

count = 0;
while ~feof(fid)
    line = fgetl(fid);
    count = count + 1;

    if contains(line,'SW version: ')
        idx = strfind(line,'SW version: ');
        space_idx = strfind(line,' ');
        s_idx = idx+12;
        e_idx = space_idx(find(space_idx > s_idx,1,'first'))-2;
        sw_ver = line(s_idx:e_idx);
    end

    if contains(line,'Time')
        header_names = strsplit(line,', ');
        header_line_number = count;
    end

end
approx_n_entries = (count - 8)/2;
frewind(fid);

date_vec = cell(approx_n_entries,1);
time_vec = cell(approx_n_entries,1);
data_mat = nan(approx_n_entries,length(header_names)-2);    % two less entries because the date and time will be separated 

count = 0;
entry_count = 0;
while ~feof(fid)
    line = fgetl(fid);
    count = count + 1;
    
    if count > header_line_number + 1 && ~isempty(line)
        entry_count = entry_count + 1;
        vals = strsplit(line,',');
        for i = 1:length(vals)
            if i == 1
                temp = strsplit(vals{i},'/');
                date_vec{entry_count,1} = datetime(year(datetime),str2double(temp{1}),str2double(temp{2}));
            elseif i == 2
                temp = strsplit(vals{i},':');
                if length(temp{1}) == 1
                    temp{1} = ['0',temp{1}];
                end
                
                if length(temp{3}) == 1
                    temp{3} = ['0',temp{3}];
                end
                time_vec{entry_count,1} = temp{1} + ":" + temp{2} + ":" + temp{3};
            else
                data_mat(entry_count,i-2) = str2double(vals{i});
            end
        end
    end
    
end

% Create date time variables as it will be useful when plotting
time = createDateTime(date_vec,time_vec);
datetime_vec = datetime(time,'InputFormat','dd-MMM-yyyy HH:mm:ss.SSS');
datetime_vec.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';

% Create output table 
out_table = array2table(data_mat);
out_table = [array2table(time),array2table(datetime_vec),out_table];
out_table.Properties.VariableNames = {'RawStringDateTime','DateTime',header_names{3:end}};

[save_file, save_path] = uiputfile('*.*');
if ~isempty(save_file) && ~isempty(save_path)
    writetable(out_table,[fullfile(save_path,save_file),'.csv']);
    save([fullfile(save_path,save_file),'.mat'],'out_table');
    
end

end

function time_out = createDateTime(date_vec,time_vec)
fs = 100;
time_out = repmat("",length(time_vec),1);

if length(date_vec) ~= length(time_vec)
    error('Number of dates and times do not match.');
end

chunk_time = time_vec{1};
chunk_date = date_vec{1};
start_chunk_ind = 1;
for i = 2:length(time_vec)
    if time_vec{i} == chunk_time(1) && date_vec{i} == chunk_date(1) % Add to chunk of data
        chunk_time = [chunk_time;time_vec{i}];
        chunk_date = [chunk_date;date_vec{i}];

    else    % End of a chunk, create time strings
        % Check to see there are the correct amount of samples for the
        % sampling rate. Occasionally will have more samples than the
        % sample rate for some reason.
        if length(chunk_date) <= 100
            temp = strsplit(num2str(990-(length(chunk_date)-1)*1000/fs:1000/fs:990,'% 04i')," ")'; % Need one more in field width to account for space
            mills = "." + temp;
            time_out(start_chunk_ind:i-1) = datestr(chunk_date)+ " " + chunk_time + mills;

            start_chunk_ind = i;
            chunk_time = time_vec{i};
            chunk_date = date_vec{i};
        else
            temp = strsplit(num2str(990-(100-1)*1000/fs:1000/fs:990,'% 04i')," ")'; % Need one more in field width to account for space
            mills = "." + temp;
            time_out(start_chunk_ind:start_chunk_ind+99) = datestr(chunk_date(1:100))+ " " + chunk_time(1:100) + mills;
            
            start_chunk_ind = i;
            chunk_time = time_vec{i};
            chunk_date = date_vec{i};
        end
        
    end
end

if i == length(time_vec) && length(temp) > length(chunk_time)
    temp = strsplit(num2str(0:1000/fs:(length(chunk_date)-1)*10,'% 04i')," ")'; % Need one more in field width to account for space
    mills = "." + temp;
    time_out(start_chunk_ind:i) = datestr(chunk_date)+ " " + chunk_time + mills;
end

end