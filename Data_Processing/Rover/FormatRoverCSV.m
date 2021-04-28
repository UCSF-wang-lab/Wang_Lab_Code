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

sw_ver = [];
header_names = {};
date_vec = [];
time_vec = [];
data_mat = [];

while ~feof(fid)
    line = fgetl(fid);
    
    if contains(line,'SW version: ')
        idx = strfind(line,'SW version: ');
        space_idx = strfind(line,' ');
        s_idx = idx+12;
        e_idx = space_idx(find(space_idx > s_idx,1,'first'))-2;
        sw_ver = line(s_idx:e_idx);
    end
    
    if contains(line,'Time')
        header_names = strsplit(line,', ');
    end
    
    if ~isempty(sw_ver) && ~isempty(header_names) && ~isempty(line) && ~contains(line,'Time')
        vals = strsplit(line,',');
        for i = 1:length(vals)
            if i == 1
                temp = strsplit(vals{i},'/');
                date_vec{end+1,1} = datetime(year(datetime),str2double(temp{1}),str2double(temp{2}));
            elseif i == 2
                temp = strsplit(vals{i},':');
                if length(temp{1}) == 1
                    temp{1} = ['0',temp{1}];
                end
                
                if length(temp{3}) == 1
                    temp{3} = ['0',temp{3}];
                end
                time_vec{end+1,1} = temp{1} + ":" + temp{2} + ":" + temp{3};
            else
                data_mat(size(date_vec,1),i-2) = str2double(vals{i});
            end
        end
    end
    
end

time = createDateTime(date_vec,time_vec);
datetime_vec = datetime(time,'InputFormat','dd-MMM-yyyy HH:mm:ss.SSS');
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
time_out = [];

if length(date_vec) ~= length(time_vec)
    error('Number of dates and times do not match.');
end

chunk_time = time_vec{1};
chunk_date = date_vec{1};
for i = 2:length(time_vec)
    if time_vec{i} == chunk_time(1) && date_vec{i} == chunk_date(1)
        chunk_time = [chunk_time;time_vec{i}];
        chunk_date = [chunk_date;date_vec{i}];
    else
        temp = strsplit(num2str(990-(length(chunk_date)-1)*1000/fs:1000/fs:990,'% 04i')," ")'; % Need one more in field width to account for space
        mills = "." + temp;
        time_out = [time_out; datestr(chunk_date)+ " " + chunk_time + mills];
        
        chunk_time = time_vec{i};
        chunk_date = date_vec{i};
    end
end

if (length(time_out) ~= length(time_vec))
    temp = strsplit(num2str(990-(length(chunk_date)-1)*1000/fs:1000/fs:990,'% 04i')," ")'; % Need one more in field width to account for space
        mills = "." + temp;
        time_out = [time_out; datestr(chunk_date)+ " " + chunk_time + mills];
end

end