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
                if length(temp{3})
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
out_table = array2table(data_mat);
out_table = [time,out_table];
out_table.Properties.VariableNames = header_names{2:end};

end

function time_out = createDateTime(date_vec,time_vec)
time_out = [];

end