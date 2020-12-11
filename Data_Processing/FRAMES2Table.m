function FRAMES2Table(filename,savepath)
% FRAMES2Table: used to convert Vicon force plate data into a csv format.
% Author:   Kenneth Louie
% Date:     12/09/2020
[col_names,data_blocks] = getColNames(filename);

data_table = array2table(readmatrix(filename,'Range',[data_blocks(1,1),1,data_blocks(1,2),length(col_names{1})]));
data_table.Properties.VariableNames = col_names{1};

data_table_per_second = array2table(readmatrix(filename,'Range',[data_blocks(2,1),1,data_blocks(2,2),length(col_names{2})]));
data_table_per_second.Properties.VariableNames = col_names{2};

[dir,file,ext] = fileparts(filename);

if ~exist('savepath','var') || isempty(savepath)
    save_name = fullfile(dir,['ForcePlateData',ext]);
    save_name_per_second = fullfile(dir,['ForcePlateData_perSecond',ext]);
else
    save_name = fullfile(savepath,[file,'_ForcePlateData',ext]);
    save_name_per_second = fullfile(savepath,[file,'_ForcePlateData_perSecond',ext]);
end
writetable(data_table,save_name);
writetable(data_table_per_second,save_name_per_second);
end

function [col_names,data_blocks] = getColNames(filename)
fid = fopen(filename);

count = 1;
col_names = [];
data_blocks = [];
units = [];
while ~feof(fid)
    cur_line = fgetl(fid);
    splits = strsplit(cur_line,',','CollapseDelimiters',false);
    if ~isempty(splits{1}) && strcmp(splits{1},'Frame')
        % Column names
        col_names{end+1} = splits;
        count = count + 1;
        
        % Units
        cur_line = fgetl(fid);
        splits = strsplit(cur_line,',','CollapseDelimiters',false);
        units{end+1} = splits;
        
        data_blocks(end+1,1) = count + 1;
        
    elseif sum(cellfun(@(x) isempty(x),splits)) == length(splits)
        data_blocks(end,2) = count-1;
    elseif feof(fid)
        data_blocks(end,2) = count;
    end
    count = count + 1;
end
fclose(fid);

for i = 1:length(col_names)
    for j = 1:length(col_names{i})
        if ~isempty(units{i}{j})
            col_names{i}{j} = [col_names{i}{j},'_',units{i}{j}];
        end
    end
end

col_names{1} = cellfun(@(x) strrep(x,' ',''),col_names{1},'UniformOutput',false);
col_names{1} = cellfun(@(x) strrep(x,'.','_'),col_names{1},'UniformOutput',false);

col_names{2} = cellfun(@(x) strrep(x,' ',''),col_names{2},'UniformOutput',false);
col_names{2} = cellfun(@(x) strrep(x,'.','_'),col_names{2},'UniformOutput',false);
col_names{2} = cellfun(@(x) strrep(x,'/','_'),col_names{2},'UniformOutput',false);
col_names{2} = cellfun(@(x) strrep(x,'''',''),col_names{2},'UniformOutput',false);
end