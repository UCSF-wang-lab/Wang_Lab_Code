function FRAMES2Table(filename,savepath)
% FRAMES2Table: used to convert Vicon force plate data into a csv format.
%
% INPUTS:   filename    [=] Vicon csv file to load in.
%           savepath    [=] Path to save converted csv file. Optional.
%
% OUTPUTS:
%
% Author:   Kenneth Louie
% Date:     12/16/2020
[col_names,data_blocks] = getColNames(filename);

data_table = array2table(readmatrix(filename,'Range',[data_blocks(1,1),1,data_blocks(1,2),length(col_names{1})]));
data_table.Properties.VariableNames = col_names{1};

% Rotate coordinates so it is in the right directions and magnitude
if sum(contains(col_names{1},'FP')) == 0
    % Single force plate
    temp = data_table.Cx_mm;
    data_table.Cx_mm = data_table.Cy_mm;
    data_table.Cy_mm = -temp;
    data_table.Cz_mm = -data_table.Cz_mm;
    
    temp = data_table.Fx_N;
    data_table.Fx_N = data_table.Fy_N;
    data_table.Fy_N = -temp;
    data_table.Fz_N = -data_table.Fz_N;
end

[dir,file,ext] = fileparts(filename);

if ~exist('savepath','var') || isempty(savepath)
    save_name = fullfile(dir,['ForcePlateData',ext]);
else
    save_name = fullfile(savepath,[file,'_ForcePlateData',ext]);
end
writetable(data_table,save_name);

% if per second variables exist
if length(col_names) == 2
    data_table_per_second = array2table(readmatrix(filename,'Range',[data_blocks(2,1),1,data_blocks(2,2),length(col_names{2})]));
    data_table_per_second.Properties.VariableNames = col_names{2};
        
    temp = data_table_per_second.Cx_mm_s;
    data_table_per_second.Cx_mm_s = data_table_per_second.Cy_mm_s;
    data_table_per_second.Cy_mm_s = -temp;
    data_table_per_second.Cz_mm_s = -data_table_per_second.Cz_mm_s;
    
    temp = data_table_per_second.Fx_N_s;
    data_table_per_second.Fx_N_s = data_table_per_second.Fy_N_s;
    data_table_per_second.Fy_N_s = -temp;
    data_table_per_second.Fz_N_s = -data_table_per_second.Fz_N_s;
    
    if ~exist('savepath','var') || isempty(savepath)
        save_name_per_second = fullfile(dir,['ForcePlateData_perSecond',ext]);
    else
        save_name_per_second = fullfile(savepath,[file,'_ForcePlateData_perSecond',ext]);
    end
    
    writetable(data_table_per_second,save_name_per_second);
end

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
for i = 1:length(col_names{1})
    repeats =  cellfun(@(x) strcmp(x,col_names{1}{i}),col_names{1});
    if sum(repeats) > 1
        rep_inds = find(repeats);
        for j = 1:length(rep_inds)
            col_names{1}{rep_inds(j)} = [col_names{1}{rep_inds(j)},'_FP',num2str(j)];
        end
    end
end

if length(col_names) > 1
    col_names{2} = cellfun(@(x) strrep(x,' ',''),col_names{2},'UniformOutput',false);
    col_names{2} = cellfun(@(x) strrep(x,'.','_'),col_names{2},'UniformOutput',false);
    col_names{2} = cellfun(@(x) strrep(x,'/','_'),col_names{2},'UniformOutput',false);
    col_names{2} = cellfun(@(x) strrep(x,'''',''),col_names{2},'UniformOutput',false);
end
end