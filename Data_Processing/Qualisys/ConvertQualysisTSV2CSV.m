%% ConvertQualisysTSV2CSV
%   This file takes in a TAB separated file (.tsv) from Qualysis software
%   and converts it to a comma separated file (.csv) after removing force
%   place calibration values and extra fluff that is not pure data.
%
%   Inputs:     filename    [=] The full path to the .tsv file to convert.
%
%   Outputs:    Data_table  [=] The raw data from the .tsv file.
%
%   Author: Kenneth Louie, PhD
%   Email:  kenneth.louie@ucsf.edu; ken.louie9@gmail.com
%   Date:   05/15/2023
function Data_table = ConvertQualysisTSV2CSV(filename)

if ~exist('filename','var') || isempty(filename)
    [file,path] = uigetfile('*.tsv');
    filename = fullfile(path,file);
end

try
    fid = fopen(filename);
catch
    error('Full file path should be passed in');
end

% Estimate number of values and get extra info
n_samples = [];
header_names = {};

count = 0;
while ~feof(fid)
    line = fgetl(fid);
    count = count + 1;

    if contains(line,'NO_OF_SAMPLES')
        str_parts = strsplit(line,'\t');
        n_samples = str2double(str_parts{2});
    end

    if contains(line,'COP')
        header_names = strsplit(line,'\t');
        n_headers = length(header_names);
        header_line_number = count;
    end

    if (~isempty(n_samples) && ~isempty(header_names)) && ~exist('Data_mat','var')
        Data_mat = nan(n_samples,length(header_names));
    end

    if exist('header_line_number','var') && (count > header_line_number)
        str_parts = strsplit(line,'\t');
        Data_mat(count-header_line_number,:) = cellfun(@str2double,str_parts(1:n_headers));
    end

end

Data_table = array2table(Data_mat,'VariableNames',header_names);
end