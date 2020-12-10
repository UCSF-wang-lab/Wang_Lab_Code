function varargout = processTargetStepInfo(filename,savename)
% Removes any erroneous steps (multiple hits of same target) and saves as a
% .mat file containing the original table, the original table sorted by
% target number, and a filtered sorted table.
%
% Inputs:   file_name   [=] File of the extracted target information
%           save_path   [=] Save name and path of the .mat file
% Outputs:  varargout{1}[=] Resultant sorted and filtered table of targets
%           varargout{2}[=] Original table but sorted by target number
%           varargout{3}[=] The original table
%
% Author:   Kenneth Louie
% Date:     10/12/20

data_table = readtable(filename);
data_table_sorted = sortrows(data_table);

removeInds = [];
for i = 2:height(data_table_sorted)
    if data_table_sorted.TargetNumber(i-1) == data_table_sorted.TargetNumber(i) ||...
            data_table_sorted.TargetNumber(i-1) > data_table_sorted.TargetNumber(i)
        removeInds(end+1) = i;
    end
end

data_table_sorted_filt = data_table_sorted;
data_table_sorted_filt(removeInds,:) = [];

% Save data
if ~exist('savename','var')
    [A,B,~] = fileparts(filename);
    savename = fullfile(A,B)+"_target.csv";
end

try
    writetable(data_table_sorted_filt,savename);
catch
    error('Unable to write table to file');
end

% Prepare output of function
varargout{1} = data_table_sorted_filt;
varargout{2} = data_table_sorted;
varargout{3} = data_table;
end