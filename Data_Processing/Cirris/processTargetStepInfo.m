function varargout = processTargetStepInfo(filename,savepath)
% Removes any erroneous steps (multiple hits of same target) and saves as a
% .mat file containing the original table, the original table sorted by
% target number, and a filtered sorted table.
%
% Inputs:       filename[=]     File of the extracted target information
%               savepath[=]     Save name and path of the .mat file
%
% Outputs:  varargout{1}[=]     Resultant sorted and filtered table of targets
%           varargout{2}[=]     Original table but sorted by target number
%           varargout{3}[=]     The original table
%
% Author:   Kenneth Louie
% Contributor: Eleni Patelaki
% Date:     10/12/20

if isempty(filename)
    [filename,filepath] = uigetfile('*all_targets_aligned.csv');
else
    if ispc
        fileparts = strsplit(filename,'\');
    elseif ismac
        fileparts = strsplit(filename,'/');
    else
        error('Platform not currently supported.');
    end
    filepath  = fullfile(fileparts{1:end-1});
    filename =  fullfile(fileparts{end});
end

data_table = readtable(fullfile(filepath,filename));
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
if ~exist('savepath','var')||isempty(savepath)
    % [~,fname,~] = fileparts(filename);
    fname= extractBefore(filename,'.');
    fname = strrep(fname,'all_targets','processed_targets');
    savepath = fullfile(filepath,strcat(fname,'.csv'));
end

try
    writetable(data_table_sorted_filt,savepath);
catch
    error('Unable to write table to file');
end

% Prepare output of function
varargout{1} = data_table_sorted_filt;
varargout{2} = data_table_sorted;
varargout{3} = data_table;
end