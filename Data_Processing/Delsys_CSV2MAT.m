function out_struct = Delsys_CSV2MAT(filename,savepath)
%{ DelsysCSV2Mat
% Converts exported Delsys CSV file to a struct. The struct has four
% fields: Chan_names, srates, Data, and Time. "Chan_names" are the names
% given to the signal in the Delsys data collection software. "srates" are
% the sampling rates for each recorded channel. "Data" is a struct with
% each field the name of the channel, and the signal data is truncated to
% only the number of points Delsys said it has collected. The number of
% points collected can be found in the raw Delsys CSV file. "Time" is also
% a struct and is laid out exactly like "Data," but with time instead.
%
% INPUTS:   filename    [=] The filename of the Delsys CSV file to convert.
%                           If the current matlab directory is not in the
%                           same folder as the file, the entire file path
%                           must be supplied
%
% OUTPUTS:  out_struct  [=] A structure with fields: Chan_names, srates,
%                           Data, and Time
% Author:   Kenneth Louie
% Date:     11/17/2020
%}

stopWatch = tic;

[out_struct,data_start_ind] = readHeader(filename);

raw_data = readtable(filename,'HeaderLines',data_start_ind+1);

for i = 1:length(out_struct.Chan_names)
    nDataPointsSignal = length(out_struct.Time.(out_struct.Chan_names{i}));
    out_struct.Data.(out_struct.Chan_names{i}) = raw_data{1:nDataPointsSignal,i*2};
end

[default_path,default_filename,default_ext] = fileparts(filename);
if exist('savepath','var') && ~isempty(savepath)
    savename = fullfile(savepath,[default_filename,'.mat']);
else
    savename = fullfile(default_path,[default_filename,'.mat']);
end
save(savename,'out_struct');
end

function [out_struct,data_start_ind] = readHeader(filename)
fid = fopen(filename);

stopFlag = 0;
out_struct.srates = [];
out_struct.Chan_names = {};
count = 1;
while stopFlag == 0
   cur_line = fgetl(fid);
   inds = extractPropertyIndices(cur_line);
   if sum(cellfun(@isempty,inds)) == 0
       % variable name
       var_name = replaceInvalidVarName(cur_line(inds{1}+7:inds{2}-2));
       var_name = strrep(var_name,' ','_');

       out_struct.Chan_names{end+1,1} = var_name;
       
       % sampling frequency
       Fs = str2double(cur_line(inds{2}+20:inds{3}-2));
       out_struct.srates(end+1,1) = Fs;
       
       % empty time and data vectors
       nPoints = round(str2double(cur_line(inds{3}+18:inds{4}-2)));
       out_struct.Data.(var_name) = zeros(1,nPoints);
       out_struct.Time.(var_name) = (0:nPoints-1).*(1/Fs);
       
       count = count + 1;
   end
   
   if strcmp(cur_line(1),'X')
       stopFlag = 1;
       data_start_ind = count;
       fclose(fid);
   end
end

end

function indices = extractPropertyIndices(cur_line)
f = @(x) strfind(cur_line,x);
indices = cellfun(f,{'Label','Sampling frequency','Number of points','start'},'UniformOutput',false);
end

function new_name = replaceInvalidVarName(cur_name)
new_name = cur_name;

reserved_char = [':','(',')','.','^','+','-','*','/','>','<','[',']','{','}'];

% Get rid of any unique characters reserved for matlab
for i = 1:length(reserved_char)
    new_name = strrep(new_name,reserved_char(i),'');
end

end