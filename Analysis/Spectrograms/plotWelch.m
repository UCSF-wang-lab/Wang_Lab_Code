function plotWelch(varargin)
% Calculates the grand average gait cycle average coherogram from all the
% files passed into this function. Multiple optional inputs allows
% customization of the computation including resolution of the gait cycle,
% and normalization strategy.
%
% INPUTS:  Optional
%             data_fullpath[=] Full path to the aligned datafile to be
%                               analyzed
%
%                   savePlot[=] Boolean option to save the resulting plots.
%                               Default is false.
%
% Date:     05/19/2024
% Author: Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
for i = 1:2:nargin
    switch varargin{i}
        case 'data_fullpath'
            data_fullpath = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

%% Determine parent directory
correct_path = false;
while ~correct_path
    parent_dir = uigetdir();
    if ismac
        split_parent_dir = strsplit(parent_dir,'/');
    elseif ispc
        split_parent_dir = strsplit(parent_dir,'\');
    else
        error('Platform not supported');
    end

    if strcmp(split_parent_dir{end},'Data')
        correct_path = true;
    else
        warning('Please select the folder called "Data"');
    end
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    plot_save_path = fullfile(parent_dir,'Figures','GSLT','WelchPSD',condition);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end


%% Load the aligned datafile
if ~exist('data_fullpath','var') || isempty(data_fullpath)
    [aligned_fname,data_folder] = uigetfile('*w_Gait_Events_Eleni_Adaptive.mat','Select the file'); 
    load(fullfile(data_folder,aligned_fname));
else
    load(fullfile(data_fullpath));
end

%% Compute the Welch periodogram
fs = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);

[pxx,f] = pwelch(aligned_data.left_LFP_table.key3,[],[],[],fs);

figure; plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

end