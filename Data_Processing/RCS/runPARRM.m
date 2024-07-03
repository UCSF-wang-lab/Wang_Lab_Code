function runPARRM(filename,varargin)
% Identifies and removes the adaptive stimulation artifacts using ICA and
% filtering
%
% INPUTS: Required
%           filename[=]     The name of the aligned datafile. The full path
%                           needs to be provided.
%         Optional
%
%         show_plots[=]     Flag for generating plots (1=yes,0=no).
%                           Default = 0
%                         
% OUTPUTS:       NONE (for now)
% Example call: removeAdaptiveStimArtifacts(filename)
%
% Author:      Eleni Patelaki
% Date:     4/8/23

%% Option variables
for i = 1:2:nargin-1
    switch varargin{i}
        case 'show_plots'
            show_plots = varargin{i+1};
    end
end

%% Set default values for optional parameters
if ~exist('show_plots','var') || isempty(show_plots)
    show_plots = 0;
end

%% Add the PARRM toolbox to the MATLAB path
addpath('C:\Users\el_pa\PARRM');

%% Load the mat file containing the aligned data
if isempty(filename)
    [filename,filepath] = uigetfile('*.mat');
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
lfp_struct = load(fullfile(filepath,filename));

%% Extract the LFP data
t_lfp = lfp_struct.aligned_data.left_taxis;
lfp_data_key0 = lfp_struct.aligned_data.left_LFP_table.key0;
lfp_data_key1 = lfp_struct.aligned_data.left_LFP_table.key1;
fs = lfp_struct.aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);

% t_lfp = aligned_data.left_taxis;
% lfp_data_key0 = aligned_data.left_LFP_table.key0;
% lfp_data_key1 = aligned_data.left_LFP_table.key1;
% fs = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);

%% PARRM
winSize=5000; % Width of the window in sample space for PARRM filter
skipSize=0; % Number of samples to ignore in each window in sample space
winDir='both'; % Filter using samples from the past and future
guessPeriod=fs/130; % Starting point for period grid search

% Find the period of stimulation in simulated data
Period=FindPeriodLFP(lfp_data_key0,[1,length(lfp_data_key0)-1],guessPeriod);
perDist=0.01; % Window in period space for which samples will be averaged
% Create the linear filter
PARRM=PeriodicFilter(Period,winSize,perDist,skipSize,winDir);
% Filter using the linear filter and remove edge effects
Filtered=((filter2(PARRM.',lfp_data_key0','same')-lfp_data_key0')./(1-filter2(PARRM.',ones(size(lfp_data_key0')),'same'))+lfp_data_key0')';
d = filter2(PARRM.',lfp_data_key0')-lfp_data_key0';

% %%
% % Plot the components to visually inspect them and select the neural one
% figure; 
% subplot(2,1,1);
% plot(t_lfp,ica_activations(1,:));
% title('Component 1');
% subplot(2,1,2);
% plot(t_lfp,ica_activations(2,:));
% title('Component 2');
% 
% % Request input from the user to select the neural component to be kept
% prompt = "Which is the neural component? ";
% comps = 1:ica_mdl.NumLearnedFeatures;
% comp2keep= input(prompt);
% if ~ismember(comp2keep,comps)
%     error('The component number provided is not valid');
% end
% icaclean_signals = invWeights(:,comp2keep)*ica_activations(comp2keep,:);
% 
% % Plot the ica-clean signals
% figure;
% subplot(2,1,1);
% plot(icaclean_signals(1,:));
% title('1st ICA-cleaned subcortical signal');
% subplot(2,1,2);
% plot(icaclean_signals(2,:));
% title('2nd ICA-cleaned subcortical signal');
% 
% % Request input from the user to select the ICA-cleaned signal that is most
% % contaminated with stim artifacts (normally it's the stim
% % contact, but not always...). Usually it's the one that has substantially
% % lp-reater amplitude.
% prompt = "Which ICA-cleaned signal is the most contaminated? ";
% sigs = 1:ica_mdl.NumPredictors;
% sig2keep= input(prompt);
% if ~ismember(sig2keep,sigs)
%     error('The signal number provided is not valid');
% end
% 
% %% Filter (after ICA)
% lfp_data_noisierKey_ica_filt = filtfilt(hpf.sosMatrix,hpf.scaleValues,icaclean_signals(sig2keep,:)');
% lfp_data_noisierKey_ica_filt = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_noisierKey_ica_filt);
% 
% %% Plot the results
% if show_plots
%     figure; 
%     plot(t_lfp,lfp_data_key0_filt);
%     hold on; plot(t_lfp,lfp_data_key1_filt);
%     hold on; plot(t_lfp,lfp_data_noisierKey_ica_filt);
%     legend({'filt noisier key','filt cleaner key', 'ica + filt noisier key'});
% end

end