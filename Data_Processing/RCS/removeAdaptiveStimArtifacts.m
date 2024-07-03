function removeAdaptiveStimArtifacts(filename,varargin)
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
%           Fstop_hp[=]     Frequency at the end of the stopband
%                           (highpass filter). Default = 1.5
%
%           Fpass_hp[=]     Frequency at the start of the passband
%                           (highpass filter). Default = 2.5
%
%           Fpass_lp[=]     Frequency at the end of the passband band
%                           (lowpass filter). Default = 50
%
%           Fstop_lp[=]     Frequency at the start of the stopband
%                           (lowpass filter). Default = 60
%         
%              Apass[=]     Amount of ripples in the passband (in dB).
%                           Defualt = 1
%
%              Astop[=]     Attenuation in the stop band (in dB).
%                           Defualt = 65
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
        case 'Fstop_hp'
            Fstop_hp = varargin{i+1};
        case 'Fpass_hp'
            Fpass_hp = varargin{i+1};
        case 'Fpass_lp'
            Fpass_lp = varargin{i+1};
        case 'Fstop_lp'
            Fstop_lp = varargin{i+1};
        case 'Apass'
            Apass = varargin{i+1};
        case 'Astop'
            Astop = varargin{i+1};
    end
end

%% Set default values for optional parameters
if ~exist('show_plots','var') || isempty(show_plots)
    show_plots = 0;
end

if ~exist('Fstop_hp','var') || isempty(Fstop_hp)
    Fstop_hp = 1.5;
end

if ~exist('Fpass_hp','var') || isempty(Fpass_hp)
    Fpass_hp = 2.5;
end

if ~exist('Fpass_lp','var') || isempty(Fpass_lp)
    Fpass_lp = 50;
end

if ~exist('Fstop_lp','var') || isempty(Fstop_lp)
    Fstop_lp = 60;
end

if ~exist('Apass','var') || isempty(Apass)
    Apass = 1;
end

if ~exist('Astop','var') || isempty(Astop)
    Astop = 60;
end

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

%% Extract adaptive stimulation current waveform
val_stimArtifact_left = lfp_struct.aligned_data.left_adaptive_table.CurrentProgramAmplitudesInMilliamps(:,1);
t_stimArtifact_left = lfp_struct.aligned_data.left_adaptive_taxis;


%% Filter (before ICA)
% Generate high/low-pass filters     
hpSpecs = fdesign.highpass(Fstop_hp,Fpass_hp,Astop,Apass,fs);
hpf = design(hpSpecs,'cheby2','MatchExactly','stopband'); 
lpSpecs = fdesign.lowpass(Fpass_lp,Fstop_lp,Apass,Astop,fs);
lpf = design(lpSpecs,'cheby2','MatchExactly','stopband');

% Apply the filters
lfp_data_key0_filt = filtfilt(hpf.sosMatrix,hpf.scaleValues,lfp_data_key0);
lfp_data_key0_filt = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_key0_filt);

lfp_data_key1_filt = filtfilt(hpf.sosMatrix,hpf.scaleValues,lfp_data_key1);
lfp_data_key1_filt = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_key1_filt);

%% ICA
raw_lfps = [lfp_data_key0,lfp_data_key1];
ica_mdl = rica(raw_lfps,2);
invWeights = pinv(ica_mdl.TransformWeights);
ica_activations = ica_mdl.TransformWeights*raw_lfps';

% Plot the components to visually inspect them and select the neural one
figure; 
subplot(2,1,1);
plot(t_lfp,ica_activations(1,:));
title('Component 1');
subplot(2,1,2);
plot(t_lfp,ica_activations(2,:));
title('Component 2');

% Request input from the user to select the neural component to be kept
prompt = "Which is the neural component? ";
comps = 1:ica_mdl.NumLearnedFeatures;
comp2keep= input(prompt);
if ~ismember(comp2keep,comps)
    error('The component number provided is not valid');
end
icaclean_signals = invWeights(:,comp2keep)*ica_activations(comp2keep,:);

% Plot the noise component along with the stimulation current waveform to
% esnure that they look similar
if show_plots
    figure; 
    subplot(2,1,1);
    comp2rej = comps(comps~=comp2keep);
    t_start = max([t_lfp(1),t_stimArtifact_left(1)]);
    t_end = min([t_lfp(end),t_stimArtifact_left(end)]);
    
    plot(t_lfp,ica_activations(comp2rej,:));
    title('Noise component');
    xlim([t_start,t_end]);
    subplot(2,1,2);
    plot(t_stimArtifact_left,val_stimArtifact_left);
    title('Stimulation waveform');
    xlim([t_start,t_end]);
end 

% Plot the ica-clean signals
figure;
subplot(2,1,1);
plot(icaclean_signals(1,:));
title('1st ICA-cleaned subcortical signal');
subplot(2,1,2);
plot(icaclean_signals(2,:));
title('2nd ICA-cleaned subcortical signal');

% Request input from the user to select the ICA-cleaned signal that is most
% contaminated with stim adaptive stim artifacts (normally it's the stim
% contact, but not always...). Usually it's the one that has substantially
% lp-reater amplitude.
prompt = "Which ICA-cleaned signal is the most contaminated? ";
sigs = 1:ica_mdl.NumPredictors;
sig2keep= input(prompt);
if ~ismember(sig2keep,sigs)
    error('The signal number provided is not valid');
end

%% Filter (after ICA)
lfp_data_key1_ica_filt = filtfilt(hpf.sosMatrix,hpf.scaleValues,icaclean_signals(sig2keep,:)');
lfp_data_key1_ica_filt = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_key1_ica_filt);

%% Plot the results
if show_plots
    figure; 
    plot(t_lfp,lfp_data_key0_filt);
    hold on; plot(t_lfp,lfp_data_key1_ica_filt);
    legend({'filt key 0', 'ica + filt key 1'});
end

end