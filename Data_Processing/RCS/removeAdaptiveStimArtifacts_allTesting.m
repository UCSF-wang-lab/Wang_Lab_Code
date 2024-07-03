function removeAdaptiveStimArtifacts_allTesting(filename,varargin)
% Computes the freeze index and wavelet coherence for a desired time
% interval of the walking task
%
% INPUTS: Required
%           filename[=]     The name of the aligned datafile. The full path
%                           needs to be provided.
%         Optional
%          save_plot[=]     Flag for spectrogram plotting (1=yes,0=no).
%                           Default = 0
%                           
% OUTPUTS:       NONE (for now)
% Example call: removeAdaptiveStimArtifacts(filename)
%
% Author:      Eleni Patelaki
% Date:     12/7/23

%% Option variables
for i = 1:2:nargin-1
    switch varargin{i}
        case 'save_plot'
            save_plot = varargin{i+1};
    end
end

%% Set default values for optional parameters
if ~exist('plot_spec','var') || isempty(save_plot)
    save_plot = 0;
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
lfp_data_keyA = lfp_struct.aligned_data.left_LFP_table.key0;
lfp_data_keyB = lfp_struct.aligned_data.left_LFP_table.key1;
% lfp_data_keyC = lfp_struct.aligned_data.left_LFP_table.key2;
% lfp_data_keyD = lfp_struct.aligned_data.left_LFP_table.key3;
fs = lfp_struct.aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);

%% Extract adaptive stimulation data
val_stimArtifact = lfp_struct.aligned_data.left_adaptive_table.CurrentProgramAmplitudesInMilliamps(:,1);
t_stimArtifact_left = lfp_struct.aligned_data.left_adaptive_taxis;

%% Adaptive Filtering
% First, resample the signals
[val_stimArtifact_res,t_stimArtifact_res] = resample(val_stimArtifact,t_stimArtifact_left,500);
[lfp_data_keyA_res,t_lfp_res] = resample(lfp_data_keyA,t_lfp,500);
[lfp_data_keyB_res,~] = resample(lfp_data_keyB,t_lfp,500);

% Truncate to the same size
minLen = min(length(val_stimArtifact_res),length(lfp_data_keyA_res));
val_stimArtifact_res = val_stimArtifact_res(1:minLen);
t_stimArtifact_res = t_stimArtifact_res(1:minLen);
lfp_data_keyA_res = lfp_data_keyA_res(1:minLen);
lfp_data_keyB_res = lfp_data_keyB_res(1:minLen);
t_lfp_res = t_lfp_res(1:minLen);

% Filtering
filt.Numerator = fir1(51,0.2,'bandpass');
coeffs = (filt.Numerator).'-0.01; % Set the filter initial conditions.
mu = 0.0001; % Set the step size for algorithm updating.
lms = dsp.LMSFilter(10,'Method','LMS','StepSize',mu,'InitialConditions',0);
% lms = dsp.LMSFilter(12,'Method','Sign-Data LMS','StepSize',mu);

[~,lfp_data_keyA_filt_adaptive] = lms(lfp_data_keyA_res,val_stimArtifact_res);


%% Classic filtering
Fstop1 = 4;
Fpass1 = 5;
Fpass2 = 50;
Fstop2 = 55;

% Gains
Astop1 = 65;
Astop2 = 65;
Apass = 1;

% Generate high/low-pass filters     
hpSpecs = fdesign.highpass(Fstop1,Fpass1,Astop1,Apass,fs);
hpf = design(hpSpecs,'cheby2','MatchExactly','stopband'); 
lpSpecs = fdesign.lowpass(Fpass2,Fstop2,Apass,Astop2,fs);
lpf = design(lpSpecs,'cheby2','MatchExactly','stopband'); 
lfp_data_keyA_filt_classic = filtfilt(hpf.sosMatrix,hpf.scaleValues,lfp_data_keyA_filt_adaptive);
lfp_data_keyA_filt_classic = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_keyA_filt_classic);

lfp_data_keyA_filt_classic_raw = filtfilt(hpf.sosMatrix,hpf.scaleValues,lfp_data_keyA);
lfp_data_keyA_filt_classic_raw = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_keyA_filt_classic_raw);

lfp_data_keyB_filt_classic_raw = filtfilt(hpf.sosMatrix,hpf.scaleValues,lfp_data_keyB);
lfp_data_keyB_filt_classic_raw = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_keyB_filt_classic_raw);


%% Polynomial fit and detrending
% opol = 15;
% [p,~,mu] = polyfit(t_lfp,lfp_data_keyA,opol);
% f_y = polyval(p,t_lfp,[],mu);
% 
% lfp_data_keyA_filt_detrended = lfp_data_keyA_filt_classic - f_y;

%% Wavelet denoising
lfp_data_keyA_filt_wden = wdenoise(lfp_data_keyA,'DenoisingMethod','SURE','ThresholdRule','Hard','NoiseEstimate','LevelDependent');


%% ICA
% [Zica, W, T, mu] = fastICA([lfp_data_keyA(1:10000) lfp_data_keyB(1:10000)],2);
% noisy_signals = [lfp_data_keyA,lfp_data_keyB,lfp_data_keyC,lfp_data_keyD];
noisy_signals = [lfp_data_keyA,lfp_data_keyB];
Mdl = rica(noisy_signals,2);
ica_activations = Mdl.TransformWeights*noisy_signals';
figure; 
subplot(2,1,1);
plot(ica_activations(1,:));
title('Component 1');
subplot(2,1,2);
plot(ica_activations(2,:));
title('Component 2');

invWeights = pinv(Mdl.TransformWeights);
prompt = "Which is the neural component? ";
comp2keep= input(prompt);
clean_signals = invWeights(:,comp2keep)*ica_activations(comp2keep,:);

figure;
subplot(2,1,1);
plot(clean_signals(1,:));
title('ICA-cleaned subcortical 1st key');
subplot(2,1,2);
plot(clean_signals(2,:));
title('ICA-cleaned subcortical 2nd key');

lfp_data_keyA_filt_classic_ica = filtfilt(hpf.sosMatrix,hpf.scaleValues,clean_signals(2,:)');
lfp_data_keyA_filt_classic_ica = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_keyA_filt_classic_ica);

lfp_data_keyB_filt_classic_ica = filtfilt(hpf.sosMatrix,hpf.scaleValues,clean_signals(1,:)');
lfp_data_keyB_filt_classic_ica = filtfilt(lpf.sosMatrix,lpf.scaleValues,lfp_data_keyB_filt_classic_ica);

%% Plotting
figure; 
% plot(t_lfp,lfp_data_keyB); 
% hold on; plot(t_lfp_res,lfp_data_keyA_filt_adaptive); 
% hold on; 
% plot(t_lfp_res,lfp_data_keyA_filt_classic);
% plot(t_lfp,ica_activations(1,:)')
hold on; plot(t_lfp,lfp_data_keyA_filt_classic_raw);
% hold on; plot(t_lfp,lfp_data_keyB_filt_classic_raw);
% hold on; plot(t_lfp,clean_signals(2,:)');
hold on; plot(t_lfp, lfp_data_keyA_filt_classic_ica);
hold on; plot(t_lfp, lfp_data_keyB_filt_classic_ica);

% hold on; plot(t_lfp, lfp_data_keyA_filt_detrended);
% hold on; plot(t_lfp, lfp_data_keyA_filt_wden);
% legend({'raw','adaptive filtered','classic filtered','matched'});
% legend({'raw','adaptive filtered','classic filtered'});
% legend({'classic adaptive','classic raw', 'classic ica'});

% legend({'raw key 1','filt key 0', 'filt key 1', 'ica + filt key 1'});
% legend({'filt key 0', 'filt key 1', 'ica + filt key 0', 'ica + filt key 1'});

legend({'filt key 0', 'ica + filt key 0', 'ica + filt key 1'});


end