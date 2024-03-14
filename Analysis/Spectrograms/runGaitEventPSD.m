% Load aligned data
[filename_aligned_data,path] = uigetfile('*.mat');
load(fullfile(path,filename_aligned_data));
% fs_right = aligned_data.DeviceSettings.Right.fftSettings.TDsampleRates;
fs_left = aligned_data.DeviceSettings.Left.fftSettings.TDsampleRates;
% assert(fs_left==fs_right);
fs = fs_left;

% Calculate the STFT
windowLength = 1;
win = round(fs*windowLength);
win = win + mod(win,2);
percentOverlap = 0.9;
noverlap = round(win*percentOverlap);
nfft = 2^nextpow2(fs);
gapFillType = [];
STFT = calcRCS_STFT(aligned_data,gapFillType,windowLength,percentOverlap,nfft);

% Calculate & plot the PSD
GaitEventPSD(aligned_data,STFT)