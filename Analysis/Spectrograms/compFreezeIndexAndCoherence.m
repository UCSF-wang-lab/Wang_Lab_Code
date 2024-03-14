function compFreezeIndexAndCoherence(filename,varargin)
% Computes the freeze index and wavelet coherence for a desired time
% interval of the walking task
%
% INPUTS: Required
%           filename[=]     The name of the aligned datafile. The full path
%                           needs to be provided.
%         Optional
%             system[=]     The name of the system that the selected data
%                           was redorded with. Acceptable values are
%                           'xsens' and 'delsys'. For now, this code works
%                           only for accelerometry data.
%
%           landmark[=]     The anatomical location of the sensor of the
%                           selected data.
%                           For Delsys, acceptable values are 'sol' for soleus,
%                           and 'ta' for tibialis anterior
%                           For Xsens, acceptable values are 'lleg' for
%                           lower leg and 'uleg' for upper leg.
%
%               axis[=]     The motion axis to be analyzed. Acceptable
%                           values are 'x','y','z'.
%
%         start_time[=]     Start time of the time interval to be analyzed.
%                           Default = start of the entire walking period.
%
%           end_time[=]     End time of the time interval to be analyzed.
%                           Default = end of the entire walking period.
%
%          fi_method[=]     Method to use to compute the freeze index
%                           (FI). Acceptable values are 'wav' for wavelet
%                           transform spectrogram and 'four' Fourier transform 
%                           spectrogram. Default = 'wav'
%
%          plot_spec[=]     Flag for spectrogram plotting (1=yes,0=no).
%                           Default = 0
%
%           coh_pair[=]     Pair of contacts to compute coherence between.
%                           Example (and default): {'key0','key2'}
%                        
%                           
% OUTPUTS:       NONE (for now)
% Example call: compFreezeIndexAndCoherence(filename,'system','delsys','landmark','sol','axis','y','start_time',70,'end_time',100,'coh_pair',{'key0','key3'})
%
% Author:      Eleni Patelaki
% Date:     12/7/23

%% Option variables
for i = 1:2:nargin-1
    switch varargin{i}
        case 'system'
            sys = varargin{i+1};
        case 'landmark'
            loc = varargin{i+1};
        case 'axis'
            axis = varargin{i+1};
        case 'start_time'
            t_start = varargin{i+1};
        case 'end_time'
            t_end = varargin{i+1};
        case 'fi_method'
            fi_method = varargin{i+1};
        case 'plot_spec'
            plot_spec = varargin{i+1};
        case 'coh_pair'
            coh_pair = varargin{i+1};
    end
end

%% Set default options if not passed in by user
if ~exist('sys','var') || isempty(sys)
    sys = 'delsys';
end

if ~exist('loc','var') || isempty(loc)
    loc = 'sol';
end

if ~exist('axis','var') || isempty(axis)
    axis = 'y';
end

if ~exist('fi_method','var') || isempty(fi_method)
    fi_method = 'wav';
end

if ~exist('plot_spec','var') || isempty(plot_spec)
    plot_spec = 0;
end

if ~exist('coh_pair','var') || isempty(coh_pair)
    coh_pair = {'key0','key2'};
end

%% Run the necessary checks to make sure that the inputs are valid
if ~strcmpi(sys,'delsys') && ~strcmpi(sys,'xsens')
    error('Invalid system type.');
end

if strcmpi(sys,'delsys')
    if ~strcmpi(loc,'sol') && ~strcmpi(loc,'ta')
        error('Invalid anatomical landmark for Delsys.');
    end
end

if strcmpi(sys,'xsens')
    if ~strcmpi(loc,'lleg') && ~strcmpi(loc,'uleg')
        error('Invalid anatomical landmark for Xsens.');
    end
end

if strcmpi(sys,'xsens')
    if ~strcmpi(loc,'lleg') && ~strcmpi(loc,'uleg')
        error('Invalid anatomical landmark for Xsens.');
    end
end

if ~strcmpi(fi_method,'wav') && ~strcmpi(fi_method,'four')
    error('Invalid freeze index computation method.');
end

keys = {'key0','key1','key2','key3'};
if length(coh_pair)~=2 || ~matches(coh_pair{1},keys) || ~matches(coh_pair{2},keys)
    error('Invalid contact pair for coherence computation.');
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
load(fullfile(filepath,filename));

%% Extract the accelerometry data
if strcmpi(sys,'delsys')
    if strcmpi(loc,'sol')
        datastr = strcat('Soleus_Right_8_ACC',upper(axis),'_8');
    elseif strcmpi(loc,'ta')
        datastr = strcat('Tibialis_Anterior_Right_6_Acc_6',upper(axis),'_IM');
    end
    accel_time = getfield(aligned_data.Delsys.Time,datastr);
    accel_data = getfield(aligned_data.Delsys.Data,datastr);
    fs = aligned_data.Delsys.srates(find(strcmp(aligned_data.Delsys.Chan_names,datastr)==1));
elseif strcmpi(sys,'xsens')
    if strcmpi(loc,'lleg')
        datastr = strcat('RightLowerLeg_Acc',upper(axis));
    elseif strcmpi(loc,'uleg')
        datastr = strcat('RightUpperLeg_Acc',upper(axis));
    end
    accel_time = getfield(aligned_data.Xsens,'Time');
    accel_data = getfield(aligned_data.Xsens,datastr);
    fs = 60;
end

%% Extract the LFP data
lfp_time = aligned_data.left_taxis;
lfp_data_keyA = aligned_data.left_LFP_table.key1;
lfp_data_keyB = aligned_data.left_LFP_table.key2;
sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);

%%  Define necessary variables
% Start & end index of the entire walking period
walking_start_ind = 1;
walking_end_ind = length(lfp_time);
% walking_start_ind = find(lfp_time >= t_start,1,'first');
% walking_end_ind = find(lfp_time <= t_end,1,'last');

% Start & the end time of the desired time interval
if ~exist('t_start','var') || isempty(t_start)
    t_start = lfp_time(1);
end

if ~exist('t_end','var') || isempty(t_end)
    t_end = lfp_time(end);
end


%% Compute the normalized magnitude-squared LFP coherence 
% Compute the magnitude-squared coherence of the entire walking period
[wcoh,~,freq_spec] = wcoherence(lfp_data_keyA,lfp_data_keyB,sr,'VoicesPerOctave',12);

% Truncate portions of matrices x,y that correspond to frequencues lower than 0.1Hz and higher than 140Hz
wcoh = wcoh(freq_spec>=0.1&freq_spec<=140,:);
freq_spec = freq_spec(freq_spec>=0.1&freq_spec<=140,:);

%  Define the coherence normalization vector based on the entire walking
%  period
norm_vec = [mean(wcoh(:,walking_start_ind:walking_end_ind),2,'omitnan'),...
            std(wcoh(:,walking_start_ind:walking_end_ind),0,2,'omitnan'),...
            median(wcoh(:,walking_start_ind:walking_end_ind),2,'omitnan')];

% Normalize magnitude-squared Coherence
norm_wcoh= (wcoh-norm_vec(:,1))./norm_vec(:,2); % z-scoring
% norm_wcoh= (wcoh-norm_vec(:,1)); % baseline subtraction

%% Compute the freeze Index (FI)
% Compute the power spectrum
if strcmpi(fi_method,'wav')
    % Wavelets
    [wt,freq_spec] = cwt(accel_data,fs);
    spec = abs(wt).^2;
    time_spec = accel_time; 
elseif strcmpi(fi_method,'four')
    % Fourier
    windowLength = 3;
    win = round(fs*windowLength);
    win = win + mod(win,2);
    percentOverlap = 0.9;
    noverlap = round(win*percentOverlap);
    nfft = 2^nextpow2(fs);
    
    [ft,freq_spec,time_spec,~]=spectrogram(accel_data,win,noverlap,nfft,fs,'psd');
    spec = abs(ft).^2;
end

% Plot the spectrogram for the selected time interval (optionally)
if plot_spec
    figure;
    ax = pcolor(time_spec,log2(freq_spec),spec);
    ticks = logspace(log10(0.5),log10(8),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(0.5),log2(8)]);
    xlim([t_start t_end]); 
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title('Acceleration Power Spectrum');
end

% Calculate the delta power
pow_delta = spec(freq_spec>=0.5&freq_spec<=3,:);
pow_delta_avg = mean(pow_delta,1);

% Calculate the theta power
pow_theta = spec(freq_spec>=3&freq_spec<=8,:);
pow_theta_avg = mean(pow_theta,1);

% Calculate the FI as the theta/delta ratio
fi = pow_theta_avg./pow_delta_avg;

%% Plot the acceleration, FI and wavelet coherence during the desired time interval
figure; 

subplot(3,1,1); plot(accel_time,accel_data,'k');
xlim([t_start t_end]);
set(gca,'FontSize',12);
ylabel('Acceleration (g)');

subplot(3,1,2); plot(accel_time,fi,'r','LineWidth',1.2);
xlim([t_start t_end]); 
set(gca,'FontSize',12);
ylabel('Freeze Index');

subplot(3,1,3); 
ax = pcolor(lfp_time,log2(freq_spec),norm_wcoh);
ticks = logspace(log10(2.5),log10(50),5);
ax.Parent.YTick = log2(ticks);
ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
ylim([log2(2.5),log2(50)]);
xlim([t_start t_end]);
shading interp;
colormap jet;
max_abs=max(norm_wcoh,[],'all');
clim([-max_abs,max_abs]);
set(gca,'FontSize',12);
ylabel('Frequency (Hz)');
xlabel('Time (s)');
subtitle('Coherence: GPi-M1 (Left)');

end