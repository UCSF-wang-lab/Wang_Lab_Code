function shift_t = calcGLSTAlignment(varargin)
% Calculates the Cirris-Xsens alignment using cross-correlation
%
% Inputs:       plot_flag[=]     0 = does not save the acceleration and speed plots, 
%                                1 = saves the acceleration and speed plots
%                                Default is 0.
%
% Author:   Eleni Patelaki
% Date:     3/14/24

%% Parse optional arguments
for i = 1:2:nargin
    switch varargin{i}
        case 'plot_flag'
            plot_flag = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('plot_flag','var') || isempty(plot_flag)
    plot_flag = 0;
end

%% Read XSens kinematic data
[xsens_fname, xsens_fpath] = uigetfile('*.csv','Select Xsens csv');
table = readtable(fullfile(xsens_fpath,xsens_fname));

posX_R_Xsens = table.RightFoot_PosX;
posY_R_Xsens = table.RightFoot_PosY;
posZ_R_Xsens = table.RightFoot_PosZ;
velX_R_Xsens = table.RightFoot_VelX;
velY_R_Xsens = table.RightFoot_VelY;
velZ_R_Xsens = table.RightFoot_VelZ;
accX_R_Xsens = table.RightFoot_AccX;
accY_R_Xsens = table.RightFoot_AccY;
accZ_R_Xsens = table.RightFoot_AccZ;
posX_L_Xsens = table.RightFoot_PosX;
posY_L_Xsens = table.RightFoot_PosY;
posZ_L_Xsens = table.RightFoot_PosZ;
velX_L_Xsens = table.LeftFoot_VelX;
velY_L_Xsens = table.LeftFoot_VelY;
velZ_L_Xsens = table.LeftFoot_VelZ;
accX_L_Xsens = table.LeftFoot_AccX;
accY_L_Xsens = table.LeftFoot_AccY;
accZ_L_Xsens = table.LeftFoot_AccZ;
time = table.Time;

%% Read the Delsys data
[delsys_fname, delsys_fpath] = uigetfile('*.mat','Select Delsys');
delsys_file = load(fullfile(delsys_fpath,delsys_fname));
% chan_name = 'FSR_adapter_15_Left_FSRA_15';
chan_name = 'FSR_adapter_16_Right_FSRA_16';
delsys_hs_data = delsys_file.out_struct.Data.(chan_name);
delsys_hs_time = delsys_file.out_struct.Time.(chan_name);

%% Read Cirris foot kinematic data
[cirris_rfoot,cirris_path] = uigetfile('g*_RFoot.csv','Select right foot Cirris csv');
rfootData = readtable(fullfile(cirris_path,cirris_rfoot));
[csv_lfoot,cirris_path] = uigetfile('g*_LFoot.csv','Select left foot Cirris csv');
lfootData = readtable(fullfile(cirris_path,csv_lfoot));

timeR=rfootData.Time;
posX_R_Cirris=rfootData.PosX;
posY_R_Cirris=rfootData.PosY;
posZ_R_Cirris=rfootData.PosZ;
speedX_R_Cirris=rfootData.SpeedX;
speedY_R_Cirris=rfootData.SpeedY;
speedZ_R_Cirris=rfootData.SpeedZ;
timeL=lfootData.Time;
posX_L_Cirris=lfootData.PosX;
posY_L_Cirris=lfootData.PosY;
posZ_L_Cirris=lfootData.PosZ;
speedX_L_Cirris=lfootData.SpeedX;
speedY_L_Cirris=lfootData.SpeedY;
speedZ_L_Cirris=lfootData.SpeedZ;
accelX_L_Cirris=diff(speedX_L_Cirris)./diff(timeL);
accelY_L_Cirris=diff(speedY_L_Cirris)./diff(timeL);
accelZ_L_Cirris=diff(speedZ_L_Cirris)./diff(timeL);

%% Read the csv file containing the Cirris hit/miss info
[csv_hits,cirris_path] = uigetfile('g*_processed_targets.csv','Select Cirris performance csv');
hitData = readtable(fullfile(cirris_path,csv_hits));
left_hits = hitData(strcmp(hitData.Side,'Left'),:);
right_hits = hitData(strcmp(hitData.Side,'Right'),:);

%% Identify the time lag between the Cirris  Xsens time series using cross-correlation
% Cirris signal
% Left
% sigA = speedZ_L_Cirris;
% sigA = speedY_L_Cirris;
% sigA = speedX_L_Cirris;

% Xsens signal
% sigB = velX_L_Xsens;
% sigB = velZ_L_Xsens;
% sigB = velY_L_Xsens;

% Right
% sigA = speedZ_R_Cirris;
% sigA = speedY_R_Cirris;
sigA = speedX_R_Cirris;

% Xsens signal
% sigB = velX_R_Xsens;
% sigB = velZ_R_Xsens;
sigB = velY_R_Xsens;


% Resample the Cirris signal to a fixed 60Hz rate (same as Xsens)
[sigA_res,~] = resample(sigA,timeL,60);

% Compute the cross-correlation
[r,lags] = xcorr(sigA_res,sigB);

% Indentify the lag (where cross-correlation is maximal)
[~,tidx_max]=max(abs(r));
t_max = lags(tidx_max)/60;

% Plot the cross-correlation and the indenitied lag
figure; plot(lags/60,r); hold on; plot(t_max,r(tidx_max),'r*','MarkerSize',10);
title('Cross correlation');

% Calculate the shift_t by which the Cirris time series has to shift to be
% aligned to the Xsens one
shift_t = timeL(1)+t_max;

%% Plotiing
if plot_flag
    % Plot speeds
    sigA = speedZ_R_Cirris;
    sigB = velX_R_Xsens;
    figure; plot(timeL-shift_t,(sigA-mean(sigA))./std(sigA),'b');hold on; plot(time,(sigB-mean(sigB,'omitnan'))./std(sigB,'omitnan'),'r');
    hold on; plot(delsys_hs_time,(delsys_hs_data-mean(delsys_hs_data,'omitnan'))./std(delsys_hs_data,'omitnan'));
    hold on; scatter(right_hits.TaskTimer-shift_t,0,'g','LineWidth',2);
    
    sigA = speedX_L_Cirris;
    sigB = velY_L_Xsens;
    figure; plot(timeL-shift_t,(sigA-mean(sigA))./std(sigA),'b');hold on; plot(time,(sigB-mean(sigB,'omitnan'))./std(sigB,'omitnan'),'r');
    hold on; scatter(hitData.TaskTimer-shift_t,0,'g','LineWidth',2);
    
    sigA = speedY_L_Cirris;
    sigB = velZ_L_Xsens;
    figure; plot(timeL-shift_t,(sigA-mean(sigA))./std(sigA),'b');hold on; plot(time,(sigB-mean(sigB,'omitnan'))./std(sigB,'omitnan'),'r');
    hold on; scatter(hitData.TaskTimer-shift_t,0,'g','LineWidth',2);
    
    % Plot accelerations
    sigA = accelZ_L_Cirris;
    sigB = accX_L_Xsens;
    figure; plot(timeL(1:end-1)-shift_t,(sigA-mean(sigA))./std(sigA),'b');hold on; plot(time,(sigB-mean(sigB,'omitnan'))./std(sigB,'omitnan'),'r');
    hold on; scatter(hitData.TaskTimer-shift_t,0,'g','LineWidth',2);
    
    sigA = accelX_L_Cirris;
    sigB = accY_L_Xsens;
    figure; plot(timeL(1:end-1)-shift_t,(sigA-mean(sigA))./std(sigA),'b');hold on; plot(time,(sigB-mean(sigB,'omitnan'))./std(sigB,'omitnan'),'r');
    hold on; scatter(hitData.TaskTimer-shift_t,0,'g','LineWidth',2);
    
    sigA = accelY_L_Cirris;
    sigB = accZ_L_Xsens;
    figure; plot(timeL(1:end-1)-shift_t,(sigA-mean(sigA))./std(sigA),'b');hold on; plot(time,(sigB-mean(sigB,'omitnan'))./std(sigB,'omitnan'),'r');
    hold on; scatter(hitData.TaskTimer-shift_t,0,'g','LineWidth',2);
end
end