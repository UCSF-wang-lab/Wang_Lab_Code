function calcRCSAlignment(varargin)
% Calculates the RCS-Delsys accelerometry alignment using cross-correlation
%
% Author:   Eleni Patelaki
% Date:     3/14/24

%% Read Delsys mat
data_str = 'DBS_5hz_Left_1_Acc_1Y_IM';
[delsys_fname, delsys_fpath] = uigetfile('*.mat','Select Delsys mat file');
delsys_mat = load(fullfile(delsys_fpath,delsys_fname));

delsys_time = getfield(delsys_mat.out_struct.Time,data_str); 
delsys_accel = getfield(delsys_mat.out_struct.Data,data_str);
% delsys_fs = delsys_mat.out_struct.srates(find(strcmp(delsys_mat.out_struct.Chan_names,data_str)==1));

%% Read accelelation RCS data
[rcs_fname,rcs_path] = uigetfile('RawDataAccel.mat');
assert(contains(rcs_path,'Left INS'));
rcs_mat = load(fullfile(rcs_path,rcs_fname));

rcs_time = rcs_mat.accelDataTable.DerivedTime/1000;
rcs_time = rcs_time - rcs_time(1);
rcs_accel = rcs_mat.accelDataTable.XSamples;

%% Resample both acceleration time series to a fixed 60-Hz rate
[rcs_accel_res,~] = resample(rcs_accel,rcs_time,60);
[delsys_accel_res,~] = resample(delsys_accel,delsys_time,60);

%% Compute cross-correlation
[r,lags] = xcorr(rcs_accel_res,delsys_accel_res);

%% Indentify the lag (where cross-correlation is maximal)
[~,tidx_max]=max(abs(r));
t_max = lags(tidx_max)/60;

%% Plot the cross-correlation and the indenitied lag
figure; plot(lags/60,r); hold on; plot(t_max,r(tidx_max),'r*','MarkerSize',10);
title('Cross correlation');

end
