function play_sound(trial_type)
%%
% play_sound(...)
%
% Author: Kenneth Louie
% Project: General Gait Protocol
% Date: 05/20/2021
% Version: 1.0
%
% Description:
%   Plays a 50ms sound that has a different frequency depending on the
%   trial type
%
% Inputs:   trial_type  [=] integer (1-4)
%
% Outputs:  NONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent soundVector

if isempty(soundVector)
    amp = 5;                % amplitude of sound wave
    fs = 48000;             % sampling frequency
    duration = 0.05;        % in seconds (50ms currently)
    freq=250:250:1000;      % Different tones for different keys        
    values=0:1/fs:duration;     % time vector to produce noise
    soundVector=amp*sin( 2*pi.*freq'.*values);
end

switch trial_type
    case 0 % sound test
        sound(soundVector(1,:),48000,16);
    case 1 % cued
        sound(soundVector(2,:),48000,16);
    case 2 % self-initiated
        sound(soundVector(3,:),48000,16);
    case 3 % prepare-and-go
        sound(soundVector(4,:),48000,16);
end

end