clc
clear all
close all
% gain = 0.5;
% wind = (0.001:0.005:0.026); % [m/s]
gain = 0.9;
wind = (0.031:0.005:0.056); % [m/s]
% wind = (0.0145:0.0005:0.0165); % [m/s]
% wind = (0.0170:0.0005:0.0190); % [m/s]
% gain = 0.8;
% wind = (0.0195:0.0005:0.0215); % [m/s]
% gain = 0.9;
% wind = (0.0245:0.0005:0.0260); % [m/s]
% wind = (0.027:0.001:0.074); % [m/s]


n = length(wind);
% wind = [0.001, 0.004]
fprintf('Number of wind speeds: %d\n', n);
for i = 1:n
    fprintf('Running simulation for wind speed: %.6f m/s, iteration: %d/%d\n', wind(i), i, n);
    PSF_R_AOSsim_spec(wind(i), gain, i+1) % +1 to avoid overwriting the first values
end