clc
clear all
close all

wind = (0.0095:0.0005:0.0115); % [m/s]
n = length(wind);
% wind = [0.001, 0.004]
fprintf('Number of wind speeds: %d\n', n);
for i = 1:n
    fprintf('Running simulation for wind speed: %.6f m/s, iteration: %d/%d\n', wind(i), i, n);
    PSF_R_AOSsim_spec(wind(i), i+1) % +1 to avoid overwriting the first values
end