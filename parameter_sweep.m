clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainScript = 'PSF_R_AOSsim_oversampled.m';

% Parameter to sweep
ToSweep = [0, 2, 4, 6, 8]; % NGS magnitude values to sweep through
sweptParameter = 'magnitude'; % automatic replacement of the correct value not implemented

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% helper function to have variable separated from this script.
function run_one_case(mainScript)
    run(mainScript);
end
%%%



nRuns = numel(ToSweep);
h = waitbar(0, 'Starting parameter sweep...');
pos = get(h, 'Position');      % [left bottom width height]
pos(4) = 120;                  % make it a bit taller
set(h, 'Position', pos);

t0 = tic;
for k = 1: nRuns
    fprintf('Running simulation with %s = %g\n', sweptParameter, ToSweep(k));
    
    sweptValue = ToSweep(k);
    
    if k == 1
        msg = sprintf(['Run %d/%d\n' ...
                   'Value = %g\n' ...
                   'Runs left = %d\n' ...
                   'Est. time left = too early to say...'], ...
                   k, nRuns, sweptValue, nRuns-k+1);
        waitbar((k-1)/nRuns, h, msg);
    else
        msg = sprintf(['Finished run %d/%d\n' ...
                   'Runs left = %d\n' ...
                   'Est. time left = %.1f s'], ...
                   k, nRuns, nRuns-k, remainingTime);
        if ~isgraphics(h)
            h = waitbar((k-1)/nRuns, 'Starting parameter sweep...');
        else
            waitbar((k-1)/nRuns, h, msg);
        end
    end

    runTag{k} = sprintf('%s_%0.4g', sweptParameter, sweptValue);
    runTag{k} = strrep(runTag{k}, '.', 'p');

    outDir = fullfile('output', runTag{k});
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fid = fopen('ao_inputs.txt', 'w');

    fprintf(fid, '# AO Simulation Parameters\n');
    fprintf(fid, '# Units are indicated in comments; edit values only\n\n');

    
    fprintf(fid, 'SEED = 42                       # integer or None\n\n');

    fprintf(fid, '# Atmosphere\n');
    fprintf(fid, 'r0           = 3.75e-3          # [m] (double) Fried parameter\n');
    fprintf(fid, 'L0           = 30               # [m] (double) Outer scale\n');
    fprintf(fid, 'Asl          = [0.02]           # [m] (array of double) Layer altitude(s)\n');
    fprintf(fid, 'wind         = [0.04]           # [m/s] (array of double) Wind speed(s)\n');
    fprintf(fid, 'windDir      = [355/113]        # [rad] (array of double) Wind direction(s)\n\n');

    fprintf(fid, 'Sensors\n');
    fprintf(fid, 'sensor_type = double             # [format] Sensor resolution: uint8, uint16, single, double. Doubles are required for simulations with noise. The value will be overwritten. There is no maping of the values for an equivalent full well capacity. The magnitude or exposure time has to be adjusted accordingly. Bothe the PSF and SHWFS sensor must be of the same type.\n\n');

    fprintf(fid, '# NGS\n');
    fprintf(fid, 'NGSmagnitude = %g               # [m] (double) Natural guide star magnitude\n\n', sweptValue);

    fprintf(fid, '# DM / WFS geometry\n');
    fprintf(fid, 'oversampling = 1                # [-] (uint)\n');
    fprintf(fid, 'nActWSF      = 10               # [-] (uint)\n');
    fprintf(fid, 'edge_act     = 0.5              # [-] (uint)\n');
    fprintf(fid, 'nPx          = 27               # [-] (uint)\n');
    fprintf(fid, 'D            = 0.0175           # [m] (double)\n');
    fprintf(fid, 'dmStroke     = 10e-9            # [m] (double)\n\n');

    fprintf(fid, '# WFS camera\n');
    fprintf(fid, 'SH_ill_thresh = 0.5             # [-] (double)\n');
    fprintf(fid, 'photonNoise  = True             # [-] (bool)\n');
    fprintf(fid, 'readOutNoise = 13               # [e-/px] (double)\n\n');

    fprintf(fid, '# Timing\n');
    fprintf(fid, 'samplingFreq = 500              # [Hz] (double)\n');
    fprintf(fid, 'lag_c          = 1              # [cycles] lag = lag_c/samplingFreq [s]\n\n');

    fprintf(fid, '# Data storage\n');
    fprintf(fid, 'chunksize    = 200e6             # [B] (double)\n');
    fprintf(fid, 'exposureTime = 2000            # [iterations]\n');
    fprintf(fid, 'startDelay   = 10              # [iterations]\n');
    fprintf(fid, 'gain_cl      = 0.5              # [-] Integrator gain\n\n');

    fprintf(fid, '# Log\n');
    fprintf(fid, 'SAVEWF          = false\n');
    fprintf(fid, 'SAVESLOPES      = false\n');
    fprintf(fid, 'SAVELIGHTFIELD  = false\n');
    fprintf(fid, 'SAVEDM          = false\n');
    fprintf(fid, 'SAVEPSF         = true\n');
    fprintf(fid, 'SAVERWFE        = true\n');
    fprintf(fid, 'SAVEDIFFLIMITED = true\n');
    fprintf(fid, 'SAVEINSTANTDIFFLIMITED = true\n\n');


    fprintf(fid, '# Output file prefixes\n');
    fprintf(fid, 'outputDir           = %s\n', outDir);
    fprintf(fid, 'fileID_WF           = ao_WF_\n');
    fprintf(fid, 'fileID_WFS          = ao_WFS_\n');
    fprintf(fid, 'fileID_lightfield   = ao_lightfield_\n');
    fprintf(fid, 'fileID_DM           = ao_DM_\n');
    fprintf(fid, 'fileID_ipsf          = ao_ipsf_\n');
    fprintf(fid, 'fileID_ipsf_diff_lim= ao_ipsf_difflim_\n');
    fprintf(fid, 'fileID_rwfe         = ao_rwfe_\n');
    fprintf(fid, 'fileID_diff_limited = ao_diff_limited_\n');
    fprintf(fid, 'fileID_metadata     = ao_metadata\n');

    fclose(fid);

    run_one_case(mainScript);
    elapsed = toc(t0);
    avgTime = elapsed / k;
    remainingTime = avgTime * (nRuns - k);

    msg = sprintf(['Finished run %d/%d\n' ...
                   'Runs left = %d\n' ...
                   'Est. time left = %.1f s'], ...
                   k, nRuns, nRuns-k, remainingTime);

    waitbar(k/nRuns, h, msg);
end
msg = sprintf('Finished\n\n');
waitbar(1, h, msg);
