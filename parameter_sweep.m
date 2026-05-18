% helper funcion to have variavle separeded from this script.
function run_one_case(mainScript)
    run(mainScript);
end

mainScript = 'PSF_R_AOSsim_oversampled.m';


% Parameter to sweep
ToSweep = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]; % LGS magnitude values to sweep through
sweptParameter = 'magnitude';

for k = 1:numel(ToSweep)
    fprintf('Running simulation with %s = %g\n', sweptParameter, ToSweep(k));
    sweptValue = ToSweep(k);

    runTag{k} = sprintf('%s_%0.4g', sweptParameter, sweptValue);
    runTag{k} = strrep(runTag{k}, '.', 'p');

    outDir = fullfile('output', runTag{k});
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fid = fopen('ao_inputs.txt', 'w');

    fprintf(fid, '# AO Simulation Parameters\n');
    fprintf(fid, '# Units are indicated in comments; edit values only\n\n');

    fprintf(fid, '# Atmosphere\n');
    fprintf(fid, 'r0           = 3.75e-3          # [m] (double) Fried parameter\n');
    fprintf(fid, 'L0           = 30               # [m] (double) Outer scale\n');
    fprintf(fid, 'Asl          = [0.02]           # [m] (array of double) Layer altitude(s)\n');
    fprintf(fid, 'wind         = [0.04]           # [m/s] (array of double) Wind speed(s)\n');
    fprintf(fid, 'windDir      = [355/113]        # [rad] (array of double) Wind direction(s)\n\n');

    fprintf(fid, '# LGS\n');
    fprintf(fid, 'NGSmagnitude = %g                # [m] (double) Laser guide star magnitude\n\n', sweptValue);

    fprintf(fid, '# DM / WFS geometry\n');
    fprintf(fid, 'oversampling = 1                # [-] (uint)\n');
    fprintf(fid, 'nActWSF      = 10               # [-] (uint)\n');
    fprintf(fid, 'edge_act     = 0.5              # [-] (uint)\n');
    fprintf(fid, 'nPx          = 27               # [-] (uint)\n');
    fprintf(fid, 'D            = 0.0175           # [m] (double)\n');
    fprintf(fid, 'dmStroke     = 10e-9            # [m] (double)\n\n');

    fprintf(fid, '# WFS camera\n');
    fprintf(fid, 'SH_ill_thresh = 0.005           # [-] (double)\n');
    fprintf(fid, 'photonNoise  = false            # [-] (bool)\n');
    fprintf(fid, 'readOutNoise = 0                # [-] (double)\n\n');

    fprintf(fid, '# Timing\n');
    fprintf(fid, 'samplingFreq = 500              # [Hz] (double)\n\n');

    fprintf(fid, '# Data storage\n');
    fprintf(fid, 'chunksize    = 10e6           # [B] (double)\n');
    fprintf(fid, 'exposureTime = 10000            # [iterations]\n');
    fprintf(fid, 'startDelay   = 100              # [iterations]\n');
    fprintf(fid, 'gain_cl      = 0.5              # [-] Integrator gain\n\n');

    fprintf(fid, '# Log\n');
    fprintf(fid, 'SAVEWF          = false\n');
    fprintf(fid, 'SAVESLOPES      = false\n');
    fprintf(fid, 'SAVELIGHTFIELD  = false\n');
    fprintf(fid, 'SAVEDM          = false\n');
    fprintf(fid, 'SAVEPSF         = false\n');
    fprintf(fid, 'SAVERWFE        = false\n');
    fprintf(fid, 'SAVEDIFFLIMITED = false\n\n');

    fprintf(fid, '# Output file prefixes\n');
    fprintf(fid, 'outputDir           = %s\n', outDir);
    fprintf(fid, 'fileID_WF           = ao_WF_\n');
    fprintf(fid, 'fileID_WFS          = ao_WFS_\n');
    fprintf(fid, 'fileID_lightfield   = ao_lightfield_\n');
    fprintf(fid, 'fileID_DM           = ao_DM_\n');
    fprintf(fid, 'fileID_psf          = ao_psf_\n');
    fprintf(fid, 'fileID_rwfe         = ao_rwfe_\n');
    fprintf(fid, 'fileID_diff_limited = ao_diff_limited_\n');
    fprintf(fid, 'fileID_metadata     = ao_metadata\n');

    fclose(fid);

    run_one_case(mainScript);
end
