function PSF_R_AOSsim(wind, i)
    close all
    addpath('OOMAO')

    % fixed parameters
    ngs = source;
    r0 = 3.75e-3; %[m]
    L0 = 30; % [m]
    Asl = [0.02]; % [m]
    windDir = [pi]; % [rad]
    nAct = 11; % number of actuators across the pupil, including the ones outside the pupil
    nL   = nAct-1;

    nPx  = 17;
    nRes = nL*nPx;
    D    = 0.0195;
    d    = D/nL; % lenslet pitch
    samplingFreq = 10;  %[Hz]



    fprintf('Running for wind = %.4f m/s\n', wind);
    atm = atmosphere(photometry.HeNe,r0,L0,'fractionnalR0',[1],'altitude',Asl,'windSpeed',wind,'windDirection',windDir, 'logging',false);
    tel = telescope(D,'resolution',nRes,'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);
    wfs = shackHartmann(nL,nRes,0.85);
    ngs = source;
    ngs = ngs.*tel*wfs;
    wfs.INIT
    +wfs;
    % avoid displaying every iteration for speed; comment out if you want visuals
    % figure
    % imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
    % slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))
    wfs.camera.frameListener.Enabled = false;
    wfs.slopesListener.Enabled = false;
    wfs.pointingDirection = zeros(2,1);
    pixelScale = ngs.wavelength/(2*d*wfs.lenslets.nyquistSampling);
    tipStep = pixelScale/2;
    nStep   = floor(nPx/3)*2;
    sx      = zeros(1,nStep+1);
    u       = 0:nStep;
    wfs.camera.frameListener.Enabled = false;
    wfs.slopesListener.Enabled = false;
    warning('off','oomao:shackHartmann:relay')
    for kStep=u
        ngs.zenith = -tipStep*kStep;
        +ngs;
        sx(kStep+1) = median(wfs.slopes(1:end/2));
    end
    Ox_in  = u*tipStep*constants.radian2arcsec;
    Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
    slopesLinCoef = polyfit(Ox_in,Ox_out,1);
    wfs.slopesUnits = 1/slopesLinCoef(1);
    ngs.zenith = 0;
    wfs.pointingDirection = [];
    bifa = influenceFunction('monotonic',0.75);
    dm = deformableMirror(nL+1,'modes',bifa,...
        'resolution',tel.resolution,...
        'validActuator',wfs.validActuator);
    wfs.camera.frameListener.Enabled = false;
    wfs.slopesListener.Enabled = false;
    ngs = ngs.*tel;
    calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2,'noshow', true);
    tel = tel + atm;
    ngs = ngs.*tel*wfs;
    %% Diffraction limited performance
    cam = imager();
    instantCam = imager();
    ngs = source('zenith',0,'azimuth',0,'magnitude',8);     % AO source
    science = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',8);    % long psf source
    instantScience = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',8); %instantaneous psf source
    tel = tel - atm;
    science = science.*tel*cam;
    +science
    cam.referenceFrame = cam.frame;
    +science
    instantScience = instantScience.*tel*instantCam;
    +instantScience
    instantCam.referenceFrame = instantCam.frame;
    +instantScience
    fprintf('Long PSF Strehl ratio ref: %4.1f\n',cam.strehl);
    fprintf('Ipsf Strehl ratio ref: %4.1f\n',instantCam.strehl);
    tel = tel + atm;
    +science
    +instantScience
    fprintf('Long PSF Strehl ratio single frame init: %4.1f\n',cam.strehl);
    fprintf('Ipsf Strehl ratio single frame init: %4.1f\n',instantCam.strehl);
    % Setting the the actual paths
    science = science.*tel*dm*cam;
    instantScience = instantScience.*tel*dm*instantCam;
    ngs = ngs.*tel*dm*wfs;

    fprintf('Long PSF Strehl ratio single frame init: %4.1f\n',cam.strehl);
    fprintf('Ipsf Strehl ratio single frame init: %4.1f\n',instantCam.strehl);
    % Setting the the actual paths

    science = science.*tel*dm*cam;
    instantScience = instantScience.*tel*dm*instantCam;
    ngs = ngs.*tel*dm*wfs;

    cam.clockRate    = 1;
    instantCam.clockRate    = 1;
    exposureTime     = 100;
    cam.exposureTime = exposureTime;
    instantCam.exposureTime = 1;
    startDelay       = 10;
    
    gain_cl  = 0.5 % integrator gain
    
    wfs.camera.photonNoise = false;
    wfs.camera.readOutNoise = 0;
    
    
    dm.coefs = zeros(dm.nValidActuator,1);
    cam.startDelay = startDelay;
    nIteration = startDelay + exposureTime;


    psfHistory = zeros(size(instantCam.frame, 1), size(instantCam.frame, 2),nIteration);
    rwfe_waves_history = zeros(nIteration,1);

    cam.frameListener.Enabled = false;
    instantCam.frameListener.Enabled = false;
    flush(cam)
    for k=1:nIteration
        % Objects update
        flush(instantCam)
        +tel;
        +ngs;
        +science;
        +instantScience;
        % Closed-loop controller
        dm.coefs = dm.coefs - gain_cl*calibDm.M*wfs.slopes;
        dm.coefs = min(max(dm.coefs, -1), 1);
        % record
        rwfe_waves_history(k) = sqrt(var(ngs))./2/pi; % [waves]
        psfHistory(:,:,k) = instantCam.frame;
    end

    psf_sum = sum(psfHistory(:,:,startDelay+1:end), 3);

    % psf_sum_flux = sum(psf_sum(:));
    % long_psf_flux = sum(cam.frame(:));

    strehl_long = cam.strehl;
    % strehl_inst = instantCam.strehl;

    if i == 1
        fid = fopen('strehl_long.csv','w');
        if fid == -1
            error('Could not open strehl_long.csv for writing.');
        else
        fprintf(fid, 'wind,strehl_long\n'); % Write header
        end
    end
    fid = fopen('strehl_long.csv','a');
    if fid == -1
        error('Could not open strehl_long.csv for appending.');
    end
    fprintf(fid, '%.6f,%.6f\n', wind, strehl_long);
    fclose(fid);
    


    if i == 1
        fid = fopen('rwfe_waves_history.csv','w');
        if fid == -1
            error('Could not open rwfe_waves_history.csv for writing.');
        else

        fprintf(fid, '\n');
        end
    end

    fid = fopen('rwfe_waves_history.csv','a');
    if fid == -1
        error('Could not open rwfe_waves_history.csv for appending.');
    end
    fprintf(fid, '%.6f', wind);
    for k = 1:numel(rwfe_waves_history)
        fprintf(fid, ',%.6f', rwfe_waves_history);
    end
    fprintf(fid, '\n');
    fclose(fid);
    
clear all
end    