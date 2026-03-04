clc
clear all
close all
addpath('OOMAO')

ngs =source;
r0 = 1.5e-3; %[m]
L0 = 30; % [m]
atm = atmosphere(photometry.HeNe,r0,L0,'fractionnalR0',[1],'altitude',[.5],'windSpeed',[.5],'windDirection',[pi]);

nL   = 10;
nPx  = 17;
nRes = nL*nPx;
D    = 0.015;
d    = D/nL; % lenslet pitch
samplingFreq = 500;

tel = telescope(D,'resolution',nRes,...
    'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);

wfs = shackHartmann(nL,nRes,0.85);

ngs = ngs.*tel*wfs;

wfs.INIT

+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
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
    % drawnow
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
%figure
%plot(Ox_in,Ox_out)
%grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);
ngs.zenith = 0;
wfs.pointingDirection = [];
%


bifa = influenceFunction('monotonic',0.75);
dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

ngs = ngs.*tel;
calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);

tel = tel + atm;
% figure
% imagesc(tel)
ngs = ngs.*tel*wfs;

%% Diffraction limited performance
cam = imager();
instantCam = imager();

science = source('wavelength',photometry.HeNe);

tel = tel - atm;
science = science.*tel*cam;
% figure(31416)
% imagesc(cam,'parent',subplot(2,1,1))


cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)
% Atmospheric turbulence performance

tel = tel + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)
% Regulation ?

% ngsCombo = source('zenith',zeros(1,2),'azimuth',zeros(1,2),'magnitude',8);
ngsCombo = source('zenith',0,'azimuth',0,'magnitude',8);

ngsCombo = ngsCombo.*tel*dm*wfs;
% scienceCombo = source('zenith',zeros(1,2),'azimuth',zeros(1,2),'wavelength',photometry.HeNe,'magnitude',1);
scienceCombo = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',8);
instantScience = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',8);

scienceCombo = scienceCombo.*tel*dm*cam;
instantScience = instantScience.*tel*dm*instantCam;

%%

flush(cam)
cam.clockRate    = 1;
exposureTime     = 100;
cam.exposureTime = exposureTime;
startDelay       = 20;

instantCam.clockRate    = 1;
instantCam.exposureTime = 1;

figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
% subplot(2,1,2)
% h = imagesc(catMeanRmPhase(scienceCombo));
% axis xy equal tight
% colorbar


%% Regulation
gain_cl  = 0.5;
% dm.coefs = zeros(dm.nValidActuator,2);
dm.coefs = zeros(dm.nValidActuator,1);

% set(scienceCombo, 'logging', true);  
% set(scienceCombo, 'phaseVar', []);  
flush(cam)

cam.clockRate    = 1;
instantCam.clockRate    = 1;
exposureTime     = 1000;
cam.exposureTime = exposureTime;
instantCam.exposureTime = 1;
startDelay       = 20;


cam.startDelay = startDelay;
psf_short = zeros(size(instantCam.frame,1), size(instantCam.frame,2), exposureTime);
nIteration = startDelay + exposureTime;

% the start delay could be implemented using 2 loops. the 1st is a startup to stabilise the regulator, and the 2nd is the main loop to collect data.
for k=1:nIteration
    % Objects update
    flush(instantCam)
    +tel;
    +ngsCombo;
    +scienceCombo;
    +instantScience;
    % Closed-loop controller
    dm.coefs = dm.coefs - gain_cl*calibDm.M*wfs.slopes;
    % rwfe_history(k) = scienceCombo.meanRmOpd
    % rwfe_history(k) = ngsCombo.meanRmOpd;  % 

end
psf_sum = sum(psf_short(:,:,startDelay+1:end), 3);
% psf_sum = sum(psf_short, 3);

% size_variable = size(scienceCombo.meanRmO);
% fprintf('Size of variable: [%s]\n', mat2str(size_variable));
% fprintf('Strehl ratio: %4.1f\n',cam.strehl)
% fprintf('Strehl ratio: %4.1f\n',instantCam.strehl);
figure;
imagesc(cam,'parent',subplot(3,1,1));
title('Long Exposure PSF', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;
imagesc(instantCam,'parent',subplot(3,1,2));
title('Instantaneous PSF', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;
imagesc(psf_sum,'parent',subplot(3,1,3));
title('Long psf from instantaneous PSF', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;

sgtitle(sprintf('AO Strehl: Long=%.2f, Instant=%.2f', cam.strehl, instantCam.strehl));


% figure;
% semilogy(phi, 'o-');
% xlabel('Iteration'); ylabel('Residual WF RMS (waves)');
% title('AO Loop Convergence');
% grid on;


% phaseEstResRms = ngs.opdRms;  
% fprintf('Residual wavefront error: %4.2fnm\n', 1e9*phaseEstResRms/ngs.waveNumber)


% telLowRes = telescope(tel.D,'resolution',nL+1,'fieldOfViewInArcsec',30,'samplingTime',1/500);
% telLowRes.pupil = wfs.validActuator;

% telLowRes= telLowRes + atm;
% ngs = ngs.*telLowRes;
% phase = ngs.meanRmOpd;

% phaseEst = tools.meanSub( wfs.finiteDifferenceWavefront*ngs.wavelength ,...
% wfs.validActuator);

% ngs = ngs.*telLowRes*{wfs.validActuator,-2*pi*wfs.finiteDifferenceWavefront};
% phaseEstRes = ngs.meanRmOpd;
% phaseEstResRms = ngs.opdRms;
% % phaseEst = wfs.finiteDifferenceWavefront;  
% % % Compute residual by propagating through estimated phase  
% % ngs = ngs.*telLowRes*{wfs.validActuator, -2*pi*phaseEst};  
% % residualRms = ngs.opdRms;

% fprintf(phaseEstResRms)