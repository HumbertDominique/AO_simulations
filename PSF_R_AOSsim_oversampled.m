clc
clear all
close all
addpath('OOMAO')



ngs = source;

% TODO: put these into a txt file for input
r0 = 1.5e-3; %[m]
L0 = 30; % [m]
nAct = 11; % number of actuators across the pupil, including the ones outside the pupil
nL   = nAct*4;
nPx  = 27;

nRes = nL*nPx;
D    = 0.0195;
d    = D/nL; % lenslet pitch
samplingFreq = 500;

atm = atmosphere(photometry.HeNe,r0,L0,'fractionnalR0',[1],'altitude',[.5],'windSpeed',[.12],'windDirection',[pi]);

tel = telescope(D,'resolution',nRes,'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);

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
dm = deformableMirror(nAct,'modes',bifa, 'resolution',tel.resolution);
% dm = deformableMirror(nAct,'modes',bifa, 'resolution',tel.resolution, 'validActuator', wfs.validActuator);


calibDm = calibration(dm,wfs,ngs,ngs.wavelength/40);




wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

ngs = ngs.*tel;

tel = tel + atm;
% figure
% imagesc(tel)
ngs = ngs.*tel*wfs;

%% Diffraction limited performance
cam = imager();
instantCam = imager();

ngs = source('zenith',0,'azimuth',0,'magnitude',8);     % AO source
science = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',8);    % long psf source, Magniture is arbitrary
instantScience = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',8); %instantaneous psf source. could be the same as the long one.

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

%% Regulation settings

% TODO: put these into a txt file for input
cam.clockRate    = 1;
instantCam.clockRate    = 1;
exposureTime     = 100;
cam.exposureTime = exposureTime;
instantCam.exposureTime = 1;
startDelay       = 20;

gain_cl  = 0.5 % integrator gain

wfs.camera.photonNoise = false;
wfs.camera.readOutNoise = 0;


dm.coefs = zeros(dm.nValidActuator,1);
cam.startDelay = startDelay;
nIteration = startDelay + exposureTime;

% Performance history

% TODO: put these into a txt file for input
fileID_WF = 'ao_WF.h5';
fileID_WFS = 'ao_WFS.h5';
fileID_lightfield = 'ao_lightfield.h5';
fileID_DM = 'ao_DM.h5';
fileID_psf = 'ao_psf.h5';
fileID_rwfe = 'ao_rwfe.h5';


WFHistory = zeros(size(ngs.meanRmPhase,1),size(ngs.meanRmPhase,2),nIteration);
WFSHistory = zeros(length(wfs.slopes),nIteration);
lightfieldHistory = zeros(size(wfs.camera.frame, 1), size(wfs.camera.frame, 2),nIteration);
dmCommandsHistory = zeros(length(dm.coefs),nIteration);
psfHistory = zeros(size(instantCam.frame, 1), size(instantCam.frame, 2),nIteration);
rwfe_history = zeros(nIteration,1);
if exist(fileID_WF, 'file'), delete(fileID_WF); end
if exist(fileID_WFS, 'file'), delete(fileID_WFS); end
if exist(fileID_DM, 'file'), delete(fileID_DM); end
if exist(fileID_lightfield, 'file'), delete(fileID_lightfield); end
if exist(fileID_psf, 'file'), delete(fileID_psf); end
if exist(fileID_rwfe, 'file'), delete(fileID_rwfe); end


%% Regulation

flush(cam)
% the start delay could be implemented using 2 loops. the 1st is a startup to stabilise the regulator, and the 2nd is the main loop to collect data.
for k=1:nIteration
    % Objects update
    flush(instantCam)
    +tel;
    +ngs;
    +science;
    +instantScience;
    % Closed-loop controller
    % dm.coefs = dm.coefs - gain_cl*calibDm.M*wfs.slopes;
    dm.coefs = dm.coefs - gain_cl*dmCalib.M*wfs.slopes;
    dm.coefs = min(max(dm.coefs, -1), 1);
    % local log
    WFHistory(:,:,k) = ngs.meanRmPhase;
    WFSHistory(:,k) = wfs.slopes;
    lightfieldHistory(:,:,k) = wfs.camera.frame;
    dmCommandsHistory(:, k) = dm.coefs;
    rwfe_waves_history(k) = sqrt(var(ngs))./2/pi; % [waves]
    psfHistory(:,:,k) = instantCam.frame;
end

psf_sum = sum(psfHistory(:,:,startDelay+1:end), 3);

fprintf('Strehl ratio: %4.1f\n',cam.strehl);
fprintf('Strehl ratio: %4.1f\n',instantCam.strehl);

figure;
subplot(2,1,1);
plot(rwfe_waves_history, 'o-');
xlabel('Iteration'); ylabel('Residual WF RMS (waves)');
title('AO Loop Convergence (Linear Scale)');
grid on;
subplot(2,1,2);
semilogy(rwfe_waves_history, 'o-');
xlabel('Iteration'); ylabel('Residual WF RMS (waves)');
title('AO Loop Convergence (Logarithmic Scale)');
grid on;


figure;
imagesc(cam,'parent',subplot(2,2,1));
title('Long Exposure PSF', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;
imagesc(instantCam,'parent',subplot(2,2,2));
title('Instantaneous PSF', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;
imagesc(psf_sum,'parent',subplot(2,2,3));
title('Long psf from instantaneous PSF', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;
imagesc(cam.frame-psf_sum,'parent',subplot(2,2,4));
title('Long psf - iPsfSum', 'FontSize', 12, 'FontWeight', 'bold');
colorbar; axis image;
sgtitle(sprintf('AO Strehl: Long=%.2f, Instant=%.2f', cam.strehl, instantCam.strehl));


psf_sum_flux = sum(psf_sum(:));
long_psf_flux = sum(cam.frame(:));

fprintf('Flux in long exposure PSF: %.2e\n', long_psf_flux);
fprintf('Flux in sum of instantaneous PSFs: %.2e\n', psf_sum_flux);

flux_ratio = psf_sum_flux / long_psf_flux;
fprintf('Flux ratio (iPsfSum / Long PSF): %.3f\n', flux_ratio);


h5create(fileID_WF, '/wf', size(WFHistory));
h5write(fileID_WF, '/wf', WFHistory);
h5create(fileID_WFS, '/wfs', size(WFSHistory));
h5write(fileID_WFS, '/wfs', WFSHistory);
h5create(fileID_lightfield, '/wf_lightfield', size(lightfieldHistory));
h5write(fileID_lightfield, '/wf_lightfield', lightfieldHistory);
h5create(fileID_DM, '/dm_commands', size(dmCommandsHistory));
h5write(fileID_DM, '/dm_commands', dmCommandsHistory);
h5create(fileID_psf, '/psf_history', size(psfHistory));
h5write(fileID_psf, '/psf_history', psfHistory);
h5create(fileID_rwfe, '/rwfe_waves_history', size(rwfe_waves_history));
h5write(fileID_rwfe, '/rwfe_waves_history', rwfe_waves_history);