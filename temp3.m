clc
clear all
close all
addpath('OOMAO')

ngs =source;
r0 = 50e-2; %[m]
L0 = 30; % [m]
atm = atmosphere(photometry.HeNe,r0,L0,'fractionnalR0',[1],'altitude',[12e3],'windSpeed',[20],'windDirection',[pi]);

nL   = 10;
nPx  = 17;
nRes = nL*nPx;
D    = 10;
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
% scienceCombo = source('zenith',zeros(1,2),'azimuth',zeros(1,2),'wavelength',photometry.HeNe);
scienceCombo = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',1);

scienceCombo = scienceCombo.*tel*dm*cam;
%%

flush(cam)
cam.clockRate    = 1;
exposureTime     = 100;
cam.exposureTime = exposureTime;
startDelay       = 1;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(scienceCombo));
axis xy equal tight
colorbar


%%

gain_cl  = 0.3;
dm.coefs = zeros(dm.nValidActuator,2);
flush(cam)
cam.clockRate    = 1;
exposureTime     = 200;
cam.exposureTime = exposureTime;
startDelay       = 200;
nIteration = startDelay + exposureTime;
for k=1:nIteration
    % Objects update
    +tel;
    +ngsCombo;
    +scienceCombo;
    % Closed-loop controller
    dm.coefs = dm.coefs - gain_cl*calibDm.M*wfs.slopes;    
end
figure;
imagesc(cam)
fprintf('Strehl ratio: %4.1f\n',cam.strehl)
