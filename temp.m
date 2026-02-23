clc
clear all
close all
addpath('OOMAO')

ngs =source('wavelength', photometry.HeNe);
atm = atmosphere(photometry.HeNe,15e-2,30,'altitude',5e3,'windSpeed',10,'windDirection',pi/3);
cam = imager();

nL = 60;            % [-] Number of lenslet accross the pupil 1mmu lens/actuator
nPx = 10;           % [-] Px per lenslet
nRes = nL*nPx;      % [-] Resolution
D = 25;             % [-] Telescope diameter
d = D/nL;           % [-] Lenslet pitch
samplingFreq = 500; % [Hz] Sampling frequency

tel = telescope(D,'resolution', nRes, 'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);
wfs = shackHartmann(nL,nRes, 0.85);

cam.referenceFrame = cam.frame;  % perfect PSF as ref

ngs = ngs.*tel*wfs;

wfs.INIT;
+wfs;
% figure
% imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
% slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

wfs.pointingDirection = zeros(2,1);

pixelScale = ngs.wavelength/(2*d*wfs.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep = floor(nPx/3)*2;
sx = zeros(1,nStep+1);
u = 0:nStep;

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
warning('off', 'oomao:shackHartmann:relay')

% moves the natural guide star in the field of view and record the median slope value for each step
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end

warning('on', 'oomao:shackHartmann:relay')

Ox_in = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;

% figure
% plot(Ox_in, Ox_out)
% grid on
slopeLinCoef = polyfit(Ox_in, Ox_out, 1);
wfs.slopesUnits = 1/slopeLinCoef(1); %This is the slope unit conversion factor

% resets the star position
ngs.zenith = 0;
% has the wfs always aligned with the star
wfs.pointingDirection = [];

% DM, a custom influence function can be used
bifa = influenceFunction('monotonic', 0.75);
dm = deformableMirror(nL+1, 'modes', bifa, 'resolution', tel.resolution, 'validActuator', wfs.validActuator);


%interaciton matrix
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

ngs.*tel; % propagation ut to, but not including DM ans wfs
dm.coefs = zeros(dm.nValidActuator,1);  % start flat
CalibDm = calibration(dm,wfs,ngs,ngs.wavelength); %this is not the cevtorised implementation.
CalibDm.threshold = 1e6;
disp(CalibDm)


tel = tel + atm;
figure
imagesc(tel)


%% regulation loop
nIter = 10;  % number of loop iterations


slopesHistory = zeros(nIter, length(wfs.slopes));
dmCommandsHistory = zeros(nIter, length(dm.coefs));
fileID = 'ao_data.h5';
if exist(fileID, 'file'), delete(fileID); end  % Fresh run

ngs = ngs.*tel;
ngs = ngs.*tel*dm*wfs;       


wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
gain = 0.2;    % integrator gain

M = CalibDm.M;

science = source('wavelength', photometry.HeNe);  % Science channel
% cam = imager();
science = science.*tel*cam;  % Corrected path (tel includes dm)

for k = 1:nIter
    fprintf('Iteration %d\n', k)
    +tel;          % update atmosphere phase screen on telescope
    +ngs;          % propagate source through current optical path (tel*dm*wfs)
    dc = -gain * (M * wfs.slopes);  % DM command increment (minus sign for correction)
    dm.coefs = dm.coefs + dc;       % integrator controller
    +science;       % propagate science source through corrected path (tel*dm*cam)
    slopesHistory(k,:) = wfs.slopes;
    dmCommandsHistory(k,:) = dm.coefs;
    psfHistory(:,:,k) = cam.frame;

end

figure
imagesc(psfHistory(:,:,5));


h5create(fileID, '/wf_slopes', size(slopesHistory));
h5write(fileID, '/wf_slopes', slopesHistory);
h5create(fileID, '/dm_commands', size(dmCommandsHistory));
h5write(fileID, '/dm_commands', dmCommandsHistory);
h5create(fileID, '/psf_history', size(psfHistory));
h5write(fileID, '/psf_history', psfHistory);
