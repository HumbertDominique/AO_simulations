clc
clear all
close all
addpath('OOMAO')

fileID_WF = 'ao_WF.h5';
fileID_lightfield = 'ao_lightfield.h5';
fileID_WFS = 'ao_WFS.h5';
fileID_DM = 'ao_DM.h5';
fileID_psf = 'ao_psf.h5';


ngs =source('wavelength', photometry.HeNe);
atm = atmosphere(photometry.HeNe,15e-3,30,'altitude',5e3,'windSpeed',1,'windDirection',pi/3);

nL = 60;            % [-] Number of lenslet accross the pupil 1mmu lens/actuator
nPx = 10;           % [-] Px per lenslet
nRes = nL*nPx;      % [-] Resolution
D = 25;             % [-] Telescope diameter
d = D/nL;           % [-] Lenslet pitch
samplingFreq = 500; % [Hz] Sampling frequency

tel = telescope(D,'resolution', nRes, 'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);
wfs = shackHartmann(nL,nRes, 0.85);

% cam.referenceFrame = cam.frame;  % perfect PSF as ref


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

ngs = ngs.*tel; % propagation ut to, but not including DM ans wfs
dm.coefs = zeros(dm.nValidActuator,1);  % start flat
CalibDm = calibration(dm,wfs,ngs,ngs.wavelength); %this is not the cevtorised implementation.
% CalibDm  = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);
CalibDm.threshold = 1e6;
M = CalibDm.M;
disp(CalibDm)

%% regulation loop

nIter = 3;  % number of loop iterations
gain = 0.2;    % integrator gain


slopesHistory = zeros(nIter, length(wfs.slopes));
dmCommandsHistory = zeros(nIter, length(dm.coefs));
% lightfieldHistory = zeros(nIter, size(wfs.camera.frame));
% psfHistory = zeros(nIter, size(camH.frame));
if exist(fileID_WF, 'file'), delete(fileID_WF); end  % Fresh run
if exist(fileID_WFS, 'file'), delete(fileID_WFS); end  % Fresh run
if exist(fileID_DM, 'file'), delete(fileID_DM); end  % Fresh run
if exist(fileID_psf, 'file'), delete(fileID_psf); end  % Fresh run
if exist(fileID_lightfield, 'file'), delete(fileID_lightfield); end  % Fresh run
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

tel = tel + atm;
% tel = tel;
% figure
% imagesc(tel) 
tel = tel + atm;                      % bind atmosphere to telescope
ngs = source('wavelength', photometry.HeNe);
ngs = ngs.*tel*dm*wfs;    

camH = imager('diameter',25, 'fieldStopSize',30,'nyquistSampling',8*2);

scienceH = source('wavelength', photometry.HeNe);  % Science channel
scienceH = scienceH.*tel*dm*camH;  % Corrected path (tel includes dm)




for k = 1:nIter
    fprintf('Iteration %d\n', k)

    +tel;          % update atmosphere phase screen on telescope
    +ngs;          % propagate source through current optical path (tel*dm*wfs)
    +scienceH
    % Regulator
    dc = -gain * (M * wfs.slopes);  % DM command increment (minus sign for correction)
    dm.coefs = dm.coefs + dc;       % integrator controller
    % dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calibDm.M*wfs.slopes(:,1);
    % dm.coefs(:,2) = (1-gain_pol)*dm.coefs(:,2) + ...
% gain_pol*iF*( slmmse*( wfs.slopes(:,2) - calibDm.D*dm.coefs(:,2) ) );
    % +scienceH;       % propagate science source through corrected path (tel*dm*cam)
    slopesHistory(k,:) = wfs.slopes;
    lightfieldHistory(k,:) = wfs.camera;
    dmCommandsHistory(k,:) = dm.coefs;
    psfHistory(:,:,k) = camH.frame;

end

figure
imagesc(psfHistory(:,:,2));
% imagesc(wfs.camera);
% imagesc(det.frame);

% h5create(fileID_lightfield, '/wf_slopes', size(slopesHistory));
% h5write(fileID_lightfield, '/wf_slopes', slopesHistory);
% h5create(fileID_WF, '/wf_slopes', size(slopesHistory));
% h5write(fileID_WF, '/wf_slopes', slopesHistory);
% h5create(fileID_WFS, '/wfs_slopes', size(slopesHistory));
% h5write(fileID_WFS, '/wfs_slopes', slopesHistory);
% h5create(fileID_DM, '/dm_commands', size(dmCommandsHistory));
% h5write(fileID_DM, '/dm_commands', dmCommandsHistory);
% h5create(fileID_psf, '/psf_history', size(psfHistory));
% h5write(fileID_psf, '/psf_history', psfHistory);
