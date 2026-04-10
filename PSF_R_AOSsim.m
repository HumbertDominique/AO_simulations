clc
clear all
close all
addpath('OOMAO')

cfg = readConfig('ao_inputs.txt');
r0           = cfg.r0;
L0           = cfg.L0;
Asl          = cfg.Asl;
wind         = cfg.wind;
windDir      = cfg.windDir;
nAct         = cfg.nAct;
oversampling  = cfg.oversampling;
if oversampling ~= 1
    nAct = nAct * oversampling;
else
    nL           = nAct - 1;
end
nPx          = cfg.nPx;
nRes         = nL * nPx;
D            = cfg.D;
d            = D / nL;
samplingFreq = cfg.samplingFreq;
chunksize    = cfg.chunksize;
exposureTime = cfg.exposureTime;
startDelay   = cfg.startDelay;
gain_cl      = cfg.gain_cl;
SAVEWF       = cfg.SAVEWF;
fileID_WF           = cfg.fileID_WF;
fileID_WFS          = cfg.fileID_WFS;
fileID_lightfield   = cfg.fileID_lightfield;
fileID_DM           = cfg.fileID_DM;
fileID_psf          = cfg.fileID_psf;
fileID_rwfe         = cfg.fileID_rwfe;
fileID_diff_limited = cfg.fileID_diff_limited;
fileID_metadata     = cfg.fileID_metadata;

ngs = source;

atm = atmosphere(photometry.HeNe,r0,L0,'fractionnalR0',[1],'altitude',Asl,'windSpeed',wind,'windDirection',windDir);
tel = telescope(D,'resolution',nRes,'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);

wfs = shackHartmann(nL,nRes,0.50);
wfs.camera.photonNoise  = cfg.photonNoise;
wfs.camera.readOutNoise = cfg.readOutNoise;

ngs = ngs.*tel*wfs;

wfs.INIT

+wfs;
figure
% imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
% slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))
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
% figure
% plot(Ox_in,Ox_out)
% grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);
ngs.zenith = 0;
wfs.pointingDirection = [];
%


bifa = influenceFunction('monotonic',0.75);
dm = deformableMirror(nL+2,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

ngs = ngs.*tel;
calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);

tel = tel + atm;
figure
imagesc(tel)
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

diff_limited = cam.frame;

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
cam.exposureTime = exposureTime;
instantCam.exposureTime = 1;

dm.coefs = zeros(dm.nValidActuator,1);
cam.startDelay = startDelay;
nIteration = startDelay + exposureTime;

% Performance history

nBatch = ceil(nIteration*(size(wfs.camera.frame, 1)* size(wfs.camera.frame, 2))/(chunksize));
batchItSize = floor(nIteration/(nIteration*(size(wfs.camera.frame, 1)* size(wfs.camera.frame, 2))/(chunksize)));
LastBatchItSize = nIteration - batchItSize*(nBatch-1);
fprintf('----------------------------------------------------------------\n');
fprintf('Total number of iterations to store: %d\n',nIteration);
fprintf('Qty of batches: %d\n',nBatch);
fprintf('Batch size [iterations]: %d\n', batchItSize);
fprintf('Last batch size [iterations]: %d\n', LastBatchItSize);
fprintf('----------------------------------------------------------------\n');


for i = 1 : batchItSize
    if SAVEWF
         if exist(fileID_WF+string(i)+".h5", 'file'), delete(fileID_WF+string(i)+".h5"); end
    end
    if exist(fileID_WFS+string(i)+".h5", 'file'), delete(fileID_WFS+string(i)+".h5"); end
    if exist(fileID_DM+string(i)+".h5", 'file'), delete(fileID_DM+string(i)+".h5"); end
    if exist(fileID_lightfield+string(i)+".h5", 'file'), delete(fileID_lightfield+string(i)+".h5"); end
    if exist(fileID_psf+string(i)+".h5", 'file'), delete(fileID_psf+string(i)+".h5"); end
    if exist(fileID_rwfe+string(i)+".h5", 'file'), delete(fileID_rwfe+string(i)+".h5"); end
    if exist(fileID_diff_limited+string(i)+".h5", 'file'), delete(fileID_diff_limited+string(i)+".h5"); end
end
if exist(fileID_metadata+".txt", 'file'), delete(fileID_metadata+".txt"); end

% WFHistory = zeros(size(ngs.meanRmPhase,1),size(ngs.meanRmPhase,2),batchItSize);
WFSHistory = zeros(length(wfs.slopes),batchItSize);
lightfieldHistory = zeros(size(wfs.camera.frame, 1), size(wfs.camera.frame, 2),batchItSize, 'uint8');
dmCommandsHistory = zeros(length(dm.coefs),batchItSize);
psfHistory = zeros(size(instantCam.frame, 1), size(instantCam.frame, 2),batchItSize, 'uint8');
rwfe_waves_history = zeros(batchItSize,1);


%% Regulation
% gain_cl = .9;
flush(cam)
flush(instantCam)

indexInBatch = 0;

for k=1:nIteration
    indexInBatch = indexInBatch + 1;
    % Objects update
    % flush(instantCam)
    instantCam.frame = [];
    instantCam.frameCount = 0;
    +tel;
    +ngs;
    +science;
    +instantScience;
    % Closed-loop controller
    dm.coefs = dm.coefs - gain_cl*calibDm.M*wfs.slopes;
    dm.coefs = min(max(dm.coefs, -1), 1);
    % local log
    WFHistory(:,:,indexInBatch) = ngs.meanRmPhase;
    WFSHistory(:,indexInBatch) = wfs.slopes;
    lightfieldHistory(:,:,indexInBatch) = uint8(floor(wfs.camera.frame*255));
    dmCommandsHistory(:, indexInBatch) = dm.coefs;
    rwfe_waves_history(indexInBatch) = sqrt(var(ngs))./2/pi; % [waves]
    psfHistory(:,:,indexInBatch) = uint8(floor(instantCam.frame*255));
    
    if mod(k-1, round(nIteration/50)) == 0 || k == nIteration
        fprintf('Progress: %d%% done\n', round(100*k/nIteration));
    end
    % fprintf('Current iteration: %d\n', mod(k, batchItSize))
    if mod(k, batchItSize) == 0 || k == nIteration
        fprintf('Saving batch %d/%d to disk...\n', ceil(k/batchItSize), nBatch);
        batchIndex = ceil(k/batchItSize);
        
        if SAVEWF
            sz = size(WFHistory);
            h5create(fileID_WF+string(batchIndex)+".h5", '/wf', sz, 'ChunkSize', [sz(1) sz(2) 1], 'DataType', 'double');
            h5write(fileID_WF+string(batchIndex)+".h5", '/wf', WFHistory);
        end
        sz = size(WFSHistory);
        h5create(fileID_WFS+string(batchIndex)+".h5", '/wfs', sz, 'ChunkSize', [sz(1) 1], 'DataType', 'double');
        h5write(fileID_WFS+string(batchIndex)+".h5", '/wfs', WFSHistory);
        WFSHistory = zeros(length(wfs.slopes),batchItSize);

        sz = size(lightfieldHistory);
        totBytes = prod(sz);   % for uint8, this is roughly the expected file size
        fprintf('Dataset size: [%d,%d,%d] -> %d bytes\n', sz, totBytes);
        h5create(fileID_lightfield+string(batchIndex)+".h5", '/wf_lightfield', sz, 'ChunkSize', [sz(1) sz(2) 1],'DataType', 'uint8');
        h5write(fileID_lightfield+string(batchIndex)+".h5", '/wf_lightfield', lightfieldHistory);
        lightfieldHistory = zeros(size(wfs.camera.frame, 1), size(wfs.camera.frame, 2),batchItSize);

        sz = size(dmCommandsHistory);
        h5create(fileID_DM+string(batchIndex)+".h5", '/dm_commands', sz, 'ChunkSize', [sz(1) 1], 'DataType', 'double');
        h5write(fileID_DM+string(batchIndex)+".h5", '/dm_commands', dmCommandsHistory);
        dmCommandsHistory = zeros(length(dm.coefs),batchItSize);

        sz = size(psfHistory);
        h5create(fileID_psf+string(batchIndex)+".h5", '/psf_history', sz, 'ChunkSize', [sz(1) sz(2) 1], 'DataType', 'uint8');
        h5write(fileID_psf+string(batchIndex)+".h5", '/psf_history', psfHistory);
        psfHistory = zeros(size(instantCam.frame, 1), size(instantCam.frame, 2),batchItSize);

        sz = size(rwfe_waves_history);
        h5create(fileID_rwfe+string(batchIndex)+".h5", '/rwfe_waves_history', sz, 'ChunkSize', [sz(1) 1], 'DataType', 'double');
        h5write(fileID_rwfe+string(batchIndex)+".h5", '/rwfe_waves_history', rwfe_waves_history);
        rwfe_waves_history = zeros(batchItSize,1);

        indexInBatch = 0;
    end
end
sz = size(diff_limited);
h5create(fileID_diff_limited+".h5", '/diff_limited', sz, 'ChunkSize', [sz(1) 1]);
h5write(fileID_diff_limited+".h5", '/diff_limited', diff_limited);
% h5create(fileID_WF, '/wf', size(WFHistory));
% h5write(fileID_WF, '/wf', WFHistory);
% h5create(fileID_WFS, '/wfs', size(WFSHistory));   
% h5write(fileID_WFS, '/wfs', WFSHistory);
% h5create(fileID_lightfield, '/wf_lightfield', size(lightfieldHistory));
% h5write(fileID_lightfield, '/wf_lightfield', lightfieldHistory);
% h5create(fileID_DM, '/dm_commands', size(dmCommandsHistory));
% h5write(fileID_DM, '/dm_commands', dmCommandsHistory);
% h5create(fileID_psf, '/psf_history', size(psfHistory));
% h5write(fileID_psf, '/psf_history', psfHistory);
% h5create(fileID_rwfe, '/rwfe_waves_history', size(rwfe_waves_history));
% h5write(fileID_rwfe, '/rwfe_waves_history', rwfe_waves_history);
% h5create(fileID_diff_limited, '/diff_limited', size(diff_limited));
% h5write(fileID_diff_limited, '/diff_limited', diff_limited);

rowNames = {'D';'r0';'L0';'Asl';'wind';'windDir';'Exposure time';'nIteration';'gain_cl';'batchItSize';'nBatch';'LastBatchItSize'};
values =    [D;  r0;  L0;  Asl;  wind;  windDir;  exposureTime;   nIteration;  gain_cl;  batchItSize;  nBatch;  LastBatchItSize];
T = table(values,'RowNames',rowNames);
writetable(T,fileID_metadata+".txt",'Delimiter','\t','WriteRowNames',true);

%% s
maxValue = max(psfHistory, [], 'all');
fprintf('Maximum value in frame 30: %f\n', maxValue);
% figure;
fprintf('max long psf value: %f\n', sum(cam.frame(:)))

figure;
imshow(psfHistory(:,:,21), []);
%%
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


clear WFHistory WFSHistory lightfieldHistory dmCommandsHistory psfHistory rwfe_history