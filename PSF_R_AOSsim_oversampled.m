% clc
clear all
close all
addpath('OOMAO')
input_file = "ao_inputs.txt"

cfg = readConfig('ao_inputs.txt');
r0           = cfg.r0;
L0           = cfg.L0;
Asl          = cfg.Asl;
wind         = cfg.wind;
windDir      = cfg.windDir;
nActWSF         = cfg.nActWSF;
oversampling  = cfg.oversampling;
edge_act     = cfg.edge_act;

NGSmagnitude = cfg.NGSmagnitude;

nL = nActWSF * oversampling;

nPx          = cfg.nPx;
nRes         = nL * nPx;
D            = cfg.D;
d            = D / nL;
dmStroke     = cfg.dmStroke;

samplingFreq = cfg.samplingFreq;
lag_c        = cfg.lag_c;

chunksize    = cfg.chunksize;
exposureTime = cfg.exposureTime;
startDelay   = cfg.startDelay;
gain_cl      = cfg.gain_cl;
SH_ill_thresh = cfg.SH_ill_thresh;
photonNoise = cfg.photonNoise;
readOutNoise = cfg.readOutNoise;

if readOutNoise == 1
    sensor_type = 'double'
else
    sensor_type = cfg.sensor_type;
end

SAVEWF       = cfg.SAVEWF;
SAVESLOPES      = cfg.SAVESLOPES;
SAVELIGHTFIELD  = cfg.SAVELIGHTFIELD;
SAVEDM          = cfg.SAVEDM;
SAVEPSF         = cfg.SAVEPSF;
SAVERWFE        = cfg.SAVERWFE;
SAVEDIFFLIMITED = cfg.SAVEDIFFLIMITED;
SAVEINSTANTDIFFLIMITED = cfg.SAVEINSTANTDIFFLIMITED;

outputDir           = cfg.outputDir;
fileID_WF           = cfg.fileID_WF;
fileID_WFS          = cfg.fileID_WFS;
fileID_lightfield   = cfg.fileID_lightfield;
fileID_DM           = cfg.fileID_DM;
fileID_ipsf          = cfg.fileID_ipsf;
fileID_ipsf_diff_lim= cfg.fileID_ipsf_diff_lim;
fileID_rwfe         = cfg.fileID_rwfe;
fileID_diff_limited = cfg.fileID_diff_limited;
fileID_metadata     = cfg.fileID_metadata;

metadataFile = outputDir + "/metadata.txt";

%% code

ngs = source;

atm = atmosphere(photometry.HeNe,r0,L0,'fractionnalR0',[1],'altitude',Asl,'windSpeed',wind,'windDirection',windDir);
tel = telescope(D,'resolution',nRes,'fieldOfViewInArcsec',3,'samplingTime',1/samplingFreq);

wfs = shackHartmann(nL,nRes,SH_ill_thresh);

ngs = ngs.*tel*wfs;
wfs.INIT

+wfs;
% figure
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

%figure
%plot(Ox_in,Ox_out)
%grid

slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);
ngs.zenith = 0;
wfs.pointingDirection = [];


bifa = influenceFunction('monotonic',0.47);
act_tot = nActWSF + 2*edge_act;
% Create a circular mask for the DM actuators to only allow the pupil+1 actuators to be active.
[x, y] = meshgrid(1:act_tot, 1:act_tot);
c = (act_tot + 1) / 2;
r = act_tot / 2;
DM_MASK = (x - c).^2 + (y - c).^2 <= r^2;

dm = deformableMirror(act_tot,'modes',bifa, 'resolution',tel.resolution, 'validActuator', DM_MASK); % valid actuators is used to ensure proper calibration matrix
calibDm = calibration(dm,wfs,ngs,dmStroke);

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

wfs.camera.photonNoise = photonNoise;
wfs.camera.readOutNoise = readOutNoise;

ngs = ngs.*tel;

tel = tel + atm;
% figure
% imagesc(tel)

ngs = ngs.*tel*wfs;

%% temp

% figure;
% imagesc(dm.validActuator);
% axis square tight;
% title('validActuator: true=inside pupil, false=outside pupil');
% colorbar;

%% Diffraction limited performance
cam = imager();
instantCam = imager();

camera.photonNoise = photonNoise;
cam.readOutNoise = readOutNoise;

instantCam.photonNoise = photonNoise;
instantCam.readOutNoise = readOutNoise;


ngs = source('zenith',0,'azimuth',0,'magnitude',NGSmagnitude);     % AO source
ngs.log.verbose = false;
science = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',NGSmagnitude);    % long psf source, Magniture is arbitrary
science.log.verbose = false;
instantScience = source('zenith',0,'azimuth',0,'wavelength',photometry.HeNe,'magnitude',NGSmagnitude); % instantaneous psf source. could be the same as the long one.
instantScience.log.verbose = false;
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

% fprintf('Long PSF Strehl ratio single frame init: %4.1f\n',cam.strehl);
% fprintf('Ipsf Strehl ratio single frame init: %4.1f\n',instantCam.strehl);

% Setting the the actual paths

science = science.*tel*dm*cam;
if SAVEINSTANTDIFFLIMITED
    if exist(outputDir+"\ipsf_difflim.h5", 'file'), delete(outputDir+"\ipsf_difflim.h5"); end
    tel = tel-atm
    instantScience = instantScience.*tel*dm*instantCam;
    +instantScience
    iPSF_strehl = type_cast(instantCam.frame,sensor_type);
    sz = size(iPSF_strehl);
            h5create(outputDir+"\ipsf_difflim.h5", '/ipsf_difflim', sz, 'ChunkSize', [sz(1) sz(2)], 'DataType', sensor_type);
            h5write(outputDir+"\ipsf_difflim.h5", '/ipsf_difflim', iPSF_strehl);
    clear iPSF_strehl
    tel = tel+atm
end

instantScience = instantScience.*tel*dm*instantCam;

ngs = ngs.*tel*dm*wfs;

%% Regulation settings

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

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for i = 1 : batchItSize
    if SAVEWF
         if exist(outputDir+"\"+fileID_WF+string(i)+".h5", 'file'), delete(outputDir+"\"+fileID_WF+string(i)+".h5"); end
        WFHistory = zeros(size(ngs.meanRmPhase,1),size(ngs.meanRmPhase,2),batchItSize);
    end
    if SAVESLOPES
        if exist(outputDir+"\"+fileID_WFS+string(i)+".h5", 'file'), delete(outputDir+"\"+fileID_WFS+string(i)+".h5"); end
        WFSHistory = zeros(length(wfs.slopes),batchItSize);
    end
    if SAVEDM
        if exist(outputDir+"\"+fileID_DM+string(i)+".h5", 'file'), delete(outputDir+"\"+fileID_DM+string(i)+".h5"); end
        dmCommandsHistory = zeros(length(dm.coefs),batchItSize);
    end
    if SAVELIGHTFIELD
        if exist(outputDir+"\"+fileID_lightfield+string(i)+".h5", 'file'), delete(outputDir+"\"+fileID_lightfield+string(i)+".h5"); end
        lightfieldHistory = zeros(size(wfs.camera.frame, 1), size(wfs.camera.frame, 2),batchItSize, sensor_type);
    end
    if SAVEPSF
        if exist(outputDir+"\"+fileID_ipsf+string(i)+".h5", 'file'), delete(outputDir+"\"+fileID_ipsf+string(i)+".h5"); end
        psfHistory = zeros(size(instantCam.frame, 1), size(instantCam.frame, 2),batchItSize, sensor_type);
    end
    if SAVERWFE
        if exist(outputDir+"\"+fileID_rwfe+string(i)+".h5", 'file'), delete(outputDir+"\"+fileID_rwfe+string(i)+".h5"); end
        % rwfe_waves_history = zeros(batchItSize,1);
    end
    if SAVEDIFFLIMITED
        if exist(outputDir+"\"+fileID_diff_limited+".h5", 'file'), delete(outputDir+"\"+fileID_diff_limited+".h5"); end
        diff_limited = zeros(size(cam.frame));
    end
end
if exist(outputDir+"\"+fileID_metadata+".txt", 'file'), delete(outputDir+"\"+fileID_metadata+".txt"); end

rwfe_waves_history = zeros(batchItSize,1);

%% Regulation

lag_buffer = zeros(length(dm.coefs),lag_c+1);   % +1 becaus Matlab arrays start at 1. Also it allows for a lag of 0 to work without special case handling.
flush(cam)
flush(instantCam)
indexInBatch = 0;
for k=1:nIteration
    indexInBatch = indexInBatch + 1;
    % Objects update
    instantCam.frame = [];
    instantCam.frameCount = 0;
    +tel;
    +ngs;
    +science;
    +instantScience;

    % Closed-loop controller
    lag_buffer(:,lag_c+1) = min(max(dm.coefs - gain_cl*calibDm.M*wfs.slopes, -1), 1);
    % Moving the last element of the buffer to the first position and shifting the rest to the right
    for j = 1: lag_c
        lag_buffer(:,j) = lag_buffer(:,j+1);
    end

    dm.coefs = lag_buffer(:,1);

        % local log
    if SAVEWF
        WFHistory(:,:,indexInBatch) = ngs.meanRmPhase;
    end
    if SAVESLOPES
        WFSHistory(:,indexInBatch) = wfs.slopes;
    end
    if SAVELIGHTFIELD
        lightfieldHistory(:,:,indexInBatch) = type_cast(wfs.camera.frame, sensor_type);
    end
    if SAVEDM
        dmCommandsHistory(:, indexInBatch) = dm.coefs;
    end
    if SAVERWFE
        rwfe_waves_history(indexInBatch) = sqrt(var(ngs))./2/pi; % [waves]
    end
    if SAVEPSF
        psfHistory(:,:,indexInBatch) = type_cast(instantCam.frame, sensor_type);
    end
    if mod(k-1, round(nIteration/50)) == 0 || k == nIteration
        fprintf('Progress: %d%% done\n', round(100*k/nIteration));
    end
    % fprintf('Current iteration: %d\n', mod(k, batchItSize))
    if mod(k, batchItSize) == 0 || k == nIteration
        fprintf('Saving batch %d/%d to disk...\n', ceil(k/batchItSize), nBatch);
        batchIndex = ceil(k/batchItSize);

        if SAVEWF
            sz = size(WFHistory);
            h5create(outputDir+"\"+fileID_WF+string(batchIndex)+".h5", '/wf', sz, 'ChunkSize', [sz(1) sz(2) 1], 'DataType', 'double');
            h5write(outputDir+"\"+fileID_WF+string(batchIndex)+".h5", '/wf', WFHistory);
        end
        if SAVESLOPES
            sz = size(WFSHistory);
            h5create(outputDir+"\"+fileID_WFS+string(batchIndex)+".h5", '/wfs', sz, 'ChunkSize', [sz(1) 1], 'DataType', 'double');
            h5write(outputDir+"\"+fileID_WFS+string(batchIndex)+".h5", '/wfs', WFSHistory);
            WFSHistory = zeros(length(wfs.slopes),batchItSize);
        end
        if SAVELIGHTFIELD
            sz = size(lightfieldHistory);
            totBytes = prod(sz);   % for uint8, this is roughly the expected file size
            % fprintf('Dataset size: [%d,%d,%d] -> %d bytes\n', sz, totBytes);
            h5create(outputDir+"\"+fileID_lightfield+string(batchIndex)+".h5", '/wf_lightfield', sz, 'ChunkSize', [sz(1) sz(2) 1],'DataType', sensor_type);
            h5write(outputDir+"\"+fileID_lightfield+string(batchIndex)+".h5", '/wf_lightfield', lightfieldHistory);
            lightfieldHistory = zeros(size(wfs.camera.frame, 1), size(wfs.camera.frame, 2),batchItSize, sensor_type);
        end
        if SAVEDM
            sz = size(dmCommandsHistory);
            h5create(outputDir+"\"+fileID_DM+string(batchIndex)+".h5", '/dm_commands', sz, 'ChunkSize', [sz(1) 1], 'DataType', 'double');
            h5write(outputDir+"\"+fileID_DM+string(batchIndex)+".h5", '/dm_commands', dmCommandsHistory);
            dmCommandsHistory = zeros(length(dm.coefs),batchItSize);
        end
        if SAVEPSF
            sz = size(psfHistory);
            h5create(outputDir+"\"+fileID_ipsf+string(batchIndex)+".h5", '/psf_history', sz, 'ChunkSize', [sz(1) sz(2) 1], 'DataType', sensor_type);
            h5write(outputDir+"\"+fileID_ipsf+string(batchIndex)+".h5", '/psf_history', psfHistory);
            psfHistory = zeros(size(instantCam.frame, 1), size(instantCam.frame, 2),batchItSize);
        end
        if SAVERWFE
            sz = size(rwfe_waves_history);
            h5create(outputDir+"\"+fileID_rwfe+string(batchIndex)+".h5", '/rwfe_waves_history', sz, 'ChunkSize', [sz(1) 1], 'DataType', 'double');
            h5write(outputDir+"\"+fileID_rwfe+string(batchIndex)+".h5", '/rwfe_waves_history', rwfe_waves_history);
            rwfe_waves_history = zeros(batchItSize,1);
        end
        indexInBatch = 0;
    end
end
%%
if SAVEDIFFLIMITED
    sz = size(diff_limited);
    h5create(outputDir+"\"+fileID_diff_limited+".h5", '/diff_limited', sz, 'ChunkSize', [sz(1) 1]);
    h5write(outputDir+"\"+fileID_diff_limited+".h5", '/diff_limited', diff_limited);
end
%%
%%
% rowNames = {'D';'r0';'L0';'Asl';'wind';'windDir';'Exposure time';'nIteration';'gain_cl';'batchItSize';'nBatch';'LastBatchItSize'; 'oversampling'; 'nActWSF';'edge_act'; 'startDelay'; 'longStrehl'; 'magnitude'};
% values =    [D;  r0;  L0;  Asl;  wind;  windDir;  exposureTime;   nIteration;  gain_cl;  batchItSize;  nBatch;  LastBatchItSize; oversampling; nActWSF; edge_act; startDelay; cam.strehl; NGSmagnitude];
% T = table(values,'RowNames',rowNames);
% writetable(T,outputDir+"\"+fileID_metadata+".txt",'Delimiter','\t','WriteRowNames',true);

if exist(metadataFile, 'file')
    delete(metadataFile)
end
fidIn = fopen(input_file, 'r');
fidMeta = fopen(metadataFile, 'w');
while ~feof(fidIn)
    line = fgetl(fidIn);
    if ischar(line)
        fprintf(fidMeta, '%s\n', line);
    end
end

fprintf(fidMeta, '\n\n');
fprintf(fidMeta, '---------------------- OUTPUTS----------------------\n\n');
fprintf(fidMeta, 'batchItSize = %d\n', batchItSize);
fprintf(fidMeta, 'nBatch = %d\n', nBatch);
fprintf(fidMeta, 'LastBatchItSize = %d\n', LastBatchItSize);
fprintf(fidMeta, 'LongStrehl = %d\n', cam.strehl);

fclose(fidIn);
fclose(fidMeta);
%% s
% maxValue = max(psfHistory, [], 'all');
% fprintf('Maximum value in frame 30: %f\n', maxValue);
% % figure;
% fprintf('max long psf value: %f\n', sum(cam.frame(:)))

% figure;
% imshow(psfHistory(:,:,21), []);
%%
% psf_sum = sum(psfHistory(:,:,startDelay+1:end), 3);   

% fprintf('Strehl ratio: %4.1f\n',cam.strehl);
% fprintf('Strehl ratio: %4.1f\n',instantCam.strehl);

% figure;
% subplot(2,1,1);
% plot(rwfe_waves_history, 'o-');
% xlabel('Iteration'); ylabel('Residual WF RMS (waves)');
% title('AO Loop Convergence (Linear Scale)');
% grid on;
% subplot(2,1,2);
% semilogy(rwfe_waves_history, 'o-');
% xlabel('Iteration'); ylabel('Residual WF RMS (waves)');
% title('AO Loop Convergence (Logarithmic Scale)');
% grid on;


% figure;
% imagesc(cam,'parent',subplot(2,2,1));
% title('Long Exposure PSF', 'FontSize', 12, 'FontWeight', 'bold');
% colorbar; axis image;
% imagesc(instantCam,'parent',subplot(2,2,2));
% title('Instantaneous PSF', 'FontSize', 12, 'FontWeight', 'bold');
% colorbar; axis image;
% imagesc(psf_sum,'parent',subplot(2,2,3));
% title('Long psf from instantaneous PSF', 'FontSize', 12, 'FontWeight', 'bold');
% colorbar; axis image;
% imagesc(cam.frame-psf_sum,'parent',subplot(2,2,4));
% title('Long psf - iPsfSum', 'FontSize', 12, 'FontWeight', 'bold');
% colorbar; axis image;
% sgtitle(sprintf('AO Strehl: Long=%.2f, Instant=%.2f', cam.strehl, instantCam.strehl));


% psf_sum_flux = sum(psf_sum(:));
% long_psf_flux = sum(cam.frame(:));

% fprintf('Flux in long exposure PSF: %.2e\n', long_psf_flux);
% fprintf('Flux in sum of instantaneous PSFs: %.2e\n', psf_sum_flux);

% flux_ratio = psf_sum_flux / long_psf_flux;
% fprintf('Flux ratio (iPsfSum / Long PSF): %.3f\n', flux_ratio);


%% clear

if exist('WFHistory', 'var'), clear WFHistory; end
if exist('WFSHistory', 'var'), clear WFSHistory; end
if exist('lightfieldHistory', 'var'), clear lightfieldHistory; end
if exist('dmCommandsHistory', 'var'), clear dmCommandsHistory; end
if exist('psfHistory', 'var'), clear psfHistory; end
if exist('rwfe_waves_history', 'var'), clear rwfe_waves_history; end