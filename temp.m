clc
clear all
addpath('OOMAO')


tel = telescope(1,'resolution', 100, 'fieldOfViewInArcsec',30,'samplingTime',1/500);


wfs = shackHartmann(10,100, 0.5);


ngs = source('wavelength',photometry.HeNe);
science = source('wavelength',photometry.HeNe); % same source for experiment


bifa = influenceFunction('monotonic', 0.75);
dm = deformableMirror(11, 'modes', bifa, 'resolution', tel.resolution, 'validActuator', wfs.validActuator);


atm = atmosphere(photometry.HeNe, 15e-2, 30, 'altitude', 5e3, 'windSpeed', 10, 'windDirection', pi/3);


cam = imager();


ngs = ngs.*tel*wfs
science = science.*tel*cam

wfs.INIT

calibDm = calibration(dm, wfs, ngs, ngs.wavelength);
calibDm.threshold = 1e6;
disp(calibDm)



%% closed-loop adaptive optics
gain = 0.5;
dm.coefs = 0;
ngs.logging = true;
ngs = ngs.*tel*dm*wfs;
figure
h = imagesc(ngs.meanRmOpd*1e6);
colorbar
cam.exposureTime = 150;
cam.clockRate = 1;
science = science.*tel*dm*cam;
dmCoefs = size(dm.coefs,2);
pause(1)
for k=1:150
    +tel
    +ngs
    +science
    dm.coefs = dm.coefs - gain*calibDm.M*wfs.slopes;
%     dmCoefs(:,2) = dmCoefs(:,1) - gain*calibDm.M*wfs.slopes;
%     dm.coefs = dmCoefs(:,1);
%     dmCoefs(:,1) = dmCoefs(:,2);
    set(h,'Cdata',ngs.meanRmOpd*1e6)
    drawnow
end

%% reporting performance
hf = figure;
plot(1e6*sqrt(ngs.phaseVar(1:150*2))/ngs.waveNumber,'.')