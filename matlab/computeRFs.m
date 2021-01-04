% computeRFs.m
%
% This script will loop over a directory of data, and calculate receiver
% functions. As part of the algorithm, it will remove the instrument
% response from each trace using a corresponding poles/zeros file

% Overview: Should be a function which takes in a directory as input,
% as well as some parameters controlling the resulting receiver functions

%--------------------------------------------------------------------------
% Last updated 1/4/2020 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% -------------------------------------------------------------------------
% Read in seismic data
% -------------------------------------------------------------------------
dataDir = '/Users/aburky/PycharmProjects/bermudaRFs/matlab/';

% First, try for files with 'BHE' or 'BHN' naming convention
eData = dir(fullfile(dataDir,'*BHE*SAC'));
nData = dir(fullfile(dataDir,'*BHN*SAC'));
% If no data read in, assume files use 'BH1' or 'BH2' convention
if isempty(eData)
    eData = dir(fullfile(dataDir,'*BH1*SAC'));
elseif isempty(nData)
    nData = dir(fullfile(dataDir,'*BH2*SAC'));
end
zData = dir(fullfile(dataDir,'*BHZ*SAC'));

nFiles = [length(eData), length(nData), length(zData)];

% Poles and Zeros data
pzData = dir(fullfile(dataDir,'SAC_PZs*'));

% Read the data and save it into a structure, 'sac'
for i = 1:max(nFiles)
    try
        [sac{i}.te,sac{i}.de,sac{i}.he] = fread_sac(fullfile(dataDir,...
                                            eData(i).name));
    catch
        disp('Reached the limit for E channel')
    end
    try
        [sac{i}.tn,sac{i}.dn,sac{i}.hn] = fread_sac(fullfile(dataDir,...
                                            nData(i).name));
    catch
        disp('Reached the limit for N channel')
    end
    try
        [sac{i}.tz,sac{i}.dz,sac{i}.hz] = fread_sac(fullfile(dataDir,...
                                            zData(i).name));
    catch
        disp('Reached the limit for Z channel')
    end
end

% -------------------------------------------------------------------------
% Remove the instrument response for each channel
% -------------------------------------------------------------------------

% Instrument response removal parameters
r = 0.1;
freqlims = [0.002 0.004 5 10];

% East Channel
for i = 1:length(eData)
    sacIdx = extractBefore(eData(i).name,'.SAC');
    sacIdx = sacIdx(end);
    % Pre-process the data
    sac{i}.de = sac{i}.de - mean(sac{i}.de);
    sac{i}.de = detrend(sac{i}.de);
    sac{i}.de = sac{i}.de.*tukeywin(sac{i}.he.npts,r);
    for j = 1:length(pzData)
        pzFile = pzData(j).name;
        pzIdx = pzFile(end);
        pzFile = fullfile(dataDir,pzFile);
        % If the SAC data file index matches, remove the response!
        if strcmp(sacIdx,pzIdx)
            sac{i}.de = transfer(sac{i}.de,sac{i}.he.delta,freqlims,...
                                 'velocity',pzFile);
        end
    end
end

% North Channel
for i = 1:length(nData)
    sacIdx = extractBefore(nData(i).name,'.SAC');
    sacIdx = sacIdx(end);
    % Pre-process the data
    sac{i}.dn = sac{i}.dn - mean(sac{i}.dn);
    sac{i}.dn = detrend(sac{i}.dn);
    sac{i}.dn = sac{i}.dn.*tukeywin(sac{i}.hn.npts,r);
    for j = 1:length(pzData)
        pzFile = pzData(j).name;
        pzIdx = pzFile(end);
        pzFile = fullfile(dataDir,pzFile);
        % If the SAC data file index matches, remove the response!
        if strcmp(sacIdx,pzIdx)
            sac{i}.dn = transfer(sac{i}.dn,sac{i}.hn.delta,freqlims,...
                                 'velocity',pzFile);
        end
    end
end

% Vertical Channel
for i = 1:length(zData)
    sacIdx = extractBefore(zData(i).name,'.SAC');
    sacIdx = sacIdx(end);
    % Pre-process the data
    sac{i}.dz = sac{i}.dz - mean(sac{i}.dz);
    sac{i}.dz = detrend(sac{i}.dz);
    sac{i}.dz = sac{i}.dz.*tukeywin(sac{i}.hz.npts,r);
    for j = 1:length(pzData)
        pzFile = pzData(j).name;
        pzIdx = pzFile(end);
        pzFile = fullfile(dataDir,pzFile);
        % If the SAC data file index matches, remove the response!
        if strcmp(sacIdx,pzIdx)
            sac{i}.dz = transfer(sac{i}.dz,sac{i}.hz.delta,freqlims,...
                                 'velocity',pzFile);
        end
    end
end

%% ------------------------------------------------------------------------
% Rotate Horizontal Components
% -------------------------------------------------------------------------

% Iterate over the longer of the two horizontal file lists
% nFiles = [length(eData), length(nData)];

% for i = 1:max(nFiles)
disp('Rotation loop execution time:')
tic
for i = 1:length(eData)
    % Check if the two files correspond to the same event
    for j = 1:length(nData)
        if sac{i}.he.nzjday == sac{j}.hn.nzjday
            if sac{i}.he.nzhour == sac{j}.hn.nzhour
                [sac{j}.dn,sac{i}.de] = seisne(sac{j}.dn,sac{i}.de,...
                                                sac{j}.hn.cmpaz);
                [sac{j}.dr,sac{i}.dt] = seisrt(sac{j}.dn,sac{i}.de,...
                                                sac{j}.hn.baz);
                break
            end
        end
    end
end
toc

%% ------------------------------------------------------------------------
% Pre-Process Data: Cut, Taper, and Optionally Filter
% -------------------------------------------------------------------------

% First, optionally filter data
fc = [0.02 0.2];
for i = 1:length(nData)
    fs = 1/sac{i}.hn.delta;
    [b,a] = butter(3,fc/(fs/2),'bandpass');
    % sac{i}.dr = filtfilt(b,a,sac{i}.dr);
    sac{i}.dr = filter(b,a,sac{i}.dr);
end

for i = 1:length(zData)
    fs = 1/sac{i}.hz.delta;
    [b,a] = butter(3,fc/(fs/2),'bandpass');
    % sac{i}.dz = filtfilt(b,a,sac{i}.dz);
    sac{i}.dz = filter(b,a,sac{i}.dz);
end

% Second, cut and taper the data
cut_b = 30;
cut_e = 90;
taperw = 0.25;

for i = 1:length(nData)
    % Get P-wave arrival time from 'T0' header
    pidx = fix(sac{i}.hn.t(1)/sac{i}.hn.delta);
    bidx = pidx - fix(cut_b/sac{i}.hn.delta);
    eidx = pidx + fix(cut_e/sac{i}.hn.delta);
    
    sac{i}.drc = sac{i}.dr(bidx:eidx);
    sac{i}.drc = sac{i}.drc.*tukeywin(length(sac{i}.drc),taperw);
end

for i = 1:length(zData)
    pidx = fix(sac{i}.hz.t(1)/sac{i}.hz.delta);
    bidx = pidx - fix(cut_b/sac{i}.hz.delta);
    eidx = pidx + fix(cut_e/sac{i}.hz.delta);
    
    sac{i}.dzc = sac{i}.dz(bidx:eidx);
    sac{i}.dzc = sac{i}.dzc.*tukeywin(length(sac{i}.dzc),taperw);
end

%% ------------------------------------------------------------------------
% Calculate receiver functions!
% -------------------------------------------------------------------------

% In this section, need to check that the vertical and radial component
% correspond to the same event before making receiver function!

gw = 1.0;
tshift = 10;
itmax = 1000;
tol = 0.001;

for i = 1:length(nData)
    for j = 1:length(zData)
        if sac{i}.hn.nzjday == sac{j}.hz.nzjday
            if sac{i}.hn.nzhour == sac{j}.hz.nzhour
                npts = length(sac{i}.drc);
                [rf{i}.d,rf{i}.rms] = makeRFitdecon_la(sac{i}.drc,...
                    sac{j}.dzc,sac{i}.hn.delta,npts,tshift,gw,itmax,tol);
                rf{i}.t = 0:sac{i}.hn.delta:(length(rf{i}.d)-1)*...
                    sac{i}.hn.delta;
                rf{i}.h = sac{j}.hz;
            end
        end
    end
end

%% ------------------------------------------------------------------------
% Save the resulting receiver function data to SAC files
% -------------------------------------------------------------------------

% Write a separate function - saveRF.m - to do this...
% Model it after 'saveSAC.m'

for i = 1:length(rf)
    % saveRF.m
end

%% Compare to results using Python...

pyDir = '/Users/aburky/IFILES/NETWORKS_TEST/TA/N61A/NULL/RFUNCS_VEL/FILTERED_0.02_0.2/GW10/';
pyFile = '2014.04.18.14.27.24.TA.N61A.NULL.RF.SAC';

[py.t,py.d,py.h] = fread_sac(fullfile(pyDir,pyFile));

plot(py.t,py.d,'k')
hold on
plot(rf{3}.t,rf{3}.d,'r')
