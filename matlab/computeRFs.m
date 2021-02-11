% computeRFs.m
%
% This script will loop over a directory of data, and calculate receiver
% functions. As part of the algorithm, it will remove the instrument
% response from each trace using a corresponding poles/zeros file

% Overview: Should be a function which takes in a directory as input,
% as well as some parameters controlling the resulting receiver functions

%--------------------------------------------------------------------------
% Last updated 2/11/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

tic

% -------------------------------------------------------------------------
% Read in seismic data
% -------------------------------------------------------------------------

% Location of seismic data
dataDir = '/Users/aburky/IFILES/NETWORKS/II/SACV/00/RFQUAKES/';

% Network (for filename convention)
network = 'II';

% Directory where receiever functions will be saved
rfDir = '/Users/aburky/IFILES/NETWORKS/II/SACV/00/RFUNCS_VEL/';
if ~exist(rfDir,'dir')
    mkdir(rfDir)
elseif exist(rfDir,'dir') == 7
    disp('Receiver function directory already exists! Deleting it...')
    rmdir(rfDir,'s')
    mkdir(rfDir)
end

% Files with 'BHE' or 'BHN' naming convention
eData = dir(fullfile(dataDir,'*BHE*SAC'));
nData = dir(fullfile(dataDir,'*BHN*SAC'));

% Files with 'BH1' or 'BH2' naming convention
nData1 = dir(fullfile(dataDir,'*BH1*SAC'));
eData2 = dir(fullfile(dataDir,'*BH2*SAC'));

% Vertical component files
zData = dir(fullfile(dataDir,'*BHZ*SAC'));

% Determine which component has the most data
nFiles = [length(eData)+length(eData2), length(nData)+length(nData1), ...
          length(zData)];    
[~,idx] = max(nFiles);

% Poles and Zeros data
pzData = dir(fullfile(dataDir,'SAC_PZs*'));

% Read in events for which all 3 components exist
j = 1;
for i = 1:max(nFiles)
    % 'BHE' channel has most data
    if idx == 1
        if i <= length(eData)
            evh = strsplit(eData(i).name,network);
            traces = dir(fullfile(dataDir,[evh{1},'*SAC']));
            if length(traces) == 3
                [sac{j}.te,sac{j}.de,sac{j}.he] = fread_sac(fullfile(...
                    dataDir,traces(1).name));
                sac{j}.efile = traces(1).name;
                [sac{j}.tn,sac{j}.dn,sac{j}.hn] = fread_sac(fullfile(...
                    dataDir,traces(2).name));
                sac{j}.nfile = traces(2).name;
                [sac{j}.tz,sac{j}.dz,sac{j}.hz] = fread_sac(fullfile(...
                    dataDir,traces(3).name));
                sac{j}.zfile = traces(3).name;
                j = j + 1;
            end
        else
            evh = strsplit(eData2(i - length(eData)).name,network);
            traces = dir(fullfile(dataDir,[evh{1},'*SAC']));
            if length(traces) == 3
                [sac{j}.te,sac{j}.de,sac{j}.he] = fread_sac(fullfile(...
                    dataDir,traces(2).name));
                sac{j}.efile = traces(2).name;
                [sac{j}.tn,sac{j}.dn,sac{j}.hn] = fread_sac(fullfile(...
                    dataDir,traces(1).name));
                sac{j}.nfile = traces(1).name;
                [sac{j}.tz,sac{j}.dz,sac{j}.hz] = fread_sac(fullfile(...
                    dataDir,traces(3).name));
                sac{j}.zfile = traces(3).name;
                j = j + 1;
            end
        end
    % 'BHN' channel has most data
    elseif idx == 2
        if i <= length(nData)
            evh = strsplit(nData(i).name,network);
            traces = dir(fullfile(dataDir,[evh{1},'*SAC']));
            if length(traces) == 3
                [sac{j}.te,sac{j}.de,sac{j}.he] = fread_sac(fullfile(...
                    dataDir,traces(1).name));
                sac{j}.efile = traces(1).name;
                [sac{j}.tn,sac{j}.dn,sac{j}.hn] = fread_sac(fullfile(...
                    dataDir,traces(2).name));
                sac{j}.nfile = traces(2).name;
                [sac{j}.tz,sac{j}.dz,sac{j}.hz] = fread_sac(fullfile(...
                    dataDir,traces(3).name));
                sac{j}.zfile = traces(3).name;
                j = j + 1;
            end
        else
            evh = strsplit(nData1(i - length(nData)).name,network);
            traces = dir(fullfile(dataDir,[evh{1},'*SAC']));
            if length(traces) == 3
                [sac{j}.te,sac{j}.de,sac{j}.he] = fread_sac(fullfile(...
                    dataDir,traces(2).name));
                sac{j}.efile = traces(2).name;
                [sac{j}.tn,sac{j}.dn,sac{j}.hn] = fread_sac(fullfile(...
                    dataDir,traces(1).name));
                sac{j}.nfile = traces(1).name;
                [sac{j}.tz,sac{j}.dz,sac{j}.hz] = fread_sac(fullfile(...
                    dataDir,traces(3).name));
                sac{j}.zfile = traces(3).name;
                j = j + 1;
            end
            
        end
    % 'BHZ' channel has most data
    elseif idx == 3
        evh = strsplit(zData(i).name,network);
        traces = dir(fullfile(dataDir,[evh{1},'*SAC']));
        if length(traces) == 3
            if contains(traces(1).name,'BHE')
                [sac{j}.te,sac{j}.de,sac{j}.he] = fread_sac(fullfile(...
                    dataDir,traces(1).name));
                sac{j}.efile = traces(1).name;
            elseif contains(traces(1).name,'BH1')
                [sac{j}.tn,sac{j}.dn,sac{j}.hn] = fread_sac(fullfile(...
                    dataDir,traces(1).name));
                sac{j}.nfile = traces(1).name;
            elseif contains(traces(2).name,'BHN')
                [sac{j}.tn,sac{j}.dn,sac{j}.hn] = fread_sac(fullfile(...
                    dataDir,traces(2).name));
                sac{j}.nfile = traces(2).name;
            elseif contains(traces(2).name,'BH2')
                [sac{j}.te,sac{j}.de,sac{j}.he] = fread_sac(fullfile(...
                    dataDir,traces(2).name));
                sac{j}.efile = traces(2).name;
            end
            [sac{j}.tz,sac{j}.dz,sac{j}.hz] = fread_sac(fullfile(...
                dataDir,traces(3).name));
            sac{j}.zfile = traces(3).name;
            j = j + 1;
        end
    end
end

% -------------------------------------------------------------------------
% Remove the instrument response for each channel
% -------------------------------------------------------------------------

% Instrument response removal parameters
r = 0.1;
freqlims = [0.002 0.004 5 10];

for i = 1:length(sac)
    % Pre-process East Channel
    eIdx = extractBefore(sac{i}.efile,'.SAC');
    eIdx = eIdx(end);
    sac{i}.de = sac{i}.de - mean(sac{i}.de);
    sac{i}.de = detrend(sac{i}.de);
    sac{i}.de = sac{i}.de.*tukeywin(sac{i}.he.npts,r);
    
    % Pre-process North Channel
    nIdx = extractBefore(sac{i}.nfile,'.SAC');
    nIdx = nIdx(end);
    sac{i}.dn = sac{i}.dn - mean(sac{i}.dn);
    sac{i}.dn = detrend(sac{i}.dn);
    sac{i}.dn = sac{i}.dn.*tukeywin(sac{i}.hn.npts,r);
    
    % Pre-process Vertical Channel
    zIdx = extractBefore(sac{i}.zfile,'.SAC');
    zIdx = zIdx(end);
    sac{i}.dz = sac{i}.dz - mean(sac{i}.dz);
    sac{i}.dz = detrend(sac{i}.dz);
    sac{i}.dz = sac{i}.dz.*tukeywin(sac{i}.hz.npts,r);
    
    % Remove instrument response from each channnel
    for j = 1:length(pzData)
        pzFile = pzData(j).name;
        pzIdx = pzFile(end);
        pzFile = fullfile(dataDir,pzFile);
        if strcmp(eIdx,pzIdx)
            sac{i}.de = transfer(sac{i}.de,sac{i}.he.delta,freqlims,...
                'velocity',pzFile);
        elseif strcmp(nIdx,pzIdx)
            sac{i}.dn = transfer(sac{i}.dn,sac{i}.hn.delta,freqlims,...
                'velocity',pzFile);
        elseif strcmp(zIdx,pzIdx)
            sac{i}.dz = transfer(sac{i}.dz,sac{i}.hz.delta,freqlims,...
                'velocity',pzFile);
        end
    end
end

% -------------------------------------------------------------------------
% Rotate Horizontal Components
% -------------------------------------------------------------------------

% Iterate over the longer of the two horizontal file lists
% nFiles = [length(eData), length(nData)];

k = 1;
for i = 1:length(sac)
    % Make sure lengths of data vectors are equal before rotating
    if isequal(length(sac{i}.dn),length(sac{i}.de))
        [sac{i}.dn,sac{i}.de] = seisne(sac{i}.dn,sac{i}.de,...
            sac{i}.hn.cmpaz);
        [sac{i}.dr,sac{i}.dt] = seisrt(sac{i}.dn,sac{i}.de,...
            sac{i}.hn.baz);
        rGood(k) = i;
        tGood(k) = i;
        k = k + 1;
    end
end

% -------------------------------------------------------------------------
% Pre-Process Data: Cut, Taper, and Optionally Filter
% -------------------------------------------------------------------------

% First, optionally filter data
fc = [0.02 0.2];

% Radial component
for i = 1:length(rGood)
    fs = 1/sac{rGood(i)}.hn.delta;
    [b,a] = butter(3,fc/(fs/2),'bandpass');
    % sac{i}.dr = filtfilt(b,a,sac{i}.dr);
    sac{rGood(i)}.dr = filter(b,a,sac{rGood(i)}.dr);
end

% Vertical component
for i = 1:length(sac)
    fs = 1/sac{i}.hz.delta;
    [b,a] = butter(3,fc/(fs/2),'bandpass');
    % sac{i}.dz = filtfilt(b,a,sac{i}.dz);
    sac{i}.dz = filter(b,a,sac{i}.dz);
end

% Second, cut and taper the data
cut_b = 30;
cut_e = 90;
taperw = 0.25;

% Radial component
for i = 1:length(rGood)
    % Get P-wave arrival time from 'T0' header
    pidx = fix(sac{rGood(i)}.hn.t(1)/sac{rGood(i)}.hn.delta);
    bidx = pidx - fix(cut_b/sac{rGood(i)}.hn.delta);
    eidx = pidx + fix(cut_e/sac{rGood(i)}.hn.delta);
    
    sac{rGood(i)}.drc = sac{rGood(i)}.dr(bidx:eidx);
    sac{rGood(i)}.drc = sac{rGood(i)}.drc.*...
                        tukeywin(length(sac{rGood(i)}.drc),taperw);
end

% Vertical component
for i = 1:length(sac)
    pidx = fix(sac{i}.hz.t(1)/sac{i}.hz.delta);
    bidx = pidx - fix(cut_b/sac{i}.hz.delta);
    eidx = pidx + fix(cut_e/sac{i}.hz.delta);
    
    sac{i}.dzc = sac{i}.dz(bidx:eidx);
    sac{i}.dzc = sac{i}.dzc.*tukeywin(length(sac{i}.dzc),taperw);
end

% -------------------------------------------------------------------------
% Calculate receiver functions!
% -------------------------------------------------------------------------

% In this section, need to check that the vertical and radial component
% correspond to the same event before making receiver function!

gw = 1.0;
tshift = 10;
itmax = 1000;
tol = 0.001;

j = 1;
for i = 1:length(rGood)
    if isequal(length(sac{rGood(i)}.dz),length(sac{rGood(i)}.dr))
        npts = length(sac{rGood(i)}.drc);
        [rf{j}.d,rf{j}.rms] = makeRFitdecon_la(sac{rGood(i)}.drc,...
            sac{rGood(i)}.dzc,sac{rGood(i)}.hz.delta,npts,tshift,gw,...
            itmax,tol);
        rf{j}.t = 0:sac{rGood(i)}.hz.delta:(length(rf{j}.d)-1)*...
            sac{rGood(i)}.hz.delta;
        rf{j}.h = sac{rGood(i)}.hz;
        % Get indices for calculating signal to noise ratio
        % (using method of Gao and Liu, 2014)
        noise_b = round((cut_b - 20)*(1/sac{rGood(i)}.hz.delta))-1;
        noise_e = round((cut_b - 10)*(1/sac{rGood(i)}.hz.delta));
        signal_b = round((cut_b - 8)*(1/sac{rGood(i)}.hz.delta))-1;
        signal_e = round((cut_b + 12)*(1/sac{rGood(i)}.hz.delta));
        % Calculate vertical component SNR
        vn = abs(mean(sac{rGood(i)}.dzc(noise_b:noise_e)));
        vs = max(abs(sac{rGood(i)}.dzc(signal_b:signal_e)));
        rf{j}.vsnr = vs/vn;
        % Calculate radial component SNR
        rn = abs(mean(sac{rGood(i)}.drc(noise_b:noise_e)));
        rs = max(abs(sac{rGood(i)}.drc(signal_b:signal_e)));
        rf{j}.rsnr = rs/rn;
        j = j + 1;
    end
end   

% -------------------------------------------------------------------------
% Save the resulting receiver function data to SAC files
% -------------------------------------------------------------------------

for i = 1:length(rf)
    saveRF(rf{i},rfDir);
end

disp('Execution time:')
toc

%% Compare to results using Python...

% pyDir = '/Users/aburky/IFILES/NETWORKS_TEST/TA/N61A/NULL/RFUNCS_VEL/FILTERED_0.02_0.2/GW10/';
% pyFile = '2014.04.18.14.27.24.TA.N61A.NULL.RF.SAC';

% [py.t,py.d,py.h] = fread_sac(fullfile(pyDir,pyFile));

% plot(py.t,py.d,'k')
% hold on
% plot(rf{3}.t,rf{3}.d,'r')
