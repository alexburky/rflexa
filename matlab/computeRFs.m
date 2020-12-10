% computeRFs.m
%
% This script will loop over a directory of data, and calculate receiver
% functions. As part of the algorithm, it will remove the instrument
% response from each trace using a corresponding poles/zeros file

% Overview: Should be a function which takes in a directory as input,
% as well as some parameters controlling the resulting receiver functions

%--------------------------------------------------------------------------
% Last updated 12/10/2020 by aburky@princeton.edu
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
