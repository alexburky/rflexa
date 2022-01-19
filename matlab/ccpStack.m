% ccpStack.m
%
% This script is used to make CCP stacks of the TA dataset (US Array) with
% a focus on the mantle transition zone. It is assumed that the data are
% structured in a specific directory structure, and that the paths to
% relevant piercing point data are specified. The output of this script is
% a .mat file with the CCP stack.
%
%--------------------------------------------------------------------------
% Last updated 1/19/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% Path to piercing point data
pcDir{1} = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/piercePoints/',...
            'pierceMoho/'];
pcDir{2} = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/piercePoints/',...
            'pierce210/'];
pcDir{3} = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/piercePoints/',...
            'pierce410/'];
pcDir{4} = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/piercePoints/',...
            'pierce660/'];

% Path to receiver function data        
rfDir = '/Users/aburky/IFILES/NETWORKS/TA/';

% Velocity model to use
vModel = 'LLNL';

if strcmp(vModel,'iasp91')
    rfDet = '/NULL/RFUNCS_VEL/FILTERED_0.02_0.2/GW10/';
elseif strcmp(vModel,'LLNL')
    rfDet = '/NULL/RFUNCS_VEL/FILTERED_0.02_0.2/GW10/LLNL3D/';
end
        
%--------------------------------------------------------------------------
% Load all piercing point data
%--------------------------------------------------------------------------        
     
% Piercing point files
for i = 1:length(pcDir)
    pcFile{i} = dir(fullfile(pcDir{i},'*.dat'));
end

% Read in all piercing point data
for j = 1:length(pcDir)
    
    pcNam = [];
    pcLat = [];
    pcLon = [];
    for i = 1:length(pcFile{j})
        % pcData{i} = dlmread(fullfile(pcDir,pcFile(i).name),'',0,1);
        % pcData{i} = importdata(fullfile(pcDir,pcFile(i).name));
        pcData = readtable(fullfile(pcDir{j},pcFile{j}(i).name));
        
        % Only save variables if data isn't empty
        if ~isempty(pcData)
            pcNam{i} = pcData.Var1;
            pcLat{i} = pcData.Var2;
            pcLon{i} = pcData.Var3;
        end
    end
    
    file{j} = [];
    lat{j} = [];
    lon{j} = [];
    % Concatenate data?
    for i = 1:length(pcFile{j})
        file{j} = [file{j}; pcNam{i}];
        lat{j} = [lat{j}; pcLat{i}];
        lon{j} = [lon{j}; pcLon{i}];
    end
    
end

%--------------------------------------------------------------------------
% Bin the data
%--------------------------------------------------------------------------

% Define lat/lon limits
lonBin = -98.5:0.5:-63.5;
latBin = 22.5:0.5:51.5;

for i = 1:length(pcDir)
    [N{i},~,~,binX{i},binY{i}] = histcounts2(lon{i},lat{i},lonBin,latBin);
end

%--------------------------------------------------------------------------
% Calculate the CCP stack!
%--------------------------------------------------------------------------

% Define depth increment
dz = 0.1;
z = 0:dz:800;

% Stack 410 bits
nRows = size(N{3},1);
nCols = size(N{3},2);
ccp410 = cell(nRows,nCols);
for i = 1:nRows
    for j = 1:nCols
        sbin = [i j];
        [tf,loc] = ismember([binX{3},binY{3}],sbin,'row');
        idx = find(loc);
        
        % Read in all of the receiver functions in this bin
        rf = {};
        if strcmp(vModel,'iasp91')
            ccp410{i,j} = zeros(8001,1);
        elseif strcmp(vModel,'LLNL')
            ccp410{i,j} = zeros(1,7501);
        end
        % fprintf('[%i %i] N = %i\n',i,j,length(idx))
        for k = 1:length(idx)
            rfData = file{3}{idx(k)};
            tmp = strsplit(rfData,'TA.');
            stat = extractBefore(tmp{2},'.');
            
            if strcmp(vModel,'iasp91')
                % Read in a receiver function
                rfPath = fullfile(rfDir,stat,rfDet,rfData);
                [rf{k}.t,rf{k}.d,rf{k}.h] = fread_sac(rfPath);

                % Get depth conversion data
                P = rf{k}.h.user(10)/deg2km(1);
                dt = rf{k}.h.delta;

                % Depth convert the receiver function
                rf{k}.dep = rfDepcon(rf{k}.d,dt,P,dz,'iasp91','true');

                % Add the depth converted receiver function to the stack
                ccp410{i,j} = ccp410{i,j} + rf{k}.dep;
            elseif strcmp(vModel,'LLNL')
                % Path to stacked LLNL data
                rfData = [extractBefore(rfData,'SAC'),'LLNL3D.mat'];
                
                % Read in the stacked receiver function
                rfPath = fullfile(rfDir,stat,rfDet,rfData);
                load(rfPath);
                
                ccp410{i,j} = ccp410{i,j} + rfDepth(1:7501);
            end
            
        end
        
        % Normalize the stack by the number of RFs in the bin
        ccp410{i,j} = ccp410{i,j}/length(idx);
    end
end

% Stack 660 bits
nRows = size(N{4},1);
nCols = size(N{4},2);
ccp660 = cell(nRows,nCols);
for i = 1:nRows
    for j = 1:nCols
        sbin = [i j];
        [tf,loc] = ismember([binX{4},binY{4}],sbin,'row');
        idx = find(loc);
        
        % Read in all of the receiver functions in this bin
        rf = {};
        if strcmp(vModel,'iasp91')
            ccp660{i,j} = zeros(8001,1);
        elseif strcmp(vModel,'LLNL')
            ccp660{i,j} = zeros(1,7501);
        end
        % fprintf('[%i %i] N = %i\n',i,j,length(idx))
        for k = 1:length(idx)
            rfData = file{4}{idx(k)};
            tmp = strsplit(rfData,'TA.');
            stat = extractBefore(tmp{2},'.');
            
            if strcmp(vModel,'iasp91')
                % Read in a receiver function
                rfPath = fullfile(rfDir,stat,rfDet,rfData);
                [rf{k}.t,rf{k}.d,rf{k}.h] = fread_sac(rfPath);

                % Get depth conversion data
                P = rf{k}.h.user(10)/deg2km(1);
                dt = rf{k}.h.delta;

                % Depth convert the receiver function
                rf{k}.dep = rfDepcon(rf{k}.d,dt,P,dz,'iasp91','true');

                % Add the depth converted receiver function to the stack
                ccp660{i,j} = ccp660{i,j} + rf{k}.dep;
                
            elseif strcmp(vModel,'LLNL')
                % Path to stacked LLNL data
                rfData = [extractBefore(rfData,'SAC'),'LLNL3D.mat'];
                
                % Read in the stacked receiver function
                rfPath = fullfile(rfDir,stat,rfDet,rfData);
                load(rfPath);
                
                ccp660{i,j} = ccp660{i,j} + rfDepth(1:7501);
                
            end
            
        end
        
        % Normalize the stack by the number of RFs in the bin
        ccp660{i,j} = ccp660{i,j}/length(idx);
    end
end

% Stitch together the 410 and 660 portions
for i = 1:nCols
    for j = 1:nRows
        if strcmp(vModel,'iasp91')
            CCP(i,j,:) = zeros(8001,1);
        elseif strcmp(vModel,'LLNL')
            CCP(i,j,:) = zeros(1,7501);
        end
        CCP(i,j,2851:5350) = ccp410{j,i}(2851:5350);
        CCP(i,j,5351:7850) = ccp660{j,i}(5351:7850);
    end
end

%--------------------------------------------------------------------------
% Save the result
%--------------------------------------------------------------------------

save('CCP_Stack_Filt0.02_0.2_GW10.mat','CCP','ccp410','ccp660',...
        'latBin','lonBin');
