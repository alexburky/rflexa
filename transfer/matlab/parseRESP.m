function [z,p,k] = parseRESP(respFile)
% PARSERESP Parse a SAC style RESP file and get the poles, zeros, and
%           constant for removing the instrument response
%
% >> [z,p,k] = parseResp(respFile)
%
%---Input Variables--------------------------------------------------------
% respFile - Full path to the RESP file
% 
%---Output Variables-------------------------------------------------------
% z      - Vector containing the zeros
% p      - Vector containing the poles
% k      - Value of the constant in the transfer function
% 
%--------------------------------------------------------------------------
% Last updated 3/2/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Find lines delimiting zeros, poles, and constant entries
F = regexp(fileread(respFile),'\n','split');
nZero = find(contains(F,'B053F09'));
zLine = find(contains(F,'B053F10-13'));
nPole = find(contains(F,'B053F14'));
pLine = find(contains(F,'B053F15-18'));
aLine = find(contains(F,'B053F07'));
sLine = find(contains(F,'Sensitivity:'));

% Initiliaze iterators
lineNum = 0;
nzidx = 1;
zidx = 1;
npidx = 1;
pidx = 1;

fileID = fopen(respFile,'r');
% First, determine how many poles and zeros there are
while ~feof(fileID)
    lineNum = lineNum + 1;
    line = fgetl(fileID);
    if ismember(lineNum,nZero)
        nzl = split(line);
        nz(nzidx) = str2double(nzl{5});
        nzidx = nzidx + 1;
    elseif ismember(lineNum,nPole)
        npl = split(line);
        np(npidx) = str2double(npl{5});
        npidx = npidx + 1;
    end
end

nzeros = max(nz);
npoles = max(np);

lineNum = 0;
fileID = fopen(respFile,'r');
% Next, read in the poles, zeros, A0, and sensitivity
while ~feof(fileID)
    lineNum = lineNum + 1;
    line = fgetl(fileID);
    if ismember(lineNum,zLine)
        zn = split(line);
        z(zidx) = complex(str2double(zn{3}),str2double(zn{4}));
        zidx = zidx + 1;
    elseif ismember(lineNum,pLine)
        pn = split(line);
        p(pidx) = complex(str2double(pn{3}),str2double(pn{4}));
        pidx = pidx + 1;
    % This currently takes the first value of A0 as the correct one, check
    % and make sure this is always the case!
    elseif ismember(lineNum,aLine(1))
        an = split(line);
        a0 = str2double(an{5});
    elseif ismember(lineNum,sLine)
        sn = split(line);
        sens = str2double(sn{3});
    end
end

% Zeros, poles, and constant
z = z(1:nzeros);
z = transpose(z);
p = p(1:npoles);
p = transpose(p);
k = a0*sens;
