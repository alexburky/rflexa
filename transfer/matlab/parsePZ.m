function [z,p,k] = parsePZ(pzfile)
% PARSEPZ Parse a SAC style Pole Zero file and get the poles, zeros, and
%         constant for removing the instrument response
%
% >> [z,p,k] = parsePZ(pzfile)
%
%---Input Variables--------------------------------------------------------
% pzfile - Full path to the SAC_PZs file
% 
%---Output Variables-------------------------------------------------------
% z      - Vector containing the zeros
% p      - Vector containing the poles
% k      - Value of the constant in the transfer function
% 
%--------------------------------------------------------------------------
% Last updated 10/27/2020 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Find lines delimiting zeros, poles, and constant entries
F = regexp(fileread(pzfile),'\n','split');
zLine = find(contains(F,'ZEROS'));
pLine = find(contains(F,'POLES'));
cLine = find(contains(F,'CONSTANT'));

lineNum = 0;
zidx = 1;
pidx = 1;
fileID = fopen(pzfile,'r');

% Maybe add a check here - does the number of Zeros on the line
% that says 'ZEROS' match the length of the zeros array? If not,
% need to add the zeros?

% Loop over lines of file and fill vectors
while ~feof(fileID)
    lineNum = lineNum + 1;
    line = fgetl(fileID);
    if (lineNum > zLine && lineNum < pLine)
        zn = textscan(line,'%f');
        z(zidx) = complex(zn{1}(1),zn{1}(2));
        zidx = zidx + 1;
    elseif (lineNum > pLine && lineNum < cLine)
        pn = textscan(line,'%f');
        p(pidx) = complex(pn{1}(1),pn{1}(2));
        pidx = pidx + 1;
    elseif (lineNum == cLine)
        cn = textscan(line,'%s %f');
        k = cn{2};
    end
end

z = transpose(z);
p = transpose(p);

fclose(fileID);