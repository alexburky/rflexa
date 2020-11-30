function savePZ(channel,directory,varargin)
% SAVEPZ Constructs and saves a SAC_PZs poles and zeros file from a
%        channel structure.
%
% >> savePZ(channel,directory)
%
%---Input Variables--------------------------------------------------------
% channel   - channel structure containing all of the relevant poles and
%             zeros information (specifically, one output by 
%             irisFetch.Channels('RESPONSE'))
% directory - string of the folder where you would like to save the file
%
% Optional Variables:
% pz        - index to append to PZ filename to indicate relationship to
%             corresponding SAC files
%
%--------------------------------------------------------------------------
% Last updated 11/30/2019 by aburky@princeton.edu
%--------------------------------------------------------------------------

% To do: Add some error handling...

% Get the poles, zeros, and constant
z = channel.Response.Stage(1).PolesZeros.Zero;
p = channel.Response.Stage(1).PolesZeros.Pole;
k = double(channel.Response.Stage(1).PolesZeros.NormalizationFactor)*...
    double(channel.Response.InstrumentSensitivity.Value);

% Format the PZ file name
if isempty(channel.LocationCode)
    pzFile = sprintf('SAC_PZs_%s_%s_%s',channel.NetworkCode,...
             channel.StationCode,channel.ChannelCode);
else
    pzFile = sprintf('SAC_PZs_%s_%s_%s_%s',channel.NetworkCode,...
             channel.StationCode,channel.ChannelCode,channel.LocationCode);
end
if nargin > 2
    if strcmp(varargin{1},'pz')
        pzIndex = varargin{2};
        pzFile = sprintf('%s.%i',pzFile,pzIndex);
    end
end

% Save the data to the PZ file
fID = fopen(fullfile(directory,pzFile),'w');

% Write the start and end dates
fprintf(fID,sprintf('* Start date: %s\n',channel.StartDate));
fprintf(fID,sprintf('* End date: %s\n',channel.EndDate));

% Save the poles, zeros, and constant
fprintf(fID,sprintf('ZEROS %d\n',length(z)));
for i = 1:length(z)
    fprintf(fID,sprintf('%+e %+e\n',real(z(i)),imag(z(i))));
end
fprintf(fID,sprintf('POLES %d\n',length(p)));
for i = 1:length(p)
    fprintf(fID,sprintf('%+e %+e\n',real(p(i)),imag(p(i))));
end
fprintf(fID,sprintf('CONSTANT %e',k));
fclose(fID);
