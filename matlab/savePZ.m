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
% Last updated 12/01/2020 by aburky@princeton.edu
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

% Format time strings
tStart = datestr(datetime(channel.StartDate),'yyyy-mm-ddTHH:MM:SS');
if ~isempty(channel.EndDate)
    tEnd = datestr(datetime(channel.EndDate),'yyyy-mm-ddTHH:MM:SS');
else
    tEnd = '';
end
dtLocal = datetime('now','TimeZone','Local');
dtUTC = datestr(datetime(dtLocal,'TimeZone','Z'),'yyyy-mm-ddTHH:MM:SS');

% Write the header (emulate RDSEED v5.2)
fprintf(fID,sprintf('* **********************************\n'));
fprintf(fID,sprintf('* NETWORK   (KNETWK): %s\n',channel.NetworkCode));
fprintf(fID,sprintf('* STATION    (KSTNM): %s\n',channel.StationCode));
fprintf(fID,sprintf('* LOCATION   (KHOLE): %s\n',channel.LocationCode));
fprintf(fID,sprintf('* CHANNEL   (KCMPNM): %s\n',channel.ChannelCode));
fprintf(fID,sprintf('* CREATED           : %s\n',dtUTC));
fprintf(fID,sprintf('* START             : %s\n',tStart));
fprintf(fID,sprintf('* END               : %s\n',tEnd));
fprintf(fID,sprintf('* DESCRIPTION       : %s\n',channel.StationName));
fprintf(fID,sprintf('* LATITUDE          : %0.6f\n',channel.Latitude));
fprintf(fID,sprintf('* LONGITUDE         : %0.6f\n',channel.Longitude));
fprintf(fID,sprintf('* ELEVATION         : %0.1f\n',channel.Elevation));
fprintf(fID,sprintf('* DEPTH             : %0.1f\n',channel.Depth));
fprintf(fID,sprintf('* DIP               : %0.1f\n',abs(channel.Dip)-90));
fprintf(fID,sprintf('* AZIMUTH           : %0.1f\n',channel.Azimuth));
fprintf(fID,sprintf('* SAMPLE RATE       : %0.1f\n',channel.SampleRate));
fprintf(fID,sprintf('* INPUT UNIT        : M\n'));
fprintf(fID,sprintf('* OUTPUT UNIT       : COUNTS\n'));
fprintf(fID,sprintf('* INSTTYPE          : %s\n',...
    channel.Sensor.Description));
fprintf(fID,sprintf('* INSTGAIN          : %0.6e (M/S)\n',...
    channel.Response.Stage(1).StageGain.Value));
fprintf(fID,sprintf('* COMMENT           : \n'));
fprintf(fID,sprintf('* SENSITIVITY       : %0.6e (M/S)\n',...
    channel.Response.InstrumentSensitivity.Value));
fprintf(fID,sprintf('* A0                : %0.6e \n',...
    channel.Response.Stage(1).PolesZeros.NormalizationFactor));
fprintf(fID,sprintf('* **********************************\n'));

% Save the poles, zeros, and constant
nzeros = 3;
z = nonzeros(z);
% Use the IRIS convention of not saving zeros equal to (0,0)
fprintf(fID,sprintf('ZEROS %d\n',length(z) + nzeros));
fprintf(fID,sprintf('%+e %+e\n',0.0,0.0));
fprintf(fID,sprintf('%+e %+e\n',0.0,0.0));
fprintf(fID,sprintf('%+e %+e\n',0.0,0.0));
for i = 1:length(z)
    fprintf(fID,sprintf('%+e %+e\n',real(z(i)),imag(z(i))));
end
fprintf(fID,sprintf('POLES %d\n',length(p)));
for i = 1:length(p)
    fprintf(fID,sprintf('%+e %+e\n',real(p(i)),imag(p(i))));
end
fprintf(fID,sprintf('CONSTANT %e',k));
fclose(fID);
