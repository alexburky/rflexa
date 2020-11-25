function saveSAC(trace,startTime,directory,varargin)
% SAVESAC Saves trace data to a SAC file and fills all of the appropriate
%         header variables.
%
% >> saveSAC(trace,startTime,directory)
%
%---Input Variables--------------------------------------------------------
% trace     - trace structure containing data and station information
%             (specifcally, the output of irisFetch.Traces)
% startTime - date string indicating the start time of the data
%             formatted as: 'yyyy-mm-dd HH:MM:SS.FFF'
% directory - string of the folder where you would like to save the data
% 
% Optional Variables:
% event     - event structure containing information about the
%             earthquake that seismic data potentially pertains to
%
%--------------------------------------------------------------------------
% Last updated 11/25/2019 by aburky@princeton.edu
%--------------------------------------------------------------------------

% SAC timeseries data
sacmat.data = trace.data;
sacmat.hdr.delta = 1/trace.sampleRate;

% SAC station data
sacmat.hdr.cmpaz = trace.azimuth;
sacmat.hdr.cmpinc = trace.dip;
sacmat.hdr.stla = trace.latitude;
sacmat.hdr.stlo = trace.longitude;
sacmat.hdr.stel = trace.elevation;
sacmat.hdr.stdp = trace.depth;
sacmat.hdr.kcmpnm = trace.channel;
sacmat.hdr.knetwk = trace.network;
sacmat.hdr.kstnm = trace.station;
sacmat.hdr.khole = trace.location;

% SAC data calendar time information
sacmat.hdr.nzyear = year(datetime(startTime));
sacmat.hdr.nzjday = day(datetime(startTime),'dayofyear');
sacmat.hdr.nzhour = hour(datetime(startTime));
sacmat.hdr.nzmin = minute(datetime(startTime));
dseconds = second(datetime(startTime));
if dseconds == 0
    sacmat.hdr.nzsec = 0;
    sacmat.hdr.nzmsec = 0;
else
    sacmat.hdr.nzsec = str2double(extractBefore(num2str(dseconds),'.'));
    sacmat.hdr.nzmsec = str2double(extractAfter(num2str(dseconds),'.'));
end

% The user has entered optional variables
if nargin > 3
    for i = 1:length(varargin)
        % User has input an event data structure
        if strcmp(varargin{i},'event')
            event = varargin{i+1};
            sacmat.hdr.mag = event.PreferredMagnitudeValue;
            sacmat.hdr.evla = event.PreferredLatitude;
            sacmat.hdr.evlo = event.PreferredLongitude;
            sacmat.hdr.evdp = event.PreferredDepth;
            % Event - station derived data
            sacmat.hdr.baz = azimuth(sacmat.hdr.stla,sacmat.hdr.stlo,...
            sacmat.hdr.evla,sacmat.hdr.evlo);
            [sacmat.hdr.gcarc,sacmat.hdr.az] = distance(sacmat.hdr.evla,...
            sacmat.hdr.evlo,sacmat.hdr.stla,sacmat.hdr.stlo);
            sacmat.hdr.dist = deg2km(sacmat.hdr.gcarc);
        % User wants to save PZ Index
        elseif strcmp(varargin{i},'pz')
            pzIndex = varargin{i+1};
        end
    end
end

% Format the SAC filename
fname = sprintf('%i.%02d.%02d.%02d.%02d.%02d.%s.%s.%s.%s.%i.SAC',...
                sacmat.hdr.nzyear,month(datetime(startTime)),...
                day(datetime(startTime)),sacmat.hdr.nzhour,...
                sacmat.hdr.nzmin,sacmat.hdr.nzsec,trace.network,...
                trace.station,trace.location,trace.channel,pzIndex);
   
% Save the data to a SAC file
fwrite_sac(sacmat,fullfile(directory,fname));
    