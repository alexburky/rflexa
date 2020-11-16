function saveSAC(trace,startTime,directory)
% SAVESAC Saves trace data to a SAC file and fills all of the appropriate
%         header variables.
%
% >> saveSAC(trace,startTime,directory)
%
%---Input Variables--------------------------------------------------------
% trace     - trace structure containing data and station information
%             (specifcally, the output of irisFetch.Traces)
% startTime - date string indicated the start time of the data
%             formatted as: 'yyyy-mm-dd HH:MM:SS.FFF'
% directory - string of the folder where you would like to save the data
%
%--------------------------------------------------------------------------
% Last updated 11/15/2019 by aburky@princeton.edu
%--------------------------------------------------------------------------

% TO DO: add the option to also save event information a la receiver
%        function workflow

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

% Format the SAC filename
fname = sprintf('%i.%02d.%02d.%02d.%02d.%02d.%s.%s.%s.%s.SAC',...
                sacmat.hdr.nzyear,month(datetime(startTime)),...
                day(datetime(startTime)),sacmat.hdr.nzhour,...
                sacmat.hdr.nzmin,sacmat.hdr.nzsec,trace.network,...
                trace.station,trace.location,trace.channel);
   
% Save the data to a SAC file
fwrite_sac(sacmat,fullfile(directory,fname));
    