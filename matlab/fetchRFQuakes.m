% fetchRFQuakes.m

% This script fetches earthquake data for receiver function analysis,
% and additionally fetches poles and zeros data for the selected
% stations and channels. The earthquake data are saved as .SAC files,
% and their filenames are formatted so that they can be easily linked
% to their corresponding poles and zeros data.
%
%--------------------------------------------------------------------------
% Last updated 12/05/2020 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Next step: Write a function which removes the instrument response
%            entirely within MATLAB. Inject this function into a receiver
%            function workflow.

% Can only get the true SACPZ data from a trace object, so...
% For each channel, retrieve a single trace object and construct
% the pz file from the trace object AND the channel object

clear,clc

% Define the directory where you would like to save the data
sacDir = '/Users/aburky/PycharmProjects/bermudaRFs/matlab/';

% Define the network, station that you would like to fetch data for
network = 'TA';
station = 'O61A';
location = '*';
channel = 'BH*';

% network = 'II';
% station = 'SACV';
% location = '10';
% channel = 'BHZ';

% Define desired earthquake parameters
minMag = 7.0;
maxMag = 7.5;
minRad = 30;
maxRad = 90;

ch = irisFetch.Channels('RESPONSE',network,station,location,channel);

timeFormat = 'yyyy-mm-dd HH:MM:SS.FFF';

% Loop over each channel and make the PZ file
for i = 1:length(ch)
    % Get the channel start and end dates
    t1 = ch(i).StartDate;
    t2 = ch(i).EndDate;
    if isempty(t2)
        dtLocal = datetime('now','TimeZone','Local');
        t2 = datestr(datetime(dtLocal,'TimeZone','Z'),...
            'yyyy-mm-dd HH:MM:SS.FFF');
    end
    
    % Get the PZ data by fetching a trace object
    
    % Need to recursively call a function here, check if there is trace
    % data within the channel operation period, if not, iterate time
    startChan = t1;
    endChan = datetime(t1) + hours(0.5);
    endChan = datestr(endChan,timeFormat);
    looping = 1;
    while looping == 1
        tr = irisFetch.Traces(network,station,location,...
                ch(i).ChannelCode,startChan,endChan,'includePZ');
        if isempty(tr)
            startChan = datetime(startChan) + days(1);
            startChan = datestr(startChan,timeFormat);
            endChan = datetime(startChan) + hours(0.5);
            endChan = datestr(endChan,timeFormat);
        else
            pz(i) = tr.sacpz;
            looping = 0;
        end
    end

    % pz(i) = tr.sacpz;
    
    % Save the PZ file
    savePZ(ch(i),pz(i),sacDir,'pzindex',i);
    % While we are iterating over the channel, check for earthquakes
    % that meet our search criterion during its operation
    donut = [ch(i).Latitude,ch(i).Longitude,maxRad,minRad];
    ev = irisFetch.Events('MinimumMagnitude',minMag,'MaximumMagnitude',...
            maxMag,'radialcoordinates',donut,'startTime',t1,...
            'endTime',t2);
        
    % Loop over each event, get the trace for the channel, and save
    for j = 1:length(ev)
        % Get start time of event
        ev_start = ev(j).PreferredTime;
        ev_end = datetime(ev_start) + hours(1);
        ev_end = datestr(ev_end,timeFormat);
        % Fetch trace data from the current channel
        tr = irisFetch.Traces(network,station,location,...
                ch(i).ChannelCode,ev_start,ev_end);
        % Save the trace data to a SAC file
        if ~isempty(tr)
            saveSAC(tr,ev_start,sacDir,'event',ev(j),'pz',i);
        end
    end
end

% Options to add:
% - RESP files? (saveRESP function)
% Note:
% Working with the RESP files is a bit hairy, for negligible difference
% in results. For the time being, do everything with PZ files.

% % Loop over each channel and get the RESP file
% for i = 1:length(ch)
%     t1 = ch(i).StartDate;
%     t2 = datetime(t1) + milliseconds(1);
%     formatOut = 'yyyy-mm-dd HH:MM:SS.FFF';
%     t2 = datestr(t2,formatOut);
%     % Fetch the RESP file
%     re = irisFetch.Resp(network,station,location,ch(i).ChannelCode,t1,t2);
%     % Format the RESP filename
%     respFile = sprintf('RESP.%s.%s.%s.%s.%d',network,station,...
%                ch(i).LocationCode,ch(i).ChannelCode,i);
%     % Save the RESP file
%     fID = fopen(respFile,'w');
%     fprintf(fID,re);
%     fclose(fID);
% end

%%

t1 = ch(1).StartDate;
% Use the start date plus a very short time as the time window
t2 = '2014-08-28 20:00:00.000';
% t2 = ch(1).EndDate;

channel = 'BHZ';

re = irisFetch.Resp(network,station,location,channel,t1,t2);

% Loop over all channels and save a RESP file for each one!

R = regexp(re,'\n','split');
% Get the number of poles and the number of zeroes
nzLine = find(contains(R,'Number of zeroes:'));
tmp = split(R{nzLine});
nz = str2double(tmp{end});
npLine = find(contains(R,'Number of poles:'));
tmp = split(R{npLine});
np = str2double(tmp{end});

% Get the line where the zeroes begin
zLine = find(contains(R,'Complex zeroes'));
for i = 2:(nz + 1)
    tmp = split(R{zLine + i});
    z(i-1) = complex(str2double(tmp{3}),str2double(tmp{4}));
end

% Get the line where the poles begin
pLine = find(contains(R,'Complex poles'));
for i = 2:(np + 1)
    tmp = split(R{pLine + i});
    p(i-1) = complex(str2double(tmp{3}),str2double(tmp{4}));
end

% Get the poles and zeroes constant (A0 * sensitivity)
a0Line = find(contains(R,'A0 normalization factor'));
tmp = split(R{a0Line});
a0 = str2double(tmp{end});
sensLine = find(contains(R,'Sensitivity'));
tmp = split(R{sensLine(end)});
sens = str2double(tmp{end});
k = a0 * sens;

% Try saving the instrument response information to a file
fID = fopen('RESP.TEST','w');
fprintf(fID,re);
fclose(fID);

% End date of one channel is equal to the start date of the next.
% what to do about this?

%% Make Bode plots

[z,p,k] = parsePZ('/Users/aburky/PycharmProjects/bermudaRFs/matlab/SAC_PZs_TA_O61A_BHZ.6');

npts = 10000000;
nfft = 2^nextpow2(npts);
nfreq = (nfft / 2) + 1;
f = linspace(0,40,nfreq);

[b,a] = zp2tf(z,p,k);
[h,w] = freqs(b,a,2*pi*f);

% Amplitude response plot
figure(1)
loglog(f,abs(h),'r','linewidth',1)
hold on
plot([20 20],[1e0 1e10],'k--')
xlim([1e-5 40])
ylim([5e3 3e9])
grid on
title('TA.O61A..BHZ: Amplitude Response')
ylabel('Gain (Counts / m/s)')
xlabel('Frequency (Hz)')
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
set(gcf,'Position',[0 0 600 400])

% Phase response plot
figure(2)
semilogx(f,((180/pi) * angle(h)),'b','linewidth',1)
hold on
plot([20 20],[-200 200],'k--')
xlim([1e-5 40])
ylim([-181 180])
grid on
title('TA.061A..BHZ: Phase Response')
ylabel('Phase ($^{\circ}$)')
xlabel('Frequency (Hz)')
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [-180, -135, -90, -45, 0, 45, 90, 135, 180];
set(gcf,'Position',[650 0 600 400])