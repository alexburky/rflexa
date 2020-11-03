% fetchRFQuakes.m

% This script fetches earthquake data for receiver function analysis,
% and additionally fetches poles and zeros data for the selected
% stations and channels. The earthquake data are saved as .SAC files,
% and their filenames are formatted so that they can be easily linked
% to their corresponding poles and zeros data.
%
%--------------------------------------------------------------------------
% Last updated 11/03/2020 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% Define the network, station that you would like to fetch data for
network = 'TA';
station = 'O61A';
% network = 'II';
% station = 'SACV';
location = '*';
channel = 'BH*';

ch = irisFetch.Channels('RESPONSE',network,station,location,channel);

% Options to add:
% - Directory where you would like to save the output data
% - RESP, or PZ files?

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

% Loop over each channel and make the PZ file
for i = 1:length(ch)
    % Get the channel start and end dates
    t1 = ch(i).StartDate;
    t2 = ch(i).EndDate;
    % Get the poles, zeros, and constant
    z = ch(i).Response.Stage(1).PolesZeros.Zero;
    p = ch(i).Response.Stage(1).PolesZeros.Pole;
    k = double(ch(i).Response.Stage(1).PolesZeros.NormalizationFactor)*...
        double(ch(i).Response.InstrumentSensitivity.Value);
    % Format the PZ filename
    if isempty(ch(i).LocationCode)
        pzFile = sprintf('SAC_PZs_%s_%s_%s.%d',network,station,...
                 ch(i).ChannelCode,i);
    else
        pzFile = sprintf('SAC_PZs_%s_%s_%s_%s.%d',network,station,...
                 ch(i).ChannelCode,ch(i).LocationCode,i);
    end
    % Save data to the PZ file
    fID = fopen(pzFile,'w');
    % Write the start and end dates
    fprintf(fID,sprintf('* Start date: %s\n',t1));
    fprintf(fID,sprintf('* End date:   %s\n',t2));
    % Save the poles, zeros, and constant
    fprintf(fID,sprintf('ZEROS %d\n',length(z)));
    for j = 1:length(z)
        fprintf(fID,sprintf('%+e %+e\n',real(z(j)),imag(z(j))));
    end
    fprintf(fID,sprintf('POLES %d\n',length(p)));
    for j = 1:length(p)
        fprintf(fID,sprintf('%+e %+e\n',real(p(j)),imag(p(j))));
    end
    fprintf(fID,sprintf('CONSTANT %e',k));
    fclose(fID);
end

% Once you have the instrument response files, get some earthquake data!


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

%% Sort out some date time logic here

t2 = datetime(t1) + milliseconds(1);
formatOut = 'yyyy-mm-dd HH:MM:SS.FFF';
t2 = datestr(t2,formatOut);

%% Calculate the poles zeros constant?

z = ch(1).Response.Stage(1).PolesZeros.Zero;
p = ch(1).Response.Stage(1).PolesZeros.Pole;
k = double(ch(1).Response.Stage(1).PolesZeros.NormalizationFactor) * ...
    double(ch(1).Response.InstrumentSensitivity.Value);

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