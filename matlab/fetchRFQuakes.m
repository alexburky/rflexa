% fetchRFQuakes.m

clear,clc

% Define the network, station that you would like to fetch data for
network = 'TA';
station = 'O61A';
% network = 'IU';
% station = 'BBSR';
location = '*';
channel = 'BH*';

ch = irisFetch.Channels('RESPONSE',network,station,location,channel);

% So, it looks like you can't get a RESP file for a station that has
% been updated without getting the RESP information for all prior
% periods that station was operational.
% Need to come up with a way to deal with this - looks like we will need
% to use some datestr manipulation!
% Just save the one RESP file that has it all...

% Loop over all channels, and then find the latest operational date for
% each channel (BHE,BHN,BHZ,etc.) - thenn fetch and save the RESP file
% for this time period!

% Loop over each channel and get the RESP file
for i = 1:length(ch)
    t1 = ch(i).StartDate;
    t2 = datetime(t1) + milliseconds(1);
    formatOut = 'yyyy-mm-dd HH:MM:SS.FFF';
    t2 = datestr(t2,formatOut);
    % Fetch the RESP file
    re = irisFetch.Resp(network,station,location,ch(i).ChannelCode,t1,t2);
    % Format the RESP filename
    respFile = sprintf('RESP.%s.%s.%s.%s.%d',network,station,...
               ch(i).LocationCode,ch(i).ChannelCode,i);
    % Save the RESP file
    fID = fopen(respFile,'w');
    fprintf(fID,re);
    fclose(fID);
end

% Loop over each channel and make the PZ file
for i = 1:length(ch)
    % t1 = ch(i).StartDate;
    % t2 = datetime(t1) + milliseconds(1);
    % formatOut = 'yyyy-mm-dd HH:MM:SS.FFF';
    % t2 = datestr(t2,formatOut);
    % Get the poles, zeros, and constant
    z = ch(i).Response.Stage(1).PolesZeros.Zero;
    p = ch(i).Response.Stage(1).PolesZeros.Pole;
    k = double(ch(i).Response.Stage(1).PolesZeros.NormalizationFactor)*...
        double(ch(i).Response.InstrumentSensitivity.Value);
    % Format the PZ filename
    pzFile = sprintf('SAC_PZs_%s_%s_%s.%d',network,station,...
             ch(i).ChannelCode,i);
    % Save the PZ file
    fID = fopen(pzFile,'w');
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
