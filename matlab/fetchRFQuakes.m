% fetchRFQuakes.m

clear,clc

% Define the network, station that you would like to fetch data for
network = 'TA';
station = 'O61A';
location = '*';
channel = 'BH*';

ch = irisFetch.Channels('RESPONSE',network,station,location,channel);

t1 = ch(1).StartDate;
% Use the start date plus a very short time as the time window
t2 = '2014-08-28 20:00:00.000';
% t2 = ch(1).EndDate;

channel = 'BHZ';

re = irisFetch.Resp(network,station,location,channel,t1,t2);

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

