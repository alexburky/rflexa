% figure4.m
%
% This code is meant to serve as a companion to the 2021 EduQuakes paper,
% 'Instrument response removal and the 2020 M3.1 Marlboro, New Jersey,
% earthquake.' It produces the fourth figure of the paper.
%
% External dependencies:
% - fread_sac
% - parseRESP
% - parsePZ
% - transfer
% - transfer_evalresp
%
%--------------------------------------------------------------------------
% Last updated 6/8/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% Specify the locations of the seismic data (traces)
sacDir = '../data/traces/';
data{1} = 'PP.S0001.00.HHZ.D.2020.253.060000.SAC';

% Specify the locations of the response files
respDir = '../data/responses/';
resp{1} = 'PP.S0001.00.HHZ.resp';

% Read in the seismograms and responses
for i = 1:length(data)
    [s{i}.t,s{i}.d,s{i}.h] = fread_sac(fullfile(sacDir,data{i}));
    [z{i},p{i},k{i}] = parseRESP(fullfile(respDir,resp{i}));
end

% Read in data that has been corrected in SAC
sac1 = 'PP.S0001.00.HHZ.D.2020.253.060000.VEL.SACPZ.SAC';
sac2 = 'PP.S0001.00.HHZ.D.2020.253.060000.VEL.EVALRESP.SAC';

[sac{1}.t,sac{1}.d,sac{1}.h] = fread_sac(fullfile(sacDir,sac1));
[sac{2}.t,sac{2}.d,sac{2}.h] = fread_sac(fullfile(sacDir,sac2));

% Path to files output by EVALRESP
respDir = '../data/responses/';
respAmp = 'AMP.PP.S0001.00.HHZ';
respPhase = 'PHASE.PP.S0001.00.HHZ';
aFile = fullfile(respDir,respAmp);
pFile = fullfile(respDir,respPhase);

%--------------------------------------------------------------------------
% Instrument response removal
%--------------------------------------------------------------------------

% Transfer parameters
delta = 0.01;
flims = [0.1 0.2 10 20];

% Desired output units
units = 'velocity';

% Remove the instrument response from S0001 (SACPZ)
s1resp = fullfile(respDir,resp{1});
s{1}.d_sacpz = s{1}.d - mean(s{1}.d);
s{1}.d_sacpz = transfer(s{1}.d_sacpz,delta,flims,units,s1resp,'resp');

% Remove the instrument response from S0001 (EVALRESP)
s{1}.d_evresp = s{1}.d - mean(s{1}.d);
s{1}.d_evresp = transfer_evalresp(s{1}.d_evresp,1/100,flims,aFile,pFile);

%--------------------------------------------------------------------------
% FIGURE 4
%--------------------------------------------------------------------------

% Plot parameters
xmin = 10;
xmax = 60;
ymin = -7e4;
ymax = 7e4;

color{1} = [0.25 0.5 0.15];
color{2} = [0 0 0];

% Subplot 1: SACPZ Records
subplot(4,1,1)
plot(s{1}.t,s{1}.d_sacpz,'Color',color{1})
hold on
plot(sac{1}.t,sac{1}.d,'Color',color{2})
legend('rflexa: \texttt{transfer from polezero}',...
    'SAC: \hspace{0.1cm} \texttt{transfer from polezero}')
grid on
xlim([xmin xmax])
ylim([-2e-4 2e-4])
ax1 = gca;
ax1.TickDir = 'out';
ax1.XTickLabel = {};
ax1.XTick = 10:2.5:60;
ax1.YTick = [-2e-4 -1e-4 0 1e-4 2e-4];
ax1.YTickLabel = {'-2e-4','-1e-4','0','1e-4','2e-4'};
rectangle('Position',[10.6 -1.85e-4 2.15 1.35e-4],'FaceColor',[1 1 1])
text(10.85,-1.2e-4,'(a)','FontSize',12)

% Subplot 2: SACPZ Difference
subplot(4,1,2)
plot(s{1}.t,s{1}.d_sacpz - sac{1}.d,'Color',[1 0 0])
legend('Difference')
grid on
xlim([xmin xmax])
ylim([-8e-11 8e-11])
ax2 = gca;
ax2.TickDir = 'out';
ax2.XTickLabel = {};
ax2.XTick = 10:2.5:60;
ax2.YTick = [-8e-11 -4e-11 0 4e-11 8e-11];
ax2.YTickLabel = {'-8e-11','-4e-11','0','4e-11','8e-11'};
ylabel('Velocity (m/s)')
ax2.YLabel.Position(2) = -1.2e-10;
rectangle('Position',[10.6 -7e-11 2.15 6e-11],'FaceColor',[1 1 1])
text(10.85,-4.25e-11,'(b)','FontSize',12)

% Subplot 3: EVALRESP Records
subplot(4,1,3)
plot(s{1}.t,s{1}.d_evresp,'Color',color{1})
hold on
plot(sac{2}.t,sac{2}.d./1e9,'Color',color{2})
legend('rflexa: \texttt{transfer from evalresp}',...
    'SAC: \hspace{0.1cm} \texttt{transfer from evalresp}')
grid on
xlim([xmin xmax])
ylim([-2e-4 2e-4])
ax3 = gca;
ax3.TickDir = 'out';
ax3.XTickLabel = {};
ax3.XTick = 10:2.5:60;
ax3.YTick = [-2e-4 -1e-4 0 1e-4 2e-4];
ax3.YTickLabel = {'-2e-4','-1e-4','0','1e-4','2e-4'};
rectangle('Position',[10.6 -1.85e-4 2.15 1.35e-4],'FaceColor',[1 1 1])
text(10.85,-1.2e-4,'(c)','FontSize',12)

% Subplot 4: EVALRESP Difference
subplot(4,1,4)
plot(s{1}.t,s{1}.d_evresp - transpose(sac{2}.d)./1e9,'Color',[1 0 0])
legend('Difference')
grid on
xlim([xmin xmax])
ylim([-8e-12 8e-12])
ax4 = gca;
ax4.TickDir = 'out';
ax4.YTick = [-8e-12 -4e-12 0 4e-12 8e-12];
ax4.YTickLabel = {'-8e-12','-4e-12','0','4e-12','8e-12'};
ax4.XTick = 10:2.5:60;
ax4.XTickLabel = {'','','06:00:15','','','','06:00:25','','','',...
                  '06:00:35','','','','06:00:45','','','','06:00:55',...
                  '',''};
xlabel('Time (UTC)')
rectangle('Position',[10.6 -7e-12 2.15 6e-12],'FaceColor',[1 1 1])
text(10.85,-4.5e-12,'(d)','FontSize',12)

% Final figure formatting
set(gcf,'Position',[0 0 600 300]);
sgtitle('rflexa \texttt{transfer} vs. SAC \texttt{transfer}')

