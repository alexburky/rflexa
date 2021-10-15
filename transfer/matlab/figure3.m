% figure3.m
%
% This code is meant to serve as a companion to the 2021 EduQuakes paper,
% 'Instrument response removal and the 2020 M3.1 Marlboro, New Jersey,
% earthquake.' It produces the third figure of the paper.
%
% External dependencies:
% - fread_sac
% - parseRESP
% - parsePZ
% - transfer
%
%--------------------------------------------------------------------------
% Last updated 10/14/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% Option: display boxes with figure labels (a,b,c,etc.)
boxOn = false;

% Specify the locations of the seismic data (traces)
sacDir = '../data/traces/';
data{1} = 'PP.S0001.00.HHZ.D.2020.253.060000.SAC';
data{2} = 'PP.S0002.00.HHZ.D.2020.253.060000.SAC';
data{3} = 'PP.S0002.10.HNZ.D.2020.253.060000.SAC';
data{4} = 'AM.R36A4.00.EHZ.2020.09.09.06.00.00.SAC';

% Specify the locations of the response files
respDir = '../data/responses/';
resp{1} = 'PP.S0001.00.HHZ.resp';
resp{2} = 'PP.S0002.00.HHZ.resp';
resp{3} = 'PP.S0002.10.HNZ.resp';
resp{4} = 'SAC_PZs_R36A4_00';

% Read in the seismograms and responses
for i = 1:length(data)
    [s{i}.t,s{i}.d,s{i}.h] = fread_sac(fullfile(sacDir,data{i}));
    if i < 4
        [z{i},p{i},k{i}] = parseRESP(fullfile(respDir,resp{i}));
    else
        [z{i},p{i},k{i}] = parsePZ(fullfile(respDir,resp{i}));
    end
end

%--------------------------------------------------------------------------
% Instrument response removal
%--------------------------------------------------------------------------

% Transfer parameters
delta = 0.01;
flims = [0.1 0.2 10 20];

% Desired output units
units = 'velocity';

% Remove the instrument response from S0001
s1resp = fullfile(respDir,resp{1});
s{1}.d = s{1}.d - mean(s{1}.d);
s{1}.d = transfer(s{1}.d,delta,flims,units,s1resp,'resp');

% Remove the instrument response from S0002
s2resp = fullfile(respDir,resp{2});
s{2}.d = s{2}.d - mean(s{2}.d);
s{2}.d = transfer(s{2}.d,delta,flims,units,s2resp,'resp');

% Remove the instrument response from S0002 (Accelerometer - special case)
s3resp = fullfile(respDir,resp{3});
s{3}.d = s{3}.d - mean(s{3}.d);
if strcmp(units,'displacement')
    s{3}.d = transfer(s{3}.d,delta,flims,'velocity',s3resp,'resp');
elseif strcmp(units,'velocity')
    s{3}.d = transfer(s{3}.d,delta,flims,'acceleration',s3resp,'resp');
elseif strcmp(units,'acceleration')
    s{3}.d = transfer(s{3}.d,delta,flims,'test',s3resp,'resp');
end

% Remove the instrument response from the Raspberry Shake (special case)
s{4}.d = s{4}.d - mean(s{4}.d);
raspPZ = fullfile(respDir,resp{4});
if strcmp(units,'displacement')
    s{4}.d = transfer(s{4}.d,delta,flims,'raspdisp',raspPZ,'sacpz');
elseif strcmp(units,'velocity')
    s{4}.d = transfer(s{4}.d,delta,flims,'displacement',raspPZ,'sacpz');
elseif strcmp(units,'acceleration')
    s{4}.d = transfer(s{4}.d,delta,flims,'velocity',raspPZ,'sacpz');
end

%--------------------------------------------------------------------------
% FIGURE 3
%--------------------------------------------------------------------------

% Plot parameters
xmin = 10;
xmax = 60;
ymin = -7e4;
ymax = 7e4;

color{1} = [0.25 0.5 0.15];
color{2} = [0 0.9 1];
color{3} = [0 0 0];
color{4} = [0.89 0.043 0.365];

ylims = [-2e-4, -1e-4, 0, 1e-4, 2e-4];
[rdx,ex] = radexp(ylims);
ytl = cell(size(ylims));
for i = 1:length(ylims)
    if rdx(i) == 0
        ytl{i} = '0';
    else
        ytl{i} = sprintf('$%i{\\times}10^{%i}$',rdx(i),ex(i));
    end
end

% Raspberry Shake Data
subplot(4,1,1)
plot(s{4}.t,s{4}.d,'Color',color{4})
xlim([xmin xmax])
ax1 = gca;
ax1.TickDir = 'out';
grid on
ax1.XTick = 10:2.5:60;
ax1.XTickLabel = {''};
% Displacement
if strcmp(units,'displacement')
    ylim([-7e-6 7e-6])
    ax1.YTick = [-6e-6, -3e-6, 0, 3e-6, 6e-6];
    ax1.YTickLabel = {'-6e-6','-3e-6','0','3e-6','6e-6'};
    ylabel('Displacement (m)')
    text(3.5,-5.5e-6,'\textbf{R36A4.00.EHZ}','Rotation',90)
% Velocity
elseif strcmp(units,'velocity')
    ylim([-2e-4 2e-4])
    ax1.YTick = ylims;
    ax1.YTickLabel = ytl;
    ylabel('Velocity (m/s)')
    % text(3.5,-1.6e-4,'\textbf{R36A4.00.EHZ}','Rotation',90)
    text(1.5,-1.6e-4,'\textbf{R36A4.00.EHZ}','Rotation',90)
    ax1.Position(1) = 0.15;
% Acceleration
elseif strcmp(units,'acceleration')
    ylim([-0.02 0.02])
    ax1.YTick = [-0.02 -0.01 0 0.01 0.02];
    ax1.YTickLabel = {'-0.02','-0.01','0','0.01','0.02'};
    ylabel('Acceleration (m/s$^{2}$)')
    text(4,-0.016,'\textbf{R36A4.00.EHZ}','Rotation',90)
end
rectangle('Position',[18.5 -0.5e-4 1.5 1e-4])
if boxOn == true
    rectangle('Position',[10.6 -1.85e-4 2.15 0.75e-4],'FaceColor',[1 1 1])
    text(10.85,-1.5e-4,'(a)','FontSize',12)
    title('2020-09-09 06:00:13 $M_{Lg}$ = 3.1 Marlboro, New Jersey Earthquake')
end

ax1.Title.FontSize = 15.0;
ax1.Title.Position(2) = 2.25e-4;

% S0001 Data
subplot(4,1,2)
plot(s{1}.t,s{1}.d,'Color',color{1})
xlim([xmin xmax])
ax2 = gca;
ax2.TickDir = 'out';
grid on
ax2.XTick = 10:2.5:60;
ax2.XTickLabel = {''};
% Displacement
if strcmp(units,'displacement')
    ylim([-7e-6 7e-6])
    ax2.YTick = [-6e-6, -3e-6, 0, 3e-6, 6e-6];
    ax2.YTickLabel = {'-6e-6','-3e-6','0','3e-6','6e-6'};
    ylabel('Displacement (m)')
    text(3.5,-5.25e-6,'\textbf{S0001.00.HHZ}','Rotation',90)
% Velocity
elseif strcmp(units,'velocity')
    ylim([-2e-4 2e-4])
    ax2.YTick = ylims;
    ax2.YTickLabel = ytl;
    ylabel('Velocity (m/s)')
    ylabel('Velocity (m/s)')
    text(1.5,-1.6e-4,'\textbf{S0001.00.HHZ}','Rotation',90)
    ax2.Position(1) = 0.15;
% Acceleration
elseif strcmp(units,'acceleration')
    ylim([-0.02 0.02])
    ax2.YTick = [-0.02 -0.01 0 0.01 0.02];
    ax2.YTickLabel = {'-0.02','-0.01','0','0.01','0.02'};
    ylabel('Acceleration (m/s$^{2}$)')
    text(4,-0.015,'\textbf{S0001.00.HHZ}','Rotation',90)
end
rectangle('Position',[18.5 -0.5e-4 1.5 1e-4])
if boxOn == true
    rectangle('Position',[10.6 -1.85e-4 2.15 0.75e-4],'FaceColor',[1 1 1])
    text(10.85,-1.5e-4,'(b)','FontSize',12)
end

% S0002 Data
subplot(4,1,3)
plot(s{2}.t,s{2}.d,'Color',color{2})
xlim([xmin xmax])
ax3 = gca;
ax3.TickDir = 'out';
grid on
ax3.XTick = 10:2.5:60;
ax3.XTickLabel = {''};
% Displacement
if strcmp(units,'displacement')
    ylim([-7e-6 7e-6])
    ax3.YTick = [-6e-6, -3e-6, 0, 3e-6, 6e-6];
    ax3.YTickLabel = {'-6e-6','-3e-6','0','3e-6','6e-6'};
    ylabel('Displacement (m)')
    text(3.5,-5.25e-6,'\textbf{S0002.00.HHZ}','Rotation',90)
% Velocity
elseif strcmp(units,'velocity')
    ylim([-2e-4 2e-4])
    ax3.YTick = ylims;
    ax3.YTickLabel = ytl;
    ylabel('Velocity (m/s)')
    text(1.5,-1.6e-4,'\textbf{S0002.00.HHZ}','Rotation',90)
    ax3.Position(1) = 0.15;
% Acceleration
elseif strcmp(units,'acceleration')
    ylim([-0.02 0.02])
    ax3.YTick = [-0.02 -0.01 0 0.01 0.02];
    ax3.YTickLabel = {'-0.02','-0.01','0','0.01','0.02'};
    ylabel('Acceleration (m/s$^{2}$)')
    text(4,-0.015,'\textbf{S0002.00.HHZ}','Rotation',90)
end
rectangle('Position',[18.5 -0.5e-4 1.5 1e-4])
if boxOn == true
    rectangle('Position',[10.6 -1.85e-4 2.15 0.75e-4],'FaceColor',[1 1 1])
    text(10.85,-1.5e-4,'(c)','FontSize',12)
end

% S0002 Data (Accelerometer)
subplot(4,1,4)
plot(s{3}.t,s{3}.d,'Color',color{3})
xlim([xmin xmax])
ax4 = gca;
ax4.TickDir = 'out';
grid on
ax4.XTick = 10:2.5:60;
ax4.XTickLabel = {'','','06:00:15','','','','06:00:25','','','',...
                  '06:00:35','','','','06:00:45','','','','06:00:55',...
                  '',''};
% Displacement
if strcmp(units,'displacement')
    ylim([-7e-6 7e-6])
    ax4.YTick = [-6e-6, -3e-6, 0, 3e-6, 6e-6];
    ax4.YTickLabel = {'-6e-6','-3e-6','0','3e-6','6e-6'};
    ylabel('Displacement (m)')
    text(3.5,-5.25e-6,'\textbf{S0002.10.HNZ}','Rotation',90)
% Velocity
elseif strcmp(units,'velocity')
    ylim([-2e-4 2e-4])
    ax4.YTick = ylims;
    ax4.YTickLabel = ytl; 
    ylabel('Velocity (m/s)')
    text(1.5,-1.6e-4,'\textbf{S0002.10.HNZ}','Rotation',90)
    ax4.Position(1) = 0.15;
% Acceleration
elseif strcmp(units,'acceleration')
    ylim([-0.02 0.02])
    ax4.YTick = [-0.02 -0.01 0 0.01 0.02];
    ax4.YTickLabel = {'-0.02','-0.01','0','0.01','0.02'};
    ylabel('Acceleration (m/s$^{2}$)')
    text(4,-0.015,'\textbf{S0002.10.HNZ}','Rotation',90)
end
rectangle('Position',[18.5 -0.5e-4 1.5 1e-4])
if boxOn == true
    rectangle('Position',[10.6 -1.85e-4 2.15 0.75e-4],'FaceColor',[1 1 1])
    text(10.85,-1.5e-4,'(d)','FontSize',12)
    xlabel('Time (UTC)')
end

% Make an inset showing a zoomed in section
axes('Position',[0.735 0.773 0.21 0.156])
box on
plot(s{4}.t,s{4}.d,'Color',color{4})
hold on
scatter(19.011,1.039e-5,20,'Filled','MarkerFaceColor',color{4})
text(18.55,4e-5,'6.011 s, $1.039{\times}10^{-5}$ m/s','FontSize',11)
grid on
xlim([18.5 20])
ylim([-0.5e-4 0.5e-4])
ax5 = gca;
ax5.XTick = [];
ax5.YTick = [];
ax5.XTickLabel = {};
ax5.YTickLabel = {};

axes('Position',[0.735 0.553 0.21 0.156])
box on
plot(s{1}.t,s{1}.d,'Color',color{1})
hold on
scatter(18.98,8.159e-6,20,'Filled','MarkerFaceColor',color{1})
text(18.55,4e-5,'5.98 s, $8.159{\times}10^{-6}$ m/s','FontSize',11)
grid on
xlim([18.5 20])
ylim([-0.5e-4 0.5e-4])
ax6 = gca;
ax6.XTick = [];
ax6.YTick = [];
ax6.XTickLabel = {};
ax6.YTickLabel = {};

axes('Position',[0.735 0.334 0.21 0.156])
box on
plot(s{2}.t,s{2}.d,'Color',color{2})
hold on
scatter(18.98,8.236e-6,20,'Filled','MarkerFaceColor',color{2})
text(18.55,4e-5,'5.98 s, $8.236{\times}10^{-6}$ m/s','FontSize',11)
grid on
xlim([18.5 20])
ylim([-0.5e-4 0.5e-4])
ax6 = gca;
ax6.XTick = [];
ax6.YTick = [];
ax6.XTickLabel = {};
ax6.YTickLabel = {};

axes('Position',[0.735 0.114 0.21 0.156])
box on
plot(s{3}.t,s{3}.d,'Color',color{3})
hold on
scatter(18.98,7.398e-6,20,'Filled','MarkerFaceColor',color{3})
text(18.55,4e-5,'5.98 s, $7.398{\times}10^{-6}$ m/s','FontSize',11)
grid on
xlim([18.5 20])
ylim([-0.5e-4 0.5e-4])
ax6 = gca;
ax6.XTick = [];
ax6.YTick = [];
ax6.XTickLabel = {};
ax6.YTickLabel = {};

% Print the figure to a PDF and open it with Preview
set(gcf,'Position',[0 0 600 600])
print(gcf,'-dpdf','-r600','Figure3');
close;

system('open Figure3.pdf &');

