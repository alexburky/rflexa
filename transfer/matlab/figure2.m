% figure2.m
%
% This code is meant to serve as a companion to the 2021 EduQuakes paper,
% 'Instrument response removal and the 2020 M3.1 Marlboro, New Jersey,
% earthquake.' It produces the second figure of the paper.
%
% External dependencies:
% - fread_sac
% - parseRESP
% - parsePZ
%
%--------------------------------------------------------------------------
% Last updated 10/13/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% Option: display boxes with figure labels (a,b,c,etc.)
boxOn = false;

% Option: display inset map
inset = false;

% Option: ONLY display inset map
insetOnly = false;

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
% FIGURE 2
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

% Y Values of labels
ylims = [-6e4 -3e4 0 3e4 6e4];
% Construct the Y Axis Label strings
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
if insetOnly == false
subplot(4,1,1)
plot(s{4}.t,s{4}.d - mean(s{4}.d),'Color',color{4})
xlim([xmin xmax])
ylim([ymin ymax])
ax1 = gca;
ax1.TickDir = 'out';
grid on
ax1.YTick = ylims;
ax1.YTickLabel = ytl;
ax1.XTick = 10:2.5:60;
ax1.XTickLabel = {''};
ylabel('Counts')
text(2.5,-5.75e4,'\textbf{R36A4.00.EHZ}','Rotation',90)
if boxOn == true
    rectangle('Position',[10.6 -6.5e4 2.15 2.75e4],'FaceColor',[1 1 1])
    text(10.85,-5.25e4,'(a)','FontSize',12)
    title(['2020-09-09 06:00:13 $M_{Lg}$ = 3.1 Marlboro ',...
        'New Jersey Earthquake'])
end
ax1.Position(1) = 0.14;
ax1.Title.FontSize = 15.0;
ax1.Title.Position(2) = 77500;

% Station S0001 Data
subplot(4,1,2)
plot(s{1}.t,s{1}.d - mean(s{1}.d),'Color',color{1})
xlim([xmin xmax])
ylim([ymin ymax])
ax2 = gca;
ax2.TickDir = 'out';
grid on
ax2.YTick = ylims;
ax2.YTickLabel = ytl;
ax2.XTick = 10:2.5:60;
ax2.XTickLabel = {''};
ylabel('Counts')
text(2.5,-5.5e4,'\textbf{S0001.00.HHZ}','Rotation',90)
if boxOn == true
    rectangle('Position',[10.6 -6.5e4 2.15 2.75e4],'FaceColor',[1 1 1])
    text(10.85,-5.25e4,'(b)','FontSize',12)
else
    ax2.XTick = 10:2.5:60;
    ax2.XTickLabel = {'','','06:00:15','','','','06:00:25','','','',...
                  '06:00:35','','','','06:00:45','','','','06:00:55',...
                  '',''};
end
ax2.Position(1) = 0.14;

% Station S0002 Data
subplot(4,1,3)
plot(s{2}.t,s{2}.d - mean(s{2}.d),'Color',color{2})
xlim([xmin xmax])
ylim([ymin ymax])
ax3 = gca;
ax3.TickDir = 'out';
grid on
ax3.YTick = ylims;
ax3.YTickLabel = ytl;
ax3.XTick = 10:2.5:60;
ax3.XTickLabel = {''};
ylabel('Counts')
text(2.5,-5.5e4,'\textbf{S0002.00.HHZ}','Rotation',90)
if boxOn == true
    rectangle('Position',[10.6 -6.5e4 2.15 2.75e4],'FaceColor',[1 1 1])
    text(10.85,-5.25e4,'(c)','FontSize',12)
end
ax3.Position(1) = 0.14;

% Accelerometer Data
subplot(4,1,4)
plot(s{3}.t,s{3}.d - mean(s{3}.d),'Color',color{3})
xlim([xmin xmax])
ylim([ymin ymax])
ax4 = gca;
ax4.TickDir = 'out';
grid on
ax4.YTick = ylims;
ax4.YTickLabel = ytl;
ax4.XTick = 10:2.5:60;
ax4.XTickLabel = {'','','06:00:15','','','','06:00:25','','','',...
                  '06:00:35','','','','06:00:45','','','','06:00:55',...
                  '',''};
ylabel('Counts')
text(2.5,-5.5e4,'\textbf{S0002.10.HNZ}','Rotation',90)
if boxOn == true
    rectangle('Position',[10.6 -6.5e4 2.15 2.75e4],'FaceColor',[1 1 1])
    text(10.85,-5.25e4,'(d)','FontSize',12)
    xlabel('Time (UTC)')
end
ax4.Position(1) = 0.14;
end

% Inset Map
if inset == true
    axes('Position',[0.535 0.11 0.38 0.3765])
    box on
    ax6 = gca;
    ax6.XLim = [0 1];
    ax6.YLim = [0 1];
    ax6.XTick = [];
    ax6.YTick = [];
    if boxOn == true
        rectangle('Position',[0.025 0.025 0.1 0.1],'FaceColor',[1 1 1])
        text(0.04,0.07,'(e)','FontSize',12)
    end

    axes('Position',[0.605 0.1865 0.3 0.275])
    geoshow('usastatehi.shp','FaceColor',[0.15 0.5 0.15])
    hold on
    geoshow(40.3460,-74.655,'DisplayType','Point','Marker','^',...
        'MarkerFaceColor',color{2},'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',8)
    geoshow(40.3325,-74.6572,'DisplayType','Point','Marker','^',...
        'MarkerFaceColor',color{4},'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',8)
    geoshow(40.302,-74.289,'DisplayType','Point','Marker','p',...
        'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',12)
    box on
    ax5 = gca;
    ax5.TickDir = 'out';
    ax5.TickLength = [0.025 0.025];
    ax5.XLabel.String = 'Longitude ($^{\circ}$)';
    ax5.YLabel.String = 'Latitude ($^{\circ}$)';
    ax5.XLim = [-75.5 -73.5];
    ax5.YLim = [39.5 41.5];
end

% Print the figure to a PDF and open it with Preview
set(gcf,'Position',[0 0 600 600])
if boxOn == true
    print(gcf,'-dpdf','-r600','Figure2');
else
    print(gcf,'-dpdf','-r600','Figure2_NoBox');
end
close;

% system('open Figure2.pdf &');
