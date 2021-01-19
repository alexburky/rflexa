% rfCrossSection.m
%
% This script produces cross sectional plots through a volume of
% receiver function data for the Eastern US.
%
% Inputs: Start and end coordinates of cross section?
%
%------------------------------------------------------------------
% Last updated 1/19/2021 by aburky@princeton.edu
%------------------------------------------------------------------

clear,clc

% Location of Depth Converted Stacks
stackDir = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/',...
            'Stacks_Filt_0.02_0.2_GW10/'];
stacks = dir(fullfile(stackDir,'*.mat'));

% Location of station latitude abd longitude information
metaDir = '/Users/aburky/GFILES/TA_Website';
metaFile = 'TA_Locations.dat';

metaData = readtable(fullfile(metaDir,metaFile));

for i = 1:length(stacks)
    load(fullfile(stackDir,stacks(i).name));
    station = extractBefore(stacks(i).name,'_');
    for j = 1:height(metaData)
        if strcmp(metaData(j,1).Var1{1},station) == 1
            rf{i}.code = station;
            rf{i}.lat = metaData(j,3).Var3;
            rf{i}.lon = metaData(j,2).Var2;
        end
    end
    rf{i}.d = stack;
end

%% Plot a cross section!

% Location to save the resulting cross sectional image
imgDir = '/Users/aburky/IFILES/NETWORKS/TA_Analysis/TA_XSections/';

% Do this for stations D through Z

% Look for all stations starting with station code 'T'
% (This plots a line of station with roughly equal latitude)
stat = 'C';
p1 = subplot(1,2,1);
scale = 10;
minLon = 180;
maxLon = -180;
for i = 1:length(rf)
    if strncmpi(rf{i}.code,stat,1) && strncmpi('KMSC',rf{i}.code,4) == 0 && strncmpi('SPMN',rf{i}.code,4) == 0 && strncmpi('SFIN',rf{i}.code,4) == 0 && strncmpi('TIGA',rf{i}.code,4) == 0
        pos = scale*rf{i}.d + rf{i}.lon;
        pos(pos < rf{i}.lon) = rf{i}.lon;
        xPatchPos = [rf{i}.lon * ones(size(pos)), fliplr(pos)];
        yPatchPos = [z, fliplr(z)];
        fill(xPatchPos,yPatchPos,'r','LineStyle','None')
        hold on
        neg = scale*rf{i}.d + rf{i}.lon;
        neg(neg > rf{i}.lon) = rf{i}.lon;
        xPatchNeg = [rf{i}.lon * ones(size(neg)), fliplr(neg)];
        yPatchNeg = [z, fliplr(z)];
        fill(xPatchNeg,yPatchNeg,'b','LineStyle','None')
        plot(scale*rf{i}.d + rf{i}.lon,z,'k')
        ylim([300 750])
        
        if rf{i}.lon < minLon
            minLon = rf{i}.lon;
        end
        
        if rf{i}.lon > maxLon
            maxLon = rf{i}.lon;
        end
    end
    ax1 = gca;
    ax1.YDir = 'reverse';
    ax1.TickDir = 'out';
    xlabel('Longitude ($^{\circ}$)')
    ylabel('Depth (km)')
    title(sprintf('TA Cross Section: Station Code %s*',stat))
    ax1.FontSize = 12;
    plot([-200 200],[410 410],'k--')
    plot([-200 200],[660 660],'k--')
end
xlim([minLon - 1 maxLon + 1])
set(p1,'Position',[0.1 0.15 0.4 0.75]);

% Plot a map above the cross section to give geographic context
axes('Position',[0.55 0.5 0.45 0.45]);
landareas = shaperead('landareas.shp','UseGeoCoords',true);
axesm('mercator','MapLatLimit',[20 50],'MapLonLimit',[-110 -60]);
geoshow(landareas,'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.1);
axis off
framem
gridm
setm(gca,'parallellabel','on')
setm(gca,'meridianlabel','on')

for i = 1:length(rf)
    if strncmpi(rf{i}.code,stat,1) && strncmpi('KMSC',rf{i}.code,4) == 0 && strncmpi('SPMN',rf{i}.code,4) == 0 && strncmpi('SFIN',rf{i}.code,4) == 0 && strncmpi('TIGA',rf{i}.code,4) == 0
        scatterm(rf{i}.lat,rf{i}.lon,40,'r','^','filled','MarkerEdgeColor','k');
    end
end

% Save the completed figure!
fName = ['TA_',stat,'*_XS.png'];
print(fullfile(imgDir,fName),'-dpng','-r600');

