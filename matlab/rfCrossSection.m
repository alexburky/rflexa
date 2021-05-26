% rfCrossSection.m
%
% This script produces cross sectional plots through a volume of
% receiver function data for the Eastern US.
%
% Beta: Adding an option to overlay the cross section on a tomographic
%       cross section?
%
% Inputs: Start and end coordinates of cross section?
%
%------------------------------------------------------------------
% Last updated 5/26/2021 by aburky@princeton.edu
%------------------------------------------------------------------

clear,clc

% Location of Depth Converted Stacks
stackDir = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/',...
            'Stacks/LLNL_3D/FILTERED_0.02_0.2/GW10/'];
stacks = dir(fullfile(stackDir,'M*.mat'));

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
    % Option 1: 1D corrected data
    % rf{i}.d = stack;
    % rf{i}.evla = evla;
    % rf{i}.evlo = evlo;
    
    % Option 2: 3D corrected data
    rf{i}.d = mbstk;
    rf{i}.evla = 30;
    rf{i}.evlo = 30;
    z = x;
    
end

%% Create a volume from all of the stacks for slicing

% 696 * 696 * 8001 = ~4e9
% This is 4 billion floats, a float is 4 bytes
% a double is 8 bytes
% This is 16GB or 32GB
% There must be a better way?

%% Plot a cross section!

% Location to save the resulting cross sectional image
imgDir = ['/Users/aburky/IFILES/NETWORKS/TA_Analysis/',...
          'XSections/FILTERED_0.2_5/GW20/'];

% Disbale figure plotting
set(0,'DefaultFigureVisible','off');

% Crust or Mantle Section?
% section = 'crust';
section = 'mantle';

% List of stations to make cross sections through (const. latitude)
stats = {'1','2','3','4','5','C','D','E','F','G','H','I','J','K','L',...
         'M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

% Crust scale factor
if strcmp(section,'crust')
    scale = 2;
end

% MTZ scale factor
if strcmp(section,'mantle')
    scale = 10;
end

for stat = 1:length(stats)
    clf
    p1 = subplot(1,2,1);
    minLon = 180;
    maxLon = -180;
    for i = 1:length(rf)
        if strncmpi(rf{i}.code,stats{stat},1) && ...
           strncmpi('KMSC',rf{i}.code,4) == 0 && ...
           strncmpi('SPMN',rf{i}.code,4) == 0 && ...
           strncmpi('SFIN',rf{i}.code,4) == 0 && ...
           strncmpi('TIGA',rf{i}.code,4) == 0
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

            % Focus on the crust
            if strcmp(section,'crust')
                ylim([0 100])
            end

            % Focus on the transition zone
            if strcmp(section,'mantle')
                ylim([300 750])
            end

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
        title(sprintf('TA Cross Section: Station Code %s*',stats{stat}))
        ax1.FontSize = 12;
        plot([-200 200],[410 410],'k--')
        plot([-200 200],[660 660],'k--')
    end
    xlim([minLon - 1 maxLon + 1])
    set(p1,'Position',[0.1 0.15 0.4 0.75]);

    % Plot a map above the cross section to give geographic context
    axes('Position',[0.5 0.525 0.55 0.45]);
    landareas = shaperead('landareas.shp','UseGeoCoords',true);
    axesm('mercator','MapLatLimit',[20 50],'MapLonLimit',[-100 -60]);
    geoshow(landareas,'FaceColor',[0.75 0.75 0.75],'EdgeColor','k',...
        'LineWidth',0.1);
    axis off
    framem
    gridm
    gax1 = gca;
    setm(gca,'parallellabel','on')
    setm(gca,'meridianlabel','on')
    setm(gca,'MLineLocation',[-95 -90 -85 -80 -75 -70 -65])
    setm(gca,'MLabelLocation',[-95, -85, -75, -65])
    setm(gca,'PLineLocation',[25 30 35 40 45])
    setm(gca,'PLabelLocation',[25, 30, 35, 40, 45])
    setm(gax1,'FontName','CMU Serif');

    j = 0;
    stla_m = 0;
    stlo_m = 0;
    for i = 1:length(rf)
        if strncmpi(rf{i}.code,stats{stat},1) && ...
           strncmpi('KMSC',rf{i}.code,4) == 0 && ...
           strncmpi('SPMN',rf{i}.code,4) == 0 && ...
           strncmpi('SFIN',rf{i}.code,4) == 0 && ...
           strncmpi('TIGA',rf{i}.code,4) == 0
            scatterm(rf{i}.lat,rf{i}.lon,40,'r','^','filled',...
                'MarkerEdgeColor','k');
            stla_m = stla_m + rf{i}.lat;
            stlo_m = stlo_m + rf{i}.lon;
            j = j + 1;
        end
    end
    stla_m = stla_m/j;
    stlo_m = stlo_m/j;

    % Plot a map at the bottom right showing the distribution of events
    axes('Position',[0.55 0.075 0.45 0.45]);
    axesm('eqdazim','Origin',[stla_m,stlo_m]);
    geoshow(landareas,'FaceColor',[0.75 0.75 0.75]);
    axis off
    framem
    %gridm
    for i = 1:length(rf)
        if strncmpi(rf{i}.code,stats{stat},1) && ...
           strncmpi('KMSC',rf{i}.code,4) == 0 && ...
           strncmpi('SPMN',rf{i}.code,4) == 0 && ...
           strncmpi('SFIN',rf{i}.code,4) == 0 && ...
           strncmpi('TIGA',rf{i}.code,4) == 0
            scatterm(rf{i}.evla,rf{i}.evlo,20,'r','filled',...
                'MarkerEdgeColor','k');
        end
    end

    % sgtitle('Filtered 0.02 - 0.2 Hz')

    % Save the completed figure!
    if strcmp(section,'crust')
        fName = ['TA_',stats{stat},'*_CRUST_XS.png'];
        print(fullfile(imgDir,fName),'-dpng','-r600');
    end

    if strcmp(section,'mantle')
        fName = ['TA_',stats{stat},'*_MANTLE_XS.png'];
        print(fullfile(imgDir,fName),'-dpng','-r600');
    end

end

%% Load in the SEMUCB Tomography Model

% Load model SEMUCB_wm1 (do this once, it takes a long time!)
load('/Users/aburky/IFILES/MODELS/SEMUCB_WM1/SEMUCBwm1_interpolant.mat');

% Load the colormap for plotting!
load('/Users/aburky/MFILES/colormaps/vik.mat');

%% (Beta Zone) Overlay the cross section on a tomographic plot

% Get the minimum and maximum longitude of current batch of RFs
minLon = rf{1}.lon;
minLat = rf{1}.lat;
maxLon = rf{1}.lon;
maxLat = rf{1}.lat;
for i = 1:length(rf)
    if rf{i}.lon < minLon
        minLon = rf{i}.lon;
        minLat = rf{i}.lat;
    end
    if rf{i}.lon > maxLon
        maxLon = rf{i}.lon;
        maxLat = rf{i}.lat;
    end
end

%% Part 1a: Prepare tomography data for plotting

degrees = [minLat; minLon; maxLat; maxLon];

% Desired number of depth points (resolution)
n = 500;

% Depth range to plot
minDep = 300;
maxDep = 750;

% Desired depths to plot
dep = linspace(minDep, maxDep, n);

% Generate lat-lon grid
[lat,lon] = gcwaypts(degrees(1),degrees(2),...
                     degrees(3),degrees(4),n-1);
                    
% Create a grid of model values for each (lat,lon,depth)
absVs = zeros(n,n);
for j = 1:n
    [x,y,z_t] = sph2cart((pi/180)*lon(j),...
                        (pi/180)*lat(j),6371-dep);
    absVs(j,:) = Vs(x,y,z_t);
end

% Pre-allocate variables
meanVs = zeros(1, n);
delVs = zeros(n, n);

% Construct array of averages for each depth slice
for j = 1:n
    meanVs(j) = mean(absVs(:,j));
    % Replace absolute velocities with relative velocities
    delVs(:,j) = ((absVs(:,j)-meanVs(j))./meanVs(j))*100;
end

% Rotate the velocity grid to give it the right sense
flipVs = transpose(delVs);
rotVs = rot90(flipVs,2);

%% Part 1b: Plot the tomographic cross section

% Generate x and y coordinates consistent with RF data
y = linspace(minDep,maxDep,n);
x = linspace(minLon-1,maxLon+1,n);

% pcolor(x,y,rotVs);
% shading interp
% hold on
% [C,H] = contour(x,y,rotVs,0,'k');
[C,H] = contourf(x,y,fliplr(flipud(rotVs)),50,'EdgeColor','none');
ax1 = gca;
ax1.YDir = 'reverse';
set(H,'linewidth',0.01)
caxis([-2 2])
colormap(flipud(vik))
h = colorbar('location','southoutside');
set(h,'TickLabelInterpreter','latex')
% set(h,'position',[0.4 0.05 0.25 0.05])
ylabel(h,'$\delta V_{S}$ from Section 1-D Average (\%)')
h.Label.Interpreter = 'latex';
hold on

%% Part 2: Overlay the receiver functions
scale = 10;
for i = 1:length(rf)
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
end
ylim([300 750])
xlim([minLon-1 maxLon+1])
ax1 = gca;
ax1.YDir = 'reverse';
ax1.TickDir = 'out';
xlabel('Longitude ($^{\circ}$)')
ylabel('Depth (km)')
% title(sprintf('TA Cross Section: Station Code %s*',stats{stat}))
ax1.FontSize = 12;
plot([-200 200],[410 410],'k--')
plot([-200 200],[660 660],'k--')
