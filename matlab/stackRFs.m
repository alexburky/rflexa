% stackRFs.m
%
% This script takes a directory containing receiver functions and depth
% converts and stacks them. The stacks are bootstrap resampled to calculate
% standard deviations. The output of this script is a plot of the depth
% converted receiver function stack
%
%--------------------------------------------------------------------------
% Last updated 1/20/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

station = 'N41A';

% Directory containing the receiver function data
rfdir = ['/Users/aburky/IFILES/NETWORKS/TA/',station,...
         '/NULL/RFUNCS_VEL/FILTERED_0.02_0.2/GW10/'];
rfs = dir(fullfile(rfdir,'*RF.SAC'));

% Directory to save the resulting image to
imgDir = '/Users/aburky/IFILES/NETWORKS/TA_Analysis/TA_Plots/';

% Flag indicating whether PP is included or not
PP = false;

% Set QC parameters
snr_z = 2;
snr_r = 2;
nu = 0.1;
fit = 80;
if PP == true
    gcarc = 29;
else
    gcarc = 35;
end

% Read in receiver functions which match QC criteria
j = 1;
for i = 1:length(rfs)
    [rf_0{i}.t, rf_0{i}.d, rf_0{i}.h] = fread_sac(fullfile(rfdir,...
                                            rfs(i).name));
    if rf_0{i}.h.user(1) > snr_z && rf_0{i}.h.user(2) > snr_r && ...
       rf_0{i}.h.user(3) > fit && rf_0{i}.h.user(4) > nu && ...
       rf_0{i}.h.gcarc >= gcarc
        % Keep receiver functions which pass all QC
        rf{j} = rf_0{i};
        gcarcs(j) = rf{j}.h.gcarc;
        evla(j) = rf{j}.h.evla;
        evlo(j) = rf{j}.h.evlo;
        rf{j}.d = rf{j}.d/max(rf{j}.d(1:round(20/rf{j}.h.delta)));
        [~,shidx] = max(rf{j}.d(1:round(20/rf{j}.h.delta)));
        rf{j}.d = rf{j}.d(shidx:end);
        npts(j) = length(rf{j}.d);
        P(j) = rf{j}.h.user(10);
        j = j + 1;
    end
end

% clear rf_0 rfs rfdir

% Make the length of the data vectors and sample rates consistent
npmax = max(unique(npts));
for i = 1:length(rf)
    if length(rf{i}.d) == npmax
        dt = rf{i}.h.delta;
        trf = 0:dt:(npmax*dt);
        trf = trf(1:npmax);
        break
    end
end

% Upsample all data to the maximum sample rate of the dataset
for i = 1:length(rf)
    if rf{i}.h.npts ~= npmax
        rf{i}.t = 0:rf{i}.h.delta:length(rf{i}.d);
        rf{i}.t = rf{i}.t(1:length(rf{i}.d));
        rf{i}.d = interp1(rf{i}.t,rf{i}.d,trf);
        rf{i}.t = trf;
    end
end

% Depth convert and stack receiver functions!

% Depth interval
dz = 0.1;
z = 0:dz:800;

stack = zeros(size(z));

for i = 1:length(rf)
    rf{i}.depth = rfDepcon(rf{i}.d,dt,P(i)/deg2km(1),dz,'iasp91','false');
    stack = stack + rf{i}.depth;
end

stack = stack/length(rf);

% Bootstrap to determine standard deviations
bn = 1000;
bstk = zeros(bn,length(z));

% Resample
for b = 1:bn
    zran = randi(length(rf),[length(rf),1]);
    for i = 1:length(zran)
        bstk(b,:) = bstk(b,:) + rf{zran(i)}.depth;
    end
end

bstk = bstk/length(P);

% Calculate standard deviations
bdot = sum(bstk,1)/bn;
num = zeros(1,length(z));
for b = 1:bn
    num = num + sum((bstk(b,:) - bdot).^2, 1);
end
sdv = (num/(bn-1)).^(1/2);
mbstk = mean(bstk);
pos = mbstk - (2.*sdv);
pos(pos < 0) = 0;
neg = mbstk + (2.*sdv);
neg(neg > 0) = 0;

%% Make a plot of the stack and additional statistics

for i = 1:length(rf_0)
    % Save statistics as vectors
    mag(i,1) = rf_0{i}.h.mag;
    gcarc(i,1) = rf_0{i}.h.gcarc;
    nu(i,1) = rf_0{i}.h.user(4);
    fit(i,1) = rf_0{i}.h.user(3);
end

nbins = [20 20];

H1 = [fit, nu];
H2 = [gcarc, mag];

mygray = flipud(gray);
mygray(2:end,:) = mygray(2:end,:)*0.95;

% Initialize the figure
% f = figure('visible','off');
% set(0,'DefaultFigureVisible','off');

% GCARC vs. Magnitude
p1 = subplot(2,2,1);
hist3(H2,nbins,'CdataMode','auto','EdgeColor','none')
view(2)
axis square
ax1 = gca;
grid off
ax1.Box = 'on';
ax1.TickDir = 'out';
ax1.TickLength = [0.025; 0.025];
ax1.FontSize = 12;
xlabel('$\Delta (^{\circ})$')
ylabel('$M_{w}$')
ylim([5.35 8.25])
title(sprintf('TA Station: %s\n$N = %i, N'' = %i$',station,...
        length(rf_0), length(rf)))
ax1.Title.Position(1) = 115;
ax1.Title.Position(2) = 8.35;
ax1.Position(2) = 0.565;

% Fit vs. Nu
p2 = subplot(2,2,2);
hist3(H1,nbins,'CdataMode','auto','EdgeColor','none')
view(2)
axis square
ax2 = gca;
grid off
ax2.Box = 'on';
ax2.TickDir = 'out';
ax2.TickLength = [0.025; 0.025];
xlabel('$(1-\chi)\times100$')
ylabel('$\nu$')
hold on
plot3([-10 110],[0.1 0.1],[500 500],'r')
plot3([80 80],[-1 1],[500 500],'r')
xlim([-1 101])
ylim([-0.7 0.7])
ax2.YTick = [-0.5 -0.25 0 0.25 0.5];
ax2.XTick = [0 25 50 75 100];
caxis([0 10])
colormap(mygray)
ax2.FontSize = 12;
ax2.Position(2) = 0.565;

% Bootstrapped stack
p3 = subplot(2,2,[3,4]);
xPatchPos = [z, fliplr(z)];
yPatchPos = [zeros(size(mbstk)), fliplr(pos)];
fill(xPatchPos,yPatchPos,'r','LineStyle','None');
hold on
xPatchNeg = [z, fliplr(z)];
yPatchNeg = [zeros(size(mbstk)), fliplr(neg)];
fill(xPatchNeg,yPatchNeg,'b','LineStyle','None');
plot(z,stack,'k')
ax3 = gca;
xlim([300 750])
ylim([-0.1 0.1])
plot(z,mbstk + 2.*sdv,'k')
plot(z,mbstk - 2.*sdv,'k')
grid on
ylabel('$\bar{f}_{Z \rightarrow R}(p,z)$')
xlabel('Depth (km)')
title('Pre-Filtered 5 - 50s')
ax3.FontSize = 12;
ax3.TickDir = 'out';
ax3.YTick = [-0.1 -0.05 0 0.05 0.1];

% Save the completed figure!
fName = [station,'_Summary.png'];
print(fullfile(imgDir,fName),'-dpng','-r300');

% Save the depth converted stack to a .mat file
fName = [station,'_Stack.mat'];
save(fullfile(rfdir,fName),'z','stack','evla','evlo');
