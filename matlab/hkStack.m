% hkStack.m
%
% This script performs an H-k stack on a set of receiver function data,
% saved in a particular directory. This script has multiple outputs:
% - an image of the H-k stack (.png file)
% - the x,y,z data of the H-k stack (.xyz file)
% - a .mat file with the crustal thickness, Vp/Vs ratio, and station info
%
%--------------------------------------------------------------------------
% Last updated 3/11/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

% To do: Add an option to specify whether or not to Phase Weight Stack...

clear,clc

% Station code
stat = '339A';

% Location of receiver function data
rfDir = ['/Users/aburky/IFILES/NETWORKS/TA/',stat,'/NULL/RFUNCS_VEL'];
rfData = dir(fullfile(rfDir,'*SAC'));

% Location you would like to save the output data to
outDir = '/Users/aburky/IFILES/NETWORKS/TA_Analysis/hkStacks/';

% Quality control parameters
fit = 80;
nu = 0.05;

% Phase weight stack?
pwstk = 'true';

% Save the output?
saveOutput = 'false';

%--------------------------------------------------------------------------
% Read in the receiver function data
%--------------------------------------------------------------------------

j = 1;
for i = 1:length(rfData)
    [rf{i}.t,rf{i}.d,rf{i}.h] = fread_sac(fullfile(rfDir,rfData(i).name));
    % Save receiver functions which pass QC criteria
    if rf{i}.h.user(3) >= fit && rf{i}.h.user(4) >= nu
        grf{j} = rf{i};
        % Remove ten second time shift...
        shidx = round(10/grf{j}.h.delta);
        grf{j}.d = grf{j}.d(shidx:end);
        grf{j}.t = grf{j}.t(shidx:end) - 10;
        j = j + 1;
    end
end

clear rf;
rf = grf;
clear grf;

%--------------------------------------------------------------------------
% Compute the H-k stack
%--------------------------------------------------------------------------

% Change this: use the Crust1.0 Vp value nearest the station instead of
% the searching through possible Vp values -> this is how EARS does things

% H-k stack parameters
vp = 2.5:0.05:7.0;
h = 20.0:0.5:70;
k = linspace(1.3,2.3,length(h));
[H,K] = meshgrid(h,k);

PsW = 0.33333334;
PpPsW = 0.33333334;
PpSsW = 0.33333328;

for i = 1:length(vp)
    for j = 1:length(rf)
        % Ray-parameter in (s/km)
        % P(j) = rf{j}.h.user(10)/deg2km(1);
        P(j) = rf{j}.h.user(10)/111.1949;

        % Theoretical times of Pds phases
        tPs{j} = H.*(sqrt(K.^2./vp(i)^2-P(j)^2) - sqrt(vp(i)^(-2)-P(j)^2));
        tPpPs{j} = H.*(sqrt(K.^2./vp(i)^2-P(j)^2) + ...
                   sqrt(vp(i)^(-2)-P(j)^2));
        tPpSs{j} = 2*H.*(sqrt(K.^2./vp(i)^2-P(j)^2));

        Ps_idx = round(tPs{j}/rf{j}.h.delta);
        PpPs_idx = round(tPpPs{j}/rf{j}.h.delta);
        PpSs_idx = round(tPpSs{j}/rf{j}.h.delta);

        S{j} = PsW*rf{j}.d(round(tPs{j}/rf{j}.h.delta)) + ...
               PpPsW*rf{j}.d(round(tPpPs{j}/rf{j}.h.delta)) - ...
               PpSsW*rf{j}.d(round(tPpSs{j}/rf{j}.h.delta));
           
        % Optional: Phase Weight Stack
        if strcmp(pwstk,'true')
            Ps_pw = rf{j}.d(round(tPs{j}/rf{j}.h.delta)) + ...
                1i*hilbert(rf{j}.d(round(tPs{j}/rf{j}.h.delta)));
            PpPs_pw = rf{j}.d(round(tPpPs{j}/rf{j}.h.delta)) + ...
                1i*hilbert(rf{j}.d(round(tPpPs{j}/rf{j}.h.delta)));
            PpSs_pw = rf{j}.d(round(tPpSs{j}/rf{j}.h.delta)) + ...
                1i*hilbert(rf{j}.d(round(tPpSs{j}/rf{j}.h.delta)));
             
            S_pw{j} = PsW*exp(1i*angle(Ps_pw)) + ...
                PpPsW*exp(1i*angle(PpPs_pw)) - ...
                PpSsW*exp(1i*angle(PpSs_pw));
        end
    end
    
    % Optional: Phase Weight Stack
    if strcmp(pwstk,'true')
        sS_pw{i} = 0;
        for j = 1:length(rf)
            sS_pw{i} = sS_pw{i} + S_pw{j};
        end
        sS_pw{i} = sS_pw{i}/length(rf);
    end
    
    % Stack!
    sS{i} = 0;
    for j = 1:length(rf)
        % Optional: Phase Weight Stack
        if strcmp(pwstk,'true')
            sS{i} = sS{i} + S{j}.*abs(sS_pw{i});
        else
            sS{i} = sS{i} + S{j};
        end
    end
    sS{i} = sS{i}/length(rf);
    
    % Get maximum value of H-k stack for current Vp value
    [~,idx] = max(sS{i}(:));
    maxS(i) = sS{i}(idx);
end

% Get the maximum value of S(H,k) for all iterations
[~,idx] = max(maxS);

[x,y] = find(sS{idx} == max(max(sS{idx})));

% Initialize the figure
if strcmp(saveOutput,'true')
    f = figure('visible','off');
    set(0,'DefaultFigureVisible','off');
end

% Plot the maximum values of the H-k stack for each Vp value
s1 = subplot(1,2,1);
scatter(maxS,vp,30,[0.7 0.7 0.7],'Filled')
ylabel('$V_{P}$ (km/s)')
xlabel('S(H,$\kappa$)')
hold on
box on
grid on
scatter(maxS(idx),vp(idx),60,'k','Filled','MarkerEdgeColor','k')
ylim([min(vp) - 0.01*max(vp) max(vp) + 0.01*max(vp)])
xlim([min(maxS) - 0.01*max(maxS) max(maxS) + 0.01*max(maxS)])
s1.Position = [0.1 0.1 0.3 0.8];
s1.FontSize = 12;
ax1 = gca;
ax1.TickDir = 'out';
title(sprintf('$V_{P}$ = %0.2f km/s, S(H,$\\kappa$) = %0.2f',...
    vp(idx),maxS(idx)));

s2 = subplot(1,2,2);
contourf(K,H,sS{idx},400,'EdgeColor','none')
cbar = colorbar;
colormap(turbo)
cbar.Label.String = 'S(H,$\kappa$)';
cbar.Label.Interpreter = 'latex';
cbar.TickLabelInterpreter = 'latex';
cbar.FontSize = 12;
caxis([0 max(max(sS{idx}))])
xlabel('$\kappa$')
ylabel('Depth (km)')
s2.Position = [0.5 0.1 0.35 0.8];
ax2 = gca;
s2.FontSize = 12;
ax2.YDir = 'reverse';
ax2.TickDir = 'out';
ylim([20 70])
title(sprintf('H = %0.2f km, $\\kappa$ = %0.2f',H(x,y),K(x,y)))

sgtitle(sprintf('TA %s',stat))

% Optional: Save the output data and figure
if strcmp(saveOutput,'true')
    % Save the completed figure!
    fName = [stat,'_HkStack.png'];
    print(fullfile(outDir,fName),'-dpng','-r300');

    % Save the .xyz data to a .xyz file
    fName = [stat,'_hk.xyz'];
    fid = fopen(fullfile(outDir,fName),'w');
    [I,J] = size(sS{idx});
    for i = 1:I
        for j = 1:J
            fprintf(fid,sprintf('%0.6f %0.6f %0.6f\n',...
                    K(j,i),H(j,i),sS{idx}(j,i)));
        end
    end
    fclose(fid);

    % Save the H, k, and Station metadata to a .mat file
    stla = rf{1}.h.stla;
    stlo = rf{1}.h.stlo;
    h = H(x,y);
    k = K(x,y);
    fName = [stat,'_HkData.mat'];
    save(fullfile(outDir,fName),'stla','stlo','h','k');
end
