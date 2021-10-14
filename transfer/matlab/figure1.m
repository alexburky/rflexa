% figure1.m
%
% This code is meant to serve as a companion to the 2021 EduQuakes paper,
% 'Instrument response removal and the 2020 M3.1 Marlboro, New Jersey,
% earthquake.' It produces the first figure of the paper.
%
% External dependencies:
% - parseRESP
% - parsePZ
%
%--------------------------------------------------------------------------
% Last updated 10/13/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

clear,clc

% Option: display boxes with figure labels (a,b,c,etc.)
boxOn = false;

% Specify the locations of the response files
respDir = '../data/responses/';
resp{1} = 'PP.S0001.00.HHZ.resp';
resp{2} = 'PP.S0002.00.HHZ.resp';
resp{3} = 'PP.S0002.10.HNZ.resp';
resp{4} = 'SAC_PZs_R36A4_00';

% Read in the seismograms and responses
for i = 1:length(resp)
    if i < 4
        [z{i},p{i},k{i}] = parseRESP(fullfile(respDir,resp{i}));
    else
        [z{i},p{i},k{i}] = parsePZ(fullfile(respDir,resp{i}));
    end
end

% Sample rate
fs = 100;
nyq = fs/2;

% Prepare the frequency grid for Bode plot (this takes a few seconds...)
% npts = 1000;
npts = 10000000;
nfft = 2^nextpow2(npts);
nfreq = (nfft / 2) + 1;
f = linspace(0,fs,nfreq);

z1 = nonzeros(z{1});
z2 = nonzeros(z{2});
z3 = nonzeros(z{3});
z4 = nonzeros(z{4});

% Displacement
z{1} = [complex(0,0); complex(0,0); complex(0,0); z1];
z{2} = [complex(0,0); complex(0,0); complex(0,0); z2];
z{3} = [complex(0,0); complex(0,0); z3];
z{4} = [complex(0,0); complex(0,0); complex(0,0); complex(0,0); z4];

% Velocity
z{5} = [complex(0,0); complex(0,0); z1];
p{5} = p{1};
k{5} = k{1};
z{6} = [complex(0,0); complex(0,0); z2];
p{6} = p{2};
k{6} = k{2};
z{7} = [complex(0,0); z3];
p{7} = p{3};
k{7} = k{3};
z{8} = [complex(0,0); complex(0,0); complex(0,0); z4];
p{8} = p{4};
k{8} = k{4};

% Acceleration
z{9} = [complex(0,0); z1];
p{9} = p{1};
k{9} = k{1};
z{10} = [complex(0,0); z2];
p{10} = p{2};
k{10} = k{2};
z{11} = [z3];
p{11} = p{3};
k{11} = k{3};
z{12} = [complex(0,0); complex(0,0); z4];
p{12} = p{4};
k{12} = k{4};

% Construct the transfer function from poles, zeros, constant
for i = 1:12
    [b{i},a{i}] = zp2tf(z{i},p{i},k{i});
    [h{i},w{i}] = freqs(b{i},a{i},2*pi*f);
end

%--------------------------------------------------------------------------
% FIGURE 1
%--------------------------------------------------------------------------

color{1} = [0.25 0.5 0.15];
color{2} = [0 0.9 1];
color{3} = [0 0 0];
color{4} = [0.89 0.043 0.365];

% Displacement
if boxOn == true
    subplot(3,2,1)
else
    subplot(2,3,1)
end
loglog(f,abs(h{4}),'Color',color{4},'linewidth',1)
hold on
loglog(f,abs(h{1}),'Color',color{1},'linewidth',1)
loglog(f,abs(h{2}),'Color',color{2},'linewidth',1,'LineStyle','--')
loglog(f,abs(h{3}),'Color',color{3},'linewidth',1)

plot([nyq nyq],[1e0 1e14],'k--')
xlim([5e-5 fs])
ylim([5e5 1e12])
grid on
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [1e6, 1e7, 1e8 ,1e9, 1e10, 1e11, 1e12];
ax.XTick = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ylabel('Gain (Counts / m)')
if boxOn == true
    rectangle('Position',[0.85e-4,0.8e6,2.7e-4,3.5e6],'FaceColor',[1 1 1])
    text(1e-4,1.8e6,'(a)','FontSize',12)
    text(5e1,1e13,'Displacement','FontSize',14)
else
    title('Displacement')
end
legend('R36A4.00.EHZ','S0001.00.HHZ','S0002.00.HHZ','S0002.10.HNZ','Location','NorthWest')

if boxOn == true
    subplot(3,2,2)
else
    subplot(2,3,4)
end
semilogx(f,((180/pi)*angle(h{1})),'Color',color{1},'linewidth',1)
hold on
semilogx(f,((180/pi)*angle(h{2})),'Color',color{2},'linewidth',1,...
    'LineStyle','--')
semilogx(f,((180/pi)*angle(h{3})),'Color',color{3},'linewidth',1)
semilogx(f,((180/pi)*angle(h{4})),'Color',color{4},'linewidth',1)
plot([nyq nyq],[-200 200],'k--')
xlim([5e-5 fs])
ylim([-190 190])
grid on
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [-180, -135, -90, -45, 0, 45, 90, 135, 180];
ax.XTick = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ylabel('Phase ($^{\circ}$)')
if boxOn == true
    rectangle('Position',[0.85e-4,-175,2.7e-4,45],'FaceColor',[1 1 1])
    text(1e-4,-155,'(b)','FontSize',12)
else
    xlabel('Frequency (Hz)')
end

% Velocity
if boxOn == true
    subplot(3,2,3)
else
    subplot(2,3,2)
end
loglog(f,abs(h{5}),'Color',color{1},'linewidth',1)
hold on
loglog(f,abs(h{6}),'Color',color{2},'linewidth',1,'LineStyle','--')
loglog(f,abs(h{7}),'Color',color{3},'linewidth',1)
loglog(f,abs(h{8}),'Color',color{4},'linewidth',1)
plot([nyq nyq],[1e0 1e10],'k--')
xlim([5e-5 fs])
ylim([5e3 3e9])
grid on
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [1e4, 1e5, 1e6, 1e7, 1e8 ,1e9];
ax.XTick = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
if boxOn == true
    rectangle('Position',[0.85e-4,0.8e4,2.7e-4,3.5e4],'FaceColor',[1 1 1])
    text(1e-4,1.8e4,'(c)','FontSize',12)
    text(1.5e2,2e10,'Velocity','FontSize',14)
    ylabel('Gain (Counts / m/s)')
else
    title('Velocity')
end

if boxOn == true
    subplot(3,2,4)
else
    subplot(2,3,5)
end
semilogx(f,((180/pi)*angle(h{5})),'Color',color{1},'linewidth',1)
hold on
semilogx(f,((180/pi)*angle(h{6})),'Color',color{2},'linewidth',1,...
    'LineStyle','--')
semilogx(f,((180/pi)*angle(h{7})),'Color',color{3},'linewidth',1)
semilogx(f,((180/pi)*angle(h{8})),'Color',color{4},'linewidth',1)
plot([nyq nyq],[-200 200],'k--')
xlim([5e-5 fs])
ylim([-190 190])
grid on
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [-180, -135, -90, -45, 0, 45, 90, 135, 180];
ax.XTick = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
if boxOn == true
    rectangle('Position',[0.85e-4,-175,2.7e-4,45],'FaceColor',[1 1 1])
    text(1e-4,-155,'(d)','FontSize',12)
    ylabel('Phase ($^{\circ}$)')
else
    xlabel('Frequency (Hz)')
end

% Acceleration
if boxOn == true
    subplot(3,2,5)
else
    subplot(2,3,3)
end
loglog(f,abs(h{9}),'Color',color{1},'linewidth',1)
hold on
loglog(f,abs(h{10}),'Color',color{2},'linewidth',1,'LineStyle','--')
loglog(f,abs(h{11}),'Color',color{3},'linewidth',1)
loglog(f,abs(h{12}),'Color',color{4},'linewidth',1)
plot([nyq nyq],[1e0 1e10],'k--')
xlim([5e-5 fs])
ylim([5e3 1e10])
grid on
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [1e4, 1e5, 1e6, 1e7, 1e8 ,1e9];
ax.XTick = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
if boxOn == true
    rectangle('Position',[0.85e-4,0.8e4,2.7e-4,3.5e4],'FaceColor',[1 1 1])
    text(1e-4,1.8e4,'(e)','FontSize',12)
    text(0.65e2,6e10,'Acceleration','FontSize',14)
    xlabel('Frequency (Hz)')
    ylabel('Gain (Counts / (m/s$^{2}$))')
else
    title('Acceleration')
end

if boxOn == true
    subplot(3,2,6)
else
    subplot(2,3,6)
end
semilogx(f,((180/pi)*angle(h{9})),'Color',color{1},'linewidth',1)
hold on
semilogx(f,((180/pi)*angle(h{10})),'Color',color{2},'linewidth',1,...
    'LineStyle','--')
semilogx(f,((180/pi)*angle(h{11})),'Color',color{3},'linewidth',1)
semilogx(f,((180/pi)*angle(h{12})),'Color',color{4},'linewidth',1)
plot([nyq nyq],[-200 200],'k--')
xlim([5e-5 fs])
ylim([-190 190])
grid on
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [-180, -135, -90, -45, 0, 45, 90, 135, 180];
ax.XTick = [1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
xlabel('Frequency (Hz)')
if boxOn == true
    rectangle('Position',[0.85e-4,-175,2.7e-4,45],'FaceColor',[1 1 1])
    text(1e-4,-155,'(f)','FontSize',12)
    ylabel('Phase ($^{\circ}$)')
end

if boxOn == true
    set(gcf,'Position',[0 0 600 750])
else
    set(gcf,'Position',[0 0 880 475])
end

