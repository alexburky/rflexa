function bodePlot(file,fileType)
% BODEPLOT Make a Bode plot of the amplitude and phase response of a
%          seismometer from a SAC_PZs file or a RESP file
%
% >> bodePlot(file,fileType)
%
%---Input Variables--------------------------------------------------------
% file     - full path to the file containing the instrument response data
% fileType - type of file, either 'sacpz' or 'resp'
%
%--------------------------------------------------------------------------
% Last updated 3/02/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

if strcmp(fileType,'sacpz')
    [z,p,k] = parsePZ(file);
elseif strcmp(fileType,'resp')
    [z,p,k] = parseRESP(file);
else
    error(['Invalid fileType. Currently supported options are ',...
           '''sacpz'' or ''resp'''])
end

% Sample rate (get this from the data)
fs = 100;
nyq = fs/2;

npts = 10000000;
nfft = 2^nextpow2(npts);
nfreq = (nfft / 2) + 1;
f = linspace(0,fs,nfreq);

[b,a] = zp2tf(z,p,k);
[h,w] = freqs(b,a,2*pi*f);

% Amplitude response plot
figure
subplot(1,2,1)
loglog(f,abs(h),'r','linewidth',1)
hold on
plot([nyq nyq],[1e0 1e10],'k--')
% Find an appropriate xmin
[~,idx] = min(abs(abs(h)-5e3));
xlim([f(idx)*0.5 fs])
ylim([5e3 3e9])
grid on
title('Amplitude Response')
ylabel('Gain (Counts / m/s)')
xlabel('Frequency (Hz)')
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';

% Phase response plot
subplot(1,2,2)
semilogx(f,((180/pi) * angle(h)),'b','linewidth',1)
hold on
plot([nyq nyq],[-200 200],'k--')
xlim([f(idx)*0.5 fs])
ylim([-181 180])
grid on
title('Phase Response')
ylabel('Phase ($^{\circ}$)')
xlabel('Frequency (Hz)')
ax = gca;
ax.FontSize = 12;
ax.TickDir = 'out';
ax.YTick = [-180, -135, -90, -45, 0, 45, 90, 135, 180];
set(gcf,'Position',[650 0 600 250])

end

%%

% t1 = ch(1).StartDate;
% % Use the start date plus a very short time as the time window
% t2 = '2014-08-28 20:00:00.000';
% % t2 = ch(1).EndDate;
% 
% channel = 'BHZ';
% 
% re = irisFetch.Resp(network,station,location,channel,t1,t2);
% 
% % Loop over all channels and save a RESP file for each one!
% 
% R = regexp(re,'\n','split');
% % Get the number of poles and the number of zeroes
% nzLine = find(contains(R,'Number of zeroes:'));
% tmp = split(R{nzLine});
% nz = str2double(tmp{end});
% npLine = find(contains(R,'Number of poles:'));
% tmp = split(R{npLine});
% np = str2double(tmp{end});
% 
% % Get the line where the zeroes begin
% zLine = find(contains(R,'Complex zeroes'));
% for i = 2:(nz + 1)
%     tmp = split(R{zLine + i});
%     z(i-1) = complex(str2double(tmp{3}),str2double(tmp{4}));
% end
% 
% % Get the line where the poles begin
% pLine = find(contains(R,'Complex poles'));
% for i = 2:(np + 1)
%     tmp = split(R{pLine + i});
%     p(i-1) = complex(str2double(tmp{3}),str2double(tmp{4}));
% end
% 
% % Get the poles and zeroes constant (A0 * sensitivity)
% a0Line = find(contains(R,'A0 normalization factor'));
% tmp = split(R{a0Line});
% a0 = str2double(tmp{end});
% sensLine = find(contains(R,'Sensitivity'));
% tmp = split(R{sensLine(end)});
% sens = str2double(tmp{end});
% k = a0 * sens;
% 
% % Try saving the instrument response information to a file
% fID = fopen('RESP.TEST','w');
% fprintf(fID,re);
% fclose(fID);

% End date of one channel is equal to the start date of the next.
% what to do about this?
