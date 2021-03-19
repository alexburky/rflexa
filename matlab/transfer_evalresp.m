function [data] = transfer_evalresp(data,delta,freqlimits,amp,phase)
% TRANSFER_EVALRESP Remove the instrument response from a seismogram using
%                   the AMP and PHASE files output by the exteranl C
%                   program, EVALRESP. This function requires that you have
%                   already run EVALRESP to produce these files.
%
% >> [data] = transfer_evalresp(data,delta,freqlimits,amp,phase)
%
%---Input Variables--------------------------------------------------------
% data       - vector containing the input seismic data
% delta      - sample rate of the data (s)
% freqlimits - vector containing 4 frequencies in order of increasing
%              frequency. These frequencies specify the corners of a
%              cosine filter applied to the data to stabilize deconvolution
% units      - output units, 'displacement', 'velocity', or 'acceleration'
% amp        - full path to the AMP* file output by EVALRESP
% phase      - full path to the PHASE* file output by EVALRESP
%
%---Output Variables-------------------------------------------------------
% data       - vector containing the seismic data with the instrument
%              response removed
%
%--------------------------------------------------------------------------
% Last updated 3/18/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

% To do: maybe we can recast this program in a way that mimics SAC and
%        does a call to EVALRESP? Otherwise we can't have 'units' as an
%        input option

% Do some input validation...
if length(freqlimits) ~= 4
    error(['Incorrect size of parameter ''freqlimits''',newline,...
           'Should have 4 elements']);
end

fl = sort(freqlimits);

% Read in the transfer function output by EVALRESP
resp.a = readmatrix(fullfile(amp),'FileType','text');
resp.p = readmatrix(fullfile(phase),'FileType','text');

resp.f = resp.a(:,1);
resp.a = resp.a(:,2);
resp.p = resp.p(:,2);

resp.tf = complex(resp.a.*cosd(resp.p),resp.a.*sind(resp.p));
resp.tf = transpose(resp.tf);

% FFT parameters
npts = length(data);
nfft = 2^nextpow2(npts);
dfreq = 1/(nfft * delta);
nfreq = (nfft / 2) + 1;
nyq = (1/delta)*0.5;

% Initialize frequency grid
a = ones(1,nfreq);
xresp = complex(a,0);
f = linspace(0,nyq,nfreq);

% Interpolate the transfer function from EVALRESP
iSpec = interp1(transpose(resp.f),transpose(resp.tf),f);

% Take the Fourier Transform of the data
data_fft = fft(data,nfft);

% Invert the transfer function spectrum
for i = 2:nfreq
    denr = real(iSpec(i))^2 + imag(iSpec(i))^2;
    denr = 1.0/denr;
    if denr < 1e-37
        iresp(i) = complex(0,0);
    else
        iresp(i) = complex(real(iSpec(i))*denr,-imag(iSpec(i))*denr);
    end
end

% Blend the spectra
for i = 2:nfreq
    tempr = real(xresp(i))*real(iresp(i)) - imag(xresp(i))*imag(iresp(i));
    tempi = real(xresp(i))*imag(iresp(i)) + imag(xresp(i))*real(iresp(i));
    iiresp(i) = complex(tempr,tempi);
end

% Apply the taper to the transfer function spectrum
for i = 2:nfreq
    freq = (i-1) * dfreq;
    if freq < fl(1)
        fac = 0.0;
    elseif freq >= fl(1) && freq <= fl(2)
        fac = 0.5 * (1 - cos(pi*(freq - fl(1)) / (fl(2) - fl(1))));
    elseif freq >= fl(3) && freq <= fl(4)
        fac = 0.5 * (1 + cos(pi*(freq - fl(3)) / (fl(4) - fl(3))));
    elseif freq > fl(4)
        fac = 0.0;
    else
        fac = 1.0;
    end
    factor = fac;
    tresp(i) = iiresp(i) * factor;
end

% Multiply transformed data by the transfer operator
for i = 2:nfreq
    tempR = real(tresp(i)) * real(data_fft(i)) - imag(tresp(i)) * imag(data_fft(i));
    tempI = real(tresp(i)) * imag(data_fft(i)) + imag(tresp(i)) * real(data_fft(i));
    data_fft_t(i) = complex(tempR,tempI);
    
    if i < (nfreq)
        j = nfft - i + 2;
        data_fft_t(j) = complex(tempR,-tempI);
    end
end

data_fft_t(1) = complex(0,0);

data_fft_t(nfft) = complex(sqrt(real(data_fft_t(nfft))*real(data_fft_t(nfft))+...
                           imag(data_fft_t(nfft))*imag(data_fft_t(nfft))),0);
                  
data_ifft = ifft(data_fft_t,nfft);

% data = data_ifft;

data = data_ifft(1:npts);

