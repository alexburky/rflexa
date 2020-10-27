function [data] = transfer(data,delta,freqlimits,units,pzfile)
% TRANSFER Remove the instrument response from a seismogram using the
%          poles and zeros contained in an external SAC style Pole Zero
%          file. Pre-filters the data with a cosine filter defined by
%          freqlims, and outputs with user-specified units.
%
% >> [data] = transfer(data,delta,freqlimits,units,pzfile)
% 
%---Input Variables--------------------------------------------------------
% data       - vector containing the input seismic data
% delta      - sample rate of the data (s)
% freqlimits - vector containing 4 frequencies in order of increasing
%              frequency. These frequencies specify the corners of a
%              cosine filter applied to the data to stabilize deconvolution
% units      - output units, 'displacement', 'velocity', or 'acceleration'
% pzfile     - full path to the SAC_PZs file
%
%---Output Variables-------------------------------------------------------
% data       - vector containing the seismic data with the instrument
%              response removed
%
%--------------------------------------------------------------------------
% Last updated 10/27/2020 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Do some input validation...
if length(freqlimits) ~= 4
    error(['Incorrect size of parameter ''freqlimits''',newline,...
           'Should have 4 elements']);
end

fl = sort(freqlimits);

if sum(strcmp(units,{'displacement','velocity','acceleration'})) == 0
    error(['Incorrect output unit option.',newline,'Currently supported'...
          ' options are ''displacement'', ''velocity'', or '...
          '''acceleration''']);
end

% Get the poles, zeros, and constant
[z,p,k] = parsePZ(pzfile);

% Displacement option: 3 zeros
% Velocity option:     2 zeros
% Acceleration option: 1 zero

% FFT parameters
npts = length(data);
nfft = 2^nextpow2(npts);
dfreq = 1/(nfft*delta);
nfreq = (nfft/2) + 1;
nyq = (1/delta)*0.5;

% Construct the transfer function from the poles and zeros
f = linspace(0,nyq,nfreq);
[b,a] = zp2tf(z,p,k);
[h,w] = freqs(b,a,2*pi*f);

% Invert the transfer function
for i = 2:nfreq
    denr = real(h(i))^2 + imag(h(i))^2;
    denr = 1.0/denr;
    if denr < 1e-37
        h(i) = complex(0,0);
    else
        h(i) = complex(real(h(i))*denr,-imag(h(i))*denr);
    end
end

% Apply the cosine filter (make this a function?)
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
    h(i) = h(i) * fac;
end

% Fourier transform the data
data_fft = fft(data,nfft);

% Multiply Fourier transformed data by the inverted transfer function
for i = 2:nfreq
    tempR = real(h(i))*real(data_fft(i)) - imag(h(i))*imag(data_fft(i));
    tempI = real(h(i))*imag(data_fft(i)) + imag(h(i))*real(data_fft(i));
    data_fft(i) = complex(tempR,tempI);
    
    if i < (nfreq)
        j = nfft - i + 2;
        data_fft(j) = complex(tempR,-tempI);
    end
end

data_fft(1) = complex(0,0);

data_fft(nfft) = complex(sqrt(real(data_fft(nfft))*...
                 real(data_fft(nfft)) + imag(data_fft(nfft))*...
                 imag(data_fft(nfft))),0);
             
% Transform back into the time domain
data = ifft(data_fft,nfft);

data = data(1:npts);
