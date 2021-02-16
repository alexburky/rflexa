function [nlnm, nhnm] = petersonNoiseModel(units)
% PETERSONNOISEMODEL Returns structures containing the Peterson (1993) New
%                    Low and New High Noise Models, for a desired unit.
% Adapted from:
% https://github.com/ytamama/GuyotSeismology/blob/master/pnm.m
%
% >> [nlnm, nhnm] = petersonNoiseModel(units)
%
%---Input Variables--------------------------------------------------------
% units - Desired output units, 'displacement', 'velocity', or
%         'acceleration'
%
%---Output Variables-------------------------------------------------------
% nlnm  - Structure containing the New Low Noise Model
% nhnm  - Structure containing the New High Noise Model
%
%--------------------------------------------------------------------------
% Last updated 2/15/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Path to noise model CSV files
lncsv = '/Users/aburky/IFILES/MODELS/NOISE/nlnm.csv';
hncsv = '/Users/aburky/IFILES/MODELS/NOISE/nhnm.csv';

ln = readtable(lncsv,'Delimiter',',');
hn = readtable(hncsv,'Delimiter',',');

% Equation for the noise model depends on desired units
if strcmp(units,'displacement')
    % Displacement noise model: A + Blog10(P) + 20log10(P^2/(4pi^2))
    nlnm.d = ln.A+ln.B.*log10(ln.P)+20*log10((ln.P.^2)/(4*pi*pi));
    nhnm.d = hn.A+hn.B.*log10(hn.P)+20*log10((hn.P.^2)/(4*pi*pi));
    nlnm.units = 'dB (m$^{2}$/Hz)';
    nhnm.units = 'dB (m$^{2}$/Hz)';
elseif strcmp(units,'velocity')
    % Velocity noise model: A + Blog10(P) + 20log10(P/(2pi))
    nlnm.d = ln.A+ln.B.*log10(ln.P)+20*log10((ln.P)/(2*pi));
    nhnm.d = hn.A+hn.B.*log10(hn.P)+20*log10((hn.P)/(2*pi));
    nlnm.units = 'dB ((m/s)$^2$/Hz)';
    nhnm.units = 'dB ((m/s)$^2$/Hz)';
elseif strcmp(units,'acceleration')
    % Acceleration noise model: A + Blog10(P)
    nlnm.d = ln.A+ln.B.*log10(ln.P);
    nhnm.d = hn.A+hn.B.*log10(hn.P);
    nlnm.units = 'dB ((m/s$^2$)$^2$/Hz)';
    nhnm.units = 'dB ((m/s$^2$)$^2$/Hz)';
end

% X Axis Units - Frequency
nlnm.f = 1./ln.P;
nhnm.f = 1./hn.P;