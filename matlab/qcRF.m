function [rfQuality] = qcRF(rf,t)
% QCRF Calculate quality ratio for receiver function data.
%      For a complete description of the quality ratio
%      algorithm, see Burky et. al, 2020.
%
% >> [rfQuality] = QCRF(rf)
%
%---Input Variables-------------------------------------------
% rf         - Vector of receiver function data
% t          - Vector of corresponding time data (assumes RF
%              has been shifted to align the maximum at 10s)
%
%---Output Variables------------------------------------------
% rfQuality  - Scalar receiver function quality ratio
%
%-------------------------------------------------------------
% Last updated 1/5/2021 by aburky@princeton.edu
%-------------------------------------------------------------

% Take derivative of receiver function
dt = t(2) - t(1);
drf = diff(rf/dt);

% Find index of first RF value above threshold
idx1 = find(abs(rf)>1e-4,1);

% Find largest peak in first 13 seconds (main P arrival)
bidx = 1;
eidx = round(13/dt);
[pk,idx] = findpeaks(abs(rf(bidx:eidx)),'NPeaks',1,'SortStr','descend');
idx = idx + 10;

% Metric 1: Determine where RF goes negative
idx2 = find(rf(idx:end)<0,1);
idx2 = idx + idx2;

% Metric 2: Determine where first derivative changes sign
idx3 = find(drf(idx:end)>0,1);
idx3 = idx + idx3;

% Check which occurs first
if idx3 < idx2
    idx2 = idx3;
end

% Integral of main P-arrival
p_e = trapz(abs(rf(idx1:idx2)));

% Integral of entire RF
rf_e = trapz(abs(rf));

% Ratio of integrals
p_rf = p_e/rf_e;

% Check first impulse - is it negative?
[pk_r] = findpeaks(rf(bidx:eidx),'NPeaks',1,'SortStr','descend');

% Check for negative first impulse
if pk_r ~= pk
    p_rf = 0;
end

% Check for extremely small amplitudes
if pk < 1e-3
    p_rf = 0;
end

% Final Quality Ratio
rfQuality = p_rf;
