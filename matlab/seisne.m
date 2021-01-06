function [ndat,edat] = seisne(dat1,dat2,az)
% SEISNE Rotate seismograms from arbitrary components to North and
%        East components.
%
% >> [ndat,edat] = SEISNE(dat1,dat2,az)
%
%---Input Variables----------------------------------------------
% dat1   - Channel 1 of seismic data (BH1,HH1,etc.) - North
% dat2   - Channel 2 of seismic data (BH2,HH2,etc,) - East
% az     - Azimuth measured clockwise from north, use the 'CMPAZ'
%          value from the SAC header of the BH1 channel
%
%---Output Variables---------------------------------------------
% ndat   - North component of seismic data
% edat   - East component of seismic data
%
%----------------------------------------------------------------
% Last updated 3/10/2019 by aburky@princeton.edu
%----------------------------------------------------------------

m2d = [cosd(az) -sind(az); sind(az) cosd(az)];

rot = m2d*[dat1';dat2'];
ndat = rot(1,:)';
edat = rot(2,:)';