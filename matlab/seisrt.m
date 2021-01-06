function [rdat,tdat] = seisrt(ndat,edat,baz)
% SEISRT Rotate seismograms from North and East components into
%        radial and transverse components.
%
% >> [rdat,tdat] = SEISRT(ndat,edat,baz)
%
%---Input Variables----------------------------------------------
% ndat   - North component of seismic data
% edat   - East component of seismic data
% baz    - Backazimuth measured clockwise from north (degrees)
%          from the station's location to the event
%
%---Output Variables---------------------------------------------
% rdat   - Radial component of seismic data
% tdat   - Transverse component of seismic data
%
%----------------------------------------------------------------
% Last updated 3/10/2019 by aburky@princeton.edu
%----------------------------------------------------------------

if baz < 180
    az = baz + 180;
elseif baz >= 180
    az = baz - 180;
end

m2d = [cosd(az) sind(az); -sind(az) cosd(az)];

rot = m2d*[ndat';edat'];
rdat = rot(1,:)';
tdat = rot(2,:)';