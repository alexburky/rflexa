function [x,y,z] = latlon2cart(lat,lon,rad)
% LATLON2CART Converts latitude, longitude, and spherical Earth radius to
%             Cartesian coordinates. Useful for making surface plots of
%             spherical sections.
%
% >> [x,y,z] = latlon2cart(lat,lon,rad)
%
%---Input Variables--------------------------------------------------------
% lat - latitude value
% lon - longitude value
% rad - radius value
%
%---Output Variables-------------------------------------------------------
% x   - 3-D x Cartesian coordinate
% y   - 3-D y Cartesian coordinate
% z   - 3-D z Cartesian coordinate
%--------------------------------------------------------------------------
% Last updated 5/26/2021 by aburky@princeton.edu
%--------------------------------------------------------------------------

% Convert latitude, longitude to theta, phi
th = (90-lat)*(pi/180);
ph = (180+lon)*(pi/180);

x = rad.*cos(ph).*sin(th);
y = rad.*sin(ph).*sin(th);
z = rad.*cos(th);