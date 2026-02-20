
%   Computes the total solar dipole vector using the vector-sum method
%   (Tähtinen et al. 2024, 2026). Each pixel of the magnetogram is treated
%   as a vector in spherical coordinates with magnitude equal to its total
%   magnetic flux and direction given by its surface location.
%
%   Input:
%       MAP  - 3D array of size (lat, lon, time or Nmaps), containing
%              magnetic field values (gauss). Can represent a time series
%              or a collection of maps.
%
%   Output:
%       MAG   - Dipole vector magnitude for each map (units of flux).
%       THETA - Dipole latitude in degrees.
%       PHI   - Dipole longitude in degrees. Add 180 for Carrington longitude.

function [mag,theta,phi] = Calc3DVectorSum(map)

theta = asin(repmat(linspace(-1+1/size(map,1),1-1/size(map,1),size(map,1))',1,size(map,2),size(map,3)));
phi = deg2rad(repmat(linspace(-180+size(map,2)/360/2,180-size(map,2)/360/2,size(map,2)),size(map,1),1,size(map,3)));
r = map;

[x,y,z] = sph2cart(phi,theta,r);

x = reshape(x,[size(map,1)*size(map,2) size(map,3)]);
y = reshape(y,[size(map,1)*size(map,2) size(map,3)]);
z = reshape(z,[size(map,1)*size(map,2) size(map,3)]);

x = sum(x,'omitnan')';
y = sum(y,'omitnan')';
z = sum(z,'omitnan')';

[phi,theta,mag] = cart2sph(x,y,z);

mag = abs(mag);
theta = rad2deg(theta);
phi = rad2deg(phi);

idx = mag==0;
mag(idx) = NaN;
theta(idx) = NaN;

phi(idx) = NaN;
