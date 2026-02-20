
%   Sums N dipole vectors given in spherical coordinates and returns the
%   combined dipole vector in spherical form (magnitude, theta, phi).
%
%   Input:
%       VECMAT - Array of size (T, 3, N), where:
%                T = number of time steps,
%                3 = [magnitude, theta, phi],
%                N = number of dipole vectors to sum.
%                Angles must be in radians.
%
%   Output:
%       MAG   - Resulting combined dipole magnitude (T × 1).
%       THETA - Resulting dipole latitude in radians.
%       PHI   - Resulting dipole longitude in radians.

function [mag,theta,phi] = SumDipoleVectorsRad(VecMat)

if length(size(VecMat)) ~= 3 || size(VecMat,2) ~= 3
    error('Matrix must have dimensions (t,3,N), where t is time, 3 corresponds to components of the vector sum (mag,theta,phi), and N is the number of vectors.')
end

if size(VecMat,3) == 1
    error('VecMat contains only one dipole vector')
end

% Convert to cartesian
r = VecMat(:,1,:);            
theta = VecMat(:,2,:);
phi = VecMat(:,3,:);

[x,y,z] = sph2cart(phi, theta, r);

x = sum(x, 3, 'omitnan');
y = sum(y, 3, 'omitnan');
z = sum(z, 3, 'omitnan');

% Convert back to spherical
[phi, theta, mag] = cart2sph(x, y, z);

mag = abs(mag);

%%

