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