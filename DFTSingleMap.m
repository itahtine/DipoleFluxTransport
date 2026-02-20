%   Simulates the time evolution of the dipole vector for a single map
%   or a struct array of active regions using a precomputed DFT propagator.
%
%   Inputs:
%       regs       - Either:
%                    (a) struct array with fields:
%                        .vals   (flux values), 
%                        .inds   (linear indices in 180×360 map), or
%                    (b) 3D array of maps (180×360×N).
%
%       VectorMap  - Propagator matrix, size:
%                    (lat, lon, 3 components, time).
%
%       SimLength  - Number of timesteps (Carrington rotations or days,
%                    depending on propagator resolution).
%
#   Outputs:
%       vsumar - Dipole magnitude time series for each AR/map.
%       tsumar  - Dipole latitude (degrees).
%       psumar  - Dipole longitude (degrees).


function [vsumar,tsumar,psumar] = DFTSingleMap(regs,VectorMap,simLength)
% simLength in Carrington rotations

% AR vectors
len = round(simLength);
VecMap = VectorMap(:,:,:,1:len);

VecMat = permute(reshape(VecMap, [], 3, len), [3 2 1]);
VecMat(:,2,:) = deg2rad(VecMat(:,2,:));
VecMat(:,3,:) = deg2rad(VecMat(:,3,:));

[vsumar,tsumar,psumar] = deal(nan(len,1));

if isa(regs,'struct')
    for i = 1:length(regs)
        map  = regs(i).vals;
        idx = regs(i).inds;

        VMat = VecMat(:,:,idx);
        VMat(:,1,:) = VMat(:,1,:).*reshape(map, 1, 1, []);
        [vsumar(:,i),tsumar(:,i),psumar(:,i)] = SumDipoleVectorsRad(VMat);
    end
else
    for i = 1:size(regs,3)
        map = regs(:,:,i);
        idx = find(map);
        map = map(idx);

        VMat = VecMat(:,:,idx);
        VMat(:,1,:) = VMat(:,1,:).*reshape(map, 1, 1, []);
        [vsumar(:,i),tsumar(:,i),psumar(:,i)] = SumDipoleVectorsRad(VMat);
    end
end

tsumar = rad2deg(tsumar);

psumar = rad2deg(psumar);
