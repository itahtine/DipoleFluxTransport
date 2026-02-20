%   Simulates the dipole evolution due to:
%       (1) an initial synoptic map map0, and
%       (2) multiple active regions regs emerging at times tvec.
%
%   Uses precomputed DFT propagator matrices G 
%
%   Inputs:
%       regs      - Struct array with fields:
%                      .vals (flux values)
%                      .inds (pixel indices)
%       tvec      - Emergence times for each region (days).
%       map0      - Initial synoptic map (180×360).
%       G         - Propagator: G(t, 3, Npix).
%       T         - Total duration of simulation (days).
%       timestep  - Temporal resolution of G (days).
%
%   Outputs:
%       v     - Dipole magnitude time series.
%       t     - Dipole latitude (radians).
%       p     - Dipole longitude (radians).
%       vmat  - Dipole vector components for each region.

function [v,t,p,VMat] = DFTSim(regs,tvec,map0,G,T,timestep)

% AR vectors
nRegs = length(regs);
simLength = round(T/timestep);

[vsumar,tsumar,psumar] = deal(nan(simLength,nRegs));

idx1 = floor(tvec/timestep)+2;
idx1(tvec==0) = idx1(tvec==0)-1;
idx2 = simLength;

deltaT = T-tvec;
nSteps = floor(deltaT/timestep); % Number of full steps

for i = 1:nRegs
    map  = regs(i).vals;
    idx = regs(i).inds;
    
    if tvec(i) == 0
        tsteps = 0:timestep:nSteps(i)*timestep;
    else
        tsteps = 0:timestep:(nSteps(i)-1)*timestep;
        tsteps = round(tsteps + timestep-mod(tvec(i),timestep));
    end

    if nSteps(i) == 0 & deltaT(i) <= 0.5
        map = zeros(180,360);
        map(idx) =  regs(i).vals;
        
        [v,t,p] = Calc3DVectorSum(map);
        vsumar(idx2,i) = v;
        tsumar(idx2,i) = deg2rad(t);
        psumar(idx2,i) = deg2rad(p+180);
    elseif nSteps(i) == 0 & deltaT(i) > 0.5
        tsteps = round(timestep-mod(tvec(i),timestep));
        VMat = G(tsteps,:,idx);
        VMat(:,1,:) = VMat(:,1,:).*reshape(map, 1, 1, []);
        [vsumar(idx1(i):idx2,i),tsumar(idx1(i):idx2,i),psumar(idx1(i):idx2,i)] = SumDipoleVectorsRad(VMat);
    elseif tsteps(1) == 0
        VMat = G(tsteps(2:end),:,idx);
        VMat(:,1,:) = VMat(:,1,:).*reshape(map, 1, 1, []);
        [vsumar(idx1(i)+1:idx2,i),tsumar(idx1(i)+1:idx2,i),psumar(idx1(i)+1:idx2,i)] = SumDipoleVectorsRad(VMat);

        map = zeros(180,360);
        map(idx) =  regs(i).vals;
        [v,t,p] = Calc3DVectorSum(map);
        vsumar(idx1(i),i) = v;
        tsumar(idx1(i),i) = deg2rad(t);
        psumar(idx1(i),i) = deg2rad(p+180);
    elseif tsteps(1) > 0
        VMat = G(tsteps,:,idx);
        VMat(:,1,:) = VMat(:,1,:).*reshape(map, 1, 1, []);
        [vsumar(idx1(i):idx2,i),tsumar(idx1(i):idx2,i),psumar(idx1(i):idx2,i)] = SumDipoleVectorsRad(VMat);
    end
end

VMat =  reshape([vsumar;tsumar;psumar],simLength,3,[]);

% Initial vector
tsteps = round(timestep:timestep:(simLength-1)*timestep);

VMat0 = G(tsteps,:,:);
VMat0 (:,1,:) = VMat0(:,1,:).*reshape(map0, 1, 1, []);

[v0,t0,p0] = deal(zeros(simLength,1));
[v0(2:end),t0(2:end),p0(2:end)] = SumDipoleVectorsRad(VMat0);

[v,t,p] = Calc3DVectorSum(map0);
v0(1) = v;
t0(1) = deg2rad(t);
p0(1) = deg2rad(p+180);

VMat0 = reshape([v0,t0,p0],simLength,3,1);

% Sum initial vector and AR vectors

if nRegs > 0
    VMat = cat(3,VMat0,VMat);
    [v,t,p] = SumDipoleVectorsRad(VMat);
else
    VMat = VMat0;
    v = v0;
    t = t0;
    p = p0;
end



