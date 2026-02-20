%   Runs a collection of runtime tests.
%   These tests reproduce the timing comparisons presented in Table 1 of
%   the paper "Ultra-fast simulations of the solar dipole and open flux".
%
%   Included benchmarks:
%       • Single active region evolution (4 years)
%       • Daily vs Carrington resolution DFT
%       • Full Solar Cycle 24 dipole simulation
%       • Hindcast of 180 days

%% HMI synoptic maps and SHARPs
load('HMIStackFilled.mat')
load('SHARPRegsNoRepeat.mat')
regs24 = regs;
idx = [regs24.rot]' > 2097 & [regs24.rot]' < 2225;
regs24 = regs24(idx);
[~,idx] = sort([regs24.t]','ascend');
regs24 = regs24(idx);

%% Propagators
load('VectorMapDiff350Mer11.mat')
VectorMapCar = VectorMap;
load('VectorMapDailyResDiff350Mer11.mat')

%% Single AR, DFT, 54 rotations (= 4 years), Carrington res
% Type of active region data: struct

map0 = zeros(180,360);
mer = 11;
dif = 350;
regs = regs24;
rot0 = 0;
rot1 = 54;

fluxes = [];
for i = 1:length(regs)
    fluxes(i) = sum(abs(regs(i).vals));
end

% Largest active region of SC24
[~,maxidx] = max(fluxes);
reg = regs(maxidx);
reg.t = 0;

simLength = 54; % Rotations

f = @()DFTSingleMap(reg,VectorMapCar,round(simLength));
tARCar = timeit(f);
disp(seconds(tARCar))
%% Single AR, DFT, 54 rotations (= 4 years), Daily res
% Type of active region data: struct
simLength = 54*27.2753; % Days

f = @()DFTSingleMap(reg,VectorMap,round(simLength));
tARDaily = timeit(f);
disp(seconds(tARDaily))

%% SC24
map0 = HMIStack(:,:,1);
simLength = 128*27.2753;
timestep = 27.2753;
tvec = days(seconds([regs24.t]'));
idx = find(tvec > 0 & tvec < simLength);
regs = regs24(idx);
tvec = tvec(idx);

f = @() WrapperFull(regs,tvec,map0,VectorMap,simLength,timestep);
tSC24DFT = timeit(f);
disp(seconds(tSC24DFT))

%% 180 day hindcast

simLength = 180;
maps = HMIStack(:,:,1:128);

timestep = simLength;
tvec = days(seconds([regs24.t]'));
idx = find(tvec > 0 & tvec < simLength);
regs = regs24(idx);
tvec = tvec(idx);

f = @() WrapperHind(regs,tvec,maps,VectorMap,simLength,timestep);
thind=timeit(f);
disp(seconds(thind))

%% SFT AR
map0 = zeros(180,360);
mer = 11;
dif = 350;
regs = regs24;
rot0 = 0;
rot1 = 54;

% Largest active region of SC24
[~,maxidx] = max(fluxes);
reg = regs(maxidx);
reg.t = 0;

f = @()SFTWrapper(reg,dif,mer,rot0,rot1,map0);
tARSFT = timeit(f);
disp(seconds(tARSFT))

%% SFT SC24
map0 = HMIStack(:,:,1);
mer = 11;
dif = 350;
regs = regs24;
rot0 = 2097;
rot1 = 2224;

f = @()SFTWrapper(regs,dif,mer,rot0,rot1,map0);
tSC24SFT = timeit(f);
disp(seconds(tSC24SFT))

%% SFT, full cycle, individual regions, Estimate
tpercr = tSC24SFT/128;
arrots = [regs24.rot]';
T = 0;

for i = 1:length(arrots)
    T = T+(2224-arrots(i)+1)*tpercr;
end
tsft = T;

%% SFT, Hindcast, Estimate
tperday = tSC24SFT/128/27.2753;
T = 180*tperday*128;
thindsft = T;

%%
function [v,t,p,VMat] = WrapperFull(regs,tvec,map0,VectorMap,simLength,timestep)
VecMap = VectorMap(:,:,:,1:round(simLength));
VecMat = permute(reshape(VecMap, [], 3, round(simLength)), [3 2 1]);
VecMat(:,2,:) = deg2rad(VecMat(:,2,:));
VecMat(:,3,:) = deg2rad(VecMat(:,3,:));

[v,t,p,VMat] = DFTSim(regs,tvec,map0,VecMat,simLength,timestep);

t = rad2deg(t);
p = rad2deg(p);
VMat(:,2:3,:) = rad2deg(VMat(:,2:3,:));
end

function [v,t,p,VMat] = WrapperHind(regs,tvec,maps,VectorMap,simLength,timestep)
VecMap = VectorMap(:,:,:,1:round(simLength));
VecMat = permute(reshape(VecMap, [], 3, round(simLength)), [3 2 1]);
VecMat(:,2,:) = deg2rad(VecMat(:,2,:));
VecMat(:,3,:) = deg2rad(VecMat(:,3,:));

for i = 1:128
    map0 = maps(:,:,i);
    [v,t,p,VMat] = DFTSim(regs,tvec,map0,VecMat,simLength,timestep);
end

t = rad2deg(t);
p = rad2deg(p);
end

function [vsft,tsft,psft] = SFTWrapper(regs,dif,mer,rot0,rot1,map0)
[~,SFTMaps] = sft_sim_lin(regs,dif,mer,rot0,rot1,map0);
[vsft,tsft,psft] = Calc3DVectorSum(SFTMaps);

end
