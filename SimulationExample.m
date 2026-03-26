% This code should produce SimulationExample.png

load('SHARPRegsNoRepeat.mat')
load('VectorMapFullDiff350Mer11CR.mat')

dA = 4*pi*(696340*1000*100)^2/(180*360);

fluxes = zeros(length(regs),1);
for i = 1:length(fluxes)
    fluxes(i) = sum(abs(regs(i).vals));
end

% Indices from largest to smallest flux
[~,idx] = sort(fluxes,'descend');

% Simulate the evolution of the largest AR of SC24
aridx = idx(1);
ar = zeros(180,360);
ar(regs(aridx).inds) = regs(aridx).vals;

simLength = 50;
[v,t,p] = DFTSingleMap(ar,VectorMap,simLength);

% Plot
figure
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(1:360,linspace(-1,1,180),ar)
set(gca,'YDir','normal')
caxis([-50 50])
yticks(-1:0.5:1)
yticklabels(asind(yticks))
xticks(0:60:360)
ylabel('Latitude')
xlabel('Longitude')
title(sprintf('SHARP%i',regs(aridx).sharp))

nexttile
plot(v*dA)
xlabel('Time (CR)')
ylabel('Magnetic flux (G)')
title('Dipole flux')

nexttile
plot(v.*abs(sind(t))*dA)
xlabel('Time (CR)')
ylabel('Magnetic flux (G)')
title('Axial component')

nexttile
plot(v.*cosd(t)*dA)
xlabel('Time (CR)')
ylabel('Magnetic flux (G)')
title('Equatorial component')

nexttile
plot(t)
xlabel('Time (CR)')
ylabel('Latitude')
ylim([-90 90])
title('Latitude')

nexttile
scatter(1:simLength,p+180,'filled')
xlabel('Time (CR)')
ylabel('Carrington longitude')
ylim([0 360])
title('Longitude')

for i = 2:6
    nexttile(i)
    xlim([0 simLength])
end
%%