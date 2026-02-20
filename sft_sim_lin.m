%SFT simulation
%regs=active regions
%diff=diffusion coefficient
%mer=peak meridional circulation speed
%rot0=first rotation
%rot1=last rotation
%br0=initial field
%res=field at the end of last rotation
%simArray=map at the end of each rotation

function [res,simArray] = sft_sim_lin(regs,diff,mer,rot0,rot1,br0)

simArray = nan(180,360,rot1-rot0+1);
%first and last rotation

%diffusion
eta=diff/6.96e5^2;

%meridional circulation
v0=mer/6.96e8;

%create coordinate matrices
nlat=size(br0,1);
nlon=size(br0,2);
sinlatc=2./nlat*((1:nlat) - 0.5) - 1;
lonc=2.*pi/nlon*((1:nlon) - 0.5);
[~, sthc]=meshgrid(lonc,cos(asin(sinlatc)));
sthcth=[sthc,sthc(:,1)];
sinlatg=2./(nlat)*((1:(nlat+1))-1)-1;
[~, sthg]=meshgrid(lonc,cos(asin(sinlatg)));
[~, slat]=meshgrid(lonc,sinlatc);
slatcth=[slat,slat(:,1)];
latg=asin(sinlatg);
[~, latph]=meshgrid(lonc,latg);
dlon=lonc(2)-lonc(1);
dslat=sinlatc(2)-sinlatc(1);
fnet=sum(br0(:)*dlon*dslat);

%remove flux imbalance from initial map
br0=br0-fnet/(nlat*nlon*dlon*dslat);

%compute vector potential
aph=zeros(nlat+1,nlon);
ath=-sthc.*cumsum(br0,2)*dlon;
ath=[zeros(nlat,1),ath];

%differential rotation
omA=0.18/180*pi/86400; 
omB=-2.396/180*pi/86400.;
omC=-1.787/180*pi/86400.;
om=omA + omB*slatcth.^2 + omC*slatcth.^4;
vph_om = om.*sthcth;

%meridional circulation
Du = 0.041;
p = 2.33;
vth_mf = -Du.*sin(latph).*(1-sin(latph).^2).^(p/2);
vth_mf = vth_mf*v0/max(abs(vth_mf(:)));

%diffusion
hphmin=min(abs(sthc(:)*dlon));
hthmin=min(latg(2:end)-latg(1:end-1));
dt_eta=min(hphmin^2/eta,hthmin^2/eta);

%calculate time step
t_mf=abs((latph(2:end,:) - latph(1:end-1,:))./vth_mf(1:end-1,:));
dt_mf=min(t_mf(:));
dt_om=min(abs(sthcth(:)*dlon./vph_om(:)));
dt=min([dt_eta,dt_mf,dt_om])*0.2;

%fit time step to rotation length
ndt=round(27.2753*86400/dt);
dt=27.2753*86400/ndt;

brc=zeros(nlat+2,nlon+2);%radial field
nrot=rot1-rot0;%number of rotations

regcount=1;%active region counter
dtcount=0;%time step counter (for decay term)

rotCount = 0;

%run until end of last rotation with time step dt
for t=0:dt:nrot*27.2753*24*3600
    %compute radial field from vector potential
    brc(2:end-1,2:end-1)=-(sthg(2:end,:).*aph(2:end,:)-sthg(1:end-1,:)...
        .*aph(1:end-1,:))/dslat-(ath(:,2:end)-ath(:,1:end-1))/dlon./sthc;

    %insert active regions
    while regcount<=length(regs) && t>=regs(regcount).t  
        br1=brc(2:end-1,2:end-1);%take existing field
    
        % insert region
        br1(regs(regcount).inds)=br1(regs(regcount).inds)+regs(regcount).vals;

        brc(2:end-1,2:end-1)=br1;%copy back to brc
        regcount=regcount+1;%increase region counter
    end
    
    dtcount=dtcount+1;%increase time step counter
    
    %compute vector potential
    aph=zeros(nlat+1,nlon);
    ath=[zeros(nlat,1),-sthc.*cumsum(brc(2:end-1,2:end-1),2)*dlon];
    
    %deal with boundary conditions
    brc(2:end-1,1)=brc(2:end-1,end-1);
    brc(2:end-1,end)=brc(2:end-1,2);
    brc(1,2:end-1)=brc(2,2:end-1);
    brc(end,2:end-1)=brc(end-1,2:end-1);
    
    %calculate averages between pixels
    brth=0.5*(brc(2:end-1,1:end-1)+brc(2:end-1,2:end));
    brph=0.5*(brc(1:end-1,2:end-1)+brc(2:end,2:end-1));
    emfth=-vph_om.*brth;%differential rotation
    emfph=vth_mf.*brph;%meridional circulation
    
    %change of vector potential during time step
    emfth=emfth+eta./sthcth.*(brc(2:end-1,2:end)-brc(2:end-1,1:end-1))/dlon;
    emfph=emfph+eta*sthg.*(brc(2:end,2:end-1)-brc(1:end-1,2:end-1))/dslat;
    
    %compute new vector potential
    ath=ath-dt*emfth;
    aph=aph-dt*emfph;

    if t >= rotCount*27.2753*24*3600
        rotCount = rotCount + 1;
        simArray(:,:,rotCount) = brc(2:end-1,2:end-1);
    end
end

res=brc(2:end-1,2:end-1);
simArray(:,:,end) = res;
end