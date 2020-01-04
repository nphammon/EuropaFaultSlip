%function [total_slip,mean_slip,ShearHeat,maxHeat,meanStress] = FailureDepthRewrite2(tides_s,tides_n)
mu=0.1; %friction
g=1.315; %gravity
rho=780; %ice shell density
G=3.3e9; %shear modulus

%What if you use damaged modulus
G=3e8;

nz=100; %Depth nodes
dz=100; %element size
nt=100; %Time nodes

%Normal stress a function of depth and time
sigma_n=zeros(nz,nt);
sigma_n(:,1)=dz.*((1:nz))*g*rho;
for i=2:nz
    sigma_n(i,:)=sigma_n(i,1);
end

calc_local_stress=1; %Do you want to calculate local tidal stress or appproximate

if calc_local_stress
%Describe resolved on fault based on Tidal Stress global
i=139; %latitude + 90
j=232; %longitude (+West)
beta=45*(pi/180); %Degrees clockwise from north
clear Sxx3 Syy3 Sxy3 tides_s tides_n
tides_s=0; tides_n=0;

Sxx3(1,:)=Sxx(i,j,:);
Syy3(1,:)=Syy(i,j,:);
Sxy3(1,:)=Sxy(i,j,:);

tides_n=2*(0.5*(Syy3-Sxx3).*sin(2*beta)+Sxy3.*cos(2*beta));
tides_s=Sxx3.*sin(beta).^2+Syy3.*cos(beta).^2-Sxy3.*sin(2*beta);
else
    tides0=1e5;
    tides_n=tides0*sin(2*pi*(1:100)./nt);
    tides_s=tides0*cos(2*pi*(1:100)./nt);
end
%Add background stress
%tides_n=tides_n-1e6;

for j=1:nt
    sigma_n(:,j)=sigma_n(:,j)+tides_n(1,j);
end



%Calculate slip direction
slip0_dir=zeros(nz,nt);

fault_depth=zeros(1,nt);

%Slip distance
delta=zeros(nz,nt);

%Try calculating fault displacement assuming mode 3 frature solution
delta_mode3=delta;

%delta = (2*S/G)*sqrt(d^2 -x^2); 
%Where is this equation from?
%Here x is distance to fault tip, d is fault length
% S is the excess shear stress
% Excess shear stress changes with depth, so for simplicity, use excess
% shear stress calculated at surface
maxdepth=0;
for j=1:nt
    slip0=abs(tides_s(1,j))-mu.*sigma_n(:,j);
    
    %Fault not slipping
    slip0(slip0<0)=0;
    %Falt slip direction
    slip0(slip0~=0)=abs(tides_s(j))./tides_s(j);
    
    slip0_dir(:,j)=slip0;
    
    %What if actual slip direction can be calculated
    %slip_depth=1e3+(abs(tides_s(j))-tides_n(j)*mu)/(mu*rho*g);
    slip_depth=(abs(tides_s(j))-tides_n(j)*mu)/(mu*rho*g);
    slip_depth(slip_depth<0)=0;
    fault_depth(j)=slip_depth;
    
    %Funny way of calculating slip magnitude
    slip0=abs(tides_s(1,j))-mu.*sigma_n(:,j);
    %Fault not slipping
    slip0(slip0<0)=0;
    delta(:,j)=(abs(tides_s(j))/tides_s(j))*slip_depth*(2*slip0./G);
      
    %Mode III fracture calculations
    d=fault_depth(j);
    x=(1:nz)*dz;
    S=tides_s(j);
    
    delta_mode3(:,j)=(2*S/G).*sqrt(d.^2-x.^2);
    delta_mode3((d-x)<0,j)=0;
    
    if d>maxdepth
        maxdepth=d;
    end
end

maxdepth;
total_slip=2*(max(max(delta_mode3))-min(min(delta_mode3)));
total_slip
mean_slip=mean(delta_mode3(1,:));

slipWdepth=(max(delta_mode3')-min(delta_mode3'))*2;
%ShearHeat=mu*rho*g*(1:2000).*slipWdepth./(3.55*24*3600);
ShearHeat=mu*(mean(sigma_n')).*slipWdepth./(3.55*24*3600);
ShearHeat=mu*(mean(sigma_n')).*slipWdepth./(33*3600); %Enceladus

%New way to integrate shear heating
[ddelta_dt,ddelta_dz]=gradient(delta_mode3);
ShearHeat2=zeros(1,nz);
for i=1:nz
    ShearHeat2(i)=mu*sum(abs(ddelta_dt(i,:)).*sigma_n(i,:)./(dz));
end
ShearHeat2(ShearHeat2<0)=0;
ShearHeat2=ShearHeat2.*(nt/(3.55*24*3600));
ShearHeat2=ShearHeat2.*(nt/(33*3600)); %Enceladus

FricT=100+(ShearHeat+0.05).*(1:nz)./(3);

maxHeat=max(ShearHeat);
ShearHeat(ShearHeat<0)=0;

meanStress=mean(tides_s(tides_n>0));


%end
% 
depth=linspace(0,dz*nz/1000,nz);
figure

plot(ShearHeat,depth,'b-','Linewidth',2)
hold on
plot(ShearHeat2,depth,'k-','Linewidth',2)
set(gca,'Fontsize',24,'Ydir','reverse','XaxisLocation','top')
xlabel('Frictional Heating (W/m^2)')
ylabel('Depth (km)')





