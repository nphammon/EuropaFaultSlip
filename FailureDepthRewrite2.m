function [total_slip,mean_slip,ShearHeat,maxHeat,meanStress] = FailureDepthRewrite2(tides_s,tides_n)

mu=0.1; %friction
g=1.34; %gravity
rho=782; %ice shell density
G=3e9; %shear modulus

%What if you use damaged modulus
G=3e8;

nz=100; %Depth nodes
dz=20; %element size
nt=100; %Time nodes

%Normal stress a function of depth and time
sigma_n=zeros(nz,nt);
sigma_n(:,1)=dz.*((1:nz))*g*rho;
for i=2:nz
    sigma_n(i,:)=sigma_n(i,1);
end



%Add background stress
%tides_n=tides_n-5e5;

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
total_slip;
mean_slip=mean(delta_mode3(1,:));

slipWdepth=(max(delta_mode3')-min(delta_mode3'))*2;
%ShearHeat=mu*rho*g*(1:2000).*slipWdepth./(3.55*24*3600);
ShearHeat=mu*(mean(sigma_n')).*slipWdepth./(3.55*24*3600);
FricT=100+(ShearHeat+0.05).*(1:nz)./(3);

maxHeat=max(ShearHeat);

meanStress=mean(tides_s(tides_n>0));


end
% 
% depth=linspace(0,2,nz);
% figure
% hold on
% plot(ShearHeat,depth,'b-','Linewidth',2)
% set(gca,'Fontsize',24,'Ydir','reverse','XaxisLocation','top')
% xlabel('Frictional Heating (W/m^2)')
% ylabel('Depth (km)')





