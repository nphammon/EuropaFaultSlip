%function [Fs] = simpleThermal(ShearHeat)

D=dz*nz;
dz=dz;
dx=dz;
nz=round(D/dx);
nx=nz;

T=zeros(nz,nz);
Ts=100;
Tm=270; 
T(:,1)=Ts;
T(:,nz)=Tm;

for i=1:nx
    T(i,:)=linspace(Ts,Tm,nx);
end

k=3;
rho=920;
Cp=2100;
kappa=zeros(1,nz);
%Assume upper few km reduced conductivity
z_layer=1e3; 
kappa0=k./(rho*Cp);
kappa(1,:)=kappa0.*exp(-z_layer./((1:nz).*dz));
kappa(kappa<kappa0/10)=kappa0/10;

%Or have constant conductivity
%kappa(1,:)=kappa0;

dt=0.2*dx^2/max(kappa);

%Heating
Q=zeros(nz,nz);
%Grab shear heating solution from figure 5
Q(2,:)=ShearHeat(:)/dx/2;
%Q(2,1:40)=8e-4;

%Q=Q*5;

for kk=1:3e3
%Evolve temperature
Tnew=T;
for i=2:nx-1
            dTsqrd=(T(i,3:nx)+T(i,1:nx-2)+T(i-1,2:nx-1)+T(i+1,2:nx-1)-4*T(i,2:nx-1));%second gradient
            dTdt=((Q(i,2:nx-1)/(rho*Cp))+dTsqrd.*kappa(2:end-1)./(dx^2));
            
            %Hold on, we need kappa inside the derivative
            k_u=(kappa(2:end-1)+kappa(1:end-2))./2;
            k_d=(kappa(3:end)+kappa(2:end-1))./2;
            k_c=kappa(2:end-1); %Kappa in center element
                        
            dTsqrd=(T(i,3:nx).*k_d+T(i,1:nx-2).*k_u+T(i-1,2:nx-1).*k_c+T(i+1,2:nx-1).*k_c-2*T(i,2:nx-1).*k_c-T(i,2:nx-1).*k_d-T(i,2:nx-1).*k_u);%second gradient
            dTdt=((Q(i,2:nx-1)/(rho*Cp))+dTsqrd./(dx^2));
            
            Tnew(i,2:nx-1)=T(i,2:nx-1)+dt.*dTdt;        
end
Tnew(1,:)=Tnew(2,:);
Tnew(nx,:)=Tnew(nx-1,:);
T=Tnew;
T(T>Tm)=Tm;

kk_s=mean([kappa(1,1),kappa(1,2)]);


Fs=rho*Cp*kk_s*(T(1,2)-Ts)/(dx);
Fs_time=[kk,round(Fs*1000)]
end

depth=linspace(10,0,nx);
depth1=linspace(0,10,nx);
figure;
subplot(1,2,2)
contourf(depth1,depth,rot90(T))
set(gca,'Fontsize',24,'Ydir','reverse')
xlabel('Distance (km)')
ylabel('Depth (km)')
colorbar
subplot(1,2,1)
plot(kappa*rho*Cp,depth1,'k-','Linewidth',3)
set(gca,'Fontsize',24,'Ydir','reverse')
xlabel('Thermal Conductivity (W/mK)')
ylabel('Depth (km)')
%end


