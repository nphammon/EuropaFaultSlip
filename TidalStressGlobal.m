%New Program to calculate resolve stress as a function of tides
clear all
h2=1.5; M=1.9*10^27; G=3.3*10^9; v=0.33; roe=3015; a=6.7*10^8; 

A=3*h2*M*G*(1+v)/(8*pi*roe*a^3*(5+v));
%A=1e6; %Tidal stress magnitude
e=0.009;

%Add obliquity
epsilon=0*(pi/180);
w=0; %argument of periapsis

Sxx=zeros(180,360,100);
Syy=zeros(180,360,100);
Sxy=zeros(180,360,100);
Sxx2=Sxx;
Syy2=Syy;
Sxy2=Sxy;
S1=Sxx;
S2=S1;
S3=S1; 

S_mag=Sxx;

alpha=0;
theta=0;
lat=0;
lon=0;
t=0;
Stt=0;
Spp=0;
figure
for k=1:100 %time
for i=1:180  %lat
for j=1:360 %lon

t=2*pi*k/100; %time
lat=(i-90)*pi/180; %lat
lon=j*pi/180; %lon

%Fllow methods of Hurford et al., 2009, Appendix

%First calculate permanent bulge stresses
theta=acos(cos(lat)*cos(lon));
    
Stt=A*(5+3*cos(2*theta)); %simga theta theta (in direction to bulde axis
Spp=-A*(1-9*cos(2*theta)); %sigma phi phi, perpendicular

 alpha=acos(-(sin(lat)*cos(theta)/(cos(lat)*sin(theta)))); %Algle to tidal bulge axis
 alpha=real(alpha);
            %use old code lines to trouble shoot
 %           alpha=acos((-sin(lat)*cos(theta))/(cos(lat)*sin(theta)));
 %           alpha=real(alpha);

 Sxx(i,j,k)=0.5*(Spp+Stt)-0.5*(Spp-Stt)*cos(2*alpha);
 Syy(i,j,k)=0.5*(Spp+Stt)+0.5*(Spp-Stt)*cos(2*alpha);
 Sxy(i,j,k)=-0.5*(Spp-Stt)*sin(2*alpha);
 
 %Calculate stresses from time dependent bulge 
 theta=acos(sin(lat)*sin(epsilon*sin(t-w))+cos(lat)*cos(epsilon*sin(t-w))*cos(lon+(2*e*sin(t)))); %angular distance to tidal bulge axis=
 %theta=acos(cos(lat)*cos(lon));
 
 %use old code line to trouble shoot
 %theta=acos(cos(lat)*cos(lon+2*e*sin(t)));
 
 Stt2=A*(5+3*cos(2*theta))./((1-e*cos(t))^3); %simga theta theta (in direction to bulde axis
 Spp2=-A*(1-9*cos(2*theta))./((1-e*cos(t))^3); %sigma phi phi, perpendicular
 
 alpha=acos((sin(epsilon*sin(t-w))-sin(lat)*cos(theta))/(cos(lat)*sin(theta))); %Algle to tidal bulge axis
 alpha=real(alpha);
 
              %old code
            % alpha=acos((-sin(lat)*cos(theta))/(cos(lat)*sin(theta)));
            %alpha=real(alpha);
 
 Sxx2(i,j,k)=0.5*(Spp2+Stt2)-0.5*(Spp2-Stt2)*cos(2*alpha);
 Syy2(i,j,k)=0.5*(Spp2+Stt2)+0.5*(Spp2-Stt2)*cos(2*alpha);
 Sxy2(i,j,k)=-0.5*(Spp2-Stt2)*sin(2*alpha);
 
 %Replace Sxx1 with real answer
 Sxx(i,j,k)=Sxx2(i,j,k)-Sxx(i,j,k);
 Syy(i,j,k)=Syy2(i,j,k)-Syy(i,j,k);
 Sxy(i,j,k)=Sxy2(i,j,k)-Sxy(i,j,k);
 
 %Determine principal stresses
 gamma=0.5*atan(2*Sxy(i,j,k)/(Sxx(i,j,k)-Syy(i,j,k)));
 gamma=real(gamma);
 S1(i,j,k)=Sxx(i,j,k)*(cos(gamma)^2)+Syy(i,j,k)*(sin(gamma)^2) ...
    +Sxy(i,j,k)*sin(2*gamma);
 S2(i,j,k)=Sxx(i,j,k)*(sin(gamma))^2+Syy(i,j,k)*(cos(gamma))^2 ...
    -Sxy(i,j,k)*sin(2*gamma);

S1yo=max(S1(i,j,k),S2(i,j,k));
S2yo=min(S1(i,j,k),S2(i,j,k));
S1(i,j,k)=S1yo;
S2(i,j,k)=S2yo;
 

S3(i,j,k)=(S1(i,j,k)-S2(i,j,k))/2;



 S_mag(i,j,k)=(1/2)*(Sxx(i,j,k)+Syy(i,j,k));
 

end
end

stress_mag_small=S_mag(1:10:180,1:10:360,k);
[dx,dy]=gradient(stress_mag_small);


hold off
contourf(S1(:,:,k))
set(gca,'Fontsize',24)
ylabel('Lattitute')
xlabel('Longitude')    
colorbar
caxis([-1.2e5,1.2*10^5])
hold on
quiver(1:10:360,1:10:180,-dx,-dy,'k')


F(:,k)=getframe(gcf);

k
end


% %Calculate resolve stress on fault 
% beta=(1:180).*(pi/180); %Fault
% tau_fault=zeros(180,100);
% sigma_fault=zeros(180,100);
% 
% %for random location, how does stress vary on fault for a given orientation
% i=1+round(rand*179);
% j=1+round(rand*179);
% 
% i=139;
% j=232;
% i=20;
% Sxx3(1,:)=Sxx(i,j,:);
% Syy3(1,:)=Syy(i,j,:);
% Sxy3(1,:)=Sxy(i,j,:);
% 
% % for k2=1:180
% % tau_fault(k2,:)=2*(0.5*(Syy3-Sxx3).*sin(2*beta(k2))+Sxy3.*cos(2*beta(k2)));
% % sigma_fault(k2,:)=Sxx3.*sin(beta(k2)).^2+Syy3.*cos(beta(k2)).^2-Sxy3.*sin(2*beta(k2));
% % end
% 
% %Use Rhoden obliquity paper equations, 3a 3b
% %Find phase lag as function of fault orientation
% phase_lag=zeros(1,180);
% 
% for k2=1:180
%     beta2=beta(k2);
%     tau_fault(k2,:)=(-0.5*(Syy3-Sxx3).*sin(2*beta2)+Sxy3.*cos(2*beta2));
%     sigma_fault(k2,:)=0.5.*(Sxx3+Syy3)+0.5.*(Syy3-Sxx3).*cos(2*beta2)+Sxy3.*sin(2*beta2);
%     [a1,a2]=max(tau_fault(k2,:));
%     [b1,b2]=max(sigma_fault(k2,:));
%     phase_lag(k2)=(b2-a2)*3.6;
% end
% 
% figure
% plot(phase_lag);
% 
% TA=linspace(0,2*pi,100); %True Anamoly
% p=[1,30,60,90,120,150];
% 
% figure;
% %subplot(2,1,1)
% plot(TA,tau_fault(110,:)./1e3,'k-','Linewidth',3);
% hold on
% plot(TA,sigma_fault(110,:)./1e3,'k--','Linewidth',3);
% set(gca,'Fontsize',18)
% xlabel('Orbital Cycle')
% ylabel('Stress (kPa)')
% axis([0 2*pi -100 100])
% 
% 
% 
% figure
% for b=1:6
%     subplot(2,3,b)
%     plot(TA,tau_fault(p(b),:)./1e3,'k-','Linewidth',2);
%     hold on
%     plot(TA,sigma_fault(p(b),:)./1e3,'k--','Linewidth',2);
%     set(gca,'Fontsize',18)
%     xlabel('Orbital Cycle')
%     ylabel('Stress (kPa)')
%     axis([0 2*pi -200 200])
% end
% 
% maxS=0;
% maxN=0;
% max1=0;
% for b=1:180
%     max1=max(tau_fault(b,:));
%     max2=max(sigma_fault(b,:));
%     maxN=max(max2,maxN);
%     if max1>maxS
%         maxS=max1;
%         maxb=b;
%     end
% end
% maxb

% figure
% for k=1:100
% hold off
% contourf(S3(:,:,k))
% set(gca,'Fontsize',24)
% ylabel('Lattitute')
% xlabel('Longitude')    
% colorbar
% caxis([0,1.5*10^5])
% hold on
% quiver(1:10:360,1:10:180,-dx,-dy,'k')
% 
% F(:,k)=getframe(gcf);
% end