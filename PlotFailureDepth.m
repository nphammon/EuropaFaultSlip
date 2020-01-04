%Make Plot of Shear Stress and Failure Depth
%mu=0.6;
figure
subplot(2,1,2)
plot(linspace(0,2*pi,100),fault_depth(1:100)./1e3,'k-','Linewidth',3);
hold on
ar=area(linspace(0,2*pi,100),fault_depth(1:100)./1e3,'FaceColor','k');
ar.FaceAlpha = 0.1;
set(gca,'Fontsize',24,'Ydir','reverse')
ylabel('Depth (km)')
%xlabel('Tidal Cycle')
set(gca,'XTick',0:pi/2:2*pi,'XAxisLocation','bottom') 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
axis([0 2*pi 0 1.5])

 
sigma_n(:,1)=sigma_n(:,100);
subplot(2,1,1)
hold on
plot(linspace(0,2*pi,100),abs(tides_s(1:100)./1e3),'k-','Linewidth',4)
set(gca,'Fontsize',20,'XAxisLocation','bottom')
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('Shear Stress (kPa)')
xlabel('Tidal Cycle')
axis([0 2*pi -10 150])
plot(linspace(0,2*pi,100),(mu.*sigma_n(1,1:100)./1e3),'k--','Linewidth',1.5)
plot(linspace(0,2*pi,100),(mu.*sigma_n(100,1:100)./1e3),'k-.','Linewidth',1.5)
%plot(linspace(0,2*pi,100),(mu*sigma_n(100,1:100)./1e3),'k-.','Linewidth',1.5)
%plot(linspace(0,2*pi,100),(sigma_n(40,1:100)./1e3).*mu,'k-.','Linewidth',1.5)

depth=linspace(0,2,100);
figure
plot(-delta_mode3(:,1),depth,'k-','Linewidth',4)
hold on
plot(-delta_mode3(:,13),depth,'k-','Linewidth',3)
plot(-delta_mode3(:,25),depth,'k-','Linewidth',2)
plot(-delta_mode3(:,37),depth,'k--','Linewidth',2)
plot(-delta_mode3(:,50),depth,'k-.','Linewidth',2)
% plot(-delta_mode3(:,62),depth,'b-.','Linewidth',2)
% plot(-delta_mode3(:,75),depth,'b--','Linewidth',2)
% plot(-delta_mode3(:,87),depth,'b-','Linewidth',2)
axis([-0.25 0.25 0 1])
set(gca,'Fontsize',24,'Ydir','reverse','XAxisLocation','top')
xlabel('\delta Fault Displacement (m)')
ylabel('Depth (km)')

%Random table
mu=[0.1,0.12,0.15,0.2,0.3,0.4,0.5,0.6];
%WARNING, these two arrays (maxdepth and total slip), need to be recomputed
%by running this code for different coefficients of friction.
maxdepth=[0.433,0.36805,0.305,0.2423,0.1812,0.1518,0.1394,0.124];
totalslip=[0.4308,0.36,0.29,0.2202,0.1523,0.12,0.1015,0.0898];
totalslip=[0.0799,0.0667,0.0535,0.0405,0.0277,0.0215,0.018,0.0159];
maxdepth=totalslip*8;

figure
subplot(2,1,2), scatter(mu,maxdepth,150,'ks','Filled')
muline=linspace(0.1,0.6,50);
p=polyfit(mu,maxdepth,5);
y1=polyval(p,muline);
hold on
plot(muline,y1,'k--','Linewidth',3)
xlabel('\mu, Coefficient of Friction')
ylabel('Depth or Slip Distance(km)')
yyaxis right
scatter(mu,totalslip,150,'bs','Filled')
p=polyfit(mu,totalslip,5);
y2=polyval(p,muline);
plot(muline,y2,'b--','Linewidth',3)
ylabel('Total Slip per Cycle (km)')
set(gca,'Fontsize',20)

subplot(2,1,1)
plot(linspace(0,2*pi,100),tides_s(1:100)./1e3,'k-','Linewidth',4)
hold on
plot(linspace(0,2*pi,100),tides_n(1:100)./1e3,'k--','Linewidth',3)
ylabel('Tidal Stress (kPa)')
set(gca,'Fontsize',20)
xlabel('Tidal Cycle')





