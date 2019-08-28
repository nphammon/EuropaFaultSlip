%Plot max Shear Heating
maxmaxHeat=0;
maxmaxHeatmap=zeros(180,360);
maxbMap=zeros(180,360);
for i=1:10:180
    for j=1:10:360
        
        Sxx3(1,:)=Sxx(i,j,:);
        Syy3(1,:)=Syy(i,j,:);
        Sxy3(1,:)=Sxy(i,j,:);
        maxb=0;
        maxmaxHeat=0;
        for b=1:5:180
            beta=b*pi/180;
            tides_s=0.5*(Syy3-Sxx3).*sin(2*beta)+Sxy3.*cos(2*beta);
            tides_n=Sxx3.*cos(beta).^2+Syy3.*sin(beta).^2+Sxy3.*sin(2*beta);
        if j>180
            tides_s=-tides_s;
        end
        
            [totalslip,mean_slip,ShearHeat,maxHeat,meanStress]=FailureDepthRewrite2(tides_s,tides_n);

            if maxHeat>maxmaxHeat
                maxmaxHeat=maxHeat;
                maxb=b;
            end
        end
        maxmaxHeatmap(i,j)=maxmaxHeat;
        maxbMap(i,j)=maxb;
        yo=[i,j]
    end
end


figure
%contourf(1:10:360,1:10:180,maxmaxHeatmap(1:10:180,1:10:360))
maxmaxHeatmap(180,1:10:360)=maxmaxHeatmap(1,1:10:360);
contourf(1:10:360,[1:10:171,180],maxmaxHeatmap([1:10:171,180],1:10:360))
for i=1:10:180
    for j=1:10:360
        b=maxbMap(i,j)*pi/180;
        if j<180
            b=-b;
        end
        j1=j+3*sin(b); j2=j-3*sin(b); 
        j3=j+3*sin(b+pi/2); j4=j-3*sin(b+pi/2);
        
        i1=i+5*cos(b); i2=i-5*cos(b); 
        i3=i+5*cos(b+pi/2); i4=i-5*cos(b+pi/2);
        
        hold on
        plot([j1,j2],[i1,i2],'k-');
        plot([j3,j4],[i3,i4],'k-');
    end
    i
end
        

i=91;
j=106;
Sxx3(1,:)=Sxx(i,j,:);
Syy3(1,:)=Syy(i,j,:);
Sxy3(1,:)=Sxy(i,j,:);

heat_b=zeros(1,180);
for b=1:180
    beta=b*pi/180;
    tides_s=0.5*(Syy3-Sxx3).*sin(2*beta)+Sxy3.*cos(2*beta);
    tides_n=Sxx3.*cos(beta).^2+Syy3.*sin(beta).^2+Sxy3.*sin(2*beta);
    [totalslip,mean_slip,ShearHeat,maxHeat,meanStress]=FailureDepthRewrite2(tides_s,tides_n);
    heat_b(b)=maxHeat;
    %heat_b(b)=totalslip;
    b
end
figure
plot(heat_b)
% 
% i=160;
% j=120;
% Sxx3(1,:)=Sxx(i,j,:);
% Syy3(1,:)=Syy(i,j,:);
% Sxy3(1,:)=Sxy(i,j,:);
% s_max=zeros(1,180);
% sn_max=s_max;
% for b=1:180
%         beta=b*pi/180;
%         tides_s=0.5*(Syy3-Sxx3).*sin(2*beta)+Sxy3.*cos(2*beta);
%         tides_n=Sxx3.*cos(beta).^2+Syy3.*sin(beta).^2+Sxy3.*sin(2*beta);
%         s_max(b)=max(tides_s);
%         sn_max(b)=max(tides_n);
% end
% figure
% plot(s_max,'k-')
% hold on
% plot(sn_max,'k--')
%         

    