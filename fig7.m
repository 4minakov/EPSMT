clear all
close all
addpath .\data
load('davos')
load('lapaz')
load('batlow')
load(['model_parameters_ANI.mat']);
load(['Y0_ANI']);
load Data_total_ANI_MC.mat
load(['res150_ANI']);
%
D1 = inv(D);%Compliance matrixclear
% Elastic strains
ex_el  =  D1(1,1)*Sx_tot + D1(1,2)*Sy_tot + D1(1,3)*Sxy_tot;
ey_el  =  D1(2,1)*Sx_tot + D1(2,2)*Sy_tot + D1(2,3)*Sxy_tot;
exy_el = (D1(3,1)*Sx_tot + D1(3,2)*Sy_tot + D1(3,3)*Sxy_tot)/2;
% Plastic strains
ex_pl = ex_tot  - ex_el;
ey_pl = ey_tot  - ey_el;
exy_pl= exy_tot - exy_el;
%%
Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);  
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
TAU = zeros(length(M0(:)),1);
K = zeros(length(M0(:)),1);
for i = 1:length(M0(:))
    MT = [Mxx(i), Mxy(i), 0; Mxy(i), Myy(i), 0; 0 0 Mzz(i)];
    
    [TAU(i),K(i)] = MT2tauk(MT);
    
end
TAU = reshape(TAU,size(M0));
K = reshape(K,size(M0));
indmaxM0 = find(M0==max(M0(yp<-0.1)));
Mw = 2/3*log10(M0)-6;
%%
figure(11),clf
ax1 = subplot(121);
tau1 = TAU; tau1(Mw<max(Mw(:))-3)=NaN;
pcolor(xp/R0,yp/R0,tau1), hold on
        plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
        plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
shading interp, colorbar('eastoutside')
tt=title(['(a)              \tau'],'FontSize',12);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.2,'Color',[0.4 0.4 0.4])
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R0'),ylabel('y/R0')
axis([-1 1 -1 1]*2.5), %colormap(ax1,lapaz)
ax2 = subplot(122);
K1 = K; K1(Mw<max(Mw(:))-3)=NaN;
pcolor(xp/R0,yp/R0,K1), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
shading interp, colorbar('eastoutside')
tt=title(['(b)             \kappa'],'FontSize',12);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.2,'Color',[0.4 0.4 0.4])
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R0'),ylabel('y/R0')
axis([-1 1 -1 1]*2.5),
colormap(batlow)
%%
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 4]),
print(['fig\Fig7'],'-dpng','-r600')
print(['fig\Fig7'],'-painters','-depsc','-r600')
