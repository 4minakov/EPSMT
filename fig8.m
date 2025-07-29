clear all
close all
%
addpath .\data
load(['model_parameters_WET.mat']);
load(['Y0_WET']);
load(['res30_WET']);
load('davos')
load('tokyo')
load Data_total_Tens2.mat
%%
D1 = inv(D);%Compliance matrixclear
% Elastic strains
ex_el  =  D1(1,1)*Sx_tot + D1(1,2)*Sy_tot + D1(1,3)*Sxy_tot;
ey_el  =  D1(2,1)*Sx_tot + D1(2,2)*Sy_tot + D1(2,3)*Sxy_tot;
exy_el = (D1(3,1)*Sx_tot + D1(3,2)*Sy_tot + D1(3,3)*Sxy_tot)/2;
% Plastic strains
ex_pl = ex_tot  - ex_el;
ey_pl = ey_tot  - ey_el;
exy_pl= exy_tot - exy_el;

Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);
    
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
[maxM0,indmaxM0] = max(M0(:));

Mw = 2/3*log10(M0)-6;
M0(Mw<max(Mw(:))-1)=NaN;
%
figure,
subplot(121)
pcolor(xp/R0,yp/R0,1e-6*M0/R0^2), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.3,'Color',[0.4 0.4 0.4])
%Mw = 2/3*log10(sum(M0(:))) - 6;
shading interp, %title(['M_w']), 
colormap(flipud(davos)); 
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R0'),ylabel('y/R0')
caxis([0 0.03]), 
c=colorbar;%c.Label.String='$\big(\frac{M_{0}}{R_0^2}\big)\cdot 10^{-6}$';
c.Label.FontWeight='bold'; c.Label.FontSize=14;
c.Label.Interpreter='latex';c.Location='southoutside';
xlabel('x/R0'),ylabel('y/R0'),
axis([-1 1 -1 1]*2.5), 
tt = title(['(a) Seismic moment  $\big(\frac{M_{0}}{R_0^2}\big)\cdot 10^{-6}$'],...
    'Interpreter','latex'); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
subplot(122)
contourf(xp/R0,yp/R0,pf2d*1e-6), hold on
plot(xp(F_mode==1 & ~isnan(M0) )/R0,yp(F_mode==1 & ~isnan(M0) )/R0,'.','Color',[1 0.7 .7]), hold on
plot(xp(F_mode==-1 & ~isnan(M0) )/R0,yp(F_mode==-1 & ~isnan(M0) )/R0,'.y'),
shading interp, 
colormap(flipud(davos)); 
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R0'),ylabel('y/R0')
c=colorbar; c.Label.String='MPa'; c.Label.Position = [6 0 0]; c.Label.Rotation=0;
c.Location='southoutside';
c.Label.FontWeight='normal'; c.Label.FontSize=10;
text(1.1,0,'I','FontSize',12,'FontWeight','bold')
text(1.5,1.5,'II','FontSize',12,'FontWeight','bold')
xlabel('x/R0'),ylabel('y/R0'),
axis([-1 1 -1 1]*2.5), 
tt = title(['(b)  Fluid pressure'],'Interpreter','latex'); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 4]),
print('fig\Fig8','-dpng','-r600'), 
print('fig\Fig8','-painters','-depsc','-r600'), 