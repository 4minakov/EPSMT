clear all
close all
%
addpath .\data
load(['model_parameters_MC2.mat']);
load(['Y0_MC2']);
load Data_total_MC2_new.mat
load('lapaz')
load('davos')
load('batlow')
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
% Maximum shear strains
shstr_tot = sqrt((ex_tot - ey_tot ).^2/4 + exy_tot.^2);
shstr_el  = sqrt((ex_el  - ey_el  ).^2/4 + exy_el.^2);
shstr_pl  = sqrt((ex_pl  - ey_pl  ).^2/4 + exy_pl.^2);
% Bulk strain
bstr_tot = ex_tot + ey_tot;
bstr_el  = ex_el  + ey_el;
bstr_pl  = ex_pl  + ey_pl;
%
Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
[maxM0,indmaxM0] = max(M0(:));
%%
figure,  
subplot(221),
pcolor(xp/R0,yp/R0,shstr_el), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
shading interp, 
tt = title(['(a)    Shear Strain (Elastic)']); 
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
colorbar('Position',[0.47 0.58 0.02 0.18])
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.5,'Color',[0.4 0.4 0.4])
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, ylabel('y/R_0'), 
axis([-1 1 -1 1]*2.5), 
subplot(222),
pcolor(xp/R0,yp/R0,bstr_el), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
shading interp, 
tt = title(['(b)    Dilation (Elastic)']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
colorbar('Position',[0.91 0.58 0.02 0.17])
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.5,'Color',[0.4 0.4 0.4])
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, 
axis([-1 1 -1 1]*2.5), 
colormap(flipud(lapaz))
subplot(223),
pcolor(xp/R0,yp/R0,shstr_pl), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
shading interp, 
colorbar('Position',[0.48 0.11 0.02 0.17])
tt=title(['(c)  Shear Strain (Plastic)']);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.5,'Color',[0.4 0.4 0.4])
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R_0'),ylabel('y/R_0')
axis([-1 1 -1 1]*2.5), %caxis([0 10]*1e-3)
subplot(224),
pcolor(xp/R0,yp/R0,bstr_pl), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
shading interp, 
%
colorbar('Position',[0.92 0.11 0.02 0.17])
%
tt=title(['(d)  Dilation (Plastic)']);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.5,'Color',[0.4 0.4 0.4])
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R_0'),
axis([-1 1 -1 1]*2.5), 
colormap(flipud(lapaz))
%
set(gcf,'units','normalized','position',[0.1000 0.1000 0.3 0.45])
print('fig\Fig3','-dpng','-r600')
print('fig\Fig3','-painters','-depsc','-r600')