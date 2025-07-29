clear all
close all
addpath .\data
load(['model_parameters_ANI.mat']);
load(['Y0_ANI.mat']);
load Data_total_ANI_MC.mat
load(['res150_ANI']);
load davos
%
D1 = inv(D);%Compliance matrix
% Elastic strains
ex_el  =  D1(1,1)*Sx_tot + D1(1,2)*Sy_tot + D1(1,3)*Sxy_tot;
ey_el  =  D1(2,1)*Sx_tot + D1(2,2)*Sy_tot + D1(2,3)*Sxy_tot;
exy_el = (D1(3,1)*Sx_tot + D1(3,2)*Sy_tot + D1(3,3)*Sxy_tot)/2;
% Plastic strains
ex_pl = ex_tot  - ex_el;
ey_pl = ey_tot  - ey_el;
exy_pl= exy_tot - exy_el;
%
Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);
%   
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
indmaxM0 = find(M0==max(M0(yp<-0.1)));
Mw = 2/3*log10(M0)-6;
Mw(Mw<max(Mw(:))-3)=NaN;
M0(Mw<max(Mw(:))-3)=NaN;
%
figure,
pcolor(xp/R0,yp/R0,1e-5*M0/R0^2), hold on
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'k+','MarkerSize',15)
plot(xp(indmaxM0)/R0,yp(indmaxM0)/R0,'ko','MarkerSize',10)
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.1,'Color',[0.4 0.4 0.4])
shading interp,  
colormap(flipud(davos)); 
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R0'),ylabel('y/R0')
caxis([0 .75]), 
tt = title(['Seismic moment  $\big(\frac{M_{0}}{R_0^2}\big)\cdot 10^{-5}$'],...
    'Interpreter','latex'); 
c=colorbar; 
xlabel('x/R0'),ylabel('y/R0'),
axis([-1 1 -1 1]*2.5), 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4]),
%
print('fig/Fig4b','-dpng','-r600'), 
print('fig/Fig4b','-painters','-depsc','-r600'), 