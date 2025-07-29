%% Hudson plots
clear all
close all
%
addpath .\data
load Data_total_SH_MC2.mat
load('model_parameters_MC2.mat')
load('davos')
%
hudson_analyt
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
%
Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);  
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
Mw = 2/3*log10(M0)-6;
%find(Mw(Mw<-6);
TAU = zeros(length(Mw(:)),1);
K = zeros(length(M0(:)),1);
for i = 1:length(M0(:))
    MT = [Mxx(i), Mxy(i), 0; Mxy(i), Myy(i), 0; 0 0 Mzz(i)];
    
    [TAU(i),K(i)] = MT2tauk(MT);
    
end
TAU = reshape(TAU,size(M0));
K = reshape(K,size(M0));
[maxM0,inmaxM0]=max(M0(:));
M1 = log10(M0(:)/maxM0);
indmaxM0 = find(M0==max(M0(yp<-0.1)));
tau=TAU(indmaxM0);
k=K(indmaxM0);
fm = zeros(3,3); %maximum value
fm(1,1) = Mxx(indmaxM0); fm(2,2) = Myy(indmaxM0);
fm(1,2) = Mxy(indmaxM0); fm(2,1) = Mxy(indmaxM0);
fm(3,3) = Mzz(indmaxM0);
fm_ISO =fm;
%save fm_ISO fm_ISO
%
figure
% title(['P=',num2str(fix(P_out)*1e-6),'MPa, \tau = ', ...
%     num2str(fix(tau_out)*1e-6),'MPa'] )
subplot(221)
tt = title(['(a) ISO']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on,
%caxis([-3 0])
plot([-1:.1:0],1+[-1:.1:0],'k','LineWidth',2)
plot([0:.1:1],1-[0:.1:1],'k','LineWidth',2)
plot([0:.1:1],-1+[0:.1:1],'k','LineWidth',2)
plot([-1:.1:0],-1-[-1:.1:0],'k','LineWidth',2)
plot([-1,1],[0, 0],'k','LineWidth',1)
plot([0,0],[-1, 1],'k','LineWidth',1)
text(0,1.2,'\bf{\kappa}','FontSize',20),
text(-0.3,1.1,'+V','FontSize',12),
text(1.0,0.2,'\bf{\tau} ','FontSize',20)
text(1.05,0,'CLVD','FontSize',12)
text(0,-1.05,'-V','FontSize',12),
text(-1.5,-0.2,'CLVD','FontSize',12)
text(0.05,0.1,'DC','FontSize',12)
text(-1.1,0.619,'+Crack','FontSize',12)
text( 0.45,-0.619,'-Crack','FontSize',12)
text( 0.5,0.2,'I','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,0.2,'II','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,-0.2,'III','FontSize',11,'Color',[0.4 0.4 0.4])
text( 0.5,-0.2,'IV','FontSize',11,'Color',[0.4 0.4 0.4])
plot(tau_a,k_a,'-k','LineWidth',1.5,'Color',[0.4 0.4 0.4]),
scatter(TAU(Mw>-5),K(Mw>-5),1e2*exp(M1(Mw>-5)),Mw(Mw>-5),'filled')
plot(tau,k,'+k','MarkerSize',15),
plot(tau,k,'ko','MarkerSize',10)
axis([-1 1 -1 1]),axis off, axis equal
colormap(flipud(davos));
%%
load Data_total_NU_MC.mat
load(['model_parameters_NU.mat'])
%
D0 = E*(1-NU)./(1+NU)./(1-2*NU);
D11 = D0*0+1; D12 = NU./(1 -NU); D13 = D0*0+0;
D21 = NU./(1-NU); D22 = D0*0+1; D23 = D0*0+0;
D31 = D0*0+0; D32 = D0*0+0; D33 = (1-2*NU)./2./(1-NU);

S11 = zeros(size(xp(:)));S12=S11;S13=S11;
S21 = zeros(size(xp(:)));S22=S21;S23=S21;
S31 = zeros(size(xp(:)));S32=S31;S33=S31;
% 
for i=1:length(D0(:))
    
    S1 = inv( D0(i)*[D11(i) D12(i) D13(i);...
        D21(i) D22(i) D23(i);...
        D31(i) D32(i) D33(i)]);
    S11(i) = S1(1,1);S12(i) = S1(1,2);S13(i) = S1(1,3);
    S21(i) = S1(2,1);S22(i) = S1(2,2);S23(i) = S1(2,3);
    S31(i) = S1(3,1);S32(i) = S1(3,2);S33(i) = S1(3,3);
end
% 
S11 = reshape(S11,size(xp));S12 = reshape(S12,size(xp));S13 = reshape(S13,size(xp));
S21 = reshape(S21,size(xp));S22 = reshape(S22,size(xp));S23 = reshape(S23,size(xp));
S31 = reshape(S31,size(xp));S32 = reshape(S32,size(xp));S33 = reshape(S33,size(xp));
% Elastic strains
ex_el  =  S11.*Sx_tot + S12.*Sy_tot + S13.*Sxy_tot;
ey_el  =  S21.*Sx_tot + S22.*Sy_tot + S23.*Sxy_tot;
exy_el = (S31.*Sx_tot + S32.*Sy_tot + S33.*Sxy_tot)/2;
% Plastic strains
ex_pl = ex_tot  - ex_el;
ey_pl = ey_tot  - ey_el;
exy_pl= exy_tot - exy_el;
% Moment tensor
Mxx = D0.*(D11.*ex_pl + D12.*ey_pl + D13.*exy_pl).*dA;
Myy = D0.*(D21.*ex_pl + D22.*ey_pl + D23.*exy_pl).*dA;
Mxy = D0.*(D31.*ex_pl + D32.*ey_pl + D33.*exy_pl).*dA;
Mzz = NU.*(Mxx+Myy);
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
Mw = 2/3*log10(M0)-6;
%find(Mw(Mw<-6);
TAU = zeros(length(Mw(:)),1);
K = zeros(length(M0(:)),1);
for i = 1:length(M0(:))
    MT = [Mxx(i), Mxy(i), 0; Mxy(i), Myy(i), 0; 0 0 Mzz(i)];
    [TAU(i),K(i)] = MT2tauk(MT);
end
TAU = reshape(TAU,size(M0));
K = reshape(K,size(M0));
M1 = log10(M0(:)/maxM0);
%indmaxM0 = find(M0==max(M0(yp<-0.1)));
[~,indmaxM0] = max(M0(:));
tau=TAU(indmaxM0);
k=K(indmaxM0);
%
subplot(222)
tt = title(['(b) \Delta \nu']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on,
%caxis([-3 0])
plot([-1:.1:0],1+[-1:.1:0],'k','LineWidth',2)
plot([0:.1:1],1-[0:.1:1],'k','LineWidth',2)
plot([0:.1:1],-1+[0:.1:1],'k','LineWidth',2)
plot([-1:.1:0],-1-[-1:.1:0],'k','LineWidth',2)
plot([-1,1],[0, 0],'k','LineWidth',1)
plot([0,0],[-1, 1],'k','LineWidth',1)
text(0,1.2,'\bf{\kappa}','FontSize',20),
text(-0.3,1.1,'+V','FontSize',12),
text(1.0,0.2,'\bf{\tau} ','FontSize',20)
text(1.05,0,'CLVD','FontSize',12)
text(0,-1.05,'-V','FontSize',12),
text(-1.5,-0.2,'CLVD','FontSize',12)
text(0.05,0.1,'DC','FontSize',12)
text(-1.1,0.619,'+Crack','FontSize',12)
text( 0.45,-0.619,'-Crack','FontSize',12)
text( 0.5,0.2,'I','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,0.2,'II','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,-0.2,'III','FontSize',11,'Color',[0.4 0.4 0.4])
text( 0.5,-0.2,'IV','FontSize',11,'Color',[0.4 0.4 0.4])
plot(tau_a,k_a,'-k','LineWidth',1.5,'Color',[0.4 0.4 0.4]),
scatter(TAU(Mw>-5),K(Mw>-5),1e2*exp(M1(Mw>-5)),Mw(Mw>-5),'filled')
plot(tau,k,'+k','MarkerSize',15),
plot(tau,k,'ko','MarkerSize',10)
axis([-1 1 -1 1]),axis off, axis equal
colormap(flipud(davos));

%%
load Data_total_ANI_MC.mat
load('model_parameters_ANI.mat')
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
Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);  
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
Mw = 2/3*log10(M0)-6;
%find(Mw(Mw<-6);
TAU = zeros(length(M0(:)),1);
K = zeros(length(M0(:)),1);
for i = 1:length(M0(:))
    MT = [Mxx(i), Mxy(i), 0; Mxy(i), Myy(i), 0; 0 0 Mzz(i)];
    
    [TAU(i),K(i)] = MT2tauk(MT);
    
end
TAU = reshape(TAU,size(M0));
K = reshape(K,size(M0));
indmaxM0 = find(M0==max(M0(yp<-0.1)));
M1 = log10(M0(:)/maxM0);
tau=TAU(indmaxM0);
k=K(indmaxM0);
fm = zeros(3,3); %maximum value
fm(1,1) = Mxx(indmaxM0); fm(2,2) = Myy(indmaxM0);
fm(1,2) = Mxy(indmaxM0); fm(2,1) = Mxy(indmaxM0);
fm(3,3) = Mzz(indmaxM0);
fm_ANISO =fm;
%save fm_ANISO fm_ANISO
%
subplot(223)
tt = title(['(c) ANISO']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on,
%caxis([-3 0])
plot([-1:.1:0],1+[-1:.1:0],'k','LineWidth',2)
plot([0:.1:1],1-[0:.1:1],'k','LineWidth',2)
plot([0:.1:1],-1+[0:.1:1],'k','LineWidth',2)
plot([-1:.1:0],-1-[-1:.1:0],'k','LineWidth',2)
plot([-1,1],[0, 0],'k','LineWidth',1)
plot([0,0],[-1, 1],'k','LineWidth',1)
text(0,1.2,'\bf{\kappa}','FontSize',20),
text(-0.3,1.1,'+V','FontSize',12),
text(1.0,0.2,'\bf{\tau} ','FontSize',20)
text(1.05,0,'CLVD','FontSize',12)
text(0,-1.05,'-V','FontSize',12),
text(-1.5,-0.2,'CLVD','FontSize',12)
text(0.05,0.1,'DC','FontSize',12)
text(-1.1,0.619,'+Crack','FontSize',12)
text( 0.45,-0.619,'-Crack','FontSize',12)
text( 0.5,0.2,'I','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,0.2,'II','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,-0.2,'III','FontSize',11,'Color',[0.4 0.4 0.4])
text( 0.5,-0.2,'IV','FontSize',11,'Color',[0.4 0.4 0.4])
plot(tau_a,k_a,'-k','LineWidth',1.5,'Color',[0.4 0.4 0.4]),
scatter(TAU(Mw>-5),K(Mw>-5),1e2*exp(M1(Mw>-5)),Mw(Mw>-5),'filled')
plot(tau,k,'+k','MarkerSize',15),
plot(tau,k,'ko','MarkerSize',10)
axis([-1 1 -1 1]),axis off, axis equal
colormap(flipud(davos));
%%
load Data_total_Tens2.mat
load(['model_parameters_WET.mat'])
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
Mxx = (D(1,1)*ex_pl + D(1,2)*ey_pl + D(1,3)*exy_pl).*dA;
Myy = (D(2,1)*ex_pl + D(2,2)*ey_pl + D(2,3)*exy_pl).*dA;
Mxy = (D(3,1)*ex_pl + D(3,2)*ey_pl + D(3,3)*exy_pl).*dA;
Mzz = nu*(Mxx+Myy);  
M0  = 1/sqrt(2)*sqrt(Mxx.^2 + Myy.^2 + Mzz.^2 + 2*Mxy.^2);
Mw = 2/3*log10(M0)-6;
%find(Mw(Mw<-6);
TAU = zeros(length(M0(:)),1);
K = zeros(length(M0(:)),1);
for i = 1:length(M0(:))
    MT = [Mxx(i), Mxy(i), 0; Mxy(i), Myy(i), 0; 0 0 Mzz(i)];
    
    [TAU(i),K(i)] = MT2tauk(MT);
    
end
TAU = reshape(TAU,size(M0));
K = reshape(K,size(M0));
%indmaxM0 = find(M0==max(M0(yp<-0.1)));
[maxM0,indmaxM0]=max(M0(:));
M1 = log10(M0(:)/maxM0);
tau=TAU(indmaxM0);
k=K(indmaxM0);
fm = zeros(3,3); %maximum value
fm(1,1) = Mxx(indmaxM0); fm(2,2) = Myy(indmaxM0);
fm(1,2) = Mxy(indmaxM0); fm(2,1) = Mxy(indmaxM0);
fm(3,3) = Mzz(indmaxM0);
fm_TENS =fm;
%save fm_TENS fm_TENS
%
subplot(224)
tt = title(['(d) WET']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
hold on,
%caxis([-3 0])
plot([-1:.1:0],1+[-1:.1:0],'k','LineWidth',2)
plot([0:.1:1],1-[0:.1:1],'k','LineWidth',2)
plot([0:.1:1],-1+[0:.1:1],'k','LineWidth',2)
plot([-1:.1:0],-1-[-1:.1:0],'k','LineWidth',2)
plot([-1,1],[0, 0],'k','LineWidth',1)
plot([0,0],[-1, 1],'k','LineWidth',1)
text(0,1.2,'\bf{\kappa}','FontSize',20),
text(-0.3,1.1,'+V','FontSize',12),
text(1.0,0.2,'\bf{\tau} ','FontSize',20)
text(1.05,0,'CLVD','FontSize',12)
text(0,-1.05,'-V','FontSize',12),
text(-1.5,-0.2,'CLVD','FontSize',12)
text(0.05,0.1,'DC','FontSize',12)
text(-1.1,0.619,'+Crack','FontSize',12)
text( 0.45,-0.619,'-Crack','FontSize',12)
text( 0.5,0.2,'I','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,0.2,'II','FontSize',11,'Color',[0.4 0.4 0.4])
text(-0.5,-0.2,'III','FontSize',11,'Color',[0.4 0.4 0.4])
text( 0.5,-0.2,'IV','FontSize',11,'Color',[0.4 0.4 0.4])
plot(tau_a,k_a,'-k','LineWidth',1.5,'Color',[0.4 0.4 0.4]),
scatter(TAU(Mw>-5),K(Mw>-5),1e2*exp(M1(Mw>-5)),Mw(Mw>-5),'filled')
plot(tau,k,'+k','MarkerSize',15),
plot(tau,k,'ko','MarkerSize',10)
axis([-1 1 -1 1]),axis off, axis equal
colormap(flipud(davos));
%
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 19 19]),
print(['fig\Fig9'],'-dpng','-r300')
print(['fig\Fig9'],'-painters','-depsc')