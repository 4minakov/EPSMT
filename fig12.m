clear all
close all
%%
addpath .\data
load(['model_parameters_MC2.mat']);
load(['Y0_MC2']);
load Data_total_MC2_new.mat
load bplot_data
load('davos')
load('batlow')
%sim_type = 'ISO1';
tit = {'(a)','(b)','(c)'};%dt=50
tit = {'(d)','(e)','(f)'};%dt=75
tit = {'(g)','(h)','(i)'};%dt=100
tit = {'(j)','(k)','(l)'};%dt=125
tstep = 50;%50, 75, 100, 125
tstep = 75;%50, 75, 100, 125
tstep = 100;%50, 75, 100, 125
tstep = 125;%50, 75, 100, 125
load res50_MC2%50, 75, 100, 125
load res75_MC2
load res100_MC2
load res125_MC2
%%
%
D1 = inv(D);%Compliance matrixclear
% Grid setup
dth = 2*pi/(ny-1);
dr  = ((Lr)^(1/2) - R0^(1/2))/(nx-1);  % for adaptive grid
r   = R0^(1/2):dr:(Lr)^(1/2);% adaptive spatial grid
r   = r.^2;
th       = 0:dth:2*pi;
[Xr,Yth] = meshgrid(r,th(1:end-1));       % regular grid
x        = Xr.*cos(Yth);                  % deformed grid
y        = Xr.*sin(Yth);
[l2g,l2gu] = l2g4ncir(nx,ny);    % local to global numbering
dr0 = diff(r);
dA = (dth*r(1:end-1).*dr0)';%element area
dA = repmat(dA,1,ny);
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
[maxM0,indmaxM0] = max(M0(:));
Mw = 2/3*log10(M0)-6;
TAU = zeros(length(M0(:)),1);
K = zeros(length(M0(:)),1);
for n = 1:length(M0(:))
    MT = [Mxx(n), Mxy(n), 0; Mxy(n), Myy(n), 0; 0 0 Mzz(n)];
    [TAU(n),K(n)] = MT2tauk(MT);
end
TAU = reshape(TAU,size(M0));
K = reshape(K,size(M0));
%%
figure(1), clf
pcolor(xp/R0,yp/R0,1e-5*M0/R0^2), hold on
hold on, quiver(xp(1:25:end,1:15:end)/R0,yp(1:25:end,1:15:end)/R0,...
    ux(1:25:end-1,1:15:end)/R0,uy(1:25:end-1,1:15:end)/R0,0.3,'Color',[0.4 0.4 0.4])
shading interp, 
colormap(flipud(davos)); 
plot(xp(1,:)/R0,yp(1,:)/R0,'k','LineWidth',1)
axis equal tight, xlabel('x/R0'),ylabel('y/R0')
c = colorbar('position',[0.475 0.169047619047619 0.302857142857142 0.0452380952380957],'Location','southoutside');
c.Label.FontWeight='bold'; c.Label.FontSize=12;
c.Label.String='$M_0$';c.Label.Position=[0 0];
c.Label.Interpreter='latex';
xlabel('$x/R_0$','Interpreter','latex'),ylabel('$y/R_0$','Interpreter','latex'),
tt=title([tit{1},'  $\Delta t_n =$',num2str(tstep)],'Interpreter','latex','FontSize',12);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
axis([-1 1 -1 1]*2.5),
set(gcf,'units','normalized','Position',[0.0167 0.5725 0.2917 0.3500])
set(gca,'FontSize',12)
fout=['fig\Fig12',tit{1}];
print('-depsc','-painters','-r600',fout)
print('-dpng','-r600',fout)
%%
S = log10(1e-5*M0/R0^2); S=S(S>-3);
n  = 100;
Si = linspace(-3,0,n);
N  = zeros(size(Si));
for ii = 1:n-1
    N(ii) = length(find( S>Si(ii) ) );
end
%
ff = histogram(S,'Normalization','pdf','BinMethod','scott','FaceColor','k','FaceAlpha',0.2,'Visible','off');
[max_val,max_ind]=max(ff.Values);
Mc = 0.5*(ff.BinEdges(max_ind)+ff.BinEdges(max_ind+1)) ;
%
figure(2),clf% bplot
xm = Si(Si>Mc & N>1)';
ym = log10(N(Si>Mc & N>1))'; 
xm = [xm,xm*0+1];
b_ls = xm\ym;
b_mle = log10(exp(1))/(mean(xm(:,1))-Mc);
yb = interp1(Si,  log10(N), Mc) + Mc*b_mle;
plot(Si,  log10(N),'o'), hold on
plot(Si, -Si*b_mle+yb,'r','LineWidth',1.5)
ylabel('$\log N(>M)$','Interpreter','latex'),
xlabel('Magnitude ($M$)','Interpreter','latex')
text(-1.5, 2,['$b_{MLE}=$',num2str(round(b_mle,2))],...
   'Interpreter','latex','FontSize',16);
grid on
tt= title([tit{2},'     $N = $',num2str(length(S)),', $M_c =$',num2str(Mc)],...
    'Interpreter','latex','FontSize',14);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
ylim([1 3.5]),xlim([-2 0]), 
set(gca,'FontSize',14)
set(gcf,'units','normalized')
set(gcf,'Position',[0.3109 0.5567 0.2755 0.3500]);
fout=['fig\Fig12',tit{2}];
print('-depsc','-painters','-r600',fout)
print('-dpng','-r600',fout)
%%
figure(3),clf
S  = 1e-5*M0/R0^2; S = S(S>1e-2);
pd_e = fitdist(S,'Exponential');
pd_n = fitdist(S,'Normal');
pd_w = fitdist(S,'Weibull');
b=2; x=sort(S);
%power law
pdf_pow = (b-1)/min(x)*(x/min(x)).^(-b);
ff1= histogram(S,'BinMethod','fd','Normalization','pdf','FaceColor','k',...
    'FaceAlpha',0.1,'Visible','on');
pdf_emp = interp1(0.5*(ff1.BinEdges(1:end-1)+ff1.BinEdges(2:end)),ff1.Values,x);
kk = find(~isnan(pdf_emp)); 
pdf_exp = pdf(pd_e,x);
pdf_nor = pdf(pd_n,x);
pdf_wei = pdf(pd_w,x);
hold on
p1=plot(x,pdf_exp,'-k','LineWidth',1.5); Rp1 = 1-sum((pdf_exp(kk)-pdf_emp(kk)).^2)/sum(pdf_emp(kk).^2);
p2=plot(x,pdf_nor,'-g','LineWidth',1.5); Rp2 = 1-sum((pdf_nor(kk)-pdf_emp(kk)).^2)/sum(pdf_emp(kk).^2);
p3=plot(x,pdf_wei,'-b','LineWidth',1.5); Rp3 = 1-sum((pdf_wei(kk)-pdf_emp(kk)).^2)/sum(pdf_emp(kk).^2);
p4=plot(x,pdf_pow,'-r','LineWidth',1.5); Rp4 = 1-sum((pdf_pow(kk)-pdf_emp(kk)).^2)/sum(pdf_emp(kk).^2);
legend([p1,p2,p3,p4],...
    ['Expon (',num2str(round(Rp1,2)),')'],...
    ['Gauss (',num2str(round(Rp2,2)),')'],...
    ['Weib (',num2str(round(Rp3,2)),')'],...
    ['PowLaw (',num2str(round(Rp4,2)),')'] ,'FontSize',12)
xlabel('Seismic moment ($M_0$)','Interpreter','latex'), ylabel('Probability density','Interpreter','latex')
ylim([0 1.1*max(ff1.Values)])
tt=title([tit{3}, ' '],'Interpreter','latex','FontSize',12);
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
set(gca,'FontSize',12)
set(gcf,'units','normalized')
set(gcf,'Position',[0.5875 0.5567 0.2917 0.3500]);
fout=['fig\Fig12',tit{3}];
print('-depsc','-painters','-r600',fout)
print('-dpng','-r600',fout)
