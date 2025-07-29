%% theoretical Hudson diagram
lam = E*nu/(1+nu)/(1-2*nu);
% % plane strain
% D   = E*(1-nu)/(1+nu)/(1-2*nu) * [ 1         nu/(1-nu) 0;...
%     nu/(1-nu) 1         0;...
%     0         0         (1-2*nu)/2/(1-nu) ];
npnt = 1e3;
epl = [linspace(1e-6,1,npnt),linspace(-1,-1e-6,npnt)];
Mst_0 = zeros(npnt,1); tau_a = Mst_0; k_a = Mst_0;
%shear fault (x1-x3 plane, dipl in x1)
M_shea = [0 mu 0; mu 0 0; 0 0 0];
%tensile fault (x1-x3 plane, opens in x2)
M_tens = [lam 0 0; 0 lam+2*mu 0; 0 0 lam];
theta = linspace(0,pi/2,npnt);
for i=1:length(epl)
    ex  = 2*rand-1;
    ey  = 2*rand-1;
    exy = 2*rand-1;
    e_src = [ex exy 0; exy ey 0; 0 0 0];
    
    M_st = lam*eye(3)*trace(e_src) + 2*mu*e_src;
    
    Mst_0(i) = sum(sqrt(M_st(:).^2));
    [tau_a(i),k_a(i)] = MT2tauk(M_st);  
end
[~,ii]=sort(k_a);
tau_a= tau_a(ii);
k_a  = k_a(ii);
%
% figure
% %plot(tau_a,k_a,'r','LineWidth',2)
% hold on, 
% scatter(tau_a,k_a,10,log10(Mst_0),'filled'),hold on, 
% plot([-1:.1:0],1+[-1:.1:0],'k','LineWidth',2)
% plot([0:.1:1],1-[0:.1:1],'k','LineWidth',2)
% plot([0:.1:1],-1+[0:.1:1],'k','LineWidth',2)
% plot([-1:.1:0],-1-[-1:.1:0],'k','LineWidth',2)
% plot([-1,1],[0, 0],'k','LineWidth',1)
% plot([0,0],[-1, 1],'k','LineWidth',1)
% text(0,1.1,'k (0, 1)'), text(1.0,0,'\tau (0, 1)')
% text(0,-1.1,'(0, -1)'), text(-1.2,0,'(-1, 0)')
% %plot(tau,k,'+k','MarkerSize',15),
% axis([-1 1 -1 1]),axis off, axis equal
% colormap(jet)