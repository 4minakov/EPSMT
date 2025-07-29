

%Event Magntitude Scaling
Scaling_tbl = [...
     0.2    0.2e2  0.2e3  0.2e5;...%R
    -3.6   -0.8    0.45   3.1;...%M_w^{max}
    -1.6    0.8    2.0    4.7;...%M_w^{summ}
     4.1e3  6.2e7  4.7e9  5.4e13;...%M_0^{max} (N * m)
    ];


figure,

semilogx(Scaling_tbl(1,:), Scaling_tbl(2,:), 'o--', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'r', 'Color',[.3 .3 .3],'LineWidth', 1)
hold on
semilogx(Scaling_tbl(1,:), Scaling_tbl(3,:), 'd--', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'y', 'Color',[.3 .3 .3],'LineWidth', 1)
legend('M_w^{max}','M_w^{\Sigma}','Location','northwest')
xlabel('Radius (m)'), ylabel('Moment magnitude (M_w)')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4]),
print('Fig_MwScaling','-dpng','-r600'), 
print('Fig_MwScaling','-painters','-depsc','-r600'), 
%%
figure,
aa = [log(Scaling_tbl(1,:)') ones(4,1)]\log(Scaling_tbl(4,:)');
xaa = exp(linspace(-2,10,200)); 
loglog(xaa,exp(aa(2))*xaa.^(aa(1)),'Color',[.7 .7 1],'LineWidth', 1.5), axis equal tight
hold on
loglog(Scaling_tbl(1,:), Scaling_tbl(4,:), 'o', 'MarkerEdgeColor','k',...
    'MarkerFaceColor', 'r')
xlabel('R (m)'), ylabel('Scalar moment M_0 (N-m)')
legend('R^2','M_0','Location','northwest')
print('Fig_M0Scaling','-dpng','-r600'), 
print('Fig_M0Scaling','-painters','-depsc','-r600'), 

%%
