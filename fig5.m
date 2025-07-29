clear all
close all
addpath .\data
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
xlim([0.1 200])
print('fig\Fig5','-dpng','-r600'), 
print('fig\Fig5','-painters','-depsc','-r600'), 

