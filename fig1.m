clear all
close all
addpath ./data
C = 40;
m = 2;
K = (m-1)/(m+1);
K = sind(30);

mc  = @(x1, x2, x3) max(cat(4,abs(x1-x2) - K*(x1+x2),...
    abs(x2-x3)- K*(x2+x3),...
    abs(x1-x3)- K*(x1+x3),...
    (0.2*C-x3),...
    (0.2*C-x1),...
    (0.2*C-x2)...
),[],4);

x1 = C*linspace(-2, 2, 201);
x2 = C*linspace(-2, 2, 201);
x3 = C*linspace(-2, 2, 201);
[X1, X2, X3] = meshgrid(x1, x2, x3);
%% Figure 1a
figure
p = patch(isosurface(X1,X2,X3,mc(X1, X2, X3),C));
isonormals(X1,X2,X3,mc(X1, X2, X3),p)
p.FaceColor = [0 0.5 0.5];
p.EdgeColor = 'none';
p.FaceAlpha = 0.7;
daspect([1 1 1])
view(3); 
axis tight
camlight(0,30)
lighting gouraud
hold on,
plot3(x1(30:end),x2(30:end),x3(30:end),'LineWidth',2,'Color','k')
xlabel('\sigma_1 (MPa)'),ylabel('\sigma_2'),zlabel('\sigma_3 (MPa)')
grid off, box on, view(20,30)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4]),
print('fig\Fig1a','-painters','-depsc','-r600'),
%% Figure 1b
ind = find(p.Vertices(:,2)==0);
s_n = (p.Vertices(ind,1)+p.Vertices(ind,3))/2;
tau_m = (p.Vertices(ind,1)-p.Vertices(ind,3))/2;
s_n(tau_m<0)=[];
tau_m(tau_m<0)=[];
figure,
plot(s_n,tau_m,'-','Color',[0 0.5 0.5],'LineWidth',1.5), axis equal, axis([min(s_n),40,0,40]) 
xlabel('\sigma_n (MPa)'),ylabel('\tau (MPa)')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]),
print('fig\Fig1b','-depsc','-r600'),
