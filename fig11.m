clear all
close all
addpath .\data

load fm_ISO
load fm_ANISO
load fm_TENS

[dd,rr] = eig(fm_ISO);
M = diag(sort(diag(rr),'descend'));
[~,ii]=sort(diag(rr),'descend');
dd = dd(:,ii);
dd = dd';
delta = pi/200;
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : 2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta);

z1 = cos(theta);
x1 = sin(theta).*cos(phi);
y1 = sin(theta).*sin(phi);
xps = x1(1:fix(end/2),:)./(1+z1(1:fix(end/2),:));
yps = y1(1:fix(end/2),:)./(1+z1(1:fix(end/2),:));

%% P- ISO
rho_iso = size(phi(:));
rho_ani = size(phi(:));
rho_ten = size(phi(:));
kx = sin(theta).*cos(phi);
ky = sin(theta).*sin(phi);
kz = cos(theta);
for i =1:length(phi(:))
    p = [kx(i), ky(i), kz(i)]';
    q = [kx(i), ky(i), kz(i)]';
    pq = 1/2*(p.*q' + q.*p');
    rho_iso(i) = sum(fm_ISO(:).*pq(:));
    rho_ani(i) = sum(fm_ANISO(:).*pq(:));
    rho_ten(i) = sum(fm_TENS(:).*pq(:));
end
rho_iso = reshape(rho_iso,size(phi));
rho_ani = reshape(rho_ani,size(phi));
rho_ten = reshape(rho_ten,size(phi));


%P wave amplitude
figure,
subplot(121)
rr = rho_iso(1:fix(end/2),:);
rr = (rr - min(rr(:))); rr = 2*(rr/max(abs(rr(:)))-1/2);
%contourf(xps,yps,rr), shading interp,
hold on
contourf(xps,yps,rr,[-1 0 1],'k')
plot(.98*sin(0:.001:2*pi),.98*cos(0:.001:2*pi),'k','LineWidth',1.5)
axis equal tight off
T = [dd(1,1)/(1+dd(1,3)), dd(1,2)/(1+dd(1,3))];
P = [dd(3,1)/(1+dd(3,3)), dd(3,2)/(1+dd(3,3))];
B = [dd(2,1)/(1+dd(2,3)), dd(2,2)/(1+dd(2,3))];
hold on, plot(P(1),P(2),'or','MarkerSize',10),plot(T(1),T(2),'+k','MarkerSize',10)
plot(B(1),B(2),'^g','MarkerSize',10)
text(P(1)-0.1,P(2)-0.1,'P','FontSize',14,'Color','b'),text(T(1)-0.1,T(2)-0.1,'T','FontSize',14,'Color','b')
text(B(1)-0.1,B(2)-0.1,'B','FontSize',14,'Color','b')
colormap(flipud(gray(2)))
tt = title('(a)     ISO');
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
%
subplot(122)
rr = rho_ani(1:fix(end/2),:);
rr = (rr - min(rr(:))); rr = 2*(rr/max(abs(rr(:)))-1/2);
%contourf(xps,yps,rr), shading interp,
hold on
contourf(xps,yps,rr,[-1 0 1],'k')
plot(.98*sin(0:.001:2*pi),.98*cos(0:.001:2*pi),'k','LineWidth',1.5)
axis equal tight off
T = [dd(1,1)/(1+dd(1,3)), dd(1,2)/(1+dd(1,3))];
P = [dd(3,1)/(1+dd(3,3)), dd(3,2)/(1+dd(3,3))];
B = [dd(2,1)/(1+dd(2,3)), dd(2,2)/(1+dd(2,3))];
hold on, plot(P(1),P(2),'or','MarkerSize',10),plot(T(1),T(2),'+k','MarkerSize',10)
plot(B(1),B(2),'^g','MarkerSize',10)
text(P(1)-0.1,P(2)-0.1,'P','FontSize',14,'Color','b'),text(T(1)-0.1,T(2)-0.1,'T','FontSize',14,'Color','b')
text(B(1)-0.1,B(2)-0.1,'B','FontSize',14,'Color','b')
colormap(flipud(gray(2)))
tt = title('(b)     ANISO');
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';

set(gcf,'Units','inches','PaperSize',[7 4])
print('-dpng','-r300','fig\Fig11')
print('-depsc','-painters','fig\Fig11')
