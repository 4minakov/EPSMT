clear all
close all
addpath .\data
%% Radiation pattern
load fm_ISO
load fm_ANISO
%% grid and wave vectors
delta = pi/200;
theta = 0 : delta : pi; % altitude
phi = 0 : 2*delta : 2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta);
%P
kx = sin(theta).*cos(phi);
ky = sin(theta).*sin(phi);
kz = cos(theta);
%SV
% qx =  cos(theta).*cos(phi);
% qy =  cos(theta).*sin(phi);
% qz = -sin(theta);
%SH
qx = -sin(phi);
qy =  cos(phi);
qz =  0*phi;
%
figure
%% P- ISO
rho = size(phi(:));
for i =1:length(phi(:))
    p = [kx(i), ky(i), kz(i)]';
    q = [kx(i), ky(i), kz(i)]';
    pq = 1/2*(p.*q' + q.*p');
    rho(i) = sum(fm_ISO(:).*pq(:));
end
rho_iso = reshape(rho,size(phi));
%rho = abs(rho);
%Apply spherical coordinate equations
r = rho_iso.*sin(theta);
x = r.*cos(phi); % spherical coordinate equations
y = r.*sin(phi);
z = rho_iso.*cos(theta);
%
ax1= subplot(131);
surf(x,y,z,rho_iso,'EdgeColor','none'),
light,lighting phong, axis tight equal off %xlabel('x'),ylabel('y')
view(0,45),
tt =title('(a)      P-ISO');
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
xpos1 = get(ax1,'Position');
%% P- ANISO 
rho = size(phi(:));
for i =1:length(phi(:))
    p = [kx(i), ky(i), kz(i)]';
    q = [kx(i), ky(i), kz(i)]';
    pq = 1/2*(p.*q' + q.*p');
    rho(i) = sum(fm_ANISO(:).*pq(:));
end
rho_aniso = reshape(rho,size(phi));
%rho = abs(rho);
%Apply spherical coordinate equations
r = rho_aniso.*sin(theta);
x = r.*cos(phi); % spherical coordinate equations
y = r.*sin(phi);
z = rho_aniso.*cos(theta);

ax2=subplot(132);
surf(x,y,z,rho_aniso,'EdgeColor','none'),
light,lighting phong, axis tight equal off %xlabel('x'),ylabel('y')
view(0,45),
tt =title('(b)      P-ANISO');
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
xpos2 = get(ax2,'Position');
%% SH- 
rho2c = zeros(length(phi(:)),1);
for i =1:length(phi(:))
    p = [kx(i), ky(i), kz(i)]';
    q = [qx(i), qy(i), qz(i)]';
    pq = 1/2*(p.*q' + q.*p');
    rho2c(i) = sum(fm_ISO(:).*pq(:));
end
rho2c = reshape(rho2c,size(phi));
rho2 = abs(rho2c);
%Apply spherical coordinate equations
r2 = rho2.*sin(theta);
x2 = r2.*cos(phi); % spherical coordinate equations
y2 = r2.*sin(phi);
z2 = rho2.*cos(theta);

ax3 = subplot(133);
surf(x2,y2,z2,rho2c,'EdgeColor','none');
light,lighting phong, axis tight equal off %xlabel('x'),ylabel('y')
tt =title('(c)      SH');
tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left';
view(0,45)
xpos3 = get(ax3,'Position');
set(ax1,'Position',[xpos1(1),xpos1(2),xpos1(3),xpos1(4)-.1])
set(ax2,'Position',[xpos2(1),xpos2(2)+.13,xpos2(3),xpos2(4)])
set(ax3,'Position',[xpos3(1),xpos3(2)+.07,xpos3(3),xpos3(4)])
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 17 9]),
print(['fig\Fig10'],'-dpng','-r300')
print(['fig\Fig10'],'-painters','-depsc')