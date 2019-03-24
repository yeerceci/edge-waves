clear 
load('PtSalBathy.mat','xl','zl')       

zl1x = gradient(zl, xl(2)-xl(1));  % 1st order partial derivative in x 
zl2x = gradient(zl1x,xl(2)-xl(1));  % 2nd order partial derivative in x 
zl3x = gradient(zl2x,xl(2)-xl(1));

figure
% subplot(411)
plot(xl,zl) 
% subplot(412)
% plot(xl,zl1x)       
% subplot(413)
% plot(xl,zl2x)
% subplot(414)
% plot(xl,zl3x)

x2 = xl(1:548);
z2 = zl(1:548);
p2 = polyfit(x2, z2, 2);
zf2 = polyval(p2, x2);
zf2_1x = gradient(zf2, x2(2)-x2(1));  % 1st order partial derivative in x 
zf2_1x(1)-zl1x(1)


x1 = xl(549:end);
z1 = zl(549:end);
p1 = polyfit(x1, z1, 1);
zf1 = polyval(p1, x1);
zf1_1x = gradient(zf1, x1(2)-x1(1));  % 1st order partial derivative in x 
zf1_1x(end)-zl1x(end)

% subplot(411)
hold on
plot(x1,zf1,'--g')
plot(x2,zf2,'--r')