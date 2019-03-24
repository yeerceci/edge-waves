clear 
% close
iBC = 10;

deg = pi/180;
e = tan(20*deg); % coastal slope
c = 5e3;
a = e/c;

if iBC == 0.1, % Ball
    xl = 0:2e2:20e4;
    zl = c*(1 - exp(-a*xl));  
    save('PtSalBathy_Ball_test','xl','zl')
elseif iBC == 0.4,   % wedge
    xl = 0:1e3:3125e3;  
    zl = a*c*xl;
    save('PtSalBathy_wedge_test','xl','zl')
else % parabolic
    load('PtSalBathy.mat','xl','zl')    
    
    x2 = xl(1:548);
    z2 = zl(1:548);
    p2 = polyfit(x2, z2, 2);
    ratio2 = tan(20*deg)/p2(2); % scale the initial slope to 20 degree

    x1 = xl(549:end);
    z1 = zl(549:end);
    p1 = polyfit(x1, z1, 1);
    ratio1 = tan(20*deg)/p1(1); % scale the finas slope to 20 degree

    zl = zl * ratio1;
    save('PtSalBathy_Bathy_test','xl','zl')
end

figure
plot(xl,zl)
axis equal
% save('PtSalBathy_Ball_test')
