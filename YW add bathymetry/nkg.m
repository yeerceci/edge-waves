%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STRATIFIED ROTATING EDGE WAVES OVER A GENERAL SLOPE USING NEWTON-KANTOROVICH ITERATION
%
% x,z are the physical coordinates.
% y is a stretched vertical coordinate, y = F(x)*z, where F(x) = x*h'(0)/h(x).
% 
% Thus a general profile, z = -h(x), is mapped to a wedge, y = -h'(0)*x, of equal initial slope.
%
% Governing equations become:
% 
% S_xx + 2*(F'/F)*y*S_xy + [s^2*F^2 + (F'/F)^2*y^2]*S_yy + (F''/F)*y*S_y - [(s*b/2)^2 + l^2]*S = 0;
%
% F*S_y - lambda*S = 0, on y = 0;
%
% h'*S_x + [s^2*F + h'*(F'/F)*y]*S_y - [s^2*b/2 + h'*f*l/w]*S = 0, on y = h'(0)*x.
%
% - Use Evans solution S(x,y) for equivalent wedge as initial guess.
% - Perform Newton-Kantorovich iterations.
% - Compare results with asymptotics in SGLS.
%
% USING L'HOPITAL'S RULE AT r=0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
% close all
deg = pi/180;

load PtSalBathy

tic

% % physical parameters
% f = 1; % Coriolis frequency (f < B unless B = 0)
% B = 3*f; % buoyancy frequency (denoted N in the literature)
% g = 1;
% n = 2; % mode number
% 
% % Ball's profile: H(x) = c*[1 - exp(-a*x)] => H'(0) = a*c
% e = tan(20*deg); % coastal slope
% c = 1;
% a = e/c;
% 
% w = 1.38*f;

% physical parameters
f = 1e-2; % Coriolis frequency (f < B unless B = 0)
B = 3*f; % buoyancy frequency (denoted N in the literature)
g = 10;
n = 0; % mode number

% Ball's profile: H(x) = c*[1 - exp(-a*x)] => H'(0) = a*c
e = tan(20*deg); % coastal slope
c = 5e3;
a = e/c;

w = 10*f;

% associated quantities
s = sqrt((w^2-f^2)/(w^2-B^2));
b = B^2/g;
lambda = w^2/g - b/2;

% r-coordinate: [0,infty)
N = 32;
rmax = 5*g/w^2;
rpts = lagroots(N-1); % unscaled grid points
scale = rpts(end)/rmax; % scale factor to get r(end) = rmax
[r,D] = lagdif(N,2,scale);
D1r = D(:,:,1);
D2r = D(:,:,2);

% theta-coordinate: [0,-atan(e)]
M = 16;
tmax = atan(e);
[t,D] = chebdif(M,2);
t = tmax*(t-1)/2;
D1t = (2/tmax)*D(:,:,1);
D2t = (2/tmax)^2*D(:,:,2);

% set up tensor product grid
[tt,rr] = meshgrid(t,r);
[xx,yy] = pol2cart(tt,rr); % new vertical coordinate

warning off

%% Modified topography
% YW edit
zl1x = gradient(zl, xl(2)-xl(1));  % 1st order partial derivative in x 
zl2x = gradient(zl1x,xl(2)-xl(1));  % 2nd order partial derivative in x 
zl3x = gradient(zl2x,xl(2)-xl(1));  % 3rd order partial derivative in x 

h   = interp1(xl,zl, xx,'spline');  % size(h) = 32*16
h1x = interp1(xl,zl1x,xx,'spline');
h2x = interp1(xl,zl2x,xx,'spline');
h3x = interp1(xl,zl3x,xx,'spline');

% SGSL edit
% h1x = D1r*h*diag(cos(t)) - diag(1./r)*D1t*h*diag(sin(t));  % typo
% h2x = cos(tt).*(D1r*h1x) - sin(tt)./r.*(D1t*h1x)  % use rr instead of r
% h3x = cos(tt).*(D1r*h2x) - sin(tt)./r.*(D1t*h2x)  % use rr instead of r

% h1x = D1r*h*diag(cos(t)) - diag(1./r)*h*D1t*diag(sin(t)); % size(h1x) = 32*16
% h2x = cos(tt).*(D1r*h1x) - sin(tt)./rr.*(h1x*D1t)  % size(h2x) = 32*16
% h3x = cos(tt).*(D1r*h2x) - sin(tt)./rr.*(h2x*D1t)  % size(h3x) = 32*16
% h2x =s D1r*h1x*diag(cos(t))% - diag(1./r)*h1x*D1t*diag(sin(t))  % size(h2x) = 32*16
% h3x = D1r*h2x*diag(cos(t))% - diag(1./r)*h2x*D1t*diag(sin(t))  % size(h3x) = 32*16

% return

F = h1x(1)*xx./h;
F1 = 1./xx - h1x./h; % F'/F
F2 = -(2*h1x.*F1 + h2x)./h; % F''/F

%% Ball's bottom and associated functions
% h = c*(1 - exp(-a*xx));
% h1x = a*c*exp(-a*xx);
% h2x = -a*h1x;
% h3x = -a*h2x; % need this for l'Hopital's rule
% F = h1x(1)*xx./h;
% F1 = 1./xx - h1x./h; % F'/F
% F2 = -(2*h1x.*F1 + h2x)./h; % F''/F

% tanh bottom
% h = c*tanh(a*xx);
% h1x = a*c*sech(a*xx).^2;
% h2x = -2*h.*h1x;
% h3x = -2*(h1x.^2 + h.*h2x);
% F = a*c*xx./h;
% F1 = 1./xx - h1x./h; % F'/F
% F2 = -(2*h1x.*F1 + h2x)./h ; % E = G''/G

% general bottom
% h = c*xx.*(xx+a);
% h1x = c*(2*xx+a);
% h2x = 2*c;
% h3x = 0;
% F = a*c*xx./h;
% F1 = 1./xx - h1x./h; % F'/F
% F2 = -(2*h1x.*F1 + h2x)./h; % F''/F

% wedge bottom
% h = a*c*xx;
% h1x = a*c*ones(size(xx));
% h2x = zeros(size(xx));
% h3x = h2x;
% F = a*c*xx./h;
% F1 = 1./xx - h1x./h; % F'/F
% F2 = -(2*h1x.*F1 + h2x)./h; % F''/F

% wedge bottom (perturbed)
% h = a*c*xx + (a*c)^2*sin(xx);
% h1x = a*c + (a*c)^2*cos(xx);
% h2x = -(a*c)^2*sin(xx);
% h3x = -(a*c)^2*cos(xx);
% F = a*c*xx./h;
% F1 = 1./xx - h1x./h; % F'/F
% F2 = -(2*h1x.*F1 + h2x)./h; % F''/F

warning on

% Apply l'Hopital's rule at x = 0 (note: extremely tedious for E)
F(1,:) = 1;
F1(1,:) = -h2x(1)/2/h1x(1);
F2(1,:) = (3*h2x(1)-2*h1x(1)*h3x(1))/(6*h1x(1)^2);

% physical vertical coordinate
zz = yy./F;

% differentiation matrices
C = diag(cos(t));
CC = diag(cos(t).^2);
C2 = diag(cos(2*t));
S = diag(sin(t));
SS = diag(sin(t).^2);
S2 = diag(sin(2*t));
R1 = diag(r.^-1);
R2 = diag(r.^-2);

D1X = kron(C,D1r) - kron(S*D1t,R1);
D1Y = kron(S,D1r) + kron(C*D1t,R1);
D2X = kron(CC,D2r) + kron(SS,R1*D1r) + kron(SS*D2t,R2) + kron(S2*D1t,R2-R1*D1r);
D2Y = kron(SS,D2r) + kron(CC,R1*D1r) + kron(CC*D2t,R2) - kron(S2*D1t,R2-R1*D1r);
DXY = kron(S2/2,D2r-R1*D1r) - kron(S2/2*D2t,R2) + kron(C2*D1t,R1*D1r-R2);
I = eye(M*N);

% l'Hopital's rule at r=0
D1X0 = kron(C-S*D1t,D1r);
D1Y0 = kron(S+C*D1t,D1r);
D2X0 = kron(eye(M)+SS*D2t/2-S2*D1t/2,D2r);
D2Y0 = kron(eye(M)+CC*D2t/2+S2*D1t/2,D2r);
DXY0 = kron(C2*D1t/2-S2*D2t/4,D2r);
r0pts = 1:N:M*N;
D1X(r0pts,:) = D1X0(r0pts,:);
D1Y(r0pts,:) = D1Y0(r0pts,:);
D2X(r0pts,:) = D2X0(r0pts,:);
D2Y(r0pts,:) = D2Y0(r0pts,:);
DXY(r0pts,:) = DXY0(r0pts,:);

% clear redundant large variables
clear D D1t D2t D1r C CC C2 S SS S2 R1 R2 D1X0 D1Y0 D2X0 D2Y0 DXY0

% calculate parts of N-K matrices that do not change with every iteration
P = D2X + longdiag(2*F1.*yy)*DXY + longdiag(s^2*F.^2 + (F1.*yy).^2)*D2Y + longdiag(F2.*yy)*D1Y;
Q = longdiag(F)*D1Y - lambda*I;
R = longdiag(h1x)*D1X + longdiag(s^2*F + h1x.*F1.*yy)*D1Y;

% clear again
clear D1X D1Y D2X D2Y DXY

disp(['Differentiation matrices constructed in ',num2str(toc),'s.'])

% use Evans solution in y as initial guess
[lp,lm,ssp,ssm] = evans(xx,yy,atan(e),w,g,f,B,n);
% l = lp; ssguess = ssp; % choose +/- branch (coast on L/R)
l = lm; ssguess = ssm; % choose +/- branch (coast on L/R)
% if isnan(l)==1
%     error('No trapped Evans mode exists for these parameters.')
% end

% load previous.mat, ssguess = ss;

% pick normalisation point at which s is held fixed (NOT on boundary or nodal line)
j = 1;

% iterations
tic
itns = 0;
lhist = l;
tol = 1e-9; % terminate when dl/l < tol
dvec = l; % doesn't matter, just make it larger than tol*l
svec = ssguess(:);
while abs(dvec(end)/l) > tol
% left side
    L =  P - ((s*b/2)^2 + l^2)*I;
    BC1 = Q;
    L(1:N,:) = BC1(1:N,:); % theta = 0
    BC2 = R - longdiag(s^2*b/2 + f/w*l*h1x);
    L(end-N+1:end,:) = BC2(end-N+1:end,:); % theta = -beta
% right side
    rvec = -L*svec;
% add extra row and column
    L = [L,-2*l*svec];
    L(1:N,end) = 0;
    L(end-N+1:end,end) = -f/w*h1x(:,end).*svec(end-N+1:end);
    L(end+1,j) = 1;
    rvec = [rvec;0];
% solve system
    dvec = L\rvec;
    itns = itns+1; % number of iterations
    svec = svec+dvec(1:end-1);
    l = l+dvec(end);
    lhist = [lhist;l];
% break out of loop if too many iterations
    if itns == 20
        disp('Too many iterations. Something is wrong...')
        break
    end    
end
ss = reshape(svec,N,M);
disp(['Performed ',num2str(itns),' iterations in ',num2str(toc),' seconds.'])

% comapre with Ball's result
% lball = l_ball_p(w,e,c,g,n,f);
% ['Evans    ',num2str(lp,'%.9g')]
% ['Ball     ',num2str(lball,'%.9g')]
% ['Numeric  ',num2str(l,'%.9g')]

% plot initial guess
figure(1), clf
subplot(2,2,1)
pcolor(xx,zz,ssguess)
shading interp
% axis image
xlabel(['Initial S(x,z)'])
title(['\omega = ',num2str(w)])

% plot final efun
subplot(2,2,2)
pcolor(xx,zz,ss)
shading interp
% axis image
xlabel(['Final S(x,z)'])
title(['k = ',num2str(l)])

% plot eval and relative error at each iteration
subplot(2,2,3)
plot(0:itns,lhist,'.-','markersize',12)
xlabel('itns')
ylabel('k')
subplot(2,2,4)
plot(1:itns,log10(abs(lhist(1:end-1)-lhist(2:end))),'marker','.','markersize',12)
xlim([1,itns])
grid on
xlabel('itns')
ylabel('log \delta')

figure(2), clf
subplot(2,1,1)
plot(r,ssguess(:,1),':k'), hold on, plot(r,ss(:,1),'-k')
axis tight
ylabel('S(x,0)')
title(['\it l = ',num2str(l),'   rel err = ',num2str(abs(1-l/lm))])
subplot(2,1,2)
plot(r,ssguess(:,end),':k'), hold on, plot(r,ss(:,end),'-k')
axis tight
xlabel('x')
ylabel('S(x,-h)')


% flag = input('Save data? Enter 1 for Yes, 0 for No......  ');
% if flag == 1
%     save previous.mat l ss
% %    load wedge_fN_0m.mat
%     wvec(end+1) = w;
%     lvec(end+1) = l;
%     ssi(:,:,end+1) = ss;
%     ssi=ss;
%     xxi(:,:,end+1) = xx;
%     xxi=xx;
%     zzi(:,:,end+1) = zz;
%     zzi=zz;
%     lmvec(end+1) = lm;
%     ssmi(:,:,end+1) = ssm;
%     save wedge_fN_0m.mat wvec lvec ssi xxi zzi lmvec ssmi
% end