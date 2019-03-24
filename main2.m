% Find the wavenumber corresponding to a single frequency,
% and double check the residue of the integral equation.

clear
eps = tan(20/180*pi); % a slope of 20 degree
h = @(ksi) ksi; % bathymetry: h(ksi) or H(x), where ksi = eps*x 
n = 0; % Mode number
delta0 = 2.8e-2; % from Table 1
delta1 = 1.6e-2;
delta2 = 1.1e-2;

lam = 0.02; % frequency: lambda

k = @(ksi) findk(h(ksi),lam); % find kappa using eq.(4.14) given h and lambda 
ksi_c =  @(l) atan(lam/l)/l; % at the caustic, upper bound of the integral
f1 = @(ksi,l) sqrt(k(ksi).^2-l.^2);

%% Zhevandrov
fint = @(l) integral( @(ksi) f1(ksi,l), 0,ksi_c(l) ) -eps*(n+0.5)*pi; %

%% Shen & Keller
% fint = @(l) integral( @(ksi) f1(ksi,l), 0,ksi_c(l) ) -eps*(n+0.5-delta0/2)*pi; 

options = optimset('TolFun', 1e-6, 'TolX', 1e-12 );
l = fsolve(fint,0,options)  % initial guess of 0 is way more fast

%% Doublecheck the residue of the integral equation
l_list = linspace(0,0.3,31);
res_list = zeros(size(l_list));

for i = 2:length(l_list), res_list(i) = real(fint(l_list(i))); end

figure(3), hold on
plot(l_list,res_list,'linewidth',1)
plot(l,0,'x','linewidth',2,'markersize',8)

xlabel('wavenumber {\it l }')
ylabel('residue of eq.(4.24)')
title('Zhevandrov')
% legend([1,3,5,7],'\lambda = 0.1','\lambda = 0.08','\lambda = 0.06','\lambda = 0.04','\lambda = 0.02')

% plot(l_list, zeros(size(l_list)), '--')
set(gca, 'FontSize',14)