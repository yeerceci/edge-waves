% DISPERSION RELATION
% The method is to fix frequency (lambda) and seek wavenumber (l)
% Outputs: l_z: numerical solution using Zhevandrov's integral equation
%          l_s: numerical solution using Shen & Keller's integral equation
%          l_ursell: analytical solution of Ursell
%          l_zhevandrov: analytical solution of Zhevandrov
%          l_shen: analytical solution of Shen & Keller

clear
eps = tan(20/180*pi); % a slope of 20 degree
h = @(ksi) ksi; % bathymetry: h(ksi) or H(x), where ksi = eps*x 
n = 0; % Mode number
delta0 = 2.8e-2;
delta1 = 1.6e-2;
delta2 = 1.1e-2;

%% range of frequency
lam_list = linspace(0,0.1,10);
l_z = zeros(size(lam_list));
l_s = l_z;

for i = 2:length(lam_list)
    lam = lam_list(i);
    l_z(i) = dispersion(lam,h,eps,n,0);
    l_s(i) = dispersion(lam,h,eps,n,delta0);
end

%% compare with exact and Ursell's dispersion
l_ursell = lam_list/sin( (2*n+1)*atan(eps) );
l_zhevandrov = lam_list/sin( (2*n+1)*eps );
l_shen = lam_list/sin( (2*n+1-delta0)*eps);

% close
figure(1), hold on
p_z = plot(l_z,lam_list,'x','linewidth',2,'markersize',8);
p_s = plot(l_s,lam_list,'x','linewidth',2,'markersize',8);
p_ursell = plot(l_ursell, lam_list,'linewidth',1);
p_zhevandrov = plot(l_zhevandrov, lam_list,'linewidth',1);
p_shen = plot(l_shen, lam_list,'linewidth',1);


xlabel('l')
ylabel('lambda')
title('Dispersion relation')
legend([p_z,p_s,p_ursell,p_zhevandrov,p_shen], 'My numerical Zhevandrov','My numerical Shen & Keller','Ursell','Zhevandrov','Shen & Keller','location','southeast')
set(gca, 'FontSize',14)

% save('dispersion.mat')