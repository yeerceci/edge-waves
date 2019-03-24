function l = dispersion(lam,h,eps,n,delta_n)

    k = @(ksi) findk(h(ksi),lam); % find kappa using eq.(4.14) given h and lambda 
%     ksi_c =  @(l) atanh(min(lam/l,1))/l; % at the caustic, upper bound of the integral
    ksi_c =  @(l) atanh(lam/l)/l; % at the caustic, upper bound of the integral
    f1 = @(ksi,l) sqrt(k(ksi).^2-l.^2);

    fint = @(l) integral( @(ksi) f1(ksi,l), 0,ksi_c(l) ) -eps*(n+0.5-delta_n/2)*pi;
    options = optimset('TolFun',1e-6,'TolX',1e-12, 'MaxIter',200, 'MaxFunEvals',200);
%     l = fsolve(fint,1,options);
    l = fsolve(fint,1.5*lam,options); % initial guess of l=lam
    l = real(l); 