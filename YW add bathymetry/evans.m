% FUNCTION TO CALCULATE EVANS EIGENVALUES AND EIGENFUNCTIONS (NORMALISED TO UNIT MAGNITUDE AT COAST)

function [lp,lm,vp,vm] = evans(xx,zz,delta,w,g,f,N,n)

% xx,zz = coordinates for unscaled wedge in mesh format
% delta = wedge angle
% w = frequency
% g = gravity
% f = coriolis frequency
% N = buoyancy frequency
% n = mode (0,1,2,...)

% associated quantities
b = N^2/g;
s = sqrt((w^2-f^2)/(w^2-N^2));
lambda = s*(w^2/g - b/2);
beta = atan(tan(delta)/s);
yy = zz/s;

% find l+ and l- by solving quadratic equation
P = (f/w*sin(beta))^2 - sin((2*n+1)*beta)^2;
Q = f/w*sin(beta)*(s*b*cos(beta) - 2*lambda*cos((2*n+1)*beta));
R = lambda^2 - s*b*lambda*cos(beta)*cos((2*n+1)*beta) + (s*b/2)^2*(cos(beta)^2 - sin((2*n+1)*beta)^2);
lp = (-(Q/P)+sqrt((Q/P)^2-4*R/P))/2;
lm = (-(Q/P)-sqrt((Q/P)^2-4*R/P))/2;

% suppress values of l which do not satisfy dispersion relation (although they do satisfy the quadratic)
% also suppress values of l which do not satisfy the trapping condition
kp = sqrt(lp^2+(s*b/2)^2);
ap = -(f/w*lp*sin(beta) + s*b/2*cos(beta));
right_branch = abs( lambda + ap*cos((2*n+1)*beta) - sqrt(kp^2-ap^2)*sin((2*n+1)*beta) );
wrong_branch = abs( lambda + ap*cos((2*n+1)*beta) + sqrt(kp^2-ap^2)*sin((2*n+1)*beta) );
lp( find(right_branch > wrong_branch) ) = NaN;
trap = real( (ap+lambda*cos((2*n+1)*beta))/sin((2*n+1)*beta) );
for m=(2*n+1):-1:1
    trap = min( trap, real( (ap+lambda*cos(m*beta))/sin(m*beta) ) );
end
lp(find(trap<=0)) = NaN;

km = sqrt(lm^2+(s*b/2)^2);
am = -(f/w*lm*sin(beta) + s*b/2*cos(beta));
right_branch = abs( lambda + am*cos((2*n+1)*beta) - sqrt(km^2-am^2)*sin((2*n+1)*beta) );
wrong_branch = abs( lambda + am*cos((2*n+1)*beta) + sqrt(km^2-am^2)*sin((2*n+1)*beta) );
lm( find(right_branch > wrong_branch) ) = NaN;
trap = real( (am+lambda*cos((2*n+1)*beta))/sin((2*n+1)*beta) );
for m=(2*n+1):-1:1
    trap = min( trap, real( (am+lambda*cos(m*beta))/sin(m*beta) ) );
end
lm(find(trap<=0)) = NaN;

% calculate wave function using Evans expression

% l+ function (coast on R)
l = lp;
if isnan(l)
    vp = NaN;
else
	k = kp;
	alpha = ap;
	chi = asin(alpha/k);
	v(:,:,1) = exp( -k*( xx*cos(beta-chi) - yy*sin(beta-chi) ) );
%     disp(['(2n+1)beta-chi = ',num2str((2*n+1)*beta-chi)])
%     disp(['m = 0 exponent is ',num2str(-k*cos(beta-chi))])
	for m = 1:n
%        disp(['m = ',num2str(m),' exponent is ',num2str(-k*cos((2*m+1)*beta-chi))])
        v(:,:,m+1) = coeff(m,n,beta,chi)*( exp( -k*( xx*cos((2*m-1)*beta-chi) + yy*sin((2*m-1)*beta-chi) ) )...
            + tan(m*beta-chi)/tan(m*beta)*exp( -k*( xx*cos((2*m+1)*beta-chi) - yy*sin((2*m+1)*beta-chi) ) ) );
	end
	vp = sum(v,3);
	vp = vp/vp(1); % normalise to unit magnitude at origin
end

% l- function (coast on L)
l = lm;
if isnan(l)
    vm = NaN;
else
	k = km;
	alpha = am;
	chi = asin(alpha/k);
	v(:,:,1) = exp( -k*( xx*cos(beta-chi) - yy*sin(beta-chi) ) );
	for m = 1:n
        v(:,:,m+1) = coeff(m,n,beta,chi)*( exp( -k*( xx*cos((2*m-1)*beta-chi) + yy*sin((2*m-1)*beta-chi) ) )...
            + tan(m*beta-chi)/tan(m*beta)*exp( -k*( xx*cos((2*m+1)*beta-chi) - yy*sin((2*m+1)*beta-chi) ) ) );
	end
	vm = sum(v,3);
	vm = vm/vm(1); % normalise to unit magnitude at origin
end

% SUBFUNCTION: coefficient in Evans expansion
% function fun = coeff(m,n,beta,chi)
% for r = 1:m
%     u(r) = tan((n+1-r)*beta)*tan(r*beta-chi)/tan(r*beta)/tan((r+n)*beta-chi);
% end
% fun = (-1)^m*tan(m*beta)/tan(m*beta-chi)*prod(u);