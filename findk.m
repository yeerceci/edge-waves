function k = findk(h,lam)
    % initial value of k is important
    % positive k is selected
    options = optimset('Display','off');
    fk = @(k) k.*tanh(k.*h)-lam; 
    k = fsolve(fk, ones(size(h))*lam,options); % vectorized
