function [y_fit, x_fit] = datafit_smoothingSplines(x, y)
% fits smoothing splines, estimates the smoothing coefficient such that the
% R-square value => 0.90

% (c) Jiri, Dec16

assert(size(x,1) == size(y,1));

%% settings for fitting
rSq_lim = 0.50;
rSq_tol = 0.01;

%% initial fits
% upper bound
pSm_hi = 1.0;   % param for smoothing
[fitObj, gof] = fit(x, y, 'smoothingspline', 'SmoothingParam',pSm_hi);
rSq_hi = gof.rsquare;
assert(rSq_hi >= rSq_lim);

% lower bound
pSm_lo = 0.0;
[fitObj, gof] = fit(x, y, 'smoothingspline', 'SmoothingParam',pSm_lo);
rSq_lo = gof.rsquare;
assert(rSq_lo <= rSq_lim);

% init
smP_new = mean([pSm_lo,pSm_hi],2);
[fitObj, gof] = fit(x, y, 'smoothingspline', 'SmoothingParam',smP_new);
rSq_new = gof.rsquare;

%% find "OK" smoothing param
k = 1;
while abs(rSq_lim - rSq_new) > rSq_tol
    %display(['Iteration: ' num2str(k) ', smoothParam = ' num2str(smP_new) ', r-square = ' num2str(rSq_new)]);
    
    % find new bounds
    if rSq_new < rSq_lim
        smP_new = mean([smP_new,pSm_hi],2);
    elseif rSq_new > rSq_lim
        smP_new = mean([smP_new,pSm_lo],2);
    else
        error('should have converged');
    end
    
    % new fit
    [fitObj, gof] = fit(x, y, 'smoothingspline', 'SmoothingParam',smP_new);
    rSq_new = gof.rsquare;   
    k = k+1;
    if mod(k,100) == 1
        rSq_tol = rSq_tol+0.01;
        display(['WARNING: k = ' num2str(k) ', convergence problem. Increasing rSq_tol = ' num2str(rSq_tol)]);
    end
end
display(['Final round: ' num2str(k) ', smoothParam = ' num2str(smP_new) ', r-square = ' num2str(rSq_new)]);

%% interpolation
x_fit = [x(1):0.01:x(end)]';
y_fit = feval(fitObj, x_fit);
assert(size(x_fit,1) == size(y_fit,1));

        