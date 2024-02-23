function df = estFirstDerivative(f, h)
% estimate 1st derivative (velocity) using: centered difference formulas for five-point stencils
% see wikipedia

% (c) Jiri, Mar11

df = nan(size(f));

df(1,:) = (f(2,:) - f(1,:))./h;
df(2,:) = (f(3,:) - f(2,:))./h;
%c = repmat([-1, 8, -8, 1]', 1, size(f,2));      % BUG (found 21.07.2017)
c = repmat([1, -8, 8, -1]', 1, size(f,2));
for j = 3:size(f,1)-2
    df(j,:) = sum([f(j-2,:); f(j-1,:); f(j+1,:); f(j+2,:)].*c, 1)/(12*h);
end
df(end-1,:) = (f(end-1,:) - f(end-2,:))./h;
df(end,:)   = (f(end,:) - f(end-1,:))./h;

assert(isempty(find(isnan(df), 1, 'first')));