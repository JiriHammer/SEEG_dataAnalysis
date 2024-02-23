function yOut = getBinCenters(y, doLinear)
% returns centers of the bin boundaries 'y'
% if doLinear = true, then it finds the first and last bin center values
% 'yOut' based on linear extrapolation of the inner bins

% (c) Jiri, Jan12, Jan13

if nargin < 2
    doLinear = false;
end

yC = nan(1,length(y)-1);
for n = 1:length(y)-1
    yC(n) = (y(n) + y(n+1))/2;
end

% opt: linearize the outer bin values
if doLinear
    xBin = 1:length(yC);            % bin index
    b = regress(yC(2:end-1)', cat(2, ones(size(yC(2:end-1),2),1), xBin(2:end-1)'));    % first and last bins were
    dx = diff(xBin);
%     y_beg = yC(2) - dx(1)*b(2);
%     y_end = yC(end-1) + dx(end)*b(2);
%     yOut = [y_beg, yC(2:end-1), y_end];
    
    inds = 2:5;
    b = regress(yC(inds)', cat(2, ones(length(inds),1), xBin(inds)'));    % first bins
    y_beg = yC(2) - dx(1)*b(2);
    
    inds = length(yC)-3:length(yC)-1;
    b = regress(yC(inds)', cat(2, ones(length(inds),1), xBin(inds)'));    % last bins 
    y_end = yC(end-1) + dx(end)*b(2);
    
    yOut = [y_beg, yC(2:end-1), y_end];
else
    yOut = yC;
end