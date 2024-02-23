function [ylimits, outliers] = getYLims(statData, stdThr)
% ylimits: cats values in statData (3D: samples x channels x trials) together,
% fits normal distribution & sets y-limits at 6 sigma
% outliers: if a value of a channel is out at 6 sigma of the normal
% disctribution of all channels, the channel is marked as an outlier

% (c) Jiri, Mar16

% threshold
if nargin < 2
    stdThr = 6;     % in [SD]
end

% fit normal distribution
x = statData(:);
%figure; histfit(x);
x(isnan(x)) = []; assert(~isempty(x));
[muhat,sigmahat] = normfit(x);

% ylimits
ylimits = [muhat - stdThr*sigmahat, muhat + stdThr*sigmahat];
assert(ylimits(1) < ylimits(2));

% outliers
if nargout == 2
    outliers = [];
    valThr = muhat + stdThr * sigmahat;
    for clz = size(statData,3)
        for ch = 1:size(statData,2)
            if any(abs(statData(:,ch,clz)) > valThr)
               outliers = [outliers, ch];
            end
        end
    end
end