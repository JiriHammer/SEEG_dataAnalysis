function outliers = getOutliers(statData)
% statData: 3D = samples x channels x trials 
% outliers: if a value of a channel is out at 6 sigma of the normal
% disctribution of all channels, the channel is marked as an outlier
% outliers: 2D = channels x trials

% (c) Jiri, Mar16

% threshold
stdThr = 6;     % in [SD]

% fit normal distribution
x = statData(:);
%figure; histfit(x);
[muhat,sigmahat] = normfit(x);

% outliers
outliers = zeros(size(statData,2), size(statData,3));
valThr = muhat + stdThr * sigmahat;
for ch = 1:size(statData,2)
    for clz = size(statData,3)
        if any(abs(statData(:,ch,clz)) > valThr)
           outliers(ch,clz) = 1;
        end
    end
end
