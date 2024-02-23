function clims = getColorLimits(threeDimArray)
% returns color limits from 'threeDimArray'

% (c) Jiri, Jan11

yVals =[];
for ch1 = 1:size(threeDimArray)
    for ch2 = 1:size(threeDimArray,2)
        for ch3 = 1:size(threeDimArray,3)
            yVals = cat(1, yVals, threeDimArray(ch1,ch2,ch3));
        end
    end
end

absMaxVal = max([abs(min(yVals)), abs(max(yVals))]);
clims = [-absMaxVal, absMaxVal];