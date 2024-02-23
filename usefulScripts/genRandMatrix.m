function r = genRandMatrix(vals, nRows, nCols)
% given vector of values vals, generates a matrix with elements randomly
% chosen from vals

% (c) Jiri, Apr11

% indexes
inds = 1:length(vals);

% matrix
r = nan(nRows,nCols);
for c = 1:nCols
    r(:,c) = vals(round(linTransform(rand(nRows,1), [0 1], [min(inds)-0.5, max(inds)+0.5])));
end

%r = vals(round(linTransform(rand(nCols, nRows), [0 1], [min(inds)-0.5, max(inds)+0.5])));