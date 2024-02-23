function bVals = binningVals(vals, binRangeLabels)
% given a column vector 'vals' and binRangeLabels with format:
% binRangeLabels =
% 
%    -1.0000   -0.7854    0.7854
%     1.0000   -3.1416   -2.3562
%     1.0000    2.3562    3.1416
%     label    loRange   hiRange
% bins the vector into N unique bins and returns the bin numbers

% (c) Jiri, Apr11

%% calling function (ex: dir. discretization)
% % x-direction discretization
% binRangeLabels = [-1, -pi/4, pi/4;
%                   +1, -pi, -3/4*pi;
%                   +1, +3/4*pi, +pi];
% xDir = binningVals(dataSession(useInd{sess},i_dir), binRangeLabels);
% % y-direction discretization
% binRangeLabels = [-1, -3/4*pi, -pi/4;
%                   +1, pi/4, 3/4*pi];
% yDir = binningVals(dataSession(useInd{sess},i_dir), binRangeLabels);
    
    
%% init
assert(size(vals,2) == 1);
bVals = nan(size(vals,1), size(vals,2));

%% TO DO ... bin it
nLabels = unique(binRangeLabels(:,1));
for r = 1:size(binRangeLabels ,1)
end