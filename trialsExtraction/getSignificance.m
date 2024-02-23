function [pVals, hVals] = getSignificance(trials, labels, alphaVal)
% determines if data in 'trials' are statistically significantly different
% based on their 'labels'
% trials: 3D = samples x channels x trials
% labels: 1D = trials labels

% (c) Jiri, Mar16

if nargin < 3
    alphaVal = 0.01;        % significance level
end

clzs = unique(labels);
if numel(clzs) == 2
    disp(['Significance testing (Wilcoxon ranksum test + FDR), significance level = ' num2str(alphaVal) ' ...']);
else
    disp(['WARNING: cannot test for significance. Number of classes = ' num2str(numel(clzs)) '.']);
    disp(' - comparison possible only for 2 classes!');
    pVals = nan(size(trials,1),size(trials,2));
    hVals = nan(size(trials,1),size(trials,2));
    return;
end
assert(size(trials,3) == length(labels));

%% go over the trials: 1st & 2nd dim (samples + channels)
pVals = nan(size(trials,1), size(trials,2));
for ch = 1:size(trials,2)
    for t = 1:size(trials,1)
        tmp_data = cell(1,length(clzs));
        for clz = 1:length(clzs)
            i_clz = find(labels == clzs(clz));            
            tmp = trials(t,ch,i_clz);
            tmp_data{1,clz} = tmp(:);
        end
        pVals(t,ch) = ranksum(tmp_data{1}, tmp_data{2}, 'alpha',alphaVal);
    end
    if mod(ch,10) == 0
        disp([' - channel: ' num2str(ch) ' done.']);
    end    
end

%% FDR correction (False Discovery Rate)
pVals_tmp = pVals(:);
[n_signif,i_signifTest] = fdr(pVals_tmp, alphaVal, 'original');

%% significant values (1 = yes, 0 = no)
hVals_tmp = zeros(size(pVals_tmp,1),1);
if n_signif > 0
    hVals_tmp(i_signifTest) = 1;
end
hVals = reshape(hVals_tmp, size(trials,1), size(trials,2));


