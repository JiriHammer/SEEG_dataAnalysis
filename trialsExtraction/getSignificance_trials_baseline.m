function [pVals, hVals] = getSignificance_trials_baseline(trials, time, labels, baseline_period)
% determines if data in 'trials' are statistically significantly different from baseline
% based on their 'labels'
% trials: 3D = samples x channels x trials
% time: 1D = samples
% labels: 1D = trials labels
% pVals, hVals: : 3D = samples x channels x nClz

% (c) Jiri, Feb19

%% significance level
alphaVal = 0.01;        
disp(['Significance against baseline (Wilcoxon ranksum test + FDR), significance level = ' num2str(alphaVal) ' ...']);

%% dim checks
clzs = unique(labels);
assert(size(trials,3) == length(labels));
assert(size(trials,1) == size(time,1));

%% define baseline period
i_base = closestval(time,baseline_period(1)):closestval(time,baseline_period(2));

%% go over the classes & time points
pVals = nan(size(trials,1), size(trials,2), numel(clzs));
for clz = 1:numel(clzs)
    disp([' - class: ' num2str(clz) ' ...']);
    i_clz = find(labels == clzs(clz));            
    for ch = 1:size(trials,2)
        data_base = trials(i_base,ch,i_clz);
        for t = 1:size(trials,1)
            data_chnl = squeeze(trials(t,ch,i_clz));
            pVals(t,ch,clz) = ranksum(data_chnl, data_base(:), 'alpha',alphaVal);
        end
        if mod(ch,10) == 0
            disp([' -- channel: ' num2str(ch) ' done.']);
        end    
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
hVals = reshape(hVals_tmp, [size(trials,1), size(trials,2), numel(clzs)]);


