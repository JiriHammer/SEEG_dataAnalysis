function rejectData(params, varName_prefix)
% rejects outliers in data based on normal distribution 

% (c) Jiri, Apr16

% load data
varName = [varName_prefix '_prc'];
load(params.storage.sessionCacheFiles{params.thisSess}, varName);
d = eval(varName);
r = logical(zeros(size(d)));

% normal distrubution fit
sigmaThr = params.rejection.normalDistr.sigmaThr;
for ch = 1:size(d,2)
    [muhat,sigmahat] = normfit(d(:,ch));
    i_rej = d(:,ch) > sigmaThr*sigmahat;
    r(i_rej,ch) = 1;
end
r = logical(r);

% save
outputVarName = [varName_prefix '_rej'];
eval([outputVarName ' = r;' ]);
save(params.storage.sessionCacheFiles{params.thisSess}, outputVarName, '-append');
