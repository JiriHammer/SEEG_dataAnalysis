function t_SA = getSignificantActivationTimes(hVals, pVals, xVals, yVals, SA_params)
% t_SA = cell(1,2) of significant activation times
% t_SA{1} = class 1 > class 2
% t_SA{2} = class 1 < class 2
% 
% hVals = 1D: lag x 1
% pVals = 1D: lag x 1
% yVals = 2D: lag x clz, where nClz = 2 !!! works only for two classes !!!

% (c) Jiri, Feb19

method = SA_params.whichActivationTimes;
assert(size(hVals,1) == size(pVals,1));
assert(size(hVals,1) == size(yVals,1));
assert(size(hVals,1) == size(xVals,1));
assert(size(yVals,2) == 2); % !!! works only for two classes !!!

%% extract indices of significant activation times -> i_SA
if strcmp(method, 'first_significant')
    i_SA = find(hVals == 1, 1, 'first');
    
elseif strcmp(method, 'max_significant')
    i_all = find(hVals == 1);
    [m, i_SA] = min(pVals);
    assert(ismember(i_SA, i_all));
    
elseif strcmp(method, 'all_significant')
    i_SA = find(hVals == 1);
    
else
    error(['unknown significant activation method: ' method]);
end

%% assign significant activation times
t_SA = cell(1,2);
for lag = i_SA'
    if strcmp(SA_params.signifValsNames{1},'pVals')
        if yVals(lag,1) > yVals(lag,2)                                
            t_SA{1} = cat(1, t_SA{1}, xVals(lag));
        else
            t_SA{2} = cat(1, t_SA{2}, xVals(lag));
        end
    elseif strcmp(SA_params.signifValsNames{1},'pVals_base')
        c = SA_params.thisClz;
        if yVals(lag,c) > 0             % assumes that baseline = 0 !!!                                
            t_SA{1} = cat(1, t_SA{1}, xVals(lag));
        else
            t_SA{2} = cat(1, t_SA{2}, xVals(lag));
        end
    end
        
end

%% debug figure (uncomment)
% figure; hold on;
% clrs = {'b','r'};
% clz = 1;
% plot(xVals, yVals(:,clz), clrs{clz});
% clz = 2;
% plot(xVals, yVals(:,clz), clrs{clz});
% for clz = 1:size(t_SA,2)
%     plot(t_SA{clz}, 0, '*', 'Color', clrs{clz});
% end

