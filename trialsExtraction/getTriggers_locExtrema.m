function triggers = getTriggers_locExtrema(params)
% triggers on local extrema of predictor variable, defined in:
% params.locExtrema = struct(...
%     'predictor', params.response, ...
%     'minOrMax', 'both', ...                      % options: min, max, both
%     'deadPeriod', 1, ...                        % in [s], TO DO...
%     'nBins', 5 ...
%     );
% returns cell-struct: triggers{1}.sampleSessLabel

% (c) Jiri, Mar12

assert(isfield(params, 'locExtrema'));
sampleSessVals = [];

%% load predictor signal
loadOrFilter(params, params.locExtrema.predictor);
nSess = length(params.simulation.trainSession);
sessSizes = nan(nSess,1);
pred = cell(nSess,1);
for sess = 1:nSess

    % load/select predictor channel
    thisSess = params.simulation.trainSession(sess);
    params.signalType = params.locExtrema.predictor.signalType;
    params.car.rejCh = union(params.locExtrema.predictor.excludedChnls, params.car.rejCh);
    [filterSettings catFilterNames] = getFeatureMethod(params, params.locExtrema.predictor.dataProcessing);
    inputVarName = [catFilterNames params.locExtrema.predictor.signalType];
    load(params.storage.sessionCacheFiles{thisSess}, inputVarName);  % try to load the processed dataset
    assert(exist(inputVarName, 'var') == 1);
    data = eval(inputVarName);
    clear(inputVarName);
    % select channel
    predCh = 1;
    position = params.locExtrema.predictor.selectedChnls{predCh};               % TO DO: for multiple predictors...
    selCh = getDecodingChannels(params.connectionTable, params.locExtrema.predictor.signalType, position);
    info.selCh = setdiff(selCh, params.locExtrema.predictor.excludedChnls);    % expects non-empty row vector of selected channels
    info.rejCh = params.locExtrema.predictor.excludedChnls;
    info.name = params.locExtrema.predictor.signalType;
    info.chPos = position;   
    info.nGroups = 1;                                                           % TO DO: for mulitple groups...
    if info.nGroups == 1
        pred{sess} = data(:, info.selCh);
    else
        pred{sess} = squeeze(data(:, info.thisGroup, info.selCh));
    end    
    clear data;   
    sessSizes(sess) = size(pred{sess},1);

    % extract local extrema
    if strcmp(params.locExtrema.minOrMax, 'both')
        sampleValsMax = getLocalMaxima(params, pred{sess});
        sampleValsMin = getLocalMaxima(params, -pred{sess});
        sampleValsMin(:,2) = -sampleValsMin(:,2);
        sampleVals = cat(1, sampleValsMax, sampleValsMin);
    elseif strcmp(params.locExtrema.minOrMax, 'max')
        sampleVals = getLocalMaxima(params, pred{sess});
    elseif strcmp(params.locExtrema.minOrMax, 'min')
        sampleVals = getLocalMaxima(params, -pred{sess});
        sampleVals(:,2) = -sampleVals(:,2);
    end
    sampleSessVals = cat(1, sampleSessVals, cat(2, sampleVals(:,1), sess*ones(size(sampleVals,1),1), sampleVals(:,2)));
end

%% bin extrema values into nBins
yVals = sampleSessVals(:,3);
tunInfo = struct;
tunInfo.minV = min(yVals);
tunInfo.maxV = max(yVals);
tunInfo.nBins = params.locExtrema.nBins;
tunInfo.binMethod = 'constSamplesPerBin';
tunInfo.binSizes = [];                        % empty -> bin sizes get computed (if needed)
[~, tunInfo.binSizes] = binVector(tunInfo, yVals);   
[binsArray, ~] = binVector(tunInfo, yVals);
assert(size(binsArray,1) == size(yVals,1));
binValues = sort(unique(binsArray));    

%% output
triggers{predCh}.sampleSessLabel = cat(2, sampleSessVals(:,1), sampleSessVals(:,2), binsArray);     % mandatory field
triggers{predCh}.info = info;
triggers{predCh}.info.binning = tunInfo;
triggers{predCh}.info.sessSizes = sessSizes;                                                        % mandatory field
triggers{predCh}.info.dataType = 'timeSeries';                                                      % mandatory field
triggers{predCh}.info.clzColors = colorPalette(length(binValues));

%% plot detected classes, if selected
if params.triggering.plotDetection
    % define sampling rate
    load(params.storage.cacheFile, 'sampleRateAfterFiltering');
    assert(exist('sampleRateAfterFiltering','var') == 1);
    srate = sampleRateAfterFiltering;        
    for sess = 1:nSess                                                      % !!! sess-by-sess time series !!!
        t = [1:size(pred{sess},1)]./srate;
        f = figure('visible', 'on', 'Position', [1, 1, 1920, 1200]);
        set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);
        hold on;
        plot(t, pred{sess}, 'r');
        plot([min(t), max(t)], [0, 0], '--k');
        i_sess = find(triggers{predCh}.sampleSessLabel(:,2) == sess);
        trig_sampleClass = triggers{predCh}.sampleSessLabel(i_sess,:);
        for clz = binValues'
            i_clz = find(trig_sampleClass(:,3) == clz);
            clr = triggers{predCh}.info.clzColors(clz,:);
            for n = 1:length(i_clz)
                trig = trig_sampleClass(i_clz(n),1);
                plot([trig, trig]./srate, get(gca,'ylim'), 'Color',clr);        % triggers
            end   
        end
        for n = 1:length(tunInfo.binSizes)
            binS = tunInfo.binSizes(n);
            clrs = cat(1, triggers{predCh}.info.clzColors(1,:),triggers{predCh}.info.clzColors);
            plot([min(t), max(t)], [binS, binS], 'Color',clrs(n,:));
        end
        xlabel('time [s]');
        % output directory & save
        outDir = [params.storage.outputDir filesep 'triggerDetection'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end  
        figname = [info.chPos '_timeSeries_sess' num2str(sess)];
        saveas(f, [outDir filesep figname '.fig']);
        close(f);                    
    end
    f = figure('visible', 'on', 'Position', [1, 1, 1920, 1200]);                % !!! histogram of trigger values !!!
    set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);
    hold on;    
    hist(yVals, 500);
    for n = 1:length(tunInfo.binSizes)
        binS = tunInfo.binSizes(n);
        plot([binS, binS], get(gca,'ylim'), '--k');
    end        
    binCenters = getBinCenters(tunInfo.binSizes);
    for n = 1:length(binCenters)
        binS = binCenters(n);
        clrs = triggers{predCh}.info.clzColors;
        plot([binS, binS], get(gca,'ylim'),'Color',clrs(n,:), 'LineWidth',2);
    end 
    figname = [info.chPos '_histWithBinCenters'];
    saveas(f, [outDir filesep figname '.fig']);
    print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname]);        
    close(f);         
end
