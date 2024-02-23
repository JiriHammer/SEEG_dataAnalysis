function trialsData = prepare2plotTrials(params, trials)
% arranges the 'trials' structure to trialsData structure used by plotTrials.m

% (c) Jiri, Nov16

if ~isfield(trials, 'get_pValsBase')
    trials.get_pValsBase = true;
end
if ~isfield(trials, 'get_pVals')
    trials.get_pVals = true;
end
if ~isfield(params, 'save_trialsData')
    params.save_trialsData = true;
end
if ~isfield(params.triggering, 'baselineCorrection')
    params.triggering.baselineCorrection = true;
end

%% auxiliary info
trialsData = struct;
trialsData.info.verLines = params.plot_triggering.verLines;
trialsData.info.horLines = params.plot_triggering.horLines;
trialsData.info.xlabel = 'time (s)';

%% y-limits for each channel (if commented out, the program estimates common y-lims for all channels)
%trialsData.info.yLims = [-0.15, 0.15];      % defined common y-limits for all channels
%trialsData.info.yLims = [];                 % each channel has its own y-limits

%% storage info
if isfield(trials, 'outDir')
    trialsData.info.outDir = trials.outDir;
else
    %trialsData.info.outDir = [params.storage.outputDir filesep 'classMeansPlots'];  % default
    trialsData.info.outDir = [params.storage.dir_results filesep params.storage.outName];
end
if isfield(trials, 'figname')
    trialsData.info.figName = trials.figname;
else
    trialsData.info.figName = 'noNameFig';  % default
end

%% text on figure
if isfield(trials, 'selFreq')
    trialsData.info.text = ['MEAN +/- SEM: subject: ' params.storage.subjTag ...
        ', frequency band: ' num2str(trials.selFreq(1)) ' - ' num2str(trials.selFreq(2)) ...
        ', triggered: ' params.triggering.cutPoint ...
        ', classes: '];
else
    trialsData.info.text = ['MEAN +/- SEM: subject: ' params.storage.subjTag ...
        ', ERP' ...
        ', triggered: ' params.triggering.cutPoint ...
        ', classes: '];
end

%% downsample? (for ERP, to speed up significance computations)
if mean(diff(trials.time),1) < params.stft_freq.timeStep
    clear srate;
    load(params.storage.cacheFile, 'srate');
    timeStep = ceil(params.stft_freq.timeStep * srate);    % in [samples], if any downsampling occurs, 'timeStep' MUST exist
    trials.time = downsample(trials.time, timeStep);
    % decimate by timeStep
    decData = nan(size(trials.time,1),size(trials.data,2),size(trials.data,3));
    for ch = 1:size(trials.data,2)
        for tr = 1:size(trials.data,3)
            decData(:,ch,tr) = decimate(trials.data(:,ch,tr), timeStep);
        end
    end           
    trials.data = decData;
    clear decData;
%     trials.data = downsample(trials.data, timeStep);    % ??? or decimate ???
    if isfield(trials, 'rejected')
        trials.rejected = downsample(trials.rejected, timeStep);
    end
    if isfield(trials, 'sessTime')
        trials.sessTime = downsample(trials.sessTime, timeStep);
    end
end
    
%% selected time range to plot
%i_t = closestval(trials.time,params.plot_triggering.time2plot(1)):closestval(trials.time,params.plot_triggering.time2plot(2));

%% data dimesions
nSamples = size(trials.data,1);
nCh = size(trials.data,2);

%% values: class means +/- SEM
clzLabels = unique(trials.labels);
nClz = length(clzLabels);
trialsData.xVals = trials.time;
trialsData.yVals = nan(nSamples, nCh, nClz);
trialsData.yErrs = nan(nSamples, nCh, nClz);
trialsData.yVars = nan(nSamples, nCh, nClz);
trialsData.info.time2plot = params.plot_triggering.time2plot;           % selected time to plot
trialsData.info.colors = colorPalette(size(trials.clzNames,1));         % should contain all defined classes
trialsData.info.outliers = zeros(nCh, nClz);                            % uncomment to auto-detect outliers
trialsData.info.nClz = [];
trialsData.info.clzNames = cell(nClz,1);
trialsData.info.clzNumber = nan(nClz,1);
for clz = 1:nClz
    i_clz = find(trials.labels == clzLabels(clz));                      % selected trials
    if isfield(trials, 'rejected')
        inds = ~trials.rejected;                                            % not rejected indices = 1, rejected indices = 0
%         i_clz = find(trials.labels == clzLabels(clz));                      % selected trials
        i_otherTrials = setdiff(1:size(trials.data,3),i_clz);               % all other trials
        inds(:,:,i_otherTrials) = false;                                    % not selected trials are set to 0 (i.e. rejected)
        vals = trials.data;
        vals(inds == 0) = NaN;                                              % rejected indices are set to NaN
    else
        vals = trials.data(:,:,i_clz);         % no rejection
    end
    
    % baseline correction for each class (again, spectral power seemed shifted)
    if params.triggering.baselineCorrection
        assert(~isempty(params.triggering.baseline));
        i_b = closestval(trials.time,params.triggering.baseline(1)):closestval(trials.time,params.triggering.baseline(2));
        base = nanmean(nanmean(vals(i_b,:,i_clz),1),3);   % 1 x ch, baseline (mean over time & clz)
        vals(:,:,i_clz) = vals(:,:,i_clz) - repmat(base, [size(vals,1),1,size(i_clz,2)]);   % baseline subtraction
        trials.data(:,:,i_clz) = trials.data(:,:,i_clz) - repmat(base, [size(vals,1),1,size(i_clz,2)]);   % baseline subtraction
    end
    
    % >>> AVG +/- SEM <<<
    trialsData.yVals(:,:,clz) = nanmean(vals,3);                        % ~ mean (average) over non-rejected trials
    trialsData.yErrs(:,:,clz) = nanstd(vals,0,3)./length(i_clz)^0.5;    % ~ standard error of the mean (SE or SEM)
    trialsData.yVars(:,:,clz) = nanvar(vals,0,3);                       % ~ standard error of the mean (SE or SEM)

    clr = trialsData.info.colors(clz,:);
    [tf, i_c] = ismember(clzLabels(clz),cell2mat(trials.clzNames(:,2)));    % i_c = pointer to all defined classes !!! (some may not be found)
    assert(tf);   
    trialsData.info.nClz = [trialsData.info.nClz, length(i_clz)];
    trialsData.info.text = [trialsData.info.text ...
        '\color[rgb]{' num2str(clr(1)) ' ' num2str(clr(2)) ' ' num2str(clr(3)) '}' trials.clzNames{i_c,1} '(' num2str(length(i_clz)) ') '];
    trialsData.info.clzNames{clz,1} = trials.clzNames{i_c,1};               % class name from all defined classes (some may not be found)
    trialsData.info.clzNumber(clz) = i_c;                                 % i_c = pointer to all defined classes !!! (some may not be found)
end

%% test for significance of difference in distribution between 2 classes
if trials.get_pVals
    [pVals, hVals] = getSignificance(trials.data, trials.labels);
    trialsData.pVals = pVals;
    trialsData.hVals = hVals;
end

%% test for significance against baseline
if trials.get_pValsBase
    [pVals_base, hVals_base] = getSignificance_trials_baseline(trials.data, trials.time, trials.labels, params.triggering.baseline);
    trialsData.pVals_base = pVals_base;
    trialsData.hVals_base = hVals_base;
end

%% channel names
params.nCh = size(trialsData.yVals,2);
trialsData.info.chNames = getChannelNames(params, lower(params.response.signalType), 'resp');

%% paradigm times to mark (reaction times, ...)
if isfield(trials, 'info')
    if isfield(trials.info, 'paraTimes')
        trialsData.info.paraTimes = trials.info.paraTimes;  % pre-computed
    end
else
    trialsData.info.paraTimes = getParadigmTimes(params);   % compute
end

%% save trialsData to cache file
if params.save_trialsData
    outputVarName = ['trialsData_' trials.freqBandName];
    eval([outputVarName '=trialsData;' ]); 
    disp(['Saving trials data: ' trials.freqBandName ' ...']);
    save(params.storage.cacheFile, outputVarName, '-append');  
end    

