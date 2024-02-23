function triggers = getTriggers_zeroSpeed(params)
% detects events of zero speed (pauses) as continuous data payches whose
% |velocity| < thr [std]
% pauses -> trials
% binned pauses -> classes

% (c) Jiri, Dec12

%% try to load the processed dataset
if ismember('triggers', who('-file', params.storage.cacheFile))  
    display('Triggers found.');
    load(params.storage.cacheFile, 'triggers');
    return;
end

%% settings
assert(isfield(params, 'zeroSpeed'));
predCh = 1;
position = 'absVel';
signalType = 'Traj';
selCh = getDecodingChannels(params.connectionTable, signalType, position);
info.selCh = selCh;                                                     % expects non-empty row vector of selected channels
info.rejCh = [];
info.name = signalType;
info.chPos = position;   
info.nGroups = 1;                                                       % TO DO: for mulitple groups...    

%% detect zeroSpeed patches & their lengths
nSess = length(params.simulation.trainSession);
pred = cell(1,nSess);
sessSizes = nan(nSess,1);
i_trigs = [];                                                           % indices  of triggers 
t_trigs = [];                                                           % time  of triggers, in [s]
d_trigs = [];                                                           % duration of triggers, in [samples]
i_sess  = [];                                                           % session indices
for sess = 1:nSess

    % detect triggers & their duration
    thisSess = params.simulation.trainSession(sess);
%     normSpeed = getStdSpeed(params, thisSess, 'absVel');
    params.thisSess = thisSess;
    [filterSettings ~] = getFeatureMethod(params, {'raw'});                 % to reject bad recordings epochs
    normSpeed = filterData(params, getStdSpeed(params, thisSess, 'absVel'), filterSettings);
    
    sessTime = [1:size(normSpeed,1)]'./params.amp.srate;
    sessSizes(sess) = size(normSpeed,1);
    pred{sess} = normSpeed; 
    i_zero = find(normSpeed < params.zeroSpeed.threshold);              % all in [std]
    binarySpeed = zeros(size(normSpeed,1),1);
    binarySpeed(i_zero) = 1;
    [i_event, t_event] = getContValsPatches(binarySpeed, 1, 'beg');     % beg & duration of pause, in [samples]

    % pool over sessions
    i_trigs = cat(1, i_trigs, i_event);                                 % event sample index 
    t_trigs = cat(1, t_trigs, sessTime(i_event));                       % time of event, in [s]
    d_trigs = cat(1, d_trigs, t_event);                                 % event duration, in [samples]
    i_sess  = cat(1, i_sess, sess*ones(size(i_event,1),1));
end

%% bin pause durations
yVals = d_trigs./params.amp.srate;                                      % event duration, in [s]
tunInfo = struct;
tunInfo.minV = min(yVals);
tunInfo.maxV = max(yVals);
tunInfo.nBins = params.zeroSpeed.nBins;
tunInfo.binMethod = params.zeroSpeed.binMethod;
if strcmp(tunInfo.binMethod,'defBinSizes'), tunInfo.defBinSizes = params.binSizes.ampDefinition; end
tunInfo.binSizes = [];                        % empty -> bin sizes get computed (if needed)
[~, tunInfo.binSizes] = binVector(tunInfo, yVals);   
[binsArray, ~] = binVector(tunInfo, yVals);
assert(size(binsArray,1) == size(yVals,1));
tunInfo.binCenters = findBinCenters(binsArray, yVals);
binValues = sort(unique(binsArray));    

%% output
triggers{predCh}.sampleSessLabel = cat(2, i_trigs, i_sess, binsArray);     % mandatory field, triplets = [sample, sess, classLabel]
triggers{predCh}.timeSessLabel = cat(2, t_trigs, i_sess, binsArray);       % incl. time of event, in [s]
triggers{predCh}.pauseDuration = yVals;                                    % in [s]
triggers{predCh}.eventDuration = yVals;                                    % mandatory field, in [s]
triggers{predCh}.info = info;
triggers{predCh}.info.binning = tunInfo;
triggers{predCh}.info.sessSizes = sessSizes;                                                        % mandatory field
triggers{predCh}.info.dataType = 'timeSeries';                                                      % mandatory field
triggers{predCh}.info.clzColors = colorPalette(length(binValues));

save(params.storage.cacheFile, 'triggers', 'd_trigs', 'tunInfo', '-append');



%% TO DO: plot detected classes, if selected
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
        plot(t, pred{sess}, 'r');                                                                   % in [std]
        plot([min(t), max(t)], [0, 0], '--k');
        plot([min(t), max(t)], [params.zeroSpeed.threshold, params.zeroSpeed.threshold], ':k');     % in [std]
        i_sess = find(triggers{predCh}.sampleSessLabel(:,2) == sess);
        trig_sampleClass =  triggers{predCh}.sampleSessLabel(i_sess,:);                 % ~ rows of: [sample, sess, classLabel]
        trig_eventDuration =triggers{predCh}.eventDuration(i_sess)*srate;               % in [samples] 
        for clz = binValues'
            i_clz = find(trig_sampleClass(:,3) == clz);
            clr = triggers{predCh}.info.clzColors(clz,:);
            for n = 1:length(i_clz)
                trigBeg = trig_sampleClass(i_clz(n),1);                                 % in [samples]
                trigEnd = trigBeg + trig_eventDuration(i_clz(n));                       % in [samples]
                plot([trigBeg, trigBeg]./srate, get(gca,'ylim'), 'Color',clr);          % plot triggers begin
                plot([trigEnd, trigEnd]./srate, get(gca,'ylim'), 'Color',clr);          % plot triggers end
            end   
        end
%         for n = 1:length(tunInfo.binSizes)
%             binS = tunInfo.binSizes(n);
%             clrs = cat(1, triggers{predCh}.info.clzColors(1,:),triggers{predCh}.info.clzColors);
%             plot([min(t), max(t)], [binS, binS], 'Color',clrs(n,:));
%         end
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
    hist(yVals, 5000);
    for n = 1:length(tunInfo.binSizes)                  % bin boundaries
        binS = tunInfo.binSizes(n);
        plot([binS, binS], get(gca,'ylim'), '--k');
    end        
    binCenters = tunInfo.binCenters;                    % bin centers (medians)
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
