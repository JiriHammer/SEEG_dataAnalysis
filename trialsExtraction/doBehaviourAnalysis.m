function doBehaviourAnalysis(params)
% analyzes reaction times or error rates (right/wrong/missed) trials
% loads selected trials from: pred_prc (see: getTriggers_..._timesLabels.m)

% (c) Jiri, Jun20

if ~isfield(params.triggering, 'behaviour')
    return;
end
disp(' - triggering: behavioural analysis ...');

%% reaction times 
if ismember('RT', params.triggering.behaviour)
    disp(' - triggering: behavioural analysis: reaction times ...');
    
    nClz = size(params.triggering.classes,1);    % number of all classes (some may not be found)  
    RT_all = [];
    for sess = 1:size(params.storage.sessionCacheFiles,2)
        
        % load triggers
        clear triggers classNamesLabels;
        load(params.storage.sessionCacheFiles{sess}, 'triggers', 'classNamesLabels');
        if ~isempty(triggers)
            %assert(nClz == size(triggers,1));   % triggers = 2D: row1 = RT; row2 = clzLabel
            %assert(size(pred_prc,1) == size(triggers,2));

            % load predictor (trials): pred_prc (see: getTriggers_..._timesLabels.m)
            clear pred_prc
            load(params.storage.sessionCacheFiles{sess}, 'pred_prc'); 
        
            % extract reaction times from selected trials
            for tr = 1:size(pred_prc,1)
                if pred_prc(tr).wasSelected
                    if isfield(pred_prc(tr), 'reactionTime')
                        if ~isnan(pred_prc(tr).reactionTime)
                            RT_all = cat(2, RT_all, [pred_prc(tr).reactionTime; pred_prc(tr).clzLabel]);
                        end
                    end
                end
            end
        end
    end
    if isempty(RT_all)
        disp(' - no trials with reaction times found.');
        return;
    end
    
    % separate classes
    y_sem = [];
    y_avg = [];
    RT_clz = cell(nClz,1);
    for clz = 1:nClz
        i_clz = find(RT_all(2,:) == classNamesLabels{clz,2});
        RT_clz{clz} = RT_all(1,i_clz);
        y_avg = cat(2, y_avg, mean(RT_all(1,i_clz),2));
        y_sem = cat(2, y_sem, sem(RT_all(1,i_clz),2));
    end
    x = 1:size(y_avg,2);
    
    % save mean RTs to cacheFile
    RT = cat(2, y_avg', y_sem');
    save(params.storage.cacheFile, 'RT', 'RT_all', '-append');
    
    % significance test: hardcoded for 2 classes
    if ~isempty(RT_clz{1}) && ~isempty(RT_clz{2})
        P = ranksum(RT_clz{1},RT_clz{2});
    else
        disp('Empty class encountered, returning.');
        return;
    end
    
    % figure
    f = fig_make;
    subplot(1,1,1); hold on;
    barwitherr(y_sem,y_avg);
    set(gca, 'xlim', [-2, nClz+2]);
    set(gca,'XTick',x);
%     set(gca,'XTickLabel',classNamesLabels(:,1), 'FontSize',12); 
    set(gca,'XTickLabel',strrep(classNamesLabels(:,1), '_', '\_'), 'FontSize',12); 
    ylabel('reaction time (s)');
    title('Reaction times.');
    % TO DO: makeSignificanceBar_jiri(x,y,p,thr)
    outDir = [params.storage.outputDir filesep 'behaviour'];
    fig_save(f, 'RT', outDir);
    close(f);
end

%% answer ok: right-wrong-missed 
if ismember('RWM', params.triggering.behaviour)
    disp(' - triggering: behavioural analysis: answers: right-wrong-missed ...');
    
    nClz = size(params.triggering.classes,1);    % number of all classes (some may not be found)  
    RWM_all = zeros(nClz,4);            % = 4 because: right(=1), wrong(=2), missed(=3),noClue(=4)
    for sess = 1:size(params.storage.sessionCacheFiles,2)
        
        % load triggers
        clear triggers classNamesLabels;
        load(params.storage.sessionCacheFiles{sess}, 'triggers', 'classNamesLabels');
        if ~isempty(triggers)
            %assert(nClz == size(triggers,1));   % triggers = 2D: row1 = RT; row2 = clzLabel
            %assert(size(pred_prc,1) == size(triggers,2));

            % load predictor (trials): pred_prc (see: getTriggers_..._timesLabels.m)
            clear pred_prc
            load(params.storage.sessionCacheFiles{sess}, 'pred_prc');         

            % extract reaction times from selected trials
            for tr = 1:size(pred_prc,1)
                if pred_prc(tr).wasSelected
                    answer_given = [];
                    if isfield(pred_prc(tr), 'answer_ok')
                        answer_given = pred_prc(tr).answer_ok;
                    elseif isfield(pred_prc(tr), 'answer')
                        answer_given = pred_prc(tr).answer;
                    end
                    if ~isempty(answer_given)
                        if strcmp(answer_given, 'right')
                            RWM_all(pred_prc(tr).clzLabel,1) = RWM_all(pred_prc(tr).clzLabel,1)+1;
                        elseif strcmp(answer_given, 'wrong')
                            RWM_all(pred_prc(tr).clzLabel,2) = RWM_all(pred_prc(tr).clzLabel,2)+1;
                        elseif strcmp(answer_given, 'missed')
                            RWM_all(pred_prc(tr).clzLabel,3) = RWM_all(pred_prc(tr).clzLabel,3)+1;
                        elseif strcmp(answer_given, 'noClue')
                            RWM_all(pred_prc(tr).clzLabel,4) = RWM_all(pred_prc(tr).clzLabel,4)+1;    
                        else
                            disp(['WARNING:  unknown name: ' answer_given ', in sess = ' num2str(sess) ', trial = ' num2str(tr)]);
                        end
                    end
                end
            end
        end
    end
    if isempty(RWM_all)
        disp(' - no trials with right/wrong/missed responses found.');
        return;
    end
    
    % save mean RTs to cacheFile
    RWM = RWM_all;
    save(params.storage.cacheFile, 'RWM', '-append');
    
    % figure
    f = fig_make;
    subplot(1,1,1); hold on;
    bar(RWM_all);
    x = 1:size(RWM_all,1);
    set(gca, 'xlim', [0, nClz+1]);
    set(gca,'XTick',x);
%     set(gca,'XTickLabel',classNamesLabels(:,1), 'FontSize',12);  
    set(gca,'XTickLabel',strrep(classNamesLabels(:,1), '_', '\_'), 'FontSize',12);
    ylabel('N trials');
    title('answers: right-wrong-missed');
    legend({'right','wrong','missed','noClue'});
    outDir = [params.storage.outputDir filesep 'behaviour'];
    fig_save(f, 'RWM', outDir);
    close(f);
end

%% button pressed: red - green
if ismember('BP_RedGrn', params.triggering.behaviour)
    disp(' - triggering: behavioural analysis: answers: yes-no buttons ...');
    
    nClz = size(params.triggering.classes,1);    % number of all classes (some may not be found)  
    RWM_all = zeros(nClz,2);
    for sess = 1:size(params.storage.sessionCacheFiles,2)
        
        % load triggers
        clear triggers classNamesLabels;
        load(params.storage.sessionCacheFiles{sess}, 'triggers', 'classNamesLabels');
        if ~isempty(triggers)
            %assert(nClz == size(triggers,1));   % triggers = 2D: row1 = RT; row2 = clzLabel
            %assert(size(pred_prc,1) == size(triggers,2));
           
            % load predictor (trials): pred_prc (see: getTriggers_..._timesLabels.m)
            clear pred_prc
            load(params.storage.sessionCacheFiles{sess}, 'pred_prc'); 
        
            % extract reaction times from selected trials
            for tr = 1:size(pred_prc,1)
                if pred_prc(tr).wasSelected
                    if isfield(pred_prc(tr), 'button_pressed')
                        if ~isempty(pred_prc(tr).button_pressed)
                            if strcmp(pred_prc(tr).button_pressed, 'red') || strcmp(pred_prc(tr).button_pressed, 'L')
                                RWM_all(pred_prc(tr).clzLabel,1) = RWM_all(pred_prc(tr).clzLabel,1)+1;
                            elseif strcmp(pred_prc(tr).button_pressed, 'green') || strcmp(pred_prc(tr).button_pressed, 'R')
                                RWM_all(pred_prc(tr).clzLabel,2) = RWM_all(pred_prc(tr).clzLabel,2)+1;
                            else
                                disp(['WARNING:  unknown name: ' pred_prc(tr).button_pressed ', in sess = ' num2str(sess) ', trial = ' num2str(tr)]);
                            end
                        end
                    end
                end
            end
        end
    end
    if isempty(RWM_all)
        disp(' - no trials with red/green buttons press responses found.');
        return;
    end
    
    % save mean RTs to cacheFile
    YNB = RWM_all;
    save(params.storage.cacheFile, 'YNB', '-append');
    
    % figure
    f = fig_make;
    subplot(1,1,1); hold on;
    bar(RWM_all);
    x = 1:size(RWM_all,1);
    set(gca, 'xlim', [0, nClz+1]);
    set(gca,'XTick',x);
    set(gca,'XTickLabel',classNamesLabels(:,1), 'FontSize',12);  
    title('answers: green-red buttons pressed');
    legend({'red','green'});
    outDir = [params.storage.outputDir filesep 'behaviour'];
    fig_save(f, 'BP_RedGrn', outDir);
    close(f);
end

