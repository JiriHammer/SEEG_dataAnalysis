function poolSubjects_anatomyAreas(params, AA_type)
% sorts channels into defined anatomical areas, defined by 'AA_type'
% pools together subjects and frequency bands processed by a given job
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Jan18

disp('Pooling results over all subjects ...');

%% default anatomy area
if nargin < 2 || isempty(AA_type)
    AA_type = 'neurologists_anatomy';   
end

%% define anatomy areas
AA_list = anatomicalAreas_getList(params, AA_type);

%% plot frequency bands + ERP (over subjects and anatomy areas)
freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands
for freq = 1:size(freqBands,1)
    varName = ['trialsData_' freqBands{freq,1}]; 
    d = cell(size(AA_list,1),1);            % data
    sbj = cell(size(AA_list,1),1);          % subjects
    
    % go thru all subjects
    for subj = 1:size(params.storage.subjList,1)
        subjTag = params.storage.subjList{subj,1};

        % cache file
        cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
        assert(exist(cacheFile,'file') == 2);
        
        % load structure H & selected channels
        clear H selCh_H_resp;
        load(cacheFile, 'H', 'selCh_H_resp');
        selCh_H = selCh_H_resp;
        
        % try to load freq. band activation from cache file
        if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
            disp(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
        else
            clear(varName);
            load(cacheFile, varName);
            assert(exist(varName,'var') == 1);
            eval(['trialsData = ' varName ';']); 
            clear(varName);
            assert(size(trialsData.yVals,2) == size(selCh_H_resp,2));
            assert(size(trialsData.yVals,3) == size(trialsData.info.clzNames,1));
            
            % allocate data for all defined classes (some classes may be missing!) (e.g. in case subject did not do a freeflow paradigm)
            nClz_allDef = size(params.triggering.classes,1);
            y = nan(size(trialsData.yVals,1),size(trialsData.yVals,2),nClz_allDef);
            for clz = 1:size(trialsData.yVals,3)
                [tf, i_c] = ismember(trialsData.info.clzNames{clz},params.triggering.classes(:,1));
                assert(tf);
                y(:,:,i_c) = trialsData.yVals(:,:,clz);
            end

            % sort used channels to AA
            for ch = 1:size(selCh_H,2)
                
                % get anatomy area name of each channel
                thisAA = [];
                for n = 1:size(selCh_H,1)
                    thisCh = selCh_H(n,ch);                               % chnl index in H  
                    thisAA = cat(2, thisAA, anatomicalArea_getName(AA_type, AA_list, H.channels(thisCh)));
                end

                % was assigned?
                wasAssigned = false;
                for a = 1:size(thisAA,2)
                    [aa_found,i_aa] = ismember(lower(thisAA{a}), lower(AA_list(:,1)));
                    if aa_found
                        d{i_aa} = cat(2, d{i_aa}, y(:,ch,:));
                        sbj{i_aa} = cat(2, sbj{i_aa}, subj);
                        wasAssigned = true;
                    end
                end
                if ~wasAssigned
                    for n = 1:size(selCh_H,1)
                        thisCh = selCh_H(n,ch);                               % chnl index in H 
                        disp([subjTag ', ch = ' H.channels(thisCh).name '(' num2str(thisCh) ') - not assigned: neurologist = ' H.channels(thisCh).neurologyLabel ', cytoarch = ' H.channels(thisCh).ass_cytoarchMap ', toolbox = ' H.channels(thisCh).ass_brainAtlas]);
                    end
                end
            end      
        end
    end
    
    % rewrite trialsData structure with values in 'd'
    y_vals = [];
    y_errs = [];
    chNames = [];
    for aa = 1:size(d,1)
        if size(d{aa},2) >= 2    % at least 2 different channels contribute
            y_vals = cat(2, y_vals, nanmean(d{aa},2));
            y_errs = cat(2, y_errs, sem(d{aa},2));
            chNames{1,end+1} = {AA_list{aa}, ['nCh = ' num2str(size(d{aa},2)) ', nSubj = ' num2str(numel(unique(sbj{aa})))]};
        end
    end
    y_vals(isnan(y_vals)) = 0;  % eliminate NaNs for plotting (non-existent data will be zero!)
    y_errs(isnan(y_errs)) = 0;
    trialsData.yVals = y_vals;
    trialsData.yErrs = y_errs;
    trialsData.info.chNames = chNames;
    trialsData.info.outliers = zeros(size(chNames,2),1);
    assert(size(trialsData.info.colors,1) >= size(params.triggering.classes,1));
    trialsData.info.text = ['MEAN +/- SEM: ALL subjects, ' freqBands{freq,1} ', triggered: ' params.triggering.cutPoint ', classes: '];
    for clz = 1:size(params.triggering.classes,1)
        clr = trialsData.info.colors(clz,:);
        trialsData.info.text = [trialsData.info.text ...
            '\color[rgb]{' num2str(clr(1)) ' ' num2str(clr(2)) ' ' num2str(clr(3)) '}' params.triggering.classes{clz,1} ' '];
    end
    if isfield(trialsData, 'hVals'), trialsData = rmfield(trialsData, 'hVals'); trialsData = rmfield(trialsData, 'pVals'); end
    if isfield(trialsData, 'hVals_base'), trialsData = rmfield(trialsData, 'hVals_base'); trialsData = rmfield(trialsData, 'pVals_base'); end
    
    % plot
    trialsData.info.outDir = [params.storage.dir_results filesep params.storage.outName];
    trialsData.info.figName = [trialsData.info.figName '_ALL'];    
    plotTrials(params, trialsData);
    
end
        
%% plot spectra (over subjects and anatomy areas)
for clz = 1:size(params.triggering.classes,1)
    s = cell(size(AA_list,1),1);        % data
    clear spectralData;
    for subj = 1:size(params.storage.subjList,1)
        subjTag = params.storage.subjList{subj,1};

        % cache file
        cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
        assert(exist(cacheFile,'file') == 2);

        % load structure H & selected channels
        load(cacheFile, 'H', 'selCh_H_resp');
        selCh_H = selCh_H_resp;

        % load spectra from cache file
        varName = ['spectralData_clz' num2str(clz)]; 
        if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
            disp(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
        else        
            load(cacheFile, varName);
            assert(exist(varName,'var') == 1);
            eval(['spectralData = ' varName ';']); 
            clear(varName);
            assert(size(spectralData.cVals,3) == size(selCh_H,2));

            % sort used channels to AA
            for ch = 1:size(selCh_H_resp,2)
                % get anatomy area name of each channel
                thisAA = [];
                for n = 1:size(selCh_H,1)
                    thisCh = selCh_H(n,ch);                               % chnl index in H  
                    thisAA = cat(2, thisAA, anatomicalArea_getName(AA_type, AA_list, H.channels(thisCh)));
                end

                % was assigned?
                for a = 1:size(thisAA,2)
                    [aa_found,i_aa] = ismember(lower(thisAA{a}), lower(AA_list(:,1)));
                    if aa_found
                        s{i_aa} = cat(3, s{i_aa}, spectralData.cVals(:,:,ch));
                    end
                end
            end  
        end
    end
    
    % rewrite spectralData structure with values in 's'
    if exist('spectralData','var')
        c_vals = [];
        chNames = [];
        for aa = 1:size(s,1)
            if size(s{aa},3) >= 2    % at least 2 different channels contribute
                c_vals = cat(3, c_vals, mean(s{aa},3));
                chNames{1,end+1} = {AA_list{aa}, ['nCh = ' num2str(size(s{aa},3)) ', nSubj = ' num2str(numel(unique(sbj{aa})))]};
            end
        end
        spectralData.cVals = c_vals;
        spectralData.info.chNames = chNames;
        spectralData.info.text = ['ALL subjects, triggered: ' params.triggering.cutPoint ', class: ' params.triggering.classes{clz,1}];

        % plot
        spectralData.info.outDir = [params.storage.dir_results filesep params.storage.outName];
        spectralData.info.figName = ['spectra_AA_clz' num2str(clz)];    
        plotSpectra(params, spectralData);
    end

end



