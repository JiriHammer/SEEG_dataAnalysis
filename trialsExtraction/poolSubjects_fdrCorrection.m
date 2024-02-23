function poolSubjects_fdrCorrection(params, pLevel, fields2eval, freqBands)
% performs FDR (false discovery rate) correction for multiple testing
% pools together subject processed by a given job
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Feb19

%% defaults
if nargin < 2
    pLevel = 0.05;
end

if nargin < 3
    fields2eval = {'pVals','hVals'};
end

% activations (ERP + frequency bands)
if nargin < 4
    freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands
end

disp(['FDR correction: ' fields2eval{1} ', significance level = ' (num2str(pLevel))]);

%% pool p-values: load activations (over all subjects)
P_all = [];
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    disp([' - loading data, subj: ' subjTag]);
    
    % cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    
    % load freq. band activations
    for freq = 1:size(freqBands,1)
        
        % load freq. band activation from cache file
        varName = ['trialsData_' freqBands{freq,1}];  
        load(cacheFile, varName);
        if exist(varName,'var') == 1
            eval(['trialsData = ' varName ';']); 
            clear(varName);        
        else
            disp(['WARNING: subj = ' subjTag ' - trialsData_' freqBands{freq,1} ' not found']);
            break;
        end
        
        % only two classes selected? (only then are p-values computed)
        if size(trialsData.yVals,3) ~= 2 && strcmp(fields2eval{1},'pVals')
            disp('WARNING: p-values are computed only for 2 classes');
            disp('WARNING: No FDR correction computed.');
            return;
        else
            assert(isfield(trialsData, fields2eval{1}));
            tmp = trialsData.(fields2eval{1});
            assert(~any(isnan(tmp(:))));       % trialsData.pVals = 2D: time x chnls
            P_all = cat(1, P_all, tmp(:));     % cats p-vals columns-wise (for chnls along time dim)
        end
    end       
end

%% FDR correction
disp(' - computing FDR correction ...');
[n_signif,i_signifTest] = fdr(P_all, pLevel, 'original');
H_all = zeros(size(P_all));         % default: 0 -> not significant
if n_signif > 0
    H_all(i_signifTest) = 1;        % set those who passed the FDR correction to 1 -> was significant
end

% H_all2 = fdr_bh(P_all);  % different implementation of FDR (gave same results as above)

%% update h-values (if they survived the FDR correction)
n = 1;
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    disp([' - updating significances, subj: ' subjTag]);
    
    % cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    
    % load freq. band activations
    for freq = 1:size(freqBands,1)
        
        % load freq. band activation from cache file
        varName = ['trialsData_' freqBands{freq,1}];  
        load(cacheFile, varName);
        if exist(varName,'var') == 1
            eval(['trialsData = ' varName ';']); 
            clear(varName);        
        else
            disp(['WARNING: subj = ' subjTag ' - trialsData_' freqBands{freq,1} ' not found']);
            break;
        end        
        
        % update hVals
        %tmp = trialsData.pVals(:);
        tmp = trialsData.(fields2eval{1})(:);
        inds = n:n+size(tmp,1)-1;
        H_subj_vector = H_all(inds);
        H_subj_matrix = reshape(H_subj_vector, size(trialsData.(fields2eval{1})));
        trialsData.(fields2eval{2}) = H_subj_matrix;       % trialsData.hVals = 2D: time x chnls (x clz)
        
        % save trialsData
        outputVarName = ['trialsData_' freqBands{freq,1}];
        eval([outputVarName '=trialsData;' ]); 
        disp(['  -- saving trials data: ' freqBands{freq,1} ' ...']);
        save(cacheFile, outputVarName, '-append');  
        clear(outputVarName); 
        
        % update counter
        n = n + size(tmp,1);
    end       
end

