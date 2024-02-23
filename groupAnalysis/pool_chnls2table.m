function pool_chnls2table(params, groupInfo)
% writes a list of iEEG into Excel table
% 1 row = 1 channel

% (c) Jiri, Oct23

%% check if output dir exists
if exist([params.storage.dir_results filesep params.storage.outName],'dir') ~= 7
    error(['output directory for results not found, dir = ' params.storage.dir_results filesep params.storage.outName]);
end
disp('Pooling analysis: creating channels table ...');

%% freq. bands
freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands

%% create table
TBL_all = [];
TBL_header = [];
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    TBL_subj = [];

    % cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    
    % load structure H & selected channels
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');   
    selCh_H = selCh_H_resp; 
    
    % >>> channel info <<<
    for ch = 1:size(selCh_H,2)
        thisCh = selCh_H(ch);
        m = 1;
        TBL_subj{ch,m} = subjTag;                             TBL_header{1,m} = 'subj';   m = m+1;    
        TBL_subj{ch,m} = H.channels(thisCh).name;             TBL_header{1,m} = 'chnl';   m = m+1; 
        TBL_subj{ch,m} = H.channels(thisCh).ass_marsLat_name(1);TBL_header{1,m} = 'lat.'; m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).neurologyLabel;   TBL_header{1,m} = 'nrlg';   m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).ass_yeo7_name;    TBL_header{1,m} = 'Yeo7';   m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).ass_yeo7_dist;    TBL_header{1,m} = 'd';      m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).ass_yeo17_name;   TBL_header{1,m} = 'Yeo17';  m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).ass_yeo17_dist;   TBL_header{1,m} = 'd';      m = m+1; 
        TBL_subj{ch,m} = H.channels(thisCh).ass_mars_name;    TBL_header{1,m} = 'Mars';   m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).ass_mars_dist;    TBL_header{1,m} = 'd';      m = m+1;      
        TBL_subj{ch,m} = H.channels(thisCh).MNI_x;            TBL_header{1,m} = 'x-MNI';  m = m+1;
        TBL_subj{ch,m} = H.channels(thisCh).MNI_y;            TBL_header{1,m} = 'y-MNI';  m = m+1; 
        TBL_subj{ch,m} = H.channels(thisCh).MNI_z;            TBL_header{1,m} = 'y-MNI';  m = m+1;
    end
        
    % significance of freq. band activations -> SGNF_fb
    SGNF_fb = zeros(size(selCh_H,2),size(freqBands,1));
    for freq = 1:size(freqBands,1)
        TBL_header{1,m} = ['S ' freqBands{freq,1}];  m = m+1;
    end
    for freq = 1:size(freqBands,1)
        
        % is freq. band activation (e.g. trialsData_beta) in cache file?
        varName = ['trialsData_' freqBands{freq,1}];  
        if ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists

            % load freq. band activations
            load(cacheFile, varName);
            trialsData = eval(varName);
            assert(size(trialsData.yVals,2) == size(selCh_H,2));
            
            for ch = 1:size(trialsData.yVals,2)
                SGNF_fb(ch,freq) = ch2roi_getSignificance(groupInfo, 'was_sgnf', trialsData.hVals(:,ch), trialsData.xVals);
            end
            
            clear trialsData;
            clear(varName);
        end
    end
    TBL_subj = cat(2, TBL_subj, num2cell(SGNF_fb));
    
    % cat to all subj
    TBL_all = cat(1, TBL_all, TBL_subj);
end

%% include TBL_header
TBL_all = cat(1, TBL_header, TBL_all);

%% write table: export to excell
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'chnlsTable'];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end 
fileName = ['chnlsTable.xlsx'];
xlswrite([outDir filesep fileName],TBL_all);
disp(['Table saved into: ' [outDir filesep fileName]]);


