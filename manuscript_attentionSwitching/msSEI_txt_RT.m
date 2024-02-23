%% ms SEI text: reaction times (RT_all) & right-wrong answers (RWM)
% (c) Jiri, Nov23

%% settings (v30)
% subjList = msSEI_getSubjList('SEI_withRS');
% pathBeg = 'F:\dox\ms_switch_EI\data\v30_stft_baseRS_bip\switchin_EI_IE_bip';

%% settings (v31)
subjList = msSEI_getSubjList('SEI_all');
pathBeg = 'F:\dox\ms_switch_EI\data\v31_stft_session_bip\switchin_EI_IE_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)

%% get RT & RWM
RT_avg = [];
RWM_avg = [];
nTr_avg = [];
for subj = 1:size(subjList,1)
    subjTag = subjList{subj,1};
    
    % cache file
    cacheFile = [pathBeg filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);    
    
    % load
    clear RT_all RWM params
    load(cacheFile, 'RT_all','RWM','params');
    
    % class labels
    clzLabels = unique(RT_all(2,:));
    
    % mean RT
    RT_subj = nan(1,length(clzLabels));
    for clz = 1:length(clzLabels)
        i_clz = RT_all(2,:) == clzLabels(clz);
        RT_subj(clz) = mean(RT_all(1,i_clz));
    end
    RT_avg = cat(1, RT_avg, RT_subj);
    
    % mean RWM
    for clz = 1:length(clzLabels)
        i_clz = RT_all(2,:) == clzLabels(clz);  
        RWM(clz,:) = RWM(clz,:)./sum(i_clz,2);
    end
    RWM_avg = cat(3, RWM_avg, RWM);
    
    % number of trials
    clear trialsData_alpha;
    load(cacheFile, 'trialsData_alpha');
    nTr_avg = cat(1, nTr_avg, trialsData_alpha.info.nClz);
end

%% display result RT
for clz = 1:size(RT_avg,2)
    disp(['clz = ' num2str(clz) ' = ' params.triggering.classes{clz,1} ', RT: mean +/- SEM = ' num2str(mean(RT_avg(:,clz),1)) ' +/- ' num2str(sem(RT_avg(:,clz),1)) ' s']);
    disp(['clz = ' num2str(clz) ' = ' params.triggering.classes{clz,1} ', RT: range = ' num2str(min(RT_avg(:,clz))) ' - ' num2str(max(RT_avg(:,clz))) ' s']);
end
P_RT = ranksum(RT_avg(:,1),RT_avg(:,2), 'alpha',0.05);
disp(['RT: P = ' num2str(P_RT)]);
disp('-------------------------')
disp('');

%% display result RWM
for clz = 1:size(RWM_avg,1)
    for rw = 1:2
        if rw == 1
            rw_str = 'right';
        else
            rw_str = 'wrong';
        end
        disp(['clz = ' num2str(clz) ' = ' params.triggering.classes{clz,1} ', answer = ' rw_str ': mean +/- SEM = ' num2str(mean(RWM_avg(clz,rw,:),3)) ' +/- ' num2str(sem(RWM_avg(clz,rw,:),3)) ' %']);
        disp(['clz = ' num2str(clz) ' = ' params.triggering.classes{clz,1} ', answer = ' rw_str ': range = ' num2str(min(RWM_avg(clz,rw,:))) ' - ' num2str(max(RWM_avg(clz,rw,:))) ' %']);
    end
end
disp('-------------------------')
disp('');

%% number of trials
for clz = 1:size(RT_avg,2)
    disp(['clz = ' num2str(clz) ' = ' params.triggering.classes{clz,1} ', nTrials: mean +/- SEM = ' num2str(mean(nTr_avg(:,clz),1)) ' +/- ' num2str(sem(nTr_avg(:,clz),1))]);
end
disp('-------------------------')
disp('');

%% age
disp('-------------------------')
age = [
31
26
38
41
15
16
44
48
12
15
43
46
29
35
55
46
37
44
36
38
32
28
24
49
30
];
disp(['age: mean +/- SEM = ' num2str(mean(age,1)) ' +/- ' num2str(sem(age,1)) ' years']);
disp(['age: range = ' num2str(min(age)) ' - ' num2str(max(age)) ' years']);
