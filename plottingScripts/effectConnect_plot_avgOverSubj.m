%% script that loads data from effective connectivity (dDTF) from all subj
% warning: harcoded paths, ...
% (c) Jiri, Sep22

loadDataAgain = false;


% set subjects
params_default;
params.paradigm.usedParadigm = 'switchEI';
params.paradigm.specificType = 'switchEI';
params.storage.dir_shareData = 'G';
[params.storage.pathBeg, params.storage.subjList] = get_subjectList(params);
subjList = params.storage.subjList;
% subjList = { 
%     '20_PR3'; ...  
% };

if loadDataAgain
    % plot
    FB_subj = [];
    S_subj = [];
    for subj = 1:size(subjList,1)
        subjTag = subjList{subj,1};
    %     effectConnect_plot;
        [FB_y, S, freqAxis, timeAxis] = effectConnect_plot(subjTag);
    
        % cat FB
        if ~isempty(FB_y)
            FB_freq = [];
            for freq = 1:size(FB_y,1)
                FB_group = [];
                for g = 1:size(FB_y,2)
                    FB_group = cat(3, FB_group, FB_y{freq,g});  % 3D = time x clz x group
                end
                FB_freq = cat(4, FB_freq, FB_group);  % 4D = time x clz x group x freq
            end
            FB_subj = cat(5, FB_subj, FB_freq); % 5D = time x clz x group x freq x subj
        end
    
        % cat spectra S
        if ~isempty(S)
            S_group = [];
            for g = 1:size(S,1)
                S_clz = [];
                for clz = 1:size(S,2)
                    S_clz = cat(3, S_clz, S{g,clz});  % 3D = freq x time x clz
                end
                S_group = cat(4, S_group, S_clz);  % 4D = freq x time x clz x group
            end
            S_subj = cat(5, S_subj, S_group); % 5D = freq x time x clz x group x subj
        end
    end
    
    fileName = 'G:\dox\ms_switch_EI\data\v17_dtf_test\switchin_EI_IE_car\DTF_data.mat';
    save(fileName, 'S_subj', 'FB_subj');
else
    fileName = 'G:\dox\ms_switch_EI\data\v17_dtf_test\switchin_EI_IE_car\DTF_data.mat';
    assert(exist(fileName,'file') == 2);
    load(fileName, 'S_subj', 'FB_subj');
end

%% PLOT: SPECTRA from S_subj (avg over subjects)
f = fig_make;
nCols = size(S_subj,4);  % each group to each group
nRows = size(S_subj,3);  % classes below each other
nPlot = 1;
cLims = [-0.1, 0.1];
sel_ROIs = {'DMN -> DMN','DMN -> DAN','DAN -> DMN','DAN -> DAN'};
sel_Clz = {'E-I', 'I-E'};

% group-by-group connectivity -> C
for clz = 1:size(S_subj,3)
    for g = 1:size(S_subj,4)

        % data
        C = nanmean(S_subj(:,:,clz,g,:),5);     % 2D: f x t

        % plot
        subplot(nRows, nCols, nPlot);
        hold on;
        h = imagesc(timeAxis, freqAxis, C, cLims);
        colormap(gca,brewermap(256,'*RdBu'));   
        axis tight;
        title(['dDTF: ' sel_ROIs{g} ', clz = ' sel_Clz{clz}]);
        xlabel('time (s)');
        ylabel('freq (Hz)');
        xticks = [-2:2];
        for t = xticks
            plot([t t], ylim, '--k');   % time ticks
        end
        yticks = [10, 30, 50, 100, 120];
        for t = yticks
            plot(xlim, [t t], '--k');   % freq ticks
        end
        nPlot = nPlot+1;
    end
end

% save figure
figname = 'avgSubj_dDTF_chGroups_spectra';
outDir = 'G:\dox\ms_switch_EI\data\v17_dtf_test\switchin_EI_IE_car\coherency_channelGroups';                
fig_save(f, figname, outDir);
close(f); 

%% PLOT: FB from FB_subj (avg over subjects)
params_default;
params.triggering.freqBands = {...
    'delta',[0, 3];
    'theta',[4, 7];
    'alpha',[8, 12];
    'beta',[13, 30];
    'loGamma',[30, 45];
    'hiGamma',[55, 120];
    };
sel_ROIs = {'DMN -> DMN','DMN -> DAN','DAN -> DMN','DAN -> DAN'};
sel_Clz = {'E-I', 'I-E'};
clrs = {'b','r'};
f = fig_make;
nRows = size(FB_subj,4);    % freq. bands
nCols = size(FB_subj,3);                      % each group to each group
nPlot = 1;

% group-by-group connectivity -> C
for freq = 1:size(FB_subj,4)
    for g = 1:size(FB_subj,3)
        selFreq = params.triggering.freqBands{freq,2};

        % data
        y_avg = nanmean(FB_subj(:,:,g,freq,:),5);       % 2D: T x clz
        y_sem = sem(FB_subj(:,:,g,freq,:),5);           % 2D: T x clz (SEM over subj!)

        % plot
        subplot(nRows, nCols, nPlot);
        hold on;
        for clz = 1:size(FB_subj,2)
            h_ax = plotband(timeAxis, y_avg(:,clz), y_sem(:,clz), clrs{clz});
        end
        axis tight;
        title(['dDTF: ' sel_ROIs{g} ', freq = ' params.triggering.freqBands{freq,1}])
        xlabel('time (s)');
        ylabel('dDTF');
        xticks = [-2:2];
        for t = xticks
            plot([t t], ylim, '--k');   % time ticks
        end
        nPlot = nPlot+1;

    end
end

% save figure
figname = 'avgSubj_dDTF_chGroups_freqBands';
outDir = 'G:\dox\ms_switch_EI\data\v17_dtf_test\switchin_EI_IE_car\coherency_channelGroups'; 
fig_save(f, figname, outDir);
close(f);    
