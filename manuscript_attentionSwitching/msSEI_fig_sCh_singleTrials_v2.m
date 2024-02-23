%% single trials: raw traces
% - load trials
% - select channels
% - plot raw traces

% (c) Jiri, Feb24

%% settings
subjTag = '20_PR11';
selected_channels = {'Rd1-Rd2', 'Ld6-Ld7'};
selected_networks = {'DMN', 'DAN'};
pathBeg = 'F:\dox\ms_switch_EI\data\v31_stft_session_bip\switchin_EI_IE_bip';

%% load trials
params.storage.subjTag = subjTag;

% --- cache file
cacheFile = [pathBeg filesep subjTag filesep 'cacheFile.mat'];
assert(exist(cacheFile,'file') == 2);
% params.storage.cacheFile = cacheFile;               % used in connectivityAnalysis.m

% load raw trials
clear trials
load(cacheFile, 'trials');  

% load raw trials
clear H selCh_H_resp
load(cacheFile, 'H', 'selCh_H_resp');  
selCh_H = selCh_H_resp;
assert(size(selCh_H,2) == size(trials.data,2));

%% select channels
i_ch = [];
for c = 1:size(selected_channels,2)
    for ch = 1:size(selCh_H,2)
        if strcmp(H.channels(selCh_H(ch)).name, selected_channels{c})
            i_ch = cat(2, i_ch, ch);
        end
    end
end

%% figures
time2plot = [-2.3, 2.3];
trials2plot = [50, 94];
clz2plot = [2, 1];  % 2 = I-E, 1 = E-I
trialType = {'I-E','E-I'};

marg_h = [0.35 0.05];
marg_w = [0.15 0.02];
gap = [0.01, 0.01]; 

for tr = 1:size(trials2plot,2)
    for ch = 1:size(selected_channels,2)
        % figure
        f = fig_make('fig_W',10, 'fig_H',6);
%         hold on;

        % data
        i_t = closestval(trials.time,time2plot(1)):closestval(trials.time,time2plot(2));
        y = trials.data(i_t,i_ch(ch),trials2plot(tr));              % add trial below spectrum
        t = trials.time(i_t);

        % axes
        ax = subtightplot(1, 1, 1, gap, marg_h, marg_w);
        hold on;
        
        % plot
%         plot(t, y);
        msSEI_plotband(t, y, zeros(size(y)), clz2plot(tr));

        % labels        
        axis tight;
        ylim([-85, 85]);
        xlabel('Time (s)');
        ylabel('iEEG (\muV)');
        grid on;
        box on;
        
        % save fig
        figName = ['v4_tr' num2str(trials2plot(tr)) '_' trialType{tr} '_ch' num2str(i_ch(ch)) '_' num2str(selected_networks{ch})];
        figDir = [pathBeg filesep 'msSEI_fig_sCh_singleTrials' filesep 'rawTrials'];
        fig_save(f, figName, figDir, 'res',600);
        % close(f);  
        
    end
end

    

    