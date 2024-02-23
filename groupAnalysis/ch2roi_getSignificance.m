function useCh = ch2roi_getSignificance(groupInfo, signif_criteria, ch_hVals, ch_xVals, method)
% determines if to use the channel based on 'signif_criteria'
% if the channel is significantly activated for longer than a threshold
% value, useCh = true (else false)

if nargin < 5
    method = 'contEpochs';
end

%% selected time of interest = groupInfo.time2plot
assert(size(ch_hVals,1) == size(ch_xVals,1));
if ~isfield(groupInfo, 'time2plot')
    groupInfo.time2plot = [ch_xVals(1), ch_xVals(end)];
end

i_t = closestval(ch_xVals,groupInfo.time2plot(1)):closestval(ch_xVals,groupInfo.time2plot(2));
t_step = mean(diff(ch_xVals));

%% init: default
useCh = 0;    % if to include the channel based on significance criteria
t_sig = 0;

%% time of significance
if strcmp(method, 'sampleNumber')
    t_sig = t_step * sum(ch_hVals(i_t) == 1,1);
end

%% max time duration of continuous significance "epoch" (epoch = part of data)
if strcmp(method, 'contEpochs')
    ch_hVals = ch_hVals(i_t);       % only selected samples
    ch_xVals = ch_xVals(i_t);
    i_sel = find(ch_hVals > 0);     % indices (samples) higher than threshold
    i_event = find(diff(cat(1, -666, i_sel, 666666666666)) > 1); % inds of continuous events
    t_sig = 0;
    for n = 2:length(i_event)
        % samples of the event (cont. data chunks > thr)
        i_cont = i_sel(i_event(n-1)):i_sel(i_event(n)-1);       % indices of one continuous "event"
        t_tmp = t_step * length(i_cont);                        % time of the significant cont. epoch
        if t_tmp > t_sig
            t_sig = t_tmp;
        end
    end
end

%% significance time threshold
if ~isfield(groupInfo, 'significanceTimeThr')
    groupInfo.significanceTimeThr = 0.25;
end
t_thr = groupInfo.significanceTimeThr;

%% is the channel significantly activated for longer than a threshold?
if strcmp(signif_criteria,'was_sgnf') && (t_sig >= t_thr)
    useCh = 1;    % -> 'was_sgnf';
elseif strcmp(signif_criteria,'not_sgnf') && (t_sig < t_thr)
    useCh = 1;    % -> 'not_sgnf'; 
elseif strcmp(signif_criteria,'any_sgnf')
    useCh = 1;    % -> 'all';
end


    
%% TBD, old implementation (now = Sep23)
% useCh = 0; 
% if strcmp(signif_criteria,'was_sgnf') && any(ch_hVals == 1)
%     useCh = 1;    % 'was_sgnf';
% elseif strcmp(signif_criteria,'not_sgnf') && all(ch_hVals ~= 1)
%     useCh = 1;    % 'not_sgnf'; 
% elseif strcmp(signif_criteria,'any_sgnf')
%     useCh = 1;    % 'all';
% end
% 
% % from ch2roi_load
%                 useCh_significance = 0;    % if to include the channel based on significance criteria
%                 if strcmp(list_anatomy_signif{1,3},'was_sgnf') && any(trialsData.hVals(:,ch) == 1)
%                     useCh_significance = 1;    % 'was_sgnf';
%                 elseif strcmp(list_anatomy_signif{1,3},'not_sgnf') && all(trialsData.hVals(:,ch) ~= 1)
%                     useCh_significance = 1;    % 'not_sgnf'; 
%                 elseif strcmp(list_anatomy_signif{1,3},'any_sgnf')
%                     useCh_significance = 1;    % 'all';
%                 end