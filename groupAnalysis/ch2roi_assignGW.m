function useCh_GW = ch2roi_assignGW(gm_wm_selection, H_channel, gm_thr)
% was channel assigned to: GREY matter or WHITE matter or CSF or BONE ?
% uses channel structure H_channel to find if 
%   (1) it matches the gm_wm_selection
%   (2) the probablity is below threshold gm_thr
% input:
%   H_channel = channel structure from H.channels
%        - required field: H_channel.ass_gmwm_name
%        - required field: H_channel.ass_gmwm_dist
%   gm_wm_selection = string: 'gm','wm','csf','bone', 'any'
%   gm_thr = probablity threshold in [%]
% returns:
% useCh_GW = true;  if channel should be used (grouped)
% useCh_GW = false; if channel should NOT be used (grouped)

% (c) Jiri, Nov22

if ~isfield(H_channel, 'ass_gmwm_name')
    disp('WARNING: grey / white matter not assigned.');
    useCh_GW = true;
    return;
end

%% gm_wm_selection = any (no matter if GM, WM, ...)
if strcmp(gm_wm_selection, 'any')
    useCh_GW = true;
    return;
end

%% GM / WM / CSF / BONE ?
[tf, i_col] = ismember(gm_wm_selection,H_channel.ass_gmwm_name); 
assert(tf);
assert(strcmp(H_channel.ass_gmwm_name{i_col}, gm_wm_selection));
if H_channel.ass_gmwm_dist(i_col) >= gm_thr/100
    useCh_GW = true;
else
    useCh_GW = false;
end

