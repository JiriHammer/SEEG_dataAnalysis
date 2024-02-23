function xyz = getSelChnlsMNI(params, dataTag)
% loads header H and finds the used channels and returns their MNI coors
% in case of bipolar reference, returns mean of the 2 bipolar channels MNI

% (c) Jiri, Feb17, Nov18

if nargin < 2
    dataTag = 'resp';
end

%% load header structure 'H'
load(params.storage.cacheFile, 'H');

%% load selCh_H
varName = ['selCh_H_' dataTag];
load(params.storage.cacheFile, varName);
selCh_H = eval(varName);

%% return channels info
xyz = nan(3,size(selCh_H,2));
for ch = 1:size(selCh_H,2)
    tmp = [];
    for n = 1:size(selCh_H,1)
        thisCh = selCh_H(n,ch);
        tmp = cat(2, tmp, [H.channels(thisCh).MNI_x; H.channels(thisCh).MNI_y; H.channels(thisCh).MNI_z]);
    end
    xyz(:,ch) = mean(tmp,2);
end
