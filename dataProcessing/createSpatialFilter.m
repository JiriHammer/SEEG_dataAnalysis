function filterMatrix = createSpatialFilter(params, N_inputCh, filterSettings)
% creates a spatial filter for (intracranial) EEG data
% input vars:
%   N_inputCh: number of input channels
%   selCh: vector of selected channel indices
%   filterSettings: struct, CAR, BIP, ...
% output var:
%   filterMatrix: N_inputCh x N_outputCh, where N_outputCh = numel(outputCh)
% example: 
%   S_car = createSpatialFilter(110, 65:78, 'car');
% usage:
%   X_car = X_ref * S_car;   where:
%       X_ref = [samples x N_inputCh] data matrix
%       X_car = [samples x N_outputCh] data matrix
%       S_car = [N_inputCh x N_outputCh] spatial filter (CAR) matrix

% (c) Jiri, Feb 2016


%  init
filterFound = false;
           
%% CAR: common average reference
if strfind(filterSettings.name, 'car')        % ~ common average re-reference (CAR)
    filterFound = true;
    
    % define channel groups
    chGroups = getChannelGroups(params, filterSettings.chGroups);
    
    % design filter
    filterMatrix = zeros(N_inputCh);                        % init
    for grp = 1:size(chGroups,2)
        selCh = chGroups{grp};  
        assert(~isempty(selCh));
        numCh = size(selCh,2);
        filterMatrix(selCh,selCh) = eye(numCh) - 1/numCh.*ones(numCh);      % set weights for CAR channels
    end  
end

%% BIP: bipolar reference
if strcmp(filterSettings.name, 'bip')
    filterFound = true;
    
    % define channel groups & re-define channel pointers to struct H
    chGroups = getChannelGroups(params, 'bip');
    
    % design filter
    filterMatrix = zeros(N_inputCh, size(chGroups,2));      % init
    for grp = 1:size(chGroups,2)
        selCh = chGroups{grp};
        assert(~isempty(selCh));
        filterMatrix(selCh(1),grp) = 1;                     % set weights for BIP channels
        filterMatrix(selCh(2),grp) = -1;                    % set weights for BIP channels
    end  
end

%% BIP EOG: (hardcoded for FR EEG experiments)
if strcmp(filterSettings.name, 'bip_eog')
    filterFound = true;
    
    % design filter
    filterMatrix = [1, -1, 0,  0;       % hardcoded: EOG_L - EOG_R
                    0,  0, 1, -1]';     % hardcoded: EOG_U - EOG_D 
                
    % load 'H' & selCh
    name_selCh = ['selCh_H_' params.thisData];
    load(params.storage.cacheFile, 'H', name_selCh);
    selCh_H = eval(name_selCh);    
    selCh_H = selCh_H([1,3]);
    
    % update channel selection (change of data size)
    if params.thisSess == size(params.storage.sessionCacheFiles,2)
        name_selCh = ['selCh_H_' params.thisData];
        eval([name_selCh '=selCh_H;' ]); 
        save(params.storage.cacheFile, name_selCh, '-append');
    end
end

%% NAN: no spatial filter
if strcmp(filterSettings.name, 'nan')
    filterFound = true;
    filterMatrix = eye(N_inputCh);
end

assert(filterFound);
