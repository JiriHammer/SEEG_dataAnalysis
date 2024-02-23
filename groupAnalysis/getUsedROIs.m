function [N, S] = getUsedROIs(GRUP, minSubj)
% returns N = vector(freq), number of used ROIs in trialsData of each frequency band (FB)
% returns S = number of used ROIs in spectra (significance taken from the last FB)
% used = more than minSubj

% (c) Jiri, Jul22

%% number of used ROIs in trialsData of each frequency band -> N
N = zeros(size(GRUP.trials,2),1);
for freq = 1:size(GRUP.trials,2)
    for aa = 1:size(GRUP.trials,1)
        if ~isempty(GRUP.trials{aa,freq})   % must be not empty
            if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= minSubj  % must contain at least minSubj
                N(freq,1) = N(freq,1) + 1;
            end
        end
    end
end

%% number of used ROIs spectra -> S
S = 0;
for aa = 1:size(GRUP.spectra_data,1)
    if ~isempty(GRUP.spectra_data{aa})   % must be not empty
        if size(unique(GRUP.nSubj_roi{aa,end,1}),1) > minSubj  % must contain at least minSubj, taken from the last FB !!!
            S = S + 1;
        end
    end
end
