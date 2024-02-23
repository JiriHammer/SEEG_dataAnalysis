function msSEI_FC_PLV_impl(params, trials, fileName_ec)
% computes PLV 

% (c) Jiri, Dec23

nCh = size(trials.data,2);
nClz = length(unique(trials.labels));
FC_freqBands = params.connectivity.freqBands;
M = nan(nCh,nCh,size(FC_freqBands,1),size(trials.data,1),nClz);     % 5D: [ch,ch,f,t,clz]

%% -----------PLV: for FB -> M matrix (saved)--------
for freq = 1:size(FC_freqBands,1)
    disp([' - computing DTF: freq band = ' num2str(freq) '/' num2str(size(FC_freqBands,1))]);
 
    % filter data: bp, HT, ampEnv -> trials_filt
    params.triggering.trialsFiltering = {'bandPass', 'hilbPhaseAngle', 'singlePrecision'};
    params.bp_freq.freqBand = FC_freqBands(freq,:);
    trials_plv = trialsFiltering(params, trials);
    C = trials_plv.data;
    
    % compute PLV as a sum of hilbPhases
    selLabels = unique(trials.labels);
    for clz = 1:nClz
        i_tr = find(trials.labels == selLabels(clz));
        for ch_i = 1:nCh
            for ch_j = ch_i:nCh
                phi_diff = C(:,ch_i,i_tr) - C(:,ch_j,i_tr);     % angles [-pi,pi], phase difference in each trial
                c_phi = exp(1i.*phi_diff);                      % complex number on a unit circle of the phase difference in each trial
                v_plv = abs(sum(c_phi,3));                      % abs of sum over trials
                M(ch_i,ch_j,freq,:,clz) = v_plv./size(i_tr,2);  % mean over trials
                M(ch_j,ch_i,freq,:,clz) = M(ch_i,ch_j,freq,:,clz);  % symmetric FC
            end
        end
    end
end % of freq

M = single(M); % save memory

%% save PLV
disp(['saving results from: ' params.connectivity.connectivityMeasure{1} ' ...']);
EC_info = struct;
EC_info.time = trials.time;
EC_info.freq = mean(FC_freqBands,2);
EC_info.clzLabels = unique(trials.labels);
% EC_info.selCh = selCh_all;
% EC_info.selCh_groups = selCh_groups;
EC_info.methodName = params.connectivity.connectivityMeasure{1};    % !!! hardcoded for 1st only !!!
EC_info.info = trials.info;
if exist(fileName_ec,'file') == 2
    save(fileName_ec, 'M', 'EC_info', '-append');
else
    save(fileName_ec, 'M', 'EC_info', '-v7.3');
end
disp(['saved to: ' fileName_ec]);
    


% save to PLV to filename
