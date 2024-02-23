% function [ROL]= ROLbootstrap_NC(data,ROL_params)
% INPUTS:
% data: epoched data
%       .wave: actual data
%       .time: timepoint corresponding to each data pt
%       .fsample: sampling rate
% params    (see genROLparams_chao.m)

% Chao Zhang - LBCN: https://github.com/LBCN-Stanford/lbcn_preproc/blob/master/ROL/ROLbootstrap_NC.m
% Updates for ver3, BLF - LBCN 2014.
% Cleaned up code, added comments
% Modified moving window, larger with more overlap
% modified by Amy Daitch- 2016

% Questions Jiri
% - ROL on single trials???
% - when does the bootstrap kick in ?

%% ROL params
% see: https://github.com/LBCN-Stanford/lbcn_preproc/blob/master/ROL/genROLParams_NC.m
ROL_params = struct;
ROL_params.thr_value = 1; % How many SD over baseline period
ROL_params.thr_value_counter = 25; % How many bins must surpass threshold, TO DO: in [s]
% ROL_params.bin_timewin = 0.002; % How much non-overlap between bins, 2 [ms] ?
ROL_params.bin_timewin = 0.078125;   % Jiri
ROL_params.times_nonover = 15; % 30ms bins with 28ms overlap, ???

ROL_params.blc = 1; % chao
ROL_params.power = false;
ROL_params.smooth = false;
ROL_params.deactive = false;    % ???

% Set range parameters
% ROL_params.end_value = .5;    % 500 ms
ROL_params.end_value = 1;       % in [s] ???
ROL_params.start_value = 0;     % stim_onset, in [s] ???
ROL_params.bas_value = 0;       % end of baseline, in [s] ?
ROL_params.baseline = [-0.2 0]; % time window of baseline

ROL_params.pre_event = -0.2;    % baseline period used for baseline correction
ROL_params.dur = 1;             % 1000ms, ?
        
% --- Jiri ---
t_base = [-0.9, -0.1];
i_base = closestval(trials.time, t_base(1)):closestval(trials.time, t_base(2));
t_stim = [0.0, 1.0];
i_stim = closestval(trials.time, t_stim(1)):closestval(trials.time, t_stim(2));

windSize = 0.078125;    % in [s], = 10/128 (10 samples @fs = 128 Hz)
timeStep = 0.015625;    % in [s], =  2/128 (2  samples @fs = 128 Hz)
fs = 1/mean(diff(trials.time)); % in [Hz]

%% ROL algorithm
% alldata = data.wave';   % 2D: time x trials

% --- Jiri ---
ch = 1;         % TO DO
alldata = squeeze(trials.data(:,ch,:));   % 2D: time x trials

% Normalize trial
m = max(alldata(i_stim,:),[],1);
alldata = alldata./repmat(m,size(alldata,1),1);

% Calculate baseline mean and standard deviation
baseline_signal = mean(alldata(i_base,:),2);    % mean over trials
mean_bl = mean(baseline_signal);        % mean over baseline time
std_bl = std(baseline_signal);          % std  over baseline time

thr_pos = mean_bl + ROL_params.thr_value*std_bl;
thr_neg = mean_bl - ROL_params.thr_value*std_bl;

%% number of bins
nTS = 0;        % number of time steps
t_sel = [t_stim(1),t_stim(1)+windSize];  % 1st time window
tStep_samples = timeStep*fs;
i_t = closestval(trials.time, t_sel(1)):closestval(trials.time, t_sel(2));
while i_t(end) <= i_stim(end)
    nTS = nTS+1;
    i_t = i_t + tStep_samples;
end

%% Loop through trials
peaks_pos = nan(1,size(alldata,2));     % trial vector for positive activations
onsets_pos= nan(1,size(alldata,2));     % trial vector
peaks_neg = nan(1,size(alldata,2));     % trial vector for negative deactivations
onsets_neg= nan(1,size(alldata,2));     % trial vector
for tr = 1:size(alldata,2)
    
    % data for specific trial
    mdata = alldata(:,tr)';     % 1D: 1 x time
    
    % Get number of bins
%     numbins = (((ROL_params.end_value-ROL_params.start_value)/ROL_params.bin_timewin)-ROL_params.times_nonover)+1;
    numbins = nTS;
    
    % Create cell structure to hold bin data
    bin_data = cell(floor(numbins),6); % bin_data will contain (data; averages; slopes; whether average surpasses threshold; index of first bin value; index of last bin value)
    
    %% Loop through bins and create data structure
    t_window = [t_stim(1), t_stim(1)+windSize];  % set 1st time window
    for b = 1:numbins
        i_samples = [closestval(trials.time,t_window(1)):closestval(trials.time,t_window(2))];
        bin_data(b,1) = {mdata(i_samples)};     % all samples
        bin_data(b,2) = {mean(bin_data{b,1})};  % mean over samples
        % Get slope
        coefficients = polyfit(1:length(bin_data{b,1}), bin_data{b,1}, 1);
        slope = coefficients(1);
        bin_data(b,3) = {slope};        
        
        % Check if bin mean exceeds thershold (1 = yes, 0 = no)
        if bin_data{b,2} > thr_pos
            bin_data(b,4) = {1};  
        elseif bin_data{b,2} < thr_neg
            bin_data(b,4) = {-1};  
        else
            bin_data(b,4) = {0};  
        end  
        bin_data(b,5) = {i_samples(1)};         % Get index of first value
        bin_data(b,6) = {i_samples(end)};       % Get index of last value
            
        t_window = t_window+timeStep;   % shift the time window
    end     % of bins
    
    %% ROL for POS: loop through bin data to find N consequtive bins meeting threshold criterion
    counter = 0;    % counter to track slope and threshold for activated trials
    peak_bin_index = 0; 
    for b = 1:length(bin_data) % would be better with 'while'                  
        % Check if bin is above threshold, with a positive slope
        if bin_data{b,4} == 1                      
            counter = counter + 1;                     
            % Check if counter exceeds window threshold
            if counter ==  ROL_params.thr_value_counter
                % save the value of tenth bin
                peak_bin_index = b;
                break
            end
        else
            counter = 0;
        end
    end
    
    % Initialize ROL for POS
    ROL_Bin = (peak_bin_index-ROL_params.thr_value_counter)+1;
    if (peak_bin_index-ROL_params.thr_value_counter) >= 0

        onset_index = bin_data{ROL_Bin,5};
        onset = trials.time(onset_index);

        % POS: Estimate signal's peak from bins
        [~, Peak_Bin]=max([bin_data{:,2}]);
        peak_index = bin_data{Peak_Bin,5};
        peak = trials.time(peak_index);

        peaks_pos(tr) = peak;                                   
        onsets_pos(tr) = onset;
    end

    %% ROL for NEG: loop through bin data to find N consequtive bins meeting threshold criterion
    counter = 0;    % counter to track slope and threshold for activated trials
    peak_bin_index_neg = 0; 
    for b = 1:length(bin_data) % would be better with 'while'                  
        % Check if bin is above threshold, with a positive slope
        if bin_data{b,4} == -1                      
            counter = counter + 1;                     
            % Check if counter exceeds window threshold
            if counter ==  ROL_params.thr_value_counter
                % save the value of tenth bin
                peak_bin_index_neg = b;
                break
            end
        else
            counter = 0;
        end
    end
    
    % Initialize ROL for NEG
    ROL_Bin = (peak_bin_index-ROL_params.thr_value_counter)+1;
    if (peak_bin_index-ROL_params.thr_value_counter) >= 0

        onset_index = bin_data{ROL_Bin,5};
        onset = trials.time(onset_index);

        % NEG: Estimate signal's peak from bins
        [~, Peak_Bin]=min([bin_data{:,2}]);
        peak_index = bin_data{Peak_Bin,5};
        peak = trials.time(peak_index);

        peaks_neg(tr) = peak;                                   
        onsets_neg(tr) = onset;
    end   
    
    %% degug figure
    if ~isnan(peaks_pos(tr)) || ~isnan(peaks_neg(tr))
        f = fig_make;
        plot(trials.time, mdata);
        hold on;
        axis tight;
        box on; grid on;
        plot(xlim, [thr_pos, thr_pos], '--k');
        plot(xlim, [thr_neg, thr_neg], '--k');
        plot([peaks_pos(tr),peaks_pos(tr)],ylim, 'r');
        plot([peaks_neg(tr),peaks_neg(tr)],ylim, 'm');
        why;
    end
end
    
ROL.peaks_pos= peaks_pos;
ROL.onsets_pos = onsets_pos;
ROL.peaks_neg= peaks_neg;
ROL.onsets_neg = onsets_neg;
