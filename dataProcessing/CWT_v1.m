function [PWD, F, TW]=CWT_v1(signal, time, labels, setting)
% Function computes power of continuous wavelet transfor (P-CWT) in dB from trials data for
% each channel. Single precision used to save memory and computational time
% INPUTS:
% trials_data... EEG (time X channel X trials)
% trials_time... time vector in seconds (time X 1 X 1)
% setting... structure of individual settings
%   setting.freq_lim... frequency boundaries for CWT [start stop] ([0 fs/2] DEFAULT)
% OUTPUTS==================================================================
% Px... power scalegram (freq. X time X channel)
% freq_axis... frequency axis (freq. X 1)
% time_axis... time axis (1 X time)
% CWT_setting... structure of used settings

% setting -----------------------------------------------------------------
baseline_time=[-4.5 -3]; % reflects trials_time in seconds [start stop] 
fs=1/median(diff(time)); % sampling frequency

if setting.w_freg_lim(2)>fs/2
    flim=[setting.w_freg_lim(1) fs/2]; % frequency boundaries for CWT [start stop]
else
    flim=setting.w_freg_lim;
end
TW=time(:)';

disp(['CWT: fs=' num2str(fs) ' Hz; band=[' num2str(flim(1)) ':' num2str(flim(2)) '] Hz']);

fprintf(1,'CWT: ')
for ch=1:size(signal,2) % for each channel separately
    word=[num2str(100*ch/size(signal,2),'%.0f') '%%'];
    fprintf(1,word);
    for TR=1:size(signal,3)

        %  >>> CWT <<< in single precision to speed-up
        [wx,F]=cwt(single(signal(:,ch,TR)),'amor',fs,'FrequencyLimits',flim);
        
        % allocation of memory (after first computation)
        if ch==1 && TR==1 % declaration of "power-wavelet density (PWD)" matrix
            PWD=zeros(length(F),size(signal,2),size(wx,2),length(labels),'single'); % (F x CH x T x trials)
        end

        PWD(:,ch,:,TR) = abs(wx).^2; % Power of CWT = 4D: F x CH x T x trials
    end
    fprintf(1,repmat('\b',[1 length(word)-1]))
end
PWD=permute(PWD,[1 3 2 4]); % Power of CWT = 4D: F x T x CH x trials
disp('--- CWT done ---')

% dB ------------------------------------------
% PWD=10*log10(PWD);

% % average baseline ------------------------------------------
% bs_idx=time_axis<baseline_time(2) & time_axis>baseline_time(1);
% w2x_dB_baseline=mean(Px(:,bs_idx,:),2);
% 
% % whitening ------------------------------------------
% Px=Px-repmat(w2x_dB_baseline,[1 size(wx,2) 1 1]);

