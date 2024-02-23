function [Px, F, T] = PSD_v1(signal, time, labels, setting)
% Function computes power spectral density (PSD) in dB from trials data for
% each channel.
% INPUTS===================================================================
% signal... EEG (time X channel X trials)
% time... time vector in seconds (time X 1 X 1)
% labels... (1 x trials) defines trials to avaraging
% setting... structure of individual settings


% OUTPUTS==================================================================
% Px_avr... average power spectrum in dB (freq. X time X channel X label_type)
% F... frequency axis (freq. X 1)
% T... time axis (1 X time)
% setting... structure of used settings
% setting.ww... segmentation window (sec.)
% setting.nn... negative overlap x100% 0.9~90%


% setting -----------------------------------------------------------------
fs=1/median(diff(time)); % sampling frequency from time vector
disp(['fs=' num2str(fs) 'Hz'])

ww=round(setting.ww*fs);  % window size (samples)
nn=round(setting.nn*ww); % noverlap (samples)
nfft=2^nextpow2(ww); % FFT Zero padding (n^2) 
win=hann(ww); % hanning window

index=1:ww-nn:(size(signal,1)-ww)+1; % index of first sample of segment window (T)
nseg=length(index);

% PSD
Px=zeros(nfft/2,nseg,size(signal,2),size(signal,3)); % freq. X time X channels X trials
fprintf(1,'Spectrogram: ')
for TR=1:size(signal,3) % for each trials TR
    % verbose --------
    word=[num2str(100*TR/size(signal,3),'%.0f') '%%'];
    fprintf(1,word);
    % ----------------
    
    % ======= PSD spectrogram ===========
    X=coherence_v2(signal(:,:,TR),win,nn,nfft,fs); % fast FFT implementation of PSD (F x CH x T)
    % =======================
    
    X=2*X(1:end/2,:,:); % one-side spectrum, second half is added to first (2x)
    X=permute(X,[1 3 2]); % (F x T x CH)
    
    Px(:,:,:,TR)=(1/nfft)*abs(X.^2); % PSD (F x T x CH x TR)
    
    fprintf(1,repmat('\b',[1 length(word)-1]))
end
disp('--- PSD done ---')


T=time(index); % time of segments (start of window)
T=T(:);
F=linspace(0,fs-fs/nfft,nfft)'; % frequency axis
F=F(1:end/2); % one-side spectrum

% (uncommeted) average trials PSD
% Px_var=zeros(size(Px,1),size(Px,2),size(Px,1),length(unique(labels)));
% for L=unique(labels)
%     Px_avr(:,:,:,L)=mean(Px(:,:,:,labels==L),4);
% end

% (uncommeted) dB ------------------------------------------
% Px=10*log10(Px); % power in dB

% (uncommeted) average baseline ------------------------------------------
% bs_idx=T<baseline_time(2) & T>baseline_time(1);
% pxdb_baseline=mean(Px(:,bs_idx,:,:),2);
% 
% % whitening ------------------------------------------
% Px=Px-repmat(pxdb_baseline,[1 size(Px,2) 1 1]);





