function [Cxy, F, T, ICxy] = COH_v2(signal, time, labels, setting)
% Function computes coherence and im. coherence.
% INPUTS===================================================================
% signal... EEG (time X channel X trials)
% time... time vector in seconds (time X 1 X 1)
% labels... (1 x trials) defines trials to avaraging
% setting... structure of individual settings
% setting.ww... segmentation window (sec.)
% setting.nn... negative overlap x100% 0.9~90%
% setting.int_ww... internal coh. segmentation window (sec.)
% setting.int_nn... internal coh. negative overlap x100% 0.9~90%

% OUTPUTS==================================================================
% Cxy... average coherence (freq. X time X channel X label_type)
% ICxy... average imaginary coherence (freq. X time X channel X label_type)
% F... frequency axis (freq. X 1)
% T... time axis (1 X time)


% setting -----------------------------------------------------------------
fs=1/median(diff(time)); % sampling frequency from time vector
disp(['fs=' num2str(fs) 'Hz'])

% major segmentation
ww=round(setting.ww*fs);  % window size (samples)
nn=round(setting.nn*ww); % noverlap (samples)
nfft=2^nextpow2(ww); % FFT Zero padding (n^2) 
win=hann(ww); % hanning window

idx=1:ww-nn:size(signal,1)-ww+1;

Cxy=zeros(size(signal,2),size(signal,2),nfft/2,length(idx),length(unique(labels)));
ICxy=zeros(size(signal,2),size(signal,2),nfft/2,length(idx),length(unique(labels)));
for T=1:length(idx)
    X=signal(idx(T):idx(T)+ww-1,:,:);
    X=X.*repmat(win,[1 size(X,2) size(X,3)]); % hanning window
    
    X=fft(X,nfft,1);
    
    U=(win(:)'*win(:))*(length(win)/2);  % normalization in accordance with cpsd() results
    
    Pxy=zeros(size(signal,2),size(signal,2),size(X,1));

    for L=1:length(unique(labels)) % for each label separately
        % averaging throught trials is used as variance for COH computation
        for ch=1:size(signal,2)
            % Cross-PSD (ch X ch X F) 
            % mean trought trials
            Pxy (ch,:,:)=permute(mean( (repmat(X(:,ch,labels==L),[1 size(signal,2) 1]).*conj(X(:,:,L==labels)))/U,3),[3 2 1]);
        end
        
        Pxx=repmat( reshape(Pxy(logical(repmat(eye(size(Pxy,1)),[1 1 size(Pxy,3)]))), [1 size(Pxy,2) size(Pxy,3)] ), [size(Pxy,1) 1 1]);
        Pyy=permute(Pxx,[2 1 3]);
        
        Cxy(:,:,:,T,L)=(Pxy(:,:,1:end/2).*conj(Pxy(:,:,1:end/2)))./(Pxx(:,:,1:end/2).*Pyy(:,:,1:end/2));
        ICxy(:,:,:,T,L)=Pxy(:,:,1:end/2)./sqrt(Pxx(:,:,1:end/2).*Pyy(:,:,1:end/2));
        ICxy(:,:,:,T,L)=(imag(ICxy(:,:,:,T,L)).^2)./(ones(size(ICxy(:,:,:,T,L)))-real(ICxy(:,:,:,T,L)).^2);
    end
end

T=time(idx); % time of segments (start of window)
T=T(:);
F=linspace(0,fs-fs/nfft,nfft)'; % frequency axis
F=F(1:end/2); % one-side spectrum




