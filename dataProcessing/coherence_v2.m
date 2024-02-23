function [X,Pxy,Cxy,F,TFxy]=coherence_v2(x,winsize,noverlap,nfft,fs)
% The function estimates two side Cross Power Spectral Density (CPSD) and Magnitude
% Squared Coherence (MSCOHERE) of multi-channel signal x each with each
% other. Estimation is performed using averaging of segments. The algorithm
% is optimized by matrix operation to fastes run.
% WARNING: large consumption memory! 
% VERSION: 2.00
% PROGRAMED BY: Radek Janca, CTU in Prague, Czech Republic
%
% ------------------------------------------------------------------------
% INPUTS:
% x ... data matrix, channels of signal is in columns [time x channel]
% winsize ... length of segmentation window (default hamming(winsize))
%         ... or window mask, example: hann(256), blackman(512) etc.
%         ... [] ~ round(size(x,1)/8)
% noverlap ...  the number of overlapping segments
%          ... [] ~ round(0.5*winsize)
% nfft ... length of each segment FFT: nfft=winsize+ added_nulls
%      ... 2^n is recommended
%      ... [] ~ without nfft
% fs ... sampling frequency or omega
%    ... [] ~ normalized frequency
% ------------------------------------------------------------------------
% OUTPUTS:
% X ... Complex spectrogram
% Pxy ... Cross Power Spectral Density like pwelch(CPSD): 3D matrix [X x Y x F]
%     ... Pxx is in P(X,X,:), X is channel index
% Cxy ... Magnitude Squared Coherence (MSCOHERE): 3D matrix [X x Y x F]
% F ... vector contains corresponding values of frequency
% TFxy ... Transfer Function: 3D matrix [X x Y x F]
%
% ------------------------------------------------------------------------
% Examples:
% [X,Pxy,Cxy,F,TFxy]=coherence(x); (default: 8 segmnets, 50% overlap, nearest 2^n nfft, normalized frequency)
% [X,Pxy,Cxy,F,TFxy]=coherence(x,winsize); (default: 50% overlap, nearest 2^n nfft, normalized frequency)
% [X,Pxy,Cxy,F,TFxy]=coherence(x,winsize,noverlap); (default: nearest 2^n nfft, normalized frequency)
% [X,Pxy,Cxy,F,TFxy]=coherence(x,winsize,noverlap,nfft); (default: normalized frequency)
% [X,Pxy,Cxy,F,TFxy]=coherence(x,winsize,noverlap,nfft,fs);


%%
if nargin==1
   winsize=round(size(x,1)/8);
   noverlap=round(winsize/2);
   nfft=2^nextpow2(winsize);
   fs=2;
end

if isempty(winsize); winsize=round(size(x,1)/8); end

if max(size(winsize))==1
    win=hamming(round(winsize));
else
    win=winsize(:);
    winsize=length(win);
end

if nargin==2
   noverlap=round(winsize/2);
   nfft=2^nextpow2(winsize);
   fs=2;
end

if nargin==3
    nfft=2^nextpow2(winsize);
    fs=2;
end

if isempty(noverlap); noverlap=round(winsize/2); end

if nargin==4
    fs=2;
end

if isempty(nfft); nfft=winsize; end
if isempty(fs); fs=2; end


%%
% mask of index
l=1:1:size(x,1);

L=buffer(l,winsize,noverlap,'nodelay');

if mod(size(x,1)-winsize,(winsize-noverlap))~=0
   L=L(:,1:end-1); 
end

L=repmat(L,[1 1 size(x,2)]);
L=permute(L,[1 3 2]);
CH=repmat((0:size(x,2)-1)*size(x,1),[size(L,1) 1 size(L,3)]);

X=x(L+CH); % matrix of segments [time x chanel x segment]

%%

X=X.*repmat(win,[1 size(X,2) size(X,3)]);

X=fft(X,nfft,1);
if nargout==1; return; end

% % The window is convolved with every power spectrum peak, therefore
% % compensate for the DC value squared to obtain correct peak heights.
% U=sum(win)^2;

% compensates for the power of the window.
U=(win(:)'*win(:))*(length(win)/2);  % normalization in accordance with cpsd() results

Pxy=zeros(size(x,2),size(x,2),size(X,1));

for ch=1:size(x,2)
   Pxy (ch,:,:)=permute(mean( (repmat(X(:,ch,:),[1 size(x,2) 1]).*conj(X))/U,3),[3 2 1]);
end

if nargout==2; return; end

F = fs*linspace(0,1-1/nfft,nfft);
%% coherence

Pxx=repmat( reshape(Pxy(logical(repmat(eye(size(Pxy,1)),[1 1 size(Pxy,3)]))), [1 size(Pxy,2) size(Pxy,3)] ), [size(Pxy,1) 1 1]);
Pyy=permute(Pxx,[2 1 3]);

Cxy=(abs(Pxy).^2)./(Pxx.*Pyy);


if nargout>3
   TFxy = Pxy./Pxx;
end

