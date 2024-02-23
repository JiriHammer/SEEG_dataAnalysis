function [xs, x_2S] = FT_ampPhaseShuffle(x, dataInfo)
% manipulate amplitudes (magnitudes) or phases of a signal at specific freq
% - input: real, time-domain signal
% - computes FFT (on columns of x)
% - shuffles phases
% - iFFT (inverse FT)
% - output: modified, real, time-domain signal

% (c) Jiri, Jun21

%% settings
if ~isfield(dataInfo, 'shufflePerChnl')
    shufflePerChnl = true;      % default: diff for all chnls
else
    shufflePerChnl = dataInfo.shufflePerChnl;
end
srate = dataInfo.srate;     % sampling rate, in [Hz]
selFreq = dataInfo.selFreq;     % selected frequencies, in [Hz]

%% ========== FT ===================
X= fft(x);      % complex-valued signal
Ax= abs(X);     % amps
PHIx= angle(X); % phases

%% --- modify phases
if ~shufflePerChnl
    % shuffle phases of different frequencies 
%     randphase= randperm(size(PHIx,1));
%     while any(1:size(PHIx,1)==randphase)
%         randphase= randperm(size(PHIx,1));
%     end
%     PHIx_rand= PHIx(randphase,:);           % same for all chnls, random phases, 2D: freq x ch

    PHIx_rand_ch = rand(size(PHIx,1),1)*2*pi-pi;            % random phases for all freq
    PHIx_rand= repmat(PHIx_rand_ch, [1, size(PHIx,2)]);     % same for all chnls, random phases, 2D: freq x ch
else
    PHIx_rand= rand(size(PHIx))*2*pi-pi;    % diff for all chnls, random phases <-pi,pi>, 2D: freq x ch
end

%% ================== shuffle one-sided FFT at selected freq ==============
% --- freq. vector
L = size(x,1);
F = srate*(0:(L/2))/L;

% --- selected freq: one sided
i_freq = closestval(F, selFreq(1)):closestval(F, selFreq(2));

% ---- new data with shuffled phases at selected freq: single sided
XS = X;
XS(i_freq,:)= Ax(i_freq,:) .* exp(1i*PHIx_rand(i_freq,:));

% --- iFT (real part)
xs= real(ifft(XS));

%% ===== new data with shuffled phases at ALL freq ====
XS_surr= Ax .* exp(1i*PHIx_rand);

% --- iFT (real part)
x_surr= real(ifft(XS_surr));

%% --- FFT again (do not have the same spectra ???)
% X_surr= fft(x_surr);
% A_surr= abs(X_surr);
% 
% figure; hold on;
% ch = 21;
% plot(Ax(:,ch), 'r');
% plot(A_surr(:,ch), 'b');
% M_allFreq = corrcoef(Ax(:,ch), A_surr(:,ch))


%% ================== shuffle two-sided FFT at selected freq ==============
% --- fftshift -> center DC to middle of the spectrum
% !!! The two-sided FFT is as discussed above from DC-->Fmax-->DC with negative frequency in the first location in the returned vector.
Y = fftshift(X);
A_2S= abs(Y);     % amps

% --- freq. axis: two sided freq spectrum
F_Nyuq = srate/2;
F_2S = linspace(-F_Nyuq, F_Nyuq, size(X,1));

% ---- selected freq: two sided
i_freq_N = closestval(F_2S, -selFreq(2)):closestval(F_2S, -selFreq(1)); % negative freq
i_freq_P = closestval(F_2S, selFreq(1)):closestval(F_2S, selFreq(2));   % positive freq

% ---- new data with shuffled phases at selected freq: single sided
Y_2S = Y;
Y_2S(i_freq_N,:)= A_2S(i_freq_N,:) .* exp(1i*PHIx_rand(i_freq_N,:));
Y_2S(i_freq_P,:)= A_2S(i_freq_P,:) .* exp(1i*PHIx_rand(i_freq_P,:));   % complex conjugate???

% --- iFT (2 sided)
X_2S = ifftshift(Y_2S);
x_2S = real(ifft(X_2S));


1+1;

%% uncomment to plot
% t = [1:size(x,1)]'./srate;
% figure; hold on;
% ch = 1;
% plot(t,x(:,ch), 'r');
% plot(t,xs(:,ch), 'b');
% plot(t,x_surr(:,ch), 'm');
% plot(t,x_2S(:,ch), 'k');
% 
% 1+1;

