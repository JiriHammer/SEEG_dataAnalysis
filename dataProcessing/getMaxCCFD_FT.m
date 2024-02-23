function fd_ft = getMaxCCFD_FT(predTraj, respTraj, realTraj, parms)
% computes windowed CC between predTraj and respTraj, both have format:
% predTraj: samples x FD groups (freq)
% and returns its maximum freq. -> ccMax
% computes windowed FT, smoothes and returns its maximum freq. -> ftMax

% (c) Jiri, Jul12

nSamples = size(respTraj,1);
assert(size(predTraj,1) == nSamples);
assert(size(predTraj,2) == size(respTraj,2));
assert(size(realTraj,1) == nSamples);

%% settings
wS_samples = ceil(parms.wSize * parms.srate);
tS_samples = ceil(parms.tStep * parms.srate);

%% time-resolved windowing, FT max
doMaxFT = false;
if doMaxFT
    ccMax = [];
    ftMax = [];
    currWindow = 1:wS_samples;
    while currWindow(end) <= nSamples

        % CC between predTraj and respTraj, returns its maximum freq. -> ccMax
        ccFD = nan(size(predTraj,2),1);
        for group = 1:size(predTraj,2)
            cc = corrcoef(predTraj(currWindow,group), respTraj(currWindow,group));
            ccFD(group) = cc(1,2);
        end
        [m, i_m] = max(ccFD);
        ccMax = cat(1, ccMax, parms.freqFD(i_m));

        % FT, smoothes and returns its maximum freq. -> ftMax
        y = realTraj(currWindow);
        NFFT = 2^nextpow2(size(y, 1));                  % number of FFT points
        fftData = fft(y, NFFT)/(0.5*size(y, 1));        % normalizes to windowSize = all samples
        yFT = abs(fftData(1:NFFT/2+1,:,:));
        freqAxis = parms.srate/2 * linspace(0,1,NFFT/2+1);

    %     Wn = 5/(parms.srate/2);                         % smooth with 5 Hz low pass
    %     [b,a] = butter(3, Wn, 'low');
    %     bpSmoothFFT = filtfilt(b,a, yFT);  
    %     f = figure; hold on;
    %     plot(freqAxis, filtFT, 'b');
    %     plot(freqAxis, bpSmoothFFT, 'r');
    %     close(f);
    %     yFT = bpSmoothFFT;

        % select maximum
        i_f = closestval(freqAxis, parms.freqFD(1)):closestval(freqAxis, parms.freqFD(end));
        yFT_selFreq = yFT(i_f);
        frAx = freqAxis(i_f);
        [m, i_m] = max(yFT_selFreq);
        ftMax = cat(1, ftMax, frAx(i_m));

        % shift window by time step
        currWindow = currWindow + tS_samples;
    end
    fd_ft = [ccMax, ftMax];
end

%% maximum variance explained by traj. FC
doMaxVarExplained = false;
if doMaxVarExplained
    ccMax = [];
    ftMax = [];    
    currWindow = 1:wS_samples;
    while currWindow(end) <= nSamples

        % CC between predTraj and respTraj, returns its maximum freq. -> ccMax
        ccFD = nan(size(predTraj,2),1);
        for group = 1:size(predTraj,2)
            cc = corrcoef(predTraj(currWindow,group), respTraj(currWindow,group));
            ccFD(group) = cc(1,2);
        end
        [m, i_m] = max(ccFD);
        ccMax = cat(1, ccMax, parms.freqFD(i_m));

        % FT, smoothes and returns its maximum freq. -> ftMax
        y = realTraj(currWindow);
        NFFT = size(y, 1);                  % number of FFT points
        freqAxis = parms.srate/2 * linspace(0,1,NFFT/2+1);

        fcUsed = closestval(freqAxis, parms.freqFD(1)):closestval(freqAxis, parms.freqFD(end));
        nFC = closestval(freqAxis, parms.freqFD(end));
        
        fftData = fft(y);           % !! FT !!

        currTraj = y;
        currVar = nan(nFC+1,1);
        currVar(1) = var(y, 0, 1);

        % !!!! reconstruct traj from individual freq. components: fcTraj !!!!
        fc = 1;             % DC component
        i_f = fc;
        fcData = zeros(length(fftData),1);
        fcData(i_f) = fftData(i_f);
        fcTraj = ifft(fcData);
        currTraj = currTraj - fcTraj;
        currVar(2) = var(currTraj, 0, 1);

        % other freq. components
        for fc = 2:nFC
            i_f = [fc, length(fftData)-fc+2];
            fcData = zeros(length(fftData),1);
            assert(fftData(i_f(1)) == conj(fftData(i_f(2))));
            fcData(i_f) = fftData(i_f);                     % iFFT vector
            fcTraj = ifft(fcData);                          % !! iFT !!
            assert(isreal(fcTraj));

            currTraj = currTraj - fcTraj;
            currVar(fc+1) = var(currTraj, 0, 1);
        end
        fcVarExplain = diff(currVar)./currVar(1)*(-1);      % percentage of variance explained

        % smooth FFT result
%         Wn = 1/10;
%         [b,a] = butter(3, Wn, 'low');    
%         yFT = filtfilt(b,a, fcVarExplain);
        yFT = fcVarExplain;

        %figure; hold on; plot(freqAxis(1:nFC), fcVarExplain, 'b'); plot(freqAxis(1:nFC),yFT, 'r');
        
        % select maximum
        i_f = closestval(freqAxis, parms.freqFD(1)):closestval(freqAxis, parms.freqFD(end));
        yFT_selFreq = yFT(i_f);
        frAx = freqAxis(i_f);
        [m, i_m] = max(yFT_selFreq);
        ftMax = cat(1, ftMax, frAx(i_m));

        % shift window by time step
        currWindow = currWindow + tS_samples;
    end
    fd_ft = [ccMax, ftMax];        
end

%% matrix product of ccFD & yFT at each time step
doMatrixProduct = true;
if doMatrixProduct
    fd_ft = [];
    currWindow = 1:wS_samples;
    while currWindow(end) <= nSamples

        % CC between predTraj and respTraj, returns its maximum freq. -> ccMax
        ccFD = nan(size(predTraj,2),1);
        for group = 1:size(predTraj,2)
            cc = corrcoef(predTraj(currWindow,group), respTraj(currWindow,group));
            ccFD(group) = cc(1,2);
        end

        % FT, smoothes and returns its maximum freq. -> ftMax
        y = realTraj(currWindow);
        NFFT = size(y, 1);                  % number of FFT points
        freqAxis = parms.srate/2 * linspace(0,1,NFFT/2+1);
        nFC = closestval(freqAxis, parms.freqFD(end));
        
        fftData = fft(y);           % !! FT !!

        currTraj = y;
        currVar = nan(nFC+1,1);
        currVar(1) = var(y, 0, 1);

        % !!!! reconstruct traj from individual freq. components: fcTraj !!!!
        fc = 1;             % DC component
        i_f = fc;
        fcData = zeros(length(fftData),1);
        fcData(i_f) = fftData(i_f);
        fcTraj = ifft(fcData);
        currTraj = currTraj - fcTraj;
        currVar(2) = var(currTraj, 0, 1);

        % other freq. components
        for fc = 2:nFC
            i_f = [fc, length(fftData)-fc+2];
            fcData = zeros(length(fftData),1);
            assert(fftData(i_f(1)) == conj(fftData(i_f(2))));
            fcData(i_f) = fftData(i_f);                     % iFFT vector
            fcTraj = ifft(fcData);                          % !! iFT @ single freq !!
            assert(isreal(fcTraj));

            currTraj = currTraj - fcTraj;
            currVar(fc+1) = var(currTraj, 0, 1);
        end
        yFT = diff(currVar)./currVar(1)*(-1);      % percentage of variance explained
        
        fd_ft_step = ccFD * yFT';                   % matrix product
        %figure; hold on; imagesc(fd_ft_step);

        % shift window by time step
        fd_ft = cat(3, fd_ft, fd_ft_step);
        currWindow = currWindow + tS_samples;
    end
end

