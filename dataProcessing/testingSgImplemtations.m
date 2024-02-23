x = rawData(:,1);

%% sgolay
wSize = 251;
ord = 2;
tic; 
sgData = sgolayfilt(x, ord, wSize);
tSg=toc;

%% asg
wSize = 251;
ord = 2;
rOffset = 0.1*500;
tic;
%smoothKernel = SavitzkyGolay([-(wSize-rOffset-1),rOffset],'PolynomialOrder',ord);  !!! WRONG, has to be in reversed order???
smoothKernel = flipdim(SavitzkyGolay([-(wSize-rOffset-1), rOffset],'PolynomialOrder',ord), 2);
tmp = cat(1, flipdim(x(1:floor(wSize/2),:),1), x, flipdim(x(end-floor(wSize/2):end,:),1));
asgData = nan(size(x));
for t = 1:size(x,1)
    i_t = [1:wSize] + (t-1);
    asgData(t,:) = sum(tmp(i_t,:).*repmat(smoothKernel', 1,size(x,2)), 1);
end
tAsg=toc;

%% fast asg
tic;
fasgData = fsgwindow(x,ord,wSize,rOffset); % (c) S.Rotter
tFasg = toc;

%zeroLagFasgData = fsgwindow(x,ord,wSize,0); % (c) S.Rotter
halfWLagFasg = fsgwindow(x,ord,wSize,ceil(wSize/2)); % (c) S.Rotter

%% plots of different sg implementations
figure;
hold on;

plot(x, 'c')
plot(sgData, 'b');
plot(asgData, 'r');
plot(fasgData, 'g');
plot(halfWLagFasg, 'm');
legend({'raw', 'sg','asg','fasg', 'halfWLagFasg'});

%plot(zeroLagFasgData, 'c');
% legend({'sg','asg','fasg', 'zeroLagFasg', 'halfWLagFasg'});

%% test different freq. response of low/hi pass filters
loF = 2.0;
freqNyquist = 250;
Wn = loF/freqNyquist;           % normalized bandpass frequencies
n = [2:7];                      % butterworth order

a = cell(1,length(n));
b = cell(1,length(n));
for nn = 1:length(n)
    thisN = n(nn);
    [b{nn},a{nn}] = butter(thisN, Wn, 'low');  
end


fhandle = fvtool(b{1},a{1}, b{2},a{2}, b{3},a{3}, b{4},a{4}, b{5},a{5}, b{6},a{6}, 'Analysis','magnitude');
set(fhandle, 'Fs',params.amp.srate, 'NormalizedFrequency','off');
set(gca, 'xlim',[0 10]);
legend(fhandle, 'n 2', 'n 3', 'n 4', 'n 5', 'n 6', 'n 7');

