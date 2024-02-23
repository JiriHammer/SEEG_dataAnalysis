function xs = surrogate_data(x,shuffle,smoothspect)
% function surrogate_data - produces random data with a specific auto-correlation
% 
% syntax: xs = surrogate_data(x)
%         xs = surrogate_data(x,shuffle,smoothspect)
%
% Input:
% x:  template data
% xs: surrogate data, produced by phase-shuffling in x, which has the
%     same autocorrelation (i.e. the same power-spectrum) as x.
% shuffle: if set to 1, the phases are actually shuffled (default, i.e. if
%          left out), if set to 0, the new phase-values are just random values
%          between -pi and pi.
% smoothspect: if set to 1 (default is 0), amplitude spectrum is smoothed,
%              in order to only reflect the more general properties of x in xs.
% xs is a matrix of the same size as x, where x is interpreted column-wise 
% (if x is not a vector).

% (c) Tobias Pistohl

if nargin<3
    smoothspect= 0;
end

if nargin<2 || isempty(shuffle)
    shuffle= 1;
end

rowinput= false;
if isvector(x)
    if size(x,1)==1
        rowinput= true;
    end
    x= x(:);
end

% determine the mean and standard-deviation of x to scale the output accordingly
stdx= std(x,[],1);
meanx= mean(x,1);
% Fourier transform of template data in x
X= fft(x);
Ax= abs(X); PHIx= angle(X);
clear X
% produce random phase
if shuffle
    randphase= randperm(size(PHIx,1));
    while any(1:size(PHIx,1)==randphase)
        randphase= randperm(size(PHIx,1));
    end
    PHIx= PHIx(randphase,:);
else
    PHIx= rand(size(PHIx))*2*pi-pi;
end
if smoothspect
    figure
    semilogy(Ax(:,1),'b')
    [b,a]= butter(1,0.2);
    Ax= [Ax(1,:);filtfilt(b,a,Ax(2:end,:))];
%     Ax= [Ax(1,:);sgolayfilt(Ax(2:end,:),3,floor(size(Ax,1)/400)*2+1)];
    hold on
    semilogy(Ax(:,1),'r')
end
XS= Ax .* exp(1i*PHIx);
clear Ax PHIx
xs= real(ifft(XS));

% correct for variance and mean
stdxs= std(xs,[],1);
xs= detrend(xs.*repmat(stdx./stdxs,size(x,1),1),'constant')+repmat(meanx,size(x,1),1);

if rowinput
    xs= xs';
end
