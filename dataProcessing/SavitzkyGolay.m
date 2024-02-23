function fir = SavitzkyGolay(n,varargin)

% Version 1.4 - 31/03/00 - rotter@biologie.uni-freiburg.de
%
% fir = SavitzkyGolay([nl,nr]) gives a Savitzky-Golay smoothing filter
% of total width nr-nl+1. It extends -nl data points to the left and nr
% data points to the right of the point of operation. Works by least
% squares fitting of a polynomial to nr-nl+1 consecutive data points. A
% differentiating filter is obtained by specifying the <option, value>
% pair <'DerivativeOrder', ld> where ld=1 for the first derivative, ld=2
% for the second derivative, and so on. It is obtained by evaluating the
% ld-th derivative of the fitted polynomial at the position of the
% (1-nl)-th data point. The default is ld=0 for smoothing. The option
% <'TimeStep', h> can be used to account for any sampling rate. This
% effectively leads to a scaling of all filter coefficients by
% (1/h)^ld. The order of the fitted polynomial can be specified by
% <'PolynomialOrder', m>, the default is m=2 corresponding to a
% parabola. It is expected that m+1<=nr-nl+1. Data windowing can be
% enforced by <'WindowFunction', wf> where wf='Uniform' is the
% default, 'Welch' can also be specified. See Numerical Recipes in C
% (Second Edition) for further details.

%
% Set default values for all options.
%

if length(n) == 1
  nl = -n;
  nr = n;
else
  nl = n(1);
  nr = n(2);
end


h  = 1;
ld = 0;
m  = 2;
wf = 'Uniform';

for i = 1:2:length(varargin)
  switch varargin{i}
    case 'TimeStep'
      h = varargin{i+1};
    case 'DerivativeOrder'
      ld = varargin{i+1};
    case 'PolynomialOrder'
      m = varargin{i+1};
    case 'WindowFunction'
      wf = varargin{i+1};
    otherwise
      disp(['Unknown option ', varargin{i}, '.']);
  end
end

%
% Check parameters to avoid an underdetermined situation
%

if nr-nl+1 < m+1
  error('Value for PolynomialOrder is too high.');
end

%
% Generate list of weights.
%

switch wf
  case 'Uniform'
    w = ones(1,nr-nl+1);
  case 'Welch'
    w = 1 - [(nl:-1)/(nl-1) 0 (1:nr)/(nr+1)].^2;
  otherwise
    disp('Unknown value for option WindowFunction');
end

%
% This is the (weighted) design matrix of the corresponding
% least-squares problem.
%

A = zeros(nr-nl+1,m+1);

for j = 0 : m
  for i = nl : nr
    A(i-nl+1,j+1) = i^j * w(i-nl+1);
  end
end

%
% Finally, compute the filter coefficients.
%

A = pinv(A);
fir = prod(1:ld) / h^ld * A(ld+1,:) .* w;