function colors = colorPalette(n)
% defines colors for n elements as rgb tripplets

% (c) Jiri, May12

if n < 4
    colors = [0 0 1;        % b
              1 0 0;        % r
              0 1 0];       % g
else
    func = @(x) colorspace('RGB->Lab',x);
    colors = distinguishable_colors(n,'w',func);
end