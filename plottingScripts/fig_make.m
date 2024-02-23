function f = fig_make(varargin)
% creates figure based on aguments in varargin{1}, ..., varargin{end}
% returns its handle 'f'
% default (no arguments): full screen figure
% usage examples: 
% f = fig_make('fig_W',18.3, 'fig_H',14);
% f = fig_make;

% (c) Jiri, Mar20

%% settings: backgroud color (bkgclr)
bkgclr = 'w';   % default
% bkgclr = 'none';   % default
for n = 1:size(varargin,2)
    if strcmp(varargin{n},'bkgclr')
        bkgclr = varargin{n+1};
    end
end

%% settings: units
fig_units = 'centimeters';   % default
for n = 1:size(varargin,2)
    if strcmp(varargin{n},'units')
        fig_units = varargin{n+1};
    end
end

%% settings: figure height (fig_H)
fig_H = 12;   % default, in [cm]
for n = 1:size(varargin,2)
    if strcmp(varargin{n},'fig_H')
        fig_H = varargin{n+1};
    end
end

%% settings: figure width (fig_W)
fig_W = 18;   % default, in [cm]
for n = 1:size(varargin,2)
    if strcmp(varargin{n},'fig_W')
        fig_W = varargin{n+1};
    end
end

%% >>> create figure <<<
if size(varargin,1) == 0
    f = figure('units','normalized','outerposition',[0 0 1 1],'color',bkgclr);
else
    f = figure('units',fig_units,'outerposition',[1 1 fig_W fig_H], 'color',bkgclr);
end

% preserve black background ?
if strcmp(bkgclr, 'k')
    set(f, 'InvertHardcopy','off');                 % preserves black background
end
