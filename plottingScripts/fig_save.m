function fig_save(varargin)
% saves figure
% varargin{1} = figure handle
% varargin{2} = figure name
% varargin{3} = output directory
% further input arguments in varargin are tuples: attribute,value
% usage:
%   fig_save(f, 'figname', outDir, 'format',{'png','fig'}, 'res',600);

% (c) Jiri, Mar20

assert(size(varargin,2) >= 3);

%% figure handle
f = varargin{1};
assert(ishandle(f));

%% figure name
figname = varargin{2};
assert(ischar(figname));

%% output directory
outdir = varargin{3};
if ~exist(outdir, 'dir')
    mkdir(outdir);
end 

%% settings: resolution (res)
res = '-r0';   % default
for n = 4:size(varargin,2)
    if strcmp(varargin{n},'res')
        assert(isnumeric(varargin{n+1}));
        res = ['-r' num2str(varargin{n+1})];     % ~ '-r600'
    end
end

%% settings: output format (format)
outform = {'fig', 'png'};   % default
for n = 4:size(varargin,2)
    if strcmp(varargin{n},'format')
        outform = varargin{n+1};
    end
end

%% >>> save <<<
set(f, 'PaperPositionMode','auto');

% save as .fig
if ismember('fig',outform)
    saveas(f, [outdir filesep figname '.fig']);
end

% save as .png
if ismember('png',outform)
    print(f, '-dpng',res, [outdir filesep figname '.png']);
end

% save as .tiff
if ismember('tiff',outform)
    print(f, '-dtiff',res, [outdir filesep figname '.tiff']);
end

disp(['Figure: ' figname ', stored in: ' outdir]);
