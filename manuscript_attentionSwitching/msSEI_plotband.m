function msSEI_plotband(x, y_avg, y_sem, clz)
% plotband of SEI results
% E->I = b->m
% I->E = r->c

% (c) Jiri, Nov23

%% colors
% clrs = {
%     'b','m';
%     'r','c';
%     };

% clrs = {
%     'b','r';
%     'r','b';
%     };

% clrs = {
%     [0 0 255]./255, [255 0 255]./255; % b -> m
%     [255 0 0]./255, [255,140,0]./255; % r -> o
% };

clrs = {
    [0 0 255]./255, [255 0 255]./255;   % dark b -> m
    [255 0 0]./255, [81,174,255]./255;  % r -> light b
};

%% index at t = 0
i_zero = closestval(x, 0.0);

%% plotband (change color at t = 0)
% plot from beg -> 0 (= task switch)
plotband(x(1:i_zero), y_avg(1:i_zero), y_sem(1:i_zero), clrs{clz,1});

% plot from 0 -> end
plotband(x(i_zero:end), y_avg(i_zero:end), y_sem(i_zero:end), clrs{clz,2});

% plot from 0 -> end
% if i_zero == size(x,1)
%     return;
% end
% plotband(x(i_zero+1:end), y_avg(i_zero+1:end), y_sem(i_zero+1:end), clrs{clz,2});


%% colors as characters or RGB code?
% if ischar(clrs{1,1})
%     % plot from beg -> 0 (= task switch)
%     plotband(x(1:i_zero), y_avg(1:i_zero), y_sem(1:i_zero), clrs{clz,1});
% 
%     % plot from 0 -> end
%     if i_zero == size(x,1)
%         return;
%     end
%     plotband(x(i_zero+1:end), y_avg(i_zero+1:end), y_sem(i_zero+1:end), clrs{clz,2});
% else
%     % plot from beg -> 0 (= task switch)
%     plotband(x(1:i_zero), y_avg(1:i_zero), y_sem(1:i_zero), clrs{clz,1});
% 
%     % plot from 0 -> end
%     if i_zero == size(x,1)
%         return;
%     end
%     plotband(x(i_zero+1:end), y_avg(i_zero+1:end), y_sem(i_zero+1:end), clrs{clz,2});
% end
