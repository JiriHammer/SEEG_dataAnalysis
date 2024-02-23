function plotEC_plotband_2NN(FC, EC_info)
% plotband of effective coinnectivity (e.g. DTF) for different classes
% mean fover subj
% FC = 3D: [t, clz, subj], may inlcude NaNs

% (c) Jiri, Oct23

% plotband subj mean
for clz = 1:size(FC,2)
    plotband(EC_info.time, nanmean(FC(:,clz,:),3), sem(FC(:,clz,:),3), EC_info.info.colors(clz,:));
end
axis tight;

plot([0 0], ylim, '--k');

% significance
vals = reshape(FC, [size(FC,1), 1, size(FC,2)*size(FC,3)]); % 3D: [t, 1, subj-clz], where 1 = "ch", subj-clz = "trials"
labels = repmat([1,2], [1, size(FC,3)]);
[P_vals, H_vals] = getSignificance(vals, labels, EC_info.P_level);
plot2axes_signif_filled(EC_info.time, H_vals, P_vals);

% labels
xlabel('time (s)');
ylabel('DTF (z-score)');
title(EC_info.str_title);

grid on;
box on;
