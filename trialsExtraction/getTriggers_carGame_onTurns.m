function getTriggers_carGame_onTurns(params)
% extracts discrete turns from car game paradigm

% (c) Jiri, Oct16

for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load road
    load(params.storage.sessionCacheFiles{sess}, 'road','time_raw');
    assert(size(road,1) == size(time_raw,1));

    % road derivative: "instrcuted velocity"
    road_raw = road(:,1);
    params.thisSess = sess;
    params.sgolay.wSize = 0.25;
    road_der = filterData(params, road_raw, {'sgolay','derivative'}); 
    
    % extract local maxima => turns
    road_der_abs = abs(road_der);
    thr = 200;      % in [pixels/s], threshold defining a turn
    i_epochs = findContPatches(find(road_der_abs > thr));
    t_turn_instructed = [];
    v_turn_instructed = [];
    for ev = 1:size(i_epochs,1)
        absRoadDer_epoch = road_der_abs(i_epochs(ev,1):i_epochs(ev,2));
        [m, i_m] = max(absRoadDer_epoch);
        i_t = i_epochs(ev,1)+i_m-1;
        t_turn_instructed = cat(1, t_turn_instructed, time_raw(i_t));
        v_turn_instructed = cat(1, v_turn_instructed, m);
    end
      
    % find "good trials"
end
    

%% debug plots
figure; hold on;
plot(time_raw, road_raw, 'b');
plot(time_raw, road_der_abs, 'r');
plot(get(gca,'xlim'), [thr, thr], '--k');
for ev = 1:size(i_epochs)
    plot([time_raw(i_epochs(ev,1)),time_raw(i_epochs(ev,1))], get(gca,'ylim'), '--k');
    plot([time_raw(i_epochs(ev,2)),time_raw(i_epochs(ev,2))], get(gca,'ylim'), 'k');
end
