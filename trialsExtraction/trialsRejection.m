function trials_out = trialsRejection(trials_in)
% rejects (leaves out) trials marked in

% (c) Jiri, Nov22

%% copy the entire structure
trials_out = trials_in;

%% time of interest (TO DO: pass from params)
toi = [-2,2];
i_t = closestval(trials_in.time,toi(1)):closestval(trials_in.time,toi(end));

%% modify trials.data
nTr_ok = 0;
d_ok = [];
r_ok = [];
c_ok = [];
for tr = 1:size(trials_in.data,3)
    r_tr = trials_in.rejected(i_t,:,tr);
    if ~any(r_tr(:))
        d_ok = cat(3, d_ok, trials_in.data(:,:,tr));
        r_ok = cat(3, r_ok, trials_in.rejected(:,:,tr));
        c_ok = cat(2, c_ok, trials_in.labels(tr));
        nTr_ok = nTr_ok + 1;
    end
end

%% overwrite trials
trials_out.data = d_ok;
trials_out.rejected = r_ok;
trials_out.clzNames = c_ok;

        
        
    
