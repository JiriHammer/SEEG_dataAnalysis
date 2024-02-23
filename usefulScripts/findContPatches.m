function epochs = findContPatches(dataInds)
% localizes beginnings and ends of continouus patches in vector of data indices

% (c) Jiri, Oct16

% include discrete edges
dataInds = dataInds(:);     % column vector
d = cat(1, -666, dataInds, max(dataInds)+10);

% determine jumps in data index continuity
i_jump = find(diff(d)>1);

% return beginning / end of each conti. segment
epochs(:,1) = d(i_jump(1:end-1)+1);
epochs(:,2) = d(i_jump(2:end));

