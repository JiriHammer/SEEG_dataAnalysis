function loadAndFilter(params, dataStruct)  
% loads & filters data stored in session cache files
% required structure 'dataStruct' with fields (ex.):
% dataStruct =  
%            D_field: 'tracker'
%         signalType: 'xVel'
%      excludedChnls: []
%     dataProcessing: {'bp_loPass'  'norm2std'}
%           saveName: 'resp'
    
% (c) Jiri, Feb16

%% load & filter each session
disp(['Filtering of: ' dataStruct.saveName ' ...']);
inputVarName = [dataStruct.saveName '_raw'];
for sess = 1:size(params.storage.sessionCacheFiles,2)
    disp(['Processing session: ' num2str(sess) ' ...']);
    
    % load
    load(params.storage.sessionCacheFiles{sess}, inputVarName);
    assert(exist(inputVarName,'var') == 1);
    rawData = eval(inputVarName);
    clear(inputVarName);

    % filter: rawData -> filteredData
    params.usedSignalType = dataStruct.signalType;
    params.thisSess = sess;
    params.thisData = dataStruct.saveName;
    filteredData = filterData(params, rawData, dataStruct.dataProcessing, params.init_srate);

    % save
    outputVarName = [dataStruct.saveName '_prc'];
    eval([outputVarName '=filteredData;' ]); 
    disp([' - session: ' num2str(sess) ' done. Saving: ' outputVarName]);
    save(params.storage.sessionCacheFiles{sess}, outputVarName, '-append');      
    clear rawData outputVarName filteredData;
end    
