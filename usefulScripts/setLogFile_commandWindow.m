function setLogFile_commandWindow(params)
% sets log file for output that appears in command window
% file = 'log_commandWindow.txt' in output directory
% copies existing (empty) 'log_commandWindow.txt' from usefulScripts

% (c) Jiri, Oct21

%% source file
diaryFileName_src = [params.storage.dir_neuroHammer filesep 'usefulScripts' filesep 'log_commandWindow.txt'];
assert(exist(diaryFileName_src,'file') == 2);

%% destination file
diaryDirName_dst = [params.storage.dir_results filesep params.storage.outName];
if ~exist(diaryDirName_dst,'dir')
    mkdir(diaryDirName_dst);
end
diaryFileName_dst = [diaryDirName_dst filesep 'log_commandWindow.txt'];

%% copy
copyfile(diaryFileName_src, diaryFileName_dst);
assert(exist(diaryFileName_dst,'file') == 2);

%% set logging on
diary(diaryFileName_dst);
if strcmp(get(0,'Diary'),'off')
    diary on;
end
