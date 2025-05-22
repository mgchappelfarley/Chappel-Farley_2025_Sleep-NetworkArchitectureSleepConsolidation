%% Combine rsfMRI connectivity matrices for different projects (i.e., 3-min batch1 scans and 10-min batch2 whole brain)

%% Step 1 - set directories
HomeDir = '/name/of/dir/';
CONNDir_scansBatch1 = '/path/to/conn/first-level/results/batch1/'
CONNDir_scansBatch2 = '/path/to/conn/first-level/results/batch2/' 
CONNDir_scansBatch3 = '/path/to/conn/first-level/results/batch3/'
%% Step 2 - Load in files and save subjects matrices with identifier for scan type

load([CONNDir_scansBatch1,'/resultsROI_Condition001.mat']); % FC matrix here
scansBatch1ubjects = Z; 
batch1DOF = DOF;
batch1SE = SE;
clear DOF names names2 regressors SE xyz Z; % keep workspace tidy

load([CONNDir_scansBatch2,'/resultsROI_Condition001.mat']); 
batch2DOF = DOF;
batch2SE = SE;
scansBatch2ubjects = Z; % subjects with the batch2 scan resting state
clear SE DOF; 

load([CONNDir_scansBatch3,'/resultsROI_Condition001.mat']); 
batch3DOF = DOF;
batch3SE = SE;
scansBatch3ubjects = Z; % subjects with the batch2 scan resting state
clear SE DOF; 

%% Step 3 - combine the three 3D matrices containing the connectivity matrices for each scan type
CombinedScans = cat(3, scansBatch1ubjects, scansBatch2ubjects, scansBatch3ubjects);

%% Step 4 - create subject list
cd(HomeDir);
sub_info = readtable('list_of_subject_ids.csv');
sublist = sub_info.beacon_id;

%% Step5 - Save combined file
cd(HomeDir);
save('resultsROI_Condition001_n40_Combinedbatch1batch2scansBatch3.mat','CombinedScans', 'sublist', 'batch1DOF', 'batch1SE','batch2DOF','batch2SE','batch3DOF', 'batch3SE','names', 'names2', 'regressors', 'xyz');

