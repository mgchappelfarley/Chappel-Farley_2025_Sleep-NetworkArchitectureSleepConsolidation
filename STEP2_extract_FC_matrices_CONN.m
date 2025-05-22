%% Step 1. Extract FC matrix from CONN  - make sure to run STEP 1 prior to this
clear all; clc;
HomeDir = '/home/dir/';
CONNDir = '/dir/with/combined/conn/matrix'; 


fprintf('\n Step 1 Complete- Home and Conn output directories set \n');

%% Step 2. Load CONN ROI-to-ROI matrix (3D)
load([CONNDir,'/resultsROI_Condition001_n40_Combinedbatch1batch2scansBatch3.mat']); % FC matrix here
% matrix with ROI-to-ROI results is "Z" - number of ROIs (names) x number
% of ROIs  (names2) x number of subjects

numROIs = length(names); %Determine the number of ROIs
clear names2 regressors soyunDOF batch2SE xyz batch1DOF batch1SE batch3DOF batch3SE

fprintf('\n Step 2 Complete- Conn ROI file loaded \n');
    
%% Step 3. Extract ROI-to-ROI FC (3D matrix with all subjects)
    
%FCmatrix = Z(1:numROIs,1:numROIs,:); % for one scan type
FCmatrix = CombinedScans(1:numROIs,1:numROIs,:); % for combined scan types

% this pulls out just the ROIs that you want for all the subjects (last :)
% each layer of FCmatrix is another subject ROI-ROI correlation matrix for a specific
% subject 

% save files    
savDir = [HomeDir,'/FC_matrix_bivariate/'];
%mkdir(savDir);
save([savDir,'FCmatrix_151ROIs_N40.mat'],'FCmatrix','names', 'sublist', 'numROIs');
clear FCmatrix;

fprintf('\n Step 3 Complete- ROI-ROI FC matrix extracted and saved \n');

%% Step 4. Extract 2D matrix per subject
% this slices up the big matrix and extracts subject level data and saves it
for i = 1:length(sublist)
    subject = num2str(sublist(i));
    disp(['Running ', subject]);
    
    FCmatrix = CombinedScans(1:numROIs,1:numROIs,i); % all ROIs
    
   
    save([savDir,subject,'_FCmatrix_151ROIs.mat'],'FCmatrix','names', 'numROIs') % change number of ROIs in mat file name to save
    clear FCmatrix
  
end

%FCmatrix_XROIs_Nnumberof subs.mat is the FC matrix for all X subjects. IF
%you wanted average across everyone you would average across the third
%demention. The other outputs are the individual slices for each subjects.


fprintf('\n Step 4 Complete- 2D Matrix/ Subject has been extracted and saved \n');