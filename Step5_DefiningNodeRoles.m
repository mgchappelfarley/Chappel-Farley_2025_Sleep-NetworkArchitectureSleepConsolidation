% Miranda Chappel-Farley
% The purpose of this script is to calculate several node-specific
% mesasures including:

% Node Strength

% Eigenvector Centrality 

% Participation coefficient 

% Betweenness Centrality

%% Load in directories
clear all; clc;

HomeDir = '/home/dir/';
FCDir = [HomeDir,'/FC_matrix_bivariate'];
BinDir = [FCDir, '/WeightedMatrix/'];
CommDir = [BinDir, 'BCT_output_Step3_NetworkModularity_Weighted/'];

cd(HomeDir);
sub_info = readtable('list_of_subject_ids.csv');
subsToExclude = sub_info.to_exclude;

cd(FCDir);
sublist = load('FCmatrix_151ROIs_N40.mat', 'sublist');
sublist = sublist.sublist;
numSubs = length(sublist);

% Select diretory with input data (symmetric matrices)
sym_matrix_dir = BinDir;
cd(BinDir);
% Set generic file path
sym_file_path = fullfile(sym_matrix_dir, '*BCTweighted.mat');
% Create object with all the matrix files
files = dir(sym_file_path);
clear sym_file_path

% Load in the first file to automatically set the number of nodes 
load(files(1).name, 'A_sym_norm');
adj_matrix = A_sym_norm;
clear A_sym_norm;
numNodes = size(adj_matrix, 1);
clear adj_matrix ;
%% Calculate Node Strengths
% -----
% Node strength (weighted network) - the sum of connectivity weights of the
% edges attached to each node - can compute positive and negative strength
% Node strength is the sum of weights of links connected to the node.
% -----
cd(BinDir);
TotalNodeWeight = zeros(numNodes,2,numSubs); 
PosStrength = zeros(numNodes, numSubs);
NegStrength = zeros(numNodes, numSubs);
NormStrength = zeros(numNodes, numSubs);



for i = 1:length(files)
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating node strengths for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A = A_sym_norm; 
    [Spos, Sneg, vpos, vneg] = strengths_und_sign(A);
    PosStrength(:,i) = Spos;
    NegStrength(:,i) = Sneg;
    Snorm = Spos - ((Sneg/Spos + Sneg).*Sneg);
    % normalized node strength, Rubinov & Sporns 2011
    % s* = s+ - (s-/s+ + s-)s-
    NodeStrengths(:,1, i) = vpos;
    NodeStrengths(:,2, i) = vneg;
    clear A 
end

fprintf('Node strengths calculated and compiled for all subjects: saved in NodeStrengths \n');
%% Plot the distributions of node strengths
figure
subplot(3,1,1)
histogram(PosStrength(:,:))
xlabel('Nodal Strength for Positive Weights')
ylabel('Count')
subplot(3,1,2)
histogram(NegStrength(:,:))
xlabel('Nodal Strength for Negative Weights')
ylabel('Count')
subplot(3,1,3)
histogram(NormStrength(:,2))
xlabel('Normalized Nodal Strength')
ylabel('Count')

figure
subplot(2,1,1)
histogram(TotalNodeWeight(:,1,:))
xlabel('Total Positive Weight')
ylabel('Count')
subplot(2,1,2)
histogram(TotalNodeWeight(:,2,:))
xlabel('Total Negative Weight')
ylabel('Count')


%% Eigenvector Centrality
%Eigenvector centrality accounts for the quantity and quality of a nodes
%degree -- can compute with weighted networks, either add 1 to negative
%weights, but then strong negative correlation will be given lower weight
%than a weak pos corr. If you take abs val. a strong neg correlation will
%have equal weighting as a strong pos. 

% choosing to remap by adding 1, since the meaning behind neg correlation
% is unclear, but don't want to disregard them entirely or give them teh
% same weight. (Lohmann 2010: Eigenvector centrality mapping for analyzing
% connectivity patterns in fmri data of the human brain)

evCentrality = zeros(numSubs, numNodes);
B = ones(numNodes, numNodes); % create a matrix of ones with the same size as our adj matrix
for i = 1:length(files)
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating eigenvector centrality for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A = A_sym_norm; 
    Apos = A + B; % add 1 to everything so no negatives, but not entirely ignoring them
    v = eigenvector_centrality_und(Apos);
    evCentrality(i,:) = v;
end

histogram(evCentrality(:,:)) % check distribution
xlabel('Eigenvector Centrality')
ylabel('Count')
fprintf('Eigenvector Centrality calculated and compiled for all subjects: saved in evCentrality \n');
%% Participation coefficient - requires module assignment from modularity maximization (Step 3)
mod_dir = CommDir;
% Set file path to the partitions
cd(CommDir);
mod_file_path = fullfile(mod_dir, '*FinalLouvainPartitions.mat');
% Create object with all the matrix files
mod_files = dir(mod_file_path);
load(mod_files(1).name, 'gamma');
% pre-alllocate matrix to hold participation coefficients for each level of
% gamma
levelsOfgamma = length(gamma);

Pos_pCoef = zeros(numNodes,levelsOfgamma, numSubs);
Neg_pCoef = zeros(numNodes,levelsOfgamma, numSubs);

for i = 1:length(mod_files)
    fprintf(1, 'Calculating Participation Coefficient for %s\n', mod_files(i).name);
    load(mod_files(i).name); %load in the final partition matrix for each level of gamma
    cd(BinDir);
    load(files(i).name) % load in the weighted connectvity matrix
    cd(CommDir);
    communities = Dfinal; % the community affiliation vectors at each level of gamma
    W = A_sym_norm; % undirected connection matrix with pos and neg weights
    for i_gamma = 1:levelsOfgamma
        Ci = communities(:,i_gamma); % community affiliation vector
        [Ppos, Pneg] = participation_coef_sign(W,Ci);
        Pos_pCoef(:,i_gamma,i) = Ppos;
        Neg_pCoef(:,i_gamma,i) = Pneg;
    end
    
end

figure
subplot(2,1,1)
histogram(Pos_pCoef(:,1));
xlabel('Participation Coeff for Positive Weights')
ylabel('Count')
subplot(2,1,2)
histogram(Neg_pCoef(:,1)); % chceking distribution
xlabel('Participation Coeff for Negative Weights')
ylabel('Count')



fprintf('Positive and Negative Participation Coeffs calculated and compiled for all subjects: saved in Pos_pCoef and Neg_pCoef \n');
 
%% Load in ROI labels for analyses
maskDir = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/AnalysisMask/n40_AnalysisMask'
cd(maskDir);
ROIlabels = readtable('n40_finalMask_ROINumsAndLabels.csv');
%% SAVE to mat file 
save([BinDir,'NodeRoles_n40.mat'],'numNodes','sublist', 'numSubs', 'PosStrength', 'NegStrength', 'NormStrength', 'TotalNodeWeight', 'evCentrality', 'Pos_pCoef', 'Neg_pCoef' , 'ROIlabels');

%save recalculated pos participation coefficient - EC did not need to be
%recalc with new range of gamma
%save([BinDir,'ParticipationCoefficients_n40.mat'],'numNodes','sublist', 'numSubs', 'gamma', 'levelsOfgamma','Pos_pCoef', 'Neg_pCoef' , 'ROIlabels');
fprintf('Participation Coefficient .mat file generated');
%% Pull values for specific ROIs
% Subcortical ROIs in mask
hipp = find(contains(ROIlabels.regionLabel, 'Hipp'));
amg = find(contains(ROIlabels.regionLabel, 'Amyg'));


%vmPFC ROIs in mask
%A10
a10 = find(contains(ROIlabels.regionLabel, 'A10'));

%Control ROIs in mask: Occipital lobe
vmPos = find(contains(ROIlabels.regionLabel, 'vmPOS'));

%Entorhinal Cortex and temporal lobe control (parahippocampal)
a35 = find(contains(ROIlabels.regionLabel, 'A35')); %parahippocampal

%Compile all indices of ROIS
ROIS = vertcat(hipp, amg,vmPos, a10, a35);
ROIS = sort(ROIS); 

centrality_ROIS = evCentrality(:,ROIS);
Spos_ROIS = PosStrength(ROIS,:)';
Pcoef_ROIS = Pos_pCoef(ROIS,:,:);

%grab the least similar calc participation coefficients according to NMI
Pcoef_ROIS_gamma1 = squeeze(Pcoef_ROIS(:,1,:))'; 
Pcoef_ROIS_gamma4 = squeeze(Pcoef_ROIS(:,4,:))'; 

%pick median gamma value
Pcoef_ROIS_medianGamma = squeeze(Pcoef_ROIS(:,3,:))';
selectedROIS = ROIlabels(ROIS,:); %table of selected ROIS with labels
%% SAVE to mat file 

labels = selectedROIS.regionLabel'; % pull labels

centrality = array2table(centrality_ROIS);
centrality.Properties.VariableNames = labels; 
centrality.beacon_id = sublist;

Spos = array2table(Spos_ROIS);
Spos.Properties.VariableNames = labels;
Spos.beacon_id = sublist; 

Pcoef_mg = array2table(Pcoef_ROIS_medianGamma);
Pcoef_mg.Properties.VariableNames = labels;
Pcoef_mg.beacon_id = sublist;

% Pcoef_g1 = array2table(Pcoef_ROIS_gamma1);
% Pcoef_g1.Properties.VariableNames = labels;
% Pcoef_g1.beacon_id = sublist;
% 
% Pcoef_g4 = array2table(Pcoef_ROIS_gamma4);
% Pcoef_g4.Properties.VariableNames = labels;
% Pcoef_g4.beacon_id = sublist;

%To save just participation coeff
%save([BinDir,'ParticipationCoeff_37ROIs_n40.mat'],'sublist', 'Pcoef_mg', 'selectedROIS');

save([BinDir,'NodeRoles_37ROIs_n40.mat'],'sublist', 'numSubs','centrality', 'Spos', 'Pcoef_mg', 'selectedROIS');
fprintf('Node Roles for selected ROIS saved in .mat file: NodeRoles_37ROIs_n40.mat ');
%% Write to csv files for analyses
cd(HomeDir);
% convert total group Eff values to csv
writetable(centrality, 'EigenVectorCentrality_n40_37ROIs.csv');
writetable(Spos, 'PositiveNodeStrength_n40_37ROIs.csv');
writetable(Pcoef_mg, 'ParticipationCoefficient_medianGamma_n40_37ROIs.csv');

%% Betweenness Centrality -- a measure of information flow. Determines how central a vertex is for information propogation
% Data must be converted to connection-length matrix. Nodes with a high
% betweenness centrality participate in a large number of shortest paths.
% 
% The input matrix must be a connection-length matrix, typically
% obtained via a mapping from weight to length. For instance, in a
% weighted correlation network higher correlations are more naturally
% interpreted as shorter distances and the input matrix should
% consequently be some inverse of the connectivity matrix. 
% 
% Betweenness centrality may be normalised to the range [0,1] as
% BC/[(N-1)(N-2)], where N is the number of nodes in the network.

betweenness_centrality = zeros(numSubs, numNodes); % pre-allocate a matrix for our metrics
cd(BinDir);
for i = 1:length(files)
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating betweenness centrality for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A_sym_pos = max(A_sym_norm,0); % keeps everything zero and above (i.e., drops negative corrs)
    A_sym_pos_len = weight_conversion(A_sym_pos, "lengths"); %convert corrs to lengths. high corr = short length, low corr = long length
    bc = betweenness_wei(A_sym_pos_len); % calculate betweenness centrality on the connection-length metrix
    bc_norm = bc/((numNodes-1)*(numNodes-2)); %normalize to [0,1]
    betweenness_centrality(i,:) = bc_norm;
end
clear i;


betweenness_ROIs = betweenness_centrality(:,ROIS); %pull the ROIS of interest
bc_table = array2table(betweenness_ROIs); %make it a table
bc_table.Properties.VariableNames = labels;
bc_table.beacon_id = sublist; 

%add labels and sub list
histogram(betweenness_centrality) % chceking distribution
xlabel('Betweenness Centrality')
ylabel('Count')

cd(HomeDir);
writetable(bc_table, 'BetweennessCentrality_n40_37ROIs.csv');

