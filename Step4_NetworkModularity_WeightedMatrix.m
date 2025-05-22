%% Calculate network modularity using weighted signed networks across a range of resolution parameters

% This script is for for Weighted Adjacency Matrices using a
% uniform null model that is run over a series of resolution parameters.

% Script is broadly based on methods from Lacichinetti & Fortuanto 2012 and Cohen & D'Espocito 2016. 

%% Step 1- Set directories
clear all; clc;

HomeDir = '/home/dir';
FCDir = [HomeDir,'/FC_matrix_bivariate'];
BinDir = [FCDir, '/WeightedMatrix/'];
OutputDir = [BinDir,'BCT_output_Step3_NetworkModularity_Weighted/'];
% Create folder to store node partitions for each subject*condition  matrix
partitionDir = OutputDir;

cd(HomeDir);
sub_info = readtable('list_of_subject_ids.csv');
subsToExclude = sub_info.to_exclude;

cd(FCDir);
sublist = load('FCmatrix_151ROIs_N40.mat', 'sublist');
sublist = sublist.sublist;
numSubs = length(sublist);
cd(BinDir);

fprintf('\n Step 1 complete- necessary directories set and subject list created. \n NOTE: Adjust free parameters in Step 3 before running. Also make sure to specify which subjects to exclude \n');

%% Step 2- Find Gamma range based on 50th percentile of all positive functional associations across participants
% (Kennett, Betzel & Beaty, 2020). anything greater than or equal to the
% median connection weight

% Select diretory with input data (symmetric matrices)
sym_matrix_dir = BinDir;
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

% initialize empty matrix to contain all of the processed matrices for each
% subject
allProcessedMatrices = zeros(size(adj_matrix,1),size(adj_matrix,1),length(files));


for i = 1:length(files)
    fprintf(1, 'Compiling all subject matrices: now reading %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load each file, in particular the weighted matrix
    allProcessedMatrices(:,:,i) = A_sym_norm; % put the matrix into the respective subject slice
end
clear i;

% Check that median weights are not calculated using any subjects you want
% to exclude. Use this section below for manual exclusion. 
sub_info.to_exclude(INDEX) = 1; % SALSA057 doesnt have a scan so we need to exclude them

subsToInclude = sub_info.sub_id(sub_info.to_exclude ~= 1) % process all subjects who are being kept

allProcessedMatrices = allProcessedMatrices(:,:,length(subsToInclude)); % conditional indexing based on the subject IDs to keep

B=permute(allProcessedMatrices,[1:length(subsToInclude)]); % go in order of first subject to last of those IDs we keep
out=B(:); % combine all edge weights into one massive vector
%histogram(out);
% calculate STDEV from entire distribution to help with setting gamma range
out_std = std(out);
out_mean = mean(out);
% Callculate median edge weight of positives only - 2STD for
% distribution of gammas
out_pos = out(out>0);
medianWeight = median(out_pos);% median edge weight, corresponding to 50th percentile
histogram(out_pos);
xline(medianWeight, '--r', 'Median Connection Weight');
maxWeight = max(out_pos); % max edge weight, 100th percentile
gamma_range_max = out_mean + (2*out_std)

clear allProcessedMatrices;
fprintf('\n Step 2 complete- adjust number of partitions and gamma range below if desired');

%% Step 3- Load in files, select total number of partitions

%  -------FREE PARAMETERS- ADJUST ACCORDINGLY -------
% Set number of Louvain partitions to create - currently this is based on Cohen & D'Espocito
numPartitions = 150; % select  number of desired partitions
gamma = 0:0.05:gamma_range_max; % 5 possible values of gamma, ranging 0, to the median, and up to max edge weight in 0.05 steps
% b/c we are using the uniform null model, our values of gamma cannot
% exceed the max correlation value (gamma == expected correlation
% coefficient)
%  -------------------------------------------------

% Pre-allocate empty cells with proper dimensions to hold each subjects
% partitions and Q values across costs while varying the level of gamma
levelsOfGamma = length(gamma);

% For a given matrix, generate an empty matrix P of possible Louvain node partitions 
% P will have dimensions = (node# x partition# xgammaVal)
P = zeros(numNodes, numPartitions, levelsOfGamma);

% For a given matrix, generate an empty matrix Qstat of possible Q values across partitions
% Q stat will have dimensions = #partitions
Qstat = zeros(numPartitions,levelsOfGamma);
clear adj_matrix;  % we can clear adj_matrix since we load each file later in the loop in step 3

% Pre-allocate empty matrix with #ROIx Levels of Gamma dimensions 
% Pre-allocate cells with proper dimensions to hold the agreement & consensus matrices for each subject for each level of gamma
% Pre-allocate matrix for final consensus clusters 
% Pre-allocate matrix for final modularity values
D = zeros(numNodes, numNodes, levelsOfGamma);
Dp = zeros(numNodes, numNodes, levelsOfGamma);
Dfinal = zeros(numNodes, levelsOfGamma);
Qfinal = zeros(2, levelsOfGamma);
Qfinal(1,:) = gamma; 
% Pre-allocate empty matrix to hold final Q values for entire group
GroupQ = zeros(numSubs,levelsOfGamma+1);
GroupQ(:,1)=sublist;
fprintf('\n Step 3 Complete- Files loaded, parameters set, cells and matrices pre-allocated for step 3 \n');
 
%% Step 4- Loop through the files and calculate Modularity and generate module structure for each partition at each level of gamma

% Cycle through each file
for index = 1:length(files)
  % Assign file name to variable  
  fileName = files(index).name;
  % Attach path to file name
  fullFileName = fullfile(sym_matrix_dir, fileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Load file into matrix A
  load(fullFileName,'A_sym_norm');
  A = A_sym_norm;
  subject = num2str(sublist(index)); 
  
  % General info - primarily for binary networks
  % We need to iterate through different resolution parameters. When
  % gamma (resolution parameter) < 1, larger communities are resolved,
  % whereas gamma > 1 yeilds more communities containing fewer nodes. When γ < 1, 
  % modularity maximization detects larger (and fewer) communities, resulting in 
  % a higher value for Q; when γ > 1 the algorithm produces several smaller communities
  % and a smaller Q statistic. One method for selecting the resolution parameter is to report 
  % community structure at the value of gamma at which partitions are most similar to
  % each other (Bassett et al. 2013 & Sporns & Betzel 2016). This can be
  % computed using normalized mutual information (use the z-score). As a
  % benchmarch, Cohen & D'Esposito 2016 used gamma = 1.25. According to
  % Fortunato & Barthelemy 2007 " a value of Q larger than 0.3–0.4 is a clear 
  % indication that the subgraphs of the corresponding partition are
  % modules." From Rubinov & Sporns 2011: "The assumption that positive and 
  % negative weights are equally important is neurobiologically problematic,
  % because the role and importance of positive and negative weights in functional
  % networks is intrinsically unequal. Positive weights associate nodes with modules
  % explicitly, while negative weights associate nodes with modules implicitly, by 
  % dissociating nodes from other modules." 
  
  % ----------------------- 
  % CHOICE AND REASONING FOR UNIFORM NULL MODEL IN COMMUNITY_LOUVAIN

  % WHY IS THE DEFAULT NULL IN BCT SUBOPTIMAL?
  % Because correlational networks (correlation coefficient = edges) are
  % inherently non-independent and more clustered than random networks,
  % they are more likely to be small-world networks. This presents an issue
  % when using a null model that only preseveres degree distribution, so
  % they need to be benchmarked against null models that respect the
  % topological structure induced by the correlational nature or transivity 
  % (This is due to the transitive nature of the correlation coefficient, which effectively 
  % means indirect paths and their corresponding direct connections are more likely to 
  % coexist in correlation networks compared to random networks. 
  % Otherwise, the extent of small-world structure will be over estimated with full correlation or under
  % estimated with partial correlation. (Zalesky, Fornito, Bullmore, 2012).
  
  % WHY IS U NULL MODEL BETTER? ﻿More recently, another study pointed out that the configuration
  % model has a complicated interpretation in the context of correlation networks (appropriateness aside) 
  % ( Bazzi et al., 2016; MacMahon, Gar- laschelli ). As an alternative, the authors suggested that a uniform 
  % null model may be appropriate for correlation matrices. The uniform null model is one in which it is assumed 
  % that every element in the network is mutually correlated with the same magnitude. That is, P=1*gamma, where X is a 
  % matrix whose elements are all 1 and gamma is the magnitude of mutual correlation. This model performed well 
  % in benchmark analyses and has the added advantage of not suffering from the resolution limit ( Traag et al., 2011 ).
  % Notably, this model also has a convenient interpretation. Communities discovered at a given gamma value correspond to 
  % groups of nodes whose mean connectivity to one another exceeds a value of gamma on average. This truism facilitates 
  % a clearer interpretation of the relationship between communities and
  % the resolution parameter (Esfahlani...Betzel, 2014, NeuroImage)
  % -----------------------
  
 

  for i = 1:levelsOfGamma % for each level of gamma
     for iter = 1:numPartitions % for each individual partition
        [M,Q1] = community_louvain(ones(numNodes),[],[],((A - gamma(i)).*~eye(length(A)))); % calc modularity and partition using uniform null
        P(:,iter,i) = M; % save partition into P for partition # at a given level of gamma
        Qstat(iter,i) = Q1; % save Q into Qstat for each partition at a given level of gamma
     end 
  end
  
  clear i iter ;
  
   %[M,Q1] = community_louvain(A(:,:),gamma(i),[],'negative_asym'); calculate partition and modularity at a given level of gamma using neg_asynm 
  
  
  
  % Now need to generate an agreement matrix at a given value of gamma -
  % this tells us the number of times that a node was assigned to the same
  % module across the 150 partitions. 
 
  for i = 1:levelsOfGamma
      D(:,:,i) = agreement(P(:,:,i)); % put agreement matrix for given gamma in D
  end 
  
  clear i; 
  
  % Now, need to create a consensus matrix across each layer of 150 partitions.
  % Then, determine the optimal cost and level of gamma based on the
  % conensus matrix at each cost. Cohen & D'Espocito thresholded each cell
  % in the consensus matrix to 0.5, so that the agreement of all pairs that
  % were not assigned to the same matrix at least 50% of the time were set
  % to 0. Then, they created a consensus partition by running the Louvain
  % algorithm on the agree ment matrix 100 times until a single
  % representative partition was obtained
  
  % Convert agreement matrix D to proportion of agreement matrix Dp .
  % This line converts counts in the agreement matrix to probabilities. 
  % Each element in Dp will have a value equal to or less than the total 
  % number of partitions. We do this by dividing matrix D by the scalar
  % total number of partitions. This converts each element into a proportion.
  for i = 1:levelsOfGamma
      Dp(:, :, i) = D(:,:,i)/numPartitions;  
  end
   
  clear i; 
  
  %Use BCT consensus function to run the agreement matrix through the
  %louvain algorithm again. Set threshold to .3 Lancichinetti & Fortuando 
  %2012 recommend a low threshold --.3 or below
  
  % Current tau/threhsold = .3, 
  
 for i = 1:levelsOfGamma % for each level of gamma
     Dfinal(:,i) = consensus_und(Dp(:,:,i),.3,numPartitions); % determine consensus clustering of each slice, threshold at .3, and run 150 times
 end 
 

 % Get the final Q values from each consensus partition
 % Qfinal = zeros(numSubs, levelsOfGamma);
 
 for i = 1:levelsOfGamma
     % ((A - gamma(i)).*~eye(length(A)))
     [finalM, FinalQ]=community_louvain(A,gamma(i),Dfinal(:,i),'negative_asym');
     Qfinal(2,i) = FinalQ; % first row is the level of gamma, second if final Q from consensus partition
 end
 
 GroupQ(index,2:end)= Qfinal(2,:); 
 
  % need to save for each subject 
  save([OutputDir,subject,'_FinalLouvainPartitions.mat'],'Dfinal','Qfinal','P','D','Dp','Qstat','sublist','numNodes','numPartitions', 'gamma')
  save([OutputDir, 'GroupFinalQValues.mat'], 'GroupQ', 'gamma');
  % Dfinal is the optimal clustering based on consensus clustering
  % Qfinal holds the final Q values based on consensus partitions
  % P is each partition (150)
  % D is agreement matrix, 
  % Dp is percentage of of modules assigned to same module across particition
  % QStat holds the Q calues for each of the 150 partitions across each
  % level of gamma
  % sublist is the suject list
  % numNodes is tGroupQ(subjectsToExclude,:) = NaN; he number of nodes (ROIs)
  % numPartitions tells us the number of partitions were ran and how many
  % times wer ran the consensus clustering 
  % gamma tells us our range of gamma values

  subject = num2str(sublist(index));
  fprintf(1, 'Step 4 Complete- Final Partitions and Modularity Values saved for subject %s\n', subject);

 
end

fprintf(' \n Step 4 completed \n');
%% Save final values 
FinalGroupQ = GroupQ % Copy GroupQ to hold final Q values for all subjects across all levels of gamma
GroupQ(:,9) = subsToExclude; % add a column where a '1' specifies that they should be excluded, and '0' means to keep
FinalGroupQ((GroupQ(:,9) == 1),2:end) = NaN; %wherever there is a 0, set the Q values across gamma levels to NaN

% convert total group Q values to csv
cd(OutputDir);
csvwrite('GroupFinalQValues_n40.csv', FinalGroupQ);

