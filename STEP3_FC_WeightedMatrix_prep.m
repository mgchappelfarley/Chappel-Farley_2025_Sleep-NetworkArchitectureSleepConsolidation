%% Format Weighted Matrices for BCT Graph Theory Analyses

%% Step 1. Set directories and load subject list mat file.
clear all; clc;
HomeDir =  '/home/dir/';
FCDir = [HomeDir, '/FC_matrix_bivariate/'];
cd(FCDir);
FCmatrix = load('FCmatrix_151ROIs_N40.mat');

sublist = FCmatrix.sublist;

fprintf('Step 1 Complete- directories set and subject list loaded \n');

%% Step 2. Initialize zeros matrix to hold information about negative edge weights, symmetrize matrix, save.
negatives_summary = zeros(length(sublist),2); % in case later we want to know how many negative edges we had

for i = 1:length(sublist)
      subject = num2str(sublist(i));
      disp(['Running ', subject]);
    
    %% Symmetrize/Fix Matrix
      % averaging corresponding lower and upper diagonal elements in the matrices (semipartial values are not symmetric)
     
      % Load matrix output from CONN
      load([FCDir,subject,'_FCmatrix_151ROIs.mat'])
      A = FCmatrix; % make A into our functional connectivity matrix
     
      %Make matrix symmetric by summing it with its transpose and dividing all elements in half. This is essentially taking the average of the
      %corresponding elements across the diagonal, and setting each element to the average. The result is a symmetric matrix, where each element in the
      %average of the two beta-interaction terms for eachs et of  ROI-pairs
      A_sym = (A + A') / 2;
      
      %This line requires BCT toolbox to be installed. It cleans up symmetric matrices in numerous ways, including setting the diagonal to zeros 
      A_sym = weight_conversion(A_sym, 'autofix'); % sets NaNs to zero (& along diagnonal), ensures symmetry, etc
     
      %Normalization rescales all weight magnitudes to the range [0,1].Optimal for weighted analyses
      A_sym_norm = weight_conversion(A_sym, 'normalize'); 

      % count how many negative values, in case we want to know later. 
      negatives_orig = (length(nonzeros(A_sym(A_sym<0)))/2);
      negatives_orig_perc = negatives_orig/((length(A)*length(A))/2);
      negatives_summary(i,1) = negatives_orig; %saves as a group summary stat
      negatives_summary(i,2) = negatives_orig_perc;
      
      % IF YOU WANT ONLY POSITIVE MATRIX, LINE BELOW WILL REMOVE NEGATIVE VALUES (for no-cost weighted matrices)
      % A_sym_norm_pos = max(A_sym_norm,0); % to just look at positives but still weighted

      %Save Weighted Matrix, Consisting of both positive and negative
      %values.
      
      save([FCDir,subject,'_FCmatrix_BCTweighted.mat'],'A_sym_norm', 'A_sym','names','negatives_orig_perc')
      
      %clear A* negatives_orig* FCmatrix
     
      
    %% STEP 3: Save relevant variables
    
     save([FCDir,subject,'_allROIs_FCmatrix_processed.mat'],'A_sym','A_sym_norm','negatives_orig','negatives_orig_perc',...
         'names')
    

     clear *sym* allROIs* conn* A cost 
     
     
    disp(['Completed ', subject]);
end % i

fprintf('Complete- Symmetrical Weighted Matrix Saved \n');

%% FOR GROUP SUMMARY STATS
    data = [sublist, negatives_summary];
    varnames = {'subject','num_neg','neg_pnct'};
    save([FCDir,'negative_values_bivariate_summary_N40_151ROIs.mat'],'data','varnames')

    fprintf(' Symmetrical Weighted Matrix Saved \n');

%     varnames = {'subjects','negatives_orig','negatives_orig_perc','neg_cost05','neg_cost10','neg_cost15','neg_cost20',...
%         'neg_cost25','neg_cost30','neg_cost35','neg_cost40','neg_cost45','neg_cost50'};