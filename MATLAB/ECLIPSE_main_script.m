%setting for the figures
set(0,'Defaulttextfontsize',20); set(0,'Defaultaxesfontsize',20);set(0,'Defaultlinemarkersize',12);set(0,'Defaultlinelinewidth',2);set(0,'Defaultaxesfontweight','bold')
load blue_map2; load red_map2; load red_map; load blue_map_intense; load mycmap; load white_jet

%% SECTION A import/load data
% to import data from fsc, lmd, txt, csv, format. 
% Choose import if you need to import data, if you have already the .mat struct, load it and go to the next section

[S_original, VariableNames_original] = import_or_load_data([]);

% %% SECTION B (Optional)Outlier removal based on 95% quantiles of each variable
% % Remove the outlier cell caused by instrumental problems
% % Option for remove outlier only within controls, patients or both groups
% 
% [~,S_original,percentage_retained_idv,percentage_retained_total] =rem_outliers_class(S_original,5, 'control');
% % class='control'
% % class='patient'
% % class='all'
% % output:   S_original struct of data without extreme outliers within the
% %           class chosen
%% SECTION B data selection
% select the variables you wish to use for the analysi
% paired data question
paired_data = paired_data_question;
% paired data=1 if the two groups consists of the same individual (e.g. before and after treament)
% paired data=0 if the two groups NOT consists of the same individual
selection_mode = [];
[S, VariableNames, selection_mode] = data_selection(S_original, VariableNames_original,selection_mode);
%selection_mode{1}=[];          % variables used
%selection_mode{2}=0;           % labels control
%selection_mode{3}=1;           % labels patient(and/or treated)
%% SECTION C preprocessing (Step 0a flowchart)

% pre-process: different options are possible (see the Supplementary Material I of the manuscript)

 [S_scaled, pre_process_modes] = Pre_process_MFC_takingintoaccountthenumberofcells(S, paired_data, []);

% example, pre_process_modes:[2 1 2]:  arcsinh transform; mean center over all controls; scaling over
%           all the controls
% []:       leave empty if you wish to change the chosen pre-processing
%% Sub-SECTION C
% visualize the effect of the pre-processing
hist_markers(S_scaled, VariableNames);
% per class: red histograms for responder/diseased, blue histograms for
% controls
%% SECTION D - Simultaneous Component Analysis (Step 1 and Step 3 flowchart)
[Pr] = PRESSplot_control(S_scaled);
pc_used=input('Which PCs should be used for the control model?'); % used the first two PCs; change it if you want to check other PCs combination
[controlmodel, S_scores, S_error, S_scores_PC_used, ~]=ControlPCA(S_scaled, pc_used);

%% SECTION E - Compute densities (Step 2A and 2B flowchart)
% density estimation are normalized over the number of cells per each set

[DensityControl, DensityResponse, X, Y, bandwidth1, bandwidth2, MIN_XY, MAX_XY ] = ComputeWeightedDensities2D(S_scores_PC_used);
[DensityControl, DensityResponse, DensityDiff, ~, ~] = ProcessDensities(DensityControl, DensityResponse);

%% SECTION Fa- Visualization of the estimated probability densities per class
figure; PlotDens(DensityControl,X,Y,'imagesc', pc_used,controlmodel.VAR,controlmodel.loadings,VariableNames);colormap(blue_map_intense); title('Control Density')
figure; PlotDens(DensityResponse,X,Y,'imagesc', pc_used,controlmodel.VAR,controlmodel.loadings,VariableNames);colormap(red_map); title('Response Density')
figure; PlotDiff(DensityDiff, X, Y,'imagesc', pc_used,controlmodel.VAR,controlmodel.loadings,VariableNames);title('Difference between Control Densiy (blue) and Response Density (red)')
tileFigs

%% SECTION Fb - Visualization of Density of Patients and Difference between Densities

IDs=question_ID(S_scaled);
PlotDensInd (S_scores_PC_used,S_scaled,MIN_XY,MAX_XY,IDs,pc_used, controlmodel, VariableNames)
PlotDensDiffIndividuals (S_scores_PC_used,S_scaled,MIN_XY, MAX_XY, DensityControl, IDs, pc_used,controlmodel,VariableNames);tileFigs

%% SECTION G - Indices to extract 'abnormal' cells 
% build index to select cells that exceed the SumPredictionError limit

% build indices to select cells lying in the Density region specific to the
% response
[inoroutDbD, ~, ~]=DbDlimit_linear(S_scores_PC_used);

[inoroutSPE, ~, count_SPE, SPElim95]=SPElimitfunc(S_scaled, S_scores, controlmodel , pc_used);
%% SECTION H: Elimination Step
% Define abnormal cells cells that exceed the SPE limit
[S_in_original, percentage_eclipse]=Identify_response_cells(S_original, inoroutSPE, inoroutDbD);
open percentage_eclipse
%% SECTION J - Select again the variables for S_in (after removal of normal cells) 
[S_in, ~, ~]=data_selection(S_in_original, VariableNames_original, selection_mode); % selection mode defined in SECTION C

%% SECTION K - SCA on the patient/response group
% option:
% general model: to build the SCA model on the all the individuals
% partial model: to build the SCA model only on a sbusets of the data


[modelECLIPSE, DensityECLIPSEmodel, S_scaled_ECLIPSE, S_scores_PC_used_eclipse,X,Y,MIN_XY,MAX_XY, pre_process_modes, pc_used, bins]=ECLIPSE_space_general_partial(S_in, paired_data, DensityControl, VariableNames);
 
%% Plotting the density of the model obtained in section K
% options:
% visualize all the individuals separately 
% visualize only few individuals
plotting_ECLIPSE(S_scores_PC_used_eclipse,S_scaled_ECLIPSE,MIN_XY,MAX_XY, modelECLIPSE, VariableNames)

%% SECTION L - Build individual ECLIPSE model
% Note: Drawback is that the personal models cannot be compared between individuals
ECLIPSE_model_perindividual(S_in,DensityECLIPSEmodel, paired_data, pre_process_modes, pc_used, VariableNames)
