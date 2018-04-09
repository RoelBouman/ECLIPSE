% Create PCA model on control data with number of PCs given

%controlmodel, function that builds controlmodel space on control
%individual
% INPUT:                        S_scaled (structure)
% OUTPUT: controlmodel(struct)  Tc (scores of controls)
%                               Pc (loadings of controls)
%                               Tpc (scores of patients)
%                               VAR (variance explained by PCs)

%% Dependecies:
% PCAecon
% blockscale_struct
% NDhist
%%

function [controlmodel, S_scores, S_error, S_scores_PC_used, Hist_control]=ControlPCA(S_scaled, pc_used)
%% initialize starting values:
Labels=vertcat(S_scaled.Labels);
number_bins = 250000;
N_PCs = length(pc_used); 
minmaxcutoff = 0.001; %part not taken for determining the edges, standard: 0.01  = 1%
edgesfactor = 1; %Increase the edgesfactor.
smooth_factor = 5; %smooth_factor  = lambda (article), 1-5 normal, 0 for no smoothing.
if isfield(S_scaled(1), 'train')
    trainset = logical(vertcat(S_scaled.train));
else
    trainset = true(length(S_scaled),1);
end

%% Create PCA model based on controls which belong to the training set
[control_prep_blockscaled] = blockscale_struct(S_scaled(Labels == 0 & trainset));
[~, LDS, VAR] = PCAecon(control_prep_blockscaled);
for L1=1:(min(length(VAR),7))
    if sum(L1 == pc_used) == 1
        disp([ '  Used | '  ' | Expained Variance - PC ' num2str(L1) ' = ' num2str(VAR(L1)) ' %']);
    else
        disp([ 'Unused | '  ' | Expained Variance - PC ' num2str(L1) ' = ' num2str(VAR(L1)) ' %']);
    end
end

%% projection of data in the control space

S_scores = S_scaled;
for l1 = 1:length(S_scaled)
    S_scores(l1).Data = S_scaled(l1).Data*LDS;
end

S_scores_PC_used=S_scores;
for l1 = 1:length(S_scaled)
    S_scores_PC_used(l1).Data = S_scaled(l1).Data*LDS(:,pc_used);
end 

S_error=S_scores;
for l1 = 1:length(S_scaled)
    S_error(l1).Data = S_scaled(l1).Data-S_scores(l1).Data(:,pc_used)*LDS(:,pc_used)';
end



%% Create binedges
X = vertcat(S_scores(trainset).Data);
if N_PCs == 1
    binsize = 1000;
else
    binsize = repmat(round(number_bins^(1/N_PCs)),1,N_PCs);
end

edges = binedges(X(:,pc_used), binsize, minmaxcutoff, edgesfactor);

%% Create Histograms
% Time: 90 seconds

Hist_control = zeros([length(S_scores) binsize]);
for l1 = 1:length(S_scores)
    Hist_control(l1,:,:,:,:,:,:) = NDhist(S_scores(l1).Data(:,pc_used), edges, smooth_factor);
end

%% Calculate benchmark
CI=0.8;   
ControlCH_control = con_hul2(permute(Hist_control(Labels == 0 & trainset,:,:), [2 3 1]), CI*100);
ControlCH_patient = con_hul2(permute(Hist_control(Labels == 1 & trainset,:,:), [2 3 1]), CI*100);
%% Create output struct
controlmodel.loadings = LDS;
controlmodel.VAR = VAR;
controlmodel.control_error=vertcat(S_error(Labels==0));
controlmodel.patient_error=vertcat(S_error(Labels==1));
controlmodel.Ec=vertcat(controlmodel.control_error.Data);
controlmodel.Epc=vertcat(controlmodel.patient_error.Data);
controlmodel.edges = edges; 
controlmodel.CH_control = ControlCH_control;
controlmodel.CH_patient = ControlCH_patient; 
controlmodel.Labels = Labels; 
controlmodel.ID = vertcat(S_scaled.ID); 
controlmodel.CI = CI; 
end

