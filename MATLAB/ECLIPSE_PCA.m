% Function to build SCA model on the remaining cells after removal of normal cells

% Dependecies: 
% blockscale_struct
% PCAecon

%% Written by Rita Folcarelli
%%
function [model_eclipse, S_scores_eclipse,S_scores_PC_used_eclipse,S_error_eclipse]=ECLIPSE_PCA(S_scaled_eclipse, pc_used_eclipse)


Labels=vertcat(S_scaled_eclipse.Labels);
[patient_prep_blockscaled] = blockscale_struct(S_scaled_eclipse(vertcat(S_scaled_eclipse.Labels) == 1));

[~, LDS, VAR] = PCAecon(patient_prep_blockscaled);

% projection of data in the ECLIPSE space

S_scores_eclipse = S_scaled_eclipse;
for l1 = 1:length(S_scaled_eclipse)
    S_scores_eclipse(l1).Data = S_scaled_eclipse(l1).Data*LDS;
end

S_scores_PC_used_eclipse=S_scores_eclipse;
for l1 = 1:length(S_scaled_eclipse)
    S_scores_PC_used_eclipse(l1).Data = S_scaled_eclipse(l1).Data*LDS(:,pc_used_eclipse);
end 

% estimation of residuals in the control space

S_error_eclipse=S_scores_eclipse;
for l1 = 1:length(S_scaled_eclipse)
    S_error_eclipse(l1).Data = S_scaled_eclipse(l1).Data-S_scores_PC_used_eclipse(l1).Data*LDS(:,pc_used_eclipse)';
end


control_scores=vertcat(S_scores_eclipse(Labels==0));
patient_scores=vertcat(S_scores_eclipse(Labels==1));
Tc_eclipse = vertcat(control_scores.Data);
Tp_eclipse= vertcat(patient_scores.Data);


control_error_eclipse=vertcat(S_error_eclipse(Labels==0));
patient_error_eclipse=vertcat(S_error_eclipse(Labels==1));
Ec_eclipse=vertcat(control_error_eclipse.Data);
Epc_eclipse=vertcat(patient_error_eclipse.Data);



model_eclipse=struct('Tc',Tc_eclipse , 'VAR', VAR, 'loadings', LDS, 'Tp',Tp_eclipse, 'Ec',Ec_eclipse,'Ep', Epc_eclipse, 'pc_used', pc_used_eclipse);

end

