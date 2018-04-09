% Funtion to build SCA model on individual sample and visualize the results
% Drawback is that the personal models cannot be compared between individuals

% INPUT:  S_in contained the cells after removal of normal cells
% OUTPUT: plot

%% Dependencies:
% Pre_process_MFC_takingintoaccountthenumberofcells_ECLIPSE
% kde2dimproved
% ProcessDensity
% PlotDens
%
%% Written by R.Folcarelli (based on a script written by R. Bouman)


function [modelECLIPSE, DensityECLIPSEmodel, S_scaled_ECLIPSE, S_scores_PC_used_eclipse,X,Y,MIN_XY,MAX_XY,pre_process_modes, pc_used_eclipse, bins]=ECLIPSE_space_general_partial(S_in, paired_data, DensityControl, VariableNames)
load('mycmap')

invalidinput = 1;

    
    indQ = input('Do you want to build a partial response space or the general response? (partial/general): ', 's');
    
    pc_used_eclipse=input('Which PCs do you want to use?');
    
    if(strcmp(indQ,'general')||strcmp(indQ,'General'))
        
        [S_scaled_ECLIPSE, pre_process_modes] = Pre_process_MFC_takingintoaccountthenumberofcells_ECLIPSE(S_in, paired_data, [1 5 1]);
        
        
        
        [modelECLIPSE, ~,S_scores_PC_used_eclipse,~]=ECLIPSE_PCA(S_scaled_ECLIPSE,pc_used_eclipse);
        
        [~,DensityECLIPSEmodel,X,Y,MIN_XY,MAX_XY,bins] = kde2dimproved(modelECLIPSE.Tp(:,pc_used_eclipse), size(DensityControl,1));
        [DensityECLIPSEmodel] = ProcessDensity(DensityECLIPSEmodel);
        
        figure
        PlotDens(DensityECLIPSEmodel,X,Y,'imagesc', pc_used_eclipse,modelECLIPSE.VAR,modelECLIPSE.loadings,VariableNames);
        title('Density of the Response model')
        colormap(mycmap)
        
        
        
    elseif(strcmp(indQ,'partial')||strcmp(indQ,'Partial'))
        
        disp('The following IDs belong to Responders:');
        allPIDs = vertcat(S_in(vertcat(S_in.Labels)==1).ID);
        
        for i = 1:length(allPIDs)
            disp(strcat(int2str(i),': ID: ',int2str(allPIDs(i))));
        end
        
        RMIDs = input('Enter the IDs of the individuals you want to base the model on. ([1 3 5], 1:5, etc.) : ');
        
        S_temp = S_in(ismember(vertcat(S_in.ID),RMIDs));
        
        [S_scaled_ECLIPSE, pre_process_modes] = Pre_process_MFC_takingintoaccountthenumberofcells_ECLIPSE(S_temp, paired_data, [1 5 1]);
        
        [modelECLIPSE, S_scores_ECLIPSE,S_scores_PC_used_eclipse,~] = ECLIPSE_PCA(S_scaled_ECLIPSE, pc_used_eclipse);
        
        [~,DensityECLIPSEmodel,X,Y,MIN_XY,MAX_XY,bins] = kde2dimproved(modelECLIPSE.Tp(:,pc_used_eclipse), size(DensityControl,1));
        [DensityECLIPSEmodel] = ProcessDensity(DensityECLIPSEmodel);
        
        figure
        PlotDens(DensityECLIPSEmodel,X,Y,'imagesc', pc_used_eclipse,modelECLIPSE.VAR,modelECLIPSE.loadings,VariableNames);
        title('Density of the partial Response model')
        colormap(mycmap)
        
    else 
        
    end

end 

