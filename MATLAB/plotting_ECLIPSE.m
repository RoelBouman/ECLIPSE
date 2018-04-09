
% Funtion to visualized partial or general SCA model built with the
% function ECLIPSE_space_general_partial


% INPUT:  
% - S_scores_PC_used_eclipse: struct containing the scores of the two PCs
% chosen in the model
% - S_scaled_ECLIPSE:   struct containing the data trasformed
% - MIN_XY MIN_XY:      axis limit necessary for the figure
% - modelECLIPSE       
% - VariableNames
% OUTPUT: plot

%% Dependencies:
% PlotDensInd
%% Written by R.Folcarelli (based on a script by R. Bouman)

%%
function plotting_ECLIPSE(S_scores_PC_used_eclipse,S_scaled_ECLIPSE,MIN_XY,MAX_XY, modelECLIPSE, VariableNames)


indQ2 = input('Do You want to visualize some, none, or all responding individuals in the general response space using ECLIPSE? (some/none/all) : ', 's');

if(strcmp(indQ2,'some')||strcmp(indQ2,'Some'))
    
    disp('The following IDs can be visualized:');
    allPIDs = vertcat(S_scores_PC_used_eclipse.ID);
    
    for i = 1:length(allPIDs)
        disp(strcat(int2str(i),': ID: ',int2str(allPIDs(i))));
    end
    
    IDs = input('Enter the IDs of these individuals which you want to visualize. ([1 3 5], 1:5, etc.) : ');
    
elseif(strcmp(indQ2,'all')||strcmp(indQ2,'All'))
    
    IDs = vertcat(S_scores_PC_used_eclipse.ID);
    
elseif(strcmp(indQ2,'none')||strcmp(indQ2,'None'))
    
else
    disp('Input was invalid! Please re-enter');
    
end

PlotDensInd (S_scores_PC_used_eclipse,S_scaled_ECLIPSE,MIN_XY,MAX_XY,IDs,modelECLIPSE.pc_used, modelECLIPSE, VariableNames)

end
