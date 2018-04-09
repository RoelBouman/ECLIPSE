% Funtion to build SCA model on individual sample and visualize the results
% Drawback is that the personal models cannot be compared between individuals

% INPUT:  S_in contained the cells after removal of normal cells
% OUTPUT: plot


%% Written by R.Folcarelli (based on script written by R. Bouman)
%%

function ECLIPSE_model_perindividual(S_in, DensityECLIPSEmodel, paired_data, pre_process_modes, pc_used, VariableNames)


[S_scaled_ECLIPSE, ~] = Pre_process_MFC_takingintoaccountthenumberofcells(S_in, paired_data, pre_process_modes);

indQ = input('Do You want to visualize some, or all responding individuals using individual ECLIPSE model? (some/all) : ', 's');

if(strcmp(indQ,'some')||strcmp(indQ,'Some'))
    
    allPIDs = vertcat(S_in.ID);
    for i = 1:length(allPIDs)
        disp(strcat(int2str(i),': ID: ',int2str(allPIDs(i))));
    end
    IDs = input('Enter the IDs of the individuals you want to visualize. ([1 3 5], 1:5, etc.) : ');
    
elseif(strcmp(indQ,'all')||strcmp(indQ,'All'))
    
    IDs = vertcat(S_scaled_ECLIPSE.ID);
    
end

for l1 = 1:length(IDs)
    
    [~, LDS, VAR] = PCAecon(S_scaled_ECLIPSE(vertcat(S_scaled_ECLIPSE.ID)==IDs(l1)).Data);
    
    REC =S_scaled_ECLIPSE(vertcat(S_scaled_ECLIPSE.ID)==IDs(l1)).Data*LDS;
    
    [~,DensityResponse,X,Y,~,~,~] = kde2dimproved(REC(:,pc_used), size(DensityECLIPSEmodel,1));
    
    [DensityResponse] = ProcessDensity(DensityResponse);
    
    figure
    PlotDens(DensityResponse,X,Y,'imagesc', pc_used,VAR,LDS,VariableNames);
    title(strcat('Density of individual: ', int2str(IDs(l1)), ' using ECLIPSE per individual'))
    
end
end