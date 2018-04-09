function PlotDensDiffIndividuals (S_scores,S_scaled,MIN_XY, MAX_XY, DensityControl, IDs, pc_used,controlmodel,VariableNames)
%% PlotDensDiffIndividuals plots one or more individuals against the density of a control group
%
%% Dependencies:
%  kde2dimproved.m
%  ProcessDensities.m
%  PlotDiff.m
%
%% Input:
%  Data_Struct, a Structure containing processed data in 2 dimensions in
%  the .Data entry, a numerical identifier in the .ID entry, and a label in
%  the .Labels entry.
%
%  MIN_XY and MAX_XY, indicating the minimum and maximum in the entirety of
%  the DensityControl
%
%  DensityControl, a grid containing the Density evaluated at every point
%  in the space. DensityControl must be a n-by-n grid where n is a power of
%  2.
%
%  IDs, a vector containing the IDs, which are matched with the .ID
%  category of the Data_struct, which signify the ID of all individuals to
%  be compared.
%
%  pc_used, a 2-by-1 vector containing the PCs which were used. They are
%  used for setting the labels on the X and Y axis for entries 1 and 2
%  respectively.
%
%  controlmodel, a structure containing the processed data of the control
%  group. .Pc contains the Scores of the control group, and .VAR contains
%  the variances accounted for by each pc.
%
%  VariableNames (optional), a vector containing the names of the used
%  variables. When it is used, a loadingsplot is plotted over the density
%  difference.
%
%% Output:
%  None
%
%% Author: Roel Bouman (roelbouman@gmail.com). Modified by R. Folcarelli to allow the plot of paired data.
%% Remarks:
%  The function ends with hold off after plotting loadings. Consider this
%  when using it.
Labels=vertcat(S_scores.Labels);
ID=vertcat(S_scores.ID);
ID2=vertcat(S_scaled.ID);
N_PCs = length(pc_used);

if isfield(controlmodel(1), 'Pp')
    LDS=controlmodel.Pp;
elseif isfield(controlmodel(1), 'Pc')
    LDS=controlmodel.Pc;
elseif    isfield(controlmodel(1), 'loadings')
    LDS=controlmodel.loadings;
end
IDs=unique(IDs);
two_plots = false;
for l1 = 1:length(IDs)
    S_scoresID=vertcat(S_scores(ID==IDs(l1)));
    S_scaledID=vertcat(S_scaled(ID2==IDs(l1)));
    if length(S_scoresID)==2
        %     if vertcat(S_scores_ECLIPSE.ID)==IDs(l1)
        S_scoresID=vertcat(S_scores(ID==IDs(l1)& Labels==0));
        S_scaledID=vertcat(S_scaled(ID2==IDs(l1)& Labels==0));
        two_plots = true;
    end
    [~,DensityInd,X,Y,~,~,~] = kde2dimproved(S_scoresID.Data,1024,MIN_XY,MAX_XY);
    [~,~,DensityDiff,~,~,~] = ProcessDensities(DensityControl, DensityInd);
    
    figure
    
    if(nargin == 9)
        PlotDiff(DensityDiff,X,Y,'imagesc',pc_used,controlmodel.VAR,controlmodel.loadings,VariableNames);
    else
        PlotDiff(DensityDiff,X,Y,'imagesc',pc_used,controlmodel.VAR);
    end
    
    
    title(strcat('Individual with ID: ',int2str(IDs(l1)), ' compared to the control model'),'FontSize',20);
    VAR = zeros(1, N_PCs);
    for ii = 1:N_PCs
        VAR(1,ii) = 1 - sum(var(S_scaledID.Data - S_scoresID.Data(:,ii)*LDS(:,pc_used(ii))'))/sum(var(S_scaledID.Data));
    end
    VAR = round((VAR*100));
    xlabel(strcat('PC ', int2str(pc_used(1)),': Individual VAR: ', num2str(VAR(pc_used(1))),'%'));
    ylabel(strcat('PC ', int2str(pc_used(2)),': individual VAR: ', num2str(VAR(pc_used(2))),'%'));
    
    if (nargin == 9)
        hold on
        DA_FLplotloadings(LDS(:,pc_used),VariableNames);
        hold off
    end
    
    
    %
    if two_plots
        S_scoresID=vertcat(S_scores(ID==IDs(l1)& Labels==1));
        S_scaledID=vertcat(S_scaled(ID2==IDs(l1)& Labels==1));
        
        [~,DensityInd,X,Y,~,~,~] = kde2dimproved(S_scoresID.Data, 1024,MIN_XY,MAX_XY);
        [~,~,DensityDiff,~,~,~] = ProcessDensities(DensityControl, DensityInd);
        
        Xisc = X(1,:);
        Yisc = Y(:,1);
        
        figure
        
        if(nargin == 9)
            PlotDiff(DensityDiff,X,Y,'imagesc',pc_used,controlmodel.VAR,controlmodel.loadings,VariableNames);
        else
            PlotDiff(DensityDiff,X,Y,'imagesc',pc_used,controlmodel.VAR);
        end
        
        
        title(strcat('Individual with ID: ',int2str(IDs(l1)), ' compared to the control model'),'FontSize',20);
        VAR = zeros(1, N_PCs);
        for ii = 1:N_PCs
            VAR(1,ii) = 1 - sum(var(S_scaledID.Data - S_scoresID.Data(:,ii)*LDS(:,pc_used(ii))'))/sum(var(S_scaledID.Data));
        end
        VAR = round((VAR*100));
        xlabel(strcat('PC ', int2str(pc_used(1)),': Individual VAR: ', num2str(VAR(pc_used(1))),'%'));
        ylabel(strcat('PC ', int2str(pc_used(2)),': individual VAR: ', num2str(VAR(pc_used(2))),'%'));
        
        if (nargin == 9)
            hold on
            DA_FLplotloadings(LDS(:,pc_used),VariableNames);
            hold off
        end
        
        
    end
    
    
    
end