function PlotDensInd (S_scores,S_scaled,MIN_XY,MAX_XY,IDs,PCs, controlmodel, VariableNames)
%% PlotDens plots the density as either an imagesc or surf plot.
% Variation of PlotDens  accomodating for plotting the Ind variance for
% each individual. Replace it instead of PlotDens whenever you want to plot
% Individual density.

%% Dependencies:
%  Da_FLplotloadings.m
%  individual_VAR
%
%% Input:
%
%  Density is a n-by-n matrix containing the local densities evaluated at
%  meshgrid coordinates given by X and Y.
%
%  X must either be a n-by-n matrix with identical elements on each column
%
%  Y must either be a n-by-n matrix with identical elements on each row
%
%  plotMode is a string argument with value 'imagesc' or 'surf' for imagesc
%  and surface plots respectively.
%
%  PCs indicate the labels, or indices of the used Principal Components,
%  this must be a 2-by-1 vector
%
%  VAR is a 1-by-p vector, where p is the total number of PC's, indicating
%  the variances accounted for by each PC.
%
%  LDS are the Loadings of the control group constructed with PCA
%
%  VariableNames are the name of the used variables, given in a cell vector.
%
%% Output:
%  None
%
%% Author: Roel Bouman (roelbouman@gmail.com). Date: 3-05-2016. Modified by R. Folcarelli to allow the plotting of paired data.

%% Remarks:
%  The function ends with hold off after plotting loadings. Consider this
%  when using it.
Labels=vertcat(S_scores.Labels);
ID=vertcat(S_scores.ID);
ID2=vertcat(S_scaled.ID);
N_PCs = length(PCs);

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
    [DensityInd] = ProcessDensity(DensityInd);
    
    Xisc = X(1,:);
    Yisc = Y(:,1);
    
    figure
    imagesc(Xisc(1,:), Yisc(:,1), DensityInd,[0 .5]);
    set(gca,'Ydir','normal')
    colorbar;
    
    title(strcat('Density of individual: ', int2str(IDs(l1)), ' Label ',num2str(S_scoresID.Labels)),'FontSize',20)
    VAR = zeros(1, N_PCs);
    for ii = 1:N_PCs
        VAR(1,ii) = 1 - sum(var(S_scaledID.Data - S_scoresID.Data(:,ii)*LDS(:,PCs(ii))'))/sum(var(S_scaledID.Data));
    end
    VAR = round((VAR*100));
    xlabel(strcat('PC ', int2str(PCs(1)),': Individual VAR: ', num2str(VAR(PCs(1))),'%'));
    ylabel(strcat('PC ', int2str(PCs(2)),': individual VAR: ', num2str(VAR(PCs(2))),'%'));
    
    if (nargin == 8)
        hold on
        DA_FLplotloadings(LDS(:,PCs),VariableNames);
        hold off
    end
    
    load('mycmap')
    colormap(mycmap)
    
    
    
    % 
    if two_plots
    S_scoresID=vertcat(S_scores(ID==IDs(l1)& Labels==1));
    S_scaledID=vertcat(S_scaled(ID2==IDs(l1)& Labels==1));
    
    [~,DensityInd,X,Y,~,~,~] = kde2dimproved(S_scoresID.Data, 1024,MIN_XY,MAX_XY);
    [DensityInd] = ProcessDensity(DensityInd);
    
    Xisc = X(1,:);
    Yisc = Y(:,1);
    
    figure
    imagesc(Xisc(1,:), Yisc(:,1), DensityInd,[0 .5]);
    set(gca,'Ydir','normal')
    colorbar;
    
    title(strcat('Density of individual: ', int2str(IDs(l1)), ' Group 1' ),'FontSize',20)
    VAR = zeros(1, N_PCs);
    for ii = 1:N_PCs
        VAR(1,ii) = 1 - sum(var(S_scaledID.Data - S_scoresID.Data(:,ii)*LDS(:,PCs(ii))'))/sum(var(S_scaledID.Data));
    end
    VAR = round((VAR*100));
    xlabel(strcat('PC ', int2str(PCs(1)),': Individual VAR: ', num2str(VAR(PCs(1))),'%'));
    ylabel(strcat('PC ', int2str(PCs(2)),': individual VAR: ', num2str(VAR(PCs(2))),'%'));
    
    if (nargin == 8)
        hold on
        DA_FLplotloadings(LDS(:,PCs),VariableNames);
        hold off
    end
    
    load('mycmap')
    colormap(mycmap)
    
    end
        
        
    
end
