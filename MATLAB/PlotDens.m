function PlotDens(Density,X,Y,plotMode, PCs, VAR, LDS, VariableNames)
%% PlotDens plots the density as either an imagesc or surf plot.
%
%% Dependencies:
%  Da_FLplotloadings.m
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
%% Remarks:
%  The function ends with hold off after plotting loadings. Consider this
%  when using it.
%% Author: Roel Bouman (roelbouman@gmail.com).
 
%% Plotting
load mycmap.mat

if (nargin > 2)
    
    if (strcmp(plotMode, 'imagesc')||nargin < 4)    
        
        % ensuring proper form of X and Y
        
        Xisc = X(1,:);
        
        Yisc = Y(:,1);
        
        imagesc(Xisc(1,:), Yisc(:,1), Density);
        set(gca,'Ydir','normal')
        colorbar;
        
        if (nargin > 4)
            xlabel(strcat('PC ', int2str(PCs(1)),': Variance explained: ', int2str(VAR(PCs(1))),'%'));
            ylabel(strcat('PC ', int2str(PCs(2)),': Variance explained: ', int2str(VAR(PCs(2))),'%'));
        end
        
        if (nargin == 8)
            hold on
            DA_FLplotloadings(LDS(:,PCs),VariableNames);
            hold off
        end
        
    elseif(strcmp(plotMode, 'surf')) 
        
        h = surf(X,Y,Density);
        set(h,'LineStyle','none');
        
    end
    
else
    
    imagesc(Density)
    set(gca,'Ydir','normal')
    colorbar;
    colormap(mycmap)
end

end
