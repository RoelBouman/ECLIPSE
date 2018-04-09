function [ Clim ] = PlotDiff( DensityDiff, X, Y,plotMode, PCs, VAR, LDS, VariableNames)
%% PlotDiff plots the differences between 2 densities described by
%  DensityDiff. The plots may additionally be improved upon by entering X
%  and Y for correct axes, PCs and VAR to add labels to the axes, and LDS 
%  and VariableNames to also plot loadings. Supported plot options are
%  surface and imagesc plots.
%
%% Dependencies:
%  Da_FLplotloadings.m
%
%
%% Input:
%
%  DensityDiff is a n-by-n matrix containing the local densities at each
%  point. For proper visualization, they must be centered around 0.
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
%  Clim, the computed limits based on the difference plot. (optional)
%
%% Remarks:
%  The function ends with hold off after plotting loadings. Consider this
%  when using it.
%
%% Author: Roel Bouman (roelbouman@gmail.com). Date: 3-05-2016

%% calculation of Clims
maxdiff = max(max(DensityDiff));
mindiff = min(min(DensityDiff));
climtemp = max([maxdiff abs(mindiff)]);
Clim = [-climtemp climtemp];

%% Loading cmap

%load('cmap_bluered');
% load('cmap2'); %Uses an uneven length colormap
load('blue_red_intense');
%% Plotting
if (nargin > 3)
    
    if (strcmp(plotMode, 'imagesc'))    
        
        % ensuring proper form of X and Y
        
        Xisc = X(1,:);
        
        Yisc = Y(:,1);
        
        imagesc(Xisc(1,:), Yisc(:,1), DensityDiff,Clim);
        set(gca, 'Ydir','normal');
        colormap(blue_red_intense);
%         set(gca,'XLim',[-8 10])
        colorbar;
        
        if (nargin > 4)
            xlabel(strcat('PC ', int2str(PCs(1)),': Variance explained: ', int2str(VAR(PCs(1))),'%'));
            ylabel(strcat('PC ', int2str(PCs(2)),': Variance explained: ', int2str(VAR(PCs(2))),'%'));
        end
        
        if (nargin == 8)
            hold on
            DA_FLplotloadings(LDS,VariableNames);
            hold off
        end
        
    elseif(strcmp(plotMode, 'surf')) 
        
        % based on http://nl.mathworks.com/help/matlab/ref/caxis.html
        [m,~] = size(blue_red_intense);        
        
        index = fix((DensityDiff-Clim(1,1))/(Clim(1,2)-Clim(1,1))*m)+1;
        %Clamp values outside the range [1 m]
        index(index<1) = 1;
        index(index>m) = m;
        
        cmapmod = blue_red_intense(min(min(index)):max(max(index)),:);
        
        %index = index-min(min(index))+1;        
        
        %set opacity
        alphamap = 1-double((index == round(length(blue_red_intense)/2)));
        
        h = surf(X,Y,DensityDiff,index);
        colormap(cmapmod);
        set(h,'LineStyle','none');
        alpha(alphamap)
        
    end
    
else
    
    imagesc(DensityDiff,Clim)
    colormap(blue_red_intense);
    colorbar;
   
    
end

end