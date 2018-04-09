function [densities] = LocalDensity(densitydiff,locs, X, Y)
%% LocalDensity calculates the local value of the density, or density diff, based on the location of point on a grid specified by vector X and Y. This implementation is 2D.
%
%% Dependencies:
%  binlocs2d.m
%
%% Input:
%  densitydiff, a m-by-n matrix containing densities, or density difference values evaluated at every point in the matrix.
%  locs, a p-by-2 vector indicating the location in 2D space of the points where the density should be evaluated.
%  X, the m-by-1 vector specifying the grid on the x-axis
%  Y, the n-by-1 vector specifying the grid on the y-axis
%
%% Output:
%  densities, a p-by-1 vector containing the densities evaluated at the point in locs.
%% Author: Roel Bouman (roelbouman@gmail.com).


binindex = binlocs2d(locs,X,Y);
densities = densitydiff(sub2ind(size(densitydiff),binindex(:,1),binindex(:,2)));  
    


end

