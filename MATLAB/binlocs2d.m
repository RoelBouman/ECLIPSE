function [binindex]=binlocs2d(data,X,Y)
%% binlocs2d bins data based on a meshgrid indicated by X and Y vectors to determine 
%  the locations of the data on the grid. This implements this method in 2d space.
%
%% Dependencies:
%  none
%
%% Input:
%  data must be a n-by-2 matrix containing the coordinates in 2d space
%  X must be a m-by-1 vector containing the centers of the meshgrid on the x axis
%  Y must be a p-by-1 vector containing the centers of the meshgrid on the y axis
%
%% Output:
%  binindex, a n-by-2 vector containing the locations of the coordinates in data on the %  grid. Every row has coordinates structured as follows: [ygridcoord, xgridcoord]
%
%% Author: Roel Bouman (roelbouman@gmail.com)
%% Date: 19-05-2016
% X = X(:);
% Y = Y(:);

binindex = zeros(length(data),2);

[~,binindex(:,1)] = histc(data(:,1),Y);
[~,binindex(:,2)] = histc(data(:,2),X);

binindex(binindex < 1) = 1;
binindex(binindex > length(X)) = length(X);




end

