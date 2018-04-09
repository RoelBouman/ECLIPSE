function [labels] = DensDiffThreshold(densities,threshold)
%% DensDiffThreshold calculates labels indicating whether or not the densities of coordinates are above or below a certain threshold. 1 indicates that the density value at a certain point is higher than or equal to the density, 0 indicates that a certain density at a point is lower than the threshold.
%
%% Dependencies:
%  none
%
%% Input:
%  Densities, a n-by-1 vector containing the densities, most often based on a density difference plot.
%  Threshold, a constant indicating the threshold to compare to.
%
%% Output:
%  labels, a n-by-1 logical vector containing the result of comparison to the threshold.
%% Author: Roel Bouman (roelbouman@gmail.com).
%% Date: 19-05-2016


labels = (densities > threshold);


end


