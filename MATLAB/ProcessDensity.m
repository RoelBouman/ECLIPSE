function [Density, LogDensity ] = ProcessDensity( Density )
%% ProcessDensity returns the Density with small values below machine precision error, and negative values, removed. In addition, the log scaled Density is also returned.
%
%% Dependencies:
%  None
%
%% Input:
%  Density, a m-by-n matrix containing the densities evaluated at every point.
%
%% Output:
%  Density, a m-by-n matrix containing the densities, but with small values below machine precision error, and negative values removed.
%  logDensity, a m-by-n matrix containing the log scaled processed Density.
%
%% Author: Roel Bouman (roelbouman@gmail.com)
%% Date: 3-05-2016

Density(Density <10^-5) = 0;
LogDensity = log(Density);

end

