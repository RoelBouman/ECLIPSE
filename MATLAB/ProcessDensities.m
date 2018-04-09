function [Density1, Density2, DensityDiff, iclDensityDiff, LogDensity1, LogDensity2] = ProcessDensities(Density1, Density2)
%% ProcessDensities processes densities computed using ComputeDensities2D
%
%% Dependencies:
% ProcessDensity.m
%
%% Input:
%  Density1 and Density2, both m-by-n matrices of densities computed upon the same grid
%
%% Output:
%  Density1 and Density2, both processed to exclude small negative values and values below machine precision.
%  Densitydiff, the difference between Density1 and Density2 by subtraction of Density2 from density1.
%  iclDensityDiff, the inverse cutoff log scaled DensityDiff, allowing for log scaling around zero, whilst retaining useful information.
%  LogDensity1 and Logdensity2, the log scaled Density1 and Density2.
%
%% Author: Roel Bouman (roelbouman@gmail.com). Date: 3-05-2016


[Density1, LogDensity1 ] = ProcessDensity( Density1 );
[Density2, LogDensity2 ] = ProcessDensity( Density2 );


DensityDiff = Density1 - Density2;

%inverse cutoff log scaling
iclDensityDiff = zeros(size(DensityDiff));
[m,n] = size(DensityDiff);
for i = 1:m
    for j = 1:n
        if(DensityDiff(i,j)>0)
            
            iclDensityDiff(i,j) = -1/log(DensityDiff(i,j));
            
        elseif(DensityDiff(i,j)<0)
            
            iclDensityDiff(i,j) = 1/(log(-DensityDiff(i,j)));
            
        end
           
    end
end

end

%%

function [Density, LogDensity ] = ProcessDensity( Density )
% ProcessDensity returns the Density with small values below machine precision error, and negative values, removed. In addition, the log scaled Density is also returned.
%
% Dependencies:
%  None
%
% Input:
%  Density, a m-by-n matrix containing the densities evaluated at every point.
%
% Output:
%  Density, a m-by-n matrix containing the densities, but with small values below machine precision error, and negative values removed.
%  logDensity, a m-by-n matrix containing the log scaled processed Density.

Density(Density <10^-5) = 0;
LogDensity = log(Density);

end
