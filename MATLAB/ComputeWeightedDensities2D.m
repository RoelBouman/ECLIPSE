function [Density1, Density2, X, Y, bandwidth1, bandwidth2, MIN_XY, MAX_XY] = ComputeWeightedDensities2D(Data_Struct, Nr_Of_Bins)
%% ComputeWeightedDensities2D computes 2D weighted Densities based on an individual scaling basis, instead of a cell based approach.
%
%% Dependencies: (concatenated at the end)
%  Weighted2DHistogram.m
%  HistKDE2D.m
%
%% Input:
%  Data_Struct, a Structure containing at least .Data and .Labels entries.
%  The .Data section should contain 2D coordinates of each cell. If more
%  than 2 dimensions are used as input, only the first 2 will be used in
%  the computation of the histogram.
%
%  Nr_Of_Bins (optional), indicating the number of bins in each dimension.
%  Due to the nature of how KDE is calculated, this should be a power of 2.
%  If Nr_Of_Bins isn't specified the default value of 2^10 will be used.
%
%% Output:
%  Density1, Density2, the densities of data1 and data2 projected onto the same meshgrid.
%  
%  X, Y, both m-by-n matrices indicating the centers of the meshgrid onto which Density1 and Density2 have been computed.
%  
%  bandwidth1 and bandwidth2, the bandwidths used for the computation of the density of Density1 and density2 respectively.
%  
%  MIN_XY and MAX_XY, indicating the minimum and maximum in the entirety of the given coordinates in both of the data matrices.
%
%% Author: Roel Bouman (roelbouman@gmail.com). Date: 25-05-2016


if nargin < 2
   Nr_Of_Bins = 2^10;    
end

[WHistControl, WHistPatient, MAX_XY, MIN_XY, NC, NP] = Weighted2DHistogram( Data_Struct,Nr_Of_Bins );

[bandwidth1,Density1,X,Y] = HistKDE2D(WHistControl,MIN_XY,MAX_XY,NC);
[bandwidth2,Density2,~,~] = HistKDE2D(WHistPatient,MIN_XY,MAX_XY,NP);



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

function [WHistControl, WHistPatient, MAX_XY, MIN_XY, NC, NP] = Weighted2DHistogram( Data_Struct,Nr_Of_Bins )
%% Weighted2DHistogram computes histograms, weighted for every individual, based on a FCS Struct.
%  The resulting histograms are specific to the Control and Patient groups,
%  as indicated by the labels.
%
%% Dependencies:
%  None
%
%% Input:
%  Data_Struct, a Structure containing at least .Data and .Labels entries.
%  The .Data section should contain 2D coordinates of each cell. If more
%  than 2 dimensions are used as input, only the first 2 will be used in
%  the computation of the histogram.
%
%  Nr_Of_Bins (optional), indicating the number of bins in each dimension.
%  Due to the nature of how KDE is calculated, this should be a power of 2.
%  If Nr_Of_Bins isn't specified the default value of 2^10 will be used.
%
%% Output:
%  WHistControl, the weighted histogram of the control group.
%
%  WHistPatient, the weighted histogram of the patient group.
%
%  MAX_XY, a 1-by-2 vector containing the maxima on the X and Y axis on
%  indices 1 and 2 respectively
%
%  MIN_XY, a 1-by-2 vector containing the minima on the X and Y axis on
%  indices 1 and 2 respectively 
%
%  NC, the number of datapoints used in the calculation of the control
%  histogram
%
%  NP, the number of datapoints used in the calculation of the patient
%  histogram
%
%% Remarks:
%  The histograms are weighted in such a manner that sum(sum(histogram)) =
%  1.
%
%  PC's are expected to be already selected, otherwise, the first 2 pc's will
%  be selected.
%
%  The labels for control and patients groups are 0 and 1 respectively
%
%% Date: 25-05-2016



if nargin<4
    Nr_Of_Bins=2^10;
end
Nr_Of_Bins=2^ceil(log2(Nr_Of_Bins)); % round up n to the next power of 2;

for i = 1:length(Data_Struct)
    
    if(size(Data_Struct(i).Data,2) > 2)
        warning('Data with more than 2 dimensions was used as input, only the first 2 dimensions are used')
        Data_Struct(i).Data = Data_Struct(i).Data(:,1:2);
    end
    
end

Data_Struct_Control = Data_Struct(vertcat(Data_Struct.Labels) == 0);
    
Data_Struct_Patient = Data_Struct(vertcat(Data_Struct.Labels) == 1);

Max_Control = zeros(length(Data_Struct_Control),2);
Max_Patient = zeros(length(Data_Struct_Patient),2);


Min_Control = zeros(length(Data_Struct_Control),2);
Min_Patient = zeros(length(Data_Struct_Patient),2);

for i = 1:length(Data_Struct_Control)
   
    Max_Control(i,:) = max(Data_Struct_Control(i).Data,[],1);
    Min_Control(i,:) = min(Data_Struct_Control(i).Data,[],1);
    
end


for i = 1:length(Data_Struct_Patient)
   
    Max_Patient(i,:) = max(Data_Struct_Patient(i).Data,[],1);
    Min_Patient(i,:) = min(Data_Struct_Patient(i).Data,[],1);
    
end

MAX_XY = max([Max_Control;Max_Patient],[],1);
MIN_XY = min([Min_Control;Min_Patient],[],1);

Range=MAX_XY-MIN_XY;

MAX_XY=MAX_XY+Range/4; 
MIN_XY=MIN_XY-Range/4;

% Transformation for regular binning, code from: kde2d.m by Z. I. Botev, J. F. Grotowski and D. P. Kroese
scaling=MAX_XY-MIN_XY;

NC = 0;
NP = 0;

if(~isempty(Data_Struct_Control))
    
    WHistControl = zeros(Nr_Of_Bins);
    
    for i = 1:length(Data_Struct_Control)
        
        N=size(Data_Struct_Control(i).Data,1);
        NC = NC + N;
        transformed_data_Control=(Data_Struct_Control(i).Data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);
        %bin the data uniformly using regular grid;
        WHistControl = WHistControl + ndhist(transformed_data_Control,Nr_Of_Bins);
        
    end
    
else
    
    WHistControl = [];
    
end

WHistControl = WHistControl/length(Data_Struct_Control);

if(~isempty(Data_Struct_Patient))
    
    WHistPatient = zeros(Nr_Of_Bins);
    
    for i = 1:length(Data_Struct_Patient)
        
        N=size(Data_Struct_Patient(i).Data,1);
        NP = NP + N;
        transformed_data_Patient=(Data_Struct_Patient(i).Data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);
        %bin the data uniformly using regular grid;
        WHistPatient = WHistPatient + ndhist(transformed_data_Patient,Nr_Of_Bins);
        
    end
    
    WHistPatient = WHistPatient/length(Data_Struct_Patient);
    
else
    
    WHistPatient = [];
    
end

end

% %% From kde2d.m by Z. I. Botev, J. F. Grotowski and D. P. Kroese
% function binned_data=ndhist(data,M)
% % this function computes the histogram
% % of an n-dimensional data set;
% % 'data' is nrows by n columns
% % M is the number of bins used in each dimension
% % so that 'binned_data' is a hypercube with
% % size length equal to M;
% [nrows,ncols]=size(data);
% bins=zeros(nrows,ncols);
% for i=1:ncols
%     [~,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
%     bins(:,i) = min(bins(:,i),M);
% end
% % Combine the  vectors of 1D bin counts into a grid of nD bin
% % counts.
% binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
% 
% 
% 
% end

%%

%%%%%%%

function [bandwidth,density,X,Y]=HistKDE2D(binned_data,MIN_XY,MAX_XY, N)
%% HistKDE2D computes the 2D KDE based on a histogram, its scale, and the number of points it was based on.
%
%% Dependencies:
%  None
%
%% Input:
%  Binned Data, a n-by-n histogram containing bin counts. The data is
%  scaled to have a total sum of 1 across the entire histogram for
%  normalization purposes.
%
%  MIN_XY, a 1-by-2 vector containing the minima on the X and Y axis on
%  indices 1 and 2 respectively 
%
%  MAX_XY, a 1-by-2 vector containing the maxima on the X and Y axis on
%  indices 1 and 2 respectively
%
%  N, the number of datapoints used to calculate the histogram.
%
%% Output:
%  Bandwidth, the bandwidth used for the KDE estimation.
%
%  Density, the KDE based probability density across the histogram grid.
%
%  X, a n-by-n matrix containing the X coordinates at a density evaluation
%  point.
%  
%  Y, a n-by-n matrix containing the Y coordinates at a density evaluation
%  point.
%
%% Remarks:
%  The KDE implementation is based upon the original KDE2D method described
%  in:
%  Reference: Z. I. Botev, J. F. Grotowski and D. P. Kroese
%             "KERNEL DENSITY ESTIMATION VIA DIFFUSION" ,Submitted to the
%             Annals of Statistics, 2009
%
%  The largest part is a direct adaptation of this original implementation.
%
%% Author: Roel Bouman (roelbouman@gmail.com)
%% Date: 25-05-2016

[n1, n2] = size(binned_data);

scaling=MAX_XY-MIN_XY;
    
if(n1 ~= n2)
    error('binned_data must be a square matrix');
elseif(n1 ~= 2^ceil(log2(n1)))
    error('the size along both dimensions must be a power of 2')
elseif(n1 == n2 && n1 == 2^ceil(log2(n1)))
    n = n1;
end

% Ensure the binned_data is scaled properly
initial_data = binned_data./(sum(sum(binned_data)));

% discrete cosine transform of initial data
a= dct2d(initial_data);

% now compute the optimal bandwidth^2
I=(0:n-1).^2; A2=a.^2;

test = @(t)(t-evolve(t, N, I, A2));

t_star=fzero( test,[0,0.1]);

p_02=func([0,2],t_star, N, I, A2);p_20=func([2,0],t_star, N, I, A2); p_11=func([1,1],t_star, N, I, A2);
t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);

% smooth the discrete cosine transform of initial data using t_star
a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 

% now apply the inverse discrete cosine transform
if nargout>1
    density=idct2d(a_t)*(numel(a_t)/prod(scaling));
    [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
end

bandwidth=sqrt([t_x,t_y]).*scaling; 

end

%%
function  [out,time]=evolve(t, N, I, A2)

Sum_func = func([0,2],t, N, I, A2) + func([2,0],t, N, I, A2) + 2*func([1,1],t, N, I, A2);

time=(2*pi*N*Sum_func)^(-1/3);

out=(t-time)/time;

end

%%
function out=func(s ,t , N, I, A2)

if (sum(s)<=4)
    
    Sum_func=func([s(1)+1,s(2)],t, N, I, A2)+func([s(1),s(2)+1],t, N, I, A2); const=(1+1/2^(sum(s)+1))/3;
    time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
    out=psi(s,time, I, A2);
    
else
    out=psi(s,t, I, A2);
end

end

%%
function out=psi(s,Time, I, A2)
% s is a vector

w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];

wx=w.*(I.^s(1));
wy=w.*(I.^s(2));

out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));

end

%%
function out=K(s)

out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);

end

%%
function data=dct2d(data)

% computes the 2 dimensional discrete cosine transform of data
% data is an nd cube
[nrows,ncols]= size(data);
if nrows~=ncols
    error('data is not a square array!')
end
% Compute weights to multiply DFT coefficients
w = [1;2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
weight=w(:,ones(1,ncols));
data=dct1d(dct1d(data, weight)', weight)';


end

%%
function transform1d=dct1d(x, weight)

% Re-order the elements of the columns of x
x = [ x(1:2:end,:); x(end:-2:2,:) ];

% Multiply FFT by weights:
transform1d = real(weight.* fft(x));

end

%%
function data = idct2d(data)
% computes the 2 dimensional inverse discrete cosine transform

[nrows,ncols]=size(data);

% Compute weights
w = exp(1i*(0:nrows-1)*pi/(2*nrows)).';

weights=w(:,ones(1,ncols));

data=idct1d(idct1d(data, weights, nrows, ncols)', weights, nrows, ncols);

end

%%
function out=idct1d(x, weights, nrows, ncols)

y = real(ifft(weights.*x));

out = zeros(nrows,ncols);

out(1:2:nrows,:) = y(1:nrows/2,:);
out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);

end

%%
function [binned_data,bins]=ndhist(data,M)
% this function computes the histogram
% of an n-dimensional data set;
% 'data' is nrows by n columns
% M is the number of bins used in each dimension
% so that 'binned_data' is a hypercube with
% size length equal to M;

[nrows,ncols]=size(data);

bins=zeros(nrows,ncols);

for i=1:ncols
    [~,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
    bins(:,i) = min(bins(:,i),M);
end

% Combine the  vectors of 1D bin counts into a grid of nD bin counts.
binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));

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

