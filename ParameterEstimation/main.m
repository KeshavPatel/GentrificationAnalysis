clc; clear; close all;

%% Fake Data
% Use this to create artificial wealth and amenity data. Change the
% variable "omega" to define valid and invalid regions (1 is valid, 0 is
% invalid). Change the variable "map" to define a probability density map
% of wealth and amenities. The points used in the MLE are randomly sampled
% using this distribution.
%{
meshSize = 250;
boundaryDim1 = 1/10;
boundaryDim2 = 1/10;
[mx, my] = meshgrid(0:meshSize-1);
PDEToDatax = @(x) 1+floor((x+S/2)*(meshSize-1)/S);

omega = sqrt((mx-meshSize/2).^2+(my-meshSize/2).^2)>=10 & sqrt((mx-meshSize/2).^2+(my-meshSize/2).^2)<=meshSize/2;
nzNeighbor = omega(1:end-2,1:end-2)+omega(2:end-1,1:end-2)+omega(3:end,1:end-2)+omega(1:end-2,2:end-1)+ ...
    omega(3:end,2:end-1)+omega(1:end-2,3:end)+omega(2:end-1,3:end)+omega(3:end,3:end);
boundary = zeros(meshSize);
boundary(2:end-1,2:end-1) = (nzNeighbor ~= 8 & nzNeighbor ~= 0) & ...
    omega(2:end-1,2:end-1)==0;

map = abs(150-sqrt((mx-125).^2+(my-125).^2)) .* omega;
map = map/sum(sum(map));

sampleInd = datasample(1:meshSize^2, 1000, 'Weights', map(:));
data = zeros(meshSize);
data(sampleInd) = 1;

Wmap = data;
Amap = data;

map = map/sum(sum((S/meshSize)^2*map));
%}

%% Real Data
%%{
meshSize = 100;
boundaryDim1 = 1/10;
boundaryDim2 = 1/10;
[mx, my] = meshgrid(0:meshSize-1);
PDEToDatax = @(x) 1+floor((x+S/2)*(meshSize-1)/S);
omega = ones(meshSize);
nzNeighbor = omega(1:end-2,1:end-2)+omega(2:end-1,1:end-2)+omega(3:end,1:end-2)+omega(1:end-2,2:end-1)+ ...
    omega(3:end,2:end-1)+omega(1:end-2,3:end)+omega(2:end-1,3:end)+omega(3:end,3:end);
boundary = zeros(meshSize);
boundary(2:end-1,2:end-1) = (nzNeighbor ~= 8 & nzNeighbor ~= 0) & ...
    omega(2:end-1,2:end-1)==0;

%{
Wmap = double(rgb2gray(imread('6-17_chicagoMedPropVal.png')));
Wmap = Wmap(1:2:2*meshSize,50:2:49+2*meshSize);
Amap = double(rgb2gray(imread('6-17_chicago_Amenities.png')));
Amap = Amap(1:2:2*meshSize,50:2:49+2*meshSize);
%}
Wmap = double(rgb2gray(imread('6-17_chicagoMedPropVal.png')));
Wmap = Wmap(floor(linspace(1,size(Wmap,1),meshSize)),floor(linspace(50,550,meshSize)));
Amap = double(rgb2gray(imread('6-17_chicago_Amenities.png')));
Amap = Amap(floor(linspace(1,size(Wmap,1),meshSize)),floor(linspace(50,550,meshSize)));
    
map = Wmap;
%}

Rrange = [1 4-2*sqrt(2)];
Erange = [0 .5];
Brange = [0 .01];
Arange = [1 2];

timeSteps = 8000;

fun = @(X) likelihood(X(1), X(2), X(3), X(4), timeSteps, Wmap, Amap);

Y = fmincon(fun,[1.01 .05 .002 1.1],[],[],[],[], ...
    [Rrange(1) Erange(1) Brange(1) Arange(1)],[Rrange(2) Erange(2) Brange(2) Arange(2)]);

%%
[result,model,d,ICVec] = continuumSimulation(Y(1), Y(2), Y(3), 0.5, timeSteps, 1, Y(4), []);
u = result.NodalSolution;

%%
figure(1)
subplot(1,2,1)
pdeplot(model,'XYData',u(:,1,end),'ColorBar','off');
colorbar('northoutside')
colormap jet
title({'Basic MLE', ['R=' num2str(Y(1)) '; E=' num2str(Y(2)) ...
    '; B=' num2str(Y(3)) '; A=' num2str(Y(4))]})

subplot(1,2,2)
plotMap = map; 
plotMap(~omega) = NaN;
h = pcolor(plotMap(end:-1:1,:));
colorbar('northoutside')
colormap jet
caxis([0 max(max(plotMap))])
set(h, 'EdgeColor', 'none');
title('Test Data PMF')
