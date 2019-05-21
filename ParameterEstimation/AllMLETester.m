%% Full MLE Tester
%{
Author:         Keshav Patel
Advisor:        Nancy Rodriguez-Bunn
last modified:  Keshav Patel
Date modified:  5/20/2019

This file is the main file for running our various MLE schemes. The
determination of which tests are run is given by the variable "testWhich",
and the best solution for each scheme is given in a single figure.

The most recent update adds support for using data taken from
DataCollection as the initial condition of the PDE
%}

close all; clc; clear;

%% Select Tests
% Input into "testWhich" the indices of the strings given in "tests" to use
% this test (i.e. add 2 if you want to test the change in PDE geometry)
tests = {'Basic MLE'; 'Changing PDE Geometry'; 'Penalty Term'; 'Spatially Dependent Constants'};
testWhich = [1 2 3 4];

%% Parameter Ranges
% The script will loop through every combination of the following 4 values
% when computing the solution to the PDE system

%R = $\rho$, the ratio of investment to decay
Rrange = 1:.01:4-2*sqrt(2); %chosen st hotspot region lies within
%E = $\eta$, the strength of neighborhood effects
Erange = 0:.01:.2;
%B = $\gamma$, homogeneous rate of investment in amenities
Brange = 0:.001:.005;
%A = reaction scaling term (effectively increases domain size)
Arange = 1:.5:3;

%% Constants
dt = .5; %Temporal step size
tTotal = 200; %total time
S = 1; %size of geometry - effective size is A*S

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
meshSize = 250;
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

Wmap = double(rgb2gray(imread('6-17_chicagoMedPropVal.png')));
Wmap = Wmap(1:2:2*meshSize,50:2:49+2*meshSize);
Amap = double(rgb2gray(imread('6-17_chicago_Amenities.png')));
Amap = Amap(1:2:2*meshSize,50:2:49+2*meshSize);

map = Wmap;
%}

%% Load Data as Initial Conditions
%This is a predictive test. Uncomment this to use a predefined image as the
%initial conditions for wealth and amenities. Use this to predict what will
%happen to wealth and amenities given a set of parameters defined above
%{
IC1 = double(rgb2gray(imread('6-17_chicagoMedPropVal.png')));
IC1 = IC1(1:250, 51:300);
IC1(IC1==221)=0;
IC2 = double(rgb2gray(imread('6-17_chicago_Amenities.png')));
IC2 = IC2(1:250, 51:300);
IC2(IC2==221)=0;
initialConditionCell = {IC1; IC2};
%}

%% Run Tests
if find(testWhich == 1)
    %Normal MLE
    disp('Basic MLE')
    run MLETester
    bestModelBasic = bestModel;
    bestuBasic = bestu;
    bestParamBasic = bestParam;
    besttBasic = bestt;
    LSEBasic = [WbestLSE AbestLSE];
end

if find(testWhich == 2)
    %Changing PDE Geometry
    disp('Testing Geometry')
    run MLETester_geometry
    bestModelGeom = bestModel;
    bestuGeom = bestu;
    bestParamGeom = bestParam;
    besttGeom = bestt;
    LSEGeom = [WbestLSE AbestLSE];
end

if find(testWhich == 3)
    %Penalty MLE
    disp('Testing Penalty')
    run MPLETester
    bestModelPenalty = bestModel;
    bestuPenalty = bestu;
    bestParamPenalty = bestParam;
    besttPenalty = bestt;
    LSEPenalty = [WbestLSE AbestLSE];
end

if find(testWhich == 4)
    %Spatially Dependent Const
    disp('Testing Spatially Dependent Constants')
    run MLETesterHeteroConst
    bestModelSpatDepn = bestModel;
    bestuSpatDepn = bestu;
    bestParamSpatDepn = bestParam;
    besttSpatDepn = bestt;
    LSESpatDepn = [WbestLSE AbestLSE];
end
%clearvars -except map omega bestModelGeom bestuGeom bestParamGeom besttGeom LSEGeom bestModelPenalty bestuPenalty bestParamPenalty besttPenalty LSEPenalty bestModelSpatDepn bestuSpatDepn bestParamSpatDepn besttSpatDepn LSESpatDepn

%% Plotting
%%{
today = datestr(now, 'mm-dd');
figure('pos',[1 387 1279 311])
currSubplot = 1;
if find(testWhich == 1)
    subplot(1,length(testWhich)+2,currSubplot)
    pdeplot(bestModelBasic,'XYData',bestuBasic(:,1,besttBasic),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuBasic(:,1,besttBasic)))])
    title({'Basic MLE', ['LSE=' num2str(LSEBasic(1))], ['R=' num2str(bestParamBasic(1)) ...
        '; E=' num2str(bestParamBasic(2)) '; B=' num2str(bestParamBasic(3)) '; A=' ...
        num2str(bestParamBasic(4)) '; t=' num2str(besttBasic*dt)]})
    currSubplot = currSubplot+1;
end

if find(testWhich == 2)
    subplot(1,length(testWhich)+2,currSubplot)
    pdeplot(bestModelGeom,'XYData',bestuGeom(:,1,besttGeom),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuGeom(:,1,besttGeom)))])
    title({'Changing PDE Geometry', ['LSE=' num2str(LSEGeom(1))], ['R=' num2str(bestParamGeom(1)) ...
        '; E=' num2str(bestParamGeom(2)) '; B=' num2str(bestParamGeom(3)) '; A=' ...
        num2str(bestParamGeom(4)) '; t=' num2str(besttGeom*dt)]})
    currSubplot = currSubplot+1;
end

if find(testWhich == 3)
    subplot(1,length(testWhich)+2,currSubplot)
    pdeplot(bestModelPenalty,'XYData',bestuPenalty(:,1,besttPenalty),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuPenalty(:,1,besttPenalty)))])
    title({'Penalty Term', ['LSE=' num2str(LSEPenalty(1))], ['R=' num2str(bestParamPenalty(1)) ...
        '; E=' num2str(bestParamPenalty(2)) '; B=' num2str(bestParamPenalty(3)) '; A=' ...
        num2str(bestParamPenalty(4)) '; t=' num2str(besttPenalty*dt)]})
    currSubplot = currSubplot+1;
end

if find(testWhich == 4)
    subplot(1,length(testWhich)+2,currSubplot)
    pdeplot(bestModelSpatDepn,'XYData',bestuSpatDepn(:,1,besttSpatDepn),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuSpatDepn(:,1,besttSpatDepn)))])
    title({'Spatially Dependent Constants', ['LSE=' num2str(LSESpatDepn(1))], ['R=' num2str(bestParamSpatDepn(1)) ...
        '; E=' num2str(bestParamSpatDepn(2)) '; B=' num2str(bestParamSpatDepn(3)) '; A=' ...
        num2str(bestParamSpatDepn(4)) '; t=' num2str(besttSpatDepn*dt)]})
    currSubplot = currSubplot+1;
end

%data plot
subplot(1,length(testWhich)+2,currSubplot)
plotMap = map;
plotMap(~omega) = NaN;
h = pcolor(plotMap(end:-1:1,:));
caxis([0 max(max(plotMap))])
set(h, 'EdgeColor', 'none');
title('Test Data PMF')
currSubplot = currSubplot+1;

%IC plot
subplot(1,length(testWhich)+2,currSubplot)
pdeplot(bestModelBasic,'XYData',bestuBasic(:,1,1),'ColorBar','off');
caxis([0 max(max(bestuBasic(:,1,1)))])
title('IC')

cusColorMap = [linspace(0,1,100)', zeros(100,1), zeros(100,1)];
colormap(cusColorMap)

saveas(gcf, [today '_MLEMethods_temp.png'])
%}
