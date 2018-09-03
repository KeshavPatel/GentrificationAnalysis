close all; clc; clear;

%% Parameter Ranges
Rrange = 1.05:.01:1.15;%4-2*sqrt(2); %chosen st hotspot region lies within
Erange = 0.01:.01:.07;
Brange = 0;
Arange = 1:.5:3;

%% Constants
dt = .5;
tTotal = 200;      % total time
S = 1;
tests = {'Basic MLE'; 'Changing PDE Geometry'; 'Penalty Term'; 'Spatially Dependent Constants'};
testWhich = [1 2 3 4];

%% Fake Data
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

map = abs(150-sqrt((mx-125).^2+(my-125).^2));
map = map/sum(sum(map));

sampleInd = datasample(1:meshSize^2, 1000, 'Weights', map(:));
data = zeros(meshSize);
data(sampleInd) = 1;

Wmap = data;
Amap = data;

map = map/sum(sum((S/meshSize)^2*map));


%% Load Data
%{
Wmap = double(rgb2gray(imread('06-27_chicagoMedPropx2.png')));
Wmap = Wmap(1:525, 1:675);
Wmap(Wmap==221)=0;
Wmap(Wmap~=0)=1000;

Amap = rgb2gray(imread('06-27_chicagoAmenityPosition.png'));
Amap = Amap(1:525, 1:675);
Amap(Amap==221)=0;
Amap(Amap~=0)=1000;
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
    subplot(1,length(testWhich)+1,currSubplot)
    pdeplot(bestModelBasic,'XYData',bestuBasic(:,1,besttBasic),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuBasic(:,1,besttBasic)))])
    title({'Basic MLE', ['LSE=' num2str(LSEBasic(1))], ['R=' num2str(bestParamBasic(1)) ...
        '; E=' num2str(bestParamBasic(2)) '; B=' num2str(bestParamBasic(3)) '; A=' ...
        num2str(bestParamBasic(4)) '; t=' num2str(besttBasic*dt)]})
    currSubplot = currSubplot+1;
end

if find(testWhich == 2)
    subplot(1,length(testWhich)+1,currSubplot)
    pdeplot(bestModelGeom,'XYData',bestuGeom(:,1,besttGeom),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuGeom(:,1,besttGeom)))])
    title({'Changing PDE Geometry', ['LSE=' num2str(LSEGeom(1))], ['R=' num2str(bestParamGeom(1)) ...
        '; E=' num2str(bestParamGeom(2)) '; B=' num2str(bestParamGeom(3)) '; A=' ...
        num2str(bestParamGeom(4)) '; t=' num2str(besttGeom*dt)]})
    currSubplot = currSubplot+1;
end

if find(testWhich == 3)
    subplot(1,length(testWhich)+1,currSubplot)
    pdeplot(bestModelPenalty,'XYData',bestuPenalty(:,1,besttPenalty),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuPenalty(:,1,besttPenalty)))])
    title({'Penalty Term', ['LSE=' num2str(LSEPenalty(1))], ['R=' num2str(bestParamPenalty(1)) ...
        '; E=' num2str(bestParamPenalty(2)) '; B=' num2str(bestParamPenalty(3)) '; A=' ...
        num2str(bestParamPenalty(4)) '; t=' num2str(besttPenalty*dt)]})
    currSubplot = currSubplot+1;
end

if find(testWhich == 4)
    subplot(1,length(testWhich)+1,currSubplot)
    pdeplot(bestModelSpatDepn,'XYData',bestuSpatDepn(:,1,besttSpatDepn),'ColorBar','off');
    %colorbar('westoutside')
    caxis([0 max(max(bestuSpatDepn(:,1,besttSpatDepn)))])
    title({'Spatially Dependent Constants', ['LSE=' num2str(LSESpatDepn(1))], ['R=' num2str(bestParamSpatDepn(1)) ...
        '; E=' num2str(bestParamSpatDepn(2)) '; B=' num2str(bestParamSpatDepn(3)) '; A=' ...
        num2str(bestParamSpatDepn(4)) '; t=' num2str(besttSpatDepn*dt)]})
    currSubplot = currSubplot+1;
end

subplot(1,length(testWhich)+1,currSubplot)
plotMap = map;
plotMap(~omega) = NaN;
h = pcolor(plotMap(end:-1:1,:));
caxis([0 max(max(plotMap))])
set(h, 'EdgeColor', 'none');
title('Test Data PMF')

cusColorMap = [linspace(0,1,100)', zeros(100,1), zeros(100,1)];
colormap(cusColorMap)

saveas(gcf, [today '_MLEMethods_temp.png'])
%}
