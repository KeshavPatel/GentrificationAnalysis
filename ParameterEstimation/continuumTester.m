%% PDE Simulation Tester
%{
Author:         Avishai Halev
Advisor:        Nancy Rodriguez-Bunn
Last modified:  Keshav Patel
Date modified:  09/04/2018

This file runs the "continuumSimulation.m" file and plots the solution. A
movie is saved at completion of the script
%}

close all
clearvars -except oldResult pResult


%% Define parameters
dt = 0.5;
tTotal = 200;      % total time
R = 1.1;
E = 0.0;
B = 0.001;
A = 2;
S = 1;

%% Load Data as Initial Conditions
%%{
IC1 = double(rgb2gray(imread('6-17_chicagoMedPropVal.png')));
IC1 = IC1(1:250, 51:300);
IC1(IC1==221)=0;
IC2 = double(rgb2gray(imread('6-17_chicago_Amenities.png')));
IC2 = IC2(1:250, 51:300);
IC2(IC2==221)=0;
initialConditionCell = {IC1; IC2};
%}

%% Model Condiitons
meshSize = 250;
[mx, my] = meshgrid(1:250);
PDEToDatax = @(x) 1+floor((x+S/2)*(meshSize-1)/S);
omega = (mx-meshSize/2).^2+(my-meshSize/2).^2>=25;
boundaryDim1 = 1/10;
boundaryDim2 = 1/10;
rho = @(x,y) R;
eta = @(x,y) E;
alpha = @(x,y) A;
beta = @(x,y) B;

timeSteps = tTotal/dt+1; %Defines timeSteps based off tTotal and dt above

%% Run and time simulation

tic
if ~exist('initialConditionCell','var')
    fprintf('Using default initial condition\n')
    try
        [result,model,d,ICVec,time] = continuumSimulation(rho, eta, beta, dt, ...
            timeSteps, S, alpha, []);
    catch
        
    end
else
    fprintf('Using old result as initial condition\n')
    [result,model,d,ICVec,time] = continuumSimulation(rho, eta, beta, dt, ...
        timeSteps, S, alpha, initialConditionCell);
end
toc

%% If running on the cluster
%save('gent1Step.mat','result','model','d','timeSteps','R','A','S','E','gamma','-v7.3')


%% Getting results
u = result.NodalSolution;
j = 1;


%% This is for debugging or trying to figure out WTF is happening
%for k=floor(0.2*timeSteps):(0.5*timeSteps)
%    fprintf('Absolute difference %4.2e relative difference %4.2f mean %4.2e\n',...
%        max(u(:,1,k))-min(u(:,1,k)),max(u(:,1,k))/min(u(:,1,k)),mean(u(:,1,k)));
    %disp(max(u(:,1,k))-min(u(:,1,k)));
%end

%% Plotting

%Plots every step timeSteps
step = 1;
fig = figure('pos',[1 387 1279 311]);
for i=1:step:time
    subplot(1,2,1)
    pdeplot(model,'XYData',u(:,1,i),'ColorMap','jet','ColorBar','off');
    colorbar('westoutside')
    title({'Wealth'; ['R = ' num2str(R) '; E = ' num2str(E) ...
        '; B = ' num2str(B) '; A = ' num2str(A)]})
    xlabel 'X-coordinate, meters'
    ylabel 'Y-coordinate, meters'
    ax = gca;                   % fixing color limits.
    ax.CLim = [0 .2];
    xlim(d(1,:))
    ylim(d(2,:))
    zlim([0 .2])
    
    subplot(1,2,2)
    pdeplot(model,'XYData',u(:,2,i),'ColorMap','jet','ColorBar','off');
    colorbar('westoutside')
    title({'Amenities'; ['R = ' num2str(R) '; E = ' num2str(E) ...
        '; B = ' num2str(B) '; A = ' num2str(A)]})
    xlabel 'X-coordinate, meters'
    ylabel 'Y-coordinate, meters'
    ax = gca;                   % fixing color limits.
    ax.CLim = [0 .2];
    xlim(d(1,:))
    ylim(d(2,:))
    zlim([0 .2])
    Mv(j)=getframe(fig);
    j = j+1;
    pause(.1)
end

%% Saving a movie. This can take a long time


v = VideoWriter(['../../Fall2018/November/6Nov/movie_R_' num2str(round(1000*R)) ...
    '_E_' num2str(round(10000*E)) '_B_' num2str(round(1000*B)) '_A_' num2str(A)]);
v.FrameRate = 5;
v.Quality = 100;
open(v)
writeVideo(v,Mv)
close(v)


