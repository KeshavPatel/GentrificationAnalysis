%% Energy Minimization via fmincon
%{
Last Modified: 7/30/2020

Description: This script uses fmincon to compute an energy functional
decribing the least-squares difference between the steady-state PDE
solution and a given set of data. The energy functional can be thought of
as a function of four parameters from the PDE. The output of the script is
the steady-state PDE solution that mininimizes the energy and a figure
depicting the solution with the given data.
%}

clc; clear; close all;

%% Setup

%Parameter Domain
Rrange = [1 4-2*sqrt(2)];
Erange = [0 .5];
Brange = [0 .1];
Arange = [1 10000];

%Numerical Parameters
S = 1;
timeDilate = 0.00001;

%load data
load('R_1pt010029333613098_1pt108595084546357_E_1em4_B_1em8_A_10.mat')

%Find minima and maxima of dataset
Wset = griddedInterpolant(xmesh, result(1,:));
Aset = griddedInterpolant(xmesh, result(2,:));
xData = []; WData = []; AData = [];

xIter = linspace(-.5,.5,10001);
sgn = Wset(xIter(2))-Wset(xIter(1)) > 0;
for j = 2:10000
    sgnNow = Wset(xIter(j+1))-Wset(xIter(j)) > 0; %1 if slope of Wset is positive, 0 otherwise
    if sgnNow ~= sgn
        sgn = sgnNow;
        xData = [xData; xIter(j)];
        WData = [WData; Wset(xIter(j))];
        AData = [AData; Aset(xIter(j))];
    end
end

%% START Steady State PDE Solver
if isfile('temp.mat')
    delete('temp.mat')
end

%as of 7/30/2020 this program computes a 2D minimization. The definition of
%'X' should be updated here and in the next section to add/remove/change
%which parameters undergo minimization
fun = @(X) likelihood_1DSS_LS([X(1) X(2) 1e-4 1e-4 1e-8 1e-8 10 10], [xData WData AData]);
opts = optimoptions('fmincon', 'PlotFcn', 'optimplotfval', 'algorithm', 'sqp', 'DiffMaxChange', .01);

[Y, fval] = fmincon(fun, [1 1], [],[],[],[], [Rrange(1) Rrange(1)], [Rrange(2) Rrange(2)], [],opts);

%% Compute Minimizer Solution
X = [Y(1) Y(2) 1e-4 1e-4 1e-8 1e-8 10 10];

rho=@(x) X(1)+(X(2)-X(1))*heaviside(x); eta=@(x) X(3)/timeDilate+(X(4)-X(3))/timeDilate*heaviside(x);
gamma=@(x) X(5)+(X(6)-X(5))*heaviside(x); alpha=@(x) X(7)/timeDilate+(X(8)-X(7))/timeDilate*heaviside(x);
[xmesh, result_dt] = continuumSimulation_1D(rho,eta,gamma,S,alpha,10001,[]);
eta=@(x) X(3)+(X(4)-X(3))*heaviside(x); alpha=@(x) X(7)+(X(8)-X(7))*heaviside(x);
result = steadyStateSimulation_1D(rho,eta,gamma,S,alpha,{xmesh, squeeze(result_dt(end,:,:))},2);


%% Plot
figure(2)
subplot(2,1,1)
plot(result.x, result.y(1,:))
hold on
scatter(xData,WData)
hold off
set(gca, 'fontsize', 18); set(gca, 'tickLabelInterpreter', 'latex');
xlabel('$x$', 'fontsize', 20, 'interpreter', 'latex')
ylabel('$\bar W$', 'fontsize', 20, 'interpreter', 'latex')
title({['Energy: ' num2str(fval)]; ['initial pos: $\rho_{\rm left}=1.1$; $\rho_{\rm right}=1.02$'];
    ['fmincon soln: $\rho_{\rm left}=$' num2str(Y(1),10) '; $\rho_{\rm right}=$' num2str(Y(2),10)];
    ['Data: $\rho_{\rm left}=$1.010029333613098; $\rho_{\rm right}=$1.108595084546357; $\eta=$1e-4' ...
    '; $\gamma=$1e-8; $\alpha=$10']}, ...
    'fontsize', 20, 'interpreter', 'latex')
ylim([0 1.1*max([result.y(1,:) WData'])])

subplot(2,1,2)
plot(result.x, result.y(3,:))
hold on
scatter(xData,AData)
hold off
set(gca, 'fontsize', 18); set(gca, 'tickLabelInterpreter', 'latex');
xlabel('$x$', 'fontsize', 20, 'interpreter', 'latex')
ylabel('$\bar A$', 'fontsize', 20, 'interpreter', 'latex')
ylim([0 1.1*max([result.y(3,:)./Aintegral AData'./ADataIntegral])])
saveas(gcf,['figures/2D_minimization_MaxMinDataPts_temp.png'])


