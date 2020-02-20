clc; clear; close all;

%% Fake Data
Wset = [.4];
Aset = Wset;

Rrange = [1 4-2*sqrt(2)];
Erange = [0 .5];
Brange = [0 .1];
Arange = [0 6];

timeSteps = 5000;

%% FULL START
%%{
fun = @(X) likelihood_1D(X(1), X(2), X(3), X(4), timeSteps, Wset, Aset);
opts = optimoptions('fmincon', 'PlotFcn', 'optimplotfval');

Y = fmincon(fun,[1.04 0.1 0.0001 4], ...
    [],[],[],[], [Rrange(1) Erange(1) Brange(1) Arange(1)],[Rrange(2) Erange(2) Brange(2) Arange(2)], ...
    [],opts);

%%
result = continuumSimulation_1D(Y(1), Y(2), Y(3), 1, timeSteps*10, 1, Y(4), []);

%%
figure(2)
subplot(1,2,1)
Wintegral = 1/200*sum(result(end,:,1));
plot(linspace(-.5,.5,200), result(end,:,1)./Wintegral)
title(['R=' num2str(Y(1)) '; E=' num2str(Y(2)) ...
    '; B=' num2str(Y(3)) '; A=' num2str(Y(4))])
%ylim([0 .2])

subplot(1,2,2)
scatter(Wset,ones(size(Wset)))
xlim([-.5 .5])
ylim([0 1.1])
title('Test Data PMF')
%}

%% 1D Parameter Search START
%{
fun = @(X) likelihood_1D(X, .2, 0, 1, timeSteps, Wset, Aset);
opts = optimoptions('fmincon', 'PlotFcn', 'optimplotfval');

Y = fmincon(fun,1.01, ...
    [],[],[],[], Rrange(1),Rrange(2), ...
    [],opts);

%%
result = continuumSimulation_1D(Y, .2, 0, 1, timeSteps*10, 1, 1, []);

%%
figure(2)
subplot(1,2,1)
Wintegral = 1/200*sum(result(end,:,1));
plot(linspace(-.5,.5,200), result(end,:,1)./Wintegral)
title(['R=' num2str(Y(1)) '; E=' num2str(.2) ...
    '; B=' num2str(0) '; A=' num2str(1)])
%ylim([0 .2])

subplot(1,2,2)
scatter(Wset,ones(size(Wset)))
xlim([-.5 .5])
ylim([0 1.1])
title('Test Data PMF')
%}