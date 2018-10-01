%% Basic MLE Maximizer - Heterogeneous Parameters
%{
Author:         Keshav Patel
Advisor:        Nancy Rodriguez-Bunn
Last modified:  Keshav Patel
Date modified:  09/04/2018

This file scans over the parameter space and computes an MLE between the
solution and a given dataset. This file tests the case where parameters are
set to 0 in the invalid regions defined by the data. The best model and 
parameters are stored as variables at the end of the program.
%}

%% Ranges to calculate PDE with
clc; clear; close all;
Rrange = [0 .9:.005:1.17];
Erange = 0.001%:.005:.5;
Brange = 0:.001:.05;
Arange = 2:.1:4;
tTotal = 200;
dt = 0.5;
S = 1;

meshSize = 250;
[mx, my] = meshgrid(0:meshSize-1);
PDEToDatax = @(x) 1+floor((x+S/2)*(meshSize-1)/S);
omega = sqrt((mx-meshSize/2).^2+(my-meshSize/2).^2)>=10;
nzNeighbor = omega(1:end-2,1:end-2) + omega(2:end-1,1:end-2) + omega(3:end,1:end-2) ... 
    + omega(1:end-2,2:end-1) + omega(3:end,2:end-1) + omega(1:end-2,3:end) ...
    + omega(2:end-1,3:end)+omega(3:end,3:end);
boundary = zeros(meshSize);
boundary(2:end-1,2:end-1) = (nzNeighbor ~= 8 & nzNeighbor ~= 0) & ...
    omega(2:end-1,2:end-1)==0;

%% Cell arrays to contain results of every run
models = cell(length(Rrange),length(Erange),length(Brange),length(Arange));
us = cell(length(Rrange),length(Erange),length(Brange),length(Arange));

%% Start
for Rind = 1:length(Rrange)
    for Eind = 1:length(Erange)
        for Bind = 1:length(Brange)
            for Aind = 1:length(Arange)
                R = Rrange(Rind);
                E = Erange(Eind);
                B = Brange(Bind);
                A = Arange(Aind);
                disp(['R=' num2str(R) '; E=' num2str(E) '; B=' num2str(B) '; A=' num2str(A)])
                
                %% Heterogeneous constants defn
                rho = @(x,y) R*double(omega(PDEToDatax(y), PDEToDatax(x)))+ ...
                             R/2*double(boundary(PDEToDatax(y), PDEToDatax(x)));
                eta = @(x,y) E*double(omega(PDEToDatax(y), PDEToDatax(x)))+ ...
                             E/2*double(boundary(PDEToDatax(y), PDEToDatax(x)));
                alpha = @(x,y) A;
                beta = @(x,y) B*double(omega(PDEToDatax(y), PDEToDatax(x)))+ ...
                              B/2*double(boundary(PDEToDatax(y), PDEToDatax(x)));
                
                timeSteps = tTotal/dt+1; %Defines timeSteps based off tTotal and dt above
                
                %% Run and time simulation
                
                tic
                if ~exist('oldResult','var')
                    fprintf('Using default initial condition\n')
                    [result,model,d] = continuumSimulation(rho, eta, beta, dt, ...
                        timeSteps, S, alpha, []);
                else
                    fprintf('Using old result as initial condition\n')
                    [result,model,d] = continuumSimulation(rho, eta, beta, dt, ...
                        timeSteps, S, alpha, oldResult);
                end
                toc
                
                %% If running on the cluster
                %save('gent1Step.mat','result','model','d','timeSteps','R','A','S','E','gamma','-v7.3')
                
                
                %% Getting results
                models{Rind,Eind,Bind,Aind} = model;
                us{Rind,Eind,Bind,Aind} = result.NodalSolution;
                if size(sum(result.NodalSolution)) == [1,2,2]
                    models{Rind,Eind,Bind,Aind} = NaN;
                    us{Rind,Eind,Bind,Aind} = NaN;
                end
                %pdeplot(model,'XYData',u(:,1,end),'ColorMap','jet','ColorBar','on');
                %drawnow
            end
        end
    end
end


%% Save everything
clearvars -except models us Rrange Erange Brange Arange tTotal dt
today = datestr(now, 'mm-dd');
save([today '_PDEResults_SpatHeteroConst_squareGeom_removedCircleCenter_noFluxBoundary.mat'])

