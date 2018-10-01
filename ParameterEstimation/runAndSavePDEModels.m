%% Basic MLE Maximizer
%{
Author:         Keshav Patel
Advisor:        Nancy Rodriguez-Bunn
Last modified:  Keshav Patel
Date modified:  09/04/2018

This file scans over the parameter space and computes an MLE between the
solution and a given dataset. The best model and parameters are stored as
variables at the end of the program.
%}

%% Ranges to calculate PDE with
Rrange = [0 .9:.005:1.17];
Erange = 0:.005:.5;
Brange = 0:.001:.05;
Arange = .5:.1:4;
tTotal = 200;
dt = 0.5;
S = 1;

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
                rho = @(x,y) R;
                eta = @(x,y) E;
                alpha = @(x,y) A;
                beta = @(x,y) B;
                
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
                
            end
        end
    end
end

%% Save everything
clearvars -except models us Rrange Erange Brange Arange tTotal dt
today = datestr(now, 'mm-dd');
save([today '_PDEResults_squareGeom_noFluxBoundary.mat'])

