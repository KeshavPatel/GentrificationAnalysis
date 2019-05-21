%% Basic MLE Maximizer - Changing Geometry
%{
Author:         Keshav Patel
Advisor:        Nancy Rodriguez-Bunn
Last modified:  Keshav Patel
Date modified:  09/04/2018

This file scans over the parameter space and computes an MLE between the
solution and a given dataset. This file tests the case where the geometry
of the space is determined by invalid regions in the data. The best model 
and parameters are stored as variables at the end of the program.
%}

%% MLE Initializations
wealthMaxLogLikeFunc = -Inf;
amenMaxLogLikeFunc = -Inf;
bestModel = 0;
bestu = 0;
bestt = 0;
bestParam = [];
WbestLSE = 0;
AbestLSE = 0;

%% Start

for R = Rrange
    for E = Erange
        for B = Brange
            for A = Arange
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
                    [result,model,d] = continuumSimulationGeometry(rho, eta, beta, dt, ...
                        timeSteps, S, alpha, []);
                else
                    fprintf('Using old result as initial condition\n')
                    [result,model,d] = continuumSimulationGeometry(rho, eta, beta, dt, ...
                        timeSteps, S, alpha, oldResult);
                end
                toc
                
                %% If running on the cluster
                %save('gent1Step.mat','result','model','d','timeSteps','R','A','S','E','gamma','-v7.3')
                
                
                %% Getting results
                u = result.NodalSolution;
                %pdeplot(model,'XYData',u(:,1,end),'ColorMap','jet','ColorBar','on');
                %drawnow
                
                %% Compare PDE and data
                dx = S/250;
                [meshx, meshy] = meshgrid(linspace(-S/2+dx,S/2-dx,250));
                
                %Interpolate the PDE solution to match the probability map
                %that contains the wealth and amenity data. Later we will
                %use this to normalize the solution as well
                Wmesh = interpolateSolution(result,[meshx(:) meshy(:)]', 1, 1:timeSteps);
                Wmesh(isnan(Wmesh)) = 0;
                Amesh = interpolateSolution(result,[meshx(:) meshy(:)]', 2, 1:timeSteps);
                Amesh(isnan(Amesh)) = 0;
                
                %we also need the data from the probability map
                %Wvalues is the number of times a value has been randomly
                %sampled
                %Wintrp is the set of interpolated solutions at the places
                %the map was sampled
                mapPosX = linspace(-S/2,S/2,size(Wmap,2));
                mapPosY = linspace(S/2,-S/2,size(Wmap,1));
                [nonzeroY, nonzeroX] = find(Wmap~=0);
                queryX = mapPosX(nonzeroX);
                queryY = mapPosY(nonzeroY);
                Wintrp = interpolateSolution(result,[queryX; queryY], 1, 1:timeSteps);
                %Wintrp = Wintrp(~isnan(Wintrp));
                
                %Do the same thing for the amenity map
                mapPosX = linspace(-S/2,S/2,size(Amap,2));
                mapPosY = linspace(S/2,-S/2,size(Amap,1));
                [nonzeroY, nonzeroX] = find(Amap~=0);
                queryX = mapPosX(nonzeroX);
                queryY = mapPosY(nonzeroY);
                Aintrp = interpolateSolution(result,[queryX; queryY], 2, 1:timeSteps);
                %Aintrp = Aintrp(~isnan(Aintrp));
                
                %check the solution at each time point and compare to the
                %data
                for j=floor(10/dt):timeSteps
                    %normalize PDE solution by dividing by integral
                    Wintegral = sum(sum(dx^2*reshape(Wmesh(:,1,j),[250 250])));
                    Aintegral = sum(sum(dx^2*reshape(Amesh(:,1,j),[250 250])));
                    WintrpTimej = Wintrp(:,1,j)./Wintegral;
                    AintrpTimej = Aintrp(:,1,j)./Aintegral;
                    
                    %Use this to compare maps by Least Squares
                    %(this is just an alternative metric for determining
                    % how close the solution and data are; this does nothing
                    % for the algorithm itself)
                    Wuse = ~isnan(WintrpTimej);
                    Ause = ~isnan(AintrpTimej);
                    WmeshTimej = reshape(Wmesh(:,1,j)./Wintegral,[250,250]);
                    AmeshTimej = reshape(Amesh(:,1,j)./Aintegral,[250,250]);
                    
                    %calculate log likelihood
                    wealthLLF = sum(log(WintrpTimej(Wuse)));
                    amenLLF = sum(log(AintrpTimej(Ause)));
                    
                    % Check if these values are the new minimum
                    if wealthLLF > wealthMaxLogLikeFunc && amenLLF > amenMaxLogLikeFunc
                        wealthMaxLogLikeFunc = wealthLLF;
                        amenMaxLogLikeFunc = amenLLF;
                        bestModel = model;
                        bestu = u;
                        bestt = j;
                        bestParam = [R E B A];
                        WbestLSE = sqrt(sum(sum((WmeshTimej-map).^2)));
                        AbestLSE = sqrt(sum(sum((AmeshTimej-map).^2)));
                    end
                end
            end
        end
    end
end
%% This is for debugging or trying to figure out WTF is happening
%for k=floor(0.2*timeSteps):(0.5*timeSteps)
%    fprintf('Absolute difference %4.2e relative difference %4.2f mean %4.2e\n', ...
%        max(u(:,1,k))-min(u(:,1,k)),max(u(:,1,k))/min(u(:,1,k)),mean(u(:,1,k)));
%disp(max(u(:,1,k))-min(u(:,1,k)));
%end

