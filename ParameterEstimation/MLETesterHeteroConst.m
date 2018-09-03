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
                u = result.NodalSolution;
                %pdeplot(model,'XYData',u(:,1,end),'ColorMap','jet','ColorBar','on');
                %drawnow
                
                %% Compare PDE and data
                dx = S/250;
                [meshx, meshy] = meshgrid(linspace(-S/2+dx,S/2-dx,250));
                Wmesh = interpolateSolution(result,[meshx(:) meshy(:)]', 1, 1:timeSteps);
                Wmesh(isnan(Wmesh)) = 0;
                %Wintegral = sum(Wmesh*dx^2);
                Amesh = interpolateSolution(result,[meshx(:) meshy(:)]', 2, 1:timeSteps);
                Amesh(isnan(Amesh)) = 0;
                %Aintegral = sum(Amesh*dx^2);
                
                mapPosX = linspace(-S/2,S/2,size(Wmap,2));
                mapPosY = linspace(S/2,-S/2,size(Wmap,1));
                [nonzeroY, nonzeroX] = find(Wmap~=0);
                Wvalues = Wmap(nonzeroY+250*(nonzeroX-1));
                queryX = mapPosX(nonzeroX);
                queryY = mapPosY(nonzeroY);
                Wintrp = interpolateSolution(result,[queryX; queryY], 1, 1:timeSteps);
                Wintrp(isnan(Wintrp)) = 0;
                
                mapPosX = linspace(-S/2,S/2,size(Amap,2));
                mapPosY = linspace(S/2,-S/2,size(Amap,1));
                [nonzeroY, nonzeroX] = find(Amap~=0);
                Avalues = Amap(nonzeroY+250*(nonzeroX-1));
                queryX = mapPosX(nonzeroX);
                queryY = mapPosY(nonzeroY);
                Aintrp = interpolateSolution(result,[queryX; queryY], 2, 1:timeSteps);
                Aintrp(isnan(Aintrp)) = 0;
                
                for j=floor(10/dt):timeSteps
                    Wintegral = sum(sum(dx^2*reshape(Wmesh(:,1,j),[250 250])));
                    Aintegral = sum(sum(dx^2*reshape(Amesh(:,1,j),[250 250])));
                    WintrpTimej = Wintrp(:,1,j)./Wintegral;
                    AintrpTimej = Aintrp(:,1,j)./Aintegral;
                    WmeshTimej = reshape(Wmesh(:,1,j)./Wintegral,[250,250]);
                    AmeshTimej = reshape(Amesh(:,1,j)./Aintegral,[250,250]);
                    wealthLLF = sum(log(WintrpTimej.^Wvalues));
                    amenLLF = sum(log(AintrpTimej.^Avalues));
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

