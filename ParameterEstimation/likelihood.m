function lf_mag = likelihood(rho, eta, gamma, alpha, timeSteps, Wmap, Amap)
S = 1;

disp(['rho = ' num2str(rho) '; eta = ' num2str(eta) ...
    '; gamma = ' num2str(gamma) '; alpha = ' num2str(alpha)])

if ~exist('initialConditionCell','var')
    fprintf('Using default initial condition\n')
    try
        [result,model,d,ICVec] = continuumSimulation(rho, eta, gamma, ...
            0.5, timeSteps, 1, alpha, []);
    catch
    end
else
    fprintf('Using old result as initial condition\n')
    [result,model,d,ICVec] = continuumSteadyState(rho, eta, gamma, ...
        0.5, timeSteps, 1, alpha, initialConditionCell);
end

dx = 1/100;
[meshx, meshy] = meshgrid(linspace(-S/2+dx,S/2-dx,100));

%Interpolate the PDE solution to match the probability map
%that contains the wealth and amenity data. Later we will
%use this to normalize the solution as well
Wmesh = reshape(interpolateSolution(result,[meshx(:) meshy(:)]', 1, timeSteps),size(Wmap));
Wmesh(isnan(Wmesh)) = 0;
Amesh = reshape(interpolateSolution(result,[meshx(:) meshy(:)]', 2, timeSteps),size(Amap));
Amesh(isnan(Amesh)) = 0;
% 
% %we also need the data from the probability map
% %Wvalues is the number of times a value has been randomly
% %sampled
% %Wintrp is the set of interpolated solutions at the places
% %the map was sampled
% mapPosX = linspace(-S/2,S/2,size(Wmap,2));
% mapPosY = linspace(S/2,-S/2,size(Wmap,1));
% [nonzeroY, nonzeroX] = find(Wmap~=0);
% Wvalues = Wmap(nonzeroY+100*(nonzeroX-1));
% queryX = mapPosX(nonzeroX);
% queryY = mapPosY(nonzeroY);
% Wintrp = interpolateSolution(result,[queryX; queryY], 1, timeSteps);
% Wintrp(isnan(Wintrp)) = 0;
% 
% %Do the same thing for the amenity map
% mapPosX = linspace(-S/2,S/2,size(Amap,2));
% mapPosY = linspace(S/2,-S/2,size(Amap,1));
% [nonzeroY, nonzeroX] = find(Amap~=0);
% Avalues = Amap(nonzeroY+100*(nonzeroX-1));
% queryX = mapPosX(nonzeroX);
% queryY = mapPosY(nonzeroY);
% Aintrp = interpolateSolution(result,[queryX; queryY], 2, timeSteps);
% Aintrp(isnan(Aintrp)) = 0;

%normalize PDE solution by dividing by integral
Wintegral = sum(sum(dx^2*Wmesh(:,1,end)));
Aintegral = sum(sum(dx^2*Amesh(:,1,end)));
WintrpTimej = Wmesh./Wintegral;
AintrpTimej = Amesh./Aintegral;

%calculate log likelihood
wealthLF = sum(sum(log(WintrpTimej.^Wmap(end:-1:1,:))));
amenLF = sum(sum(log(AintrpTimej.^Amap(end:-1:1,:))));
format long
lf_mag = -sqrt(wealthLF^2+amenLF^2);
stdWmesh = sqrt(sum(sum((Wmesh - mean(mean(Wmesh))).^2))/numel(Wmesh));

disp(['likelihood = ' num2str(lf_mag)])
disp(['std of solution = ' num2str(stdWmesh)])


end