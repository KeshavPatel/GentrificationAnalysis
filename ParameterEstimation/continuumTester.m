close all
clearvars -except oldResult pResult


%% Define parameters
dt = 0.5;
tTotal = 100;      % total time
R = 1.15;
E = 0.01;
S = 1;
A = 4;
B = 0.0;

meshSize = 250;
[mx, my] = meshgrid(1:250);
PDEToDatax = @(x) 1+floor((x+S/2)*(meshSize-1)/S);
omega = (mx-meshSize/2).^2+(my-meshSize/2).^2>=25;
boundaryDim1 = 1/10;
boundaryDim2 = 1/10;
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
    [result,model,d] = continuumSimulation(rho, eta, beta, dt, timeSteps, S, alpha, []);
else
    fprintf('Using old result as initial condition\n')
    [result,model,d] = continuumSimulation(rho, eta, beta, dt, timeSteps, S, alpha, oldResult);
end
toc

%% If running on the cluster
%save('gent1Step.mat','result','model','d','timeSteps','R','A','S','E','gamma','-v7.3')


%% Getting results
u = result.NodalSolution;
j = 1;


%% This is for debugging or trying to figure out WTF is happening
for k=floor(0.2*timeSteps):(0.5*timeSteps)
    fprintf('Absolute difference %4.2e relative difference %4.2f mean %4.2e\n',max(u(:,1,k))-min(u(:,1,k)),max(u(:,1,k))/min(u(:,1,k)),mean(u(:,1,k)));
    %disp(max(u(:,1,k))-min(u(:,1,k)));
end

%% Plotting

%Plots every step timeSteps
step = 10;

for i=1:step:timeSteps
    pdeplot(model,'XYData',u(:,1,i),'ColorMap','jet','ColorBar','off');
    colorbar('westoutside')
    title('')
    xlabel 'X-coordinate, meters'
    ylabel 'Y-coordinate, meters'
    ax = gca;                   % fixing color limits.
    ax.CLim = [0 0.2];
    xlim(d(1,:))
    ylim(d(2,:))
    zlim([0 0.1])
    Mv(j)=getframe;
    j = j+1;
    pause(.1)
end

%% Saving a movie. This can take a long time


v = VideoWriter('movie1');
v.FrameRate = 5;
v.Quality = 100;
open(v)
writeVideo(v,Mv)
close(v)


