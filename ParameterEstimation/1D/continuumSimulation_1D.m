%% PDE Toolbox Implementation for PDE System
%{
Author:         Avishai Halev
Advisor:        Nancy Rodriguez-Bunn
Last modified:  Keshav Patel
Date modified:  07/30/2020

This file implements MATLAB's PDE Toolbox package as a solver for the PDE
system given in the group. The scripts within this package implement this
file for each set of variables.
%}

function [xmesh, result] = continuumSimulation_1D(R,E,B,S,A,meshPoints,initialConditionCell)
%% CONTINUUMSIMULATION executes one run of the PDE system through time (timeSteps-1)*dt
%      Parameters:
%      R: rho
%      E: eta
%      B: gamma
%      dt: time increment, 0.5 is a good place to start
%      timeSteps: self-explanatory. See figures in write-up for ideas on
%       what dynamics should look like at various time scales
%      S: size of domain. Can be specified in two ways:
%       1) single number. In this case, domain is an SxS square grid,
%           centered at the origin.
%       2) 2x1 array. In this case, domain is an S(1)xS(2) rectangular
%           grid, centered at the origin.
%      A: Numerical parameter. Can be used to scale domain to mimic larger
%       grids. (See Murray Chapter 2). Generally, use 1.
%      initialCondition: IC, specified as a PDERESULT or as []. If empty,
%       IC is specified in IC local function (see below).

%This for loop functions to decrease the time step if a solution is not found
for i = 0:2:12 
    disp(['Timestep = ' num2str(1/10^(i+1))])
    warning('')
    
    %% create mesh
    % determine size of space
    xmesh = linspace(-S/2,S/2,meshPoints);
    xmesh = sort(unique([xmesh linspace(-1/(meshPoints-1),1/(meshPoints-1),51)]));
    
    %% BCs
    bc = @bcfun;
    
    %% ICs
    ic = @(x) icfun(x,initialConditionCell,xmesh);
    
    %% PDE to Solve
    pde = @(x,t,u,dudx) pdefun(x,t,u,dudx,R(x),E(x),B(x),A(x));
    
    %% Solve
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    result = pdepe(0,pde,ic,bc,xmesh,linspace(0,1,10^(i+1)),opts);
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
        disp('Time Dependent Run Failed! Retrying with finer timestep')
    else
        break;
    end
end



%% End of main function
end

%% Helper functions

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
pl = [0; 0];
ql = [1; 1];
pr = [0; 0];
qr = [1; 1];
end

function u0 = icfun(x,icCell, xmesh)
u0 = zeros(2,length(x));

%n defines initial condition.
%   n = 1 is little lumpy sinusoids
%   n = 2 is Gaussian
%   n = 3 is random initial data
%   default is little lumpy sinusoids
n = 1;

%num is the number of little lumpy sinusoids in each direction. Only set up
%for square grid for now.
num = 16;

% Leave this
k = pi*num;

% c and factor scale the initial conditions.
c=0.01;
factor = 1/c*0.2;

switch n
    case 2
        u0(1,:) = c*exp(-1*x.^2);
        u0(2,:) = c*exp(-1*x.^2);
        
    case 3
        sSize = size(u0);
        u0 = c+factor*c*rand(sSize);
        
    otherwise
        u0(1,:) = c*sin(k*x).^2+factor*c;
        u0(2,:) =  -c.*sin(k*x).^2+factor*c;
end

if ~isempty(icCell)
    F1 = griddedInterpolant(xmesh,icCell(end,:,1));
    F2 = griddedInterpolant(xmesh,icCell(end,:,2));
    u0(1,:) = F1(x);
    u0(2,:) = F2(x);
end
end

function [c,f,s] = pdefun(x,t,u,dudx,R,E,B,A)
c = [1; 1];
f = E*[dudx(1)-2*u(1).*dudx(2)./u(2); dudx(2)];
s = A*[R*(1-u(1)).*u(1)-u(1); -u(2)+u(1)+B];
end


