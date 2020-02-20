%% PDE Toolbox Implementation for PDE System
%{
Author:         Avishai Halev
Advisor:        Nancy Rodriguez-Bunn
Last modified:  Keshav Patel
Date modified:  09/04/2018

This file implements MATLAB's PDE Toolbox package as a solver for the PDE
system given in the group. The scripts within this package implement this
file for each set of variables.
%}

function [result] = continuumSimulation_1D(R,E,B, dt,...
    timeSteps,S,A,initialConditionCell)
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

%% create mesh
% determine size of space

xmesh = linspace(-S/2,S/2,200);

%% BCs
% default is zero neumann, leave unspecified.

bc = @bcfun;

%% ICs
ic = @(x) icfun(x,initialConditionCell,xmesh);

%% PDE to Solve
pde = @(x,t,u,dudx) pdefun(x,t,u,dudx,R,E,B,A);

%% Solve
result = pdepe(0,pde,ic,bc,xmesh,linspace(0,timeSteps*dt,timeSteps));


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


