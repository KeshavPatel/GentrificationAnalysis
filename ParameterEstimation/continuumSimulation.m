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

function [ result ,model,dimensions, ICVec] = continuumSimulation(R,E,B, dt,...
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


model = createpde(2);


%% create mesh
% determine size of space

if size(S)==[1 1]
    SS = [S S];
elseif size(S)== [1 2] || size(S)==[2 1]
    SS = S;
else
    error('Size dimension mismatch')
end

s1 = SS(1)/2;
s2 = SS(2)/2;
s3 = s1*4/5;
s4 = s2*4/5;

%For a rectangle, the first row contains 3, and the second row contains 4.
%    The next four rows contain the x-coordinates of the starting points
%    of the edges, and the four rows after that contain the y-coordinates
%    of the starting points of the edges.
R1 = [3,4,-s1 s1 s1 -s1 -s2 -s2 s2 s2]';
%R2 = [3,4,-s1+.1 -.1 -.1 -s1+.1 .1 .1 s2-.1 s2-.1]';
%R1 = [1, 0, 0, s1]';
%R2 = [1, 0, 0, s1/5]';
gm = [R1];% R2];
sf = 'R1';%-R2';
ns = char('R1');
%ns2 = char('R2');
ns = [ns'];% ns2'];

%Creates the domain on which the PDE is run. Here, it is a rectangle centered at the origin
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);

dimensions = [-s1 s1 ; -s2 s2];

%% to see geometry
%pdegplot(model,'EdgeLabels','on') % 'EdgeLabels' for 2-D

%% BCs
% default is zero neumann, leave unspecified.



%% generate mesh
hmax = 0.1*max(SS); % max element size.
msh = generateMesh(model,'Hmax',hmax);
%figure;
%pdeplot(model);
%axis equal
%title 'Plate With Triangular Element Mesh'
%xlabel 'X-coordinate, meters'
%ylabel 'Y-coordinate, meters'
p = msh.Nodes;


% Set solver options. If you're running into trouble w/R<1, try
% decreasing tolerances.
model.SolverOptions.RelativeTolerance = 1.0e-3;
model.SolverOptions.AbsoluteTolerance = 1.0e-5;
%model.SolverOptions.MinStep = 1e-7;
%model.SolverOptions.MaxIterations = 5000;
%model.SolverOptions.ReportStatistics = 'on';


%% Coefficients
% ccoeff and fcoeff are helper functions defined below.

a = [0 ; 0];

% RHS
cFn = @(region,state) ccoeff(region,state,R,E);
fFn = @(region,state) fcoeff(region,state,R,A,B);

D = [1;1];     % LHS
specifyCoefficients(model,'m',0,'f',fFn,'a',a,'d',D,'c',cFn);

%% ICs
if isempty(initialConditionCell)
    ICVec = @IC;
    setInitialConditions(model,ICVec);
else
    ICVec = @(locations) ICData(locations, initialConditionCell);
    setInitialConditions(model,ICVec);
end

%% Solve
result = solvepde(model,0:dt:timeSteps);


%% End of main function
end

%% Helper functions

function fmatrix = fcoeff(region,state,rho,alpha,beta)

% Everything that does not involve derivatives goes here.
% u(1,:) is wealth, u(2,:) amenities


n1 = 2;
nr = numel(region.x);
fmatrix = zeros(n1,nr);

fmatrix(1,:) = alpha*(rho.*(1-state.u(1,:)).*(state.u(1,:))-state.u(1,:));
fmatrix(2,:) = alpha*(state.u(1,:)-state.u(2,:)+beta);


end



function u0 = ICData(locations, result)

u0 = interpolateSolution(result,[locations.x;locations.y],1:2)';

%u0 = 3e-3*(rand(size(u0))-.5)+u0;

end



function u0 = IC(locations)

% Initial condition


u0 = zeros(2,length(locations.x));

%n defines initial condition.

n = 1;

%   n = 1 is little lumpy sinusoids
%   n = 2 is Gaussian
%   n = 3 is random initial data
%   default is little lumpy sinusoids

%num is the number of little lumpy sinusoids in each direction. Only set up
%for square grid for now.
num = 16;

% Leave this

l = max(locations.x)-min(locations.x);
k = pi*num/l;


% c and factor scale the initial conditions.
c=0.01;
factor = 1/c*0.2;

switch n
    case 1
        u0(1,:) = c*sin(k*locations.x).*sin(k*locations.y)+factor*c;
        u0(2,:) =  -c.*sin(k*locations.x).*sin(k*locations.y)+factor*c;
        
    case 2
        u0(1,:) = c*exp(-1*((locations.x).^2+(locations.y).^2));
        u0(2,:) = c*exp(-1*((locations.x).^2+(locations.y).^2));
        
    case 3
        sSize = size(u0);
        u0 = c+factor*c*rand(sSize);
        
    otherwise
        u0(1,:) = c*sin(k*locations.x).*sin(k*locations.y)+factor*c;
        u0(2,:) =  -c.*sin(k*locations.x).*sin(k*locations.y)+factor*c;
        
end

end



function cmatrix = ccoeff(region,state,~,eta)

% Diffusive terms


n1 = 16;
nr = numel(region.x);

cmatrix = zeros(n1,nr);

cmatrix(1,:) = eta;
cmatrix(4,:) = eta;

cmatrix(9,:) = -2*eta*state.u(1,:)./state.u(2,:);
cmatrix(12,:) = -2*eta*state.u(1,:)./state.u(2,:);

cmatrix(13,:) = eta;
cmatrix(16,:) = eta;

end