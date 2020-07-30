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

function result = steadyStateSimulation_1D(R,E,B,S,A,initialGuess)

% BCs
bc = @(ya,yb) bcfun(ya,yb);

% ICs
yinit = @(x) initGuess(x,initialGuess);
solinit = bvpinit(initialGuess{1},yinit);

% PDE to Solve
ode = @(x,y) odefun(x,y,R(x),E(x),B(x),A(x));

% Solve
options = bvpset('Nmax', 1000000);

result = bvp5c(ode,bc,solinit,options);

end


%% Helper Functions
function dydt = odefun(x,y,R,E,B,A)

dydt = zeros(4,1);

dydt(1) = y(2) + 2*y(1)*y(4)/y(3);
dydt(2) = -A/E*(R*(1-y(1))*y(1)-y(1));
dydt(3) = y(4);
dydt(4) = -A/E*(-y(3)+y(1)+B);

end

function res = bcfun(ya,yb)
res = zeros(4,1);

res(1) = ya(2);
res(2) = yb(2);
res(3) = ya(4);
res(4) = yb(4);

end

function solinit = initGuess(x,initialGuess)

if ~isempty(initialGuess)
    xmesh = initialGuess{1}; result_dt = initialGuess{2};
    W = griddedInterpolant(xmesh,result_dt(:,1));
    A = griddedInterpolant(xmesh,result_dt(:,2));
    Wprime = griddedInterpolant(xmesh, (W(xmesh([2:end end]))-W(xmesh([1 1:end-1])))./(xmesh([2:end end])-xmesh([1 1:end-1])));
    Aprime = griddedInterpolant(xmesh, (A(xmesh([2:end end]))-A(xmesh([1 1:end-1])))./(xmesh([2:end end])-xmesh([1 1:end-1])));
    
    solinit = [W(x);Wprime(x)-2*W(x)*Aprime(x)/A(x);A(x);Aprime(x)];
    
else
    load('temp.mat')
    W = griddedInterpolant(result.x,result.y(1,:));
    A = griddedInterpolant(result.x,result.y(3,:));
    Wprime = griddedInterpolant(result.x,result.y(2,:));
    Aprime = griddedInterpolant(result.x,result.y(4,:));
    solinit = [W(x);Wprime(x);A(x);Aprime(x)];
end
end
