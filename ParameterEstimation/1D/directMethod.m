%% Linear Sensitivity of Wealth and Amenity Steady-State Solutions
% R = rho; E = eta; B = gamma; A = alpha
clc; clear; close all;

%% Define Functions

% The ODE system given by Equations 5-8 in the 2020 Summer Updates document
dydx = @(y,R,E,B,A) [y(2)+2*y(1)*y(4)/y(3); -A/E*(R*(1-y(1))*y(1)-y(1)); y(4); -A/E*(-y(3)+y(1)+B)];

% The Jacobian operator used to compute the sensitivity DEs
%   See Equations 13-16 in the 2020 Summer Updates document
dfdy = @(y,R,E,B,A) [2*y(4)/y(3),           1,  -2*y(1)*y(4)/y(3)^2,  2*y(1)/y(3); ...
                     -A/E*(R*(1-2*y(1))-1), 0,  0,                    0;           ...
                     0,                     0,  0,                    1;           ...
                     -A/E,                  0,  A/E,                  0];

% The derivative of the ODE RHS with respect to the parameters
dfdrho = @(y,R,E,B,A)   [0; -A/E*(1-y(1))*y(1);            0;  0                   ];
dfdeta = @(y,R,E,B,A)   [0; A/E^2*(R*(1-y(1))*y(1)-y(1));  0;  A/E^2*(-y(3)+y(1)+B)];
dfdgamma = @(y,R,E,B,A) [0; 0;                             0;  -A/E                ];
dfdalpha = @(y,R,E,B,A) [0; -1/E*(R*(1-y(1))*y(1)-y(1));   0; -1/E*(-y(3)+y(1)+B)  ];

%Each row of the matrix "const" contains a row vector of parameters for
%which we would like to run; The first vector contains the smallest values
%of each parameter we would like to search; in each subsequent row, exactly
%one parameter is adjusted to be the largest value. After all 5 simulations
%are run, we can compute the linear sensitivity for the ith parameter by
%using the first simulation and the i+1 simulation.
const = [1.001 1e-4 1e-8 10];

%% Run Simulations
results = cell(4,1);
for meshPoints = round(10.^(2.5:.5:5))
    disp(['Grid Spacing = ' num2str(1/(meshPoints+1))])
    try
        
        % select the parameters needed for the current simulation
        t_c = .00001;
        R = const(1); E = const(2); B = const(3); A = const(4);
        rho=@(x) R; eta=@(x) E/t_c; gamma=@(x) B; alpha=@(x) A/t_c;
        
        % Time dependent simulation
        [xmesh, result_dt] = continuumSimulation_1D(rho,eta,gamma,1,alpha,meshPoints,[]);
        initialGuess = {xmesh, squeeze(result_dt(end,:,:))};
        
        % Boundary value simulation
        % Note that because the time dependent run is independent of the linear
        % sensitivity variables (the Z's), we can use the same time dependent
        % solution for all four boundary value runs
        disp('Time Dependent Run Finished, On to Steady-State')
        for i = 1:4
            lastwarn('', '');
            paramsCell = {'rho'; 'eta'; 'gamma'; 'alpha'};
            disp(['computing BVP with sensitivity in parameter ' paramsCell{i}])
            % Define the vector of linear sensitivity DEs based on what
            % parameter we would like to calculate the linear sensitivity for
            if i == 1
                dZdx = @(y,R,E,B,A) dfdrho(y,R,E,B,A) + dfdy(y,R,E,B,A)*y(5:end);
            elseif i == 2
                dZdx = @(y,R,E,B,A) dfdeta(y,R,E,B,A) + dfdy(y,R,E,B,A)*y(5:end);
            elseif i == 3
                dZdx = @(y,R,E,B,A) dfdgamma(y,R,E,B,A) + dfdy(y,R,E,B,A)*y(5:end);
            else
                dZdx = @(y,R,E,B,A) dfdalpha(y,R,E,B,A) + dfdy(y,R,E,B,A)*y(5:end);
            end
            
            ode = @(t,y) [dydx(y,R,E,B,A); dZdx(y,R,E,B,A)];
            bc = @bcfun;
            yinit = @(x) initGuess(x,initialGuess);
            solinit = bvpinit(initialGuess{1},yinit);
            options = bvpset('Nmax', 1000000);
            
            % Compute the time-independent solution
            result = bvp5c(ode,bc,solinit,options);
            
            % bvp5c found a fixed point, but it has negative values for y1 
            % or y3, which is not physically relevant
            if nnz(sum(sum(result.y([1 3],:)<0)))
                error('Found negative values in solution')
            end
            
            % collects any warning messages from the Boundary value run
            [warnMsg, warnId] = lastwarn();
            
            % If there was a warning (i.e. "Unable to meet tolerance"),
            % trigger the catch block to restart the simulation on a finer
            % mesh
            if(~isempty(warnId))
                error(warnMsg, warnId);
            else
                results{i} = result;
            end
        end
        
    catch
        %bvp5c ran into numerical issues
        disp('Boundary Value Run Failed! Retrying with finer grid spacing')
        continue;
    end
    break;
end

%% Calculating L
Lrho = const(1,1)*abs(results{1}.y(5,:)./results{1}.y(1,:));
Leta = const(1,2)*abs(results{2}.y(5,:)./results{2}.y(1,:));
Lgamma = const(1,3)*abs(results{3}.y(5,:)./results{3}.y(1,:));
Lalpha = const(1,4)*abs(results{4}.y(5,:)./results{4}.y(1,:));

%% Plotting
figure(1); 
semilogy(results{1}.x,Lrho, 'linewidth', 2); 
hold on; 
semilogy(results{2}.x,Leta, 'linewidth', 2); 
semilogy(results{3}.x,Lgamma, 'linewidth', 2); 
semilogy(results{4}.x,Lalpha, 'linewidth', 2);
set(gca, 'ticklabelinterpreter', 'latex');
set(gca, 'fontsize', 24);
xlabel('$x$', 'interpreter', 'latex')
ylabel('Linear Sensitivity', 'interpreter', 'latex')
legend({'$L_\rho$'; '$L_\eta$';'$L_\gamma$';'$L_\alpha$'}, 'interpreter', 'latex', 'location', 'northeastoutside')
set(gcf, 'Position',  [100, 100, 800, 800])




%% Helper Functions

function res = bcfun(ya,yb)
res = zeros(8,1);
res(1) = ya(2);
res(2) = yb(2);
res(3) = ya(4);
res(4) = yb(4);
res(5) = ya(6);
res(6) = yb(6);
res(7) = ya(8);
res(8) = yb(8);
end

function solinit = initGuess(x,initialGuess)

xmesh = initialGuess{1}; result_dt = initialGuess{2};
W = griddedInterpolant(xmesh,result_dt(:,1));
A = griddedInterpolant(xmesh,result_dt(:,2));
Wprime = griddedInterpolant(xmesh, (W(xmesh([2:end end]))-W(xmesh([1 1:end-1])))./(xmesh([2:end end])-xmesh([1 1:end-1])));
Aprime = griddedInterpolant(xmesh, (A(xmesh([2:end end]))-A(xmesh([1 1:end-1])))./(xmesh([2:end end])-xmesh([1 1:end-1])));

solinit = [W(x);Wprime(x)-2*W(x)*Aprime(x)/A(x);A(x);Aprime(x); 0; 0; 0; 0];

end