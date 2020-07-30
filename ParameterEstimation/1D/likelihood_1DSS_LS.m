%% Energy Functional Using Steady-State Solution and Least Squares Error in Data

function ls_mag = likelihood_1DSS_LS(X, data)

%numerical parameters
S = 1;
timeDilate = 0.001;

fprintf('***current position rho_left = %.7e; rho_right = %.7e; eta_left = %.7e; eta_right = %.7e; gamma_left = %.7e; gamma_right = %.7e; alpha_left %.7e; alpha_right %.7e\n', ...
    X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8));

tic;

%This for loop functions to refine the spatial grid if a solution is not found
for meshPoints = round(10.^(2.5:.5:5))
    disp(['Grid Spacing = ' num2str(1/(meshPoints+1))])
    try
        %piecewise constant parameters using time-dilating change of variables
        rho=@(x) X(1)+(X(2)-X(1))*double(x>0); eta=@(x) X(3)/timeDilate+(X(4)-X(3))/timeDilate*double(x>0);
        gamma=@(x) X(5)+(X(6)-X(5))*double(x>0); alpha=@(x) X(7)/timeDilate+(X(8)-X(7))/timeDilate*double(x>0);
        
        %run time-dependent simulation
        [xmesh, result_dt] = continuumSimulation_1D(rho,eta,gamma,S,alpha,meshPoints+1,[]);
        disp('Time Dependent Run Finished, On to Steady-State')
        
        %redefine eta and alpha without time dilation
        eta=@(x) X(3)+(X(4)-X(3))*double(x>0); alpha=@(x) X(7)+(X(8)-X(7))*double(x>0);
        result = steadyStateSimulation_1D(rho,eta,gamma,S,alpha,{xmesh, squeeze(result_dt(end,:,:))});
        
        
        if nnz(sum(sum(result.y([1 3],:)<0)))
            %bvp5c found a fixed point, but it is not physically relevant
            disp('Found negative values in solution. Retrying with finer grid spacing');
            continue;
        else
            break;
        end
    catch
        %bvp5c ran into numerical issues
        disp('Boundary Value Run Failed! Retrying with finer grid spacing')
        continue;
    end
end
disp(['Took ' num2str(toc) ' seconds to find the solution'])

%Convert Grid functions from bvp5c to piecewise continuous functions
W = griddedInterpolant(result.x,result.y(1,:));
A = griddedInterpolant(result.x,result.y(3,:));

%calculate log likelihood
wealthLS = norm(W(data(:,1)) - data(:,2),2);
amenLS = norm(A(data(:,1)) - data(:,3),2);

ls_mag = sqrt(wealthLS^2+amenLS^2);
fprintf('***Least Squares Error: %.7e\n\n', ls_mag);

end