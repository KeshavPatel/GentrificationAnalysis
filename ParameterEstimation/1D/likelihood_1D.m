function lf_mag = likelihood(rho, eta, gamma, alpha, timeSteps, Wset, Aset)

format long
[rho eta gamma alpha]

S = 1;
tol = 1e-13;

err = 2*tol;
runs = 1;
warning('')

while err > tol && runs < 300
    if ~exist('ic','var')
        ic = [];
        fprintf('Using default initial condition\n')
        try
            result = continuumSimulation_1D(rho, eta, gamma, ...
                timeSteps/3, 3, 1, alpha, []);
        catch
        end
    else
        try
            fprintf('Using old result as initial condition\n')
            result = continuumSimulation_1D(rho, eta, gamma, ...
                timeSteps/3, 3, 1, alpha, ic);
        catch
        end
    end
    [warnMsg, warnId] = lastwarn;
    %if ~isempty(warnMsg)
    %    break;
    %end
    err = sum(sum((result(end,:,:)-result(1,:,:)).^2));
    ic = result;
    runs = runs+1;
end

%normalize PDE solution by dividing by integral
Wintegral = 1/200*sum(result(end,:,1));
Aintegral = 1/200*sum(result(end,:,2));
W = griddedInterpolant(linspace(-S/2,S/2,200),result(end,:,1));
A = griddedInterpolant(linspace(-S/2,S/2,200),result(end,:,2));

%calculate log likelihood

%wealthLF = sum(log(W(Wset)));
%amenLF = sum(log(A(Aset)));
wealthLF = prod(W(Wset)./Wintegral);
amenLF = prod(A(Aset)./Aintegral);

lf_mag = -sqrt(wealthLF^2+amenLF^2);% + 100*((Wintegral-1/200*sum(W(Wset))) + (Aintegral-1/200*sum(A(Aset))));


end