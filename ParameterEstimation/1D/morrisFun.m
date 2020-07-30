function out = morrisFun(R, E, B, A, test)

fprintf('***current position %.7e %.7e %.7e %.7e\n', R, E, B, A);

%Numerical parameters
S = 1;
timeDilate = 0.0000001;

fprintf('***current position rho = %.7e; eta = %.7e; gamma = %.7e; alpha = %.7e\n', ...
    R, E, B, A);

tic;
%% Find Steady-State Solution
%This for loop functions to refine the spatial grid if a solution is not found
for meshPoints = round(10.^(2.5:.5:5))
    disp(['Grid Spacing = ' num2str(1/(meshPoints+1))])
    try
        %piecewise constant parameters using time-dilating change of variables
        rho=@(x) R;   eta=@(x) E/timeDilate;
        gamma=@(x) B; alpha=@(x) A/timeDilate;
        
        %run time-dependent simulation
        [xmesh, result_dt] = continuumSimulation_1D(rho,eta,gamma,S,alpha,meshPoints+1,[]);
        disp('Time Dependent Run Finished, On to Steady-State')
        
        %redefine eta and alpha without time dilation
        eta=@(x) E; alpha=@(x) A;
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

%% Compute Metric
if test == 2
    out = fwhm(result,S);
elseif test == 3
    out = spikeDist(result,S);
elseif test == 4
    out = spikeMin(result,S);
else
    out = spikeMax(result,S);
end

end

%% Metrics
function out = spikeMin(result,S)
whichEqn = 1;
uintegral = trapz(result.x,result.y(whichEqn,:));
u = griddedInterpolant(result.x,result.y(whichEqn,:));
out = min(u(linspace(-S/2,S/2,5000))./uintegral);

fprintf('***Amplitude: %.7e\n\n', out);

end

function out = spikeMax(result,S)
whichEqn = 1;
uintegral = trapz(result.x,result.y(whichEqn,:));
u = griddedInterpolant(result.x,result.y(whichEqn,:));
out = max(u(linspace(-S/2,S/2,5000))./uintegral);

fprintf('***Amplitude: %.7e\n\n', out);

end

function out = spikeDist(result,S)
Wintegral = trapz(result.x,result.y(1,:));
W = griddedInterpolant(result.x,result.y(1,:));
Wprime = griddedInterpolant(result.x,result.y(2,:));
s = 1;
spike_list = [];

for i = linspace(-S/2,S/2,5000)
    if sign(Wprime(i)) == -1*s
        s = sign(Wprime(i));
        if s == -1
            spike_list = [spike_list i];
        end
    end
end
out = norm(spike_list(2:end)-spike_list(1:end-1));
fprintf('***Spike Distance: %.7e\n\n', out);

end

function out = fwhm(result,S) %%TODO
whichEqn = 1;
u = griddedInterpolant(result.x,result.y(whichEqn,:));
uprime = griddedInterpolant(result.x,result.y(whichEqn+1,:));

s = sign(uprime(-S/2));
extremePrev = -inf;
fwhm_list = [];
fwhm_temp = 0;
for i = linspace(-S/2,S/2,5000)
    if sign(uprime(i)) == -1*s
        s = sign(uprime(i));
        if extremePrev == -inf
            extremePrev = i;
            continue;
        end
        if s == 1
            hm = (u(i)+u(extremePrev))/2;
            [M,k] = min(abs(u(linspace(extremePrev,i,50000))-hm));
            %k = find(abs(W(linspace(extremePrev,i,50000))-hm)<1e-12,1);
            fwhm_list = [fwhm_list extremePrev + (k-1)/50000*(i-extremePrev) - fwhm_temp];
            fwhm_temp = 0;
        else
            hm = (u(i)+u(extremePrev))/2;
            [M,k] = min(abs(u(linspace(extremePrev,i,50000))-hm));
            %k = find(abs(W(linspace(extremePrev,i,50000))-hm)<1e-12,1);
            fwhm_temp = extremePrev + (k-1)/50000*(i-extremePrev);
        end
        extremePrev = i;
    end 
end
out = norm(fwhm_list);
fprintf('***FWHM: %.7e\n\n', out);

end

