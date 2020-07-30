%% Morris Method Sensitivity Calculator
%This script uses the morris_sa.m function to create paths in parameter
% space, compute a defined metric from the steady-state solution, and
% compute the sensitivities mu and sigma.
%The scipt outputs a vector for mu and sigma and plots the sensitivity
% scatter plot.

%Current Tests:
% test = 1 --- spike maximum
% test = 2 --- spike width (full width at half-maximum)
% test = 3 --- distance between spike peaks
% test = 4 --- spike minimum

clc; clear; close all;

%%
test = 4;
fun = @(x) morrisFun(0.17*x(1)+1, 1e-2*x(2), 1e-5*x(3), 5*x(4)+1,test);
[mu, sigma] = morris_sa(fun, 4, 20);

%%
figure(1)
scatter(mu/max(mu),sigma/max(sigma))
textx = mu/max(mu)+.01; texty=sigma/max(sigma);
xlim([0 1.05])
ylim([0 1.05])
set(gca, 'fontsize', 18)
set(gca, 'tickLabelInterpreter', 'latex');
text(textx,texty,{'$\rho$', '$\eta$', '$\gamma$', '$\alpha$'}, 'interpreter', 'latex', 'fontsize', 20)
xlabel('$\mu$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$\sigma$', 'interpreter', 'latex', 'fontsize', 20)

if test == 2
    title('Morris Method --- Avg FWHM', 'interpreter', 'latex', 'fontsize', 20)
elseif test == 3
    title('Morris Method --- Avg Distance btw Spikes', 'interpreter', 'latex', 'fontsize', 20)
elseif test == 4
    title('Morris Method --- Amenity Spike Minimum', 'interpreter', 'latex', 'fontsize', 20)
else
    title('Morris Method --- Amenity Spike Amplitude', 'interpreter', 'latex', 'fontsize', 20)
end