function [ mu, order ] = morris_sa( studied_function, nfac, npoint )
% This function applies the method of Sensitivity Analysis called the
% Morris method.
%
% Sensitivity analysis is used to estimate the influence of uncertainty
% factors on the output of a function.
% The Morris method is sometimes referenced to as a qualitative method : it
% gives rough estimations with a limited number of calculations.
% The Morris method can be used to simplify a function, as a first step. It
% can identify the factors with a low influence which can be fixed.
% For further information : Saltelli, A., Tarantola, S., Campolongo, F., and Ratto, M. (2004). Sensitivity Analysis in Practice - A Guide to Assessing Scientific Models. Wiley.
%
% This function respects the recommendations from :
% Sohier, H., Farges, J. L., & Piet-Lahanier, H. (2014, August). Improvement of the Representativity of the Morris Method for Air-Launch-to-Orbit Separation. In The 19th IFAC World Congress.
% (the elementary effects are calculated with radial points sampled in a Latin hypercube and with large variations)
%
% INPUTS
% 1) studied_function : Anonymous function to be analyzed. It must be a 
% function with one multidimensional input x. x must represent the values 
% of the uncertainty factors in the quantiles hyperspace (the i-th 
% coordinate of x is not the actual value of the i-th factor, but the 
% corresponding value of its cumulative distribution function.
% To adapt your function, first calculate the actual values of 
% the factors by applying the inverse of their cumulative distribution 
% function to x; Matlab includes such inverses: 
% mathworks.com/help/stats/icdf.html ).
% Matlab functions : if the function to be analyzed is a Matlab function 
% called "test_function", input the following argument : 
% @(x)test_function(x)
% 2) nfac : Number of factors of uncertainty of the studied function
% 3) npoint : Number of points to be sampled in the quantiles hyperspace.
% More points = better accuracy = more calculations
% Number of calculations = npoint*(nfac+1)
% Recommended value = 10
%
% OUTPUTS
% 1) mu : Average of the absolute values of the elementary effects of
% each factor, sorted in descending order.
% 2) order : Indexes of the factors. Consider fixing the last factors,
% those whose "mu" is very close to zero.
%
% EXAMPLE :
% The file "test_function" is the modified Sobol test function. You can
% test the Morris method with the following line:
% [mu, order] = morris_sa(@(x)test_function(x), 20, 10)
% In each dimension, the interval [0;1] is discretized in npoint values.
% Calculation of the minimum and maximum values:
delta=1/npoint;
mini=delta/2;
maxi=mini+delta*(npoint-1);
% Sampling the points (total: npoint) in a Latin hypercube. The discretized
% values are first handled as integers, and they are then normalized. In 
% the matrix "points", each line is a point and each column is a factor.
coord = 0:npoint-1;
for i=1:nfac
    points(1:npoint,i) = coord(randperm(length(coord)));
end
points = points/(npoint-1)*(maxi-mini)+mini;
% Simulation runs
% All the simulation outputs are in table_outputs. In table_outputs, each
% line represents the outputs around a different point. The first column
% represents the outputs at the sampled points. The value at the i-th line
% and (1+j)-th column represents the output when a variation is applied to
% the j-th factor at the i-th point.
for i=1:npoint
    table_outputs(i,1) = studied_function(points(i,:)); % Output at the sampled point.
    for j=1:nfac
        if points(i,j) < 0.5 % If the coordinate is smaller than 0.5, a positive variation is applied
            table_outputs(i,1+j) = studied_function([points(i,1:j-1) points(i,j)+0.5 points(i,j+1:nfac)]); % Output after the variation of the j-th factor.
            table_ee(i,j) = (table_outputs(i,1+j)-table_outputs(i,1))/0.5; % Elementary effect of the j-th factor.
        else % If the coordinate if larger than 0.5, a negative variation is applied
            table_outputs(i,1+j) = studied_function([points(i,1:j-1) points(i,j)-0.5 points(i,j+1:nfac)]);
            table_ee(i,j) = (table_outputs(i,1+j)-table_outputs(i,1))/(-0.5);
        end
    end
end
% Estimation of the factors influence with the average of the absolute
% values of the elementary effects
for j=1:nfac
    mu_temp(j) = mean(abs(table_ee(:,j)));
end
[mu, order] = sort(mu_temp,'descend'); % Ordering.
end
