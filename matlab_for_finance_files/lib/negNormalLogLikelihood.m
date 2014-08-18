function llh = negNormalLogLikelihood(param,data)
% the function calculates the negative loglikelihood of the normal distribution
%
% Input:    param   2x1 vector of parameters mu and sigma
%           data    nx1 vector of univariate observations

nObs = numel(data);
llh = 0;
for ii=1:nObs
    llh = llh + 1/2*log(2*pi)+log(param(2))+(data(ii,1)-param(1))^2/...
        (2*param(2)^2);
end
