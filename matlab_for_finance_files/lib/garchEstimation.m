function nllh = garchEstimation(params,data,initVals)
% calculates neg. log-likelih. of GARCH model

% The function calculates the negative log-likelihood associated
% with a GARCH(1,1) model with parameters given by params and
% initial values of Y and sigmas specified by initVals.

% Inputs:
%   params      1x3 vector of GARCH model parameters
%   data        1xn vector of observed realizations
%   initVals    1x2 vector of initial values

% Output:
%   nllh        negative log-likelihood associated with given
%               input values

% first step: retrieve sigmas

% preallocate reconstructed sigmas
retrieveSigmas = zeros(numel(data),1);
retrieveSigmas(1) = initVals(2);

% reconstruction
for ii=2:numel(data)
    % current sigmas depend on last sigma and last observation
    % through parameter values of the model
    retrieveSigmas(ii) = sqrt(params(1)+...
        params(2)*retrieveSigmas(ii-1)^2 +...
        params(3)*data(ii-1)^2);
end

% calculate negative log-likelihood
nllh = sum(0.5*log(retrieveSigmas.^2*2*pi)+...
    0.5*(data.^2./retrieveSigmas.^2));