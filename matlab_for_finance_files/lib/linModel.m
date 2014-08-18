function [intcepts betas tVals] = linModel(Y,X,coeffNull)
% linModel computes parameter estimates and t-statistics

% The function computes OLS estimated coefficients in a linear
% model without intercept. Given parameter estimates,
% t-statistics for given null hypothesis coeffNull are
% calculated.

% Input:
%   Y       Nxreps matrix of reps different observed dependent
%           variable values of sample size N
%   X       Nxreps matrix of reps different explanatory variable
%           values of sample size N
%   coeffNull   1x2 vector with intercept and beta

% Output:
%   intcepts    1xreps vector of estimated intercepts
%   betas       1xreps vector of estimated coefficient values
%   tVals       2xreps vector of associated t-statistics 

[N reps] = size(Y);   % sample size / number of repetitions
intcepts = zeros(1,reps);
betas = zeros(1,reps);

for ii=1:reps
    % get current variables
    currY = Y(:,ii);
    currX = [ones(N,1) X(:,ii)];
    
    % get parameter estimates
    params = currX\currY;
    intcepts(ii) = params(1);
    betas(ii) = params(2);
    
end

% extract residuals
%resids = Y - repmat(betas,N,1).*X;
resids = Y - repmat(intcepts,N,1) - ...
    repmat(betas,N,1).*X;

% estimate standard deviation of innovations
stdDevs = sqrt(sum(resids.^2)/(N-2));

% estimate standard errors
sErr_intcepts = zeros(1,reps);
sErr_betas = zeros(1,reps);
for ii=1:reps
    sErr_intcepts(ii) = stdDevs(1,ii) * sqrt((X(:,ii)'*X(:,ii))...
        /(N*((X(:,ii)'*X(:,ii))-N*mean(X(:,ii))^2)));
    sErr_betas(ii) = stdDevs(1,ii) * sqrt((X(:,ii)'*X(:,ii)-...
        N*mean(X(:,ii))^2)^(-1));
end

tVals = (intcepts-coeffNull(1))./sErr_intcepts;
tVals = [tVals; (betas-coeffNull(2))./sErr_betas];
    
    