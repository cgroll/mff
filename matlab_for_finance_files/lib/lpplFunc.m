function [vals derivs] = lpplFunc(params)
%lpplFunc calculates function values and gradients for lppl function
%
% The function computes the values of a LPPL model as well as the
% derivate at each point. Since the LPPL model has a finite time
% singularity, the values are calculated from day 1 to one day
% before the singularity.
%
% Args:
%   params      1x8 matrix containing the parameter values 
%
% Output:
%   vals        1x(params(3)-1) vector of values
%   derivs      1x(params(3)-1) vector of derivatives

% initialize grid
grid = (1:(params(3)/params(8)-1));

% calculate values
vals = params(1) + params(2)*(params(3)-params(8)*grid)...
    .^params(4).*...
    (1+params(5)*cos(params(6)*...
    log(params(3)-params(8)*grid)+params(7)));

% calculate derivs
% derivs = params(2)*(params(3)-grid).^(params(4)-1).*...
%     (params(4)+params(4)*params(5)*...
%     cos(params(6)*log(params(3)-grid)+params(7))+...
%     params(6)*params(5)*sin(params(6)...
%     *log(params(3)-grid)+params(7)));
% derivs = -derivs;

derivs = gradient(vals);

% Note: the function returned negative derivative values in
% periodes of increase. However, the absolute value of the
% derivatives seemed correct, since I have checked them
% numerically with derivs2 = gradient(vals). Hence, I have
% entailed the last assignement. But I did not find the mistake
% so far!!