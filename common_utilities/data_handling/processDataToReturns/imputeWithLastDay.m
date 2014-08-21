function values = imputeWithLastDay(values)
%
% Inputs:
%   values  nxm matrix
%
% Outputs:
%   values  nxm matrix with imputed values

% Replace NaN with observation of last day. 
nRows = size(values, 2);

% find missing values and replace with previous observation
missingValues = isnan(values);

nansToReplace = logical([zeros(1, nRows); missingValues(2:end, :)]);
replaceWith = logical([missingValues(2:end, :); zeros(1, nRows)]);

values(nansToReplace) = values(replaceWith);

end