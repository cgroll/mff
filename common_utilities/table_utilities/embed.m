function t = embed(x, t)
% add meta-data from table t to matrix x
%
% Inputs:
%   x       nxm matrix 
%   t       nxm table
%
% Output:
%   xAsTable    table

t{:,:} = x;

end
