function xAsTable = embed(x, t)
% add meta-data from table t to matrix x
%
% Inputs:
%   x       nxm matrix 
%   t       nxm table
%
% Output:
%   xAsTable    table

xAsTable = array2table(x);
xAsTable.Properties.RowNames = t.Properties.RowNames;
xAsTable.Properties.VariableNames = t.Properties.VariableNames;
end
