function joinedTable = joinMultipleTables(cellOfTables)
%
% Input:
%   cellOfTables    1xn cell array with tables in the individual entries
%
% Output:
%   joinedTable     single table with single dates column and one column
%                   for each stock. Individual stock tables should be
%                   combined with outer join. 

nTables = length(cellOfTables);

joinedTable = cellOfTables{1, 1};
for ii=2:nTables
    joinedTable = outerjoin(joinedTable, cellOfTables{1, ii},...
        'MergeKeys',true);
end

end