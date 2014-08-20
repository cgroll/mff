function joinedTable = joinStockPriceSeries(structArr)
% 
% Input:
%   structArr       structure array as return by hist_stock_data
%
% Output:
%   joinedTable     single table with one dates column and individual
%                   stocks as columns

nStocks = length(structArr);
tableContainer = cell(1, nStocks);

for ii=1:nStocks
    tableContainer{1, ii} = singleYahooStruct2table(structArr(ii));
end

% join tables
joinedTable = joinMultipleTables(tableContainer);

% % parallel computing
% joinedTableContainer = cell(1, nCores);
% for ii=1:nCores
%    joinedTableContainer{ii} = joinMultipleTables(tableContainer(1:nQuarter));
%     
% end

end


