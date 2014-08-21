function joinedTableSorted = getPrices(dateBeg, dateEnd, tickerSymbs)
%
% Input:
%   dateBeg     same format as for hist_stock_data
%   dateEnd     same format as for hist_stock_data
%   tickerSymbs     1xn cell array of ticker symbol strings
%
% Output:
%   joinedTable     a mxn table of stock prices for multiple stocks, with
%                   all dates that occur in at least one stock price series
%                   and missing observations filled with NaNs. Dates are
%                   sorted ascending and stored as strings in the RowNames
%                   property of the table.

% download data
stockStructure = [];
for ii=1:length(tickerSymbs)
   currentStock = hist_stock_data(dateBeg, dateEnd, tickerSymbs{1, ii});
   stockStructure = [stockStructure currentStock];
end

% call joinStockPriceSeries
joinedTable = joinStockPriceSeries(stockStructure);

joinedTableSorted = sortrows(joinedTable, 1);

% get dates as row names
dats = joinedTableSorted{:, 1};
joinedTableSorted(:, 1) = [];
joinedTableSorted.Properties.RowNames = dats;

end

