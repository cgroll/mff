function joinedTable = processDataWithNaNs(dateBeg, dateEnd, tickerSymbs)
%
% Input:
%   dateBeg     same format as for hist_stock_data
%   dateEnd     same format as for hist_stock_data
%   tickerSymbs     1xn cell array of ticker symbol strings


% download data
stockStructure = [];
for ii=1:length(tickerSymbs)
   currentStock = hist_stock_data(dateBeg, dateEnd, tickerSymbs{1, ii});
   stockStructure = [stockStructure currentStock];
end

% call joinStockPriceSeries
joinedTable = joinStockPriceSeries(stockStructure);

end

