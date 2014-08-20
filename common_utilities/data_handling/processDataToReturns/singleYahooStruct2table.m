function t = singleYahooStruct2table(yahooStruct)
% 
% Input:
%   yahooStruct     single structure as returned by hist_stock_data
%
% Output:
%   t               table with Dates in first column, adjusted closing
%                   prices in second column and name of second column
%                   equals ticker symbol


% create valid column name
varname = createValidName(yahooStruct.Ticker);

t = table(yahooStruct.Date, ...
    yahooStruct.AdjClose,...
    'VariableNames', {'Dates', varname});

end