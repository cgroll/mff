function retsTable = price2retWithHolidays(prices)
%
% Input:
%   prices  nxm matrix or table of prices
%
% Output:
%   retsTable    (n-1)xm table of logarithmic returns

% get missing values
missingValues = isnan(prices{:,:});

% get log prices
logPrices = log(prices{:,:});
pricesImputed = imputeWithLastDay(logPrices);

% impute once again?
% pricesImputed = imputeWithLastDay(pricesImputed);

% calculate returns
rets = diff(pricesImputed);

% fill in NaNs again
rets(missingValues(2:end, :)) = NaN;

% embed returns in table meta-data
retsTable = embed(rets, prices(2:end, :));

end