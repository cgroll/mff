function rets = price2retWithHolidays(prices)
%
% Input:
%   prices  nxm matrix or table of prices
%
% Output:
%   rets    (n-1)xm matrix of logarithmic returns

nStocks = size(prices, 2);

% get log prices
logPrices = log(prices);

% find missing values and replace with previous observation
missingPrices = isnan(logPrices);
nansToReplace = logical([zeros(1, nStocks); missingPrices(2:end, :)]);
replaceWith = logical([missingPrices(2:end, :); zeros(1, nStocks)]);

pricesImputed = logPrices;
pricesImputed(nansToReplace) = logPrices(replaceWith);

% calculate returns
rets = diff(pricesImputed);

% fill in NaNs again
rets(missingPrices(2:end, :)) = NaN;

end