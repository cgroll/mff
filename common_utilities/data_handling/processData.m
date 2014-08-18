function [dates prices] = processData(data)
%processData processes a hist_stock_data output structure 
%
% The function receives a structure array as returned by function
% hist_stock_data as input. It then extracts adjusted closing
% prices for all dates where no asset is lacking its observation.
% These adjusted closing prices are stored in a matrix, and
% returned by the function. Moreover, the corresponding dates are
% returned in numerical date format in a column vector. Both
% prices and dates are flipped upside down, in order to list most
% recent observations at the bottom. This way, most recent
% observations will appear at the right end in graphics.
%
% Args:
%   data         structure array as returned by
%                hist_stock_data 
%
% Outputs:
%   dates        nObs x 1 vector of serial dates,
%                most recent observation last
%   prices       nObs x nAss matrix of asset prices 

% get number of assets
nAss = numel(data);

% maximal possible number of coinciding dates is smallest sample
% size - required for preallocation
numDates = ones(1, nAss);
for ii=1:nAss
    numDates(1, ii) = numel(data(1, ii).Date);
end
maxNum = min(numDates);

% each asset contains a cell array with string dates in field
% .Date that is transformed to vector of serial date numbers
for ii=1:nAss
    serialDates = datenum(data(1, ii).Date,'yyyy-mm-dd');
    data(1, ii).Date = serialDates;
end

% preallocation: init dates and prices matrices
dates = zeros(maxNum, 1);
prices = zeros(maxNum, nAss);

% start with first indices corresponding to most recent
% observations 
indCount = ones(1, nAss);

% counter for number of equal dates
nEqualDaysFound = 0;
currDates = ones(1, nAss);
currPrices = ones(1, nAss);
while(indCount <= numDates)
    
    % get horizontal vector of current date combination
    for ii=1:nAss
        currDates(1, ii) = data(1, ii).Date(indCount(1, ii));
    end
    
    % test for equal dates for current date combination
    areEqual = allVectorEntriesEqual(currDates);
    
    if (areEqual) % store data of admissable date combination
        
        % get current prices
        for ii=1:nAss
            currPrices(1, ii) = ...
                data(1, ii).AdjClose(indCount(1, ii));
        end
        
        % extract current dates and prices
        dates(nEqualDaysFound + 1, 1) = currDates(1, 1);
        prices(nEqualDaysFound + 1,:) = currPrices(1,:);
        
        % increase indices by 1
        indCount = indCount + ones(1, nAss);
        
        % increase nEqualDaysFound by 1
        nEqualDaysFound = nEqualDaysFound + 1;
    
    else  % if not all dates are equal    
        % find all assets with date larger than minimun date
        tooRecentDates = currDates > min(currDates);
        
        % decrease dates by increasing indices of dates vector
        indCount(tooRecentDates) = indCount(tooRecentDates) + 1;
    end
end
     
% skip surplus rows from preallocation
dates = dates(1:nEqualDaysFound);
prices = prices(1:nEqualDaysFound,:);

% flip data upside down to get oldest observations first
dates = flipud(dates);
prices = flipud(prices);




