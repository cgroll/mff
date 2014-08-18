function [Dates Prices] = processData(data)
% processes the data structure returned by hist_stock_data

% inputs
%   data structure of hist_stock_data

% ouputs
%   Dates       vector with serial dates, earliest first
%   Prices      associated price vector

% get number of assets
nAss = numel(data);

% determine maximum number of coinciding dates
for ii=1:nAss
    numDates(1,ii) = numel(data(1,ii).Date);
end
maxNum = min(numDates);

% transform date strings to serial date numbers
for ii=1:nAss
    SerialDates = zeros(numDates(1,ii),1);
    for jj=1:numel(SerialDates)
        SerialDates(jj,1) = datenum(data(1,ii).Date(jj,1),'yyyy-mm-dd');
    end
    data(1,ii).Date = SerialDates;
end

% init dates and Prices matrices
dates = zeros(maxNum,1);
Prices = zeros(maxNum,nAss);

% start with last observation
indCount = ones(1,nAss);

% counter for number of equal dates
counter = 0;
while(indCount <= numDates)
    
    % get current dates
    for ii=1:nAss
        currDates(1,ii) = data(1,ii).Date(indCount(1,ii));
    end
    
    % if all dates are equal
    if(sum(currDates == min(currDates))==numel(currDates))
        
        % get current Prices
        for ii=1:nAss
            currPrices(1,ii) = data(1,ii).Close(indCount(1,ii));
        end
        
        % extract current dates and prices
        dates(counter+1,1) = currDates(1,1);
        Prices(counter+1,:) = currPrices(1,:);
        
        % increase indices by 1
        indCount = indCount + ones(1,nAss);
        
        % increase counter by 1
        counter = counter + 1;
    
    % if not all dates are equal    
    else
        % find (first) asset with maximum date
        [maxVal ind] = max(currDates);
        
        % decrease date by increasing index of dates vector
        indCount(1,ind) = indCount(1,ind) + 1;
    end
end
     
% skip surplus rows
dates = dates(1:counter);
Prices = Prices(1:counter,:);

% flip upside down to get earliest observations first
Dates = flipud(dates);
Prices = flipud(Prices);








