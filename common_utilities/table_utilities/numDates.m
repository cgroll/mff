function serialDates = numDates(priceTable)
% 
% Input:
%   priceTable  nxm table with dates stored as row names
%
% Output:
%   serialDates     nx1 vector of numeric dates

serialDates = datenum(priceTable.Properties.RowNames, 'yyyy-mm-dd');

end