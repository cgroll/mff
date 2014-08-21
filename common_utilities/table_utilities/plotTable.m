function plotTable(tableToPlot)
serialDates = datenum(tableToPlot.Properties.RowNames);
plot(serialDates, tableToPlot{:, :})
datetick 'x'
end