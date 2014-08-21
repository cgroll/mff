function plotTable(tableToPlot)
serialDates = datenum(tableToPlot.Properties.RowNames);
plot(serialDates, tableToPlot{:, 1})
datetick 'x'
end