function frames = LPPLmovie(prices, dates, timeSteps)
   %lpplMovie store best fits of lppl function in frames for movie
   %
   % The function lets the user interactively choose a subperiod of
   % the given price series, which shall be analyzed stepwise with
   % fitted LPPL functions.
   %
   % The idea is, to simulate the past, in order to examine at each
   % point what kind of signal the LPPL model was giving.
   %
   % First, the user determines the overall beginning point and the
   % overall ending point. Then, the user determines the end point
   % for the first sample that is used to fit a LPPL model. This
   % input, of course, should not be chosen too closely to the
   % overall beginning. At last, the user has to determine the best
   % point, where a strategy based on LPPL models could have warned
   % the investor.
   %
   % Args:
   %  prices      nObs x 1 vector of historic prices
   %  dates       nObs x 1 vector of associated numeric dates
   %  timeSteps   the number of days the have to pass until the model
   %              will be re-estimated
   %
   % Outputs:
   %  frames      matrix of MATLAB movie frames
   %
   % The complication in this code arise from the fact that we have
   % three different scales to measure each point:
   %  - the "absolute" date of the point used in the plot
   %  - the "relative" date of each point, with respect to the
   %    beginning of the overall sample (this is equal to the index
   %    of the entry in the date vector)
   %  - the "time passed" as distance to the beginning of the chosen
   %    sample period
   
   %% Create plot for user interaction
   % transform prices to log prices
   prices = log(prices);
   
   % plot data
   plot(dates, prices)
   datetick 'x'
   title('logarithmic price series')
   
   %% User interaction: determine overall sample for examination
   % request graphical input from user
   shg
   title('Determine beginning of analysis per mouse click')
   [x(1), y(1)] = ginput(1);
   
   shg
   title('Determine end of analysis per mouse click')
   [x(2), y(2)] = ginput(1);
   
   % round values to real numbers
   dateBeg = ceil(x(1));
   dateEnd = floor(x(2));
   
   % output in words
   fprintf(['\nYou have specified the following subperiod:\n' ...
      datestr(dateBeg) ' to ' datestr(dateEnd) '\n']);
   
   % transform start and end points to indices of dates vector
   function relativeDate = indexOfDate(absDate)
      relativeDate = find(dates >= absDate, 1);
   end
   
   relBeg = indexOfDate(dateBeg);
   relEnd = indexOfDate(dateEnd);
   
   % highlight in graphic
   hold on;
   plot(dates(relBeg:relEnd), prices(relBeg:relEnd), ...
      'LineWidth', 2, 'Color', 'r');
   title('log prices with chosen subperiod');
   hold off;
   shg
   
   % get axis limits
   xLimits = get(gca, 'xLim');
   yLimits = get(gca, 'yLim');
   
   %% User interaction: determine starting point and best point
   title('Determine interval for first fit')
   [x(3), y(3)] = ginput(1);
   
   title('Determine best possible exit point')
   [x(4), y(4)] = ginput(1);
   
   % transform to indices of dates vector
   relEndPointFirstSample = indexOfDate( ceil(x(3)) );
   relBestExit = indexOfDate( floor(x(4)));
   
   function createLpplAnalysisPlot(relBeg, currRelSampleEnd)
      % get present price series
      currPrices = prices(relBeg:currRelSampleEnd);
      
      % fit LPPL function
      params = lpplFit(currPrices);
      
      % get values of LPPL fit
      [vals, ~] = lpplFunc(params);
      
      % create relative grid for values of LPPL, from beginning of
      % current sample period up to critical point
      relGridForLppl = (1:(params(3)/params(8)-1) );
      
      % plot all prices in background
      plot(dates, prices, 'Color', [1 0.9 0.9])
      hold on;
      
      % plot prices used in current sample period
      plot(dates(relBeg:currRelSampleEnd), ...
         prices(relBeg:currRelSampleEnd));
      datetick 'x'
      set(gca, 'xLim', xLimits, 'yLim', yLimits);
      
      % test, if values of LPPL fit exceed existing limit of x-axis
      isCriticalPointOutsideOfAxis = ...
         ((relBeg + relGridForLppl(end)) > numel(dates));
      
      if isCriticalPointOutsideOfAxis
         % plot values only up to end of price chart
         remainingDatesFromCurrentStartPoint = ...
            numel(dates(relBeg:end));
         
         plot(dates(relBeg:end), ...
            vals(1:remainingDatesFromCurrentStartPoint), 'g')
         datetick 'x'
         set(gca, 'xLim', xLimits, 'yLim', yLimits);
      else
         plot(dates(relBeg+relGridForLppl), vals, 'g');
         datetick 'x'
         set(gca, 'xLim', xLimits, 'yLim', yLimits);
         line(dates(relBeg + relGridForLppl(end))*[1 1],...
            yLimits, 'Color', 'k');
      end
      set(gca, 'yLim', yLimits);
      
      % calculate index of critical point
      relCriticalPoint = relBeg + floor(params(3)/params(8));
      title(['Days to critical point: '...
         num2str(relCriticalPoint - currRelSampleEnd)]);
      
      hold off;
   end
   
   
   frameCounter = 1;
   for currRelSampleEnd = relEndPointFirstSample:timeSteps:relEnd
      % create plot with analysis
      createLpplAnalysisPlot(relBeg, currRelSampleEnd)
      
      % capture frame of present plot
      frames(frameCounter) = getframe;
      
      % increase frame counter
      frameCounter = frameCounter+1;
   end
   
   %% What did LPPL model suggest at best exit point?
   createLpplAnalysisPlot(relBeg, relBestExit)
   
   frames(frameCounter) = getframe;
   fprintf('Scene is in the can\n')
end