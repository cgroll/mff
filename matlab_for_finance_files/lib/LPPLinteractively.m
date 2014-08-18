function params = LPPLinteractively(prices, dates)
   % analyses given price series with LPPL fit
   
   % The function is fitting a LPPL model to a subperiod of the
   % given price series. When executed, the user is able to decide,
   % whether he wants to specify the subperiod via graphical input
   % or by entry of strings for beginning and ending. Once the
   % subperiod is specified, an exponential function is fitted to
   % the logarithmic prices, in order to show whether the price
   % series indicates a super-exponential growth. Afterwards a LPPL
   % function is fitted to the data. Both models are compared on the
   % basis of the MSE. The LPPL model returns the date of the most
   % probable regime shift.
   
   % Input:
   %   prices      nx1 vector of normal prices
   %   dates       nx1 vector of associated serial dates
   
   % Output:
   %
   
   % transform prices to log prices
   prices = log(prices);
   
   % plot data
   clf
   plot(dates, prices)
   datetick 'x'
   title('logarithmic price series')
   shg
   
   % short break to give user time to look at price series
   pause(2)
   
   % let user choose between graphical input or dates
   commandwindow
   reply = input('\nDo you want graphical input? Y/N [Y]: ', 's');
   
   if (strcmp(reply, 'N'))
      % ask for strings to determine sample bounds
      reply1 = input('\nPlease type beginning as dd-mmm-yy:', 's');
      intS = datenum(reply1);
      reply2 = input('\nPlease type ending as dd-mmm-yy', 's');
      intE = datenum(reply2);
      
   else
      % request graphical input from user
      % plot data
      shg
      title('Click twice on price chart')
      %fprintf('\nPlease click twice with your mouse to specify period\n');
      
      [x, y] = ginput(2);
      
      % round coordinates to integer, and sort in right order
      if (x(1) < x(2))
         intS = ceil(x(1));  % coordinates could be decimal
         intE = floor(x(2));
      else
         intS = ceil(x(2));
         intE = floor(x(1));
      end
   end
   
   % transform start and end points to indices of dates vector
   indS = find(dates >= intS, 1);
   indE = find(dates <= intE, 1, 'last');
   
   % output in words
   fprintf(['\nYou have specified the following subperiod:\n' ...
      datestr(intS) ' to ' datestr(intE) '\n']);
   
   % highlight in graphic
   hold on;
   plot(dates(indS:indE), prices(indS:indE), ...
      'LineWidth', 2, 'Color', 'r');
   title('log prices with chosen subperiod');
   shg
   
   % fit regression line for exponential fit
   Y = prices(indS:indE);
   X = [ones(numel(Y), 1) (1:numel(Y))'];
   beta = X\Y;
   
   % include fit into graphics
   vals = X*beta;
   plot(dates(indS:indE), vals, 'g'); shg
   
   % calculate MSE
   mse1 = sum((Y(:)-vals(:)).^2);
   
   % fit LPPL function
   params = lpplFit(Y);
   
   % get values of LPPL fit
   [vals derivs] = lpplFunc(params);
   
   % create grid for values of LPPL up to critical point
   grid = (1:(params(3)/params(8)-1));
   
   % include in graphics
   yLimits = get(gca, 'yLim');
   
   % test, if values of LPPL fit exceed existing limit of x-axis
   if((indS+grid(end))>numel(dates))
      % plot values only up to end of price chart
      plot(dates(indS:end), vals(1:numel(dates(indS:end))), 'k')
   else
      plot(dates(indS+grid), vals, 'k');
   end
   set(gca, 'yLim', yLimits);
   
   % calculate MSE
   mse2 = sum((Y(:)'-vals(1:numel(Y))).^2);
   
   if(mse2<mse1)
      title(['The LPPL fit is better with factor '...
         num2str(mse2/mse1)])
   else
      title(['The LPPL fit is worse with factor '...
         num2str(mse2/mse1)])
   end
   shg
   
   
   % fprintf(['While the exponential function has a MSE of %3.3d, \n'...
   %     'the LPPL function has a MSE of %3.3d\n'], MSE1, MSE2);
   
   
   % calculate index of critical point
   crPoint = indS+floor(params(3)/params(8));
   
   if(crPoint < numel(prices))
      crPointString = datestr(dates(crPoint));
      fprintf(['\nThe estimated critical point was at ' ...
         crPointString ', \nwhich was ' ...
         num2str(indS+grid(end)-indE) ' business days ahead.\n']);
      line(dates(indS+grid(end))*[1 1], yLimits, 'Color', 'k');
   else
      % days until critical point
      daysIntoFuture = crPoint-numel(prices);
      
      % get sufficient future business days
      futureBusDays = ...
         busdays(dates(end), dates(end)+2*daysIntoFuture);
      crPointString = datestr(futureBusDays(1+daysIntoFuture));
      fprintf(['\nThe critical point will be at '...
         crPointString ', \nwhich is ' num2str(daysIntoFuture) ...
         ' business days ahead.\n']);
   end
   
