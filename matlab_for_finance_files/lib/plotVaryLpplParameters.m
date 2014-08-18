function plotVaryLpplParameters(origParameters, coeffNo, variations)
   %plotVaryLpplParameters plots the variation of lppl in a given
   % parameter
   %
   % Args:
   %  origParameters    1x8 vector of original parameters
   %  coeffNo           number of coefficient for variation
   %  variations        variations of coefficient
   %
   % Side Effect:
   %  the function opens a plot
   
   colors = ['b', 'k', 'g', 'y'];
   
   function grid = getLpplGrid(params)
      grid = (1:0.2:(params(3)/params(8)-1));
   end
   
   % init anonymous function to get function values
   evalLppl = @(x, params) ...
      params(1) + params(2)*(params(3)-params(8)*x).^params(4).*...
      (1+params(5)*cos(params(6).*...
      log(params(3)-params(8)*x)+params(7)));
   
   
   % get grid for original data
   gridOrig = getLpplGrid(origParameters);
   
   % get associated values
   valsOrig = evalLppl(gridOrig, origParameters);
   
   plot(gridOrig, valsOrig, '-r')
   hold on;
   
   params = origParameters;
   for ii=1:numel(variations)
      % get modified parameter vector
      params(coeffNo) = variations(ii);
      
      % get new grid
      grid = getLpplGrid(params);
      
      % get new values
      vals = evalLppl(grid, params);
      
      % include into plot
      plot(grid, vals, 'Color', colors(ii))
   
   end
   legend('original', 'lowest', 'lower', 'higher', 'highest', ...
      'location', 'NorthWest')
end

