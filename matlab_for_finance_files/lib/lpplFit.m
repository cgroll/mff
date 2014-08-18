function params = lpplFit(prices, initParams, optSettings)
   %lpplFit fit LPPL function to price series
   %
   % The function shall fit a LPPL model to the given price series.
   %
   % Args:
   %  prices         1xn vector of prices
   % Optional Args:
   %  initParams     1x8 vector of initial parameters for
   %                 optimization
   %  optSettings    structure created by optimset to change
   %                 optimization settings
   %
   % Output:
   %  params         1x8 vector of fitted parameters
   
   % specify default values for initParams
   defaultParams = [0 -0.01 (numel(prices)+100) 0.1 1 1 1 1];
   
   % set default values for initParams if function was called without
   % optional parameter input
   if ~exist('initParams', 'var')
      % initial guess for optimization routine
      initParams = defaultParams;
   end
   
   % set default values for optSettings structure, in case function
   % was called without
   if ~exist('optSettings', 'var')
      % specify options for optimization
      optSettings = optimset('display', 'off', 'TolX', 1e-18, ...
         'TolFun', 1e-14,...
         'MaxFunEvals', 4e+8, 'MaxIter', 800,...
         'algorithm', 'interior-point');
   end
   
   % make sure that initParams does not contain NANs
   if containsNans(initParams)
      fprintf('\ninit params had been adjusted\n')
      initParams = defaultParams;
   end
   
   function bol = isAdmissible(params)
      % test for singularity within sample period
      bol = (params(3)/params(8) > numel(prices));
   end
   
   % make sure that sigularity occurs only after current data sample
   if ~isAdmissible(initParams)
      initParams(3) = ceil(numel(prices)*initParams(8))+1;
   end
   
   %% specify optimization settings
   
   % anonymous function to calculate lppl values
   evalLppl = @(x, params) ...
      params(1) + params(2)*(params(3)-params(8)*x).^params(4).*...
      (1+params(5)*cos(params(6).*...
      log(params(3)-params(8)*x)+params(7)));
   
   % define objective function for minimization: mean squared error
   % function
   errFun = @(params, x, prices)...
      sum((prices(:) - evalLppl(x(:), params)).^2);
   
   % bounds on parameters:
   % second parameter must be negative to flip function upside down
   %
   % singularity may only occur after last sample point
   %
   % oscillation frequency should be reasonable small / large
   %
   % shape should be smallest linearly, and not too steep
   %
   % amplitutde should be reasonable
   %
   % shifting the phases of cosine shall only allow one full phase
   
   lb = [0.01; -inf; (numel(prices)+1); 0; -0.98; -40; 0; 0.01];
   ub = [inf; -0.01; inf; 1; 0.98; 40; 2*pi; inf];
   
   % non-linear constraint: singularity may only occur outside of
   % sample
   function [c, ceq] = nonLCFun(params, x, prices)
      c = numel(prices)+1-floor(params(3)/params(8));
      ceq = [];
   end
   fhandle = @nonLCFun;
   
   % conduct optimization
   params = fmincon(errFun, initParams, [], [], [], [],...
      lb, ub, fhandle, optSettings, (1:numel(prices)), prices);
   
end

function bol = containsNans(params)
   bol = any(isnan(params));
end



