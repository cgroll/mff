%% Working with financial data: curve fitting
% Christian Groll
%
% Chair of Financial Econometrics, Ludwig-Maximilians-University 
% Munich.  
%
% All rights reserved.

%% Stock price prediction based on curve fitting
% While the previous part was concerned with looking for an
% explanatory variable for stock returns, we now will try to find
% regularities in stock prices that allow to make predictions on
% future price movements. That is, in course of its evolution,
% any stock price seems to follow some trend at some point of
% time. Looking at charts of stock prices one usually might be
% tempted to assume that such trends could be identified in
% real-time, thereby allowing for speculative trading 
% opportunities. The idea in this chapter is to fit certain
% functions to historic stock price paths. Given that the
% function seems to be a good approximation to past prices,
% chance might be that it will still be an approximation in the
% future, so that our function could be used as stock price 
% predictor.
% However, the approach taken here is slightly different. Based
% on curve fitting tools, positive trends in stock prices shall
% be identified. But instead of trying to exactly predict future
% prices, we only try to identify points in time where the
% current dynamic changes. That is, we only try to predict
% break-offs of rising stock prices, without bothering with the
% exact type of regime evolving after the break-off.

%%
% Given that returns fluctuate around a constant positive value,
% prices should exhibit exponential growth. Such growth rates
% best can be seen on logarithmic scale, since they correspond to
% a straight line here. Hence, we first extend the data structure
% with an additional field logPrices. Visualization shows that
% DAX prices tend to exhibit super-exponential growth during
% certain periods.

% get log prices
dax.logPrices = log(dax.prices);

% specify subperiod as strings
begT = '01-Jun-1993';
endT = '29-Jul-1998';

% find indices associated with considered period
indS = find(dax.dates > datenum(begT, 'dd-mmm-yyyy'), 1);
indE = find(dax.dates > datenum(endT, 'dd-mmm-yyyy'), 1);

%%
% Note: it is not possible to access the prices with indexing
% based on the dates of the time series. Hence, dates always have
% to be converted to chronological indices first. However, the
% finance toolbox of MATLAB also includes financial time series
% objects (fints) that can be indexed by date strings. For
% example, myfts({'05/11/99', '05/21/99', '05/31/99'}) extracts
% the values of the fints object myfts at the specified dates.

% create figure window 
close
figure('Position', [50 50 1200 600])

% plot DAX prices with subperiod highlighted
ax(1) = subplot(2, 1, 1);
plot(dax.dates, dax.prices, 'Color', [1 0.8 0.8]);
hold on;
plot(dax.dates(indS:indE), dax.prices(indS:indE)); 
datetick 'x'
title('linear scale')

% plot log DAX prices with subperiod highlighted
ax(2) = subplot(2, 1, 2);
plot(dax.dates, dax.logPrices, 'Color', [1 0.8 0.8]);
hold on;
plot(dax.dates(indS:indE), dax.logPrices(indS:indE)); shg
datetick 'x'
title('logarithmic scale')

% connect axes of both graphs: zooming in applies to both plots
linkaxes([ax(1) ax(2)], 'x');

%% 
% Although it would be easier to fit a straight line to log
% prices we want to estimate to best fitting exponential growth
% for normal prices using an optimization.
% Hence, the exponentially growing function f(x)= a_1*exp(a_2*x)
% shall be fitted to the stock prices. Therefore, parameters a_1 
% and a_2 will be chosen such that the mean squared error between
% the exponential function and the historic price chart is 
% minimized.

% create new grid for subperiod, starting at 1
daysSinceBeg = 1:numel(dax.dates(indS:indE));   % stock market 
        % prices are treated as equidistant, with no distinction
        % between Friday / Monday or Monday / Tuesday

% define exponential function as anonymous function
expFun = @(x, params) params(1)*exp(x.*params(2));

% evaluating exponential function similar to normal functions
fprintf(['Calling the anonymous function according to '...
    'usual syntax\nexpFun(3,[0.5 0.5])\nreturns the value'...
    ' %1.2f.\n'], expFun(3,[0.5 0.5]))

%%

% define mean squared error function as anonymous function
errFun = @(params, x, prices)...
    sum((prices(:) - expFun(x(:), params)).^2);  % for any price 
        % series given by prices and associated x values the
        % error function computes the mean squared error between
        % exponential function with parameters params and the
        % price series

% init guess for optimization
params0 = [dax.prices(indS) ...
    log(dax.prices(indE) - dax.prices(indS))/...
    (dax.dates(indE) - dax.dates(indS))];
        % params(2) chosen so that it fulfills the equation:
        % exp((daysSinceBeg(end)-daysSinceBeg(1))*a_2) 
        %           != prices(end)-prices(1)

% specify options for optimization        
opt = optimset('display', 'off', 'TolX', 1e-18, 'TolFun', 1e-8);

% optimization
[bestParams expMSE] = fminsearch(errFun, params0, opt,...
    daysSinceBeg, dax.prices(indS:indE));

%%
% Note: since the objective function, which shall be minimized,
% also depends on the grid values of x and the given price vector
% prices, these data has to be given as fixed input into the
% optimization, since the optimization shall only be applied to
% the parameter values. 
% Therefore, the parameters of interest have to appear in the
% objective function as one vector and as first input. Additional
% inputs are included in the optimization routine fminsearch as
% additional inputs at last positions. However, this syntax is
% only allowed when the objective function is given as function
% handle to an anonymous function. An example of a similiar
% optimization task involving an already existing MATLAB function
% will be given further below.

% calculate associated exponential function values
expVals = expFun(daysSinceBeg, bestParams);

% include in given figure
subplot(2, 1, 1);
plot(dax.dates(indS+daysSinceBeg), expVals, 'r'); % Note: 
        % dax.dates(indS) + daysSinceBeg does not work, since
        % dax.dates is not numbered consecutively. dax.dates
        % refers to business days, not consecutive days!
xlabel('dates')
ylabel('prices')

subplot(2, 1, 2);
plot(dax.dates(indS+daysSinceBeg), log(expVals), 'r'); shg
xlabel('dates')
ylabel('prices')

% calculate mean squared error on logarithmic scale
mse = sum((dax.logPrices(indS+daysSinceBeg)-log(expVals(:))).^2);

% display mean squared error 
fprintf(['\nThe mean squared error between the exponential fit'...
        ' and\nthe stock price path is %3.4f.\n'], mse);

%%
% With the straight line as benchmark, one can see that the stock
% price path exhibits a convex curvature during the subperiod.
% This pattern indicates super-exponential growth rates. Such 
% growth rates usually are associated with stock market bubbles.
% Our intention now will be to identify evolving stock market
% bubbles, and try to predict the time they burst.
% According to Didier Sornette and his colleagues, stock market
% bubbles can be approximated with super-exponentially growing 
% log-periodic power law (LPPL) functions. These are 
% super-exponentially growing functions with finite-time 
% singularities and oscillating behaviour, given by the formula:
% f(x) = a_1 + a_2*(a_3-x)^(a_4)*
%       (1+a_5*cos(a_6*log(a_3-a_8*x)+a_7).
% In order to get an impression about the appropriateness of a
% LPPL function, we will fit it the subperiod and compare its
% mean squared error to the error of a simple exponential
% fucntion. Furthermore, we will examine whether the date of the
% estimated finite-time singularity could be used as indicator of
% a forthcoming change in regimes.


% fit LPPL model to subperiod
params = lpplFit(dax.logPrices(indS:indE));

% calculate approximation values to stock prices
[vals derivs] = lpplFunc(params);

% create associated grid
grid = dax.dates(indS + ( 1:(params(3)/params(8)-1) )); 
    % Note: params(3)/params(8) denotes the time in business days
    % from beginning of subperiod until finite-time singularity. 

% include in given figure
subplot(2, 1, 2);
plot(grid, vals, 'g'); shg

% include line for finite time singularity
yLimits = get(gca, 'yLim');
line(dax.dates(indS + floor(params(3)/params(8)) )*[1 1], yLimits,...
    'Color', 'k')

% calculate mean squared error on logarithmic scale
mseLppl = sum( (dax.logPrices(indS + daysSinceBeg) -...
    (vals(daysSinceBeg)')).^2);

fprintf(['\nIn contrast to the MSE of ' num2str(mse) ...
    ' obtained before,\n we now get a MSE of only '...
    num2str(mseLppl) '.\n'])

%% 
% When looking at the figure, we can see that the fitted LPPL
% model at the time of the end of the subperiod could indicates 
% an impending regime change, since the critical point given by 
% the finite-time singularity lies only days ahead.

%%
% In order to examine the validity of the LPPL model on further
% stock market indices, you can uncomment the following lines of
% code and interactively conduct experiments on historic data. As
% examples of further accurate subperiod fitting, take a look at
% Hang Seng index from 15-Dec-2004 to 21-Nov-2007, which 
% leads to an estimated regime change 52 business days ahead, or
% the German stock market index from 15-Oct-1992 to 29-Jul-1998.

%% 
% % Interactive examination of further stock market indices.
% 
% %tickerSyms = cell(8, 1);
% tickerSyms = {'^GDAXI';'^SMSI';'^SSMI';...
%     '^OMXSPI';'^NDX';'^DJI';'^HSI';'^SSEC'};
% 
% indexNames = {'DAX'; 'Madrid General';...
%     'Swiss Market'; 'Stockholm General'; 'NASDAQ'; ...
%     'Dow Jones Industrial'; 'Hang Seng';...
%     'Shanghai Composite'};
% 
% for ii=1:numel(tickerSyms)
%     fprintf(['\nIndex investigated: ' indexNames{ii} '\n'])
%     data = hist_stock_data(begT, endT, tickerSyms{ii});
%     
%     if (~isempty(data))
%        [data_dates data_prices] = processData(data);
%        params = LPPLinteractively(data_prices, data_dates)
%        title([indexNames{ii} ' -- Press key to continue'])
%        hold off
%     end
% 
%     
%     pause
% end
% 
% %% Movie
% 
% frames = LPPLmovie(data_prices, data_dates, 50);
% 
% %%
% movie(frames, 10, 2)
% 


