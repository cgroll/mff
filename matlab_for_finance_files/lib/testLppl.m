%% Explaining the parts of the lppl function
% The lppl function consists of two parts: an hyperbolic function
% that increases to infinity at a singularity, and an oscillation of
% increasing frequency.

% set up realistic parameters
params = [50.9057703950095;
   -39.2228401917215;
   1348.19815338505;
   0.0142416695966602;
   -0.00177914571357913;
   -6.22638522525240;
   0.00974118509923082;
   0.915241664436276];

% plot hyperbolic part
grid = (1:0.2:(params(3)/params(8)-1));

% init anonymous function to get function values
evalLppl = @(x, params) ...
   params(1) + params(2)*(params(3)-params(8)*x).^params(4).*...
   (1+params(5)*cos(params(6).*...
   log(params(3)-params(8)*x)+params(7)));

vals = evalLppl(grid, params);

plot(grid, vals)
%% Hyperbolic part
% If we set the fifth parameter equal to zero, the oscillating part
% of the function gets removed.

paramsHyperb = params;
paramsHyperb(5) = 0;

valsHyperb = evalLppl(grid, paramsHyperb);
plot(grid, valsHyperb)

%% Cosinus part
% If we set parameters one and four to zero, and two and five to one,
% we get oscillating part up to rescaling

paramsOsc = params;
paramsOsc([1, 4]) = 0;
paramsOsc(2) = 1;
paramsOsc(5) = 1;

valsOs = evalLppl(grid, paramsOsc);
plot(grid, valsOs);

%% Both parts together
plot(grid, valsOs + valsHyperb)

%% variation of first parameter
variations = params(1)*[0.8 0.9 1.1 1.2];
plotVaryLpplParameters(params, 1, variations);

% The first parameter shifts the function additively up and down.
% Higher values increase the level of the function.

%% variation of second parameter
clf
variations = params(2)*[0.2 0.9 1.1 1.9];
plotVaryLpplParameters(params, 2, variations);

% The second parameter linearly scales the function up and down.
% Higher values increase the level of the function.

%% variation of third parameter
clf
variations = params(3)*[0.8 0.9 1.1 1.2];
plotVaryLpplParameters(params, 3, variations);

% The third parameter additively shifts the function left and right.
% Higher values shift the function to the right.

%% variation of forth parameter
clf
variations = params(4)*[0.2 0.9 1.1 1.5];
plotVaryLpplParameters(params, 4, variations);

% The forth parameter determines the shape of the function. Hence, it
% determines how long the period of curvature persists, and hence how
% long the period of super-exponential persists.

%% variation of fifth parameter
clf
variations = params(5)*[0.01 0.5 5 20];
plotVaryLpplParameters(params, 5, variations);

% The fifth parameter determines the amplitude of the cosinus. Large
% values increase the amplitude.

%% variation of sixth parameter
clf
variations = params(6)*[0.01 0.5 5 20];
plotVaryLpplParameters(params, 6, variations);

% The sixth parameter determines the frequency of the cosinus. Large
% values increase the frequency.

%% variation of seventh parameter
clf
variations = params(7)*[-1000  -200 200 1000];
plotVaryLpplParameters(params, 7, variations);

% The seventh parameter shifts the cosinus additively left and right.
% Positive values shift to the right.

%% variation of eighth parameter
clf
variations = params(8)*[0.8 0.9 1.1 1.2];
plotVaryLpplParameters(params, 8, variations);

