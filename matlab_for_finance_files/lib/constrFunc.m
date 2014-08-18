function [c ceq] = constrFunc(x)
% this function is determining the constraints in the portfolio
% optimization

% the function is based on a global variable covMatr and a global
% variable that is determining the portfolio standard deviation
% of interest.

global COVMATR GIVENSTDDEV

ceq = x*COVMATR*x' - GIVENSTDDEV^2;
c = [];