
% B-spline basis function value B(j,n) at x.
%
% Input arguments:
% j:
%    interval index, 0 =< j < numel(t)-n
% n:
%    B-spline order (2 for linear, 3 for quadratic, etc.)
% t:
%    knot vector
% x (optional):
%    value where the basis function is to be evaluated
%
% Output arguments:
% y:
%    B-spline basis function value, nonzero for a knot span of n


% Illustrates B-spline basis functions.
clear; close all; clc

% t = [0 0 0 0.5 1 1 1];  % knot vector
t = [0 0 0 0 0.3 .6 .9 1 1 1 1];  % knot vector
n = 4;  
x = 0:0.01:1;
x2 = 0.005:0.02:1;
x = [x x2];
figure('Name', sprintf('NURBS basis functions of order %d', n));
hold on;

for j = 0 : numel(t)-n-1
    [y,x] = bspline_basis(j,n,t,x);
    plot(x, y, '.');
end
