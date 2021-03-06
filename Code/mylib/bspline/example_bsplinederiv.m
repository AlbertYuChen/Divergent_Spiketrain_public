function example_bsplinederiv
% Illustrates drawing a B-spline and its derivative.

% Copyright 2011 Levente Hunyadi

% spline order
k = 4;
% knot sequence
t = [0 0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1 1 1];
% control points
P = [ 0.1993 0.4965 0.6671 0.7085 0.6809 0.2938 0.1071 0.3929 0.5933 0.8099 0.8998 0.8906 ...
    ; 0.8377 0.8436 0.7617 0.6126 0.212 0.1067 0.3962 0.5249 0.5015 0.3991 0.6477 0.8553 ];

C = bspline_deboor(k,t,P);
[dt,dP] = bspline_deriv(k,t,P);
dC = bspline_deboor(k-1,dt,dP);

% plot control points and spline
figure;
hold all;
plot(C(1,:), C(2,:), 'b');
plot(P(1,:), P(2,:), 'kx');
plot(dC(1,:), dC(2,:), 'm');
plot(dP(1,:), dP(2,:), 'k.');
legend('original spline curve','control points','derivative spline curve','control points of derivative', ...
    'Location', 'Best');
hold off;