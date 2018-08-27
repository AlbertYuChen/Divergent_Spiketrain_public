% Author: Yu Chen
% Date: Jul 25th 2018 @ Beijing


function out = sigmoid(x)
% out = 1 ./ (1 + exp(-x));

dataClass = class(x);
lowerBnd = log(eps(dataClass)); 
upperBnd = -lowerBnd;
out = 1 ./ (1 + exp(-constrain(x,lowerBnd,upperBnd)));
end

function x = constrain(x,lower,upper)
% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;
end