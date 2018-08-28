% Date: Jan 12th 2018
% Author: Yu Chen @ MI

function [ A1 ] = fA0_Yu2( filter, ratebias, A0 )
if size(filter,1)>size(filter,2)
    filter = filter';
end

t=(1:length(filter));

% without losing generality, assume last spike is at t=0. 
for it = t
    log_lambda(it) = ratebias + filter(it) + A0*sum(filter(it+1:end));
end

lambda = exp(log_lambda);

lambda = [lambda, exp(ratebias*ones(1, 2000))];

cum_lambda = [0 lambda(1:end-1)];
S0 = exp( -cumsum(cum_lambda) );
PDF = S0.*lambda;
PDF = PDF/sum(PDF);
tao = 1:length(PDF);
A1 = 1/sum(PDF.*tao);














