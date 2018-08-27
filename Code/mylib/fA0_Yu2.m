% Date: Jan 12th 2018
% Author: Yu Chen @ MI
% This code correct the tiny problem of 
% log_lambda(it) = ratebias + filter(it) + A0*sum(filter(it+1:end)); 
% in code fA0_Yu, Qi made a mistake, which misses the integral range
% ERROR: log_lambda(it) = ratebias + filter(it) + A0*sum(filter); 

function [ A1 ] = fA0_Yu2( filter, ratebias, A0 )
if size(filter,1)>size(filter,2)
    filter = filter';
end

t=(1:length(filter));
% t=(1:3)';
% log_lambda = ratebias + filter(t) + A0*sum(filter);
% Xin Qi's method. WRONG !!!
% Should take care of the integral term. 


for it = t
    log_lambda(it) = ratebias + filter(it) + A0*sum(filter(it+1:end));
end

% filter0 = filter;
% filter0(1:2) = -1e12;
% for it = t
%     log_lambda0(it) = ratebias + filter0(it) + A0*sum(filter0(it+1:end));
% end

% % itr = 1:3;
% figure
% plot(log_lambda0, '+')
% hold on
% plot(log_lambda, 's')
% % ylim([-1000 10])

lambda = exp(log_lambda);

% figure
% plot(lambda)

lambda = [lambda, exp(ratebias*ones(1, 2000))];
% figure
% plot(lambda)

cum_lambda = [0 lambda(1:end-1)];
% S0 = exp( -cumsum(lambda) );
S0 = exp( -cumsum(cum_lambda) );
% CDF = 1 - S0;

% figure
% plot(CDF)

% P0 = S0.*lambda;
PDF = S0.*lambda;
PDF = PDF/sum(PDF);
% sum(PDF)

% figure
% plot(PDF)



tao = 1:length(PDF);
A1 = 1/sum(PDF.*tao);

% sum(PDF)
















