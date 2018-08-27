function [ A1 ] = fA0_Yu_irrg2( filter, ratebias, A0 )

t=(1:length(filter));
% log_lambda = ratebias + filter(t) + A0*sum(filter);
% log_lambda = ratebias + filter(t) + sum(A0.*filter);

for it = t
    log_lambda(it) = ratebias + filter(it) + A0*sum(filter(it+1:end));
end

lambda = exp(log_lambda);

lambda = [lambda, exp(ratebias*ones(1, 2000))];
% figure
% plot(lamda)

S0 = exp( cumsum( -lambda) );


P0 = S0.*lambda;
tao = 1:length(P0);
A1 = 1/sum(P0.*tao);

% if A1>=0.9
%     A1 = 0.9;
% end
% if isnan(A1)
%     A1 = 0.9;
% end
