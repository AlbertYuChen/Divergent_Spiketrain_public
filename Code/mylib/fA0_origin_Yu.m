function [ A1 ] = fA0_origin_Yu( filter, ratebias, A0 )

t=1:length(filter);
% lambda = exp(ratebias+filter(t) + A0*sum(exp(filter(t:end))-1) );


for it = t
    lambda(it) = exp(ratebias+filter(it) + A0*sum( exp(filter(it+1:end))-1 ) );
end


lambda = [lambda, exp(ratebias*ones(1,2000))];

% S0=[];
% for t=1:length(lamda)
%     S0 = [S0,exp(-sum(lamda(1:t)))];
%     if S0<=1e-3
%         break;
%     end
% end

cum_lambda = [0 lambda(1:end-1)];
% S0 = exp( -cumsum(lambda) );
S0 = exp( -cumsum(cum_lambda) );

P0 = S0.*lambda;
P0 = P0/sum(P0);

tao = 1:length(P0);
A1 = 1/sum(P0.*tao);


% if A1>=0.9;
%     A1 = 0.9;
% end
% if isnan(A1)
%     A1 = 0.9;
% end
