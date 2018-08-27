function [ result, y ] = whether_div_real( filter, ratebias )
% return 2 for explosion, 0 for stable
T=1e4;
nHistBins = numel(filter);
y = zeros(1 , T + nHistBins);

for t = (nHistBins+1):(T+nHistBins)
    yy = poissrnd(exp(filter * (y(t - (1:nHistBins)))' + ratebias));
    if yy ~= 0
        y(t) = 1;
    end
end

y = y(:,nHistBins+1:end);

% Qi Xin's parameter
% if length(find(y(:,end-1000:end)))>=900

if length(find(y(:,end-1000:end)))>=299 % the threshold is is 0.9 of maximum firing rate. 
    result = 2;
else
    result = 0;
end

end

