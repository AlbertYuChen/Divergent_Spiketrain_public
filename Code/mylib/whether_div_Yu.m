function [ result ] = whether_div_Yu( filter, ratebias )
% return 2 for divergent, 1 for bi-stable, 0 for stable
if size(filter,1)>size(filter,2)
    filter = filter';
end

A1 = [];
A0 = [0:1e-3:1e-1,(1e-1+1e-2):1e-2:1];
for a0 = A0
    temp = fA0_Yu2(filter, ratebias, a0);
	A1 = [A1,temp];
end

% figure
% loglog(A0, A1, 'linewidth', 4);
% hold on
% plot([1e-4 1e0], [1e-4 1e0], '--', 'Color', [.7 .7 .7]);


B1 = A1-A0;
i = 1:(length(A0)-1);
bianhao = B1(i).*B1(i+1);
bianhaodian = find(bianhao<=0);
% fixedpoint = bianhaodian(1:2:length(bianhaodian));
fixedpoint = bianhaodian;

threshold = 1/3;

if isempty(fixedpoint)
    result = 2;
else
    
    if A1(fixedpoint(end))>=threshold
        result = 1;
    else 
        result = 0;
    end
    if A1(fixedpoint(1))>=threshold
        result = result + 1;
    end

end

