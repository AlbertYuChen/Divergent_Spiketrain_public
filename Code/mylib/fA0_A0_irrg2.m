function A1 = fA0_A0_irrg2( filter, ratebias )
% return 2 for divergent, 1 for bi-stable, 0 for stable
if size(filter,1)>size(filter,2)
    filter = filter';
end

q1 = 1.4;
q2 = 1/q1;

% fA0_Yu(filter, ratebias, A0)
% fA0_Yu_irrg(filter, ratebias, A0_arry)

A1 = [];
% Aq1 = [];
% Aq2 = [];

A0 = [0:1e-4:1e-1,(1e-1+1e-2):1e-2:1];

for a0 = A0
	A1 = [A1, fA0_Yu2(filter, ratebias, a0)];
    
%     A0_arry = Calculate_Irregular_A0_arry(filter, A0, q1);
%     Aq1 = [Aq1, fA0_Yu_irrg2(filter, ratebias, A0_arry)];
%     
%     A0_arry = Calculate_Irregular_A0_arry(filter, A0, q2);
%     Aq2 = [Aq2, fA0_Yu_irrg2(filter, ratebias, A0_arry)];
end



figure('Position', [300, 300, 600, 500]);
% q0 = plot(log(A0),log(A1), 'linewidth', 4);
% hold on
% q1 = plot(log(A0),log(Aq1), '--', 'linewidth', 2, 'color', [0    0.4470    0.7410]);
% q2 = plot(log(A0),log(Aq2), ':', 'linewidth', 2, 'color', [0    0.4470    0.7410]);
% axis([-7 0 -7 0]);

q0 = loglog(A0, A1, 'linewidth', 4);
hold on
% q1 = loglog(A0, Aq1, '--', 'linewidth', 2, 'color', [0    0.4470    0.7410]);
% q2 = loglog(A0, Aq2, ':', 'linewidth', 2, 'color', [0    0.4470    0.7410]);
plot([1e-5 1e1], [1e-5 1e1], '--', 'Color', [.7 .7 .7]);
axis([1e-4 1e0 1e-4 1e0]);


xlabel('$A_0$ [spike/ms]', 'Interpreter', 'latex')
ylabel('$\mathcal{L}_h(A_0)$ [spike/ms]', 'Interpreter', 'latex')
% title('Divergent', 'interpreter', 'latex') % Stable  Divergent Fragile
% title('\textit{Human-Cortex} FLF', 'interpreter', 'latex')  
title('\textit{Monkey-PMv} FLF', 'interpreter', 'latex')  

% legend([q0 q1 q2], {'q = 1', 'q = 1.4', 'q = 0.71'}, 'Interpreter', 'latex', 'location', 'northwest')
set(gca, 'FontSize',18, 'xtick', [0.001, 0.01 0.1 1], 'TickLabelInterpreter','latex');

% detect the change point 
% A1 is the output of fA0, the fixed point is the point where f(A0) - A0 
% changes the sign. 
% B1 = A1-A0;
% i = 1:(length(A0)-1);
% bianhao = B1(i).*B1(i+1);
% bianhaodian = find(bianhao<=0);
% fixedpoint = bianhaodian(1:2:length(bianhaodian));
% if A1(fixedpoint(end))>=threshold
%     result = 1;
% else 
%     result = 0;
% end
% if A1(fixedpoint(1))>=threshold
% 	result = result + 1;
% end

end






















