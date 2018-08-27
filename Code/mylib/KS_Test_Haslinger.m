% Author: Yu Chen
% Date: May 19th 2018 @ CNBC CMU

function Tau = KS_Test_Haslinger(Tau)

rng(623)
% rng('default')

zeta_list = [];
for nn = 1:length(Tau.p_k_i)
    
    p_k_i_list = Tau.p_k_i{nn};
    q_k_i_list = -log( 1-p_k_i_list );
    zeta = sum( q_k_i_list(1:(end-1)) );
    r = rand(1);
    delta = -log( 1-r*p_k_i_list(end) );
    zeta_list = [zeta_list, zeta+delta];
end

Tau.zeta_list = zeta_list;
[eCDF, zvals] = ecdf(zeta_list);
mCDF = 1-exp(-zvals); %Model CDF at z values.

figure('Position', [300, 100, 600, 500]);
plot(mCDF, eCDF, 'Linewidth', 3)
hold on	 

% 0.95% CI
plot([0 1], [0 1]+1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]) 
h_l = plot([0 1], [0 1]-1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]); 

xlabel('Model CDF', 'interpreter', 'latex') 
ylabel('Empirical CDF', 'interpreter', 'latex')
title('Haslinger KS test')
legend(h_l, {'95\% CI'}, 'Interpreter','latex', 'location', 'southeast');

axis([0 1 0 1]);
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');












