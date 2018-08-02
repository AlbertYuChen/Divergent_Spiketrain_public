% Author: Yu Chen
% Date: 2018 @ CNBC CMU

function KS_Test(Tau)

[eCDF, zvals] = ecdf(Tau.Zscr); zvals(1)=0;
mCDF = 1-exp(-zvals); %Model CDF at z values.


figure('Position', [300, 100, 600, 500]);
plot(mCDF, eCDF, 'Linewidth', 3)
hold on	 

% 0.95% CI
plot([0 1], [0 1]+1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]) 
h_l = plot([0 1], [0 1]-1.36/sqrt(length(Tau.Zscr)), '--', 'Color', [.7 .7 .7]); 

xlabel('Model CDF', 'interpreter', 'latex') 
ylabel('Empirical CDF', 'interpreter', 'latex')
legend(h_l, {'95\% CI'}, 'Interpreter','latex', 'location', 'southeast');

axis([0 1 0 1]);
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');


















