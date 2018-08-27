clear; close all; clc

%% check working environment
% I'm swtching between PC, Mac, and CNBC cluster.
project_workspace_path = Initialization_Env_Divergent;

%% Pillow basis
nBases=10;
Baseslength=400;
tstrech=5;
raw_bs = basisFactory.makeNonlinearRaisedCos(nBases, 1, [0 Baseslength], tstrech);
basis = [raw_bs.B(58:407,6)';raw_bs.B(58:407,9)'];
basis1 = basis(1,:);
basis2 = basis(2,:);

% figure
% plot(basis1, 'linewidth', 4); hold on
% plot(basis2, 'linewidth', 4);
% legend({'Basis 1', 'Basis 2'}, 'Interpreter','latex', 'location', 'southeast'); legend boxoff  
% set(gca, 'FontSize',28, 'xtick', [], 'ytick', [], 'TickLabelInterpreter','latex');

t = (1:length(basis1))/1000;

DPI=0.2;


negx=-5:DPI:5;
y=-5:DPI:5;

% stable
baseline = -4;
x = -4;
y = -1;
h_stable = x*basis1+y*basis2;

% diverge
baseline = -4;
x = -4;
y = 2;
h_divergent = x*basis1+y*basis2;

% fragile
baseline = -4;
x = -3;
y = 0.6; % 0.6
h_fragile = x*basis1+y*basis2;

% fragile
baseline = -4;
x = -4;
y = 7;
h_test = x*basis1+y*basis2;

h = h_fragile;
% --------- for JNC paper ----------
figure('Position', [300, 300, 600, 500]);
plot(t, h_fragile, 'linewidth', 3);
hold on
plot(t, h_divergent, 'linewidth', 3);
plot(t, h_stable, 'color', [0.4660    0.6740    0.1880], 'linewidth', 3);

xlabel('time [s]', 'interpreter', 'latex')
ylabel('log firing rate', 'interpreter', 'latex')
title('Hazard function', 'interpreter', 'latex');
% legend({'Fragile', 'Divergent', 'Stable'}, 'Interpreter','latex', 'location', 'southeast');
title('Synthetic filters', 'interpreter', 'latex')
ylim([-5 3])
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');
grid on

%% k is bias term, beta0
A_stable = [];
A_divergent = [];
A_fragile = [];

A0 = 10.^(-4:0.1:0);

figure('Position', [300, 300, 600, 500]);
plot(t, h, 'linewidth', 3);
hold on
ylim([-2 5])

for a0 = A0
	A_stable = [A_stable, fA0_Yu2(h_stable, baseline, a0)]; 
    A_divergent = [A_divergent, fA0_Yu2(h_divergent, baseline, a0)]; 
    A_fragile = [A_fragile, fA0_Yu2(h_fragile, baseline, a0)]; 
end


%% plot together
figure('Position', [300, 300, 600, 500]);
loglog(A0, A_fragile, 'linewidth', 4);
hold on
loglog(A0, A_divergent, 'linewidth', 4);
loglog(A0, A_stable, 'color', [0.4660    0.6740    0.1880], 'linewidth', 4);
plot(exp(-7:1:0), exp(-7:1:0), '--', 'Color', [.7 .7 .7]);
axis(exp([-7 0 -7 0]));


xlabel('$A_0$ [spike/ms]', 'Interpreter', 'latex')
ylabel('$\mathcal{L}_h(A_0)$ [spike/ms]', 'Interpreter', 'latex')
title('Diagnosis curves', 'interpreter', 'latex')  
legend({'Fragile', 'Divergent', 'Stable'}, 'Interpreter','latex', 'location', 'southeast'); legend boxoff
set(gca, 'FontSize',18, 'xtick', [0.001, 0.01 0.1 1], 'TickLabelInterpreter','latex');


















