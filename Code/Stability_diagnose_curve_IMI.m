% Author: Yu Chen
% Date: Dec 5th 2017 @ CNBC

close all; clear; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% Load Data and parameter settings
spike_rec_res = 0.001;
stim_rec_res = 0.001;

Number_of_filters = 5; % options are 2,4,5

load([project_workspace_path 'Output/Human_IMI1_' num2str(Number_of_filters) 'h_out.mat'])

%% 
Rate = [1e-3:1e-3:1e-1,(1e-1+1e-2):1e-2:1 1:0.1:3];
% Rate = [1e-3:1e-3:1e-1,(1e-1+1e-2):1e-2:1];
% Rate = 10.^[-4:0.1:0];
ISI =  1./Rate ;
ISI = unique(round(ISI));
Rate = 1./ISI;

q1 = 2;
q2 = 0.5;

E_L0 = [];
E_L1 = [];
E_L2 = [];

for isi = ISI
E_L0 = [E_L0 f_ISI0_IMI_Yu(mod_set, isi)]; 
end
E_Rate0 = 1./E_L0;  

log_lambda_upper_bound = max(mod_set.kdc) +  max(mod_set.ih) * mod_set.SPK_Num;
lambda_upper_bound = exp(log_lambda_upper_bound);


% ----------------------------------------------
figure('Position', [300, 300, 600, 500]);

q0 = loglog(Rate, E_Rate0, 'linewidth', 4);
hold on 
loglog(exp(-10:1:5), exp(-10:1:5), '--', 'Color', [.7 .7 .7]);
loglog(exp([-10 5]), exp([log_lambda_upper_bound log_lambda_upper_bound]), '--', 'Color', [.9 .5 .5]);
leg1 = plot(exp([log_lambda_upper_bound log_lambda_upper_bound]), exp([-10 5]), '--', 'Color', [.9 .5 .5]);
axis(exp([-7 2 -7 2]));



xlabel('$A_0$ [spike/ms]', 'Interpreter', 'latex')
ylabel('$\mathcal{L}_h(A_0)$ [spike/ms]', 'Interpreter', 'latex') 
legend([leg1], {'Upper bound'}, 'Interpreter', 'latex', 'location', 'northwest')
title(['$\mathrm{FNF_{S}}$, k = ' num2str(mod_set.SPK_Num)], 'interpreter', 'latex')
set(gca, 'FontSize',18, 'xtick', [0.001, 0.01 0.1 1 10], 'ytick', [0.001, 0.01 0.1 1 10], 'TickLabelInterpreter','latex');







