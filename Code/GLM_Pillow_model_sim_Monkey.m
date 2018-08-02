% Author: Yu Chen
% Date: 2017 @ CNBC CMU

close all; clear; clc

rng('default')
%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% ========================= Raw ==========================
dataset = 'C';
stim_len = 3*1000;

load([project_workspace_path 'Data/GerhardDegerTruccolo_2017_Data/Data_Fig2' dataset '.mat'])

MC(1).spikeTimes = raster_x * 1000;
MC(1).spikeIndices = raster_y;

spike_rec_res = 0.001;
stim_rec_res = 0.001;

stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;
PSTH_win = [0 stim_len*stim_rec_res];

PSHT_time_bin_width = 0.05;
PSHT_time_bin = PSHT_time_bin_width:PSHT_time_bin_width:stim_len*stim_rec_res;
Gauss_width = 0; % should be greater equal than 2. or set 0 to disable the function

neuron_index = 1;
Trail_Index = 1:max(MC(neuron_index).spikeIndices);
N_trails = length(Trail_Index);

all_spike_train = [];
cat_trail_spike_time_stairs = [];
ISIs = [];

for ii = Trail_Index
    spike_index_array = MC(neuron_index).spikeIndices == ii;
    all_spike_train = [all_spike_train; hist(MC(neuron_index).spikeTimes(spike_index_array) * spike_rec_res, stim_time)];
    single_trail_spike_time = MC(neuron_index).spikeTimes(spike_index_array) * spike_rec_res;
    ISIs = [ISIs; diff(single_trail_spike_time)]; 
end

cat_train = reshape(all_spike_train', 1, []);
spikeindex = find(cat_train);
PSTH_reshapeRaw = plot_PSTH(all_spike_train, PSHT_time_bin, stim_time, zeros(stim_len, 0), PSTH_win, Gauss_width);

%% ========================= MOD homogeneous ==========================
% load([project_workspace_path 'Output/Fig2C_out_divergent.mat'])
load([project_workspace_path 'Output/Fig2C_out_divergent_01-Jun-2018.mat'])

N_sim = N_trails;

pseudo_all_spike_train = [];
pseudo_all_spike_train_time = {};
pseudo_all_splike_train_ISI = {};
pseudo_Var_ISI = [];
pseudo_Mean_ISI = [];
pseudo_CV_ISI = [];

for sim_ind = 1:N_sim
[tsp, sps, Itot, Istm] = simGLM_SpDi(mod_set, stim_len, sim_ind + 100);  % run model 200
pseudo_all_spike_train = [pseudo_all_spike_train; sps'];
pseudo_all_spike_train_time{end + 1} = tsp';
pseudo_all_splike_train_ISI{end + 1} = diff(tsp');
pseudo_Var_ISI = [pseudo_Var_ISI std(pseudo_all_splike_train_ISI{end})];
pseudo_Mean_ISI = [pseudo_Mean_ISI mean(pseudo_all_splike_train_ISI{end})];
pseudo_CV_ISI = [pseudo_CV_ISI pseudo_Var_ISI(end)/pseudo_Mean_ISI(end)];
end

%---------
pseudo_all_spike_trainHomo = pseudo_all_spike_train;
PSTH_reshapeHomo = plot_PSTH(pseudo_all_spike_train, PSHT_time_bin, stim_time, zeros(stim_len, 0), PSTH_win, Gauss_width);
mod_setHomo = mod_set;

%% ========================= MOD Inhomo Pillow ==========================
load([project_workspace_path 'Output/Fig2C_INTENSITY_out.mat'])
% load([project_workspace_path 'Output/Fig2C_INTENSITY_out_01-Jun-2018.mat'])

figure('Position', [300, 300, 1400, 400]);
plot(mod_set.iht, mod_set.ihbas); 
title('basis for post-spike filter h')

N_sim = N_trails;

pseudo_all_spike_train = [];
pseudo_all_spike_train_time = {};
pseudo_all_splike_train_ISI = {};
pseudo_Var_ISI = [];
pseudo_Mean_ISI = [];
pseudo_CV_ISI = [];

for sim_ind = 1:N_sim
[tsp, sps, Itot, Istm] = simGLM_INTENSITY_SpDi(mod_set, stim_len, sim_ind+15);  
pseudo_all_spike_train = [pseudo_all_spike_train; sps'];
pseudo_all_spike_train_time{end + 1} = tsp';
pseudo_all_splike_train_ISI{end + 1} = diff(tsp');
pseudo_Var_ISI = [pseudo_Var_ISI std(pseudo_all_splike_train_ISI{end})];
pseudo_Mean_ISI = [pseudo_Mean_ISI mean(pseudo_all_splike_train_ISI{end})];
pseudo_CV_ISI = [pseudo_CV_ISI pseudo_Var_ISI(end)/pseudo_Mean_ISI(end)];
end

%---------
pseudo_all_spike_trainInHomo = pseudo_all_spike_train;
PSTH_reshapeInHomo = plot_PSTH(pseudo_all_spike_train, PSHT_time_bin, stim_time, zeros(stim_len, 0), PSTH_win, Gauss_width);
mod_setInHomo = mod_set;


%% ===================== plot results =====================
PSTH_reshape1 = PSTH_reshapeRaw;

PSTH_reshape2 = PSTH_reshapeHomo;
pseudo_all_spike_train = pseudo_all_spike_trainHomo;

% PSTH_reshape2 = PSTH_reshapeInHomo;
% pseudo_all_spike_train = pseudo_all_spike_trainInHomo;

% PSTH_reshape2 = PSTH_reshapeIMI;
% pseudo_all_spike_train = pseudo_all_spike_trainIMI;

%% Raster + PSTH 
figure('Position', [300, 100, 600, 500]);
subplot(211); 
% SPT = all_spike_train;
% PSTH = PSTH_reshapeRaw;

SPT = pseudo_all_spike_train;
PSTH = PSTH_reshape2;

for ii = Trail_Index
plot(find(SPT(ii,:))*spike_rec_res, ii*ones(sum(SPT(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'markersize', 2)
hold on
end
ylim([0.2 size(SPT,1)+0.8])
xlim(PSTH_win)
ylabel('Trial \#', 'interpreter', 'latex')
title('\textit{Monkey-PMv} simulation', 'interpreter', 'latex')

set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');

subplot(212); 
plot(PSHT_time_bin, PSTH, 'LineWidth', 3)
hold on
plot(PSHT_time_bin(1:20), PSTH_reshapeRaw(1:20), 'LineWidth', 3)
legend({'Simulation', 'Real'}, 'Interpreter','latex', 'location', 'northwest');

xlim(PSTH_win)
ylim([0 0.2])
xlabel('time [s]', 'interpreter', 'latex')
ylabel('firing rate [spk/ms]' , 'interpreter', 'latex')
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex')

%% PSTH together
figure('Position', [300, 100, 600, 500]);

plot(PSHT_time_bin, PSTH, 'LineWidth', 3)
hold on
plot(PSHT_time_bin, PSTH_reshapeInHomo, 'LineWidth', 3)
plot(PSHT_time_bin(1:20), PSTH_reshapeRaw(1:20), 'LineWidth', 3, 'color', [0.4660    0.6740    0.1880])

legend({'Homogeneous FLF', 'Inhomogeneous FLF', 'Real'}, 'Interpreter','latex', 'location', 'northwest');

xlim(PSTH_win)
ylim([0 0.2])
xlabel('time [s]', 'interpreter', 'latex')
ylabel('firing rate [spike/ms]' , 'interpreter', 'latex')
title('\textit{Monkey-PMv} simulation', 'interpreter', 'latex')
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex')
























