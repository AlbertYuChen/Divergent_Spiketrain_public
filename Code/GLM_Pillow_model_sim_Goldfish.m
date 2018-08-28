% Author: Yu Chen
% Date: June 4th 2017 @ CNBC

close all; clear; clc
rng('default')
%% check working environment 
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% ========================= Raw ==========================
stim_len = 30*1000;

load([project_workspace_path 'Data/Goldfish.mat'])

MC(1).spikeTimes = raster_x * 1000;
MC(1).spikeIndices = raster_y;

spike_rec_res = 0.001;
stim_rec_res = 0.001;

stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;
PSTH_win = [0 stim_len*stim_rec_res];
% PSTH_win = [0 2];

PSHT_time_bin_width = 0.05;
PSHT_time_bin = PSHT_time_bin_width:PSHT_time_bin_width:stim_len*stim_rec_res;
Gauss_width = 0;

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

% ---------- ISI ------------
figure
histogram(ISIs, 30, 'Normalization', 'pdf');	%...compute histogram of ISIs,
title('Raw')
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');

%% ========================= MOD IMI ==========================
% load([project_workspace_path 'Output/Goldfish_IMI3_6h_out.mat'])
load([project_workspace_path 'Output/Goldfish_IMI3_5h_out.mat'])
% load([project_workspace_path 'Output/Goldfish_IMI1_h_4_out.mat'])

MC2.spikeTimes = [];
MC2.spikeIndices = [];
pseudo_all_spike_train = [];
pseudo_all_spike_train_time = {};
pseudo_all_splike_train_ISI = {};
pseudo_Var_ISI = [];
pseudo_Mean_ISI = [];
pseudo_CV_ISI = [];

for sim_ind = 1:1
[tsp, sps, Itot, Istm] = simGLM_IMI_SpDi(mod_set, stim_len, 2*sim_ind+100);  % OLD: sim_ind*2+100
pseudo_all_spike_train = [pseudo_all_spike_train; sps'];
pseudo_all_spike_train_time{end + 1} = tsp';
pseudo_all_splike_train_ISI{end + 1} = diff(tsp');
pseudo_Var_ISI = [pseudo_Var_ISI std(pseudo_all_splike_train_ISI{end})];
pseudo_Mean_ISI = [pseudo_Mean_ISI mean(pseudo_all_splike_train_ISI{end})];
pseudo_CV_ISI = [pseudo_CV_ISI pseudo_Var_ISI(end)/pseudo_Mean_ISI(end)];
%------------------ save pseudo data --------------------
MC2.spikeTimes = [MC2.spikeTimes, tsp'];
MC2.spikeIndices = [MC2.spikeIndices, sim_ind*ones(1, length(tsp))];
end

% KS_Test(mod_set.Tau)
%---------
pseudo_all_spike_trainIMI = pseudo_all_spike_train;
PSTH_reshapeIMI = plot_PSTH(pseudo_all_spike_train, PSHT_time_bin, stim_time, zeros(stim_len, 0), PSTH_win, Gauss_width);
mod_setIMI = mod_set;
IstmIMI = Istm;
ItotIMI = Itot;
 
ISIs = cell2mat(pseudo_all_splike_train_ISI)*spike_rec_res;
figure
histogram(ISIs, 30, 'Normalization', 'pdf');	%...compute histogram of ISIs,
title('IMI')
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');


%% ===================== plot results =====================
PSTH_reshape1 = PSTH_reshapeRaw;

PSTH_reshape2 = PSTH_reshapeIMI;
pseudo_all_spike_train = pseudo_all_spike_trainIMI;
%----------------------- analysis pseudo data -----------------------
figure('Position', [300, 100, 1400, 500]);
subplot(311);
for ii = Trail_Index
plot(find(all_spike_train(ii,:))*spike_rec_res, ii*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'markersize', 2)
hold on
end
xlim(PSTH_win)
ylabel('trial #')
title('Spike raster of real data');
grid on; set(gca,'FontSize',12)

subplot(312);
for ii = Trail_Index
plot(find(pseudo_all_spike_train(ii,:))*spike_rec_res, ii*ones(sum(pseudo_all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'markersize', 2)
hold on
end

xlim(PSTH_win)
ylabel('trial #')
title('Spike raster of pseudo data');
grid on; set(gca,'FontSize',12)

subplot(313); 
plot(PSHT_time_bin, PSTH_reshape1 ,'Color',[0 0.5 0], 'LineWidth', 1)
hold on
plot(PSHT_time_bin, PSTH_reshape2, 'Color',[0.5 0 0], 'LineWidth', 1)
xlim(PSTH_win)
xlabel('time [s]')
ylabel('firing rate [spikes/ms]')
legend('True', 'Model')
grid on; set(gca, 'FontSize',12)

























