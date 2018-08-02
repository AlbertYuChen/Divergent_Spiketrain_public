% Author: Yu Chen
% Date: 2017 @ CNBC CMU

close all; clear; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% Load Data and parameter settings
dataset = 'C';
stim_len = 1*1000;

%---------------- Gerhard Truccolo dataset -------------------------
load([project_workspace_path 'Data/GerhardDegerTruccolo_2017_Data/Data_Fig2' dataset '.mat'])

%-------------- convert this to my data type ------------
% MC(1).spikeTimes = raster_x' / 0.001;
MC(1).spikeTimes = round(raster_x' / 0.001);
MC(1).spikeIndices = raster_y;

%% package data
spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_time = 0.001:spike_rec_res:stim_len*stim_rec_res;

valid_point_range = 1:stim_len;
KS_valid_point_range = valid_point_range;

Trial_Index = 1:10;
N_trails = length(Trial_Index);

all_spike_train = [];
ISIs = [];

for ii = Trial_Index
    spike_index_array = MC(1).spikeIndices == ii;
    all_spike_train = [all_spike_train; hist(MC(1).spikeTimes(spike_index_array)*spike_rec_res, stim_time)];
    single_trail_spike_time = MC(1).spikeTimes(spike_index_array) * spike_rec_res;
    ISIs = [ISIs diff(single_trail_spike_time)]; 
end

cat_train = reshape(all_spike_train', 1, []);
%% visualization
%-------------- Raster plot -----------------
figure('Position', [300, 300, 800, 500]);
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2)
hold on
end
ylim([0.5 max(Trial_Index)+0.5])
grid on
%% =================== Setup fitting params ==========================
dtStim = stim_rec_res; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = spike_rec_res;

%% ===================== H =====================
% Monkey-PMv
nhbasis = 8;
hpeakFinal = 0.3; 
SPK_Num = 5; % number of spike to incorporate

gg0 = makeFittingStruct_GLM(dtStim, dtSp, 150, 0, reshape(zeros(150, 1), 150, []), nhbasis, hpeakFinal);
gg0.dtSp = spike_rec_res;
gg0.nlfun = @expfun; % default nonlinearity: exponential

gg0.irect = ones(SPK_Num, 1);
gg0.sps = cat_train';  % Insert binned spike train into fitting struct

% figure('Position', [300, 300, 1400, 400]);
% plot(gg0.iht, gg0.ihbas); 
nh = size(gg0.ihbas, 2);
Hbasis = [gg0.ihbas; zeros(15000, nh)];
%% construct design matrix
% ---- Create Spike-history Design Matrix ------------------------
Xsp = []; 
HF_u_accm = mat2cell( zeros(SPK_Num,0), ones(SPK_Num, 1) );
histH = {};

for ti = Trial_Index
    spike_index_array = MC(1).spikeIndices == ti;
    single_trail_spike_time = MC(1).spikeTimes(spike_index_array) * spike_rec_res; 
    
    Xsp_trial = 0;
%     figure; hold on
    for hi = 1:SPK_Num
        HF_u = get_IMI_input_time(single_trail_spike_time, stim_time, hi);
        ih_ind = round(HF_u/stim_rec_res);
        Xsp_irt = Hbasis(ih_ind, :);
        
        %------- not every spike start from 0 --------
        zero_area = round(single_trail_spike_time(hi)/stim_rec_res);
        Xsp_irt(1:zero_area, :) = 0;
        HF_u(1:zero_area) = NaN;
%         plot(HF_u)
        HF_u_accm{hi} = [HF_u_accm{hi}, HF_u(valid_point_range)];
        %---------------------------------------------
        Xsp_trial = Xsp_trial + Xsp_irt;       
%         plot(Xsp_irt)
    end
    Xsp = [Xsp; Xsp_trial(valid_point_range, :)];
end

% distribution of hazard function's input
HIST_bin_width = 0.005; % Monkey
figure('Position', [300, 100, 600, 500]); hold on
for hi = 1:SPK_Num
    histH{hi} = histogram(HF_u_accm{hi}, 'Normalization', 'pdf', 'BinWidth',HIST_bin_width, 'FaceAlpha', .7); 
end

leg1 = legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$','$u_6$');
set(leg1,'Interpreter','latex');
ylabel('PDF', 'interpreter', 'latex')
xlabel('time [sec]', 'interpreter', 'latex')
title('\textit{Monkey-PMv}', 'interpreter', 'latex')
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex');

%% ===================== Do ML fitting =====================
irt_ind_cum = [];
for ti = 1:N_trails  
    irt_ind_cum = [irt_ind_cum, (valid_point_range+stim_len*(ti-1))];
end

nk = 0;
b_k_ind = 1:nk;
b_dc_ind = nk+1;
b_h_ind = nk+2:nk+1+nh;

X_design = [ones(length(valid_point_range)*N_trails, 1), Xsp];

%--------- GLM ----------
% [prsML,dev,stats,loglikeli,iter,b_stack] = glmfit_YC(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');
[prsML,dev,stats] = glmfit(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');

%% synthesis output
% ---------------------------- K -------------------------
gg0.kt = prsML(1:nk);
gg0.k = gg0.ktbas*gg0.kt;
gg0.nk = nk;
%------------------------- DC + KDC ----------------------
gg0.dc = prsML(nk+1);
gg0.kdc = ones(stim_len,1)*gg0.dc;

%---------------------------- H ---------------------------
gg0.ihw = prsML(nk+1+1:nk+1+nh);
gg0.ih = gg0.ihbas*gg0.ihw;
gg0.nh = nh;
gg0.SPK_Num = SPK_Num;

%--------- for JNC paper plot ----------
figure('Position', [300, 300, 600, 500]);
plot(gg0.iht, gg0.ih, 'linewidth', 3);
hold on

xlabel('time [s]', 'interpreter', 'latex')
ylabel('log firing rate', 'interpreter', 'latex')
title('Hazard function', 'interpreter', 'latex');
title('\textit{Monkey-PMv} FNSF filter', 'interpreter', 'latex')
ylim([-1 0.5])
% xlim([0 1])
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');
grid on

% --------------------------- Misc -----------------------
gg0.valid_point_range = valid_point_range;
gg0.stim_time = stim_time;
gg0.Trial_Index = Trial_Index;

%%  GLM fitting 
lambda = glmval(prsML, X_design, 'log', 'constant', 'off');

%% goodness of fit test
Tau.Zscr = [];
Tau.p_k_i = {};
Tau.nspk = [];

for ti = 1:N_trails
    tRng = KS_valid_point_range + stim_len*(ti-1);
    spInd = find(gg0.sps(tRng)) + length(valid_point_range)*(ti-1);
    
    Tau.nspk = [Tau.nspk length(spInd)];
    for ii=2:length(spInd)
        Zscr_curr = sum(lambda((spInd(ii-1)+1):spInd(ii))); 
        Tau.p_k_i{end+1} = lambda( (spInd(ii-1)+1):spInd(ii) );
        Tau.Zscr = [Tau.Zscr Zscr_curr];

    end
end

%% append result
gg0.Tau = Tau;
mod_set = gg0;

%% Compute empirical CDF from rescaled waiting times.
KS_Test(Tau)
% KS_Test_Haslinger(Tau)
%% save output 
% save([project_workspace_path 'Output/Monkey_IMI1_5h_out_' date '.mat'], 'mod_set'); 








