% Author: Yu Chen
% Date: 2017 @ CNBC CMU

close all; clear; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% Load Data 
dataset = 'C';
stim_len = 1*1000; % 1000 ms

%---------------- Gerhard Truccolo dataset -------------------------
load([project_workspace_path 'Data/GerhardDegerTruccolo_2017_Data/Data_Fig2' dataset '.mat'])
% load([project_workspace_path 'GerhardDegerTruccolo_2017_Data/Data_Fig2E_Inhomo_Increase.mat'])

%-------------- convert this to my data type ------------
MC(1).spikeTimes = raster_x' / 0.001;
MC(1).spikeIndices = raster_y;

%% package data
spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_time = 0.001:spike_rec_res:stim_len*stim_rec_res;

valid_point_range = 1:stim_len;
KS_valid_point_range = valid_point_range;

Trial_Index = unique(MC(1).spikeIndices)';
N_trails = length(Trial_Index);

all_spike_train = [];
for ii = Trial_Index
    spike_index_array = MC(1).spikeIndices == ii;
    all_spike_train = [all_spike_train; hist(MC(1).spikeTimes(spike_index_array)*spike_rec_res, stim_time)];
end

cat_train = reshape(all_spike_train', 1, []);

% mean firing rate spk/ms
% mean(cat_train)
%% visualization
%-------------- Raster plot -----------------
figure('Position', [300, 300, 800, 500]);
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1), '.', 'Color', [0 0 0])
hold on
end
ylim([0.5 max(Trial_Index)+0.5])
title('raster plot of all spike trains')
xlabel('time [s]')
grid on

%% ========================= Setup fitting params %=========================
% this section uses piece of Jonathan Pillow's code to create basis
nhbasis = 10; 
hpeakFinal = 0.8; 

gg0 = makeFittingStruct_GLM(stim_rec_res, spike_rec_res, 10, 2, zeros(10, 1), nhbasis, hpeakFinal);
gg0.dtSp = spike_rec_res;
gg0.nlfun = @expfun; % default nonlinearity: exponential
gg0.sps = cat_train';  % Insert binned spike train into fitting struct

% post spike filter 
% figure('Position', [300, 300, 1400, 400]);
% plot(gg0.iht, gg0.ihbas); 

%% construct design matrix
% Determine # of parameters for self-coupling filter
nh = size(gg0.ihbas, 2);
% -------------- Create Spike-history Design Matrix -------------
Xsp = [];

for ti = 1:N_trails  
    trl_ind = stim_len*(ti-1) + (1:stim_len);
    Xsp_irt = hist_conv(gg0.sps(trl_ind), gg0.ihbas);
    Xsp = [Xsp; Xsp_irt(valid_point_range, :)];
end

%%  GLM fitting 
irt_ind_cum = [];
for ti = 1:N_trails  
    irt_ind_cum = [irt_ind_cum, (valid_point_range+stim_len*(ti-1))];
end

% design matrix
X_design = [ones(length(valid_point_range)*N_trails, 1), Xsp];

[prsML,dev,stats] = glmfit(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');

%% synthesize output
nk = 0;
%------------------------- DC + KDC ----------------------
gg0.dc = prsML(nk+1);
gg0.kdc = ones(stim_len,1)*gg0.dc;

%---------------------------- H ---------------------------
gg0.ihw = prsML(nk+2:nk+1+nh);
gg0.ih = gg0.ihbas*gg0.ihw;
gg0.nh = nh;

%---------
figure('Position', [300, 300, 400, 250]);
plot(gg0.iht, gg0.ih);
hold on

xlabel('time [s]')
title('post-spike filter');
axis tight;
legend('h', 'location', 'southeast');
ylim([-0.4 0.2])
grid on

% --------------------------- Misc -----------------------
gg0.valid_point_range = valid_point_range;
gg0.stim_time = stim_time;
gg0.Trial_Index = Trial_Index;

%% creat estimated firing rate lambda
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

% save([project_workspace_path 'Output/Fig2C_out_divergent_' date '.mat'], 'mod_set')










