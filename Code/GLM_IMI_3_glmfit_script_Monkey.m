% Author: Yu Chen
% Date: Nov 1st 2017 @ CNBC

close all; clear; clc

%% check working environment
% I'm swtching between PC, Mac, and CNBC cluster.
project_workspace_path = Initialization_Env_Divergent;

%% Load Data and parameter settings
dataset = 'C';
stim_len = 1*1000;

%---------------- Truccolo -------------------------
load([project_workspace_path 'GerhardDegerTruccolo_2017_Data/Data_Fig2' dataset '.mat'])
% load([project_workspace_path 'GerhardDegerTruccolo_2017_Data/Data_Fig2E_Inhomo_Increase.mat'])
% load([project_workspace_path 'GerhardDegerTruccolo_2017_Data/Data_Fig2H_Trial2trialVar_13_7.mat'])
% load([project_workspace_path 'GerhardDegerTruccolo_2017_Data/Goldfish.mat'])
% load([project_workspace_path 'GerhardDegerTruccolo_2017_Data/Burst.mat'])

% figure('Position', [300, 100, 800, 700]);
% 
% trial_ind = raster_y == 1;
% plot(raster_x(trial_ind), ones(sum(trial_ind), 1), '.')
% xlim([0 2])

%-------------- convert this to my data type ------------
MC(1).spikeTimes = round(raster_x' / 0.001);
MC(1).spikeIndices = raster_y';

%% Load Data and parameter settings
neuron_index = 1;
saveOutput = 0; % 1: SAVE  0: NOT Save

spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;

mod_set = [];
Tau_All.Zscr = [];
Tau_All.Zi1 = [];
Tau_All.Zi2 = [];

Time_fold_index = {(1:stim_len)'};  
% Time_fold_index = {(200:stim_len)'};  
% Time_fold_index = {(100:stim_len)'};  % Monkey-PMv

for vpi = 1:length(Time_fold_index)
valid_point_range = Time_fold_index{vpi};
KS_valid_point_range = valid_point_range;

% Trial_fold_Index = {1};
Trial_fold_Index = {1:10};
for iii = 1:length(Trial_fold_Index)
Trial_Index = Trial_fold_Index{iii};
N_trails = length(Trial_Index);
Nfold = length(Trial_fold_Index);

all_spike_train = [];
ISIs = [];
for ii = Trial_Index
    spike_index_array = MC(neuron_index).spikeIndices == ii;
    all_spike_train = [all_spike_train; hist(MC(neuron_index).spikeTimes(spike_index_array)*spike_rec_res, stim_time)];
    single_trail_spike_time = MC(neuron_index).spikeTimes(spike_index_array) * spike_rec_res;
    ISIs = [ISIs diff(single_trail_spike_time)]; 
end

cat_train = reshape(all_spike_train', 1, []);
spikeindex = find(cat_train);


% ------------ Preview data ------------
figure('Position', [300, 300, 800, 500]);
for ii = 1:size(all_spike_train,1)
plot(find(all_spike_train(ii,:))*spike_rec_res, Trial_Index(ii)*ones(sum(all_spike_train(ii,:)),1),...
    's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'Markersize', 2)
hold on
end
ylim([0.5 max(Trial_Index)+0.5])

% ---------- ISI ------------
figure
histogram(ISIs, 30, 'Normalization', 'pdf');	%...compute histogram of ISIs,
title('ISI distribution')
%% ===================== Setup fitting params ======================
dtStim = stim_rec_res; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = spike_rec_res;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
gg0.dtStim = dtStim;
gg0.dtSp = dtSp;
gg0.nlfun = @expfun; % default nonlinearity: exponential
gg0.sps = cat_train';  % Insert binned spike train into fitting struct

%% ===================== K =====================
nkt = 150;    % Number of time bins in stimulus filter k
gg0.ttk = dtStim*(0:nkt-1)';  % time relative to spike of stim filter taps
nk = 0;  % number of basis vectors for representing k
[~,gg0.ktbas,~] = makeKbasis(dtStim, dtSp, nkt, nk);

% ---- Convolve stimulus with spatial and temporal bases -----
Xstim = zeros(length(valid_point_range)*N_trails, nk);
if nk > 0
Xstim_One = stim_conv(stim, gg0.ktbas);
Xstim = kron(ones(N_trails, 1), Xstim_One(valid_point_range,:));
end

%% ===================== H design =====================
% SPK_Num = 6; % Goldfish
% SPK_Num = 5; % Izhikevich-burst
% SPK_Num = 7; % Monkey-PMv
SPK_Num = 5; % Monkey-PMv JCN paper

gg0.SPK_Num = SPK_Num;
HF_u_accm = mat2cell( zeros(SPK_Num,0), ones(SPK_Num, 1) );
histH = {};

basis_xlim = [0 0.5]; % Monkey-PMv
% basis_xlim = [0 8]; % Human-cortex
% basis_xlim = [0 0.3]; % Izhikevich-burst

for ti = Trial_Index
    spike_index_array = MC(neuron_index).spikeIndices == ti;
    single_trail_spike_time = MC(neuron_index).spikeTimes(spike_index_array) * spike_rec_res;
    
%     figure; hold on
%     for hi = 1:SPK_Num
    for hi = 1:min(SPK_Num, length(single_trail_spike_time))
        HF_u = get_IMI_input_time(single_trail_spike_time, stim_time, hi);
        %------- not every spike start from 0 --------
        zero_area = round(single_trail_spike_time(hi)/stim_rec_res);
        HF_u(1:zero_area) = NaN;
%         plot(HF_u)    
        HF_u_tmp = HF_u(valid_point_range);
        HF_u_tmp(isnan(HF_u_tmp)) = [];    
        HF_u_accm{hi} = [HF_u_accm{hi}, HF_u_tmp];
    end        
end

% distribution of hazard function data
% HIST_bin_width = 0.002; % Goldfish
HIST_bin_width = 0.005; % Monkey-PMv
% HIST_bin_width = 0.1; % Human-Cortex
figure('Position', [300, 100, 600, 500]); 
% for hi = [1 2 3 7]
for hi = 1:SPK_Num
    
    if hi == 10
    histH{hi} = histogram(HF_u_accm{hi}, 'Normalization', 'pdf','Facecolor', [0.4660 0.6740 0.1880], 'BinWidth', HIST_bin_width, 'FaceAlpha', .7);  hold on    
%     sum(histH{hi}.Values * histH{hi}.BinWidth)
    continue
    elseif hi ==7 
    histH{hi} = histogram(HF_u_accm{hi}, 'Normalization', 'pdf','Facecolor', [.7 .7 .7], 'BinWidth', HIST_bin_width, 'FaceAlpha', .7);  hold on   
    continue
    end
    
    histH{hi} = histogram(HF_u_accm{hi}, 'Normalization', 'pdf', 'BinWidth', HIST_bin_width, 'FaceAlpha', .7);  hold on
end


% legend('u_1','u_2','u_3','u_4')
leg1 = legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$','$u_6$');
set(leg1,'Interpreter','latex');
ylabel('PDF', 'interpreter', 'latex')
xlabel('time [sec]', 'interpreter', 'latex')
% title('\textit{Izhikevich-burst}', 'interpreter', 'latex')
title('\textit{Monkey-PMv}', 'interpreter', 'latex')
% title('\textit{Monkey-PMv} simulation', 'interpreter', 'latex')
% title('\textit{Human-Cortex}', 'interpreter', 'latex')
% title('\textit{Human-Cortex} simulation', 'interpreter', 'latex')
xlim(basis_xlim)
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex');

% ----------- calculate the quatile of u_x ----------------
% hi = 7;
% quantile(HF_u_accm{hi}, 0.25)
% quantile(HF_u_accm{hi}, 0.75)

%% ======================= H =======================
% h_num_basis = [7 6 5 4 4 4 4]; % burst
% h_num_basis = [7 5 4]; % monkey
h_num_basis = [7 6 5 4 4 4 4 4 4 4]; % Human

gg0.iht = {};
gg0.ihbas = {};
gg0.nh = {};
gg0.knots_intense = {};

for hi = 1:SPK_Num
[BasisBegin, BasisEnd, NodeArry] = NodeArray_Distribution(histH{hi}, h_num_basis(hi) );
B_s = 0;
% B_e = BasisEnd + 0.02;

if hi == 1
B_e = BasisEnd + 0.03;
elseif hi == 2
   B_e = BasisEnd + 0.02; 
else
    B_e = BasisEnd;
end

spline_order = 4;  % 4 is cubic, 3 is quadratic
if hi == 1
knots_intense = [0 0 0 0 NodeArry B_e]; 
else
knots_intense = [BasisBegin BasisBegin NodeArry B_e B_e]; 
end
XBbas = [];
for spline_index = 0 : numel(knots_intense) - spline_order - 1
    [yy, x1] = bspline_basis(spline_index, spline_order, knots_intense, (B_s:stim_rec_res:B_e)'); 
    XBbas = [XBbas yy];
end

gg0.nh{hi} = size(XBbas, 2);
gg0.iht{hi} = (B_s:stim_rec_res:B_e)';
gg0.ihbas{hi} = XBbas;
gg0.ihbasExt{hi} = [XBbas; zeros(8000, gg0.nh{hi})];
gg0.knots_intense{hi} = knots_intense;

figure('Position', [200, 200, 1200, 300]);
plot(gg0.iht{hi}, gg0.ihbas{hi})
% xlim(basis_xlim)
title(['Hazard function ' num2str(hi) ' basis'])
xlabel('time [sec]')

end


%% ------------------------ Create Spike-history Design Matrix ------------------------
Xsp = mat2cell( zeros(SPK_Num, 0), ones(SPK_Num, 1) )';
for ti = Trial_Index
    
    spike_index_array = MC(neuron_index).spikeIndices == ti;
    single_trail_spike_time = MC(neuron_index).spikeTimes(spike_index_array) * spike_rec_res;

%     figure; hold on
    for hi = 1:SPK_Num
        HF_u = get_IMI_input_time(single_trail_spike_time, stim_time, hi);
        Sample_ind = round(HF_u/stim_rec_res);
        Xsp_irt = gg0.ihbasExt{hi}(Sample_ind, :);
        %------- not every spike start from 0 --------
        if hi <= length(single_trail_spike_time)
            zero_area = round(single_trail_spike_time(hi)/stim_rec_res);
            HF_u(1:zero_area) = NaN;
            Xsp_irt(1:zero_area, :) = 0;
        else
            Xsp_irt(:) = 0;
            HF_u(:) = NaN;
        end
        %---------------------------------------------
        
        Xsp{hi} = [Xsp{hi}; Xsp_irt(valid_point_range, :)];
%         figure
%         plot(Xsp{hi})
%         plot(HF_u); xlim([0 length(HF_u)])
    end
end

%% Design matrix
% create spike train index for all concatenated spike train
irt_ind_cum = [];
for ti = 1:N_trails  
    irt_ind_cum = [irt_ind_cum; (valid_point_range+ stim_len*(ti-1))];
end

X_design = [Xstim, ones(length(valid_point_range)*N_trails, 1), cell2mat(Xsp)];

%% ============= GLM ================
[prsML,dev,stats,loglikeli,iter,b_stack] = glmfit_YC(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');

%% ============================= DC =============================
gg0.dc = prsML(nk+1);
gg0.kdc = gg0.dc*ones(stim_len, 1);

%% ============================= H =============================
gg0.ihw = {};
gg0.ih ={};
gg0.h_SE = {};

nh_arry = cell2mat(gg0.nh);
nh_arry = [0 cumsum(nh_arry)];

figure('Position', [300, 100, 600, 500]); 
for hi = 1:SPK_Num
b_h_ind = nk+1+nh_arry(hi)+1:nk+1+nh_arry(hi+1);
gg0.ihw{hi} = prsML(b_h_ind);
gg0.ih{hi} = gg0.ihbas{hi}*gg0.ihw{hi};
gg0.h_SE{hi} = sqrt(diag(gg0.ihbas{hi}*stats.covb(b_h_ind, b_h_ind)*gg0.ihbas{hi}'));
%----

x = gg0.iht{hi};
y = gg0.ih{hi};
dy = gg0.h_SE{hi};
% fill([x;flipud(x)], [y-dy; flipud(y+dy)], [.9 .85 .85], 'linestyle','none', 'facealpha',.5);

if hi == 30
plot(x, y, '-', 'linewidth', 3, 'color', [0.4660 0.6740 0.1880])   
continue
end

plot(x, y, '-', 'linewidth', 3)
hold on
end

xlabel('time [sec]', 'interpreter', 'latex')
ylabel('log firing rate', 'interpreter', 'latex')
% legend('h_{1} SE','h_{1}','h_{2} SE','h_{2}','h_{3} SE','h_{3}','h_{4} SE','h_{4}', 'location', 'southeast');
legend({'$h_{1}$', '$h_{2}$', '$h_{3}$', '$h_{4}$', '$h_{5}$'}, 'interpreter', 'latex', 'location', 'southeast');
title('\textit{Monkey-PMv} FNMF filters', 'interpreter', 'latex')
% title('\textit{Human-Cortex} filters', 'interpreter', 'latex')
% title('\textit{Goldfish} IMI-3 filters', 'interpreter', 'latex')
% title('\textit{Izhikevich-burst} IMI-3 filters', 'interpreter', 'latex')
% xlim([0 0.25])
% ylim([-2 0.5])
% ylim([-4 2])
% ylim([-20 10])
grid on
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');
%% ============================= Misc =============================
gg0.valid_point_range = valid_point_range;
gg0.stim_time = stim_time;
gg0.neuron_index = neuron_index;
gg0.Trial_Index = Trial_Index;

%% KS analysis
lambda = glmval(prsML, X_design, 'log', 'constant', 'off');

% ------------- fitting statistics ------------
log_likelihood = -lambda + gg0.sps(irt_ind_cum) .* (X_design*prsML);
log_likelihood = sum(log_likelihood);
LL = sum(log(poisspdf(gg0.sps(irt_ind_cum), lambda)));
AIC = -2*log_likelihood + 2*length(prsML);
BIC = -2*log_likelihood + length(prsML)*log( length(lambda) );
disp(['AIC: ' num2str(AIC) '  BIC: ' num2str(BIC)])
disp(['N_parameter: ' num2str(length(prsML)) ])

gg0.AIC = AIC;
gg0.BIC = BIC;

%% goodness of fit test
Tau.Zscr = [];
Tau.Zi1 = [];
Tau.Zi2 = [];

for ti = 1:N_trails
    tRng = KS_valid_point_range + stim_len*(ti-1);
    spInd = find(gg0.sps(tRng)) + length(valid_point_range)*(ti-1);
    
    Zscr_last = NaN;
    for ii=2:length(spInd)
        Zscr_curr = sum(lambda((spInd(ii-1)+1):spInd(ii))); 
        Tau.Zscr = [Tau.Zscr Zscr_curr];

        %---- for iid test ----
        if ~isnan(Zscr_last)
        Tau.Zi1 = [Tau.Zi1 Zscr_last];
        Tau.Zi2 = [Tau.Zi2 Zscr_curr];
        end
        Zscr_last = Zscr_curr;
    end

end

%% append result
gg0.Itot = X_design * prsML;
gg0.lambda = lambda;
gg0.Tau = Tau;
mod_set = [mod_set gg0];
Tau_All.Zscr = [Tau_All.Zscr Tau.Zscr];
Tau_All.Zi1 = [Tau_All.Zi1 Tau.Zi1];
Tau_All.Zi2 = [Tau_All.Zi2 Tau.Zi2];

%% Compute empirical CDF from rescaled waiting times.

KS_Test(Tau)

end % loop for folds
end % end of time window
%% Time rescaling theorem overall
% KS_Test(Tau_All)

%% save output 
if saveOutput
% save([project_workspace_path 'Output/Monkey_IMI_3h_out.mat'], 'mod_set'); 
% save([project_workspace_path 'Output/Fig2A_IMI3_5h_out.mat'], 'mod_set'); 
end
% save([project_workspace_path 'Output/Human_IMI_3h_out.mat'], 'mod_set'); 
% save([project_workspace_path 'Output/VarTrial_IMI_3h_out.mat'], 'mod_set'); 
% save([project_workspace_path 'Output/Goldfish_IMI3_6h_out.mat'], 'mod_set'); 

% save([project_workspace_path 'Output/Izhikevich_burst_IMI3_5h_out.mat'], 'mod_set'); 
% save([project_workspace_path 'Output/Monkey_IMI3_5h_out.mat'], 'mod_set'); 








