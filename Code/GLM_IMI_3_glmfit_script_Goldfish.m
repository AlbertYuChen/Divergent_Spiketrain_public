% Author: Yu Chen
% Date: Nov 1st 2017 @ CNBC

close all; clear; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% Load Data and parameter settings
dataset = 'C';
stim_len = 30*1000;

%---------------- Truccolo ------------------------- 
load([project_workspace_path 'Data/Goldfish.mat'])
%-------------- convert this to my data type ------------
MC(1).spikeTimes = round(raster_x' / 0.001);
MC(1).spikeIndices = raster_y';

%% Load Data and parameter settings
neuron_index = 1;
saveOutput = 0; % 1: SAVE  0: NOT Save

spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;

 
Time_fold_index = {(200:stim_len)'};   

valid_point_range = (200:stim_len)';
KS_valid_point_range = valid_point_range;

Trial_Index = 1;
N_trails = length(Trial_Index);

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

%% ===================== Setup fitting params ======================
dtStim = stim_rec_res; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = spike_rec_res;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
gg0.dtStim = dtStim;
gg0.dtSp = dtSp; 
gg0.sps = cat_train';  % Insert binned spike train into fitting struct

%% ===================== H design =====================
SPK_Num = 5; % Goldfish original 6 

gg0.SPK_Num = SPK_Num;
HF_u_accm = mat2cell( zeros(SPK_Num,0), ones(SPK_Num, 1) );
histH = {};

basis_xlim = [0 0.3]; % Izhikevich-burst

for ti = Trial_Index
    spike_index_array = MC(neuron_index).spikeIndices == ti;
    single_trail_spike_time = MC(neuron_index).spikeTimes(spike_index_array) * spike_rec_res;
    
    for hi = 1:min(SPK_Num, length(single_trail_spike_time))
        HF_u = get_IMI_input_time(single_trail_spike_time, stim_time, hi);
        %------- not every spike start from 0 --------
        zero_area = round(single_trail_spike_time(hi)/stim_rec_res);
        HF_u(1:zero_area) = NaN;
 
        HF_u_tmp = HF_u(valid_point_range);
        HF_u_tmp(isnan(HF_u_tmp)) = [];    
        HF_u_accm{hi} = [HF_u_accm{hi}, HF_u_tmp];
    end        
end

gg0.HF_u_accm = HF_u_accm;

% distribution of hazard function data
HIST_bin_width = 0.002; % Goldfish 
figure('Position', [300, 100, 600, 500]); 

for hi = 1:SPK_Num
    histH{hi} = histogram(HF_u_accm{hi}, 'Normalization', 'pdf', 'BinWidth', HIST_bin_width, 'FaceAlpha', .7);  hold on
end


% leg1 = legend('$u_1$','$u_2$','$u_3$','$u_7$','$u_5$','$u_6$');
% set(leg1,'Interpreter','latex');
ylabel('PDF', 'interpreter', 'latex')
xlabel('time [sec]', 'interpreter', 'latex') 
title('\textit{Monkey-PMv}', 'interpreter', 'latex') 
xlim(basis_xlim)
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex');

%% ======================= H =======================
h_num_basis = [7 6 5 4 4 4 4 4 4 4]; % burst

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
%         plot(HF_u); xlim([0 length(HF_u)])
    end
end

%% Design matrix
% create spike train index for all concatenated spike train
irt_ind_cum = [];
for ti = 1:N_trails  
    irt_ind_cum = [irt_ind_cum; (valid_point_range+ stim_len*(ti-1))];
end

X_design = [ones(length(valid_point_range)*N_trails, 1), cell2mat(Xsp)];

%% ============= GLM ================
% [prsML,dev,stats,loglikeli,iter,b_stack] = glmfit_YC(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');
[prsML,dev,stats] = glmfit(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');

%% ============================= DC =============================
gg0.dc = prsML(1);
gg0.kdc = gg0.dc*ones(stim_len, 1);

%% ============================= H =============================
gg0.ihw = {};
gg0.ih ={}; 

nh_arry = cell2mat(gg0.nh);
nh_arry = [0 cumsum(nh_arry)];

figure('Position', [300, 100, 600, 500]); hold on
for hi = 1:SPK_Num
b_h_ind = 1+nh_arry(hi)+1:1+nh_arry(hi+1);
gg0.ihw{hi} = prsML(b_h_ind);
gg0.ih{hi} = gg0.ihbas{hi}*gg0.ihw{hi}; 
%----

x = gg0.iht{hi};
y = gg0.ih{hi}; 
plot(x, y, '-', 'linewidth', 3)
end

xlabel('time [sec]', 'interpreter', 'latex')
ylabel('log firing rate', 'interpreter', 'latex') 
legend({'$h_{1}$', '$h_{2}$', '$h_{3}$', '$h_{4}$', '$h_{5}$'}, 'interpreter', 'latex', 'location', 'southeast');
title('\textit{Goldfish} $\mathrm{FNF_M}$ filters', 'interpreter', 'latex') 
xlim([0 0.4])
ylim([-3 2]) 
grid on
set(gca, 'FontSize',18, 'TickLabelInterpreter', 'latex');
%% ============================= Misc =============================
gg0.valid_point_range = valid_point_range;
gg0.stim_time = stim_time;
gg0.neuron_index = neuron_index;
gg0.Trial_Index = Trial_Index;

%% KS analysis
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

%% KS test
KS_Test(Tau)
KS_Test_Haslinger(Tau);

%% save output 
% save([project_workspace_path 'Output/Goldfish_IMI3_h_' num2str(SPK_Num) '_out.mat'], 'mod_set');  





