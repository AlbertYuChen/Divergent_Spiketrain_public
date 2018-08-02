% Author: Yu Chen
% Date: 2017 @ CNBC

close all; clear; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% Load Data and parameter settings
dataset = 'C';
stim_len = 1*1000;

%---------------- Truccolo -------------------------
load([project_workspace_path 'Data/GerhardDegerTruccolo_2017_Data/Data_Fig2' dataset '.mat'])

%-------------- convert this to my data type ------------
MC(1).spikeTimes = round(raster_x' / 0.001);
MC(1).spikeIndices = raster_y';

%% Load Data 

spike_rec_res = 0.001;
stim_rec_res = 0.001;
stim_time = 0.001:stim_rec_res:stim_len*stim_rec_res;

valid_point_range = (1:stim_len)';
KS_valid_point_range = valid_point_range;

Trial_Index = 1:10;
Nfold = length(Trial_Index);
N_trails = length(Trial_Index);

all_spike_train = [];
ISIs = [];

for ii = Trial_Index
    spike_index_array = MC(1).spikeIndices == ii;
    all_spike_train = [all_spike_train; hist(MC(1).spikeTimes(spike_index_array) * spike_rec_res, stim_time)];
    single_trail_spike_time = MC(1).spikeTimes(spike_index_array) * spike_rec_res;
    ISIs = [ISIs diff(single_trail_spike_time)]; 
end

cat_train = reshape(all_spike_train', 1, []);


%% ===================== Setup fitting params ======================
dtStim = stim_rec_res; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = spike_rec_res;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
gg0.dtStim = dtStim;
gg0.dtSp = dtSp;
gg0.nlfun = @expfun; % default nonlinearity: exponential
gg0.sps = cat_train';  % Insert binned spike train into fitting struct

%% ===================== intensity function B-spline basis =====================
B_s = valid_point_range(1)*stim_rec_res;
B_e = valid_point_range(end)*stim_rec_res; 

% Monkey
Spline_res = 0.25; % 0.05
spline_order = 4;  % 4 is cubic, 3 is quadratic

knots_intense = [B_s B_s B_s B_s:Spline_res:B_e B_e B_e B_e B_e];  % knot vector
nk = size(knots_intense, 2) - spline_order;

XBbas = [];

for spline_index = 0 : numel(knots_intense) - spline_order - 1
    [yy, x1] = bspline_basis(spline_index, spline_order, knots_intense, valid_point_range*stim_rec_res);    
    yy = [zeros((valid_point_range(1)-1), 1); yy; zeros((stim_len - valid_point_range(end)), 1)];
    XBbas = [XBbas yy];
end

gg0.ktbas = XBbas;
Xstim = kron(ones(N_trails, 1), XBbas(valid_point_range ,:));

%% ===================== H =====================
% Monkey
nh = 6; 
hpeakFinal = 0.4;

ggx = makeFittingStruct_GLM(dtStim, dtSp, 150, 2, zeros(150, 1), nh, hpeakFinal);
gg0.iht = ggx.iht;
gg0.ihbas = ggx.ihbas;

%---------------------------------
Xsp = []; % allocate
% Design matrix for self-coupling filter 

for ti = 1:N_trails  
    trl_ind = stim_len*(ti-1) + (1:stim_len);
    Xsp_irt = hist_conv(gg0.sps(trl_ind), gg0.ihbas);
    if nh>0
    Xsp = [Xsp; Xsp_irt(valid_point_range, :)];
    end
end

%% 4. Do ML fitting %=====================
irt_ind_cum = [];
for ti = 1:N_trails  
    irt_ind_cum = [irt_ind_cum; (valid_point_range+stim_len*(ti-1))];
end

b_k_ind = 1:nk;
b_dc_ind = nk+1;
b_h_ind = nk+2:nk+1+nh;

X_design = [Xstim, ones(length(Xstim),1), Xsp];

%% GLM
[prsML,dev,stats] = glmfit(X_design, gg0.sps(irt_ind_cum), 'poisson', 'constant', 'off');

%%   K  
gg0.kw = prsML(b_k_ind); 
gg0.k = gg0.ktbas*gg0.kw;

%%   DC + KDC  
gg0.dc = prsML(b_dc_ind);
gg0.kdc = gg0.k + gg0.dc;

% ------------------------------
figure('Position', [300, 100, 700, 350]);
hold on
plot(stim_time, gg0.kdc, '-', 'Color', [0.7 0 0])

xlim([valid_point_range(1) valid_point_range(end)]*stim_rec_res)
title('Intensity function');
legend('Intensity_{ML}', 'location', 'northwest');
xlabel('Time [s]')
ylabel('log firing rate')
grid on

%%   H  
gg0.ihw = prsML(b_h_ind);
gg0.ih = gg0.ihbas*gg0.ihw;

%------- History dependency filters -------
figure('Position', [300, 300, 400, 250]);
plot(gg0.iht, gg0.ih, '-', 'Color', [0.7 0 0])

title('post-spike kernel');
legend('SE', 'h_{ML}', 'location', 'southeast');
grid on
%--------------------------- Misc -----------------------
gg0.valid_point_range = valid_point_range;
gg0.stim_time = stim_time;
gg0.Trial_Index = Trial_Index;

%% Time rescaling theorem analysis for each iteration 
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
gg0.Itot = X_design * prsML;
gg0.lambda = lambda;
gg0.Tau = Tau;
mod_set = gg0;

%% Compute empirical CDF from rescaled waiting times.
KS_Test(Tau)
% KS_Test_Haslinger(Tau)

%% save output 
% save([project_workspace_path 'Output/' 'Fig2' dataset '_INTENSITY_out_' date '.mat'], 'mod_set'); 













