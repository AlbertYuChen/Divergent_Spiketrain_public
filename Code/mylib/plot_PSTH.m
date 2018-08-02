% Author: Yu Chen

function PSTH_reshape = plot_PSTH(all_spike_train, PSHT_time_bin, stim_time, stim, x_lim, Gauss_width, lambda_mod1)
plot_lambda = 0;
if nargin > 8
    plot_lambda = 1;
end

% figure('Position', [200, 300, 1500, 350]);
% subplot(211);
% plot(stim_time, stim , 'Color', [0.5 0.5 0.8])  
% xlim(x_lim)
% xlabel('Time (sec)')
% ylabel('Stimulated Spike rate')
% 
% subplot(212);
PSTH = mean(all_spike_train, 1);	%Compute histogram.
PSTH_reshape = mean(reshape(PSTH, length(PSTH)/length(PSHT_time_bin),[]), 1);

if Gauss_width > 0
PSTH_reshape = Gaussian_Filter_YC(PSTH_reshape, Gauss_width);
end 

% bar(PSHT_time_bin, PSTH, 1, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.75 0.75 0.75]);	%...and display PSTH.
% plot(PSHT_time_bin, PSTH_reshape, 'Color', [0.75 0.75 0.75]);	%...and display PSTH.
% ylabel('Spike rate (spikes/s)')
% xlim(x_lim)

if plot_lambda == 1
    
hold on;
yyaxis right
plot(stim_time, mean(lambda_mod1, 2), '.', 'Color', [0 0 0.5])  
xlim(x_lim)
% title('PSHT')
xlabel('Time (sec)')
ylabel('Stimulated Spike rate')
set(gca,'FontSize',18)
end

