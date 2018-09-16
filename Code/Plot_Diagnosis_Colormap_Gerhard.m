% Author: Yu Chen
% Date: Jan 13th 2017 @ CNBC

close all; clear; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

% I have calculated the result for you. If you want to calculate the
% performance of all combination, see code
% ExploretionRegion_diagnosis_Gerhard.m
load([project_workspace_path 'Output/Z_Diagnosis_Colormap_Gerhard.mat'])

%%
for y=1:length(b2_range)
    temp = find(z_real(:,y));
    if isempty(temp)
        edge(y) = b1_range(end)+1;
    else
        edge(y) = (temp(1)-1)*DPI_1 + b1_range(1);
    end
    
end

% ------------- Gerhard method ----------------
figure('Position', [300, 100, 600, 500]);
% imagesc(z);
z_r = rot90(z_origin,1);
imagesc(b1_range, fliplr(b2_range), z_r);
set(gca, 'Ydir','Normal')
colormap([ 0.6660 0.8740 0.2880; .2 0.6470 0.9410; .9 0.55   0.5]);
hold on

for yy = b2_range
    h4 = plot([edge(round((yy+3)/DPI_2+1)), 12],[yy, yy],':k', 'LineWidth',1.5);
end
title('Diagnostic of Gerhard et al.', 'interpreter', 'latex')
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex');

% ------------- my method ----------------
figure('Position', [300, 100, 600, 500]);
% imagesc(z);
z_r = rot90(z_myself,1);
imagesc(b1_range, fliplr(b2_range), z_r);
set(gca, 'Ydir','Normal')
colormap([ 0.6660 0.8740 0.2880; .2 0.6470 0.9410; .9 0.55   0.5]);
hold on

for yy = b2_range
    h4 = plot([edge(round((yy+3)/DPI_2+1)), 12],[yy, yy],':k', 'LineWidth',1.5);
end


% only plot for colors
h1 = fill([10 12 12 10], [1 1 2 2], [ 0.6660 0.8740 0.2880], 'linestyle','none');
h2 = fill([10 12 12 10], [1 1 2 2], [.2 0.6470 0.9410], 'linestyle','none');
h3 = fill([10 12 12 10], [1 1 2 2], [.9 0.55   0.5], 'linestyle','none');

legend([h1 h2 h3 h4], {'Stable', 'Fragile', 'Divergent', 'Unstable'},...
    'Interpreter','latex', 'location', 'southwest');

xlabel('$\beta_1$' , 'interpreter', 'latex')
ylabel('$\beta_2$', 'interpreter', 'latex')


title('Updated diagnostic', 'interpreter', 'latex')
set(gca, 'FontSize',18, 'TickLabelInterpreter','latex');























