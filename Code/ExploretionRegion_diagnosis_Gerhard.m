% Author: Yu Chen
% Date: Jan 12th 2017 @ CNBC

clear; close all; clc

%% check working environment
project_workspace_path = 'D:/Divergent_Spiketrain_public/';

%% Pillow Generate basis
nBases=10;
Baseslength=400;
tstrech=5;
raw_bs = basisFactory.makeNonlinearRaisedCos(nBases, 1, [0 Baseslength], tstrech);
basis = [raw_bs.B(58:407,6)';raw_bs.B(58:407,9)'];
basis1 = basis(1,:);
basis2 = basis(2,:);

%% -------------- Gerhard Basis --------------
t = (1:200)/1000;
basis1 = exp(-t/0.02);
basis2 = exp(-t/0.1);
dip = zeros(1, 200);
dip(1:2) = -1e12;
% ------------------------------------------ 
figure('Position', [300, 300, 600, 500]);
plot(t, basis1, 'linewidth', 3); hold on
plot(t, basis2, 'linewidth', 3);

b1 = 2;
b2 = 2;
% --------- add absolute refractory period ---------
h = b1*basis1 + b2*basis2 + dip;

%% start loop

DPI_1 = 0.2;
DPI_2 = 0.1;
b1_range = -10:DPI_1:10;
b2_range = -3:DPI_2:3;

z_real = zeros(length(b1_range),length(b2_range));
z_origin = zeros(length(b1_range),length(b2_range));
z_myself = zeros(length(b1_range),length(b2_range));

baseline = log(5/1000);
for b1=b1_range
    for b2=b2_range
        disp(['b1: ' num2str(b1) '  b2: ' num2str(b2)])      
        % --------- simulation ----------
        [result, ~] = whether_div_real(b1*basis1+b2*basis2+dip, baseline);
        z_real( round( (b1-b1_range(1) )/DPI_1+1) , round( (b2- b2_range(1) )/DPI_2+1) ) = result;
        % --------- Gerhard ----------
        result = whether_div_Yu_origin(b1*basis1+b2*basis2+dip, baseline);
        z_origin( round( (b1-b1_range(1) )/DPI_1+1) , round( (b2- b2_range(1) )/DPI_2+1) ) = result;
        % --------- Our method ----------
        result = whether_div_Yu(b1*basis1+b2*basis2+dip, baseline);
        z_myself( round( (b1-b1_range(1) )/DPI_1+1) , round( (b2- b2_range(1) )/DPI_2+1) ) = result;
    end
end

figure('Position', [300, 100, 600, 500]);
% imagesc(z);
z_r = rot90(z_real,1);
imagesc(b1_range, fliplr(b2_range), z_r);
set(gca, 'Ydir','Normal')
colormap([ 0.6660 0.8740 0.2880; .2 0.6470 0.9410; .9 0.55   0.5]);

figure('Position', [300, 100, 600, 500]);
% imagesc(z);
z_r = rot90(z_origin,1);
imagesc(b1_range, fliplr(b2_range), z_r);
set(gca, 'Ydir','Normal')
colormap([ 0.6660 0.8740 0.2880; .2 0.6470 0.9410; .9 0.55   0.5]);


figure('Position', [300, 100, 600, 500]);
% imagesc(z);
z_r = rot90(z_myself,1);
imagesc(b1_range, fliplr(b2_range), z_r);
set(gca, 'Ydir','Normal')
colormap([ 0.6660 0.8740 0.2880; .2 0.6470 0.9410; .9 0.55   0.5]);

%% save output
% save([project_workspace_path 'Output/Z_Diagnosis_Colormap_Gerhard.mat']);  
% save([project_workspace_path 'Output/Z_Diagnosis_Colormap_Gerhard_full.mat']);  


















