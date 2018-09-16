% Date: Aug 25th 2017
% Author: Yu Chen @ CNBC, CMU, Pittsburgh

function [BasisBegin, BasisEnd, NodeArry] = NodeArray_Distribution(histH, nodeNumber)

% BasisBegin = histH.BinEdges( find(histH.Values > 1, 1));
% BasisEnd = histH.BinEdges( find(histH.Values > 0.01, 1, 'last'));
BasisBegin = histH.BinLimits(1);
BasisEnd = histH.BinLimits(2);

stepWidth = 1/nodeNumber;
cdfStep = stepWidth:stepWidth:(1-stepWidth);

HistCDF = cumsum(histH.Values*histH.BinWidth);
[HistCDF, HistCDFindex] = unique(HistCDF); 
NodeArry = interp1(HistCDF, histH.BinEdges(HistCDFindex), cdfStep);

% figure('Position', [200, 200, 1200, 600]);
% plot(histH.BinEdges(HistCDFindex), HistCDF, 'b', 'LineWidth', 4); hold on
% plot(NodeArry, zeros(length(cdfStep),1), 'x', 'LineWidth', 4)
% ylabel('CDF')
% 
% yyaxis right
% plot(histH.BinEdges(HistCDFindex), histH.Values(HistCDFindex), 'g','LineWidth', 4);
% ylabel('PDF')
% % xlim([0 0.15])
% xlabel('time [sec]')
% grid on; set(gca,'FontSize',12)
% legend('CDF','B-spline nodes','PDF')








