% Author: Yu Chen
% Date: Jan 18th 2017 @ CNBC
% A self designed Gaussian filter
% smoothed_signal = Gaussian_Filter_YC(original_signal, window_width)
% change the windwo size according to your purpose.

function smoothed_signal = Gaussian_Filter_YC(original_signal, window_width)

halfWidth = window_width / 2;
gaussFilter = gausswin(window_width);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
smoothed_signal = conv(original_signal, gaussFilter);
smoothed_signal = smoothed_signal(halfWidth:end-halfWidth);

end