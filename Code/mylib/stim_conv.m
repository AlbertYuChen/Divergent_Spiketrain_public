

function Y = stim_conv(stim, K, delay)

[klen, kwid] = size(K);

Y = conv2(stim, K, 'full');

% Y = Y(1:end-klen+1, :); % old Pillow's method

% shift the stimulus to future
% Y = [Y(10:end-klen+1,:); zeros(10, kwid)];

% shift the stimulus to the past
if nargin >= 3
Y = [zeros(delay, kwid); Y];
else
Y = [zeros(0, kwid); Y];    % original as Pillow's method
% Y = [zeros(1, kwid); Y];
end



