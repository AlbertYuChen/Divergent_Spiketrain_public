

function Y = hist_conv(sps, H)

[hlen, hwid] = size(H);

Y = [];
for ii = 1:size(sps, 2)
% Do convolution and remove extra bins
y = conv2(sps(:,ii), H, 'full');
y = [zeros(1, hwid); y(1:end-hlen,:)];
Y = [Y, y];
end









