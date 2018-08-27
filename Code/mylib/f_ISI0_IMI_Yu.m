function [ Expectation_L ] = f_ISI0_IMI_Yu( mod_set, ISI )

if isfield(mod_set, 'irect')
    rect_fun = mod_set.irect;
elseif isfield(mod_set, 'SPK_Num')
    rect_fun = ones(mod_set.SPK_Num, 1);
end


if isfield(mod_set, 'ih')
hlen = length(mod_set.ih);
mod_set.ih = [mod_set.ih; zeros(10*1000, 1)];


iinxt = 1:(hlen*1.5);
log_lambda = mod_set.kdc;

% figure
% plot(log_lambda(iinxt) )
% ylim([-5 -2])
% hold on

for filter_ind = 1:mod_set.SPK_Num
    log_lambda = log_lambda(iinxt) + rect_fun(filter_ind) * mod_set.ih( iinxt + round(ISI)*(filter_ind-1) );
%     figure
%     plot(log_lambda)
end

end


lamda = exp(log_lambda);

% figure
% plot( lamda )


PDF_L = lamda .* exp( -cumsum(lamda));

% sum(PDF_L)
% figure
% plot( PDF_L )

L = (1:length(PDF_L))';

Expectation_L = sum(PDF_L.*L);

