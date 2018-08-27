
function u = get_IMI_input_time(single_trail_spike_time, stim_time, order)

X = [zeros(1, order) single_trail_spike_time (stim_time(end) + 0.1)*ones(1, order)];
u = stim_time;
base_spike_ind = 1;

for ii = 1:length(stim_time)
    
    if u(ii) - X(base_spike_ind+order) > 0.000001
        base_spike_ind = base_spike_ind + 1;
    end
    
    u(ii) = u(ii) - X(base_spike_ind);
    
    if u(ii) < 0.00001
        disp('HA')
    end
end

