function [tsp,sps,Itot,Istm] = simGLM_IMI_SpDi(glmprs, stim_len, rndseed)
rng(rndseed);
glmprs.nlfun = @exp;

% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
% dt = glmprs(1).dtSp; % bin size for simulation

slen = stim_len; % length of stimulus
rlen = slen;  % length of binned spike response


kdc = max(glmprs.kdc)*ones(stim_len, 1);
kdc(glmprs.valid_point_range) = glmprs.kdc(glmprs.valid_point_range);

Istm = kdc;
Itot = Istm; % total filter output
    
% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
sps = zeros(rlen,1); % sparse spike time matrix
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point

spk_array = zeros(1, 50);

% --------------- Append zeros ---------------------------------------
if isfield(glmprs, 'ih') && ~iscell(glmprs.ih) 
glmprs.ih = [glmprs.ih; zeros(40*1000, 1)];
else
    for filter_ind = 1:glmprs.SPK_Num    
        glmprs.ih{filter_ind} = [glmprs.ih{filter_ind}; zeros(20*1000, 1)];
    end
end

% Those cases below are left for backward compatibility
if isfield(glmprs, 'ih1')
glmprs.ih1 = [glmprs.ih1; zeros(10*1000, 1)];
end
if isfield(glmprs, 'ih2')
    glmprs.ih2 = [glmprs.ih2; zeros(10*1000, 1)];
end
if isfield(glmprs, 'ih3')
    glmprs.ih3 = [glmprs.ih3; zeros(10*1000, 1)];
end
if isfield(glmprs, 'ih4')
    glmprs.ih4 = [glmprs.ih4; zeros(300, 1)];
end
if isfield(glmprs, 'ih1h1h1')
    glmprs.ih1h1h1 = [glmprs.ih1h1h1; zeros(800, 1)];
end

if isfield(glmprs, 'irect')
    rect_fun = glmprs.irect;
elseif isfield(glmprs, 'SPK_Num')
    rect_fun = ones(glmprs.SPK_Num, 1);
end

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    if isfield(glmprs, 'ih')
        log_Firing = Itot(iinxt);
        
        if ~iscell(glmprs.ih) % for IMI-1,2: different filters. 
            for filter_ind = 1:glmprs.SPK_Num
                if spk_array(filter_ind) ~= 0
                    log_Firing = log_Firing + rect_fun(filter_ind) * glmprs.ih( iinxt - spk_array(filter_ind));
                end
            end
            rrnxt = glmprs.nlfun( log_Firing );
        else % for IMI-3: different filters. 
            for filter_ind = 1:glmprs.SPK_Num
                if spk_array(filter_ind) ~= 0
                    log_Firing = log_Firing + glmprs.ih{filter_ind}( iinxt - spk_array(filter_ind));
                end
            end
            rrnxt = glmprs.nlfun( log_Firing );
        end    
        
    % Those cases below are left for backward compatibility
    elseif isfield(glmprs, 'ih4')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1(iinxt - spk_array(1))...
             + glmprs.ih2(iinxt - spk_array(2)) + glmprs.ih3(iinxt - spk_array(3)) + glmprs.ih4(iinxt - spk_array(4)) );    
    elseif isfield(glmprs, 'ih3')
        rrnxt = glmprs.nlfun(Itot(iinxt) + ...
            glmprs.ih1(iinxt - spk_array(1)) + glmprs.ih2(iinxt - spk_array(2)) + glmprs.ih3(iinxt - spk_array(3)) ); 
    elseif isfield(glmprs, 'ih2')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1(iinxt - spk_array(1)) + glmprs.ih2(iinxt - spk_array(2))); 
    elseif isfield(glmprs, 'ih1')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1(iinxt - spk_array(1))); 
    elseif isfield(glmprs, 'ih1h1')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih12(iinxt - spk_array(1)) + glmprs.ih12(iinxt - spk_array(2))); 
    elseif isfield(glmprs, 'ih1h1h1')
        rrnxt = glmprs.nlfun(Itot(iinxt) + glmprs.ih1h1h1(iinxt - spk_array(1)) + glmprs.ih1h1h1(iinxt - spk_array(2)) ...
            + glmprs.ih1h1h1(iinxt - spk_array(3)));       
    else
        disp('does not contain history dpendency function')
    end
    
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    
    % No spike in this window
    if (tspnext >= rrcum(end)) 
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
        
    % Spike! ! !   
    else   
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin

        %---- update synthesised firing rate ----
        iiPostSpk = (spk_array(1)+1):ispk;        
        if isfield(glmprs, 'ih')
            if ~iscell(glmprs.ih) 
                for filter_ind = 1:glmprs.SPK_Num
                    Itot(iiPostSpk) = Itot(iiPostSpk) + rect_fun(filter_ind)*glmprs.ih(iiPostSpk - spk_array(filter_ind));
                end
            else
                for filter_ind = 1:glmprs.SPK_Num
                    Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih{filter_ind}(iiPostSpk - spk_array(filter_ind));
                end
            end
        elseif isfield(glmprs, 'ih3')
            Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih1(iiPostSpk - spk_array(1))... 
                + glmprs.ih2(iiPostSpk - spk_array(2)) + glmprs.ih3(iiPostSpk - spk_array(3));
        elseif isfield(glmprs, 'ih2')
            Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih1(iiPostSpk - spk_array(1))... 
                + glmprs.ih2(iiPostSpk - spk_array(2));
        elseif isfield(glmprs, 'ih1')
            Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih1(iiPostSpk - spk_array(1));
        else
            disp('does not contain history dpendency function')
        end
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
        
        spk_array = circshift(spk_array, 1);
        spk_array(1) = ispk;
            
    end
end
tsp = find(sps>0);


























