function [tsp,sps,Itot,Istm] = simGLM_INTENSITY_SpDi(glmprs, stim_len, rndseed)
%% concatenate intensity functions
kdc = max(glmprs.kdc)*ones(stim_len, 1);
kdc(glmprs.valid_point_range) = glmprs.kdc(glmprs.valid_point_range);

% figure('Position', [300, 300, 1400, 500]);
% plot(kdc)

%% start simulation
rng(rndseed);
nlfun = @exp;
glmprs.dtSp = 0.001;
% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp; % bin size for simulation

rlen = stim_len;  % length of binned spike response
hlen = size(glmprs.ih,1); % length of post-spike filter

% -------------  Compute filtered resp to signal ----------------
% Istm = glmprs(2).kdc;
Istm = kdc; 
Itot = Istm; % total filter output

% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
sps = zeros(rlen,1); % sparse spike time matrix
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = nlfun(Itot(iinxt)); %----- *dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); 
        % randomly select the boundary caused by the discrete bins.         
        tmpspk = find(rrcum>=tspnext, 1, 'first');
        sumright = rrcum(tmpspk);             
        if tmpspk > 1      
            sumleft = rrcum(tmpspk-1);
            ppleft = (sumright - tspnext)/(sumright - sumleft);         
            if ppleft > 0.5
                ispk = ispk - 1;
            end
        end
        
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
            try
            Itot(iiPostSpk) = Itot(iiPostSpk) + glmprs.ih(1:mxi-ispk);
            catch
                mxi
            end
        end
        
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
tsp = find(sps>0); %-----*dt;




















