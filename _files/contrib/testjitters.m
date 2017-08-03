% neuroelf library
n = neuroelf;

% number of tests
niter = 100;

% TR (ms)
prtr = 2000;

% first onset (ms)
firstons = 6000;

% additional time at end (ms)
plustime = 8000;

% design spec, number of conditions
numconds = 3;

% number of trials per condition (per run)
numtrialcond = 18;

% stim time (e.g. cue, pic, rating)
cstims = [2000, 6000, 3000];

% seperate predictors for those periods (per condition)?
sepstim = [false, true, false];

% jitter means (ms)
cjits = [1500, 1500, 3200];
% cjits = [1000, 1000, 1500];

% jitter distribution widths (SD, ms)
cjitw = [600, 600, 2000];
% cjitw = [400, 400, 1000];

% jitter-distribtion type (n=normal, e=exponential)
cjitd = {'n', 'n', 'e'};

% jitter max-length values (for exponential distributions, ms)
ejitm = [0, 0, 8000];
% ejitm = [0, 0, 5000];

% cell array for jitter distributions
ejits = cell(1, 3);

% contrasts of interest; depends on the "dropped" conditions!!
cons = [ ...
     0,  1,  0,  0,  0, 0; ...
     0,  0,  1,  0,  0, 0; ...
     0,  0,  0,  1,  0, 0; ...
     0,  1, -1,  0,  0, 0; ...
     0,  0,  1, -1,  0, 0]';

% compute total length (for PRT -> SDM conversion)
tlength = firstons + numconds * numtrialcond * sum(cstims + cjits) + plustime;

% target conditions
tcs = reshape(1:(numconds*numel(cstims)), numconds, numel(cstims));
for cc = 1:numel(cstims)
    if ~sepstim(cc)
        tcs(:, cc) = min(tcs(:, cc));
    end
end

% initialize iXX's
utcs = numel(unique(tcs(:)));
iXX = zeros(utcs + 1, utcs + 1, niter);
pconds = cell(1, niter);

% iterate
for ic = 1:niter
    
    % set start time
    stime = firstons;
    
    % generate PRT
    prt = xff('new:prt');
    
    % exponential distributions
    for cc = 1:numel(cstims)
        if cjitd{cc} == 'e'
            erange = log([cjits(cc) - cjitw(cc), ejitm(cc)]);
            ejits{cc} = exp(erange(1):(erange(2)-erange(1))/(numconds*numtrialcond - 1):erange(2));
            ejits{cc} = ejits{cc} - max(0, (mean(ejits{cc}) - cjits(cc)));
            ejits{cc} = ejits{cc}(randperm(numel(ejits{cc})));
        end
    end
    cjitc = [1, 1, 1];

    % add conditions
    for sc = 1:numel(cstims)
        for cc = 1:numconds
            if cc == 1 || sepstim(sc)
                prt.AddCond(sprintf('c%d_s%d', cc, sc), zeros(0, 2));
            else
                prt.AddCond(sprintf('delete%d', cc + numel(cstims) * (cc-1)), zeros(0, 2));
            end
        end
    end
    
    % generate list of items
    clist = 1 + mod(randperm(numconds * numtrialcond), numconds);

    % iterate over trials
    for tc = 1:numel(clist)
        
        % trial type
        ttype = clist(tc);
        
        % iterate over stims
        for sc = 1:numel(cstims)
            
            % add stim
            prt.Cond(tcs(ttype + (sc - 1) * numconds)).OnOffsets(end+1, :) = ...
                [stime, stime + cstims(sc)];
            
            % advance time
            if cjitd{sc} == 'n'
                stime = round(stime + cstims(sc) + cjits(sc) + cjitw(sc) * randn(1, 1));
            else
                stime = round(stime + cstims(sc) + ejits{sc}(cjitc(sc)));
                cjitc(sc) = cjitc(sc) + 1;
            end
        end
    end
    
    % remove empty conditions
    cc = prt.Cond;
    cc = {cc.OnOffsets};
    ecc = cellfun('isempty', cc);
    prt.Cond(ecc) = [];
    prt.NrOfConditions = numel(prt.Cond);
    
    % create SDM
    sdm = prt.CreateSDM(struct('prtr', prtr, 'nvol', ceil(tlength / prtr), 'rcond', []));
    
    % compute iXX
    iXX(:, :, ic) = inv(sdm.SDMMatrix' * sdm.SDMMatrix);
    
    % keep condition names (and content for later)
    cnames = sdm.PredictorNames;
    pconds{ic} = prt.Cond;
    
    % clear objects
    sdm.ClearObject;
    prt.ClearObject;
    if mod(ic, 10) == 0
        disp(ic);
        pause(0.01);
    end
end

% get contrast estimates
cest = zeros(niter, size(cons, 2));
for cc = 1:size(cons, 2)
    
    % compute contrast se terms
    cest(:, cc) = 1 ./ n.mtimesnd(n.mtimesnd(repmat(cons(:, cc)', [1, 1, niter]), iXX), repmat(cons(:, cc), [1, 1, niter]));
end

% scatter some contrast values
figure;
subplot(2, 2, 1);
scatter(cest(:, 1), cest(:, 2));
subplot(2, 2, 2);
scatter(cest(:, 1), cest(:, 3));
subplot(2, 2, 3);
scatter(cest(:, 1), cest(:, 4));
subplot(2, 2, 4);
scatter(cest(:, 4), cest(:, 5));
