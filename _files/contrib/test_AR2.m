% clear memory
clear classes;

% load GLM
glm = xff('onerun_AR2.glm');
glm.LoadTransIOData;

% get design matrix
X = glm.DesignMatrix;
iXX = inv(X' * X);

% compute semap
numtp = glm.NrOfTimePoints;
numpred = glm.NrOfPredictors;
nval = numtp - numpred - 2 * glm.NrOfStudies;
szmap = size(glm.GLMData.MCorrSS);
mrr = double(glm.GLMData.MultipleRegressionR);
mcs = double(glm.GLMData.MCorrSS);
se = sqrt((1 - mrr .^ 2) .* mcs ./ nval);
se(se == 0) = Inf;
se = se(:)';

% get beta map
b = reshape(glm.GLMData.BetaMaps(:, :, :, :), [prod(szmap), numpred])';

% clear glm
glm.ClearObject;

% mask
bnn = all(b ~= 0, 1);

% and remove voxels where at least two betas are equal (invalid contrast)
bnn(bnn) = all(diff(sort(b(:, bnn), 1)) > 0, 1);
b = b(:, bnn);
se = se(bnn);
mrr = mrr(:)';
mcs = mcs(:)';
mrr = mrr(bnn);
mcs = mcs(bnn);

% which parameter to search (or simply report average error)
sp1 = false;
sp2 = true;

% loop around lag values
l12 = zeros(420, 13);
for l1 = -.2:0.05:.75
    l1t = round(100 * l1);
    for l2 = -.5:0.05:0.5
        l2t = round(100 * l2);
        vmpfiles = findfiles(pwd, ...
            sprintf('orAR2_%03d_%03d_*.vmp', l1t, l2t));
        if isempty(vmpfiles)
            continue;
        end
        l12i = round(21 * 20 * (l1 + 0.2) + 20 * (l2 + 0.55));

        % compute initial weights
        a10 = 1 / (1 / (l1 * l1) - 1);
        a1 = [   0, (1 - l2) * (1 + a10)];
        a2 = [-a10,             1 + a10];
        a1o = a1(1) + a1(2) * l1;
        a2o = a2(1) + a2(2) * l2;
        l12(l12i, 1:4) = [l1, l2, a1o, a2o];
        l12ii = 5;

        % for each contrast with those lags
        for fc = 1:numel(vmpfiles)

            % load test VMP
            vmp = xff(vmpfiles{fc});
            vmpd = double(vmp.Map.VMPData);
            vmpd = vmpd(bnn);

            % contrast vector
            [vmpp, vmpf] = fileparts(vmpfiles{fc});
            cvec = vmpf(15:end);
            cvec(cvec == '-') = '/';
            c = [double(cvec) - 48, 0];

            % compute contrast
            iXX = inv(X' * X);
            ct = (c * b) ./ (sqrt(c * iXX * c') .* se);

            % begin with parameters and step size
            a1 = [   0, (1 - l2) * (1 + a10)];
            a2 = [-a10,             1 + a10];
            ad = 0.01;

            % until solution found
            ds = 1 + mean(abs(vmpd - ct) ./ abs(vmpd));
            dsa111 = ds; dsa112 = ds; dsa211 = ds; dsa212 = ds;
            while ds > 1e-8

                % compute difference matrix
                dX = X(3:end, :) - ...
                     (a1(1) + a1(2) * l1) .* X(2:end-1, :) - ...
                     (a2(1) + a2(2) * l2) .* X(1:end-2, :);

                % invert and compute stats
                iXX = inv(dX' * dX);

                % recompute contrast
                ct = (c * b) ./ (sqrt(c * iXX * c') .* se);

                % recompute sum
                ds = mean(abs(vmpd - ct) ./ abs(vmpd));

                % sweep parameters
                if sp1
                    dX = X(3:end, :) - ...
                        ((a1(1) - ad) + a1(2) * l1) .* X(2:end-1, :) - ...
                        ( a2(1)       + a2(2) * l2) .* X(1:end-2, :);
                    iXX = inv(dX' * dX);
                    ctt = (c * b) ./ (sqrt(c * iXX * c') .* se);
                    dsa111 = mean(abs(vmpd - ctt) ./ abs(vmpd));
                    dX = X(3:end, :) - ...
                        ((a1(1) + ad) + a1(2) * l1) .* X(2:end-1, :) - ...
                        ( a2(1) +       a2(2) * l2) .* X(1:end-2, :);
                    iXX = inv(dX' * dX);
                    ctt = (c * b) ./ (sqrt(c * iXX * c') .* se);
                    dsa112 = mean(abs(vmpd - ctt) ./ abs(vmpd));
                end
                if sp2
                    dX = X(3:end, :) - ...
                        ( a1(1)       + a1(2) * l1) .* X(2:end-1, :) - ...
                        ((a2(1) - ad) + a2(2) * l2) .* X(1:end-2, :);
                    iXX = inv(dX' * dX);
                    ctt = (c * b) ./ (sqrt(c * iXX * c') .* se);
                    dsa211 = mean(abs(vmpd - ctt) ./ abs(vmpd));
                    dX = X(3:end, :) - ...
                        ( a1(1) +       a1(2) * l1) .* X(2:end-1, :) - ...
                        ((a2(1) + ad) + a2(2) * l2) .* X(1:end-2, :);
                    iXX = inv(dX' * dX);
                    ctt = (c * b) ./ (sqrt(c * iXX * c') .* se);
                    dsa212 = mean(abs(vmpd - ctt) ./ abs(vmpd));
                end

                % which parameter worked best
                dst = [ds, dsa111, dsa112, dsa211, dsa212];
                psx = minpos(dst);

                switch (psx)

                    % nothing worked better
                    case {1}

                        % reduce ad
                        ad = 0.25 * ad;
                        if ad < 1e-8
                            break;
                        end

                    % alter parameter
                    case {2}
                        a1(1) = a1(1) - ad;
                    case {3}
                        a1(1) = a1(1) + ad;
                    case {4}
                        a2(1) = a2(1) - ad;
                    case {5}
                        a2(1) = a2(1) + ad;
                end
            end

            % show results
            l12(l12i, l12ii:l12ii+2) = [ ...
                a1(1) + a1(2) * l1, a2(1) + a2(2) * l2, ds];
            disp([sprintf('%2d,', c), ...
                sprintf('    %6.2g, %6.2g, %11g->%11g, %11g->%11g, %.10f; ...', l1, l2, ...
                a1o, l12(l12i, l12ii), a2o, l12(l12i, l12ii+1), l12(l12i, l12ii+2))]);
            pause(0.001);
            l12ii = l12ii + 3;
        end
    end
end
