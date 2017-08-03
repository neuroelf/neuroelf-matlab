function [betas, irtc, ptc, se, a] = calcbetas(rtc, tcd, dim, tol, w, ar)
% calcbetas  - perform GLM regression
%
% FORMAT:       [b, irtc, ptc, se, a] = calcbetas(X, d [, td [, tol [, w [, ar]]]])
%
% Input fields:
%
%       X           design matrix (TxR)
%       d           time course data (N-dim)
%       td          temporal dimension (default: 1)
%       tol         tolerance value to set iXX to 0 (default: 4 * eps)
%       w           regression weights (Tx1 or N-dim-x-T)
%       ar          auto-regression flag, either {0}, 1, 2
%
% Output fields:
%
%       b           betas, in time course data dimension
%       irtc        inverse design matrix
%       ptc         predicted tc
%       se          standard error
%       a           auto-regression maps
%
% Notes: the ar flag (1 or 2) produces a BrainVoyager-compatible output

% Version:  v1.0
% Build:    16010911
% Date:     Jan-09 2016, 11:47 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2013, 2016, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% argument check
if nargin < 2 || ...
   ~isa(rtc, 'double') || ...
    ndims(rtc) ~= 2 || ...
    size(rtc, 2) > size(rtc, 1) || ...
    isempty(rtc) || ...
    any(isinf(rtc(:)) | isnan(rtc(:))) || ...
   ~isnumeric(tcd) || ...
   ~any(size(tcd) == size(rtc, 1)) || ...
    isempty(tcd)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad design matrix or time course data supplied.' ...
    );
end
numrows = size(rtc, 1);
numbets = size(rtc, 2);

% further arguments
if nargin < 3 || ...
   ~isa(dim, 'double') || ...
    isempty(dim)
    rdim = findfirst(size(tcd) == numrows);
elseif numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
    fix(dim) ~= dim || ...
    dim < 1 || ...
    dim > 5 || ...
    size(tcd, real(dim)) ~= numrows
    error( ...
        'neuroelf:BadArgument', ...
        'Bad tdim argument supplied.' ...
    );
else
    rdim = real(dim);
end
if nargin < 4 || ...
   ~isa(tol, 'double') || ...
    isempty(tol)
    tol = 4 * eps;
elseif numel(tol) ~= 1 || ...
    isnan(tol) || ...
    abs(tol) > sqrt(eps)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad tol argument supplied.' ...
    );
else
    tol = abs(tol);
end
arw = 0;
if nargin == 5 && ...
    isa(w, 'double') && ...
    numel(w) == 1 && ...
   (isinf(w) || ...
    any(w == [0, 1, 2]))
    arw = w;
    w = [];
end
if nargin < 5 || ...
   ~isa(w, 'double') || ...
   (numel(w) ~= numrows && ...
    size(w, rdim) ~= numrows) || ...
   (size(w, rdim) == numrows && ...
    numel(w) ~= numrows && ...
   ~isequal(size(w), size(tcd)))
    w = [];
end
if nargin < 6 || ...
   ~isa(ar, 'double') || ...
    numel(ar) ~= 1 || ...
   (~isinf(ar) && ...
    ~any(ar == [0, 1, 2]))
    ar = arw;
end

% make data double
if ~isa(tcd, 'double')
    tcd = double(tcd);
end

% check dim
if rdim > 1
    neword = [rdim, 1:length(size(tcd))];
    newodr = find(neword == rdim);
    neword(newodr(2)) = [];
    tcd = permute(tcd, neword);
    if numel(w) == numel(tcd)
        w = permute(w, neword);
    end
    [newsrt, oldord] = sort(neword);
else
    oldord = [];
end

% reshaping data to comply
tcds = size(tcd);
numvox = prod(tcds(2:end));
tcd = reshape(tcd, [tcds(1), numvox]);
if numel(w) == numel(tcd)
    w = reshape(w, [tcds(1), numvox]);
end
tcds(1) = [];
tcdrs = tcds;
if length(tcdrs) < 2
    tcdrs(2) = 1;
end

% no weights
if isempty(w)

    % perform calculus
    if ~issparse(rtc) && ...
        sum(rtc(:) ~= 0) < (0.125 * numel(rtc)) && ...
        ar == 0
        rtc = sparse(rtc);
    end
    trtc = rtc';
    irtc = inv(trtc * rtc);
    irtc(abs(irtc) < tol) = 0;
    prtc = irtc * trtc;
    betas = prtc * tcd;

    % cov-based AR-model
    a = reshape([], [tcds, 0]);
    if isinf(ar)

        % not yet implemented
        error( ...
            'neuroelf:NotYetImplemented', ...
            'Covariance-AR modelling not yet implemented.' ...
        );

    % BVQX AR(1) model
    elseif ar == 1

        % compute residual
        ptc = rtc * betas;
        ptc = tcd - ptc;

        % compute lag map
        [acv, a] = cov_nd(ptc, ptc, -1);
        a(isinf(a) | isnan(a)) = 0;

        % more output?
        if nargout > 2
            ptc(end, :) = [];
            if nargout > 3
                se = Inf .* ones(prod(tcdrs), 1);
            end
        end

        % iterate over unique lags (in 0.05 increments)
        a1lag = 0.05 * round(20 .* a);
        a1lags = unique(a1lag);
        for a1 = a1lags(:)'

            % remove lag from design
            wrtc = rtc(2:end, :) - a1 .* rtc(1:end-1, :);
            twrtc = wrtc';
            iwrtc = inv(twrtc * wrtc);
            iwrtc(abs(iwrtc) < tol) = 0;
            pwrtc = iwrtc * twrtc;

            % find mask
            amask = (a1lag == a1);

            % remove lag from data
            wtcd = tcd(2:end, amask) - a1 .* tcd(1:end-1, amask);
            betas(:, amask) = pwrtc * wtcd;

            % further outputs?
            if nargout > 2
                ptc(:, amask) = wrtc * betas(:, amask);
                if nargout > 3
                    se(amask) = std(tcd(1:end-1, amask) - ptc(:, amask), 0) .* ...
                        sqrt((numrows - 1) / (numrows - numbets));
                end
            end
        end

        % reshape data
        betas = reshape(shiftdim(betas, 1), [tcds, numbets]);
        if nargout > 2
            ptc = reshape(ptc, [numrows - 1, tcds]);
            if nargout > 3
                se = squeeze(reshape(se, tcdrs));
                if nargout > 4
                    a = reshape(a, [tcds, ar]);
                end
            end
        end

        % return early
        return;

    % BVQX AR(2) model
    elseif ar == 2

    end

    % more output?
    if nargout > 2
        ptc = rtc * betas;
        if nargout > 3
            se = squeeze(reshape(std(tcd - ptc, 0), tcdrs)) ...
                 .* sqrt((numrows - 1) / (numrows - numbets));
            if nargout > 4
                a = reshape(a, [tcds, ar]);
            end
        end
        ptc = reshape(ptc, [numrows, tcds]);
        if ~isempty(oldord)
            ptc = permute(ptc, oldord);
        end
    end

% single list of weights
elseif numel(w) == numrows

    % compute weighted design and data
    sw = sqrt(w(:));
    ws = sum(w(:));
    wrtc = rtc .* sw(:, ones(1, numbets));

    % find useless betas
    usebets = (sum(abs(wrtc), 1) > 0);
    if sum(usebets) ~= numbets
        missbetas = true;
        wrtc = wrtc(:, usebets);
    else
        missbetas = false;
    end

    % apply weights to data (without requiring double the memory)
    for tc = 1:numrows
        tcd(tc, :) = sw(tc) .* tcd(tc, :);
    end

    % calculate betas
    irtc = inv(wrtc' * wrtc);
    irtc(abs(irtc) < tol) = 0;
    if missbetas
        betas = zeros(numbets, size(tcd, 2));
        betas(usebets, :) = irtc * wrtc' * tcd;
    else
        betas = (irtc * wrtc') * tcd;
    end

    % more output?
    if nargout > 2
        ptc = rtc * betas;
        if nargout > 3
            sef = numrows / ws;
            se = tcd;
            for tc = 1:numrows
                se(tc, :) = se(tc, :) - sw(tc) .* ptc(tc, :);
            end
            se = squeeze(reshape(sef .* sqrt(sum(se .* se) ./ (ws - (numbets + 1))), tcdrs));
        end
        ptc = reshape(ptc, [numrows, tcds]);
        if ~isempty(oldord)
            ptc = permute(ptc, oldord);
        end
    end

% full weights
else

    % no need to split into packets
    if (numvox * numrows * numbets) <= 2e7

        % compute weighted design and data
        sw = sqrt(w);
        ws = sum(w, 1);
        wrtc = repmat(rtc, [1, 1, numvox]) .* repmat(reshape(sw, [numrows, 1, numvox]), [1, numbets]);

        % apply weights to data
        tcd = sw .* tcd;

        % calculate betas
        irtc = zeros([numbets, numbets, numvox]);
        for vc = 1:numvox
            irtc(:, :, vc) = inv(wrtc(:, :, vc)' * wrtc(:, :, vc));
        end
        irtc(abs(irtc) < tol) = 0;
        betas = mtimesnd(mtimesnd(irtc, wrtc, 2), reshape(tcd, [numrows, 1, numvox]));

        % more output?
        if nargout > 2
            ptc = rtc * reshape(betas, [numbets, numvox]);
            if nargout > 3
                sef = numrows ./ ws;
                se = tcd - sw .* ptc;
                se = squeeze(reshape(sef .* sqrt(max(0, sum(se .* se) ./ (ws - (numbets + 1)))), tcdrs));
            end
            ptc = reshape(ptc, [numrows, tcds]);
            if ~isempty(oldord)
                ptc = permute(ptc, oldord);
            end
        end

    % split
    else

        % prepare outputs to receive data
        betas = zeros(numbets, numvox);
        if nargout > 1
            irtc = zeros([numbets, numbets, numvox]);
            if nargout > 2
                ptc = zeros(size(tcd));
                if nargout > 3
                    se = zeros(numvox, 1);
                end
            end
        end

        % compute size
        voxcount = 1;
        voxstep = max(1, floor(2e7 / (numbets * numrows)));
        pirtc = zeros([numbets, numbets, voxstep]);

        % loop
        while voxcount <= numvox

            % to-step
            voxend = min(numvox, voxcount + voxstep - 1);
            tnumvox = voxend + 1 - voxcount;

            % compute weighted design and data
            sw = sqrt(w(:, voxcount:voxend));
            ws = sum(w(:, voxcount:voxend), 1);
            wrtc = repmat(rtc, [1, 1, tnumvox]) .* ...
                repmat(reshape(sw, [numrows, 1, tnumvox]), [1, numbets]);

            % apply weights to data
            wtcd = sw .* tcd(:, voxcount:voxend);

            % calculate betas
            for vc = 1:tnumvox
                pirtc(:, :, vc) = inv(wrtc(:, :, vc)' * wrtc(:, :, vc));
            end
            if tnumvox < voxstep
                pirtc(:, :, tnumvox+1:end) = [];
            end
            pirtc(abs(pirtc) < tol) = 0;
            irtc(:, :, voxcount:voxend) = pirtc;
            betas(:, voxcount:voxend) = ...
                mtimesnd(mtimesnd(pirtc, wrtc, 2), reshape(wtcd, [numrows, 1, tnumvox]));

            % more output?
            if nargout > 2
                ptc(:, voxcount:voxend) = ...
                    rtc * reshape(betas(:, voxcount:voxend), [numbets, tnumvox]);
                if nargout > 3
                    sef = numrows ./ ws;
                    pse = wtcd - sw .* ptc(:, voxcount:voxend);
                    se(voxcount:voxend) = ...
                        lsqueeze(sef .* sqrt(max(0, sum(pse .* pse) ./ (ws - (numbets + 1)))));
                end
            end

            % increase counter
            voxcount = voxcount + voxstep;
        end

        % reshape outputs
        if nargout > 2
            ptc = reshape(ptc, [numrows, tcds]);
            if ~isempty(oldord)
                ptc = permute(ptc, oldord);
            end
            if nargout > 3
                se = reshape(se, [tcds, 1]);
            end
        end
    end
end

% reshape betas
betas = reshape(shiftdim(betas, 1), [tcds, numbets]);
