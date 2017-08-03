function o = poolgroups(g, r, opts)
% poolgroups  - pool several sample groups into one stats analysis
%
% FORMAT:       [o = ] poolgroups(g, r, opts)
%
% Input fields:
%
%       g           1xG cell array with groups (images/data)
%       r           regression data/model, can be either
%                   - empty or column of ones (for simple t-test)
%                   - a SxP subject-by-predictor model
%                   - a 1xG cell array each with a SxP sub model
%                     (indicating that slope-difference columns are used)
%       opts        optional settings
%        .contrasts 1xC struct array with fields
%         .outname  output filename (only used for images)
%         .type     either 't' or 'F'
%         .weights  1xG weights (only required for regression model)
%        .rank      flag, rank-transform data (only for regression models)
%        .resbox    resampling box (only used for images)
%        .resgrid   resampling grid (only used for images)
%        .resimeth  resampling interpolation method (only used for images)
%        .robust    flag, use robust stats instead of OLS
%        .slopediff flag, force slope differences for SxP model
%        .splitgrp  flag, force split groups into separate predictors
%        .thresh    relative threshold for discarding stats (default: 0.75)
%
% Output fields:
%
%       o           contrasts in the size of the data

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:08 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% first arguments check
if nargin < 1 || ...
   ~iscell(g) || ...
    isempty(g)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
ng = numel(g);
rns = [];
if nargin < 2 || ...
    isempty(r) || ...
   (~isa(r, 'double') && ...
    ~iscell(r))
    r = [];
elseif iscell(r)
    if numel(r) ~= ng
        error( ...
            'neuroelf:BadArgument', ...
            'g and r arguments must match for cell r array.' ...
        );
    end
    r = r(:)';
    rns = zeros(1, ng);
    for gc = 1:ng
        if ~isa(r{gc}, 'double') || ...
            ndims(r{gc}) ~= 2 || ...
            size(r{gc}, 1) < 2 || ...
            isempty(r{gc}) || ...
            any(isinf(r{gc}(:)) | isnan(r{gc}(:))) || ...
            any(sum(abs(r{gc}), 1) == 0) || ...
            numel(find(sum(abs(diff(r{gc}))) == 0)) > 1
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid r argument.' ...
            );
        end
        if gc == 1
            np = size(r{gc}, 2);
        elseif np ~= size(r{gc}, 2)
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid r argument.' ...
            );
        end
        rns(gc) = size(r{gc}, 1);
    end
    r = cat(1, r{:});
else
    if any(isinf(r(:)) | isnan(r(:)))
        error( ...
            'neuroelf:BadArgument', ...
            'Model must not contain Inf/NaN values.' ...
        );
    end
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isempty(rns)
    opts.slopediff = true;
end
if ~isfield(opts, 'contrasts') || ...
   ~isstruct(opts.contrasts) || ...
    isempty(opts.contrasts) || ...
   ~isfield(opts.contrasts, 'type')
    opts.contrasts = emptystruct({'type', 'weights'});
end
if ~isfield(opts, 'rank') || ...
    numel(opts.rank) ~= 1 || ...
   ~islogical(opts.rank)
    opts.rank = false;
end
if ~isfield(opts, 'resbox') || ...
   ~isa(opts.resbox, 'double') || ...
   ~isequal(opts.resbox, [2, 3]) || ...
    any(isinf(opts.resbox(:)) | isnan(opts.resbox(:))) || ...
    any(diff(opts.resbox) < 0) || ...
    any(abs(opts.resbox) > 140)
    opts.resbox = [];
end
if ~isfield(opts, 'resgrid') || ...
   ~isa(opts.resgrid, 'double') || ...
    numel(opts.resgrid) ~= 3 || ...
    any(isinf(opts.resgrid) | isnan(opts.resgrid) | opts.resgrid == 0 | abs(opts.resgrid) > 8)
    opts.resgrid = [];
else
    opts.resgrid = opts.resgrid(:)';
end
if ~isfield(opts, 'resimeth') || ...
   ~ischar(opts.resimeth) || ...
   ~any(strcmpi(opts.resimeth(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.resimeth = 'linear';
end
if ~isfield(opts, 'robust') || ...
    numel(opts.robust) ~= 1 || ...
   ~islogical(opts.robust)
    opts.robust = false;
end
if ~isfield(opts, 'slopediff') || ...
    numel(opts.slopediff) ~= 1 || ...
   ~islogical(opts.slopediff)
    opts.slopediff = false;
end
if ~isfield(opts, 'splitgrp') || ...
    numel(opts.splitgrp) ~= 1 || ...
   ~islogical(opts.splitgrp)
    opts.splitgrp = false;
end
if ~isfield(opts, 'thresh') || ...
   ~isa(opts.thresh, 'double') || ...
    numel(opts.thresh) ~= 1 || ...
    isinf(opts.thresh) || ...
    isnan(opts.thresh) || ...
    opts.thresh < 0 || ...
    opts.thresh > 1
    opts.thresh = 0.75;
end

% check group contents
hdr = [];
is = [];
ts = [];
sa = [Inf, Inf, Inf; ones(3, 3)];
ns = zeros(1, ng);
for gc = 1:ng
    if (~isa(g{gc}, 'double') && ...
        ~iscell(g{gc})) || ...
        isempty(g{gc})
        clearxffobjects(g);
        error( ...
            'neuroelf:BadArgument', ...
            'Group data must not be empty.' ...
        );
    end
    if isa(g{gc}, 'double')
        tts = size(g{gc});
        ns(gc) = tts(end);
        tts(end) = [];
        if isempty(ts)
            if ~isempty(opts.resbox) && ~isempty(opts.resgrid)
                if numel(tts) ~= 3 || ...
                    any(tts ~= [ ...
                        numel(opts.resbox(1, 1):opts.resgrid(1):opts.resbox(2, 1)), ...
                        numel(opts.resbox(1, 2):opts.resgrid(2):opts.resbox(2, 2)), ...
                        numel(opts.resbox(1, 3):opts.resgrid(3):opts.resbox(2, 3))])
                    clearxffobjects(g);
                    error( ...
                        'neuroelf:BadArgument', ...
                        'Binary data size and sampling box mismatch.' ...
                    );
                end
            end
            ts = tts;
        else
            if ~isequal(ts, tts)
                clearxffobjects(g);
                error( ...
                    'neuroelf:BadArgument', ...
                    'Numeric data must have the same D-1 sizes.' ...
                );
            end
        end
    else
        tg = g{gc};
        for ic = 1:numel(tg)
            if ~ischar(tg{ic}) || ...
                isempty(tg{ic}) || ...
                exist(tg{ic}(:)', 'file') ~= 2
                clearxffobjects(g);
                error( ...
                    'neuroelf:BadArgument', ...
                    'File not found: ''%s''.', ...
                    tg{ic}(:)' ...
                );
            end
            try
                tg{ic} = xff(strrep(tg{ic}(:)', '.img', '.hdr'));
                if numel(tg{ic}) ~= 1 || ...
                   ~isxff(tg{ic}, 'hdr')
                    error('NO_HDR');
                end
            catch ne_eo;
                clearxffobjects(g);
                error( ...
                    'neuroelf:BadArgument', ...
                    'Error loading image file: %s.', ...
                    ne_eo.message ...
                );
            end
            if isempty(is)
                if isempty(ts)
                    if isempty(opts.resbox) || ...
                        isempty(opts.resgrid)
                        cfr = tg{ic}.CoordinateFrame.Trf;
                        cfr3 = cfr(1:3, 1:3);
                        if any(sum(sum(abs(cfr3 - diag(diag(cfr3))))))
                            clearxffobjects(g);
                            error( ...
                                'neuroelf:BadArgument', ...
                                'Images must be orthogonal and standard axes.' ...
                            );
                        end
                        is = size(tg{ic}.VoxelData);
                        opts.resgrid = lsqueeze(diag(cfr3))';
                        opts.resbox = cfr(1:3, 4)' + opts.resgrid;
                        opts.resbox(2, :) = ...
                            opts.resbox + opts.resgrid .* (is - 1);
                        ts = is;
                        hdr = tg{ic};
                    else
                        is = [ ...
                            numel(opts.resbox(1, 1):opts.resgrid(1):opts.resbox(2, 1)), ...
                            numel(opts.resbox(1, 2):opts.resgrid(2):opts.resbox(2, 2)), ...
                            numel(opts.resbox(1, 3):opts.resgrid(3):opts.resbox(2, 3))];
                        ts = is;
                    end
                else
                    if isempty(opts.resbox) || ...
                        isempty(opts.resgrid)
                        clearxffobjects(g);
                        error( ...
                            'neuroelf:BadArgument', ...
                            'Mixing data and images requires resampling.' ...
                        );
                    end
                    tts = [ ...
                        numel(opts.resbox(1, 1):opts.resgrid(1):opts.resbox(2, 1)), ...
                        numel(opts.resbox(1, 2):opts.resgrid(2):opts.resbox(2, 2)), ...
                        numel(opts.resbox(1, 3):opts.resgrid(3):opts.resbox(2, 3))];
                    if ~isequal(tts, ts)
                        clearxffobjects(g);
                        error( ...
                            'neuroelf:BadArgument', ...
                            'Mixing data and images requires resampling match.' ...
                        );
                    end
                end
            end
            ns(gc) = numel(g{gc});
        end
        g{gc} = tg;
    end
end
if ~isempty(rns) && ...
    any(rns ~= ns)
    clearxffobjects(g);
    error( ...
        'neuroelf:BadArgument', ...
        'Group sizes mismatch.' ...
    );
end
nsf = 1 + [0, cumsum(ns)];
nst = cumsum(ns);
if opts.rank
    for gc = 1:ng
        r(nsf(gc):nst(gc), :) = ranktrans(r(nsf(gc):nst(gc), :), 1, ...
            struct('meancenter', true, 'nozero', true));
    end
elseif ~isempty(r)
    for gc = 1:ng
        r(nsf(gc):nst(gc), :) = ztrans(r(nsf(gc):nst(gc), :), 1);
    end
end
if ~isempty(is)
    sa([2, 4], :) = opts.resbox;
    sa(3, :) = opts.resgrid;
    if prod(is) >= 1000 || ...
        is(end) > 2
        saz = sa(2, 3):sa(3, 3):sa(4, 3);
        sa(4, 3) = sa(2, 3);
    end
    so = struct('method', opts.resimeth);
else
    saz = [];
end
if isempty(ts)
    ts = is;
end
if isempty(saz)
    saz = 1:ts(end);
end
tns = sum(ns);

% initialize output
o = zeros([ts, numel(opts.contrasts)]);
osa = struct('type', '()', 'subs', {{':'}});
osa.subs = osa.subs(ones(1, numel(ts) + 1));
dim = numel(ts) + 1;

% create mean differences model part
mp = zeros(tns, numel(g));
if ~opts.splitgrp && ...
    opts.slopediff
    mp(:, 1) = 1;
    for gc = 1:(ng - 1)
        mp(nsf(gc):nst(gc), gc + 1) = 1;
        mp(nsf(gc+1):nst(gc+1), gc + 1) = -1;
    end
else
    for gc = 1:ng
        mp(nsf(gc):nst(gc), gc) = 1;
    end
end
nsb = nsf - 1;

% remove constant from model
if ~isempty(r)
    const = (sum(abs(diff(r))) == 0);
    r(:, const) = [];
    np = size(r, 2) + 1;
else
    r = zeros(tns, 0);
    np = 1;
end

% create model
if opts.splitgrp || ...
    opts.slopediff
    m = zeros(tns, np * ng);
    for pc = 1:(np - 1)
        m(:, ((pc - 1) * ng + 1):pc*ng) = r(:, pc(ones(1, ng))) .* mp;
    end
    m(:, ((np - 1) * ng + 1):np*ng) = mp;
else
    m = [r, mp];
end
df = size(m, 1) - size(m, 2);
iXX = pinv(m' * m);

% check contrasts
if isempty(opts.contrasts)
    if opts.slopediff
        for pc = 1:np
            opts.contrasts(2*pc - 1).type = 't';
            opts.contrasts(2*pc - 1).weights = ...
                [zeros(1, (pc-1) * ng), 1, zeros(1, (1+np-pc) * ng - 1)];
            opts.contrasts(2*pc).type = 'F';
            opts.contrasts(2*pc).weights = ...
                [zeros(1, (pc-1) * ng), 0, ones(1, ng - 1), zeros(1, (np-pc) * ng)];
        end
    else
        for pc = 1:np
            opts.contrasts(pc).type = 't';
            opts.contrasts(pc).weights = ...
                [zeros(1, (pc-1) * ng), ones(1, ng), zeros(1, (np-pc) * ng)];
        end
    end
end
if ~isfield(opts.contrasts, 'outname')
    opts.contrasts(1).outname = '';
end
for cc = 1:numel(opts.contrasts)
    if ~ischar(opts.contrasts(cc).outname) || ...
        isempty(opts.contrasts(cc).outname)
        opts.contrasts(cc).outname = sprintf('poolgroups_%s%s.hdr', ...
            opts.contrasts(cc).type, sprintf('_%d', ...
            opts.contrasts(cc).weights));
    end
end

% for very small sample sizes
if prod(ts) < 1000 || ...
    ts(end) < 3

    % sample all images
    for gc = 1:gc
        if iscell(g{gc});
            for ec = 1:numel(g{gc})
                g{gc}{ec} = g{gc}{ec}.SampleData3D(sa, so);
            end
            g{gc} = cat(4, g{gc}{:});
        end
    end

    % rank ?
    if opts.rank
        for gc = 1:ng
            g{gc} = ranktrans(g{gc}, dim, ...
                struct('meancenter', true, 'nozero', true));
        end
    end
    g = cat(dim, g{:});
    g(isinf(g) | isnan(g)) = 0;
    gmsk = double(sum(g ~= 0, dim, 'double') >= (opts.thresh * tns));

    % robust ?
    if opts.robust
        [b, w] = fitrobustbisquare_img(m, g);
        for cc = 1:numel(opts.contrasts)
            osa.subs{end} = cc;
            if opts.contrasts(cc).type == 't'
                o = subsasgn(o, osa, ...
                    gmsk .* robustt(m, g, b, w, opts.contrasts(cc).weights));
            else
                o = subsasgn(o, osa, ...
                    gmsk .* robustF(m, g, b, w, opts.contrasts(cc).weights));
            end
        end

    % or OLS
    else
        g = permute(g, [dim, 1:numel(ts)]);
        g = reshape(g, tns, prod(ts));
        b = iXX * m' * g;
        rss = g - m * b;
        rss = sum(rss .* rss);
        for cc = 1:numel(opts.contrasts)
            cw = opts.contrasts(cc).weights;
            osa.subs{end} = cc;
            if opts.contrasts(cc).type == 't'
                oc = (cw * b) ./ sqrt(((cw * iXX * cw') / df) * rss);
            else
                bi = (cw ~= 0);
                rssf = sum((m(:, bi) * b(bi, :)) .^ 2);
                oc = (df / sum(bi)) * (rssf ./ rss);
            end
            oc(isinf(oc) | isnan(oc)) = 0;
            o = subsasgn(o, osa, ...
                gmsk .* permute(reshape(oc, [tns, ts(1:dim-1)]), [2:dim, 1]));
        end
    end

% larger sample
else

    % perform computation
    bsa = osa;
    bsa.subs{end-1} = 1;
    tg = ts;
    tg(end) = 1;
    tg(end + 1) = tns;
    gv = zeros(tg);

    % progress bar over slices
    try
        pbar = xprogress;
        xprogress(pbar, 'setposition', [80, 260, 640, 36]);
        xprogress(pbar, 'settitle', 'Pooling group stats...');
        xprogress(pbar, 0, ...
            sprintf('Slices completed: 0 of %d', ts(end)), 'visible', 0, ts(end));
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        pbar = [];
    end

    % iterate over last "image" dim
    for ldc = 1:ts(end)
        sa([2, 4], 3) = saz(ldc);
        bsa.subs{end-1} = ldc;

        % prepare osa and get input (with ranktrans if needed)
        osa.subs{end-1} = 1;
        for gc = 1:ng
            if isa(g{gc}, 'double')
                osa.subs{end} = nsf(gc):nst(gc);
                if opts.rank
                    gv = subsasgn(gv, osa, ...
                        ranktrans(subsref(g{gc}, bsa), dim, ...
                        struct('meancenter', true, 'nozero', true)));
                else
                    gv = subsasgn(gv, osa, subsref(g{gc}, bsa));
                end
            else
                for ec = 1:numel(g{gc})
                    osa.subs{end} = nsb(gc) + ec;
                    gv = subsasgn(gv, osa, g{gc}{ec}.SampleData3D(sa, so));
                end
                if opts.rank
                    osa.subs{end} = nsf(gc):nst(gc);
                    gv = subsasgn(gv, osa, ...
                        ranktrans(subsref(gv, osa), dim, ...
                        struct('meancenter', true, 'nozero', true)));
                end
            end
        end
        gv(isinf(gv) | isnan(gv)) = 0;
        gvmsk = double(sum(gv ~= 0, dim, 'double') >= (opts.thresh * tns));

        % perform computation
        osa.subs{end-1} = ldc;
        if opts.robust
            [b, w] = fitrobustbisquare_img(m, gv);
            for cc = 1:numel(opts.contrasts)
                osa.subs{end} = cc;
                if opts.contrasts(cc).type == 't'
                    o = subsasgn(o, osa, ...
                        gvmsk .* robustt(m, gv, b, w, opts.contrasts(cc).weights));
                else
                    o = subsasgn(o, osa, ...
                        gvmsk .* robustF(m, gv, b, w, opts.contrasts(cc).weights));
                end
            end

        % OLS case
        else
            gvr = permute(gv, [dim, 1:numel(ts)]);
            gvr = reshape(gvr, tns, prod(ts(1:dim-2)));
            b = iXX * m' * gvr;
            rss = gvr - m * b;
            rss = sum(rss .* rss);
            for cc = 1:numel(opts.contrasts)
                cw = opts.contrasts(cc).weights;
                osa.subs{end} = cc;
                if opts.contrasts(cc).type == 't'
                    oc = (cw * b) ./ sqrt(((cw * iXX * cw') / df) * rss);
                else
                    bi = (cw ~= 0);
                    rssf = sum((m(:, bi) * b(bi, :)) .^ 2);
                    oc = (df / sum(bi)) * (rssf ./ rss);
                end
                oc(isinf(oc) | isnan(oc)) = 0;
                o = subsasgn(o, osa, ...
                    gvmsk .* permute(reshape(oc, [1, ts(1:dim-2)]), [2:dim, 1]));
            end
        end

        % update progress bar
        if ~isempty(pbar)
            xprogress(pbar, ldc, ...
                sprintf('Slices completed: %d of %d', ldc, ts(end)));
        end
    end
end

% save to images?
if ~isempty(hdr)
    hdr.ImgDim.DataType = 16;
    hdr.ImgDim.BitsPerPixel = 32;
    for cc = 1:numel(opts.contrasts)
        hdr.VoxelData = single(o(:, :, :, cc));
        hdr.SaveAs(opts.contrasts(cc).outname);
    end
end

% clear progress bar
if ~isempty(pbar)
    closebar(pbar);
end

% clear objects
for gc = 1:ng
    if isa(g{gc}, 'double')
        g{gc} = [];
    end
end
clearxffobjects(g);
