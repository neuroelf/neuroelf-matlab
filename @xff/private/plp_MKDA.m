function mvmp = plp_MKDA(xo, opts)
% PLP::MKDA  - perform a multi-level kernel density analysis
%
% FORMAT:       mvmp = plp.MKDA([opts])
%
% Input fields:
%
%       opts        1x1 struct with options
%        .applymask apply mask to output results (default: false)
%        .asimiter  number of (false-positive) simulation iterations (5000)
%        .asimkeep  boolean flag, keep (average of) simulated maps (false)
%        .asimkthr  boolean flag, estimate cluster-based thresholds (true)
%        .asimmask  mask object (HDR, MSK or VMR, default: colin brain)
%        .asimsmpl  random coordinates sampling, either {'full'} or 'near'
%        .asimrthr  raw thresholds for cluster thresholds (default:
%                   [0.95, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9999,
%                    0.99995, 0.99999, 0.999995, 0.999999])
%        .bbox      bounding box (passed into newnatresvmp, default: [])
%        .cond      conditional statement, e.g.
%                   '$Study >= 1 & $Study <= 3 & $Type == 2'
%                   whereas $names are replaced by their column data
%                   and syntax as in '$Column == ''content''' is translated
%                   into a call to strcmpi (or regexpi for patterns)
%        .contcomp  contrast computation, 'diff', 'excl', or {'wexcl'}
%        .contexclw exclusion weight (default 0.5)
%        .contnames cell array with names for contrast terms
%        .contnull  null-distribution of contrasts either of
%                    'full'      - shuffle assignments overall and within
%                    'labels'    - shuffle labels within studies only
%                   {'spatial'}  - construct null distribution by random
%                                  spatial resampling within .asimmask
%                    'units'     - shuffle labels across studies/units only
%        .contrasts contrast definition cell array, default: {[1]}
%        .fixcoords fix coordinates to nearest integer (VMP space, false)
%        .grpmeth   group statistics method, one of 'ost', 'sum', {'wsum'}
%        .indivmaps keep individual (study-specific) maps (default: true)
%        .inputmask apply mask to coordinate input (PLP content, true)
%        .jbmeth    join-blobs method, either of 'max', {'rsum'}
%        .pbar      progress bar object (either xfigure or xprogress, [])
%        .res       VMP resolution (passed into newnatresvmp, default: 3)
%        .scale     scaling flag, either of 'indic', {'toone'}
%        .smkern    single or multiple smoothing kernels in mm (default: 8)
%        .smkinterp smoothing kernel interpolation (default: linear)
%        .smkmdist  distance from peak where the value is max (default: 0)
%        .smkres    smoothing kernel resolution (default: 1)
%        .sparse    perform analysis with sparse data (default: false)
%        .studycol  column name for study (statistical unit)
%        .studysel  selection of studies (default: all)
%        .stwf      per-study weight formula, e.g.
%                   'sqrt($GroupSize) .* (0.75 + 0.25 .* ($RFX==1))'
%                   whereas names are replaced by their column data
%        .stwmcc    correct for multiple maps per study (default: false)
%        .stwp      relative points-per-study weight for testing, one of
%                    'confidence'  - weigh 0.5 + gammapdf((1:maxp)./meanp)
%                    'logpoints'   - weigh by 1 + log(nr of points)
%                   {'none'}       - assume all studies contribue equally
%                    'points'      - weigh by number of points per study
%                    'sqrtpoints'  - weigh by sqrt(nr of points)
%        .stwrdiff  relative weights in diff. contrasts (default: true)
%        .stwsel    select only studies with points (default: true)
%        .unique    only use unique points (within study, default: true)
%        .usecons   flag, if given must be either the name or number of
%                   the column containing the contrast identifier
%        .usesize   flag, if given must be either the name or number of
%                   the column containing the (original) cluster size
%        .usevalue  flag, if given must be either the name of number of
%                   the column containing the (original) peak value
%
% Output fields:
%
%       mvmp        VMP container with one (set of) map(s) for each contrast
%
% Notes: all points must be in the same coordinate space (so any conversion
%        should occur prior to storing the coordinates in the PLP object)
%
%        the .asimsmpl parameter either samples the coordinates from all
%        voxels in the .asimmask or, alternatively, it adds a displacement
%        of between 0.5 and 1.5 times smpkern to each coordinate, checking they
%        remain within the mask
%
%        the .contcomp parameter defined whether contrast (within study)
%        are a straight difference between the blob maps ('diff') or
%        whether if both values of the term are ~= 0, a zero values is used
%
%        the non-spatial .contnull settings ('full' and 'labels') make
%        most sense when used with the rescaled z-Map, as FWE results are
%        likely to be too stringent
%
%        the .contrasts cell array should contain cells that give the
%        identifiers for the positive and negative terms in the contrast,
%        e.g. in a meta analysis of three contrasts, this list could be
%        {3, [3, -1], [3, -2]}, which would then compute three contrasts:
%        3 (on its own, i.e. vs. baseline), 3 > 1, and 3 > 2
%
%        the .scale parameter either creates a smooth gaussian kernel
%        around the peak, fixing the ceiling value to 1 ('toone'), or sets
%        voxels above .5 (after toone scaling) to 1 (indicator);
%        if .scale is set to 'indic' (pure indicator), the kernel size is
%        interpreted as a radius!
%
%        the .smkmdist parameter can be used to have the initial blobs be
%        either complete gaussian blobs (0, default) or flattened in a way
%        such that within a specific range the value is actually the max
%        value, after which it drops off using the gaussian smoothing;
%        this allows to create "smooth" but still "solid" blobs, e.g. by
%        setting the smkmdist to the same value(s) as smkern
%
%        the .smkres parameter can be set to a smaller resolution in case
%        the coordinates (plp.Points) contain non-rounded entries (e.g.
%        after MNI->TAL or vice versa transformation)
%
%        the .stwf (per-study weighting formula) field will default to '1'
%        unless a valid 'GroupSize' or 'N' column is found, in which case
%        the .stwf field will be set to 'sqrt($GroupSize)' or 'sqrt($N)'
%
%        the .stwp parameter sets an additional weight on studies based
%        on how many points are available for the study (at any given
%        contrast iteration); the 'confidence' setting tries to determine
%        a "sweet spot", the number of points that a "good" (average)
%        study should report, and, at most, reduces the weight by 0.5 for
%        studies that have either no points or the maximum number of points
%
% Using: bssample, clustercoordsc, dilate3d, findfirst, flexinterpn,
%        flexinterpn_method, gammapdf, histcount, limitrangec, lsqueeze,
%        minmaxmean, newnatresvmp, onesamplet, sdist, smoothkern,
%        splittocellc, varc.

% Version:  v1.1
% Build:    16040817
% Date:     Apr-08 2016, 5:05 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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

% import from neuroelf library into workspace
using(neuroelf, {'bssample', 'clustercoordsc', 'dilate3d', 'findfirst', ...
    'flexinterpn', 'flexinterpn_method', 'gammapdf', 'histcount', ...
    'limitrangec', 'lsqueeze', 'minmaxmean', 'newnatresvmp', 'onesamplet', ...
    'sdist', 'smoothkern', 'splittocellc', 'varc'});

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'plp')
    error('neuroelf:xff:BadArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
bcf = xo.F;
if isempty(bcf) && isfield(bc.RunTimeVars, 'SourceFile') && ...
   ~isempty(bc.RunTimeVars.SourceFile)
    bcf = bc.RunTimeVars.SourceFile;
end
cn = bc.ColumnNames(:);
labels = bc.Labels(:);
p = bc.Points;
if numel(cn) ~= size(p, 2)
    error('neuroelf:xff:badObject', ...
        'PLP object has an invalid ColumnNames/Points combination. Please fix.');
end
xc = findfirst(strcmpi(cn, 'x'));
if isempty(xc)
    xc = 1;
end
yc = findfirst(strcmpi(cn, 'y'));
if isempty(xc)
    yc = 2;
end
zc = findfirst(strcmpi(cn, 'z'));
if isempty(xc)
    zc = 3;
end
xyz = [xc, yc, zc];
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'applymask') || ~islogical(opts.applymask) || numel(opts.applymask) ~= 1
    opts.applymask = false;
end
if ~isfield(opts, 'asimiter') || ~isa(opts.asimiter, 'double') || numel(opts.asimiter) ~= 1 || ...
    isinf(opts.asimiter) || isnan(opts.asimiter) || opts.asimiter < 1
    asit = 5000;
else
    asit = min(1000000, max(1, opts.asimiter));
end
if ~isfield(opts, 'asimkeep') || ~islogical(opts.asimkeep) || numel(opts.asimkeep) ~= 1
    opts.asimkeep = 0;
else
    opts.asimkeep = double(opts.asimkeep);
end
if ~isfield(opts, 'asimkthr') || ~islogical(opts.asimkthr) || numel(opts.asimkthr) ~= 1
    opts.asimkthr = true;
end
if ~isfield(opts, 'asimmask') || numel(opts.asimmask) ~= 1 || ...
   ~xffisobject(opts.asimmask, true, {'hdr', 'msk', 'vmr'})
    opts.asimmask = [];
end
if ~isfield(opts, 'asimrthr') || ~isa(opts.asimrthr, 'double') || ...
    any(isinf(opts.asimrthr(:)) | isnan(opts.asimrthr(:)) | ...
        opts.asimrthr(:) < 0.5 | opts.asimrthr(:) >= 1)
    opts.asimrthr = [0.95, 0.98, 0.99, 0.995, 0.998, 0.999, ...
        0.9995, 0.9999, 0.99995, 0.99999, 0.999995, 0.999999];
else
    opts.asimrthr = opts.asimrthr(:)';
end
if ~isfield(opts, 'asimsmpl') || ~ischar(opts.asimsmpl) || ...
   ~any(strcmpi(opts.asimsmpl(:)', {'f', 'full', 'n', 'near'}))
    opts.asimsmpl = 'f';
else
    opts.asimsmpl = lower(opts.asimsmpl(1));
end
if ~isfield(opts, 'bbox') || ~isa(opts.bbox, 'double') || ~isequal(size(opts.bbox), [2, 3])
    opts.bbox = [];
end
if ~isfield(opts, 'cond') || ~ischar(opts.cond)
    opts.cond = '';
    ocond = '';
else
    opts.cond = ['(' lower(opts.cond(:)') ')'];
    ocond = opts.cond;
    cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
    while ~isempty(cregx) && ~isempty(cregx{1})
        cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*==\s*''([^'']+)''$', 'tokens');
        if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
        end
        if any(cregp{1}{2} == '*')
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                '~cellfun(''isempty'', regexpi(labels($%s), ''%s''))', cregp{1}{:}));
        else
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                'strcmpi(labels($%s), ''%s'')', cregp{1}{:}));
        end
        cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*==\s*''[^'']+'')', 'tokens');
    end
    cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
    while ~isempty(cregx) && ~isempty(cregx{1})
        cregp = regexp(cregx{1}{1}, '^\$([a-zA-Z][a-zA-Z_0-9]*)\s*~=\s*''([^'']+)''$', 'tokens');
        if ~iscell(cregp) || numel(cregp) ~= 1 || ~iscell(cregp{1}) || numel(cregp{1}) ~= 2
            error('neuroelf:xff:badArgument', 'Invalid conditional statement.');
        end
        if any(cregp{1}{2} == '*')
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                'cellfun(''isempty'', regexpi(labels($%s), ''%s''))', cregp{1}{:}));
        else
            opts.cond = strrep(opts.cond, cregx{1}{1}, sprintf( ...
                '~strcmpi(labels($%s), ''%s'')', cregp{1}{:}));
        end
        cregx = regexp(opts.cond, '(\$[a-zA-Z][a-zA-Z_0-9]*\s*~=\s*''[^'']+'')', 'tokens');
    end
end
cnl = zeros(numel(cn), 1);
for cnc = 1:numel(cn)
    cnl(cnc) = numel(cn{cnc});
end
[cnl, cni] = sort(cnl, 'descend');
cns = cn(cni);
for cnc = 1:numel(cn)
    opts.cond = strrep(opts.cond, ['$' lower(cns{cnc})], sprintf('p(:,%d)', cni(cnc)));
end
if isempty(opts.cond)
    opts.cond = '(true)';
end
if ~isfield(opts, 'contcomp') || ~ischar(opts.contcomp) || ...
   ~any(strcmpi(opts.contcomp(:)', {'d', 'diff', 'e', 'excl', 'w', 'wexl'}))
    opts.contcomp = 'w';
else
    opts.contcomp = lower(opts.contcomp(1));
end
if ~isfield(opts, 'contexclw') || ~isa(opts.contexclw, 'double') || numel(opts.contexclw) ~= 1 || ...
    isinf(opts.contexclw) || isnan(opts.contexclw)
    opts.contexclw = 0.5;
else
    opts.contexclw = max(sqrt(eps), min(1 - sqrt(eps), opts.contexclw));
end
if ~isfield(opts, 'contmode') || ~ischar(opts.contmode) || ...
   ~any(strcmpi(opts.contmode(:)', {'c', 'conj', 'm', 'mix', 'w', 'weight'}))
    opts.contmode = 'c';
else
    opts.contmode = lower(opts.contmode(1));
end
if ~isfield(opts, 'contnames') || ~iscell(opts.contnames) || isempty(opts.contnames) || ...
   ~ischar(opts.contnames{1}) || isempty(opts.contnames{1})
    opts.contnames = {};
else
    opts.contnames = opts.contnames(:);
    for cc = 1:numel(opts.contnames)
        if ~ischar(opts.contnames{cc}) || isempty(opts.contnames{cc})
            error('neuroelf:xff:badArgument', ...
                'Invalid or missing contrast name particle.');
        end
        opts.contnames{cc} = opts.contnames{cc}(:)';
    end
    if numel(unique(opts.contnames)) ~= numel(opts.contnames)
        error('neuroelf:xff:badArgument', 'Contrast name particles must be unique.');
    end
end
if ~isfield(opts, 'contnull') || ~ischar(opts.contnull) || isempty(opts.contnull) || ...
   ~any(lower(opts.contnull(1)) == 'flsu')
    opts.contnull = 's';
else
    opts.contnull = lower(opts.contnull(1));
end
if ~isfield(opts, 'contrasts') || ~iscell(opts.contrasts) || isempty(opts.contrasts)
    opts.contrasts = {1};
    if ~isempty(opts.contnames)
        opts.contnames = opts.contnames(1);
    else
        opts.contnames = {'contrast'};
    end
else
    opts.contrasts = opts.contrasts(:);
end
mxcon = 1;
pcon = [];
if isfield(opts, 'usecons') && isa(opts.usecons, 'double') && numel(opts.usecons) == 1 && ...
   ~isinf(opts.usecons) && ~isnan(opts.usecons) && any((1:numel(cn)) == opts.usecons)
    cvalues = unique(p(:, opts.usecons));
    cvalues(cvalues < 1 | cvalues > numel(labels) | cvalues ~= fix(cvalues)) = [];
    cstrings = labels(cvalues);
else
    cvalues = [];
end
for cc = 1:numel(opts.contrasts)
    if ~isempty(opts.contrasts{cc}) && ischar(opts.contrasts{cc}) && ~isempty(cvalues)
        cstring = opts.contrasts{cc}(:)';
        if sum(cstring == '>') > 1
            error('neuroelf:xff:badArgument', 'Bad contrast specification.');
        end
        cstring = splittocellc(cstring, '>');
        if numel(cstring) > 1
            cstringn = cstring{2};
        else
            cstringn = '';
        end
        cstring = cstring{1};
        cstring(cstring == '+' | cstring == ',' | cstring == ';' | cstring == '&') = ' ';
        cstringn(cstringn == '+' | cstringn == ',' | cstringn == ';' | cstringn == '&') = ' ';
        cstring = splittocellc(cstring, ' ', true);
        cstring(cellfun('isempty', cstring)) = [];
        cstringn = splittocellc(cstringn, ' ', true);
        cstringn(cellfun('isempty', cstringn)) = [];
        cvector = zeros(1, numel(cstring) + numel(cstringn));
        for cvc = 1:numel(cstring)
            cvci = find(strcmpi(cstrings, cstring{cvc}));
            if isempty(cvci) && all(cstring{cvc} >= '0' & cstring{cvc} <= '9')
                cvci = str2double(cstring{cvc});
            end
            if numel(cvci) ~= 1
                error('neuroelf:xff:badArgument', ...
                    'Bad contrast specification (term %s).', cstring{cvc});
            end
            if cvci > numel(cvalues) && any(cvci == cvalues)
                cvector(cvc) = cvci;
            elseif cvci > numel(cvalues)
                error('neuroelf:xff:badArgument', ...
                    'Bad contrast specification (term %s).', cstringn{cvc});
            else
                cvector(cvc) = cvalues(cvci);
            end
        end
        for cvc = 1:numel(cstringn)
            cvci = find(strcmpi(cstrings, cstringn{cvc}));
            if isempty(cvci) && all(cstringn{cvc} >= '0' & cstringn{cvc} <= '9')
                cvci = str2double(cstringn{cvc});
            end
            if numel(cvci) ~= 1
                error('neuroelf:xff:badArgument', ...
                    'Bad contrast specification (term %s).', cstringn{cvc});
            end
            if cvci > numel(cvalues) && any(cvci == cvalues)
                cvector(numel(cstring) + cvc) = -cvci;
            elseif cvci > numel(cvalues)
                error('neuroelf:xff:badArgument', ...
                    'Bad contrast specification (term %s).', cstringn{cvc});
            else
                cvector(numel(cstring) + cvc) = -cvalues(cvci);
            end
        end
        opts.contrasts{cc} = cvector;
    end
    if isempty(opts.contrasts{cc}) || ~isa(opts.contrasts{cc}, 'double') || ...
        numel(opts.contrasts{cc}) ~= size(opts.contrasts{cc}, 2) || ...
        any(isinf(opts.contrasts{cc}) | isnan(opts.contrasts{cc}) | ...
        opts.contrasts{cc} == 0 | opts.contrasts{cc} ~= fix(opts.contrasts{cc}))
        error('neuroelf:xff:BadArgument', 'Bad contrast specification.');
    end
    mxcon = max(mxcon, max(abs(opts.contrasts{cc})));
    pcon = unique([pcon(:); abs(opts.contrasts{cc}(:))]);
end
if numel(opts.contnames) ~= numel(opts.contrasts)
    contnames = cell(1, numel(opts.contrasts));
    for cc = 1:numel(contnames)
        cdef = opts.contrasts{cc};
        contnames{cc} = sprintf('%s + ', labels{cdef(cdef > 0)});
        contnames{cc}(end-2:end) = [];
        if any(cdef < 0)
            if sum(cdef > 0) > 1
                contnames{cc} = ['(' contnames{cc} ')'];
            end
            if sum(cdef < 0) == 1
                cminus = labels{-cdef(cdef < 0)};
            else
                cminus = ['(' sprintf('%s + ', labels{-cdef(cdef < 0)}) ')'];
                cminus(end-3:end-1) = [];
            end
            contnames{cc} = [ contnames{cc} ' > ' cminus];
        end
    end
else
    contnames = opts.contnames(:)';
end
if ~isfield(opts, 'fixcoords') || ~islogical(opts.fixcoords) || numel(opts.fixcoords) ~= 1
    opts.fixcoords = false;
end
if ~isfield(opts, 'grpmeth') || ~ischar(opts.grpmeth) || ...
   ~any(strcmpi(opts.grpmeth(:)', {'o', 'ost', 's', 'sum', 'w', 'wsum'}))
    opts.grpmeth = 'w';
else
    opts.grpmeth = lower(opts.grpmeth(1));
end
if ~isfield(opts, 'indivmaps') || ~islogical(opts.indivmaps) || numel(opts.indivmaps) ~= 1
    opts.indivmaps = true;
end
if ~isfield(opts, 'inputmask') || ~islogical(opts.inputmask) || numel(opts.inputmask) ~= 1
    opts.inputmask = true;
end
if ~isfield(opts, 'jbmeth') || ~ischar(opts.jbmeth) || ...
   ~any(strcmpi(opts.jbmeth(:)', {'m', 'max', 'r', 'rsum'}))
    opts.jbmeth = 'r';
else
    opts.jbmeth = lower(opts.jbmeth(1));
end
if ~isfield(opts, 'pbar') || numel(opts.pbar) ~= 1 || ...
   (~isa(opts.pbar, 'xfigure') && ~isa(opts.pbar, 'xprogress'))
    opts.pbar = [];
end
pb = opts.pbar;
if ~isfield(opts, 'res') || ~isa(opts.res, 'double') || numel(opts.res) ~= 1 || ...
    isinf(opts.res) || isnan(opts.res) || ~any(opts.res == (1:7))
    opts.res = 3;
end
if ~isfield(opts, 'scale') || ~ischar(opts.scale) || ...
   ~any(strcmpi(opts.scale(:)', {'1', 'i', 'indic', 'o', 'one', 't', 'toone'}))
    opts.scale = '1';
else
    opts.scale = lower(opts.scale(1));
end
if ~isfield(opts, 'smkern') || ~isa(opts.smkern, 'double') || isempty(opts.smkern) || ...
    any(isinf(opts.smkern(:)) | isnan(opts.smkern(:)) | opts.smkern(:) < 0)
    opts.smkern = 8;
else
    opts.smkern = opts.smkern(:)';
end
if ~isfield(opts, 'smkinterp') || ~ischar(opts.smkinterp) || ...
   ~any(strcmpi(opts.smkinterp(:)', {'cubic', 'lanczos3', 'linear', 'nearest'}))
    opts.smkinterp = 'linear';
else
    opts.smkinterp = lower(opts.smkinterp(:)');
end
[nsmk, opts.smkinterp] = flexinterpn_method(zeros(3, 1), [Inf; 1; 1; 3], opts.smkinterp);
if ~isfield(opts, 'smkmdist') || ~isa(opts.smkmdist, 'double') || isempty(opts.smkmdist) || ...
    any(isinf(opts.smkmdist(:)) | isnan(opts.smkmdist(:)) | opts.smkmdist(:) < 0)
    opts.smkmdist = 0;
end
if numel(opts.smkern) > 1 && numel(opts.smkmdist) > 1 && ...
    numel(opts.smkern) ~= numel(opts.smkmdist)
    error('neuroelf:xff:badArgument', 'Fields .smkern and .smkmdist must match in size.');
elseif numel(opts.smkern) > 1 && numel(opts.smkmdist) == 1
    opts.smkmdist = opts.smkmdist(1, ones(1, numel(opts.smkern)));
elseif numel(opts.smkmdist) > 1 && numel(opts.smkern) == 1
    opts.smkern = opts.smkern(1, ones(1, numel(opts.smkmdist)));
end
nsmk = numel(opts.smkern);
if ~isfield(opts, 'smkres') || ~isa(opts.smkres, 'double') || numel(opts.smkres) ~= 1 || ...
    isinf(opts.smkres) || isnan(opts.smkres)
    opts.smkres = 1;
else
    opts.smkres = max(0.25, min(1, opts.smkres));
end
resrat = opts.res ./ opts.smkres;
if ~isfield(opts, 'sparse') || ~islogical(opts.sparse) || numel(opts.sparse) ~= 1
    opts.sparse = false;
end
if ~isfield(opts, 'studycol') || ~ischar(opts.studycol) || isempty(opts.studycol) || ...
    isempty(regexpi(opts.studycol(:)', '^[a-z][a-z_0-9]+$'))
    opts.studycol = 'study';
else
    opts.studycol = lower(opts.studycol(:)');
end
if ~isfield(opts, 'studysel') || ~isa(opts.studysel, 'double') || isempty(opts.studysel) || ...
    any(isinf(opts.studysel(:)) | isnan(opts.studysel(:)) | ...
    opts.studysel(:) < 1 | opts.studysel(:) ~= fix(opts.studysel(:)))
    stcol = findfirst(strcmpi(cn, opts.studycol));
    if any(isinf(p(:, stcol)) | isnan(p(:, stcol)) | p(:, stcol) < 1 | p(:, stcol) ~= fix(p(:, stcol)))
        error('neuroelf:xff:badArgument', 'Study column not set correctly.');
    end
    opts.studysel = unique(p(:, stcol));
else
    if ~any(strcmpi(cn, opts.studycol))
        error('neuroelf:xff:badArgument', 'Study column required but not found.');
    end
    stcol = findfirst(strcmpi(cn, opts.studycol));
    if any(isinf(p(:, stcol)) | isnan(p(:, stcol)) | p(:, stcol) < 1 | ...
        p(:, stcol) ~= fix(p(:, stcol)))
        error('neuroelf:xff:badArgument', 'Study column not set correctly.');
    end
    opts.studysel = unique(opts.studysel(:));
end
if ~isfield(opts, 'stwf') || ~ischar(opts.stwf)
    opts.stwf = '';
else
    opts.stwf = ['(' lower(opts.stwf(:)') ')'];
end
for cnc = 1:numel(cn)
    opts.stwf = regexprep(opts.stwf, ['\$' lower(cns{cnc}) '\s*~=\s*''([^'']*)'''], ...
        sprintf('~strcmpi(labels(stp(:,%d)),''$1'')', cni(cnc)));
    opts.stwf = regexprep(opts.stwf, ['\$' lower(cns{cnc}) '\s*==\s*''([^'']*)'''], ...
        sprintf('strcmpi(labels(stp(:,%d)),''$1'')', cni(cnc)));
    opts.stwf = strrep(opts.stwf, ['$' lower(cns{cnc})], sprintf('stp(:,%d)', cni(cnc)));
end
if isempty(opts.stwf)
    opts.stwf = '';
end
if ~isfield(opts, 'stwmcc') || ~islogical(opts.stwmcc) || numel(opts.stwmcc) ~= 1
    opts.stwmcc = false;
end
if ~isfield(opts, 'stwp') || ~ischar(opts.stwp) || ~any(strcmpi(opts.stwp(:)', ...
    {'c', 'conf', 'confidence', 'l', 'log', 'logpoints', 'n', 'none', 'p', 'points', 's', 'sqrtpoints'}))
    opts.stwp = 'n';
else
    opts.stwp = lower(opts.stwp(1));
end
if ~isfield(opts, 'stwrdiff') || ~islogical(opts.stwrdiff) || numel(opts.stwrdiff) ~= 1
    opts.stwrdiff = true;
end
if ~isfield(opts, 'stwsel') || ~islogical(opts.stwsel) || numel(opts.stwsel) ~= 1
    opts.stwsel = true;
end
if ~isfield(opts, 'unique') || ~islogical(opts.unique) || numel(opts.unique) ~= 1
    opts.unique = true;
end
if ~isfield(opts, 'usecons') || ((~ischar(opts.usecons) || ~any(strcmpi(opts.usecons(:)', cn))) && ...
   (~isa(opts.usecons, 'double') || numel(opts.usecons) ~= 1 || isinf(opts.usecons) || ...
    isnan(opts.usecons) || opts.usecons < 1 || opts.usecons > numel(cn) || opts.usecons ~= fix(opts.usecons)))
    opts.usecons = findfirst(strcmpi(cn, 'contrast'));
elseif ischar(opts.usecons)
    opts.usecons = findfirst(strcmpi(opts.usecons(:)', cn));
end
if isempty(opts.usecons) && (numel(opts.contrasts) > 1 || ~isequal(opts.contrasts{1}, 1))
    error('neuroelf:xff:badArgument', ...
        'Computing contrasts (other than the default) requires the .usecons field.');
elseif isempty(opts.usecons)
    p(:, end + 1) = 1;
    opts.usecons = size(p, 2);
end
contcol = opts.usecons;
contcolm = max(p(:, contcol));
if ~isfield(opts, 'usesize') || ((~ischar(opts.usesize) || ~any(strcmpi(opts.usesize(:)', cn))) && ...
   (~isa(opts.usesize, 'double') || numel(opts.usesize) ~= 1 || isinf(opts.usesize) || ...
    isnan(opts.usesize) || opts.usesize < 1 || opts.usesize > numel(cn) || opts.usesize ~= fix(opts.usesize)))
    opts.usesize = findfirst(strcmpi(cn, 'size'));
    if isempty(opts.usesize)
        opts.usesize = findfirst(strcmpi(cn, 'clustersize'));
    end
elseif ischar(opts.usesize)
    opts.usesize = findfirst(strcmpi(opts.usesize(:)', cn));
end
if ~isfield(opts, 'usevalue') || ((~ischar(opts.usevalue) || ~any(strcmpi(opts.usevalue(:)', cn))) && ...
   (~isa(opts.usevalue, 'double') || numel(opts.usevalue) ~= 1 || isinf(opts.usevalue) || ...
    isnan(opts.usevalue) || opts.usevalue < 1 || opts.usevalue > numel(cn) || opts.usevalue ~= fix(opts.usevalue)))
    opts.usevalue = findfirst(strcmpi(cn, 'value'));
    if isempty(opts.usevalue)
        opts.usevalue = findfirst(strcmpi(cn, 'peakvalue'));
    end
elseif ischar(opts.usevalue)
    opts.usevalue = findfirst(strcmpi(opts.usevalue(:)', cn));
end
if opts.scale == 'i'
    opts.smkern = 2 .* opts.smkern;
end

% apply study selection
if ~isequal(opts.studysel(:), unique(p(:, stcol)))
    usest = false(max(max(opts.studysel), max(p(:, stcol))), 1);
    usest(opts.studysel) = true;
    p = p(usest(p(:, stcol)), :);
end

% apply input masking
if opts.inputmask && ~isempty(opts.asimmask)

    % remove coordinates with mask == 0 values
    p(aft_SampleData3D(opts.asimmask, p(:, xyz), struct('method', 'linear')) < 0.5, :) = [];
end

% then parse condition
if ~strcmp(opts.cond, '(true)')
    try
        pl = eval(opts.cond);
    catch xfferror
        error('neuroelf:xff:badArgument', 'Bad condition given: %s', xfferror.message);
    end
    p = p(pl, :);
end

% and then update study selection!
opts.studysel = unique(p(:, stcol));

% get coordinates (short hand)
pcrd = p(:, xyz);
ncrd = size(pcrd, 1);

% set study (pre-) weights to 1 and count unique points
numstudies = numel(opts.studysel);
stwf = ones(numstudies, 1);
stwu = zeros(numstudies, 1);
for stc = 1:numstudies
    if opts.unique
        stwu(stc) = size(unique(p(p(:, stcol) == opts.studysel(stc), xyz), 'rows'), 1);
    else
        stwu(stc) = size(p(p(:, stcol) == opts.studysel(stc), xyz), 1);
    end
end
mstwu = minmaxmean(stwu);
if opts.stwp == 'c'
    stwp = gammapdf((numel(pcon) / mstwu(3)) .* (0:mstwu(2))', 2, 1);
    stwp = (0.5 / max(stwp)) .* stwp;
end

% potentially replace the default with a better formula
if isempty(opts.stwf)
    if any(strcmpi(cn, 'groupsize')) && all(p(:, findfirst(strcmpi(cn, 'groupsize'))) > 0)
        opts.stwf = sprintf('sqrt(stp(:,%d))', findfirst(strcmpi(cn, 'groupsize')));
    elseif any(strcmpi(cn, 'n')) && all(p(:, findfirst(strcmpi(cn, 'n'))) > 0)
        opts.stwf = sprintf('sqrt(stp(:,%d))', findfirst(strcmpi(cn, 'n')));
    else
        opts.stwf = '1';
    end
end

% study weights need to be established
if ~strcmp(opts.stwf, '1')

    % evaluate preliminary study weights
    for stc = 1:numel(stwf)

        % first we find all points for that study
        stp = p(p(:, stcol) == opts.studysel(stc), :);

        % then try to compute formula
        try
            stpw = eval(opts.stwf);
        catch xfferror
            error('neuroelf:xff:badArgument', 'Bad study-weighting formula given: %s.', xfferror.message);
        end

        % the weight is then the median of those values
        stwf(stc) = median(stpw);
    end

    % remove bad studies
    badst = (isinf(stwf) | isnan(stwf) | stwf <= 0);
    if any(badst)
        opts.studysel(badst) = [];
        numstudies = numel(opts.studysel);
        usest = false(max(max(opts.studysel), max(p(:, stcol))), 1);
        usest(opts.studysel) = true;
        p = p(usest(p(:, stcol)), :);
        pcrd = p(:, xyz);
        ncrd = size(pcrd, 1);
        stwf(badst) = [];
        stwu(badst) = [];
        mstwu = minmaxmean(stwu);
        if opts.stwp == 'c'
            stwp = gammapdf((numel(pcon) / mstwu(3)) .* (0:mstwu(2))', 2, 1);
            stwp = 0.5 + (0.5 / max(stwp)) .* stwp;
        end
    end
end

% and extend for spatial transformation
pcrd(:, 4) = 1;

% and get rows that need to be shuffled
if opts.contnull == 'f'
    pconr = (p(:, contcol) > 0);
end

% get study labels
stlabels = cell(numstudies, 1);
if strcmpi(bc.ColumnNames{stcol}, 'study')
    for sc = 1:numel(stlabels)
        if numel(labels) < opts.studysel(sc)
            stlabels{sc} = sprintf('Study %d', opts.studysel(sc));
        else
            stlabels(sc) = labels(opts.studysel(sc));
        end
    end
else
    stlcol = findfirst(strcmpi(cn, 'study'));
    [stucolv, stucoli] = unique(p(:, stcol));
    stucoll = cn{stcol};
    for sc = 1:numel(stlabels)
        sci = stucoli(findfirst(stucolv == opts.studysel(sc)));
        stlabels{sc} = sprintf('%s %d (Study: %s)', stucoll, p(sci, stcol), labels{p(sci, stlcol)});
    end
end

% create vmp and get content, bounding box, and map size
mvmp = newnatresvmp(opts.bbox, opts.res, 16);
vmpb = aft_BoundingBox(mvmp);
vmpc = mvmp.C;
vmsz = size(vmpc.Map.VMPData);
nmsz = prod(vmsz);
smapp = zeros(vmsz);
smapn = zeros(vmsz);

% compute total number of maps
smaps = (numstudies + 2 + opts.asimkeep);
tmaps = smaps * numel(contnames) * nsmk;

% we need the transformation matrix
trf = vmpb.QuatT2B';

% to apply to coordinates (which are now in the maps 1-based coords)
pcrd = pcrd * trf;
pcrd(:, 4:end) = [];

% fix coordinates
if opts.fixcoords
    pcrd = round(pcrd);
end

% keep a backup
opcrd = pcrd;
if opts.contnull ~= 's'

    % and check for label across units shuffling
    if opts.contnull == 'u'

        % create array to store which study has positive/negative contrasts
        ullist = zeros(numstudies, numel(contnames));

        % iterate over contrasts
        for ctc = 1:numel(contnames)

            % get contrast particles
            cdef = opts.contrasts{ctc};
            cpos = cdef(cdef > 0);
            cneg = -cdef(cdef < 0);

            % must have positive and negative to work
            if isempty(cneg)
                error('neuroelf:xff:badArgument', ...
                    'Unit-based label shuffling requires +/- contrast.');
            end

            % units must have ONLY positive OR negative items
            for stc = 1:numstudies

                % find points for unit
                stp = find(p(:, stcol) == opts.studysel(stc));

                % ensure that points are either all positive or negative
                haspos = ~isempty(intersect(unique(p(stp, contcol)), cpos(:)));
                hasneg = ~isempty(intersect(unique(p(stp, contcol)), cneg(:)));

                % none at all or both is invalid
                if haspos && hasneg
                    error('neuroelf:xff:badArgument', ...
                        'Unit-based label shuffling requires no unit to have mixed +/- points.');
                elseif haspos
                    ullist(stc, ctc) = 1;
                elseif hasneg
                    ullist(stc, ctc) = -1;
                end
            end
        end
    end

    % then copy points
    opfull = p;
end

% interpolation grid
ipgrid = [Inf; 1; resrat; 1] * ones(1, 3);

% masking object
if isempty(opts.asimmask)
    try
        if exist([neuroelf_path '/_files/colin/colin_brain_ICBMnorm.vmr'], 'file') > 0
            asimmask = neuroelf_file('c', 'colin_brain_ICBMnorm.vmr');
        else
            asimmask = neuroelf_file('c', 'colin_brain.vmr');
        end
    catch xfferror
        neuroelf_lasterr(xfferror);
        delete(mvmp);
        error('neuroelf:xff:missingFile', ...
            'PLP::MetaVMP requires a gray-matter masking object.');
    end
else
    asimmask = opts.asimmask;
end
asmc = aft_SampleBVBox(asimmask, vmpb, 1, 'linear');

% then get mask (include one additional voxel)
asmm = (asmc > 0);
asmr = lsqueeze(dilate3d(asmm));

% and transformation matrix
asmt = vmpb.QuatB2T';

% for cluster-size based thresholding
aspdfac = 1 / 131072;
if opts.asimkthr

    % initialize cluster size arrays
    aspc = zeros(1, numel(opts.asimrthr));
    aspdist = zeros(1, 262145);
    ascsizes = zeros(asit, numel(opts.asimrthr));
end

% for voxel-wise rescaling, we need to keep ALL data for a real null
asmn = zeros(sum(asmr(:)), 1);

% clear any self-loaded mask
if isempty(opts.asimmask)

    % and then clear object
    delete(asimmask);
end

% masking coordinates (full set)
if opts.asimsmpl == 'f'

    % find voxels in mask
    [maskc, smc, stc] = ind2sub(size(asmc), find(asmc(:)));

    % compute them in the VMP space
    maskc = ([maskc, smc, stc, ones(numel(maskc), 1)] * asmt) * trf;
    maskc(:, 4:end) = [];

    % check that enough coordinates are available for sampling
    maskn = size(maskc, 1);
    if (maskn < 1000 && maskn < (0.25 * numel(asmc))) || any(varc(maskc) < max(0.5 .* opts.smkern))
        delete(mvmp);
        error('neuroelf:xff:badArgument', ...
            'Too few gray matter voxels in mask or insufficient spread.');
    end

% resample existing coordinates with replacement
else

    % compute transformation matrix from VMP to MASK
    atrf = inv(asmt * trf);
    asz = size(asmc);
    aszt = 0.5 + ones(ncrd, 1) * [asz, 1];

    % compute displacement values
    smpdisp = (1 / opts.res) .* opts.smkern;

    % and check that enough voxels are available
    if (sum(asmc(:)) / numel(asmc(:))) < 0.05 || sum(asmc(:)) < 1000
        warning('neuroelf:xff:badArgument', 'Very few gray matter voxels in mask.');
    end

    % for each smoothing kernel
    for smc = 1:nsmk

        % compute displacement values
        cdisp = smpdisp(smc) .* randn(ncrd, 3);
        dcdisp = opts.res * sqrt(sum(cdisp .* cdisp, 2));
        lcdisp = limitrangec(dcdisp, 0.5 * opts.smkern(smc), 2 * opts.smkern(smc), opts.smkern(smc));
        cdisp = cdisp .* ((lcdisp ./ dcdisp) * ones(1, 3));

        % check that for each coordinate a replacement can be found with
        pcrd = opcrd + cdisp;
        pcrd(:, 4) = 1;
        trcrd = round(pcrd * atrf);
        tpok = all(trcrd >= 0.5 & trcrd < aszt, 2);
        tpnok = find(~tpok);
        tpnok = [tpnok(:); lsqueeze(find(asmc(sub2ind(asz, ...
            trcrd(tpok, 1), trcrd(tpok, 2), trcrd(tpok, 3))) == 0))];
        tpiter = 15;
        while ~isempty(tpnok) && tpiter > 0
            cdisp = smpdisp(smc) * randn(numel(tpnok), 3);
            dcdisp = opts.res * sqrt(sum(cdisp .* cdisp, 2));
            lcdisp = limitrangec(dcdisp, ...
                0.5 * opts.smkern(smc), 2 * opts.smkern(smc), opts.smkern(smc));
            cdisp = cdisp .* ((lcdisp ./ dcdisp) * ones(1, 3));
            pcrd(tpnok, 1:3) = opcrd(tpnok, :) + cdisp;
            trcrd(tpnok, :) = round(pcrd(tpnok, :) * atrf);
            tpok = all(trcrd >= 0.5 & trcrd < aszt, 2);
            tpnok = find(~tpok);
            tpnok = [tpnok(:); lsqueeze(find(asmc(sub2ind(asz, ...
                trcrd(tpok, 1), trcrd(tpok, 2), trcrd(tpok, 3))) == 0))];
            tpiter = tpiter - 1;
        end
        if ~isempty(tpnok)
            trcrd = round([opcrd(tpnok, :), ones(numel(tpnok), 1)] * inv(trf));
            tpncrd = sprintf('#%5d: %4d,%4d,%4d\n', [tpnok(:)'; trcrd(:, 1:3)']);
            delete(mvmp);
            error('neuroelf:xff:badArgument', ...
                'Cannot find replacement coordinates for points within mask:\n%s', tpncrd);
        end
    end
end

% total number of iterations
asml = asit + 1;

% for FWE (map-wise) error correction
asmt = zeros(asit, 1);

% fill in basic stuff
vmpc.NrOfTimePoints = 0;
vmpc.Map(1).TimePointData = [];
vmpc.Map(1).RunTimeVars = struct( ...
    'MKDAMapType',    1, ...
    'PLPFile',        bcf, ...
    'PLPID',          bc.RunTimeVars.xffID, ...
    'ASimIterations', asit, ...
    'Condition',      ocond, ...
    'Contrast',       '', ...
    'ContrastColumn', contcol, ...
    'FWHM',           Inf, ...
    'FWHMMDist',      0, ...
    'Scaling',        opts.scale, ...
    'UnitColumn',     stcol, ...
    'UnitID',         0, ...
    'UnitPoints',     []);
vmpc.Map = vmpc.Map(1, ones(1, tmaps));

% initialize progress bar
simtot = asit * nsmk * numel(contnames);
tprg = 1 / simtot;
if isempty(pb)
    try
        pb = xprogress;
        xprogress(pb, 'setposition', [80, 200, 640, 36]);
        xprogress(pb, 'settitle', sprintf('Simulating MetaVMPs (%d units) ...', numstudies));
        xprogress(pb, 0, sprintf('Simulation %d/%d...', 1, simtot), 'visible', 0, 1);
    catch xfferror
        pb = [];
        neuroelf_lasterr(xfferror);
    end
end

% for each smoothing kernel
simnum = 0;
oldnow = now;
secu = 1 / 86400;
for smc = 1:nsmk

    % use smoothkern
    smk = smoothkern(opts.smkres .* opts.smkern([smc, smc, smc]));
    smks = size(smk, 1) - 1;
    smkh = 1 + 0.5 * smks;

    % what normalizing value
    smknv = flexinterpn_method(smk, [smkh - opts.smkmdist(smc), smkh, smkh], 'lanczos8');
    if smknv == 0
        smknv = min(smk(smk > 0));
    end

    % normalize, and limit to 1
    smk = limitrangec(smk ./ smknv, 0, 1, 0);

    % compute sum over all elements (that is the blob weight for size K!)
    sms = sum(smk(:)) ./ ((opts.res / opts.smkres) .^ 3);

    % now, we iterate over contrasts
    for ctc = 1:numel(contnames)

        % get the contrast
        cdef = opts.contrasts{ctc};

        % prepare shuffling
        if opts.contnull == 'u'

            % for unit shuffling, we require a random sample across units
            ullidx = find(ullist(:, ctc) ~= 0);
            ucons = permute(bssample(numel(ullidx), ...
                struct('perm', true, 'numsmp', asit)), [1, 3, 2]);
            uconl = repmat(ullist, 1, asit);
            uconl(ullidx, :) = ullist(ullidx(ucons));
        elseif opts.contnull ~= 's'
            pcons = permute(bssample(numel(pcon), ...
                struct('perm', true, 'numsmp', asit)), [1, 3, 2]);
            pconi = 1:contcolm;
        end
        cpos = cdef(cdef > 0);
        cneg = -cdef(cdef < 0);

        % and within, we start by iterating over the studies
        stsel = true(numstudies, 1);
        stidx = zeros(numstudies, 1);
        stpts = zeros(numstudies, 1);
        stpos = false(numstudies, 1);
        stneg = false(numstudies, 1);

        % keep simulation map
        if opts.asimkeep > 0
            simmap = zeros(vmsz);
        end
        asmn(:) = 0;

        % simulation iteration
        first = true;
        for asc = [asml, 1:asml]

            % for first N iterations
            if asc < asml

                % progress
                if ~isempty(pb) && (oldnow + secu) <= now
                    oldnow = now;
                    simnum = (smc - 1) * numel(contnames) * asit + (ctc - 1) * asit + asc;
                    pb.Progress(tprg * simnum, sprintf('Simulation %d/%d', simnum, simtot));
                end

                % spatial resampling
                if opts.contnull == 's'

                    % fully resample points from available space
                    if opts.asimsmpl == 'f'

                        % without restrictions
                        pcrd = maskc(ceil(maskn .* rand(ncrd, 1)), :);

                    % sample with displacement
                    else

                        % compute displacement values
                        cdisp = smpdisp(smc) .* randn(ncrd, 3);
                        dcdisp = opts.res * sqrt(sum(cdisp .* cdisp, 2));
                        lcdisp = limitrangec(dcdisp, 0.5 * opts.smkern(smc), ...
                            2 * opts.smkern(smc), opts.smkern(smc));
                        cdisp = cdisp .* ((lcdisp ./ dcdisp) * ones(1, 3));

                        % check that for each coordinate a replacement can be found with
                        pcrd = opcrd + cdisp;
                        pcrd(:, 4) = 1;
                        trcrd = round(pcrd * atrf);
                        tpok = all(trcrd >= 0.5 & trcrd < aszt, 2);
                        tpnok = find(~tpok);
                        tpnok = [tpnok(:); lsqueeze(find(asmc(sub2ind(asz, ...
                            trcrd(tpok, 1), trcrd(tpok, 2), trcrd(tpok, 3))) == 0))];
                        while ~isempty(tpnok)
                            cdisp = smpdisp(smc) * randn(numel(tpnok), 3);
                            dcdisp = opts.res * sqrt(sum(cdisp .* cdisp, 2));
                            lcdisp = limitrangec(dcdisp, 0.5 * opts.smkern(smc), ...
                                2 * opts.smkern(smc), opts.smkern(smc));
                            cdisp = cdisp .* ((lcdisp ./ dcdisp) * ones(1, 3));
                            pcrd(tpnok, 1:3) = opcrd(tpnok, :) + cdisp;
                            trcrd(tpnok, :) = round(pcrd(tpnok, :) * atrf);
                            tpok = all(trcrd >= 0.5 & trcrd < aszt, 2);
                            tpnok = find(~tpok);
                            tpnok = [tpnok(:); lsqueeze(find(asmc(sub2ind(asz, ...
                                trcrd(tpok, 1), trcrd(tpok, 2), trcrd(tpok, 3))) == 0))];
                        end

                        % remove 4th column
                        pcrd(:, 4) = [];

                        % round
                        if opts.fixcoords
                            pcrd = round(pcrd);
                        end
                    end

                % label shuffling
                else

                    % get original coords back
                    p = opfull;

                    % global shuffling
                    if opts.contnull == 'f'

                        % replace all labels
                        pconi(pcon) = pcon(pcons(:, asc));
                        p(pconr, contcol) = pconi(p(pconr, contcol));
                    end
                end

            % at last iteration
            else

                % progress
                if ~isempty(pb)
                    pb.Progress(tprg * simnum, 'Meta VMP computation');
                end

                % take actual points
                if opts.contnull ~= 's'
                    p = opfull;
                end
                pcrd = opcrd;
            end

            % then iterate over studies
            for stc = 1:numstudies

                % skip study
                if ~stsel(stc)
                    continue;
                elseif opts.contnull == 'u' && asc < asml
                    break;
                end

                % map counter (index into VMP content)
                tmc = (smc - 1) * smaps * numel(opts.contrasts) + (ctc - 1) * smaps + stc;
                stidx(stc) = tmc;

                % points counter
                stpc = 0;
                stnc = 0;

                % name the map in the VMP
                if asc == asml
                    if nsmk > 1
                        vmpc.Map(tmc).Name = sprintf('smk%.1fmm: %s (%s)', ...
                            opts.smkern(smc), stlabels{stc}, contnames{ctc});
                    else
                        vmpc.Map(tmc).Name = sprintf('%s (%s)', stlabels{stc}, contnames{ctc});
                    end
                    vmpc.Map(tmc).RunTimeVars.Contrast = contnames{ctc};
                    vmpc.Map(tmc).RunTimeVars.FWHM = opts.smkern(smc);
                    vmpc.Map(tmc).RunTimeVars.FWHMMDist = opts.smkmdist(smc);
                    vmpc.Map(tmc).RunTimeVars.UnitID = opts.studysel(stc);
                    vmpc.Map(tmc).RunTimeVars.UnitPoints = p(p(:, stcol) == opts.studysel(stc), :);

                    % set limits
                    vmpc.Map(tmc).LowerThreshold = 0.5;
                    vmpc.Map(tmc).UpperThreshold = 1;
                end

                % then we get the map (in double precision)
                smapp(:) = 0;
                if ~isempty(cneg)
                    smapn(:) = 0;
                end

                % first we find all points for that study
                stp = find(p(:, stcol) == opts.studysel(stc));

                % shuffle labels within studies
                if ~any(opts.contnull == 'su') && asc < asml

                    % sort random variable to keep balancing (stringent)
                    [rndv, rndi] = sort(rand(numel(stp), 1));
                    p(stp, contcol) = p(stp(rndi), contcol);
                end

                % then we iterate over positive particles
                for cdc = 1:numel(cpos)

                    % and file points of that particle
                    cpp = stp(p(stp, contcol) == cpos(cdc));

                    % update points
                    if isempty(cpp)
                        continue;
                    end
                    stpos(stc) = true;

                    % unique points
                    if opts.unique
                        [up, upi] = unique(opcrd(cpp, :), 'rows');
                        cpp = cpp(upi);
                    end

                    % get sub-coordinate parts
                    icrd = floor(pcrd(cpp, :));
                    scrd = pcrd(cpp, :) - icrd;

                    % compute offset (for this coordinate)
                    soff = smkh - resrat .* scrd;
                    rsoff = round(soff);
                    arsoff = (abs(soff - rsoff) < 0.001);
                    soff(arsoff) = rsoff(arsoff);

                    % now compute the first voxels to be sampled
                    smin = min(ceil(soff ./ resrat), icrd - 1);
                    sming = soff - smin .* resrat;
                    smini = icrd - smin;

                    % compute how many voxels are to be sampled
                    smax = min(ceil((smks - soff) ./ resrat), vmsz(ones(numel(cpp), 1), :) - icrd);
                    smaxg = soff + (smax + 0.1) .* resrat;
                    smaxi = icrd + smax;

                    % then we iterate over those points
                    for pc = lsqueeze(find(all(sming <= smaxg, 2)))'

                        % set as start value
                        ipgrid(2, :) = sming(pc, :);

                        % and set as end value
                        ipgrid(4, :) = smaxg(pc, :);

                        % sample kernel piece
                        smkp = flexinterpn(smk, ipgrid, opts.smkinterp{:});

                        % scaling
                        if opts.scale == 'i'
                            smkp = double(smkp >= 0.5);
                            stpc = stpc + 1;
                        else
                            stpc = stpc + sum(smkp(:)) / sms;
                        end

                        % update map
                        if opts.jbmeth == 'r'
                            smapp(smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)) = min(1, smapp( ...
                                smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)) + smkp);
                        else
                            smapp(smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)) = max(smapp( ...
                                smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)), smkp);
                        end
                    end
                end

                % then we iterate over negative particles
                for cdc = 1:numel(cneg)

                    % and file points of that particle
                    cpp = stp(p(stp, contcol) == cneg(cdc));

                    % update points
                    if isempty(cpp)
                        continue;
                    end
                    stneg(stc) = true;

                    % unique points
                    if opts.unique
                        [up, upi] = unique(opcrd(cpp, :), 'rows');
                        cpp = cpp(upi);
                    end

                    % get sub-coordinate parts
                    icrd = floor(pcrd(cpp, :));
                    scrd = pcrd(cpp, :) - icrd;

                    % compute offset (for this coordinate)
                    soff = smkh - resrat .* scrd;
                    rsoff = round(soff);
                    arsoff = (abs(soff - rsoff) < 0.001);
                    soff(arsoff) = rsoff(arsoff);

                    % now compute the first voxels to be sampled
                    smin = min(ceil(soff ./ resrat), icrd - 1);
                    sming = soff - smin .* resrat;
                    smini = icrd - smin;

                    % compute how many voxels are to be sampled
                    smax = min(ceil((smks - soff) ./ resrat), vmsz(ones(numel(cpp), 1), :) - icrd);
                    smaxg = soff + (smax + 0.1) .* resrat;
                    smaxi = icrd + smax;

                    % then we iterate over those points
                    for pc = lsqueeze(find(all(sming <= smaxg, 2)))'

                        % set as start value
                        ipgrid(2, :) = sming(pc, :);

                        % and set as end value
                        ipgrid(4, :) = smaxg(pc, :);

                        % sample kernel piece
                        smkp = flexinterpn(smk, ipgrid, opts.smkinterp{:});

                        % scaling
                        if opts.scale == 'i'
                            smkp = double(smkp >= 0.5);
                            stpc = stpc + 1;
                        else
                            stpc = stpc + sum(smkp(:)) / sms;
                        end

                        % update map
                        if opts.jbmeth == 'r'
                            smapn(smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)) = min(1, smapn( ...
                                smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)) + smkp);
                        else
                            smapn(smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)) = max(smapn( ...
                                smini(pc, 1):smaxi(pc, 1), ...
                                smini(pc, 2):smaxi(pc, 2), ...
                                smini(pc, 3):smaxi(pc, 3)), smkp);
                        end
                    end
                end

                % total points (positive and negative blobs)
                if asc == asml
                    stpts(stc) = stpc + stnc;
                    if stpts(stc) < 0.5
                        stsel(stc) = false;
                        if ~any(stsel)
                            if isempty(opts.pbar)
                                closebar(pb);
                            end
                            error('neuroelf:xff:noPointsInSet', 'No more studies with points.');
                        end
                        continue;
                    end
                    if ~opts.stwsel && opts.contnull ~= 'u'
                        stsel(:) = true;
                    end
                end

                % if necessary
                if ~isempty(cneg)

                    % compute difference with correct exclusion criterion
                    if opts.contcomp == 'w'
                        vmpc.Map(tmc).VMPData = ...
                            (smapp < opts.contexclw | smapn < opts.contexclw) .* (smapp - smapn);
                    elseif opts.contcomp == 'e'
                        vmpc.Map(tmc).VMPData = (smapp == 0 | smapn == 0) .* (smapp - smapn);
                    else
                        vmpc.Map(tmc).VMPData = smapp - smapn;
                    end

                % otherwise
                else

                    % straight map
                    vmpc.Map(tmc).VMPData = smapp;
                end
            end

            % group map
            tmc = (smc - 1) * smaps * numel(opts.contrasts) + (ctc - 1) * smaps + numstudies + 1;
            if first
                if nsmk > 1
                    vmpc.Map(tmc).Name = sprintf('smk%.1fmm: %s', opts.smkern(smc), contnames{ctc});
                else
                    vmpc.Map(tmc).Name = contnames{ctc};
                end
                vmpc.Map(tmc).RunTimeVars.FWHM = opts.smkern(smc);
                vmpc.Map(tmc).RunTimeVars.FWHMMDist = opts.smkmdist(smc);
                vmpc.Map(tmc).RunTimeVars.MKDAMapType = 2;
                vmpc.Map(tmc).RunTimeVars.Contrast = contnames{ctc};
                vmpc.Map(tmc).RunTimeVars.UnitID = opts.studysel;

                % nr-of-points weights
                switch (opts.stwp)
                    case {'c'}
                        stw = flexinterpn(stwp, 1 + stpts(stsel));
                    case {'l'}
                        stw = log(1 + stpts(stsel));
                    case {'n'}
                        stw = ones(sum(stsel), 1);
                    case {'p'}
                        stw = stpts(stsel);
                    case {'s'}
                        stw = sqrt(stpts(stsel));
                end

                % multiply with preweights, factoring in "bad points"
                if opts.stwmcc
                    stw = stw .* stwf(stsel) .* (stpts(stsel) ./ stwu(stsel));
                else
                    stw = stw .* stwf(stsel);
                end

                % normalize max to 1
                if max(stw) > 0
                    stw = stw ./ max(stw);
                end

                % relative differential weights
                if opts.stwrdiff && ~isempty(cpos) && ~isempty(cneg) && opts.contnull == 's'

                    % compute sums first
                    numstpos = 1 / sum(stw(stpos));
                    numstneg = 1 / sum(stw(stneg));
                    numstdf1 = sum(stw);
                    numstall = 1 / numstdf1;
                    numstdf1 = max(2, floor(numstdf1));

                    % reweigh studies (but not with both particles!)
                    stw( stpos & ~stneg) = stw( stpos & ~stneg) .* numstpos;
                    stw(~stpos &  stneg) = stw(~stpos &  stneg) .* numstneg;
                    stw( stpos &  stneg) = stw( stpos &  stneg) .* numstall;
                end

                % keep original weighting
                stwo = stw;
                sstw = max(2, floor(sum(stw)) - 1);
            end

            % if not unit null or first run
            if opts.contnull ~= 'u' || first

                % combine maps to 4D array
                stm = cat(4, vmpc.Map(stidx(stsel)).VMPData);

                % for unit labeling and weighting
                if opts.contnull == 'u' && opts.grpmeth == 'w'

                    % reweight positive and negative to 1
                    uconlp = (ullist(stsel) > 0);
                    uconln = (ullist(stsel) < 0);
                    uconl0 = (ullist(stsel) == 0);
                    stw(uconlp) = stw(uconlp) ./ sum(stw(uconlp));
                    stw(uconln) = stw(uconln) ./ sum(stw(uconln));
                    stw(uconl0) = 0;
                end

            % for unit null
            elseif opts.contnull == 'u' && asc < asml

                % depending on grouping mechanism
                if opts.grpmeth == 'w'
                    stw = stwo .* (1 - 2 .* (uconl(stsel, asc) ~= ullist(stsel)));

                    % reweight positive and negative to 1
                    uconlp = (uconl(stsel, asc) > 0);
                    uconln = (uconl(stsel, asc) < 0);
                    uconl0 = (uconl(stsel, asc) == 0);
                    stw(uconlp) = stw(uconlp) ./ sum(stwo(uconlp));
                    stw(uconln) = stw(uconln) ./ sum(stwo(uconln));
                    stw(uconl0) = 0;
                else
                    ulflip = (ullast(stseli) ~= uconl(stseli, asc));
                    stm(:, :, :, stseli(ulflip)) = -stm(:, :, :, stseli(ulflip));
                    ullast = uconl(:, asc);
                end

            % otherwise
            else

                % restore stw
                stw = stwo;

                % for unit labeling and weighting
                if opts.contnull == 'u' && opts.grpmeth == 'w'

                    % reweight positive and negative to 1
                    uconlp = (ullist(stsel) > 0);
                    uconln = (ullist(stsel) < 0);
                    uconl0 = (ullist(stsel) == 0);
                    stw(uconlp) = stw(uconlp) ./ sum(stw(uconlp));
                    stw(uconln) = stw(uconln) ./ sum(stw(uconln));
                    stw(uconl0) = 0;
                end
            end

            % type of grouping
            switch (opts.grpmeth)

                % one-sample t-test
                case {'o'}
                    vmpc.Map(tmc).VMPData = onesamplet(stm, 4);
                    vmpc.Map(tmc).DF1 = numel(stw) - 1;

                % sum
                case {'s'}
                    vmpc.Map(tmc).VMPData = (1 / numel(stw)) .* single(sum(stm, 4));
                    vmpc.Map(tmc).DF1 = numel(stw) - 1;

                % weighted sum
                case {'w'}
                    sstm = reshape(reshape(stm, nmsz, numel(stw)) * stw(:), vmsz);
                    if opts.stwrdiff && ~isempty(cpos) && ~isempty(cneg) && any(opts.contnull == 'su')
                        vmpc.Map(tmc).VMPData = sstm;
                        if opts.contnull == 's'
                            vmpc.Map(tmc).DF1 = numstdf1;
                        else
                            vmpc.Map(tmc).DF1 = sstw;
                        end
                    else
                        vmpc.Map(tmc).VMPData = (1 / sum(stwo)) .* sstm;
                        vmpc.Map(tmc).DF1 = sstw;
                    end
            end

            % simulation
            if asc < asml

                % keep map of simulation
                if opts.asimkeep > 0
                    simmap = simmap + vmpc.Map(tmc).VMPData;
                end

                % for cluster-size based thresholding
                if opts.asimkthr

                    % get rough idea of distribution of values
                    aspdist = aspdist + histcount( vmpc.Map(tmc).VMPData(asmm), -1, 1, aspdfac);
                    aspcdist = cumsum(aspdist(:));

                    % then determine percentile thresholds
                    for pc = 1:numel(aspc)
                        aspc(pc) = aspdfac * (findfirst(aspcdist >= ...
                            (opts.asimrthr(pc) * aspcdist(end))) - 1) - 1;
                    end
                end

                % cluster size correction
                if opts.asimkthr
                    for pc = 1:numel(aspc)
                        ascs = max(clustercoordsc(vmpc.Map(tmc).VMPData >= aspc(pc)));
                        if ~isempty(ascs)
                            ascsizes(asc, pc) = ascs;
                        end
                    end
                end

                % keep track of maximum value (for basic FWE correction)
                asmt(asc) = max(vmpc.Map(tmc).VMPData(:));

                % add to percentile map
                asnn = vmpc.Map(tmc).VMPData(asmr);
                asmn = asmn + (asnn < aspe) + 0.5 .* (asnn == aspe);

            % for real
            elseif ~first

                % if cluster size correction is requested
                if opts.asimkthr

                    % store information in RunTimeVars field
                    aspccs = aspc;
                    for pc = 1:numel(aspc)
                        asscsizes = sort(ascsizes(:, pc));
                        aspccs(pc) = 1 + asscsizes(round(0.95 * asit));
                    end

                    vmpc.Map(tmc).RunTimeVars.AlphaSimMKDA = {0.05, ...
                        [1 - opts.asimrthr(:), aspc(:), aspccs(:)]};
                end

                % for FWE-based thresholding sort values and store
                vmpc.Map(tmc).RunTimeVars.FWEMaxDist = sort(asmt);

                % first copy map
                vmpc.Map(tmc + 1) = vmpc.Map(tmc);

                % then make settings
                vmpc.Map(tmc).Type = 9;
                vmpc.Map(tmc).LowerThreshold = ...
                    vmpc.Map(tmc).RunTimeVars.FWEMaxDist(min(numel(asmt), 1 + round(0.95 * asit)));
                vmpc.Map(tmc).UpperThreshold = max(vmpc.Map(tmc).VMPData(:));
                vmpc.Map(tmc).Name = [vmpc.Map(tmc).Name ' (proportion)'];

                % then compute rescaled z-score for each value
                asmn = asmn ./ asit;
                asmn(asmn == 0) = 0.5 / asml;
                asmn(asmn == 1) = (asit + 0.5) / asml;

                % make setting on copied map
                tmc = tmc + 1;
                vmpc.Map(tmc).Type = 1;
                vmpc.Map(tmc).Name = [vmpc.Map(tmc).Name ' (z-rescaled)'];
                vmpc.Map(tmc).VMPData(~asmr) = 0;
                vmpc.Map(tmc).VMPData(asmr) = sdist('norminv', asmn, 0, 1);
                vmpc.Map(tmc).RunTimeVars.MKDAMapType = 3;
                if isempty(cneg) || isempty(cpos)
                    if isempty(cneg)
                        vmpc.Map(tmc).ShowPositiveNegativeFlag = 1;
                    else
                        vmpc.Map(tmc).ShowPositiveNegativeFlag = 2;
                    end
                else
                    vmpc.Map(tmc).ShowPositiveNegativeFlag = 3;
                end
                vmpc.Map(tmc).LowerThreshold = -sdist('tinv', 0.05, vmpc.Map(tmc).DF1);
                vmpc.Map(tmc).UpperThreshold = -sdist('tinv', 1/asit, vmpc.Map(tmc).DF1);

                % keep simulation data
                if opts.asimkeep > 0
                    vmpc.Map(tmc + 1) = vmpc.Map(tmc);
                    vmpc.Map(tmc + 1).Type = 116;
                    vmpc.Map(tmc + 1).Name = [vmpc.Map(tmc + 1).Name ' (avg. null)'];
                    vmpc.Map(tmc + 1).VMPData = (1 / asit) .* simmap;
                    ssimmap = sort(simmap(simmap ~= 0));
                    vmpc.Map(tmc + 1).LowerThreshold = (1 / asit) * ...
                        ssimmap(min(numel(ssimmap), ceil(0.95 * numel(ssimmap))));
                    vmpc.Map(tmc + 1).UpperThreshold = 2 * vmpc.Map(tmc + 1).LowerThreshold;
                    vmpc.Map(tmc + 1).RunTimeVars.MKDAMapType = 4;
                end

            % otherwise (first pass for this combination)
            else

                % for cluster-size based thresholding
                if opts.asimkthr
                    aspdist(:) = 0;
                    ascsizes(:) = 0;
                end

                % keep point estimate (to compare simulations against)
                aspe = vmpc.Map(tmc).VMPData(asmr);

                % for unit-based null distribution
                if opts.contnull == 'u'

                    % keep study selection indices and original labels
                    stseli = find(stsel);
                    ullast = ullist;
                end

                % no longer first run
                first = false;
            end
        end
    end
end

% apply mask
if opts.applymask
    asmr = ~asmr;
    for tmc = 1:numel(vmpc.Map)
        vmpc.Map(tmc).VMPData(asmr) = 0;
    end
end

% make sure all maps are single precision (memory and disk IO)
for mc = 1:numel(vmpc.Map)
    vmpc.Map(mc).VMPData = single(vmpc.Map(mc).VMPData);
end

% and make sure RunTimeVars are saved
vmpc.RunTimeVars.AutoSave = true;

% clear progress bar object
if isempty(opts.pbar)
    closebar(pb);
end

% set back into output
mvmp.C = vmpc;
