function [vb, vbv] = glm_POIBetas(xo, pfile, opts)
% GLM::POIBetas  - returns a table of POI betas (per subjects)
%
% FORMAT:       [pb, pbv] = glm.POIBetas(poi [, opts]);
%
% Input fields:
%
%       voi         VOI object
%       opts        optional settings
%        .c         CxP contrasts (defaults: each beta map on its own)
%        .pl        indices of POIs in object to sample (default: all)
%        .rmean     remove mean (of map) first (default: false)
%        .robust    perform robust mean estimation of average (for POIs
%                   with at least 4 vertices, default: false)
%        .rfx       return RFX table if possible (default: true)
%
% Output fields:
%
%       pb          SxCxV double table with data
%       pbv         SxCxV cell array with source data
%
% Using: applyspmsnc, bvcoordconv, findfirst, fitrobustbisquare,
%        flexinterpn_method, lsqueeze, makelabel, meannoinfnan.

% Version:  v1.1
% Build:    16031220
% Date:     Mar-12 2016, 8:53 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2016, Jochen Weber
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

% neuroelf library
global ne_methods;
findfirst          = ne_methods.findfirst;
fitrobustbisquare  = ne_methods.fitrobustbisquare;

% check arguments
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'glm') || ...
   ((numel(pfile) ~= 1 || ~xffisobject(pfile, true, 'poi')) && ...
    (~isa(pfile, 'double') || size(pfile, 2) ~= 1 || isempty(pfile) || ...
     any(isinf(pfile(:)) | isnan(pfile(:)) | pfile(:) < 1 | pfile(:) > xo.C.NrOfVertices)))
    error('neuroelf:xff:badArgument', 'Invalid object handle in call.');
end
bc = xo.C;
if bc.ProjectType ~= 2
    error('neuroelf:xff:badArgument', 'Only valid for MTC based GLM files.');
end
nrvert = bc.NrOfVertices;
switch (bc.ProjectTypeRFX)
    case 0
        ns = 1;
        np = bc.NrOfPredictors;
    case 1
        ns = numel(bc.GLMData.Subject);
        np = bc.NrOfSubjectPredictors;
    otherwise
        error('neuroelf:xff:badArgument', 'Invalid/unsupported ProjectTypeRFX flag in file.');
end
if ~isa(pfile, 'double')
    vc = pfile.C;
else
    vc = struct('POI', struct('Vertices', pfile));
end
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'c') || ~isa(opts.c, 'double') || isempty(opts.c) || ...
   (size(opts.c, 2) ~= np && size(opts.c, 2) ~= (np - 1)) || ...
    any(isinf(opts.c(:)) | isnan(opts.c(:))) || any(any(opts.c < 0, 2) & sum(opts.c, 2) ~= 0)
    opts.c = eye(np);
end
nc = size(opts.c, 1);
if ~isfield(opts, 'pl') || ~isa(opts.pl, 'double') || isempty(opts.pl) || ...
    any(isinf(opts.pl(:)) | isnan(opts.pl(:)))
    opts.pl = 1:numel(vc.POI);
else
    opts.pl = intersect(opts.pl(:)', 1:numel(vc.POI));
end
if ~isfield(opts, 'rmean') || ~islogical(opts.rmean) || numel(opts.rmean) ~= 1
    opts.rmean = false;
end
if ~isfield(opts, 'rfx') || ~islogical(opts.rfx) || numel(opts.rfx) ~= 1
    opts.rfx = true;
end
if ~isfield(opts, 'robust') || ~islogical(opts.robust) || numel(opts.robust) ~= 1
    opts.robust = false;
end

% prepare output
nv = numel(opts.pl);
vb = zeros(ns, nc, nv);
if nargout > 1
    vbv = cell(size(vb));
end

% check vertices
voi = cell(1, nv);
for c = 1:nv
    voi{c} = vc.POI(c).Vertices(:);
    voi{c}(isinf(voi{c}) | isnan(voi{c}) | voi{c} < 1 | voi{c} > nrvert) = [];
    voi{c} = unique(round(voi{c}));
end

% for RFX GLMs
if bc.ProjectTypeRFX == 1

    % iterate over subjects
    for sc = 1:ns

        % iterate over VOIs first
        for vc = 1:nv

            % reject empty VOIs
            if isempty(voi{vc})
                continue;
            end

            % get sampled values per contrast
            for cc = 1:nc
                conval = zeros(numel(voi{vc}), 1);
                for pc = find(opts.c(cc, :) ~= 0)
                    xval = bc.GLMData.Subject(sc).BetaMaps(voi{vc}, pc);
                    conval = conval + opts.c(cc, pc) .* xval;
                end
                if nargout > 1
                    vbv{sc, cc, vc} = conval;
                end
                conval(conval == 0) = [];
                if ~opts.robust || numel(conval) < 4
                    vb(sc, cc, vc) = sum(conval) / numel(conval);
                else
                    vb(sc, cc, vc) = fitrobustbisquare(ones(numel(conval), 1), double(conval));
                end
            end
        end
    end

% for non-RFX GLMs
else

    % iterate over VOIs first
    for vc = 1:nv

        % reject empty VOIs
        if isempty(voi{vc})
            continue;
        end

        % get sampled values per contrast
        for cc = 1:nc
            conval = zeros(numel(voi{vc}), 1);
            for pc = find(opts.c(cc, :) ~= 0)
                xval = bc.GLMData.BetaMaps(voi{vc}, pc);
                conval = conval + opts.c(cc, pc) .* xval;
            end
            if nargout > 1
                vbv{1, cc, vc} = conval;
            end
            conval(conval == 0) = [];
            if ~opts.robust || numel(conval) < 4
                vb(1, cc, vc) = sum(conval) / numel(conval);
            else
                vb(1, cc, vc) = fitrobustbisquare(ones(numel(conval), 1), double(conval));
            end
        end
    end

    % try to re-shape into RFX betas
    if opts.rfx && bc.SeparatePredictors == 2

        % create copy of vb
        vbc = vb;
        if nargout > 1
            vbvc = vbv;
        end

        % get subjects and subject predictors and full predictor names
        subjects = glm_Subjects(xo);
        subpreds = glm_SubjectPredictors(xo);
        fpnames = {bc.Predictor.Name2};

        % generate new vb
        vb = NaN .* zeros(numel(subjects), numel(subpreds), size(vbc, 3));
        if nargout > 1
            vbv = cell(size(vb));
        end

        % look up indices
        for sc = 1:numel(subjects)
            for pc = 1:numel(subpreds)

                % copy data
                vbi = findfirst(~cellfun('isempty', regexpi(fpnames, ...
                    sprintf('^subject\\s+%s\\:\\s+%s$', subjects{sc}, subpreds{pc}))));
                if ~isempty(vbi)
                    vb(sc, pc, :) = vbc(1, vbi, :);
                    if nargout > 1
                        vbv(sc, pc, :) = vbvc(1, vbi, :);
                    end
                end
            end
        end
    end
end
