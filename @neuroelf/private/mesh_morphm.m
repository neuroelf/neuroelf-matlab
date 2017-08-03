function [c, nnei, ii, jj, nm] = mesh_morphm(c, n, tri, opts, nnei, ii, jj, nm, inej)
% mesh_morphm -  morph the coordinates of a mesh (no MEX)
%
% FORMAT:       c = mesh_morphm(c, n, tri, opts)
%
% Input fields:
%
%       c           Cx3 coordinate list (double)
%       n           Cx2 neighbors list (from SRF object, 1-based !)
%       tri         Tx3 triangle list (1-based !)
%       opts        mandatory struct with settings
%        .force     morphing force (default: 0.07)
%        .niter     number of iterations (default: 1)
%                 - optionally provided settings
%        .areac     if 1x1 double := 1, keep area constant
%                   (from initial state, requires .tri to be set!)
%        .distc     if given and between 0 .. 1, perform distortion corr
%        .distw     if given and [1], perform smoothing with distance
%                   weighting (default: false)
%        .distwl    take log(dist + 1) as factor for force (false)
%        .distwrhm  remove half of minimum distance prior to computation
%        .distwsq   take dist^2 as factor for force (false)
%        .sphere    to-sphere force
%        .type      1xN char type, currently only 'smooth' supported
%
% Output fields:
%
%       c           morphed coordinates
%
% This is a Matlab (M-file) approximate implementation of mesh_morph.c

% Version:  v0.9d
% Build:    14061612
% Date:     Jun-16 2014, 12:55 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% check arguments
if nargin < 3 || ...
   ~isa(c, 'double') || ...
    size(c, 1) < 3 || ...
    ndims(c) > 2 || ...
    any(isinf(c(:)) | isnan(c(:))) || ...
   ~iscell(n) || ...
   ~isequal(size(n), [size(c, 1), 2]) || ...
   ~isa(tri, 'double') || ...
    size(tri, 2) ~= 3 || ...
    any(isinf(tri(:)) | isnan(tri(:)) | tri(:) < 1 | tri(:) > size(c, 1))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 4 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'areac') || ...
    numel(opts.areac) ~= 1 || ...
   (~islogical(opts.areac) && ...
    (~isa(opts.areac, 'double') || ...
     ~any(opts.areac == [0, 1])))
    opts.areac = false;
elseif isa(opts.areac, 'double')
    opts.areac = logical(opts.areac);
end
if ~isfield(opts, 'distw') || ...
    numel(opts.distw) ~= 1 || ...
   (~islogical(opts.distw) && ...
    (~isa(opts.distw, 'double') || ...
     ~any(opts.distw == [0, 1])))
    opts.distw = false;
elseif isa(opts.distw, 'double')
    opts.distw = logical(opts.distw);
end
if ~isfield(opts, 'distwl') || ...
    numel(opts.distwl) ~= 1 || ...
   (~islogical(opts.distwl) && ...
    (~isa(opts.distwl, 'double') || ...
     ~any(opts.distwl == [0, 1])))
    opts.distwl = false;
elseif isa(opts.distwl, 'double')
    opts.distwl = logical(opts.distwl);
end
if ~isfield(opts, 'distwrhm') || ...
    numel(opts.distwrhm) ~= 1 || ...
   (~islogical(opts.distwrhm) && ...
    (~isa(opts.distwrhm, 'double') || ...
     ~any(opts.distwrhm == [0, 1])))
    opts.distwrhm = false;
elseif isa(opts.distwrhm, 'double')
    opts.distwrhm = logical(opts.distwrhm);
end
if ~isfield(opts, 'distwsq') || ...
    numel(opts.distwsq) ~= 1 || ...
   (~islogical(opts.distwsq) && ...
    (~isa(opts.distwsq, 'double') || ...
     ~any(opts.distwsq == [0, 1])))
    opts.distwsq = false;
elseif isa(opts.distwsq, 'double')
    opts.distwsq = logical(opts.distwsq);
end
if ~isfield(opts, 'force') || ...
    numel(opts.force) ~= 1 || ...
   ~isa(opts.force, 'double') || ...
    isinf(opts.force) || ...
    isnan(opts.force) || ...
    opts.force < 0
    opts.force = 0.07;
else
    opts.force = min(1, max(0, opts.force));
end
if ~isfield(opts, 'niter') || ...
   ~isa(opts.niter, 'double') || ...
    numel(opts.niter) ~= 1 || ...
    isinf(opts.niter) || ...
    isnan(opts.niter) || ...
    opts.niter < 1
    opts.niter = 1;
else
    opts.niter = min(50000, round(opts.niter));
end

%                 - optionally provided settings
%        .distc     if given and between 0 .. 1, perform distortion corr
%        .distw     if given and [1], perform smoothing with distance
%                   weighting (default: false)
%        .sphere    to-sphere force
%        .type      1xN char type, currently only 'smooth' supported

% build sparse neighbors weighting matrix vectors
try
    nv = size(c, 1);
    if nargin < 8 || ...
       ~isa(nnei, 'double') || ...
        size(nnei, 1) ~= numel(nnei) || ...
        numel(nnei) ~= nv || ...
       ~isa(ii, 'double') || ...
       ~isa(jj, 'double') || ...
       ~isequal(size(ii), size(jj)) || ...
        size(ii, 1) ~= numel(ii) || ...
       ~issparse(nm) || ...
       ~isequal(size(nm), [nv, nv])
        [ii, jj, nnei] = mesh_neighborsarray(n, 1);
        nm = sparse(ii, jj, 1, nv, nv, numel(ii));
        [ii, jj] = find(nm);
    end
    nii = numel(ii);
    if nargin < 9 || ...
       ~islogical(inej) || ...
        numel(inej) ~= nii
        inej = (ii ~= jj);
    end
    ieqj = find(~inej);
catch ne_eo;
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid neighbors cell array.' ...
    );
end

% and keep track of mesh center
onv = ones(nv, 1);
mc = ((1 / nv) .* onv) * sum(c, 1);
c = c - mc;

% additionally, need edge index vectors to compute other things
if opts.areac
    ntri = size(tri, 1);
    ev = [tri(:, [1, 2]); tri(:, [2, 3]); tri(:, [3, 1])];

    % compute total areas
    tel = reshape(mesh_dist(c, ev), [ntri, 3]);
    a = tri_area(tel);
    ta = sum(a);
end

% create vector for all weights
w = zeros(nii, 1);

% regular force to be applied
if opts.force > 0

    % force for self
    selfforce = 1 - opts.force;

    % no other options
    if ~opts.distw

        % compute force per neighbor evenly
        fnei = opts.force ./ nnei;

        % then set
        w(inej) = fnei(ii(inej));

        % fill force for self last
        w(ieqj) = selfforce;

        % and set in sparse matrix
        setsparseval(w, nm);

        % then apply
        for nc = 1:opts.niter
            c = nm * c;
        end

    % distance weighting
    elseif opts.distw

        % loop around here
        for nc = 1:opts.niter

            % compute reweighting factor (sum)
            md = mesh_dist(c, ii, jj);

            % remove half of minimum
            if opts.distwrhm
                md = max(0, md - 0.5 * min(md(md > 0)));
            end

            % another option?
            if opts.distwl
                md = log(1 + md);
            elseif opts.distwsq
                md = md .* md;
            end
            rw = condsum(nv, ii, md);
            md(ieqj) = (selfforce / opts.force) .* rw;

            % set as default weight
            rw = opts.force ./ rw;
            w = md .* indexarray(rw, ii);

            % then compute final weights and set in sparse matrix
            setsparseval(w, nm);

            % then multiply
            c = nm * c;
        end
    end
end

% keep area constant
if opts.areac

    % compute new area
    nta = sum(tri_area(reshape(mesh_dist(c, ev), [ntri, 3])));

    % and divide coordinates
    c = sqrt(ta / nta) .* c;
end

% add back center
c = c + mc;

% more outputs
if nargout > 1
    [ii, jj] = find(nm);
end


% subfunctions


% compute area for triangles
function [a, d] = tri_area(d, tri)

% one argument, area from triangle lengths
if nargin == 1

    % sort lengths
    l3 = sort(d, 2);

    % use the stable formula of Heron
    % see http://en.wikipedia.org/wiki/Heron's_formula#Numerical_stability
    a = 0.25 .* sqrt( ...
        (l3(:, 3) + (l3(:, 2) + l3(:, 1))) .* ...
        (l3(:, 1) - (l3(:, 3) - l3(:, 2))) .* ...
        (l3(:, 1) + (l3(:, 3) - l3(:, 2))) .* ...
        (l3(:, 3) + (l3(:, 2) - l3(:, 1))));

    % ensure no error is produced
    if nargout > 1
        d = [];
    end

% two arguments
else

    % compute distance from coordinates (in l3 instead) first
    d = cat(2, ...
        sqrt(sum((d(tri(:, 1), :) - d(tri(:, 2), :)) .^ 2, 2)), ...
        sqrt(sum((d(tri(:, 2), :) - d(tri(:, 3), :)) .^ 2, 2)), ...
        sqrt(sum((d(tri(:, 3), :) - d(tri(:, 1), :)) .^ 2, 2)));
    a = tri_area(d);
end
