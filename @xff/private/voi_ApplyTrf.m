function xo = voi_ApplyTrf(xo, trf, exact, invflag)
% VOI::ApplyTrf  - apply transformations to VOI coordinates
%
% FORMAT:       [voi = ] voi.ApplyTrf(trf [,exact, invflag])
%
% Input fields:
%
%       trf         single TRF or (cell array) list of TRF (or TAL)
%       exact       flag to not round coordinates (default false)
%       invflag     either single boolean or list of boolean
%
% Output fields:
%
%       voi         VOI with transformed coordinates
%
% Using: acpc2tal, applybvtrf.

% Version:  v1.1
% Build:    16021016
% Date:     Feb-10 2016, 4:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010 - 2012, 2014, 2016, Jochen Weber
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi') || ...
   (~xffisobject(trf, true, 'trf') && ~xffisobject(trf, true, 'tal') && ...
    (~iscell(trf) || isempty(trf) || ...
     (~xffisobject(trf{1}, true, 'trf') && ~xffisobject(trf{1}, true, 'tal'))))
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if ~iscell(trf)
    trf = {trf};
else
    trf = trf(:)';
end
for tc = 2:numel(trf)
    if numel(trf{tc}) ~= 1 || (~xffisobject(trf{tc}, true, 'trf') && ~xffisobject(trf{tc}, true, 'tal'))
        error('neuroelf:xff:badArgument', 'Invalid TRF argument in cell %d.', tc);
    end
end
if nargin < 3 || ~islogical(exact) || isempty(exact)
    exact = false;
else
    exact = exact(1);
end
if nargin < 4 || ~islogical(invflag) || (numel(invflag) ~= 1 && numel(invflag) ~= numel(trf))
    invflag = false(1, numel(trf));
elseif numel(invflag) == 1
    invflag = repmat(invflag, [1, numel(trf)]);
else
    invflag = invflag(:)';
end

% get VOI array and all coordinates
voi = bc.VOI;
nvoi = numel(voi);

% for non-exact transformation, mark all "inbetween voxels" as well
if ~exact
    for vc = 1:nvoi
        vv = voi(vc).Voxels;
        nv = size(vv, 1);
        vv = [vv; intersect(vv + repmat([0.5, 0, 0], [nv, 1]), ...
            vv - repmat([0.5, 0, 0], [nv, 1]), 'rows')];
        nv = size(vv, 1);
        vv = [vv; intersect(vv + repmat([0, 0.5, 0], [nv, 1]), ...
            vv - repmat([0, 0.5, 0], [nv, 1]), 'rows')];
        nv = size(vv, 1);
        vv = [vv; intersect(vv + repmat([0, 0, 0.5], [nv, 1]), ...
            vv - repmat([0, 0, 0.5], [nv, 1]), 'rows')];
        voi(vc).Voxels = vv;
    end
end

nvox = zeros(1, nvoi);
for vc = 1:nvoi
    nvox(vc) = size(voi(vc).Voxels, 1);
end
vox = zeros(sum(nvox), 3);
tvc = 0;
for vc = 1:nvoi
    vox((tvc + 1):(tvc + nvox(vc)), :) = voi(vc).Voxels;
    tvc = tvc + nvox(vc);
end

% transform?
if strcmpi(bc.ReferenceSpace, 'tal')
    cintal = true;
    vox = 128 - vox;
else
    cintal = false;
end

% perform transformations
for tc = 1:numel(trf)

    % TRF
    if xffisobject(trf{tc}, true, 'trf')

        % check trf
        tbc = trf{tc}.C;
        if tbc.TransformationType ~= 2 || ~strcmpi(tbc.DataFormat, 'matrix') || ...
            numel(tbc.TFMatrix) ~= 16
            error('neuroelf:xff:badArgument', 'Invalid TRF object given in cell %d.', tc);
        end

        % apply
        vox = ne_methods.applybvtrf(vox, tbc.TFMatrix, ~invflag(tc));

    % TAL
    else
        vox = ne_methods.acpc2tal(vox, trf{tc}, invflag(tc));
    end
end

% transform back
if cintal
    vox = 128 - vox;
end

% round
if ~exact
    vox = round(vox);
end

% spread into VOI array
tvc = 0;
for vc = 1:nvoi
    if exact
        voi(vc).Voxels = vox((tvc + 1):(tvc + nvox(vc)), :);
    else
        voi(vc).Voxels = unique(vox((tvc + 1):(tvc + nvox(vc)), :), 'rows');
    end
    voi(vc).NrOfVoxels = size(voi(vc).Voxels, 1);
    tvc = tvc + nvox(vc);
end

% set back
bc.VOI = voi;
xo.C = bc;
