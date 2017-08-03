function v = voi_Voxels(xo, vsel, peaks)
% VOI::Voxels  - returns the voxels of all/the selected VOIs
%
% FORMAT:       voxels = voi.Voxels([vsel [, peaks]]);
%
% Input fields:
%
%       vsel        1xV VOI selection (default: all, if char regexpi match)
%       peaks       1x1 logical, peaks only (default: false)
%
% Output fields:
%
%       vnames      Nx1 list with VOI names

% Version:  v1.1
% Build:    16041811
% Date:     Apr-18 2016, 11:03 AM EST
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

% argument check
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'voi')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;

% double VOI selection
if nargin > 1 && isa(vsel, 'double') && ~isempty(vsel) && ~any(isinf(vsel(:)) | isnan(vsel(:)) | vsel(:) < 1)
    vsel = unique(min(numel(bc.VOI), round(vsel(:))));
elseif nargin > 1 && ischar(vsel) && ~isempty(vsel)
    vnames = {bc.VOI.Name};
    vsel = find(~cellfun('isempty', regexpi(vnames, vsel(:)')));
    if isempty(vsel)
        warning('neuroelf:xff:badArgument', 'VOI selection is empty.');
        v = zeros(0, 3);
        return;
    end
else
    vsel = 1:numel(bc.VOI);
end

% get all voxels
if nargin < 3 || ~islogical(peaks) || numel(peaks) ~= 1 || ~peaks
    v = unique(cat(1, bc.VOI(vsel).Voxels), 'rows');
else
    v = zeros(numel(vsel), 3);
    for vc = numel(vsel):-1:1
        if ~isempty(bc.VOI(vsel(vc)).Voxels)
            v(vc, :) = bc.VOI(vsel(vc)).Voxels(1, :);
        else
            v(vc, :) = [];
        end
    end
end
