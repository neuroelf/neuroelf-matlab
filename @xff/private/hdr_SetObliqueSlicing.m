function xo = hdr_SetObliqueSlicing(xo, osf)
% HDR::SetObliqueSlicing  - set slice obliqueness to cardinal axes
%
% FORMAT:       hdr.SetObliqueSlicing([flag]);
%
% Input fields:
%
%       flag        boolean flag, default: false
%
% No output fields.

% Version:  v1.1
% Build:    16020518
% Date:     Feb-05 2016, 6:27 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'hdr')
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 2 || ~islogical(osf) || numel(osf) ~= 1
    osf = true;
end

% get content
bc = xo.C;

% make sure RunTimeVars has DataHistBackup
if ~isfield(bc.RunTimeVars, 'DataHistBackup')
    bc.RunTimeVars.DataHistBackup = bc.DataHist;
end

% non-oblique (parallel to cardinal axes)
if ~osf

    % get size and coordinate frame
    sz = size(bc.VoxelData);
    if numel(sz) > 3
        sz(4:end) = [];
    end
    cfr = hdr_CoordinateFrame(xo);

    % factor
    f = min(256 ./ (sz .* cfr.Resolution));
    if f < 1
        f = 0.5;
    else
        f = floor(f);
    end

    % new resolution
    newres = f .* cfr.Resolution;

    % new origin
    neworg = round(-0.5 .* (newres .* sz));

    % set DataHist fields
    bc.DataHist.NIftI1.SFormCode = 2;
    bc.DataHist.NIftI1.AffineTransX = [newres(1), 0, 0, neworg(1)];
    bc.DataHist.NIftI1.AffineTransY = [0, newres(2), 0, neworg(2)];
    bc.DataHist.NIftI1.AffineTransZ = [0, 0, newres(3), neworg(3)];

    % keep track
    bc.RunTimeVars.Obliqued = false;

% or not
else

    % set back
    bc.DataHist = bc.RunTimeVars.DataHistBackup;
    bc.RunTimeVars.Obliqued = true;
end

% make sure any Trf is correct !
xo.C = bc;
cfr = hdr_CoordinateFrame(xo);
bc.RunTimeVars.Trf = inv(cfr.Trf)';
xo.C = bc;
