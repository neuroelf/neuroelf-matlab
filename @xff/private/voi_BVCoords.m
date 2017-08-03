function bvc = voi_BVCoords(xo, voi, bbox)
% VOI::BVCoords  - return BV system coords for a voi
%
% FORMAT:       bvc = voi.BVCoords(voi [, bbox])
%
% Input fields:
%
%       voi         number of name of VOI
%       bbox        optional VMR bounding box (default: 256^3 TAL VMR)
%
% Output fields:
%
%       bvc         Nx3 coords (e.g. to index into 256^3 VMR space)
%
% Using: bvcoordconv.

% Version:  v1.1
% Build:    16021016
% Date:     Feb-10 2016, 4:35 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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
   (~ischar(voi) && ~isa(voi, 'double')) || isempty(voi)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 3 || ~isstruct(bbox) || numel(bbox) ~= 1 || numel(fieldnames(bbox)) ~= 7 || ...
   ~all(strcmp(fieldnames(bbox), {'BBox'; 'DimXYZ'; 'FCube'; 'RadCnv'; 'ResXYZ'; 'QuatB2T'; 'QuatT2B'}))
    bbox = struct('BBox', [0, 0, 0; 255, 255, 255], 'DimXYZ', [256, 256, 256], ...
        'FCube', 256, 'RadCnv', 1, 'ResXYZ', [1, 1, 1], ...
        'QuatB2T', [0, 0, -1, 129; -1, 0, 0, 129; 0, -1, 0, 129; 0, 0, 0, 1], ...
        'QuatT2B', [0, -1, 0, 129; 0, 0, -1, 129; -1, 0, 0, 129; 0, 0, 0, 1]);
end
bc = xo.C;
numvois = numel(bc.VOI);
if ischar(voi)
    voiok = false;
    for vc = numvois
        if ~isempty(regexpi(bc.VOI(vc).Name, lower(voi(:)')))
            voiok = true;
            brea;
        end
    end
    if ~voiok
        error('neuroelf:xff:invalidName', 'Named VOI not found.');
    end
    voi = vc;
else
    voi = real(voi(:)');
    if any(isinf(voi) | isnan(voi) | voi < 1 | voi > numvois) || ...
        numel(unique(fix(voi))) ~= numel(voi)
        error('neuroelf:xff:badArgument', 'Invalid VOI selection.');
    end
    voi = fix(voi);
end

% get voi
voi = bc.VOI(voi);

% initialize bvc
bvc = voi.Voxels;

% transform?
if strcmpi(bc.ReferenceSpace, 'tal')
    bvc = ne_methods.bvcoordconv(bvc, 'tal2bvi', bbox);
end
