function xo2 = fsbf_SaveAsSRF(xo, srffilename, options)
% FSBF::SaveAsSRF  - convert a FreeSurfer surface into BV's SRF
%
% FORMAT:       srf = fsbf.SaveAsSRF(srffilename [, options]);
%
% Input fields:
%
%       srffilename SRF filename
%       options     if given must be 1x1 struct with optional fields
%        .MeshCenter   1x3 double, default: [128, 128, 128]
%        .ConvexRGBA   1x4 double, default: [0.333, 0.677, 0.980, 0.400]
%        .ConcaveRGBA  1x4 double, default: [0.100, 0.240, 0.333, 0.400]
%        .VertexColor  Vx4 double (see SRF file format)
%        .RefFile      char, referenced file
%
% Output fields:
%
%       srf         SRF object
%
% Using: createsrf.

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:00 AM EST
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'fsbf') || ...
   ~ischar(srffilename) || isempty(srffilename)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end

% get sizes
bc = xo.C;
vcrd = bc.VertexCoordinate;
numv = size(vcrd, 1);

% check options
if nargin < 3 || ~isstruct(options) || numel(options) ~= 1
    options = struct;
end
if isfield(options, 'MeshCenter') && isa(options.MeshCenter, 'double') && ...
    numel(options.MeshCenter) == 3 && ~any(isinf(options.MeshCenter) | isnan(options.MeshCenter))
    mshc = options.MeshCenter(:)';
else
    mshc = [128, 128, 128];
end
if isfield(options, 'ConvexRGBA') && isa(options.ConvexRGBA, 'double') && numel(options.ConvexRGBA) == 4 && ...
    ~any(isinf(options.ConvexRGBA(:)) | isnan(options.ConvexRGBA(:)) | ...
        options.ConvexRGBA(:) < 0 | options.ConvexRGBA(:) > 1)
    cvxr = options.ConvexRGBA(:)';
else
    cvxr = [0.333, 0.677, 0.98, 0.4];
end
if isfield(options, 'ConcaveRGBA') && isa(options.ConcaveRGBA, 'double') && numel(options.ConcaveRGBA) == 4 && ...
    ~any(isinf(options.ConcaveRGBA(:)) | isnan(options.ConcaveRGBA(:)) | ...
        options.ConcaveRGBA(:) < 0 | options.ConcaveRGBA(:) > 1)
    ccvr = options.ConcaveRGBA(:)';
else
    ccvr = [0.1, 0.24, 0.333, 0.4];
end
if isfield(options, 'VertexColor') && isa(options.VertexColor, 'double') && ...
    size(options.VertexColor, 1) == numv && size(options.VertexColor, 2) == 4 && ...
    numel(options.VertexColor) == (4 * numv)
    vcol = options.VertexColor;
else
    vcol = zeros(numv, 4);
end
if isfield(options, 'RefFile') && ischar(options.RefFile) && ~isempty(options.RefFile)
    rfl = options.RefFile(:)';
else
    rfl = '<none>';
end

% create SRF file in memory
xo2 = ne_methods.createsrf(vcrd, bc.TriangleVertex, struct('offset', mshc));

% set options
bc2 = xo2.C;
bc2.ConvexRGBA = cvxr;
bc2.ConcaveRGBA = ccvr;
bc2.VertexColor = vcol;
bc2.AutoLinkedSRF = rfl;
xo2.C = bc2;

% check file writability
try
    aft_SaveAs(xo2, srffilename);
catch xfferror
    warning('neuroelf:xff:fileNotWritable', ...
        'File not writable: ''%s'': %s.', srffilename, xfferror.message);
end
