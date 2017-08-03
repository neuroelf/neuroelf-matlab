function vmp = vtc_LocalCoherence(xo, opts)
% VTC::LocalCoherence  - compute local coherence map
%
% FORMAT:       vmp = vtc.LocalCoherence([opts])
%
% Input fields:
%
%       opts        options settings
%        .mask      'auto' (default: use whole brain mask) or mask object
%        .measure   one of 'meancorr', {'medcorr'}
%        .radius    search radius (used to construct search light mask, 5.5)
%        .remmedian remove median (within mask), default: true
% 
% Output fields:
%
%       vmp         VMP with a measure of local coherence
%
% Using: calcbetas, psctrans, tempfilter, ztrans.

% Version:  v1.1
% Build:    16081912
% Date:     Aug-19 2016, 12:53 PM EST
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

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vtc')
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
vtcd = double(bc.VTCData);
vtcs = size(vtcd);
vtcbb = aft_BoundingBox(xo);

% options
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end

% mask
if ~isfield(opts, 'mask') || ((numel(opts.mask) ~= 1 || ~xffisobject(opts.mask, true, 'msk')) && ...
   (~ischar(opts.mask) || isempty(opts.mask)))
    opts.mask = 'auto';
end

% load mask
clearmsk = false;
if ischar(opts.mask)

    % default mask
    if strcmpi(opts.mask, 'auto')
        opts.mask = [neuroelf_path('masks') filesep 'colin_brain_ICBMnorm_brain2mm.msk'];
    end
    try
        opts.mask = xff(opts.mask(:)');
    catch xfferror
        error('neuroelf:xff:maskNotReadable', 'Global signal mask not readable: %s', xfferror.message);
    end
    clearmsk = true;
end

% sample mask in VTC space
[msx, msy, msz] = ndgrid(1:vtcs(2), 1:vtcs(3), 1:vtcs(4));
msx = ne_methods.bvcoordconv([msx(:), msy(:), msz(:)], 'bvc2tal', vtcbb);
msz = (aft_SampleCoords(opts.mask, msx) >= 0.5);

% clear mask
if clearmsk
    aft_ClearObject(opts.mask);
end

% generate SL targets
slt = ne_methods.sltargets(reshape(msz, vtcs(2:4)));

