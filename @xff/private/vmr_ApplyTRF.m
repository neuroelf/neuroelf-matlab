function xo2 = vmr_ApplyTRF(xo, trf, opts)
% VMR::ApplyTRF  - apply transformation to VMR
%
% FORMAT:       newvmr = vmr.ApplyTRF(trf [, opts])
%
% Input fields:
%
%       trf         TRF object
%       opts        struct with optional fields
%        .asdouble  store output as double (default: false)
%        .inverse   apply inverse transformation (default: false)
%        .maxval    maximum value allowed in VMR (default: 225)
%        .method    interpolation, see flexinterpn_method (default: 'linear')
%        .v16       also apply to VMRData16 (if present, default: false)
%
% Output fields:
%
%       newvmr      transformed VMR
%
% Using: flexinterpn_method, limitrangec.

% Version:  v1.1
% Build:    16021316
% Date:     Feb-13 2016, 4:17 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2013, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'vmr') || ...
    numel(trf) ~= 1 || ~xffisobject(trf, true, 'trf')
    error('neuroelf:xff:badArgument', 'Bad or missing argument.');
end
bc = xo.C;
trfc = trf.C;
if nargin < 3 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'asdouble') || ~islogical(opts.asdouble) || numel(opts.asdouble) ~= 1
    opts.asdouble = false;
end
if ~isfield(opts, 'inverse') || ~islogical(opts.inverse) || numel(opts.inverse) ~= 1
    opts.inverse = false;
end
if ~isfield(opts, 'maxval') || ~isa(opts.maxval, 'double') || numel(opts.maxval) ~= 1 || ...
    isinf(opts.maxval) || isnan(opts.maxval) || opts.maxval <= 0
    opts.maxval = 225;
else
    opts.maxval = round(min(255, opts.maxval));
end
if ~isfield(opts, 'method') || ~ischar(opts.method) || isempty(opts.method)
    opts.method = 'linear';
else
    opts.method = lower(opts.method(:)');
end
if ~isfield(opts, 'v16') || ~islogical(opts.v16) || numel(opts.v16) ~= 1
    opts.v16 = false;
end

% the coordinate space is the full voxel space of the VMR
cs = [Inf; 1; 1; 256] * ones(1, 3);

% the basic transformations are to and from the center coordinate
fc = eye(4);
tc = eye(4);
fc(1:3, end) = -128.5;
tc(1:3, end) = 128.5;

% multiply matrices, regular transformation
if trfc.TransformationType > 1
    trfm = tc * trfc.TFMatrix * fc;

% alignment transformation
elseif trfc.AlignmentStep == 2
    trfm = tc * inv(trfc.TFMatrix) * fc;

% flip axes transform
else
    fx = zeros(4, 4);
    fx([2, 8, 9, 16]) = 1;
    trfm = tc * inv(trfc.TFMatrix) * fc * fx;
end

% inverse transform
if opts.inverse
    trfm = inv(trfm);
end

% apply transformation
if opts.v16 && isequal(size(bc.VMRData16), size(bc.VMRData))
    newvmrd16 = uint16(round(ne_methods.flexinterpn_method( ...
        bc.VMRData16(:, :, :), cs, 0, trfm, opts.method)));
else
    newvmrd16 = uint16([]);
end
newvmrd = ne_methods.flexinterpn_method(bc.VMRData(:, :, :), cs, 0, trfm, opts.method);

% copy old vmr
xo2 = aft_CopyObject(xo);
bc2 = xo2.C;

% set new data
if opts.asdouble
    bc2.VMRData = newvmrd;

% we need to set to uint8
else

    % make sure the range isn't violated
    bc2.VMRData = uint8(round(ne_methods.limitrangec(newvmrd, 0, opts.maxval, 0)));
end
bc2.VMRData16 = newvmrd16;

% add record of transformation transformed
bc2.Trf(end+1) = struct('NameOfSpatialTransformation', ...
    sprintf('NoName, %s interpolation', opts.method), ...
    'TypeOfSpatialTransformation', trfc.TransformationType, ...
    'SourceFileOfSpatialTransformation', trf.F, ...
    'NrOfSpatialTransformationValues', 16, 'TransformationValues', trfc.TFMatrix);

% set to storage
xo2.C = bc2;
