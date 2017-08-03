function [xo, par] = fmr_Preprocess(xo, opts)
% FMR::Preprocess  - preprocess FMR
%
% FORMAT:       [fmr, par = ] fmr.Preprocess([opts])
%
% Input fields:
%
%       opts        optional settings
%        .steps     1xS cell array with letters/words
%                   'i'/'intensity' for global/slicewise intensity
%                   'st'/'slicetiming' for slice-timing correction
%                   'm'/'motion' for motion detection / correction
%                   's'/'smooth' for smoothing
%                   't'/'tempfilter' for temporal filtering
%                   if not given, use default: {'st', 'm'}
%        .icorrect  correct fluctuations, either of {'no'}, 'add', 'mult'
%        .islice    perform intensity fluctuations per slice (def: false)
%        .mdetect   boolean flag, only detect motion (parameters in par)
%        .msavemean boolean flag, save mean as one-vol FMR (def: false)
%        .msavepar  boolean flag, save parameters (default: true)
%        .msmooth   smoothing kernel in mm (default: twice the voxelsize)
%        .msmpl     sampling width in mm (either 1x1 or 1x3)
%        .mtomean   two-pass realignment (if value > 2 perform N passes)
%        .mtotarget target specification, see FMR::Realign method
%        .ssmooth   smoothing kernel in mm (default: [0, 0, 0])
%        .storder   slice order, default: from FMR
%        .tempdct   DCT-based filtering (min. wavelength, default: Inf)
%        .tempdt    detrend (default: true, is overriden by dct/sc)
%        .templp    temporal lowpass (smoothing) kernel in secs (def: 0)
%        .tempsc    sin/cos set of frequencies (number of pairs, def: 0)
%
% Output fields:
%
%       fmr         FMR with preprocessed data
%       par         output parameters from steps

% Version:  v1.1
% Build:    16020311
% Date:     Feb-03 2016, 11:47 AM EST
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
if numel(xo) ~= 1 || ~xffisobject(xo, true, 'fmr')
    error('neureolf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'savemean') || ~islogical(opts.savemean) || numel(opts.savemean) ~= 1
    opts.savemean = false;
end
if ~isfield(opts, 'savepar') || ~islogical(opts.savepar) || numel(opts.savepar) ~= 1
    opts.savepar = true;
end
if ~isfield(opts, 'smooth') || ~isa(opts.smooth, 'double') || ~any(numel(opts.smooth) == [1, 3]) || ...
    any(isinf(opts.smooth) | isnan(opts.smooth) | opts.smooth <= 1)
    opts.smooth = [];
elseif numel(opts.smooth) == 1
    opts.smooth = opts.smooth([1, 1, 1]);
end
if ~isfield(opts, 'smpl') || ~isa(opts.smpl, 'double') || ~any(numel(opts.smpl) == [1, 3]) || ...
    any(isinf(opts.smpl) | isnan(opts.smpl) | opts.smpl <= 1)
    opts.smpl = [];
elseif numel(opts.smpl) == 1
    opts.smpl = opts.smpl([1, 1, 1]);
end
if ~isfield(opts, 'tomean') || (~isa(opts.tomean, 'double') && ...
    ~islogical(opts.tomean)) || numel(opts.tomean) ~= 1
    opts.tomean = 0;
else
    if islogical(opts.tomean)
        opts.tomean = double(opts.tomean);
    end
    if isinf(opts.tomean) || isnan(opts.tomean) || opts.tomean < 0  || opts.tomean > 8
        opts.tomean = 0;
    end
    if opts.tomean < 2
        opts.tomean = opts.tomean + 1;
    end
    opts.tomean = round(opts.tomean);
end
