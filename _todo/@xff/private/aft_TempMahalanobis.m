function md = aft_TempMahalanobis(xo, opts)
% AFT::TempMahalanobis  - temporal Mahalanobis Distance measure
%
% FORMAT:       md = obj.TempMahalanobis([opts]);
%
% Input fields:
%
%       opts        optional settings
%        .mthresh   implicit masking threshold, default: 1;
%                   values between 0 and 2 are relative; > 2 absolute
%        .select    selection type, one of 'hemifield', {'slice'}, 'voxel'
%                   for MTC files, this is fixed to 'voxel' (vertices)
%
% Output fields:
%
%       md          Mahalanobis Distance (across time)
%
% TYPES: FMR, HDR, HEAD, MTC, VTC

% Version:  v1.1
% Build:    16012511
% Date:     Jan-25 2016, 11:18 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, {'fmr', 'hdr', 'head', 'mtc', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
bc = xo.C;
ft = lower(xo.S.Extensions{1});
if nargin < 2 || ~isstruct(opts) || numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'mthresh') || ~isa(opts.mthresh, 'double') || numel(opts.mthresh) ~= 1 || ...
    isinf(opts.mthresh) || isnan(opts.mthresh) || opts.mthresh < 0
    opts.mthresh = 1;
end
if ~isfield(opts, 'select') || ~ischar(opts.select) || isempty(opts.select) || ...
   ~any(strcmpi(opts.select, {'h', 'hemi', 'hemifield', 's', 'slice', 'v', 'vertex', 'voxel'}))
    opts.select = 's';
else
    opts.select = lower(opts.select(1));
end

% switch on filetype
switch (ft)

    % FMR (STC)
    case 'fmr'

        % depending on FileVersion (content)
        if bc.FileVersion > 4 && numel(bc.Slice) == 1

            % get data with simple permute
            data = permute(bc.Slice.STCData(:, :, :, :), [3, 1, 2, 4]);

        % more complex storage (several files)
        else

            % get first slice data, then expand
            data = permute(bc.Slice(1).STCData(:, :, :), [3, 1, 2]);
            data(1, 1, 1, numel(bc.Slice)) = 0;

            % then get other slices
            for sc = 2:numel(bc.Slice)
                data(:, :, :, sc) = permute(bc.Slice(sc).STCData(:, :, :), [3, 1, 2]);
            end
        end

    % HDR (or NII)
    case 'hdr'

        % get data
        data = permute(bc.VoxelData(:, :, :, :), [4, 1, 2, 3]);

    % HEAD (BRIK)
    case {'head'}

        % get data
        bricks = bc.Brick;
        for c = 1:numel(bricks)
            if istransio(bricks(c).Data)
                bricks(c).Data = resolve(bricks(c).Data);
            end
        end
        data = permute(cat(4, bricks.Data), [4, 1, 2, 3]);

    % MTC
    case {'mtc'}

        % force to vertex selection
        opts.select = 'v';
        data = bc.MTCData(:, :);

    % VTC
    case {'vtc'}

        % get data
        data = bc.VTCData(:, :, :, :);
end

% implicit masking
