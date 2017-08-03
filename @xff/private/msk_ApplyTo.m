function mxo = msk_ApplyTo(xo, mxo, cpflag)
% MSK::ApplyTo  - apply MSK to object
%
% FORMAT:       [vobj = ] msk.ApplyTo(vobj [, cpflag])
%
% Input fields:
%
%       vobj        xff object to be masked (VMP, VTC, ...)
%       cpflag      copy-flag (default: false)
%
% Output fields:
%
%       vobj        reference to masked object (same or copied)

% Version:  v1.1
% Build:    16020917
% Date:     Feb-09 2016, 5:06 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2011, 2014, 2016, Jochen Weber
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
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true, 'msk') || ...
    numel(mxo) ~= 1 || ~xffisobject(mxo, true, {'ava', 'cmp', 'glm', 'vdw', 'vmp', 'vtc'})
    error('neuroelf:xff:badArgument', 'Invalid call to ''%s''.', mfilename);
end
if nargin < 3 || ~islogical(cpflag) || numel(cpflag) ~= 1
    cpflag = false;
end

% compare bounding boxes
bbo = aft_BoundingBox(xo);
mbo = aft_BoundingBox(mxo);
if ~isequal(bbo.BBox, mbo.BBox) || ~isequal(bbo.DimXYZ, mbo.DimXYZ) || ...
   ~isequal(bbo.ResXYZ, mbo.ResXYZ) || ~isequal(bbo.FCube, mbo.FCube)
    error('neuroelf:xff:badArgument', 'Spatial layout mismatch between objects.');
end

% get mask
msk = (xo.C.Mask == 0);

% get data
bc = mxo.C;

% depending on filetype
switch lower(mxo.S.Extensions{1})

    % AVA
    case 'ava'

        % for all maps
        mf = fieldnames(bc.Maps);
        for mc = 1:numel(mf)

            % number of maps (as vector of 1's)
            rep = ones(1, size(bc.Maps.(mf{mc}), 4));

            % mask = 0
            bc.Maps.(mf{mc})(msk(:, :, :, rep)) = 0;
        end

    % CMP
    case 'cmp'

        % for all maps
        for mc = 1:numel(bc.Map)

            % mask = 0;
            bc.Map(mc).CMPData(msk) = 0;
        end

    % GLM
    case 'glm'

        % for all maps (but subject)
        mf = fieldnames(bc.GLMData);
        mf(strcmpi(mf, 'subject')) = [];
        for mc = 1:numel(mf)
            if isempty(bc.GLMData.(mf{mc}))
                continue;
            end

            % number of maps (as vector of 1's)
            rep = ones(1, size(bc.GLMData.(mf{mc}), 4));

            % mask = 0
            bc.GLMData.(mf{mc})(msk(:, :, :, rep)) = 0;
        end

        % treat subject
        for mc = 1:numel(bc.GLMData.Subject)

            % number of betamaps (as vector of 1's)
            rep = ones(1, size(bc.GLMData.Subject(mc).BetaMaps, 4));

            % mask = 0
            bc.GLMData.Subject(mc).BetaMaps(msk(:, :, :, rep)) = 0;
        end

    % VDW
    case 'vdw'

        % reshape mask
        msk = reshape(msk, [1, size(msk)]);

        % along time course
        bc.VDWData(msk(ones(1, size(bc.VDWData, 1)), :, :, :)) = 0;

    % VMP
    case 'vmp'

        % for all maps
        for mc = 1:numel(bc.Map)

            % mask = 0;
            bc.Map(mc).VMPData(msk) = 0;
        end

    % VTC
    case 'vtc'

        % reshape mask
        msk = reshape(msk, [1, size(msk)]);

        % along time course
        bc.VTCData(msk(ones(1, size(bc.VTCData, 1)), :, :, :)) = 0;
end

% copy content into new object?
if cpflag
    mxo = aft_CopyObject(mxo);
end

% set content
mxo.C = bc;
