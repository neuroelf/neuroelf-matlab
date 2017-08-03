function mnames = aft_MapNames(xo, xflag)
% AFT::MapNames  - return a cell array (Nx1) of map names
%
% FORMAT:       mnames = obj.MapNames([extended]);
%
% Input fields:
%
%       extended    if given and 1x1 logical true, return extended names
%
% Output fields:
%
%       mnames      Nx1 list with map names
%
% TYPES: AVA, CMP, FSMF, GLM, HDR, HEAD, MTC, SMP, VMP, VTC

% Version:  v1.1
% Build:    16050717
% Date:     May-07 2016, 5:31 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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
if numel(xo) ~= 1 || ~xffisobject(xo, true)
    error('neuroelf:xff:badArgument', 'Invalid call to %s.', mfilename);
end
if nargin < 2 || numel(xflag) ~= 1 || ~islogical(xflag)
    xflag = false;
end
bc = xo.C;
ft = lower(xo.S.Extensions{1});
rtv = bc.RunTimeVars;
if isfield(bc, 'Map')
    map = bc.Map;
elseif isfield(rtv, 'Map')
    map = rtv.Map;
else
    map = [];
end

% switch over type
switch (ft)

    % AVA files
    case 'ava'
        am = fieldnames(bc.Maps);
        nm = zeros(1, numel(am));
        for c = 1:numel(am)
            nm(c) = size(bc.Maps.(am{c}), 4);
        end
        mnames = cell(1, sum(nm));
        mc = 1;
        for c = 1:numel(am)
            if nm(c) ~= 1
                for cc = 1:nm(c)
                    mnames{mc} = sprintf('%s %d', am{c}, cc);
                    mc = mc + 1;
                end
            else
                mnames{mc} = am{c};
                mc = mc + 1;
            end
        end

    % GLM files
    case 'glm'

        % get map names
        mnames = {bc.Predictor(:).Name2};

        % for RFX GLMs, reorder
        if bc.ProjectTypeRFX > 0
            mnameo = reshape(1:(numel(mnames) - numel(bc.GLMData.Subject)), ...
                size(bc.GLMData.Subject(1).BetaMaps, ...
                ndims(bc.GLMData.Subject(1).BetaMaps)) - 1, numel(bc.GLMData.Subject));
            mnameo(end + 1, :) = (numel(mnameo) + 1):numel(mnames);
            mnames = mnames(mnameo(:));
        end

    % HDR files
    case 'hdr'
        mnames = {bc.RunTimeVars.Map(:).Name};

    % AFNI HEAD/BRIK files
    case 'head'
        mnames = {bc.Brick(:).Label};

    % CMP/SMP/VMP files
    case {'cmp', 'fsmf', 'smp', 'vmp'}
        mnames = {bc.Map(:).Name};

    % MTC
    case 'mtc'
        if ~isfield(bc.RunTimeVars, 'AvgMTC') || ~bc.RunTimeVars.AvgMTC
            mnames = {'Time'};
        else
            mnames = bc.RunTimeVars.ConditionNames;
        end

    % VTC
    case 'vtc'
        if ~isfield(bc.RunTimeVars, 'AvgVTC') || ~bc.RunTimeVars.AvgVTC
            mnames = {'Time'};
        else
            mnames = {bc.RunTimeVars.Map.Name};
        end

    % default:
    otherwise

        % give back filename without extension
        [fn{1:3}] = fileparts(xo.F);

        % but if empty, create name in <>
        if isempty(fn{2})
            fn{2} = sprintf('<interactive %s>', ft);
        end
        mnames = fn(2);

end

% linearize
mnames = mnames(:);

% extended info
if xflag && numel(map) == numel(mnames)

    % add info
    for c = 1:numel(map)
        switch (map(c).Type)

            % t-map
            case {1}
                mnames{c} = sprintf('t[df=%d]: %s', map(c).DF1, mnames{c});

            % r-map
            case {2}
                mnames{c} = sprintf('r[df=%d]: %s', map(c).DF1, mnames{c});

            % CC-map
            case {3}
                mnames{c} = sprintf('CC-r[df=%d, %d lags]: %s', ...
                    map(c).DF1, map(c).NrOfLags, mnames{c});

            % F-map
            case {4}
                mnames{c} = sprintf('F[df=%d,%d]: %s', map(c).DF1, map(c).DF2, mnames{c});
        end
    end
end
