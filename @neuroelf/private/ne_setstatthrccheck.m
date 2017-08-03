% FUNCTION ne_setstatthrccheck: toggle .EnableClusterCheck for VMPs
function varargout = ne_setstatthrccheck(varargin)

% Version:  v0.9b
% Build:    10081109
% Date:     Aug-11 2010, 9:00 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get configuration and handles
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% and get current StatsVar and map number
stvar = cc.StatsVar;
stvix = cc.StatsVarIdx;
if ~isxff(stvar, true)
    return;
end
stvtyp = lower(stvar.Filetype);

% only adapt settings for VMPs with single map selected
if any(strcmp(stvtyp, {'cmp', 'hdr', 'head', 'vmp'})) && ...
    numel(stvix) == 1

    % get current settings (cluster check flag and cluster size)
    oldv = [stvar.Map(stvix).EnableClusterCheck, stvar.Map(stvix).ClusterSize];

    % make new settings
    stvar.Map(stvix).EnableClusterCheck = double(ch.Stats.UsekThr.Value);
    stvar.Map(stvix).ClusterSize = ...
        max(1, floor(str2double(get(ch.Stats.kThresh, 'String'))));

    % get new settings
    newv = [stvar.Map(stvix).EnableClusterCheck, stvar.Map(stvix).ClusterSize];

    % update controls
    set(ch.Stats.kThresh, 'String', sprintf('%d', newv(2)));

    % if settings changed
    if any(oldv ~= newv)

        % clear correct table
        switch (stvtyp)
            case {'cmp'}
                stvar.Map(stvix).CMPDataCT = [];
            case {'hdr'}
                stvar.VoxelDataCT{stvix} = [];
            case {'head'}
                stvar.Brick(stvix).DataCT = [];
            case {'vmp'}
                stvar.Map(stvix).VMPDataCT = [];
        end

        % re-cluster
        stvar.ClusterTable(stvix, []);
    end
end

% and update window
ne_setslicepos;
