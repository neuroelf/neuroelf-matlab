% FUNCTION ne_setstatthrtails: toggle pos/neg tail flags from controls
function varargout = ne_setstatthrtails(varargin)

% Version:  v0.9b
% Build:    13040517
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

% and currently selected stats object and map number
stvar = cc.StatsVar;
stvix = cc.StatsVarIdx;

% if this is a VMP and one map is selected
if isxff(stvar, {'cmp', 'glm', 'hdr', 'head', 'vmp', 'vtc'}) && ...
    numel(stvix) == 1

    % get current flag
    oldt = stvar.Map(stvix).ShowPositiveNegativeFlag;

    % then set new flag
    stvar.Map(stvix).ShowPositiveNegativeFlag = ...
        double(ch.Stats.PosTail.Value) + 2 * double(ch.Stats.NegTail.Value);

    % and compare value
    newt = stvar.Map(stvix).ShowPositiveNegativeFlag;

    % if not the same and clustering is requested
    if ~isequal(oldt, newt) && ...
        stvar.Map(stvix).EnableClusterCheck > 0 && ...
        isxff(stvar, 'vmp')

        % re-cluster with @VMP::ClusterTable
        stvar.Map(stvix).VMPDataCT = [];
        stvar.ClusterTable(stvix, []);
    end
end

% and update window
ne_setslicepos;
