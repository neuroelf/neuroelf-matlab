function varargout = ne_setcsrfstatmap(varargin)
% ne_setcsrfstatmap  - set current SurfStatsVarIdx SMP maps
%
% FORMAT:       ne_setcsrfstatmap(SRC, EVT [, mapidx])
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       mapidx      optional indices to set for current stats object
%
% No output fields.
%
% Example:
%
%     ne_setcsrfstatmap(0, 0, 4);
%
%     sets the 4th map in the currently selected stats object as current
%     surface map.
%
% Note: this function issues btc_meshcolor() and ne_setcsrfstatbars() calls
%       after settings the maps.

% Version:  v1.1
% Build:    16052716
% Date:     May-27 2016, 4:01 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

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

% global variable
global ne_gcfg;
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% check again
ch.MainFig.SetGroupVisible('SMPMenu', 'off');
vtype = 'none';
if ~isxff(cc.SurfVar, {'fsbf', 'srf'}) || ~isxff(cc.SurfStatsVar, true) || ...
   ~any(strcmpi(cc.SurfStatsVar.FileType, {'fsmf', 'glm', 'mtc', 'smp'})) || ...
    cc.SurfVar.NrOfVertices ~= cc.SurfStatsVar.NrOfVertices
    ch.SurfStatsVarMaps.Value = [];
    ch.SurfStatsVarMaps.String = {'none (fitting the current surface)'};
    ch.SurfStatsVarMaps.Enable = 'off';
else
    vtype = lower(cc.SurfStatsVar.FileType);
end
stvar = cc.SurfStatsVar;

% index given
stvix = [];
if nargin > 2 && isa(varargin{3}, 'double') && ~isempty(varargin{3}) && ...
   ~any(isinf(varargin{3}(:)) | isnan(varargin{3}(:)) | varargin{3}(:) < 1)
    stvix = unique(min(numel(stvar.MapNames), fix(varargin{3}(:))));
end

% get index
if isempty(stvix)
    stvix = ch.SurfStatsVarMaps.Value;
end

% update indices
if isxff(stvar, true)
    stvarmnames = stvar.MapNames;
    stvar.RunTimeVars.MapSelection = {stvarmnames(stvix), stvix(:)};
    if vtype(1) == 'f'
        vtype = 'smp';
    end
    ch.MainFig.SetGroupVisible([upper(vtype) 'Menu'], 'on');
end
ne_gcfg.fcfg.SurfStatsVarIdx = stvix;
ch.SurfStatsVarMaps.Value = stvix;

% and in reference surface
cc.SurfVar.SetHandle('Stats', {stvar, stvix});

% call colorizer
btc_meshcolor(ne_gcfg.fcfg.SurfVar, true);

% stats bars
ne_setcsrfstatbars;

% and update title?
if numel(ch.Scenery.Value) == 1 && ...
    size(ch.Scenery.UserData, 1) >= ch.Scenery.Value && ...
    isxff(ch.Scenery.UserData{ch.Scenery.Value, 4}, true) && ...
    ch.Scenery.UserData{ch.Scenery.Value, 4} == cc.SurfVar
    if isxff(stvar, true)
        [smppath, ne_gcfg.c.title{2, 2}] = fileparts(stvar.FilenameOnDisk);
        if isempty(ne_gcfg.c.title{2, 2})
            ne_gcfg.c.title{2, 2} = 'unsaved';
        end
        if numel(stvix) == 1
            ne_gcfg.c.title{2, 3} = stvar.Map(stvix).Name;
        else
            ne_gcfg.c.title{2, 3} = sprintf('(%d maps)', numel(stvix));
        end
    else
        ne_gcfg.c.title(2, 2:3) = {'', ''};
    end
    ne_updatename;
end

% set group enabled?
if numel(stvix) == 1
    ch.MainFig.SetGroupEnabled('SngSMP', 'on');
    ch = ch.SrfStats;
    cm = stvar.Map(stvix);
    ch.LThresh.String = sprintf('%.4f', cm.LowerThreshold);
    ch.UThresh.String = sprintf('%.4f', cm.UpperThreshold);
    ch.PosTail.Value = double(mod(cm.ShowPositiveNegativeFlag, 2) == 1);
    ch.NegTail.Value = double(cm.ShowPositiveNegativeFlag > 1);
    ch.UsekThr.Value = double(cm.EnableClusterCheck > 0);
    set(ch.kThresh, 'String', sprintf('%.1fmm', cm.ClusterSize));
    if cm.UseRGBColor ~= 0
        ch.UseRGB.RadioGroupSetOne;
    else
        ch.UseLUT.RadioGroupSetOne;
    end
    bcolor = min(255, max(0, cm.RGBLowerThreshPos(:)));
    ch.RGBLPos.BackgroundColor = (1 / 255) .* bcolor';
    bcolor = min(255, max(0, cm.RGBUpperThreshPos(:)));
    ch.RGBUPos.BackgroundColor = (1 / 255) .* bcolor';
    bcolor = min(255, max(0, cm.RGBLowerThreshNeg(:)));
    ch.RGBLNeg.BackgroundColor = (1 / 255) .* bcolor';
    bcolor = min(255, max(0, cm.RGBUpperThreshNeg(:)));
    ch.RGBUNeg.BackgroundColor = (1 / 255) .* bcolor';

% otherwise
else

    % disable group
    ch.MainFig.SetGroupEnabled('SngSMP', 'off');
end
