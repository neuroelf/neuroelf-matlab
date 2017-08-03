% PUBLIC FUNCTION ne_mkda_plotonsurf: plot points on surface display
function varargout = ne_mkda_plotonsurf(varargin)

% Version:  v1.1
% Build:    16032218
% Date:     Mar-22 2016, 6:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2012, 2016, Jochen Weber
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
ch = ne_gcfg.h.MKDA.h;
mh = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get content of PLP dropdown
plps = ch.PLPs;
plpud = plps.UserData;
plpid = plps.Value;
try
    plp = plpud{plpid, 3};
    if numel(plp) ~= 1 || ...
       ~isxff(plp, 'plp')
        error( ...
            'neuroelf:GUI:BadPLP', ...
            'Bad PLP object.' ...
        );
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    return;
end
rtv = plp.RunTimeVars;

% make sure *some surfaces* are loaded
if ne_gcfg.fcfg.page ~= 3
    ne_showpage(0, 0, 3);
    drawnow;
end
sval = 1;
if isempty(mh.SurfVar.UserData)
    try
        lh = [neuroelf_path('colin') , '/colin_midres_LH_SPH40k_ICBMnorm.srf'];
        rh = [neuroelf_path('colin') , '/colin_midres_RH_SPH40k_ICBMnorm.srf'];
        if exist(lh, 'file') ~= 2 || ...
            exist(rh, 'file') ~= 2
            lh = [neuroelf_path('colin') , '/colin_LH_SPH_ICBMnorm.srf'];
            rh = [neuroelf_path('colin') , '/colin_RH_SPH_ICBMnorm.srf'];
        end
        if exist(lh, 'file') > 0 && ...
            exist(rh, 'file') > 0
            lh = xff(lh);
            rh = xff(rh);
            ne_openfile(0, 0, lh, true);
            ne_openfile(0, 0, rh, true);
            sval = [1; 2];
            mh.Scenery.Value = [1; 2];
            ne_sceneryselect;
            drawnow;
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
else
    sval = mh.Scenery.Value(:);
end

% get condition particles
condstr = ch.CndParts.String;
if ~iscell(condstr)
    if isempty(condstr)
        condstr = {};
    else
        condstr = cellstr(condstr);
    end
end

% form string
if isempty(condstr)
    condstr = '';
else
    condstr = gluetostringc(condstr, ' ');
end

% extend by contrast particles
cparts = ch.Contrast.UserData;
if iscell(cparts{1})
    cparts = [cparts{1}(:); cparts{2}(:)];
else
    cparts = cparts(:);
end
ccol = ch.ContColumn.String{ch.ContColumn.Value};
ccoltext = ch.ContColumn.UserData;
if isempty(ccoltext)
    if isfield(rtv, 'ColumnIsText') && ...
        isfield(rtv.ColumnIsText, ccol)
        ccoltext = rtv.ColumnIsText.(ccol);
    end
end
for cc = 1:numel(cparts)
    if ccoltext || ...
       (isempty(regexpi(ddeblank(cparts{cc}), '^[+\-]?\d+(\.\d*)?(e[+\i]?\d+)?$')) && ...
        ~any(strcmpi(cparts{cc}, {'inf', '-inf', 'nan'})))
        cparts{cc} = sprintf(' $%s == ''%s'' ', ccol, cparts{cc});
    else
        cparts{cc} = sprintf(' $%s == %s ', ccol, cparts{cc});
    end
end
cparts = gluetostringc(cparts, '|');
if isempty(condstr)
    condstr = cparts;
else
    condstr = ['(' condstr ') & (' cparts ')'];
end

% extend by study selection
allstudies = plp.Study;
studies = unique(allstudies);
selstudies = studies(ch.Studies.Value);
if numel(selstudies) < numel(studies)
    if numel(selstudies) > (0.5 * numel(studies))
        nonselstudies = setdiff(studies(:), selstudies(:));
        spart = sprintf('$Study ~= %d & ', nonselstudies(:)');
    else
        spart = sprintf('$Study == %d | ', selstudies(:)');
    end
    if condstr(1) == '('
        condstr = [condstr ' & (' spart(1:end-3) ')'];
    else
        condstr = ['(' condstr ') & (' spart(1:end-3) ')'];
    end
end

% undock surface view
hPFig = ne_undock;
hPTag = hPFig.Tag;
hPTag = hPTag(1:8);

% create plotting config
plppcfg = struct( ...
    'axes',    ne_gcfg.cc.(hPTag).Surface, ...
    'bcolor',  rtv.Config.BorderColor, ...
    'color',   ne_gcfg.c.ini.PLPPlot.Color, ...
    'cond',    condstr, ...
    'hover',   ne_gcfg.c.ini.PLPPlot.Hover, ...
    'symsize', ne_gcfg.c.ini.PLPPlot.SymbolSize, ...
    'view',    'fr');
if ~isempty(rtv.Config.LabelColumn)
    plppcfg.label = rtv.Config.LabelColumn;
    plppcfg.labcolor = rtv.Config.LabelColor;
end

% plot points
plp.PlotOnSurface(plppcfg);
pause(0.5);
drawnow;

% re-set surfaces in scenery
mh.Scenery.Value = [1; 2];
ne_sceneryselect;
drawnow;

% output
if nargout > 0
    varargout{1} = hPFig;
end
