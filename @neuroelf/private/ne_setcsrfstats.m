% PUBLIC FUNCTION ne_setcsrfstats: set current SurfStatsVar SMP
function varargout = ne_setcsrfstats(varargin)

% Version:  v1.1
% Build:    16032218
% Date:     Mar-22 2016, 6:23 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2015, 2016, Jochen Weber
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

% get handles
ch = ne_gcfg.h;
cf = ch.MainFig;

% set menu group and also the group for loaded SRF stats to off
cf.SetGroupVisible('GLMMenu', 'off');
cf.SetGroupVisible('SMPMenu', 'off');
cf.SetGroupEnabled('SLoadSR', 'off');

% get selection and UserData (list of objects) of StatsVar
hu = ch.SurfStatsVar.UserData;
if nargin < 3 || ~isa(varargin{3}, 'double') || ...
    numel(varargin{3}) ~= 1 || isinf(varargin{3}) || isnan(varargin{3}) || ...
    varargin{3} < 1 || varargin{3} > size(hu, 1) || varargin{3} ~= fix(varargin{3})
    hi = ch.SurfStatsVar.Value;
else
    hi = varargin{3};
    ch.SurfStatsVar.Value = hi;
end

% initialize listbox of maps to no stats
ch.SurfStatsVarMaps.ListboxTop = 1;
ch.SurfStatsVarMaps.String = {'none (fitting the current surface)'};
ch.SurfStatsVarMaps.Value = [];
ch.SurfStatsVarMaps.Enable = 'off';

% and unset current configuration
ne_gcfg.fcfg.SurfStatsVar = struct('Filetype', 'NONE');
ne_gcfg.fcfg.SurfStatsVarIdx = [];
ne_gcfg.fcfg.stcvar = [];
ch.STCPlotChild.Visible = 'off';
ch.STCPlot.Visible = 'off';

% get currently selected surface
srf = ne_gcfg.fcfg.SurfVar;

% for no surface
if ~isxff(srf, {'fsbf', 'srf'})

    % return already
    return;
end

% definitely enable Scenery list
ch.Scenery.Enable = 'on';

% only for valid configurations
vtype = 'none';
if hi > 0 && hi <= size(hu, 1)

    % get object
    v = hu{hi, 4};

    % check whether isn't fitting to surface
    vselok = true;
    if v.NrOfVertices ~= srf.NrOfVertices

        % try to find a different SMP
        v = [];
        for smpc = 1:size(hu, 1)
            if isxff(hu{smpc, 4}, {'fsmf', 'smp'}) && ...
                hu{smpc, 4}.NrOfVertices == srf.NrOfVertices
                v = hu{smpc, 4};
                ch.SurfStatsVar.Value = smpc;
                vselok = false;
                break;
            end
        end
    end

    % set in Handles
    srf.SetHandle('Stats', {v, []});
    if isxff(v, true)
        v.SetHandle('Surface', srf);
    end

    % if the object has maps
    if isxff(v, {'fsmf', 'glm', 'mtc', 'smp'})
        vtype = lower(v.FileType);
        mnames = v.MapNames(ne_gcfg.c.extmapnames);
        ch.SurfStatsVarMaps.String = mnames;
        ch.SurfStatsVarMaps.Enable = 'on';
        ne_gcfg.fcfg.SurfStatsVar = v;
        vrtv = v.RunTimeVars;
        if numel(mnames) > 0 && ...
            vselok && ...
           ~isempty(vrtv.MapSelection{1})

            % map selection
            mnames = v.MapNames;
            if all(vrtv.MapSelection{2} <= numel(mnames)) && ...
                all(strcmp(vrtv.MapSelection{1}, mnames(vrtv.MapSelection{2})))
                msel = vrtv.MapSelection{2};
            else
                msel = multimatch(v.RunTimeVars.MapSelection, mnames);
                msel(msel == 0) = [];
            end

            % set listbox with map names, select first map and enable control
            ch.SurfStatsVarMaps.Value = msel;

            % also set in figure configuration
            ne_gcfg.fcfg.SurfStatsVarIdx = msel;
        end

        % also enable surface-time-course
        if strcmp(vtype, 'mtc')
            ne_gcfg.fcfg.stcvar = v;
        end
    end
end

% unload any referenced file
ch.SurfStatsVarRefs.Value = 1;
ne_openreffile(0, 0, 'SurfStatsVarRefs', 'SurfVar');

% projection button for maps and reflist only for RFX GLMs
stvar = ne_gcfg.fcfg.SurfStatsVar;
ch.SurfStatsVarRefs.Visible = 'off';

% general controls
cf.SetGroupEnabled('SLoadSR', 'on');
ch.SurfStatsVarRefs.Enable = 'off';

% so, do we have a valid xff?
if strcmp(vtype, 'glm') && ...
    stvar.ProjectType == 2

    % allow reference file calls
    ch.SurfStatsVarRefs.Visible = 'on';

    % get list of analyed files
    anafiles = stvar.Study;
    anafiles = {anafiles(:).NameOfAnalyzedFile};

    % populate StatsVarRefs list for VTCs
    anashort = anafiles;
    for ac = 1:numel(anashort)
        [anathr, anafile, anafext] = fileparts(anashort{ac});
        anashort{ac} = [anafile, anafext];
    end
    ch.SurfStatsVarRefs.String = [{'Inspect MTC file...'}; anashort(:)];
    ch.SurfStatsVarRefs.Value = 1;
    ch.SurfStatsVarRefs.UserData = [{[]}; anafiles(:)];
end

% set correct menu visible
if isxff(stvar, true) && any(strcmp(vtype, {'fsmf', 'glm', 'smp'}))
    if vtype(1) == 'f'
        vtype = 'smp';
    end
    cf.SetGroupVisible([upper(vtype) 'Menu'], 'on');
end

% show different page ?
if ne_gcfg.fcfg.page < 3
    ne_showpage(0, 0, 3);
end

% and then set configuration of currently selected map (if any)
ne_setcsrfstatmap;
