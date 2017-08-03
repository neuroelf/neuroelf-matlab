% FUNCTION ne_setcstats: set current StatsVar into figure configuration
function varargout = ne_setcstats(varargin)

% Version:  v1.1
% Build:    16052821
% Date:     May-28 2016, 9:24 PM EST
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

% global variable
global ne_gcfg;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% get handles
ch = ne_gcfg.h;
cf = ch.MainFig;

% set menu groups
cf.SetGroupVisible('GLMMenu', 'off');
cf.SetGroupVisible('VMPMenu', 'off');

% and some controls
cf.SetGroupEnabled('SLdNVMP', 'off');

% preset title
ne_gcfg.c.title(1, 2:3) = {'', ''};
ne_updatename;

% hide stats texts
set(ch.CorStatsText, 'Visible', 'off', 'String', ' ');
set(ch.ZoomStatsText, 'Visible', 'off', 'String', ' ');

% get selection and UserData (list of objects) of StatsVar
hu = ch.StatsVar.UserData;
if nargin < 3 || ~isa(varargin{3}, 'double') || numel(varargin{3}) ~= 1 || ...
    isinf(varargin{3}) || isnan(varargin{3}) || varargin{3} < 1 || ...
    varargin{3} > size(hu, 1) || varargin{3} ~= fix(varargin{3})
    hi = ch.StatsVar.Value;
else
    hi = varargin{3};
    ch.StatsVar.Value = hi;
end

% initialize listbox of maps to no stats
ch.StatsVarMaps.ListboxTop = 1;
ch.StatsVarMaps.Value = [];
ch.StatsVarMaps.String = {'none'};
ch.StatsVarMaps.Enable = 'off';
ch.TCMovieMenu.Enable = 'off';

% by default, hide TC-window if tcvar is an AvgVTC
ccf = fieldnames(ne_gcfg.cc);
if ~isempty(ccf)
    ccf(cellfun('isempty', regexp(ccf, '^TC'))) = [];
end
if ~isempty(ccf) && numel(ne_gcfg.cc.(ccf{1})) == 1 && isstruct(ne_gcfg.cc.(ccf{1})) && ...
    isfield(ne_gcfg.cc.(ccf{1}), 'Satellite') && isxfigure(ne_gcfg.cc.(ccf{1}).Satellite, true)
    tcvar = ne_gcfg.cc.(ccf{1}).Config.tcvar;
    if isxff(tcvar, 'vtc') && isfield(tcvar.RunTimeVars, 'AvgVTC') && ...
        islogical(tcvar.RunTimeVars.AvgVTC) && tcvar.RunTimeVars.AvgVTC
        ne_gcfg.cc.(ccf{1}).Satellite.Visible = 'off';
    end
end

% and unset current configuration
ne_gcfg.fcfg.StatsVar = [];
ne_gcfg.fcfg.StatsVarIdx = [];

% disable instantaneous correlations
ne_setoption(0, 0, 'instseedcorr', false);

% only for valid configurations
if hi > 0 && hi <= size(hu, 1)

    % get object
    v = hu{hi, 4};
    ne_gcfg.fcfg.StatsVar = v;

    % if the object has maps
    mnames = v.MapNames(ne_gcfg.c.extmapnames);
    ch.StatsVarMaps.String = mnames;
    ch.StatsVarMaps.Enable = 'on';
    vrtv = v.RunTimeVars;
    if numel(mnames) > 0 && ~isempty(vrtv.MapSelection{1})

        % map selection
        mnames = v.MapNames;
        if all(vrtv.MapSelection{2} <= numel(mnames)) && ...
            all(strcmp(vrtv.MapSelection{1}(:), mnames(vrtv.MapSelection{2})))
            msel = vrtv.MapSelection{2};
        else
            msel = multimatch(v.RunTimeVars.MapSelection, mnames);
            msel(msel == 0) = [];
        end

        % set listbox with map names, select first map and enable control
        ch.StatsVarMaps.Value = msel;

        % also set in figure configuration
        ne_gcfg.fcfg.StatsVarIdx = msel;
    end
else
    ne_gcfg.fcfg.StatsVar = struct('Filetype', 'NONE');
end

% projection button for maps and reflist only for RFX GLMs
stvar = ne_gcfg.fcfg.StatsVar;

% so, do we have a valid (VTC based) xff?
if isxff(stvar, 'glm') && stvar.ProjectType == 1

    % set as current object
    ne_gcfg.fcfg.glm = stvar;

    % onload any referenced file
    if numel(ne_gcfg.fcfg.StatsVarRefObj) == 1 && ...
        isxff(ne_gcfg.fcfg.StatsVarRefObj, 'glm') && ...
        ne_gcfg.fcfg.StatsVarRefObj ~= stvar
        ch.StatsVarRefs.Value = 1;
        ne_openreffile(0, 0, 'StatsVarRefs', 'SliceVar');
    end
    ne_gcfg.fcfg.StatsVarRefObj = stvar;

    % get list of analyed files
    anafiles = stvar.Study;
    anafiles = {anafiles(:).NameOfAnalyzedFile};

    % populate StatsVarRefs list for VTCs
    anashort = anafiles;
    for ac = 1:numel(anashort)
        [anathr, anafile, anafext] = fileparts(anashort{ac});
        anashort{ac} = [anafile, anafext];
    end
    ch.StatsVarRefs.String = [{'Inspect VTC file...'}; anashort(:)];
    ch.StatsVarRefs.UserData = [{[]}; anafiles(:); {stvar}];

% for average VTCs
elseif isxff(stvar, 'vtc') && isfield(vrtv, 'AvgVTC') && islogical(vrtv.AvgVTC) && vrtv.AvgVTC
    ch.TCMovieMenu.Enable = 'on';
end

% set correct menu visible
if isxff(stvar, true) && any(strcmpi(stvar.Filetype, {'glm', 'vmp'}))
    cf.SetGroupVisible([upper(stvar.Filetype) 'Menu'], 'on');
end

% also update in CM if open
if isxff(stvar, 'glm') && stvar.ProjectTypeRFX > 0 && ...
    isfield(ch, 'CM') && isstruct(ch.CM) && ...
    isfield(ne_gcfg.fcfg, 'CM') && isstruct(ne_gcfg.fcfg.CM)

    % lookup
    for glmc = 1:numel(ne_gcfg.fcfg.CM.GLMs)
        if stvar == ne_gcfg.fcfg.CM.GLMs{glmc}
            ch.CM.h.GLMs.Value = glmc;
            ne_cm_setglm;
            break;
        end
    end
end

% also update in RM if open
if isxff(stvar, 'glm') && stvar.ProjectTypeRFX > 0 && ...
    isfield(ch, 'RM') && isstruct(ch.RM) && ...
    isfield(ne_gcfg.fcfg, 'RM') && isstruct(ne_gcfg.fcfg.RM)

    % lookup
    for glmc = 1:numel(ne_gcfg.fcfg.RM.GLMs)
        if stvar == ne_gcfg.fcfg.RM.GLMs{glmc}
            ch.RM.h.GLMs.Value = glmc;
            ne_rm_setglm;
            break;
        end
    end
end

% show different page ?
if ne_gcfg.fcfg.page > 2
    ne_showpage(0, 0, 1);
end

% and then set configuration of currently selected map (if any)
ne_setcstatmap;
