function varargout = ne_closefile(varargin)
% ne_closefile  - remove file from workspace and, possibly, clear
%
% FORMAT:       ne_closefile(SRC, EVT, filespec [, finalcall])
%
% Input fields:
%
%       filespec    either of 'SliceVar', 'StatsVar', 'SurfVar', or
%                   'SurfStatsVar'; or a 1x1 (xff) object that is loaded
%       finalcall   if set to 'final' update is suppressed
%
% Example:
%
%   ne_closefile(0, 0, 'SliceVar');

% Version:  v1.1
% Build:    16031900
% Date:     Mar-19 2016, 12:37 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, 2016, Jochen Weber
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
ci = ne_gcfg.c.ini;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
    varargout{1} = false;
end

% don't bother with calls without suitable arguments
if nargin < 3
    return;
end
f = varargin{3};
if nargin > 3 && ...
    ischar(varargin{4}) && ...
    strcmpi(varargin{4}(:)', 'final')
    doupdate = false;
else
    doupdate = true;
end

% for char argument
if ischar(f)

    % if call from SliceVar
    if strcmpi(f, 'slicevar')

        % try to get object from UserData
        try
            f = ch.SliceVar.UserData{ch.SliceVar.Value, 4};
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % if call from StatsVar
    elseif strcmpi(f, 'statsvar')

        % try to get object from UserData
        try
            f = ch.StatsVar.UserData{ch.StatsVar.Value, 4};
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % if call from SurfVar
    elseif strcmpi(f, 'surfvar')

        % try to get object from UserData
        try
            f = ch.SurfVar.UserData{ch.SurfVar.Value, 4};
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end

    % if call from SurfStatsVar
    elseif strcmpi(f, 'surfstatsvar')

        % try to get object from UserData
        try
            f = ch.SurfStatsVar.UserData{ch.SurfStatsVar.Value, 4};
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
end

% check if valid object found
if numel(f) ~= 1 || ~isxff(f, true)
    return;
end
ftyp = lower(f.Filetype);

% extend filetype
if (strcmp(ftyp, 'hdr') && ...
    (~isempty(regexpi(f.DataHist.Description, ...
            'spm\{[ft]_\[\d+(\.\d+)?(,\s*[1-9][0-9\.]*)?\]\}')) || ...
      (isfield(f.RunTimeVars, 'StatsObject') && ...
       islogical(f.RunTimeVars.StatsObject) && ...
       numel(f.RunTimeVars.StatsObject) == 1 && ...
       f.RunTimeVars.StatsObject))) || ...
   (strcmp(ftyp, 'head') && ...
    ~isempty(strfind(f.TypeOfVolumes, '_FUNC')))
    ftyp = 'hfunc';

% for average VTCs
elseif strcmp(ftyp, 'vtc') && ...
    isfield(f.RunTimeVars, 'AvgVTC') && ...
    f.RunTimeVars.AvgVTC
    ftyp = 'atc';
end

% depending on filetype
fh = handles(f);
if any(strcmp(ftyp, {'dmr', 'fmr', 'hdr', 'head', 'mgh', 'msk', 'vmr', 'vtc'})) || ...
   (strcmp(ftyp, 'cmp') && ...
    f.DocumentType == 0)

    % remove from either SliceVar control and re-set current slice object
    cv_removefromlist(ch.SliceVar, f);
    if doupdate
        ne_setcvar;
    end

% or
elseif any(strcmp(ftyp, {'atc', 'ava', 'glm', 'hfunc', 'map', 'vmp'})) || ...
   (strcmp(ftyp, 'cmp') && ...
    f.DocumentType == 1)

    % from StatsVar control and re-set current stats object
    cv_removefromlist(ch.StatsVar, f);

    % for GLMs
    if strcmp(ftyp, 'glm')

        % close plot figures
        if isfield(fh, 'PlotFig') && ...
            iscell(fh.PlotFig) && ...
           ~isempty(fh.PlotFig)
            for fc = numel(fh.PlotFig):-1:1
                if isxfigure(fh.PlotFig{fc}, true)
                    try
                        ne_glmplotbetasgui(0, 0, fh.PlotFig{fc}.Tag(1:8), 'close');
                    catch ne_eo;
                        ne_gcfg.c.lasterr = ne_eo;
                    end
                end
            end
        end

        % check whether this file is the current selection in AC/CM/RM
        if isstruct(ch.AC) && ...
            isfield(ch.AC, 'ACFig') && ...
            isxfigure(ch.AC.ACFig, true)

            % this is the selection?
            if isfield(ch.AC.ACFig.UserData, 'lastglm') && ...
                isxff(ch.AC.ACFig.UserData.lastglm) && ...
                f == ch.AC.ACFig.UserData.lastglm

                % unset last GLM (no more saves to RTV!)
                ch.AC.ACFig.UserData.lastglm = [];

                % and close
                ne_ancova_ui(0, 0, 'closeui');

            % some other GLM is the selection
            else

                % find GLM in list
                for gc = 1:numel(ne_gcfg.fcfg.AC.GLMs)

                    % when found
                    if f == ne_gcfg.fcfg.AC.GLMs{gc}

                        % remove from list
                        ne_gcfg.fcfg.AC.GLMs(gc) = [];

                        % and string of control
                        gcs = ne_gcfg.h.AC.h.GLMs.String;
                        if ~iscell(gcs)
                            gcs = cellstr(gcs);
                        end
                        gcs(gc) = [];
                        if ne_gcfg.h.AC.h.GLMs.Value > gc
                            ne_gcfg.h.AC.h.GLMs.Value = ne_gcfg.h.AC.h.GLMs.Value - 1;
                        end
                        ne_gcfg.h.AC.h.GLMs.String = gcs;
                        break;
                    end
                end
            end
        end
        if isstruct(ch.CM) && ...
            isfield(ch.CM, 'CMFig') && ...
            isxfigure(ch.CM.CMFig, true)

            % this is the selection?
            if isfield(ch.CM.CMFig.UserData, 'lastglm') && ...
                isxff(ch.CM.CMFig.UserData.lastglm) && ...
                f == ch.CM.CMFig.UserData.lastglm

                % unset last GLM (no more saves to RTV!)
                ch.CM.CMFig.UserData.lastglm = [];

                % and close
                ne_cm_closeui;

            % some other GLM is the selection
            else

                % find GLM in list
                for gc = 1:numel(ne_gcfg.fcfg.CM.GLMs)

                    % when found
                    if f == ne_gcfg.fcfg.CM.GLMs{gc}

                        % remove from list
                        ne_gcfg.fcfg.CM.GLMs(gc) = [];

                        % and string of control
                        gcs = ne_gcfg.h.CM.h.GLMs.String;
                        if ~iscell(gcs)
                            gcs = cellstr(gcs);
                        end
                        gcs(gc) = [];
                        if ne_gcfg.h.CM.h.GLMs.Value > gc
                            ne_gcfg.h.CM.h.GLMs.Value = ne_gcfg.h.CM.h.GLMs.Value - 1;
                        end
                        ne_gcfg.h.CM.h.GLMs.String = gcs;
                        break;
                    end
                end
            end
        end
        if isstruct(ch.RM) && ...
            isfield(ch.RM, 'RMFig') && ...
            isxfigure(ch.RM.RMFig, true)
            if isfield(ch.RM.RMFig.UserData, 'lastglm') && ...
                isxff(ch.RM.RMFig.UserData.lastglm) && ...
                f == ch.RM.RMFig.UserData.lastglm
                ch.RM.RMFig.UserData.lastglm = [];
                ne_rm_closeui;
            else
                for gc = 1:numel(ne_gcfg.fcfg.RM.GLMs)
                    if f == ne_gcfg.fcfg.RM.GLMs{gc}
                        ne_gcfg.fcfg.RM.GLMs(gc) = [];
                        gcs = ne_gcfg.h.RM.h.GLMs.String;
                        if ~iscell(gcs)
                            gcs = cellstr(gcs);
                        end
                        gcs(gc) = [];
                        if ne_gcfg.h.RM.h.GLMs.Value > gc
                            ne_gcfg.h.RM.h.GLMs.Value = ne_gcfg.h.RM.h.GLMs.Value - 1;
                        end
                        ne_gcfg.h.RM.h.GLMs.String = gcs;
                        break;
                    end
                end
            end
        end
    end

    % update maps list, etc.
    if doupdate
        ne_setcstats;
    end

% or yet
elseif any(strcmp(ftyp, {'fsbf', 'srf'}))

    % from SurfVar control
    cv_removefromlist(ch.SurfVar, f);
    cv_removefromlist(ch.Scenery, f, 'Scene compilation (empty)');
    if isfield(fh, 'Surface') && ...
        ishandle(fh.Surface)
        delete(fh.Surface);
    end
    if isfield(fh, 'SurfaceTransform') && ...
        ishandle(fh.SurfaceTransform)
        delete(fh.SurfaceTransform);
    end
    try
        if isfield(fh, 'SurfTIO')
            delete(fh.SurfTIO);
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end
    if doupdate
        ne_setcsrf;
    end
    for dh = {'ShownInGUI', 'Surface', 'SurfStatsBars', 'SurfTIO', 'SurfaceTransform', ...
            'SUpdate', 'VertexMorphIndex', 'VertexMorphMeshes', ...
            'VertexCoordinateTal', 'VertexNormalTal'}
        try
            if isfield(fh, dh{1})
                f.DeleteHandle(dh{1});
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
    end
    if doupdate
        if ~isxff(ne_gcfg.fcfg.SurfVar, true)
            ne_showpage(0, 0, 1);
        else
            ne_setsurfpos;
        end
        if ~isempty(ch.SurfVar.UserData)
            ch.Scenery.Enable = 'on';
        end
    end

% SMP
elseif any(strcmp(ftyp, {'fsmf', 'mtc', 'smp'}))

    % from SurfStatsVar
    cv_removefromlist(ch.SurfStatsVar, f);
    if doupdate
        ne_setcsrfstats;
        if ~isempty(ch.SurfVar.UserData)
            ch.Scenery.Enable = 'on';
        end
    end
end

% clear-anyway flag
if nargin > 3 && ...
    numel(varargin{4}) == 1 && ...
    islogical(varargin{4})
    clearany = varargin{4};
else
    clearany = ci.MainFig.ClearClosed;
end

% look through workspace
wf = fieldnames(ne_gcfg.w);
for tc = 1:numel(wf)

    % get struct of object in workspace
    tstr = ne_gcfg.w.(wf{tc});

    % check handle equality
    if tstr == f

        % if loaded by tool
        if nargin < 4 && ne_gcfg.wc.(wf{tc})

            % then clear object
            clearany = true;
        end

        % then remove from workspace and control
        ne_gcfg.w = rmfield(ne_gcfg.w, wf{tc});
        ne_gcfg.wc = rmfield(ne_gcfg.wc, wf{tc});
    end
end

% clear anyway
if clearany && isxff(f, true)
    f.ClearObject;
    if nargout > 0
        varargout{1} = true;
    end
end

% update sliceing
if doupdate
    ne_setslicepos;
end
