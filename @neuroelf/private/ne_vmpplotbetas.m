% FUNCTION ne_vmpplotbetas: initialize beta plotter for referenced GLM
function varargout = ne_vmpplotbetas(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:40 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, 2014, 2015, 2016, Jochen Weber
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

% don't do anything if current var is no usable GLM
if nargin < 3 || ...
    numel(varargin{3}) ~= 1 || ...
   ~isxff(varargin{3}, 'vmp')
    stvar = cc.StatsVar;
    if numel(stvar) == 1 && ...
        isxff(stvar, 'glm') && ...
       (stvar.ProjectTypeRFX > 0 || ...
        stvar.SeparatePredictors == 2)

        % simply pass off to GLM beta plotter
        try
            if nargout > 0
                [varargout{1:nargout}] = ne_glmplotbetas;
            else
                ne_glmplotbetas;
            end
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
        end
        return;
    end
    if numel(stvar) ~= 1 || ...
       ~isxff(stvar, 'vmp')
        return;
    end
else
    stvar = varargin{3};
end
map = stvar.Map;
if ~isfield(map, 'RunTimeVars')
    return;
end
nmaps = numel(map);
if nargin < 4 || ...
    numel(varargin{4}) ~= 1 || ...
   ~isa(varargin{4}, 'double') || ...
    isinf(varargin{4}) || ...
    isnan(varargin{4}) || ...
    varargin{4} < 1 || ...
    varargin{4} > nmaps
    stvix = cc.StatsVarIdx;
    if isempty(stvix) || ...
        stvix(1) > nmaps
        return;
    end
    stvix = stvix(1);
else
    stvix = round(varargin{4});
end

% no reference data
map = map(stvix);
rtv = map.RunTimeVars;
if ~isstruct(rtv) || ...
    numel(rtv) ~= 1 || ...
   ~isfield(rtv, 'SourceGLM') || ...
    isempty(rtv.SourceGLM) || ...
   ~isfield(rtv, 'RFXGLM') || ...
   ~islogical(rtv.RFXGLM) || ...
    numel(rtv.RFXGLM) ~= 1 || ...
   ~rtv.RFXGLM
    return;
end

% look up GLM
glmh = rtv.SourceGLMID;
try
    if ischar(glmh)
        xh = xff;
        glm = xh.Document(glmh);
    else
        glm = xff(glmh);
    end
    if ~isxff(glm, 'glm') || ...
       (glm.ProjectTypeRFX <= 0 && ...
        glm.SeparatePredictors ~= 2)
        return;
    end

% GLM not loaded
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;

    % cannot be loaded
    glmh = rtv.SourceGLM;
    if ~ischar(glmh) || ...
        exist(glmh, 'file') ~= 2

        % ask
        return;
    end

    % request loading
    answer = questdlg(sprintf('GLM ''%s'' is not loaded. Load now?', glmh), ...
        'NeuroElf - user input', 'Yes', 'No', 'Yes');
    if ischar(answer) && strcmpi(answer, 'yes')
        try
            glm = xff(glmh, 't');
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            uiwait(warndlg(['Error loading GLM: ' ne_eo.message], ...
                'NeuroElf - warning', 'modal'));
            return;
        end
    else
        return;
    end
    if ~isxff(glm, 'glm') || glm.ProjectTypeRFX <= 0
        glm.ClearObject;
        uiwait(warndlg(['Not a valid RFX-GLM: ' glmh], 'NeuroElf - warning', 'modal'));
        return;
    end

    % open file
    ne_openfile(0, 0, glm, true);

    % switch back to VMP
    stvar.Browse;

    % and this time keep track of actual GLM
    stvar.Map(stvix).RunTimeVars.SourceGLMID = glm.xffID;
end

% don't plot twice!
gh = handles(glm);
if ~isfield(gh, 'PlotFig') || ~iscell(gh.PlotFig)
    stvar.SetHandle('PlotFig', {});
    stvar.SetHandle('PlotHnd', {});
    gh = handles(stvar);
end
if ~isempty(gh.PlotFig) && ...
    numel(gh.PlotFig{1}) == 1 && ...
    isxfigure(gh.PlotFig{1}, true) && ...
   ~ne_gcfg.c.ini.BetaPlot.MultiInstance
    ph = gh.PlotHnd{1}.Axes;
    tstr = get(get(ph, 'Parent'), 'Tag');
    tstr = tstr(1:8);
    pfi = 1;

% still need to be plotted
else

    % then invoke plotting dialog
    [ph, tstr, pfig] = ne_glmplotbetas(0, 0, glm);
    pfig.Visible = 'off';
    gh = handles(glm);
    pfi = numel(gh.PlotHnd);
end

% then configure plot handle
subp = glm.SubjectPredictors;
tags = ne_gcfg.cc.(tstr).Tags;
grtv = glm.RunTimeVars;
if isfield(rtv, 'SubSel') && ...
    max(rtv.SubSel) <= numel(glm.GLMData.Subject)
    ne_gcfg.cc.(tstr).Config.subsel = rtv.SubSel(:);
else
    ne_gcfg.cc.(tstr).Config.subsel = ...
        (1:numel(ne_gcfg.cc.(tstr).Config.subids))';
end
ne_gcfg.cc.(tstr).Config.gsmapx = [];
if isfield(rtv, 'GlobSigMap') && ...
    isequal(size(rtv.GlobSigMap), size(glm.GetVolume(1)))
    ne_gcfg.cc.(tstr).Config.gsmap = rtv.GlobSigMap;
    ne_gcfg.cc.(tstr).Config.reggsx = true;
    tags.RegGS.Checked = 'on';
else
    ne_gcfg.cc.(tstr).Config.gsmap = [];
    ne_gcfg.cc.(tstr).Config.reggsx = false;
    tags.RegGS.Checked = 'off';
end
if isfield(rtv, 'Regressors') && ...
   ~isempty(rtv.Regressors) && ...
    isequal(size(rtv.Regressors), [numel(glm.GLMData.Subject), 1]) && ...
    isfield(grtv, 'CovariatesData') && ...
   ~isempty(grtv.CovariatesData) && ...
    any(all(grtv.CovariatesData == ...
        rtv.Regressors(:, ones(1, size(grtv.CovariatesData, 2))))) && ...
    isfield(rtv, 'Contrast') && ...
    isfield(grtv, 'Contrasts') && ...
   ~isempty(grtv.Contrasts) && ...
    numel(rtv.Contrast) == (numel(grtv.Contrasts{1, 2}) + 1) && ...
    any(all(cat(2, grtv.Contrasts{:, 2}) == ...
        repmat(lsqueeze(rtv.Contrast(1:end-1)), 1, size(grtv.Contrasts, 1))))
    tags.Type.Value = 2;
    tags.Contrasts.Value = findfirst(all(cat(2, grtv.Contrasts{:, 2}) == ...
        repmat(lsqueeze(rtv.Contrast(1:end-1)), 1, size(grtv.Contrasts, 1))));
    tags.Covariates.Value = findfirst(all(grtv.CovariatesData == ...
        rtv.Regressors(:, ones(1, size(grtv.CovariatesData, 2)))));
    ne_glmplotbetasgui(0, 0, tstr, 'contrasts');
    ne_glmplotbetasgui(0, 0, tstr, 'covariates');
elseif isfield(rtv, 'Regressors') && ...
   ~isempty(rtv.Regressors) && ...
    isfield(rtv, 'SubSel') && ...
    max(rtv.SubSel) <= numel(glm.GLMData.Subject) && ...
    isequal(size(rtv.Regressors), [numel(rtv.SubSel), 1]) && ...
    isfield(grtv, 'CovariatesData') && ...
   ~isempty(grtv.CovariatesData) && ...
    any(all(grtv.CovariatesData(rtv.SubSel, :) == ...
        rtv.Regressors(:, ones(1, size(grtv.CovariatesData, 2))))) && ...
    isfield(rtv, 'Contrast') && ...
    isfield(grtv, 'Contrasts') && ...
   ~isempty(grtv.Contrasts) && ...
    numel(rtv.Contrast) == (numel(grtv.Contrasts{1, 2}) + 1) && ...
    any(all(cat(2, grtv.Contrasts{:, 2}) == ...
        repmat(lsqueeze(rtv.Contrast(1:end-1)), 1, size(grtv.Contrasts, 1))))
    tags.Type.Value = 2;
    tags.Contrasts.Value = findfirst(all(cat(2, grtv.Contrasts{:, 2}) == ...
        repmat(lsqueeze(rtv.Contrast(1:end-1)), 1, size(grtv.Contrasts, 1))));
    tags.Covariates.Value = findfirst(all(grtv.CovariatesData(rtv.SubSel, :) == ...
        rtv.Regressors(:, ones(1, size(grtv.CovariatesData, 2)))));
    ne_glmplotbetasgui(0, 0, tstr, 'contrasts');
    ne_glmplotbetasgui(0, 0, tstr, 'covariates');
else
    tags.Type.Value = 1;
    if isfield(rtv, 'Contrast') && ...
        isfield(rtv, 'SubPreds') && ...
        numel(rtv.Contrast) == numel(rtv.Contrast) && ...
        any([numel(subp), numel(subp) + 1] == numel(rtv.Contrast))
        tags.Conditions.Value = lsqueeze(find(rtv.Contrast ~= 0));
    else
        tags.Conditions.Value = (1:(numel(subp)-1))';
    end
    ne_glmplotbetasgui(0, 0, tstr, 'conditions');
end
if isfield(rtv, 'Robust') && ...
    islogical(rtv.Robust) && ...
    numel(rtv.Robust) == 1 && ...
    rtv.Robust
    ne_gcfg.cc.(tstr).Config.robust = true;
    tags.Robust.Checked = 'on';
else
    ne_gcfg.cc.(tstr).Config.robust = false;
    tags.Robust.Checked = 'off';
end
ne_glmplotbetasgui(0, 0, tstr, 'type');
if isfield(rtv, 'Groups') && ...
    iscell(rtv.Groups) && ...
    isfield(grtv, 'Groups') && ...
    isequal(size(rtv.Groups), size(grtv.Groups)) && ...
    all(strcmp(rtv.Groups(:, 1), grtv.Groups(:, 1)))
    tags.DoGroups.Value = 1;
    tags.Groups.Value = (1:size(rtv.Groups, 1))';
else
    tags.DoGroups.Value = 0;
end
if isfield(rtv, 'ANCOVA') && ...
    isstruct(rtv.ANCOVA) && ...
    numel(rtv.ANCOVA) == 1
    if isfield(rtv.ANCOVA, 'Subjects') && ...
        iscell(rtv.ANCOVA.Subjects) && ...
       ~isempty(rtv.ANCOVA.Subjects)
        submatch = unique(multimatch(rtv.ANCOVA.Subjects(:), glm.Subjects));
        submatch(submatch < 1) = [];
        ne_gcfg.cc.(tstr).Config.subsel = submatch;
    end
end
ne_glmplotbetasgui(0, 0, tstr, 'dogroups');

% bring figure to foreground
gh = handles(glm);
gh.PlotFig{pfi}.Visible = 'on';
figure(gh.PlotFig{pfi}.MLHandle);
drawnow;

% set outputs
if nargout > 0
    varargout{1} = ph;
    if nargout > 1
        varargout{2} = tstr;
        if nargout > 2
            varargout{3} = gh.PlotFig{pfi}.Visible;
        end
    end
end
