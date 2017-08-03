% PUBLIC FUNCTION ne_extcluster: extract cluster betas from current GLM
function varargout = ne_extcluster(varargin)

% Version:  v1.0
% Build:    15040309
% Date:     Apr-03 2015, 9:55 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, 2015, Jochen Weber
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
ci = ne_gcfg.c.ini.Statistics;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% block further calls
if ne_gcfg.c.incb
    return;
end

% test for GLM
if nargin < 3
    stvar = cc.StatsVar;
    stvix = cc.StatsVarIdx;
else
    stvar = varargin{3};
    stvix = [];
end
if numel(stvar) ~= 1 || ...
   ~isxff(stvar, 'glm') || ...
   (stvar.ProjectTypeRFX < 1 && ...
    stvar.SeparatePredictors ~= 2)
    return;
end
strtv = stvar.RunTimeVars;
ne_gcfg.h.MainFig.Pointer = 'watch';
ne_gcfg.c.incb = true;
drawnow;

% get predictor names
spred = stvar.SubjectPredictors;
if ~isfield(strtv, 'ExtractConds') || ...
   ~isa(strtv.ExtractConds, 'double')
    if numel(stvix) < 2 || ...
        any(stvix > numel(spred))
        xconds = 1:numel(spred);
    else
        xconds = stvix(:)';
    end
else
    xconds = unique(min(numel(spred), max(1, round(strtv.ExtractConds(:)))))';
end
spred = spred(xconds);
if ci.ExtractWithSubIDs
    spred = sprintf('%s\t', 'Subject', spred{:});
else
    spred = sprintf('%s\t', spred{:});
end
spred(end) = [];

% get subject names
gsub = stvar.Subjects;
slen = sprintf('%%%ds', size(char(gsub), 2));

% get cluster selection
clnam = ch.Clusters.String;
if ~iscell(clnam)
    clnam = cellstr(clnam);
end
if nargin < 4 || ...
   ~isa(varargin{4}, 'double') || ...
    isempty(varargin{4}) || ...
    any(isinf(varargin{4}(:)) | isnan(varargin{4}(:)) | varargin{4}(:) < 1)
    clsel = ch.Clusters.Value;
else
    clsel = fix(varargin{4}(:));
    clsel(clsel > numel(ne_gcfg.voi.VOI)) = [];
end
clnam = clnam(clsel);

% and then betas
if ~isempty(clsel)
    betas = stvar.VOIBetas(ne_gcfg.voi, struct('vl', clsel));
else
    betas = zeros(stvar.NrOfSubjects, stvar.NrOfSubjectPredictors, 0);
end
betas = betas(:, xconds, :);

% echo
if ne_gcfg.c.echo
    ne_echo({'betas = glm.VOIBetas(voi, struct(''vl'', %s));', any2ascii(clsel)});
end

% where to output?
if nargout < 1

    % print betas first
    edigits = ci.ExtractDigits + 1;
    eformat = sprintf('%%%dg', edigits);
    betastr = cell(size(betas));
    eseps = sprintf('%c', ci.ExtractSepChars);
    for c = 1:numel(betas)
        betastr{c} = sprintf(eformat, betas(c));
        if numel(betastr{c}) > edigits
            bf = 0.1 ^ (edigits - 2);
            rf = 10 ^ (edigits - 2);
            betastr{c} = sprintf(eformat, bf * round(rf * betas(c)));
            while numel(betastr{c}) > edigits
                bf = 10 * bf;
                rf = 0.1 * rf;
                betastr{c} = sprintf(eformat, bf * round(rf * betas(c)));
            end
        end
    end

    % strip spaces
    if ~isempty(eseps) && ...
        eseps(1) ~= ' '
        betastr = reshape(ddeblank(betastr(:)), size(betastr));
    end

    % prepare text output
    for c = 1:numel(clsel)

        % one line for each subject
        bl = cell(numel(gsub), 1);
        for sc = 1:numel(bl)
            if ci.ExtractWithSubIDs
                bl{sc} = sprintf([slen '%s\n'], gsub{sc}, ...
                    sprintf([eseps '%s'], betastr{sc, :, c}));
            else
                bl{sc} = sprintf('%s\n', sprintf(['%s' eseps], betastr{sc, :, c}));
            end
        end

        % put into cluster
        clnam{c} = sprintf('%s:\n%s', clnam{c}, sprintf('%s', bl{:}));
    end

    % put into output
    ch.ClusterTable.String = sprintf('%s\n', spred, '', clnam{:});

    % also update beta plot?
    stvh = handles(stvar);
    if numel(clsel) == 1 && ...
        isfield(stvh, 'PlotFig') && ...
        iscell(stvh.PlotFig) && ...
       ~isempty(stvh.PlotFig)
        for fc = 1:numel(stvh.PlotFig)
            if numel(stvh.PlotFig{fc}) == 1 && ...
                isxfigure(stvh.PlotFig{fc}, true)
                plottstr = stvh.PlotFig{fc}.Tag;
                plottstr = plottstr(1:8);
                if isfield(ne_gcfg.cc, plottstr) && ...
                    isfield(ne_gcfg.cc.(plottstr), 'Config') && ...
                    isstruct(ne_gcfg.cc.(plottstr).Config) && ...
                    isfield(ne_gcfg.cc.(plottstr).Config, 'upplot') && ...
                   ~ne_gcfg.cc.(plottstr).Config.upplot
                    ne_gcfg.cc.(plottstr).Data.Coords = ne_gcfg.voi.VOI(clsel).Voxels;
                    ne_glmplotbetasup(0, 0, stvar, 'fromdata', plottstr);
                end
            end
        end
    end

% simply pass on output
else
    varargout{1} = betas;
end

% unblock callbacks
ne_gcfg.c.incb = false;
ne_gcfg.h.MainFig.Pointer = 'arrow';
