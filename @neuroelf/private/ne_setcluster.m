% PUBLIC FUNCTION ne_setcluster: set the current pos (cpos) to cluster peak
function varargout = ne_setcluster(varargin)

% Version:  v1.0
% Build:    15042812
% Date:     Apr-28 2015, 12:22 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2015, Jochen Weber
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

try

    % what to do
    switch (lower(varargin{3}))

        % find nearest cluster (first coordinate)
        case {'nearest'}

            % don't do anything if empty
            voivx = ne_gcfg.voi.VOI;
            if isempty(voivx)
                return;
            end

            % build list
            voivx = {voivx(:).Voxels};
            for voivc = numel(voivx):-1:1
                if ~isempty(voivx{voivc})
                    voivx{voivc} = voivx{voivc}(1, :);
                else
                    voivx{voivc} = [1000, 1000, 1000];
                end
            end
            voivx = cat(1, voivx{:});

            % set to cluster with shortest distance
            voisel = minpos(sum( ...
                (voivx - ne_gcfg.fcfg.cpos(ones(size(voivx, 1), 1), :)) .^ 2, 2));
            ne_gcfg.h.Clusters.Value = voisel;
            ne_setcluster(0, 0, 'set', voisel);

        % set current cluster coordinate
        case {'set'}

            % get the current VOI(s)
            if nargin < 4 || ...
               ~isa(varargin{4}, 'double') || ...
                any(isinf(varargin{4}(:)) | isnan(varargin{4}(:)))
                csel = ne_gcfg.h.Clusters.Value;
            else
                csel = varargin{4}(:);
                ne_gcfg.h.Clusters.Value = csel;
            end
            cl = ne_gcfg.voi.VOI(csel);

            % output
            if nargout > 0
                varargout{1} = cl;
            end

            % update zoom value
            ne_gcfg.fcfg.strzoom = (ne_gcfg.h.ClusterZoom.Value > 0);

            % if only one VOI selected (and voxels present)
            if numel(cl) == 1 && ...
               ~isempty(cl.Voxels)

                % and update instant correlation?
                slvar = ne_gcfg.fcfg.SliceVar;
                isc = ne_gcfg.fcfg.instscorr;
                if isc && ...
                    isxff(slvar, 'vtc')
                    ne_instscorr(0, 0, slvar, cl.Voxels);
                end

                % set to peak (at least that's what comes from ClusterTable)
                ne_gcfg.fcfg.instscorr = false;
                if ~ne_gcfg.fcfg.strzoom
                    ne_gcfg.fcfg.cpos = cl.Voxels(1, :);
                    ne_gcfg.fcfg.strans = eye(4);
                    ne_setslicepos(0, 0, cl.Voxels(1, :), 'OnPosition');
                else
                    ne_gcfg.fcfg.cpos = [0, 0, 0];
                    ne_gcfg.fcfg.strans = ...
                        spmtrf(cl.Voxels(1, :), [0, 0, 0], [0.5, 0.5, 0.5]);
                    ne_setslicepos;
                end
                ne_gcfg.fcfg.instscorr = isc;

                % but possible update with TC cluster
                if isxff(slvar, {'hdr', 'head', 'vtc'}) && ...
                    slvar.NrOfVolumes > 1

                    % get handle
                    tcph = ne_gcfg.h.TCPlot;

                    % get time course (VOI avg)
                    timc = double(slvar.VOITimeCourse(ne_gcfg.voi, ...
                        struct('voisel', csel)));

                    % get minimum and maximum
                    tmmm = minmaxmean(timc, 4);

                    % and get the range
                    tmd = (tmmm(2) - tmmm(1)) + sqrt(eps);

                    % set to handle
                    set(tcph.Children, 'YData', timc);
                    ne_gcfg.fcfg.tcplotdata = timc;

                    % and adapt the limits
                    tcph.XLim = [0, numel(timc) + 1];
                    tcph.YLim = ...
                        [floor(tmmm(1) - 0.05 * tmd), ceil(tmmm(2) + 0.05 * tmd)];
                end

                % update GLM beta plot/s?
                plotc = fieldnames(ne_gcfg.cc);
                for pcc = 1:numel(plotc)
                    plotcc = ne_gcfg.cc.(plotc{pcc});
                    if isfield(plotcc, 'Config') && ...
                        isstruct(plotcc.Config) && ...
                        isfield(plotcc.Config, 'glm') && ...
                        isxff(plotcc.Config.glm, 'glm') && ...
                        isfield(plotcc.Config, 'upplot') && ...
                        plotcc.Config.upplot
                        try
                            ne_glmplotbetasup(0, 0, plotcc.Config.glm, cl.Voxels, plotc{pcc});
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    end
                end
            end

            % as well as extracting from GLM if present
            stvar = ne_gcfg.fcfg.StatsVar;
            stvix = ne_gcfg.h.StatsVarMaps.Value;
            extonselect = ne_gcfg.c.ini.Statistics.ExtractOnSelect;
            exttransio = ne_gcfg.c.ini.Statistics.ExtractTransIO;
            if numel(ne_gcfg.fcfg.plp) == 1 && ...
                isxff(ne_gcfg.fcfg.plp, 'plp') && ...
                ne_gcfg.c.ini.MKDA.LookupOnCluster && ...
                numel(ne_gcfg.h.Clusters.Value) == 1
                try
                    cpopt = struct;
                    prtv = ne_gcfg.fcfg.plp.RunTimeVars;
                    if isfield(prtv, 'MKDAAnalyses') && ...
                        isfield(ne_gcfg.h.MKDA, 'MKDAFig') && ...
                        isxfigure(ne_gcfg.h.MKDA.MKDAFig, true)
                        plph = ne_gcfg.h.MKDA.h;
                        cpopt.cond = plph.Points.UserData{1};
                        cpopt.unitcol = plph.StudyColumn.String{plph.StudyColumn.Value};
                    end
                    if numel(stvar) == 1 && ...
                        isxff(stvar, 'vmp') && ...
                        numel(stvix) == 1
                        stvmap = stvar.Map(stvix);
                        cpopt.gridres = stvar.Resolution;
                        cpopt.uvmp = stvar;
                        if isfield(stvmap, 'RunTimeVars') && ...
                            isstruct(stvmap.RunTimeVars) && ...
                            numel(stvmap.RunTimeVars) == 1 && ...
                            isfield(stvmap.RunTimeVars, 'MKDAMapType')
                            if isfield(stvmap.RunTimeVars, 'Condition') && ...
                               ~isempty(stvmap.RunTimeVars.Condition)
                                cpopt.cond = stvmap.RunTimeVars.Condition;
                                cpopt.unitcol = ...
                                    ne_gcfg.fcfg.plp.ColumnNames{stvmap.RunTimeVars.UnitColumn};
                                cpopt.uvmp = stvar;
                            end
                        end
                    end
                    ne_gcfg.h.ClusterTable.String = ...
                        ne_gcfg.fcfg.plp.CPointsAtCoord( ...
                        ne_gcfg.voi.VOI(ne_gcfg.h.Clusters.Value).Voxels, cpopt);
                catch ne_eo;
                    ne_gcfg.c.lasterr = ne_eo;
                end
            elseif ((strcmp(extonselect, 'single') && ...
                 numel(ne_gcfg.h.Clusters.Value) == 1) || ...
                strcmp(extonselect, 'multi')) && ...
                numel(stvar) == 1 && ...
                isxff(stvar, true)
                stvh = handles(stvar);
                stvt = filetype(stvar);
                if strcmpi(stvt, 'glm') && ...
                    stvar.ProjectType == 1 && ...
                    stvar.ProjectTypeRFX > 0 && ...
                   (~istransio(stvar.GLMData.Subject(1).BetaMaps) || ...
                    exttransio)
                    ne_extcluster;
                elseif strcmpi(stvt, 'vmp') && ...
                    isfield(stvh, 'SourceGLM') && ...
                    isxff(stvh.SourceGLM, 'glm')
                    if exttransio || ...
                       (stvh.SourceGLM.ProjectTypeRFX > 0 && ...
                        ~istransio(stvh.SourceGLM.GLMData.Subject(1).BetaMaps)) || ...
                       (stvh.SourceGLM.ProjectTypeRFX <= 0 && ...
                        ~istransio(stvh.SourceGLM.GLMData.BetaMaps))
                        ne_extcluster(0, 0, stvh.SourceGLM);
                    end
                elseif strcmpi(stvt, 'vmp') && ...
                    numel(stvix) == 1
                    stvr = stvar.Map(stvix).RunTimeVars;
                    if isfield(stvr, 'SourceGLM') && ...
                        ischar(stvr.SourceGLM) && ...
                       ~isempty(stvr.SourceGLM) && ...
                        isfield(stvr, 'SourceGLMID') && ...
                        ischar(stvr.SourceGLMID) && ...
                       ~isempty(stvr.SourceGLMID)
                        [mglmp, mglmf] = fileparts(stvr.SourceGLM);
                        xffroot = xff;
                        glms = xffroot.Documents('glm');
                        if isempty(glms)
                            return;
                        end
                        gmatch = [];
                        for gc = 1:numel(glms)
                            if ischar(glms{gc})
                                [glmp, glmf] = fileparts(glms{gc});
                                if strcmpi(glmf, mglmf)
                                    gmatch = xffroot.Document(glms{gc});
                                    break;
                                end
                            else
                                tglm = xffroot.Document(glms{gc});
                                if strcmpi(tglm.RunTimeVars.xffID, stvr.SourceGLMID)
                                    gmatch = tglm;
                                    break;
                                end
                            end
                        end
                        if ~isempty(gmatch)
                            if exttransio || ...
                               (gmatch.ProjectTypeRFX > 0 && ...
                                ~istransio(gmatch.GLMData.Subject(1).BetaMaps)) || ...
                               (gmatch.ProjectTypeRFX <= 0 && ...
                                ~istransio(gmatch.GLMData.BetaMaps))
                                ne_extcluster(0, 0, gmatch);
                            end
                        elseif numel(ne_gcfg.fcfg.glm) == 1 && ...
                            isxff(ne_gcfg.fcfg.glm, 'glm')
                            gmatch = ne_gcfg.fcfg.glm;
                            if exttransio || ...
                               (gmatch.ProjectTypeRFX > 0 && ...
                                ~istransio(gmatch.GLMData.Subject(1).BetaMaps)) || ...
                               (gmatch.ProjectTypeRFX <= 0 && ...
                                ~istransio(gmatch.GLMData.BetaMaps))
                                ne_extcluster(0, 0, gmatch);
                            end
                        end
                    end
                elseif numel(ne_gcfg.fcfg.glm) == 1 && ...
                    isxff(ne_gcfg.fcfg.glm, 'glm')
                    gmatch = ne_gcfg.fcfg.glm;
                    if exttransio || ...
                       (gmatch.ProjectTypeRFX > 0 && ...
                        ~istransio(gmatch.GLMData.Subject(1).BetaMaps)) || ...
                       (gmatch.ProjectTypeRFX <= 0 && ...
                        ~istransio(gmatch.GLMData.BetaMaps))
                        ne_extcluster(0, 0, gmatch);
                    end
                end
            end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
