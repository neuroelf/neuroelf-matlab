% FUNCTION ne_btdown: react on click event (down)
function ne_btdown(varargin)

% Version:  v1.1
% Build:    16042017
% Date:     Apr-20 2016, 5:02 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, 2016, Jochen Weber
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

% only if nothing is waiting
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% take note that the button is pressed!
ne_gcfg.c.btdown = gcbf;
ne_gcfg.c.btdoup = false;
if nargin < 1 || numel(varargin{1}) ~= 1 || ~ishandle(varargin{1})
    downtest = ne_gcfg.c.btdown;
else
    downtest = varargin{1};    
end

% for main window
if downtest == ne_gcfg.h.MainFigMLH

    % record where and what modifiers at the time
    ne_gcfg.fcfg.mpos.down = ne_gcfg.h.MainFig.CurrentPoint;
    ne_gcfg.fcfg.mpos.mods = ne_gcfg.fcfg.mods;

    % update StatsVar (to allow command line changes of VMPs to be added)
    if isxff(cc.StatsVar, true)
        stvmaps = ch.StatsVarMaps.String;
        if ischar(stvmaps)
            stvmaps = cellstr(stvmaps);
        end
        if numel(cc.StatsVar.Map) ~= numel(stvmaps)
            ne_setcstats;
        end
    end

    % update?
    tcplot = strcmpi(get(ch.TCPlot.MLHandle, 'Visible'), 'on');
    if cc.page < 3
        fPos = [cc.slicepos; cc.zslicepos(1, :); cc.tcpos; cc.histpos];
        nPos = ne_gcfg.h.MainFig.CurrentPoint;
        if cc.page == 1 && ...
            any(all(nPos([1, 1, 1], :) >= (fPos(1:3, 1:2) - 1), 2) & all(nPos([1, 1, 1], :) <= (fPos(1:3, 3:4) + 1), 2))
            ne_setslicepos(0, 0, [], 'OnMouse');
        elseif cc.page == 2 && ...
            all(nPos >= (fPos(4, 1:2) - 1)) && ...
            all(nPos <= (fPos(4, 3:4) + 1))
            ne_setslicepos(0, 0, [], 'OnMouse');
        elseif tcplot && ...
            isempty(cc.mods) && ...
            all(nPos >= (fPos(5, 1:2) - 1)) && ...
            all(nPos <= (fPos(5, 3:4) + 1))
            ne_setslicepos;
        elseif all(nPos >= (fPos(6, 1:2) - 1)) && ...
            all(nPos <= (fPos(6, 3:4) + 1)) && ...
            isempty(cc.mods)
            ne_setslicepos;
        else
            ne_gcfg.c.btdown = -1;
        end
    elseif cc.page == 4
        fPos = [cc.zslicepos(1, :); cc.tcpos; cc.histpos];
        nPos = ne_gcfg.h.MainFig.CurrentPoint;
        if all(nPos >= (fPos(1, 1:2) - 1)) && ...
            all(nPos <= (fPos(1, 3:4) + 1))
            ne_render_setview(0, 0, [], 'preview');
        elseif tcplot && ...
            isempty(cc.mods) && ...
            all(nPos >= (fPos(2, 1:2) - 1)) && ...
            all(nPos <= (fPos(2, 3:4) + 1))
            ne_render_setview(0, 0, [], 'preview');
        elseif all(nPos >= (fPos(3, 1:2) - 1)) && ...
            all(nPos <= (fPos(3, 3:4) + 1)) && ...
            isempty(cc.mods)
            ne_render_setview;
        else
            ne_gcfg.c.btdown = -1;
        end
    elseif cc.page ~= 3
        ne_gcfg.c.btdown = -1;
    else
        ch.STCPlot.Visible = 'off';
        ch.STCPlotChild.Visible = 'off';
    end

% or look for satellites
elseif strcmpi(get(downtest, 'Type'), 'figure')
    try
        sats = fieldnames(ne_gcfg.cc);
        for sc = 1:numel(sats)
            if ne_gcfg.c.btdown == ne_gcfg.cc.(sats{sc}).SatelliteMLH
                ne_gcfg.cc.(sats{sc}).Config.mpos.down = ...
                    ne_gcfg.cc.(sats{sc}).Satellite.CurrentPoint;
                ne_gcfg.cc.(sats{sc}).Config.mpos.mods = ne_gcfg.fcfg.mods;
                switch (ne_gcfg.cc.(sats{sc}).Config.sattype)
                    case {'slice'}
                        ne_setsatslicepos(0, 0, ...
                            ne_gcfg.cc.(sats{sc}).Config.sattag);
                    case {'surf'}
                        ne_setsurfpos(0, 0, ne_gcfg.cc.(sats{sc}).Config.sattag);
                    case {'render'}
                        ne_render_setview(0, 0, ne_gcfg.cc.(sats{sc}).Config.sattag, 'preview');
                end
                break;
            end
        end
    catch ne_eo;
        ne_gcfg.c.lasterr = ne_eo;
    end

% patches
elseif strcmpi(get(downtest, 'Type'), 'patch')
    pt = get(downtest, 'Parent');
    if ~strcmpi(get(pt, 'Type'), 'hgtransform')
        ne_gcfg.c.incb = false;
        return;
    end
    pt = get(pt, 'Matrix');
    pu = get(downtest, 'UserData');
    pv = get(downtest, 'Vertices');
    if any(pt(1:3, 4) ~= 0)
        pv(:, 4) = 1;
    else
        pt = pt(1:3, 1:3);
    end
    pv = pv * pt';
    st = get(ne_gcfg.h.SurfaceTransform, 'Matrix');
    if any(st(1:3, 4) ~= 0) && size(pv, 2) ~= 4
        pv(:, 4) = 1;
    else
        if size(pv, 2) == 4
            pv(:, 4) = [];
        end
        st = st(1:3, 1:3);
    end
    pv = pv * st';
    if size(pv, 2) == 4
        pv(:, 4) = [];
    end
    vdm = Inf;
    vnum = [];
    if isstruct(pu) && numel(pu) == 1 && isfield(pu, 'SRF') && isxff(pu.SRF, {'fsbf', 'srf'})
        hassrf = true;
        psrf = pu.SRF;
        hsrf = psrf.Handles;
        tsrf = lower(psrf.Filetype);
    else
        hassrf = false;
    end

    % only if hit
    hit = [];
    if nargin > 1 && numel(varargin{2}) == 1 && isa(varargin{2}, 'matlab.graphics.eventdata.Hit')
        if strcmpi(varargin{2}.EventName, 'hit')
            hit = varargin{2}.IntersectionPoint;
        end
    else
        try
            dp = get(downtest, 'Parent');
            while ~strcmpi(get(dp, 'Type'), 'axes')
                dp = get(dp, 'Parent');
            end
            cpoint = get(dp, 'CurrentPoint');
            cpv = pv - ones(size(pv, 1), 1) * mean(cpoint, 1);
            cpx = cpoint(2, :) - cpoint(1, :);
            cpx = cpx(:) ./ sqrt(sum(cpx .* cpx));
            cpb = cpv * cpx;
            cpv = sum((cpv - cpb * cpx') .^ 2, 2);
            vnum = find(cpv < 4);
            if ~isempty(vnum)
                hit = true;
                cpv = cpv(vnum);
                cpb = cpb(vnum);
                vdm = minpos(cpb);
                vnum = vnum(vdm);
                vdm = cpv(vdm);
            end
        catch ne_eo;
            ne_gcfg.c.lasterror = ne_eo;
        end
    end

    % valid hit
    if ~isempty(hit)

        % see if a surface is still valid
        if hassrf

            % compute closest vertex
            if isempty(vnum)
                vd = sum((pv - hit(ones(size(pv, 1), 1), :)) .^ 2, 2);
                [vdm, vnum] = min(vd);
            end

            % if close enough
            if vdm <= 4

                if isfield(hsrf, 'VertexCoordinateOrig') && ...
                    isequal(size(hsrf.VertexCoordinateOrig), size(pv))
                    srfcoord = hsrf.VertexCoordinateOrig(vnum, :);
                else
                    srfcoord = psrf.VertexCoordinate(vnum, :);
                end
                if tsrf(1) == 's'
                    srfcoord = 128 - srfcoord(:, [3, 1, 2]);
                end

                % update linked main viewer
                ne_gcfg.fcfg.spos = {psrf, vnum, srfcoord};
                cc = ne_gcfg.fcfg;
                if ne_gcfg.c.linked && ne_gcfg.fcfg.page < 3
                    ne_setslicepos(0, 0, srfcoord);
                else
                    ne_gcfg.fcfg.cpos = round(srfcoord);
                    if isfield(hsrf, 'Stats') && iscell(hsrf.Stats) && numel(hsrf.Stats) == 2 && ...
                        numel(hsrf.Stats{1}) == 1 && isxff(hsrf.Stats{1}, {'fsmf', 'glm', 'smp'}) && ...
                       ~isempty(hsrf.Stats{2}) && get(downtest, 'Parent') == ne_gcfg.h.Surface
                        srfst = hsrf.Stats{1};
                        srfsti = hsrf.Stats{2};
                        srfv = zeros(1, numel(srfsti));
                        if isxff(srfst, {'fsmf', 'smp'})
                            for cc = 1:numel(srfv)
                                srfv(cc) = srfst.Map(srfsti(cc)).SMPData(vnum);
                            end
                        else
                            ptrfx = srfst.ProjectTypeRFX;
                            if ptrfx < 1
                                for cc = 1:numel(srfv)
                                    srfv(cc) = srfst.GLMData.BetaMaps(vnum, srfsti(cc));
                                end
                            else
                                srfsp = size(srfst.GLMData.Subject(1).BetaMaps, 2);
                                for cc = 1:numel(srfv)
                                    srfspi = 1 + mod(srfsti(cc) - 1, srfsp);
                                    srfsi = 1 + round((srfsti(cc) - srfspi) / srfsp);
                                    srfv(cc) = srfst.GLMData.Subject(srfsi).BetaMaps(vnum, srfspi);
                                end
                            end
                        end
                        if numel(srfv) > 1
                            values = 's';
                        else
                            values = '';
                        end
                        srfv = strrep(any2ascii(srfv), ' ', '');
                        viewp = sprintf('vertex: %d (%.1f, %.1f, %.1f), value%s: %s', ...
                            vnum, srfcoord(1), srfcoord(2), srfcoord(3), values, srfv);
                    else
                        viewp = sprintf('vertex: %d (%.1f, %.1f, %.1f)', ...
                            vnum, srfcoord(1), srfcoord(2), srfcoord(3));
                    end
                    ne_gcfg.h.SceneryViewPoint.String = viewp;

                    % time-course display (code in line with ne_setslicepos:910+)
                    if cc.page == 3 && numel(cc.spos{1}) == 1 && isxff(cc.spos{1}, true) && ...
                        numel(cc.stcvar) == 1 && isxff(cc.stcvar, true) && ~isempty(cc.spos{2}) && ...
                        cc.spos{1}.NrOfVertices == cc.stcvar.NrOfVertices
                        stcrtv = cc.stcvar.RunTimeVars;
                        if ~isfield(stcrtv, 'AvgMTC') || ~stcrtv.AvgMTC
                            timc = meannoinfnan(cc.stcvar.MTCData(:, cc.spos{2}), 2);
                            timcx = (0.001 * cc.stcvar.TR) .* ((1:cc.stcvar.NrOfTimePoints)' - 1);
                            if ~isempty(ch.STCPlotChildren) && ishandle(ch.STCPlotChildren(1))
                                delete(ch.STCPlotChildren);
                            end
                            ne_gcfg.h.STCPlotChildren = zeros(0, 1);
                            timx = ch.STCPlot.MLHandle;
                            hold(timx, 'on');
                            tmmm = minmaxmean(timc, 4);
                            tmmm(isinf(tmmm) | isnan(tmmm)) = 0;
                            tmd = (tmmm(2) - tmmm(1)) + 0.0001;
                            set(ch.STCPlotChild, 'XData', timcx, 'YData', double(timc(:)));
                            ch.STCPlot.XLim = [timcx(1), timcx(end)];
                            ch.STCPlot.YLim = [tmmm(1) - 0.05 * tmd, tmmm(2) + 0.05 * tmd];
                            ne_gcfg.fcfg.stcplotdata = double(timc);
                            hold(timx, 'off');
                            ch.STCPlot.Visible = 'off';
                            ch.STCPlotChild.Visible = 'on';
                            ch.MainFig.SetGroupEnabled('STCPlot', 'on');
                        end

                    % disable timecourse display
                    else
                        ch.STCPlot.Visible = 'off';
                        ch.MainFig.SetGroupEnabled('STCPlot', 'off');
                    end

                end

                % any stats is surface based GLM -> beta plotter
                wf = fieldnames(ne_gcfg.w);
                for wc = 1:numel(wf)
                    glm = ne_gcfg.w.(wf{wc});
                    if isxff(glm, 'glm') && glm.ProjectType == 2 && ...
                        glm.NrOfVertices == size(pv, 1)
                        glmh = glm.Handles;
                        if isfield(glmh, 'PlotFig') && iscell(glmh.PlotFig) && ...
                           ~isempty(glmh.PlotFig) && isxfigure(glmh.PlotFig{1}, true) && ...
                            ishandle(glmh.PlotFig{1}.MLHandle)
                            ne_glmplotbetasup(0, 0, glm, vnum);
                        end
                    end
                end
            end
        end
    end
end

% unblock
ne_gcfg.c.incb = false;
