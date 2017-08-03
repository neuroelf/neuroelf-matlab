% PUBLIC FUNCTION ne_setsrfcluster: set the current pos (cpos) to POI peak
function varargout = ne_setsrfcluster(varargin)

% Version:  v1.1
% Build:    16031220
% Date:     Mar-12 2016, 8:44 PM EST
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
ch = ne_gcfg.h;

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
            poivx = ne_gcfg.poi.POI;
            if isempty(poivx)
                return;
            end

            % check last selection
            cc = ne_gcfg.fcfg;
            if numel(cc.spos{1}) == 1 && isxff(cc.spos{1}, 'srf') && ...
                numel(cc.spos{2}) == 1 && isa(cc.spos{2}, 'double') && ...
               ~isinf(cc.spos{2}) && ~isnan(cc.spos{2}) && cc.spos{2} > 0
                csrf = cc.spos{1};
                vnum = cc.spos{2};
                vpos = cc.spos{3};
            else
                return;
            end

            % build list
            poivx = {poivx(:).Vertices};
            for poivc = numel(poivx):-1:1
                if ~isempty(poivx{poivc})
                    if any(poivx{poivc} == vnum)
                        ch.ClustersSrf.Value = poivc;
                        ne_setsrfcluster(0, 0, 'set', poivc);
                        return;
                    end
                    poivx{poivc} = poivx{poivc}(1);
                else
                    poivx{poivc} = -1;
                end
            end
            poivx = cat(1, poivx{:});
            hsrf = csrf.Handles;
            if isfield(hsrf, 'VertexCoordinateOrig')
                srfcoord = 128 - hsrf.VertexCoordinateOrig(poivx, [3, 1, 2]);
            else
                srfcoord = 128 - csrf.VertexCoordinate(poivx, [3, 1, 2]);
            end
            mdist = sum((srfcoord - vpos(ones(size(srfcoord, 1), 1), :)) .^ 2, 2);

            % set to cluster with shortest distance
            poisel = minpos(mdist);
            ch.ClustersSrf.Value = poisel;
            ne_gcfg.fcfg.spos = {csrf, poivx(poisel), srfcoord(poisel, :)};
            ne_setsrfcluster(0, 0, 'set', poisel);

        % set current cluster coordinate
        case {'set'}

            % get the current VOI(s)
            cc = ne_gcfg.fcfg;
            if numel(cc.spos{1}) ~= 1 || ~isxff(cc.spos{1}, 'srf')
                return;
            end
            csrf = cc.spos{1};
            if nargin < 4 || ~isa(varargin{4}, 'double') || ...
                any(isinf(varargin{4}(:)) | isnan(varargin{4}(:)))
                csel = ne_gcfg.h.ClustersSrf.Value;
            else
                csel = varargin{4}(:);
                ne_gcfg.h.ClustersSrf.Value = csel;
            end
            cl = ne_gcfg.poi.POI(csel);

            % output
            if nargout > 0
                varargout{1} = cl;
            end

            % if only one VOI selected (and voxels present)
            if numel(cl) == 1 && ...
               ~isempty(cl.Vertices)

                % update GLM beta plot/s?
                plotc = fieldnames(ne_gcfg.cc);
                for pcc = 1:numel(plotc)
                    plotcc = ne_gcfg.cc.(plotc{pcc});
                    if isfield(plotcc, 'Config') && ...
                        isstruct(plotcc.Config) && ...
                        isfield(plotcc.Config, 'glm') && ...
                        isxff(plotcc.Config.glm, 'glm') && ...
                        plotcc.Config.glm.ProjectType == 2 && ...
                        isfield(plotcc.Config, 'upplot') && ...
                        plotcc.Config.upplot
                        try
                            ne_glmplotbetasup(0, 0, plotcc.Config.glm, cl.Vertices, plotc{pcc});
                        catch ne_eo;
                            ne_gcfg.c.lasterr = ne_eo;
                        end
                    end
                end
                hsrf = csrf.Handles;
                if isfield(hsrf, 'VertexCoordinateOrig')
                    srfcoord = 128 - hsrf.VertexCoordinateOrig(cl.Vertices(1), [3, 1, 2]);
                else
                    srfcoord = 128 - csrf.VertexCoordinate(cl.Vertices(1), [3, 1, 2]);
                end
                viewp = sprintf('vertex: %d (%.1f, %.1f, %.1f)', ...
                    cl.Vertices(1), srfcoord(1), srfcoord(2), srfcoord(3));
                ch.SceneryViewPoint.String = viewp;
            end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end
