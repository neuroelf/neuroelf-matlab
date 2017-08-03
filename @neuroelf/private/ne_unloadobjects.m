% PUBLIC FUNCTION ne_unloadobjects: remove all currently loaded vars from UI
function varargout = ne_unloadobjects(varargin)

% Version:  v0.9d
% Build:    14071116
% Date:     Jul-11 2014, 4:39 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% unload all from workspace
w = ne_gcfg.w;
oids = fieldnames(w);
if ~isempty(oids)
    for oid = oids(:)'
        ne_closefile(0, 0, w.(oid{1}));
    end
end

% close renderer
if isfield(ne_gcfg.h, 'Render') && ...
    isstruct(ne_gcfg.h.Render) && ...
    isfield(ne_gcfg.h.Render, 'RendFig') && ...
    isxfigure(ne_gcfg.h.Render.RendFig, true)
    ne_render_closeui;
end

% make sure nothing is set, regardless of workspace!
ne_gcfg.fcfg.SliceVar = struct('Filetype', 'NONE');
ne_gcfg.fcfg.StatsVar = struct('Filetype', 'NONE');
ne_gcfg.fcfg.SurfVar = struct('Filetype', 'NONE');
ne_gcfg.fcfg.SurfStatsVar = struct('Filetype', 'NONE');

% update views
ne_showpage(0, 0, 3);
ne_setsurfpos(0, 0, 1);
ne_showpage(0, 0, 1);
ne_setslicepos;
