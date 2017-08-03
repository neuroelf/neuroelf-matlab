% FUNCTION ne_render_closeui: close Render UI
function ne_render_closeui(varargin)

% Version:  v1.0
% Build:    14103017
% Date:     Oct-30 2014, 5:49 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2014, Jochen Weber
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

% force to first page
ne_showpage(0, 0, 1);

% try to store last known position
try
    ne_gcfg.c.ini.Children.RenderPosition = ...
        ne_gcfg.h.Render.RendFig.Position(1:2);
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% try to close window
try
    ne_gcfg.h.Render.RendFig.Delete;
    ne_gcfg.h.Render = [];
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% remove maps
try
    o = ne_gcfg.fcfg.Render;
    ne_gcfg.fcfg.Render = [];
    if isxff(o.stvar, true) && ...
       (~islogical(o.smstatk) || ...
        numel(o.smstatk) ~= 1 || ...
        ~o.smstatk)
        remm = setdiff(o.stvix(:)', o.stvixo(:)');
        if ~isempty(remm)
            o.stvar.Map(remm) = [];
            switch (lower(o.stvar.Filetype))
                case {'cmp', 'vmp'}
                    o.stvar.NrOfMaps = numel(o.stvar.Map);
                case {'hdr'}
                    o.stvar.VoxelData(:, :, :, remm) = [];
                    o.stvar.VoxelDataCT(remm) = [];
                case {'head'}
                    o.stvar.Brick(remm) = [];
            end
            if isxff(ne_gcfg.fcfg.StatsVar, true) && ...
                ne_gcfg.fcfg.StatsVar == o.stvar
                ne_setcstatmap(0, 0, o.stvixo(:)');
                ne_openfile(0, 0, o.stvar);
            end

        % either way
        elseif isxff(ne_gcfg.fcfg.StatsVar, true) && ...
            ne_gcfg.fcfg.StatsVar == o.stvar

            % update list!
            ne_openfile(0, 0, o.stvar);
        end
    end
    if isxff(o.slvar, true)
        slvarh = handles(o.slvar);
        for dh = {'RenderAVol', 'RenderMView', 'RenderPreview', 'RenderPView', 'RenderSVol'}
            if isfield(slvarh, dh{1})
                o.slvar.DeleteHandle(dh{1});
            end
        end
        if isfield(o.slvar.RunTimeVars, 'SliceRanges')
            o.slvar.RunTimeVars = rmfield(o.slvar.RunTimeVars, 'SliceRanges');
        end
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
end

% remove struct
ne_gcfg.fcfg.Render = [];
ne_gcfg.h.Render = [];
