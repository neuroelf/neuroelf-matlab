% FUNCTION ne_btup: take note that the mouse button is released
function ne_btup(varargin)

% Version:  v1.0
% Build:    15032017
% Date:     Mar-20 2015, 5:19 PM EST
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

% set the button state change flag
ne_gcfg.c.btdoup = true;

% return if in another callback
if ne_gcfg.c.incb
    return;
end

% main figure
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~isfield(ne_gcfg.cc, varargin{3})

    % for render
    if ne_gcfg.fcfg.page == 4 && ...
        ne_gcfg.c.rpreview
        ne_render_setview;
    end

    % then clear config
    ne_gcfg.fcfg.histset = 0;
    ne_gcfg.fcfg.mpos.ddat = {};
    ne_gcfg.fcfg.mpos.mods = ne_gcfg.fcfg.mods;

% or satellite/child
else
    if strcmpi(ne_gcfg.cc.(varargin{3}).Config.sattype, 'render')
        ne_render_setview(0, 0, varargin{3});
    end
    ne_gcfg.cc.(varargin{3}).Config.mpos.ddat = {};
end

% set the actual button flag
ne_gcfg.c.btdown = [];
