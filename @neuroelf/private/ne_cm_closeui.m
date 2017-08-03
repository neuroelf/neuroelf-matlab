% FUNCTION ne_cm_closeui: close contrast manager window
function ne_cm_closeui(varargin)

% Version:  v1.1
% Build:    16041210
% Date:     Apr-12 2016, 10:27 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2016, Jochen Weber
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
hFig = ne_gcfg.h.CM.CMFig;

% update configuration in current GLM's RunTimeVars
if numel(hFig.UserData.lastglm) == 1 && isxff(hFig.UserData.lastglm, 'glm')
    ne_cm_updatertv(0, 0, hFig.UserData.lastglm);
end

% hide figure
hFig.Visible = 'off';
mfp = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';
drawnow;

% save RunTimeVars for GLMs
cprog = ne_progress(0, 0, {true, 0, 'Saving GLM RunTimeVars...'});
glms = ne_gcfg.fcfg.CM.GLMs;
for gc = 1:numel(glms)
    if ~glms{gc}.Handles.RunTimeVarsSaved
        glms{gc}.SaveRunTimeVars;
        ch.Progress.Progress(gc/numel(glms));
    end
end

% update last known position
ne_gcfg.c.ini.Children.CMPosition = hFig.Position(1:2);

% delete figure and remove from global struct
hFig.Delete;
ne_gcfg.h.CM = [];
ne_progress(0, 0, cprog);

% release contrast manager call
ne_gcfg.c.blockcb(strcmp(ne_gcfg.c.blockcb, 'cm_open')) = [];
ch.MainFig.Pointer = mfp;
drawnow;
