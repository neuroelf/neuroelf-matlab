% FUNCTION ne_importrfxglm: call importrfxglm
function varargout = ne_importrfxglm(varargin)

% Version:  v1.1
% Build:    16020111
% Date:     Feb-01 2016, 11:34 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2016, Jochen Weber
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

% block any callbacks (doesn't update the screen, so this is OK)
if ne_gcfg.c.incb
    return;
end
ne_gcfg.c.incb = true;

% conversion
try
    glm = cell(1, 1);
    glm{1} = importrfxglmfromspms('ui');

    % re-allow callbacks
    ne_gcfg.c.incb = false;

    % if not a GLM, return
    if numel(glm{1}) ~= 1 || ~isxff(glm{1}, 'glm')
        return;
    end

    % browse
    ne_openfile(0, 0, glm{1}, true);

% on any error
catch ne_eo;

    % clear attempted object
    clearxffobjects(glm);

    % give error dialog
    errordlg(['Error converting SPM.mat files -> GLM: ' ne_eo.message], ...
        'NeuroElf GUI - error', 'modal');

    % re-allow callbacks
    ne_gcfg.c.incb = false;

    % and get out
    return;
end
