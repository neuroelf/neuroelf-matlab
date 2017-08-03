% FUNCTION ne_vmr_ihc: apply inhomogeneity correction to VMR
function varargout = ne_vmr_ihc(varargin)

% Version:  v0.9c
% Build:    11050712
% Date:     Apr-29 2011, 8:11 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, Jochen Weber
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

% initialize output
varargout = cell(1, nargout);

% test SliceVar
svar = cc.SliceVar;
if numel(svar) ~= 1 || ...
   ~isxff(svar, 'vmr')
    return;
end

% set pointer to watch
ch.MainFig.Pointer = 'watch';
drawnow;

% ensure that the UndoBuffer is set
if ~isfield(svar.RunTimeVars, 'UndoBuffer') || ...
   ~isequal(size(svar.VMRData), size(svar.RunTimeVars.UndoBuffer))
    svar.RunTimeVars.UndoBuffer = svar.VMRData(:, :, :);
end

% try to run IHC on VMR
try
    svar.InhomogeneityCorrect;

    % and also then call the SetScalingWindow function!
    svar.SetScalingWindow([0, 225], true);

    % assign to output
    varargout{1} = svar;

% handle errors
catch ne_eo;
    uiwait(warndlg(sprintf('Error using inhomogeneity correction on VMR: %s.', ...
        ne_eo.message), 'NeuroElf - error message', 'modal'));
end

% set pointer back to arrow
ch.MainFig.Pointer = 'arrow';
drawnow;
