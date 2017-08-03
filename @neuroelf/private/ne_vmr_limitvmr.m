% FUNCTION ne_vmr_limitvmr: re-compute 8-bit from 16-bit data
function varargout = ne_vmr_limitvmr(varargin)

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
   ~isxff(svar, 'vmr') || ...
    isempty(svar.VMRData16)
    return;
end

% get configuration
rangec = inputdlg({'Lower histogram boundary:', 'Upper histogram boundary:'}, ...
    'NeuroElf - user input', 1, {'  0.001', '  0.999'});
if ~iscell(rangec) || ...
    numel(rangec) ~= 2 || ...
    isempty(rangec{1}) || ...
    isempty(rangec{2})
    return;
end

% try conversion
range = [0, 1];
try
    range(1) = str2double(ddeblank(rangec{1}));
    range(2) = str2double(ddeblank(rangec{2}));
    if any(isinf(range) | isnan(range) | range < 0 | range > 1) || ...
            range(2) <= range(1)
        return;
    end
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg('Invalid user input.', 'NeuroElf - error message', 'modal'));
end

% set pointer to watch
ch.MainFig.Pointer = 'watch';
drawnow;

% ensure that the UndoBuffer is set
if ~isfield(svar.RunTimeVars, 'UndoBuffer') || ...
   ~isequal(size(svar.VMRData), size(svar.RunTimeVars.UndoBuffer))
    svar.RunTimeVars.UndoBuffer = svar.VMRData(:, :, :);
end

% try to call LimitVMR
try
    svar.LimitVMR(struct('range', range, 'recalc8b', true));

    % and also then call the SetScalingWindow function!
    svar.SetScalingWindow([0, 225], true);

    % and update cvar
    ne_setcvar;

    % assign to output
    varargout{1} = svar;

% handle errors
catch ne_eo;
    uiwait(warndlg(sprintf('Error using LimitVMR method: %s.', ...
        ne_eo.message), 'NeuroElf - error message', 'modal'));
end

% set pointer back to arrow
ch.MainFig.Pointer = 'arrow';
