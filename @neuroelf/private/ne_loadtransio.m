function varargout = ne_loadtransio(varargin)
% ne_loadtransio  - load any transio fields of chosen object into memory
%
% FORMAT:       [obj = ] ne_loadtransio(SRC, EVT, object)
%
% Input fields:
%
%       SRC, EVT    Matlab handle callback inputs (discarded)
%       object      either 'StatsVar' or xff object
%
% Output fields:
%
%       obj         xff object reference
%
% Example:
%
%   neuroelf_gui('openfile', 'results.glm');
%   glm = ne_loadtransio(0, 0, 'StatsVar');

% Version:  v1.0
% Build:    14110112
% Date:     Nov-01 2014, 12:13 PM EST
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
cc = ne_gcfg.fcfg;
ch = ne_gcfg.h;

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% only with valid input
if nargin < 3 || ...
   (~ischar(varargin{3}) && ...
    (numel(varargin{3}) ~= 1 || ...
     ~isxff(varargin{3}, true)))
    ospec = 'StatsVar';
else
    ospec = varargin{3};
    if ischar(ospec) && ...
       ~strcmpi(ospec(:)', 'statsvar')
        return;
    end
end
if ischar(ospec)
    ospec = cc.StatsVar;
    if ~isxff(ospec, true)
        return;
    end
end

% load transio
mfp = ch.MainFig.Pointer;
ch.MainFig.Pointer = 'watch';
drawnow;
try
    ospec.LoadTransIOData;
    varargout{1} = ospec;
catch ne_eo;
    ne_gcfg.c.lasterr = ne_eo;
    uiwait(warndlg(['Error loading TransIOData: ' ne_eo.message], ...
        'NeuroElf - warning', 'modal'));
end
ch.MainFig.Pointer = mfp;
drawnow;
