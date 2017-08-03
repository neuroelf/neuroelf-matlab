% FUNCTION ne_setvartrf: set TrfPlus for object
function varargout = ne_setvartrf(varargin)

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-10 2011, 4:52 PM EST
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

% preset output
if nargout > 0
    varargout = cell(1, nargout);
end

% check input
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
   ~any(strcmpi(varargin{3}(:)', {'slicevar', 'statsvar', 'surfvar'}));
    return;
end

% get configuration
cc = ne_gcfg.fcfg;

% switch
switch (lower(varargin{3}(:)'))
    case {'slicevar'}
        if numel(cc.SliceVar) ~= 1 || ...
           ~isxff(cc.SliceVar, true)
            return;
        end
        if ne_gcfg.c.echo
            ne_echo('slicevar', 'SetOrientation');
        end
        try
            otrf = cc.SliceVar.SetOrientation;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            otrf = eye(4);
        end
        ne_setslicepos;
    case {'statsvar'}
        if numel(cc.StatsVar, true) ~= 1 || ...
           ~isxff(cc.StatsVar)
            return;
        end
        if ne_gcfg.c.echo
            ne_echo('statsvar', 'SetOrientation');
        end
        try
            otrf = cc.StatsVar.SetOrientation;
        catch ne_eo;
            ne_gcfg.c.lasterr = ne_eo;
            otrf = eye(4);
        end
        ne_setslicepos;
    case {'surfvar'}
end

% output
if nargout > 0
    varargout{1} = otrf;
end
