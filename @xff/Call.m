function varargout = Call(xo, method, varargin)
%CALL  Call a method (needed for older MATLAB versions and bffio).
%   CALL(OBJ, METHOD, VARARGIN) calls method METHOD on OBJ with optional
%   arguments in VARARGIN.
%
%   [VARARGOUT] = CALL(OBJ, METHOD, VARARGIN) allows to rethrieve outputs.
%
%   See also SUBSREF

% Version:  v1.0
% Build:    16012400
% Date:     Jan-24 2016, 12:53 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2016, Jochen Weber
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

% to check availability of functions use own dir!
persistent call_methdir;
if isempty(call_methdir)
    call_methdir = [fileparts(mfilename('fullpath')) '/private/'];
end

% argument check
if nargin < 2 || numel(xo) ~= 1 || ~xffisobject(xo, true) || ~ischar(method) || isempty(method)
    error('neuroelf:xff:badArgument', 'No valid object or method given.');
end
otype = xo.S.Extensions{1};

% check all file types first
if exist([call_methdir 'aft_' method(:)' '.m'], 'file') == 2
    otype = 'aft';
end

% try call
try
    if nargout > 0
        varargout = cell(1, nargout);
        eval(['[varargout{1:nargout}]=' lower(otype) '_' method(:)' '(xo, varargin{:});']);
    else
        eval(['[varargout{1}]=' lower(otype) '_' method(:)' '(xo, varargin{:});']);
    end
catch xfferror
    rethrow(xfferror);
end
