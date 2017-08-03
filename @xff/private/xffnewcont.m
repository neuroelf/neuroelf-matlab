function bffcont = xffnewcont(type)
%XFFNEWCONT  Generate content for XFF object.
%   C = XFFNEWCONT(TYPE) generates the content for an XFF object of type
%   TYPE.
%
%   This is a private function and shouldn't be called outside of context.
%
%   It automatically sets the RunTimeVars field to a (x)struct and also
%   sets the xffID field in the RunTimeVars field to a unique identifier.
%
%   See also BFFIO, BFFPARSE, TFFIO, TFFPARSE.

% Version:  v1.1
% Build:    16031110
% Date:     Mar-11 2016, 10:08 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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

% global variables
global ne_methods xffsngl;

% check arguments
if nargin ~= 1 || ~ischar(type) || isempty(type) || ~isfield(xffsngl.EXT, lower(type(:)'))
    error('neuroelf:xff:badArgument', 'Unknown XFF filetype.');
end

% import (for use with eval) requires neuroelf functions
emptystruct    = ne_methods.emptystruct;
uint322colcode = ne_methods.uint322colcode;
unzerostring   = ne_methods.unzerostring;
zerodstring    = ne_methods.zerodstring;

% get and prepare code to be executed
fftype = lower(type(:)');
ffspec = xffsngl.EXT.(fftype);
xfft = ffspec{1}(end-2:end);
spec = xffsngl.FF.(xfft)(ffspec{2});
bffcont = struct;
newcode = spec.NewFileCode;
if lower(spec.FFTYPE(1)) == 't'
    newcode = strrep(newcode, 'tffcont', 'bffcont');
end

% execute
try
    eval(newcode);
catch xfferror
    error('neuroelf:xff:evaluationError', ...
        'Couldn''t evaluate NewFileCode snippet for type %s: %s.', ...
        fftype, xfferror.message);
end

% make sure RunTimeVars exist
if ~isfield(bffcont, 'RunTimeVars') || ~isstruct(bffcont.RunTimeVars) || ...
    numel(bffcont.RunTimeVars) ~= 1
    bffcont.RunTimeVars = struct;
end

% generate a 24-char hex string
oo = hxdouble(randn(1, 2));
oo = oo([4:15, 20:31]);

% re-generate until it's a unique one
while any(strcmpi(xffsngl.OBJS(:, 3), oo))
    oo = hxdouble(randn(1, 2));
    oo = oo([4:15, 20:31]);
end

% then set ID
bffcont.RunTimeVars.xffID = oo;
