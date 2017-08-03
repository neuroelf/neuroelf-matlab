function oo = xffmakeobject(istr)
%XFFMAKEOBJECT  Create object from struct.
%   OBJ = XFFMAKEOBJECT(ISTRUCT) creates an object out of the 1x1 input
%   struct if this is in the valid format of having a single field, with
%   the name of one of the valid types, and the content of the object.
%
%   The filename will be set to '', Handles will be reset.
%
%   If the input is invalid, the function will throw an exception.
%
%   See also XFF

% Version:  v1.1
% Build:    16012312
% Date:     Jan-23 2016, 12:18 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2016, Jochen Weber
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

% global factory
global xffsngl;

% input check
if nargin ~= 1 || ~isstruct(istr) || numel(istr) ~= 1 || numel(fieldnames(istr)) ~= 1
    error('neuroelf:xff:badArgument', 'Invalid input to xffmakeobject.');
end

% get type
ttype = fieldnames(istr);
ttype = ttype{1};

% test type
if ~isfield(xffsngl.EXT, lower(ttype))
    error('neuroelf:xff:badInputType', 'Invalid input type: %s.', upper(ttype));
end
oldcont = istr.(ttype);

% create new object of type
oo = xff(['new:' lower(ttype)]);

% not all fields the same?
if numel(fieldnames(oo.C)) ~= numel(fieldnames(oldcont)) || ...
    ~all(strcmp(fieldnames(oo.C), fieldnames(oldcont)))

    % clean up
    delete(oo);

    % then error out
    error('neuroelf:xff:badStructForClass', ...
        'Bad object type or struct given, cannot create object.');
end

% copy over
oo.C = oldcont;

% but no longer use old ID
oo.C.RunTimeVars.xffID = oo.L;
