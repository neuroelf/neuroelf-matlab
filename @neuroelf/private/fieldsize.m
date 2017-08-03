function fsz = fieldsize(str, fname)
% fieldsize  - get size (in bytes) of struct fields
%
% FORMAT:       fsz = fieldsize(str [, fname])
%
% Input fields:
%
%       str         1x1 struct to get fieldsize(s) from
%       fname       optional fieldname, if not given or not char get sizes
%                   of all fields and return as struct
%
% Output fields:
%
%       fsz         either 1x1 double with size (one field) or 1x1 struct
%                   with sizes of fields (same order)

% Version:  v0.9b
% Build:    10062423
% Date:     Jun-24 2010, 5:52 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, Jochen Weber
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

% argument check
if nargin < 1 || ...
   ~isstruct(str) || ...
    numel(str) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% all fields
if nargin < 2 || ...
   ~ischar(fname)

    % create copy
    fsz = str;

    % all fields
    ff = fieldnames(fsz);
    for fc = 1:numel(ff)

        fsz.(ff{fc}) = fsvarsize(fsz.(ff{fc}));
    end

% specific field
else

    % field doesn't exist
    if isempty(fname) || ...
       ~isfield(str, fname(:)')
        error( ...
            'neuroelf:BadArgument', ...
            'Requested field doesn''t exist.' ...
        );
    end

    % get size
    fsz = fsvarsize(str.(fname(:)'));
end


% sub function fsvarsize
function sz = fsvarsize(var)

% get space
w = whos;

% return allocation
sz = w.bytes;
