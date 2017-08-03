function s = emptystruct(fields, dim)
% emptystruct  - create an empty struct array with fields
%
% FORMAT:       s = emptystruct(fields [, dim])
%
% Input fields:
%
%       fields      fieldnames for struct s
%       dim         optional dim argument (used in repmat, default: [0, 0])
%
% Output fields:
%
%       s           struct (sized with dim) with fields
%
% Note: this uses cell2struct with a cell of size [1, 1, NUMEL(fields)]
%       to avoid a bug in Octave 3.x.

% Version:  v0.9c
% Build:    11050301
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

% argument check
if nargin < 1 || ...
   ~iscell(fields)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 2 || ...
   ~isa(dim, 'double') || ...
    numel(dim) < 2 || ...
    numel(dim) ~= size(dim, 2) || ...
    any(isinf(dim) | isnan(dim) | dim < 0 | dim >= 2^31 | dim ~= fix(dim))
    dim = [0, 0];
end
dim0 = (prod(dim) == 0);

% create struct
if numel(fields) == 0;
    s = struct;
else
    s = cell2struct(cell(numel(fields), 1), fields(:)', 1);
end
if dim0
    s(:) = [];
    if any(dim > 0)
        s = reshape(s, dim);
    end
else
    s = repmat(s, dim);
end
