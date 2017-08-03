function d = structdiff(s1, s2, af)
% structdiff  - compare two structures
%
% FORMAT:       diff = structdiff(struct1, struct2, allfields)
%
% Input fields:
%
%       struct1,2   1x1 structures
%       allfields   if given and true, diff contains all fieldnames
%
% Output fields:
%
%       diff        1x2 struct with different fields

% Version:  v0.9a
% Build:    10051716
% Date:     May-17 2010, 10:48 AM EST
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

% basic argument check
if nargin < 2 || ...
   ~isstruct(s1) || ...
    numel(s1) ~= 1 || ...
   ~isstruct(s2) || ...
    numel(s2) ~= 1
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid or missing argument(s).' ...
    );
end

% get fieldnames
f1 = fieldnames(s1);
f2 = fieldnames(s2);
if numel(f1) ~= numel(f2) || ...
   ~all(strcmp(f1, f2))
    of1 = setdiff(f1, f2);
    of2 = setdiff(f2, f1);
    fb = f1;
    for fc = 1:numel(of1)
        fb(strcmp(fb, of1{fc})) = [];
    end
else
    of1 = {};
    of2 = {};
    fb = f1;
end
ef = fb;

% create output struct
d = cell2struct(cell(2, 1, numel(of1) + numel(of2) + numel(fb)), ...
    [fb; of1; of2], 3);

% iterate over common fields
for fc = numel(fb):-1:1
    if ~isequal(s1.(fb{fc}), s2.(fb{fc}))
        d(1).fb{fc} = s1.(fb{fc});
        d(2).fb{fc} = s2.(fb{fc});
        ef(fc) = [];
    end
end

% remove equal fields
if nargin < 3 || ...
   ~islogical(af) || ...
    numel(af) ~= 1 || ...
   ~af
    d = rmfield(d, ef);
end
