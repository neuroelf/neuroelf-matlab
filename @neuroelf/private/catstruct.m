function cs = catstruct(varargin)
% catstruct  - concatenate structs with any fields
%
% FORMAT:       cs = catstruct(s1, s2 [, ... [, fields]])
%
% Input fields:
%
%       s1, s2      structs to concatenate
%       fields      list of fields (in desired order)
%
% Output fields:
%
%       cs          concatenated struct
%
% Note: if possible catstruct, by default, concatenates along the last
%       non-singleton dimension; if that is not possible the structs are
%       first squeezed into Sx1 arrays

% Version:  v0.9c
% Build:    11111615
% Date:     Nov-16 2011, 3:07 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2011, Jochen Weber
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
if nargin < 2 || ...
   ~isstruct(varargin{1}) || ...
   (~isstruct(varargin{2}) && ...
    (nargin > 2 && ...
     ~iscell(varargin{2})))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
narg = nargin;
if iscell(varargin{narg});
    fields = varargin{narg};
    narg = narg - 1;
else
    fields = fieldnames(varargin{1});
end
catd = 0;
stsize = ones(narg, 5);
fm = true;
for ac = 1:narg
    if ~isstruct(varargin{ac}) || ...
        ndims(varargin{ac}) > 5
        error( ...
            'neuroelf:BadArgument', ...
            'All arguments (but the last) must be of type struct and at most 5-D.' ...
        );
    end
    fn = fieldnames(varargin{ac});
    if numel(fn) ~= numel(fields) || ...
       ~all(strcmp(fn, fields))
        fields = uunion(fields, fn);
        fm = false;
    end
    stsize(ac, 1:ndims(varargin{ac})) = size(varargin{ac});
    if catd == 0 && ...
        max(stsize(ac, :)) > 1 && ...
        max(stsize(ac, :)) == prod(stsize(ac, :))
        catd = maxpos(stsize(ac, :));
    end
end
nf = numel(fields);
if catd == 0
    catd = findfirst(any(stsize > 1, 1), -1);
    if isempty(catd)
        catd = 1;
    end
end
if any(any(diff(stsize(:, [1:(catd-1), (catd+1):5])) ~= 0)) || ...
   (~fm && ...
    all(all(stsize(:, [1:(catd-1), (catd+1):5]) == 1)))
    lsq = true;
else
    lsq = false;
end

% if all fieldnames match, we can actually use cat!
if fm

    % with squeezing
    if lsq

        % first make a private copy
        args = varargin(1:narg);

        % squeeze them
        for ac = 1:narg
            args{ac} = args{ac}(:);
        end

        % cat
        cs = cat(1, args{:});

        % reshape?
        if catd > 1
            cs = reshape(cs, [ones(1, catd - 1), numel(cs)]);
        end

    % without squeezing
    else

        % cat in requested dimension
        cs = cat(catd, varargin{1:narg});
    end

    return;
end

% with squeezing
if lsq

    % first compute output size
    nout = sum(prod(stsize, 2));

    % then generate fields containing the data
    cs = cell(nf, nout);

    % initialize target counter
    tc = 1;

    % iterate over structs
    for ac = 1:narg

        % find the field mappings and size of struct
        fm = multimatch(fieldnames(varargin{ac}), fields);
        ns = numel(varargin{ac});

        % and assign in target
        cs(fm, tc:(tc+ns-1)) = struct2cell(varargin{ac}(:));

        % increase target counter
        tc = tc + ns;
    end

    % finally make a new struct
    cs = cell2struct(cs, fields, 1);

    % reshape
    if catd > 1
        cs = reshape(cs, [ones(1, catd - 1), numel(cs)]);
    end

% without squeezing
else

    % first compute output size
    sout = stsize(1, :);
    sout(catd) = sum(stsize(:, catd));

    % then generate fields containing the data
    cs = cell([nf, sout]);

    % initialize target counter
    tc = 1;

    % create assignment flag
    af = repmat({':'}, 1, 6);

    % iterate over structs
    for ac = 1:narg

        % find the field mappings and size of struct
        fm = multimatch(fieldnames(varargin{ac}), fields);
        ns = size(varargin{ac}, catd);

        % and assign in target
        af{1} = fm(:)';
        af{catd+1} = tc:(tc+ns-1);
        cs(af{:}) = struct2cell(varargin{ac});

        % increase target counter
        tc = tc + ns;
    end

    % finally make a new struct
    cs = cell2struct(cs, fields, 1);
end
