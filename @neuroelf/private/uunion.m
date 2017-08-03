function c = uunion(a, b, varargin)
% uunion  - unsorted set union.
%
% FORMAT:       u = uunion(a, b [, varargin])
%
% Input fields:
%
%       a, b        input vectors of same class type
%
% Output fields:
%
%       u           output union (unsorted)
%
% Note: In contrast to union that always sorts the results,
%       uunion(A,B) will keep the original order of A and add
%       any additional values of B in the order they occur in B.
%       Otherwise, the same syntax additions apply as with UNION.
%
% See also union

% Version:  v1.0
% Build:    15121018
% Date:     Dec-10 2015, 6:12 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2015, Jochen Weber
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

% work around different implementations of MATLAB's union
persistent uu_type;
if isempty(uu_type) || ...
   ~islogical(uu_type)

    % create indices of union
    [testu, tu1, tu2] = union({'abd'}, {'abc'; 'abd'; 'abe'});
    uu_type = isempty(tu1);
end

% argument check
if nargin < 2 || ...
   ~strcmp(class(a), class(b))
    error( ...
        'neuroelf:BadArgument', ...
        'Two arguments of the same class must be given for uunion.' ...
    );
end

% prepare vars
rows = false;
cac  = 0;

% check third arg
if nargin > 2
    if ischar(varargin{end}) && ...
        strcmpi(varargin{end}(:)', 'rows')
        rows = true;
        carg = nargin - 3;
    else
        rows = false;
        carg = nargin - 2;
    end
    if carg > 0
        for cac = 1:carg
            if ~strcmp(class(a), class(varargin{cac}))
                cac = cac - 1;
                break;
            end
        end
    end
end

% what dim to use for concat
if ~rows
    mdim = length(size(a));
    cdim = find(size(a) > 1);
    if isempty(cdim)
        cdim = find(size(b) > 1);
    end
    if isempty(cdim)
        cdim = 1;
    else
        cdim = cdim(1);
    end
    rdim = [ones(1, cdim - 1), 0, ones(1, mdim - cdim)];
end

% try first union call
try
    if ~rows
        a = a(:);
        b = b(:);
        if uu_type
            [c, ib, ia] = union(b, a);
        else
            [c, ia, ib] = union(a, b);
        end
    else
        if uu_type
            [c, ib, ia] = union(b, a, 'rows');
        else
            [c, ia, ib] = union(a, b, 'rows');
        end
    end
catch ne_eo;
    rethrow(ne_eo);
end

% resort array
if ~rows
    c = [a(sort(ia)); b(sort(ib))];
    try
        rdim(cdim) = numel(c);
        c = reshape(c, rdim);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
    end
else
    if ischar(a)
        c = char([cellstr(a(sort(ia), :)); cellstr(b(sort(ib), :))]);
    else
        try
            c = [a(sort(ia), :); b(sort(ib), :)];
        catch ne_eo;
            rethrow(ne_eo);
        end
    end
end

% add more arrays ?
if cac > 0
    if rows
        rows = {'rows'};
    else
        rows = {};
    end
    for cacc = 1:cac
        try
            c = uunion(c, varargin{cacc}, rows{:});
        catch ne_eo;
            warning( ...
                'neuroelf:BadArgument', ...
                'Additional argument to uunion couldn''t be added: %s.', ...
                ne_eo.message ...
            );
        end
    end
end
