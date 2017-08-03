function ss = selectstruct(si, clause, opts)
% selectstruct  - select from a struct with a clause
%
% FORMAT:       selected = selectstruct(s, clause [, opts])
%
% Input fields:
%
%       s           Nx1 structure with fields used in clause
%       clause      SQL-like selection criterion, supports
%                   'AND', 'OR', 'NOT', '=', 'LIKE'
%                   to test a char field, use single quotes!
%       opts        optional settings
%        .groupby   field (or list of fields) to group by, returns a cell
%        .print     flag, print out results in a table
%        .tstruct   equally sized Nx1 struct with fields to use in query
%
% Output fields:
%
%       ss          selected struct indices (either as struct or cell)
%
% Note: the syntax must be according to this template:
%
% 'SELECT [list-of-fields | *] FROM struct WHERE [criterion | 1];'

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

% allowed functions: SQL name -> {aggregate vs. scalar, Matlab function}
persistent ss_func;
if isempty(ss_func)
    ss_func = struct( ...
        'abs',   {{false, 'abs'}}, ...
        'avg',   {{ true, 'mean'}}, ...
        'count', {{ true, 'numel'}}, ...
        'len',   {{false, 'numel'}}, ...
        'max',   {{ true, 'max'}}, ...
        'min',   {{ true, 'min'}}, ...
        'round', {{false, 'round'}}, ...
        'sign',  {{false, 'sign'}}, ...
        'sum',   {{ true, 'sum'}} ...
    );
end

% argument check
if nargin < 2 || ...
   ~isstruct(si) || ...
   ~ischar(clause) || ...
    numel(clause) < 28
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
si = si(:);
ss = si;
cs = struct2cell(si);
ct = zeros(size(cs));
for fc = 1:numel(ct)
    if ischar(cs{fc})
        ct(fc) = 2;
    elseif ~iscell(cs{fc}) && ...
       ~isempty(cs{fc})
        ct(fc) = 1;
    end
end
if any(any(diff(double(ct ~= 1), 1, 2))) || ...
    any(any(diff(double(ct ~= 2), 1, 2)))
    error( ...
        'neuroelf:BadArgument', ...
        'Fields in structure must be consistent char or numeric.' ...
    );
end
tf = fieldnames(si);
clause = ddeblank(clause(:)');
while clause(end) == ';'
    clause(end) = '';
    clause = deblank(clause);
end
fpart = strfind(lower(clause), ' from ');
wpart = strfind(lower(clause), ' where ');
if numel(clause) < 28 || ...
    ~strcmpi(clause(1:7), 'select ') || ...
    numel(fpart) ~= 1 || ...
    numel(wpart) ~= 1 || ...
    wpart < fpart
    error( ...
        'neuroelf:BadArgument', ...
        'Bad clause argument.' ...
    );
end
if nargin < 3 || ...
   ~isstruct(opts) || ...
    numel(opts) ~= 1
    opts = struct;
end
if ~isfield(opts, 'groupby') || ...
   (~ischar(opts.groupby) && ...
    ~iscell(opts.groupby)) || ...
    isempty(opts.groupby)
    opts.groupby = [];
elseif ischar(opts.groupby)
    opts.groupby = {opts.groupby(:)'};
end
if ~isfield(opts, 'tstruct') || ...
   ~isstruct(opts.tstruct) || ...
    numel(opts.tstruct) ~= numel(si)
    opts.tstruct = [];
    ts = si;
else
    ts = opts.tstruct(:);
    cs = struct2cell(ts);
    ct = zeros(size(cs));
    for fc = 1:numel(ct)
        if ischar(cs{fc})
            ct(fc) = 2;
        elseif ~iscell(cs{fc}) && ...
           ~isempty(cs{fc})
            ct(fc) = 1;
        end
    end
    if any(diff(double(ct ~= 1), 1, 2)) || ...
        any(diff(double(ct ~= 2), 1, 2))
        error( ...
            'neuroelf:BadArgument', ...
            'Fields in structure must be consistent char or numeric.' ...
        );
    end
end
sf = fieldnames(ts);
for gc = numel(opts.groupby):-1:1
    opts.groupby{gc} = opts.groupby{gc}(:)';
    if ~any(strcmpi(opts.groupby{gc}, sf))
        opts.groupby(gc) = [];
    end
end

% get parts
spart = clause(8:fpart-1);
fpart = clause(fpart+6:wpart-1);
wpart = clause(wpart+7:end);

% check selection part, allow for '*'
spart = regexprep(spart, 'count\s*\(\s*\*\s*\)', ['count(' sf{1} ')']);
spart = strrep(spart, '*', sprintf('%s,', tf{:}));
spart = splittocellc(spart, ',', true, true);
spartas = spart;
spartfun = repmat({{}}, 1, numel(spart));
for sc = numel(spart):-1:1

    % discard empty fields
    spart{sc} = ddeblank(spart{sc});
    if isempty(spart{sc})
        spart(sc) = [];
        spartas(sc) = [];
        spartfun(sc) = [];
        continue;
    end

    % look for "as"
    [asname, aspart] = regexpi(spart{sc}, '\s+as\s+([a-z][a-z0-9_]*)$', 'tokens');
    if (any(spart{sc} == ' ') && ...
        isempty(aspart))
        error( ...
            'neuroelf:BadArgument', ...
            'Invalid syntax in fieldnames to select.' ...
        );
    end

    % split of as part and keep track
    if ~isempty(aspart)
        spart{sc} = spart{sc}(1:aspart(1)-1);
        spartas{sc} = asname{1}{1};
    end

    % check for functions
    [funpart, funmatch] = regexpi(spart{sc}, '^([a-z]+)\((.+)\)$', 'tokens');
    while ~isempty(funmatch)
        if ~any(strcmpi(funpart{1}{1}, fieldnames(ss_func)))
            error( ...
                'neuroelf:BadArgument', ...
                'Invalid function %s in select part.', ...
                funpart{1} ...
            );
        end
        spartfun{sc}{end+1} = funpart{1}{1};
        spart{sc} = ddeblank(funpart{1}{2});
        [funpart, funmatch] = regexpi(spart{sc}, '^([a-z]+)\((.+)\)$', 'tokens');
    end

    % now test whether remainder matches!
    if ~any(strcmpi(spart{sc}, sf))
        error( ...
            'neuroelf:BadArgument', ...
            'Unknown field %s or invalid syntax in select part.', ...
            spart{sc} ...
        );
    end
end
