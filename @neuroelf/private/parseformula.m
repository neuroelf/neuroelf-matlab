function f = parseformula(f, token, repl, dim, tr)
% parseformula  - parse a formula
%
% FORMAT:       f = parseformula(f, token, repl, dim, tr)
%
% Input fields:
%
%       f           formula
%       token       string to be matched
%       repl        replacement string
%       dim         concatenation dimension for output
%       tr          target indices (also used for available range)
%
% Output fields:
%
%       f           patched formula
%
% Note: this works for direct access and 1-level subfields using cat!

% Version:  v0.9b
% Build:    11050712
% Date:     Apr-09 2011, 11:08 PM EST
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
if nargin < 5 || ...
   ~ischar(f) || ...
   ~ischar(token) || ...
    isempty(token) || ...
   ~ischar(repl) || ...
    isempty(repl) || ...
    isempty(strfind(repl(:)', token(:)')) || ...
   ~isa(dim, 'double') || ...
    numel(dim) ~= 1 || ...
    isinf(dim) || ...
    isnan(dim) || ...
   ~any(1:63 == dim) || ...
   ~isa(tr, 'double') || ...
    isempty(tr) || ...
    any(isinf(tr(:)) | isnan(tr(:)) | tr(:) < 0 | tr(:) ~= fix(tr(:)))
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid combination of arguments.' ...
    );
end

% get good versions first
f = f(:)';
token = token(:)';
ftoken = reshape([token; token], 1, 2 * numel(token));
ftoken(1:2:end) = '\';
itoken = '([1-9]\d*)';
stoken = ['(\-?' itoken(2:end)];
repl = repl(:)';
tr = tr(:);
cats = sprintf('cat(%d,', dim);

% replace i1:s:i2 encounters (in reverse order!)
[rs, re, rt] = regexp(f, [ftoken itoken '\:' stoken '\:' itoken]);
for rc = numel(rt):-1:1

    % make sure all indices are valid
    i1 = str2double(f(rt{rc}(1, 1):rt{rc}(1, 2)));
    st = str2double(f(rt{rc}(2, 1):rt{rc}(2, 2)));
    i2 = str2double(f(rt{rc}(3, 1):rt{rc}(3, 2)));
    try
        ti = tr(i1:st:i2);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadArgument', ...
            'The following expression couldn''t be resolved: %s.', ...
            f(rs(rc):re(rc)) ...
        );
    end

    % and disallow empty expressions
    if isempty(ti)
        error( ...
            'neuroelf:BadArgument', ...
            'Empty expressions cannot be resolved!' ...
        );
    end

    % for one index, don't use cat
    if numel(ti) == 1
        f = [f(1:rs(rc)-1) ...
            strrep(repl, token, sprintf('%d', ti)) ...
            f(rt{rc}(end)+1:end)];
    else
        % make string
        ti = sprintf('%d,', ti);
        f = [f(1:rs(rc)-1) ...
            cats strrep(repl, token, ['[' ti(1:end-1) ']']) ')' ...
            f(rt{rc}(end)+1:end)];
    end
end

% replace i1:i2 encounters (in reverse order!)
[rs, re, rt] = regexp(f, [ftoken itoken '\:' itoken]);
for rc = numel(rt):-1:1

    % make sure all indices are valid
    i1 = str2double(f(rt{rc}(1, 1):rt{rc}(1, 2)));
    i2 = str2double(f(rt{rc}(2, 1):rt{rc}(2, 2)));

    % allow simple reversing!
    if i2 < i1
        i2o = i2;
        i2 = i1;
        i1 = i2o;
    end

    % test expression
    try
        ti = tr(i1:i2);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadArgument', ...
            'The following expression couldn''t be resolved: %s.', ...
            f(rs(rc):re(rc)) ...
        );
    end

    % and disallow empty expressions
    if isempty(ti)
        error( ...
            'neuroelf:BadArgument', ...
            'Empty expressions cannot be resolved!' ...
        );
    end

    % for one index, don't use cat
    if numel(ti) == 1
        f = [f(1:rs(rc)-1) ...
            strrep(repl, token, sprintf('%d', ti)) ...
            f(rt{rc}(end)+1:end)];
    else
        % make string
        ti = sprintf('%d,', ti);
        f = [f(1:rs(rc)-1) ...
            cats strrep(repl, token, ['[' ti(1:end-1) ']']) ')' ...
            f(rt{rc}(end)+1:end)];
    end
end

% replace i1 encounters (in reverse order!)
[rs, re, rt] = regexp(f, [ftoken itoken]);
for rc = numel(rt):-1:1

    % make sure all indices are valid
    i1 = str2double(f(rt{rc}(1, 1):rt{rc}(1, 2)));
    try
        ti = tr(i1);
    catch ne_eo;
        neuroelf_lasterr(ne_eo);
        error( ...
            'neuroelf:BadArgument', ...
            'The following expression couldn''t be resolved: %s.', ...
            f(rs(rc):re(rc)) ...
        );
    end

    % replace
    f = [f(1:rs(rc)-1) ...
        strrep(repl, token, sprintf('%d', ti)) ...
        f(rt{rc}(end)+1:end)];
end
