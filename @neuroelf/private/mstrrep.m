function [text, wasp] = mstrrep(text, pats, reps, varargin)
% mstrrep  - replace multiple patterns (multi strrep)
%
% FORMAT:       text = mstrrep(text, patterns, replaces [, rx])
%
% Input fields:
%
%       text        string to be patched
%       patterns    1xN cell array with patterns
%       replaces    1xN (or 1x1) cell array with replacement strings
%       rx          if given and true, use regexprep
%
% Output fields:
%
%       text        patched string
%       waspatched  string altered flag (convenience only, calls strcmp)
%
% Note: regexprep is called with 'ignorecase' and 'tokenize' options.
%
% See also strrep, regexprep, strcmp

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
if nargin < 3 || ...
   ~ischar(text) || ...
   ~iscell(pats) || ...
    isempty(pats) || ...
   ~iscell(reps) || ...
    isempty(reps) || ...
    all(numel(reps) ~= [1, numel(pats)])
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid combination of arguments.' ...
    );
end
text = text(:)';
pats = pats(:)';
reps = reps(:)';
npat = length(pats);
if length(reps) == 1
    reps(2:npat) = reps;
end

% keep original state for comparison only if requested!
if nargout > 1
    texo = text;
end

% strrep or regexp
userx = false;
if nargin > 3 && ...
    numel(varargin{1}) == 1 && ...
   (isnumeric(varargin{1}) || ...
    islogical(varargin{1}))
    if varargin{1}
        userx = true;
        if mainver > 6
            rxargs = {'ignorecase'};
        else
            rxargs = {'ignorecase', 'tokenize'};
        end
    end
end

% what to do
try
	if userx
        for pc = 1:npat
            text = regexprep(text, pats{pc}, reps{pc}, rxargs{:});
        end
	else
        for pc = 1:npat
            text = strrep(text, pats{pc}, reps{pc});
        end
	end
catch ne_eo;
    neuroelf_lasterr(ne_eo);
    error( ...
        'neuroelf:BadArgument', ...
        'Invalid pattern/replacement combination for pattern %d.', ...
        pc ...
    );
end

% second argout
if nargout > 1
    wasp = logical(strcmp(texo, text));
end
