function c = conjvalp(varargin)
% conjvalp  - conjugate pvalues (same sign, max value with sign)
%
% FORMAT:       c = conjvalp(v1, v2, ...)
%               c = conjvalp(nd, dim)
%
% Input fields:
%
%       v1, v2, ... values (or arrays) over which to compute conjunction
%       nd          N-dim data array
%       dim         1x1 dimension along which to run conjunction
%
% Output fields:
%
%       c           conjunction value (array)

% Version:  v0.9b
% Build:    13122114
% Date:     Apr-09 2011, 1:55 PM EST
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

% special case: array, dim
if nargin == 2 && ...
    isnumeric(varargin{1}) && ...
    isa(varargin{2}, 'double') && ...
    numel(varargin{2}) == 1 && ...
   ~isinf(varargin{2}) && ...
   ~isnan(varargin{2}) && ...
    any((1:ndims(varargin{1})) == varargin{2}) && ...
    size(varargin{1}, varargin{2}) > 1

    % deal into separate arrays
    va = cell(1, size(varargin{1}, varargin{2}));
    sra = repmat({':'}, ndims(varargin{1}));
    for vc = 1:numel(va)
        sra{varargin{2}} = vc;
        va{vc} = varargin{1}(sra{:});
    end

    % make call as usual
    c = conjvalp(va{:});

    % return
    return;
end

% argument check
if nargin < 2 || ...
   ~isnumeric(varargin{1}) || ...
   ~isnumeric(varargin{2}) || ...
    isempty(varargin{1}) || ...
    isempty(varargin{2}) || ...
   (~any(numel(varargin{1}) == [1, numel(varargin{2})]) && ...
    ~any(numel(varargin{2}) == [1, numel(varargin{1})]))
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% for multiple maps, start conjugating from the end
if nargin > 2
    xarg = cell(1, nargin);
    xarg(:) = varargin(:);
    try
        for c = nargin-1:-1:2
            xarg{c} = conjvalp(xarg{c}, xarg{c + 1});
        end
        c = conjvalp(xarg{1}, xarg{2});
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% perform conjugate (take sign of value 1 only if numel > 1!)
if numel(varargin{1}) > 1
    c = reshape( ...
        (sign(varargin{1}(:)) == sign(varargin{2}(:))) .* ...
        sign(varargin{1}(:)) .* ...
        max(abs(double(varargin{1}(:))), abs(double(varargin{2}(:)))), ...
        size(varargin{1}));
else
    c = reshape( ...
        (sign(varargin{1}(:)) == sign(varargin{2}(:))) .* ...
        sign(varargin{2}(:)) .* ...
        max(abs(double(varargin{1}(:))), abs(double(varargin{2}(:)))), ...
        size(varargin{2}));
end

% and make sure values are <= 1
c = max(-1, min(c, 1));
c(c == 0) = 1;
