function [varargout] = seqtrans(s)
%SEQTRANS  Compute the transition probability matrix of sequence.
%   T = SEQTRANS(S) computes the transition matrix (fractional occurrence)
%   of element B following element A in a square matrix, such that the
%   frequency of [A,B] occurring within S is at element T(A,B).
%
%   [T, E] = SEQTRANS(S) also returns the identities of the elements (in
%   case the elements are not ordered from 1:X).

% Version:  v1.1
% Build:    16053121
% Date:     May-31 2016, 9:41 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2016, Jochen Weber
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
if nargin < 1 || isempty(s)
    error('neuroelf:general:missingArgument', 'Missing argument: s.');
end

% parse each cell separately
if iscell(s)
    varargout{1} = cell(size(s));
    if nargout > 1
        varargout{2} = cell(size(s));
    end
    try
        if nargout == 1
            for cc = 1:numel(s)
                varargout{1}{cc} = seqtrans(s{cc});
            end
        else
            for cc = 1:numel(s)
                [varargout{1}{cc}, varargout{2}{cc}] = seqtrans(s{cc});
            end
            if all(lsqueeze(cellfun(@isequal, varargout{2}, repmat(varargout{2}(1), size(varargout{2})))))
                varargout{2} = varargout{2}{1};
            end
        end
    catch ne_eo;
        rethrow(ne_eo);
    end
    return;
end

% check data type
if ~isnumeric(s) && ~ischar(s)
    error('neuroelf:general:badArgument', 'Sequence must be numeric or char.');
end

% keep track of input class
cs = class(s);
if strcmpi(cs, 'double')
    cs = [];
else
    eval(['cs=@' cs ';']);
end

% then make sure no invalid values in input
s = double(s(:));
if any(s ~= fix(s) | isinf(s))
    error('neuroelf:general:badArgument', 'Data must be integer.');
end

% all elements within range
e = unique(s);
if nargout > 1
    varargout{2} = e;
    if ~isempty(cs)
        varargout{2} = cs(e);
    end
end
ne = numel(e);
mins = e(1);
maxs = e(end);
if ne ~= (maxs + 1 - mins)
    so = zeros(1, maxs + max(0, 1 - mins));
    if mins > 0
        so(e) = 1:ne;
        s = so(s);
    else
        so(e - (mins - 1)) = 1:ne;
        s = so(s - (mins - 1));
    end
elseif mins ~= 1
    s = s - (mins - 1);
end

% generate histogram
t = reshape(histc(s(1:end-1) + ne .* (s(2:end) - 1), 1:(ne*ne)), ne, ne);

% normalize
varargout{1} = t ./ (sum(t, 2) * ones(1, ne));
