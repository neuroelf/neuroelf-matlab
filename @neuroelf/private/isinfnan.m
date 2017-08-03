function o = isinfnan(V, fa, fi, fn)
% isinfnan  - checking whether (any) element(s) are inf/nan
%
% FORMAT:       o = isinfnan(V, fa, fi, fn)
%
% Input fields:
%
%       V           N-D single/double matrix (for other numeric arrays, false)
%       fa          "any" flag (return scalar true/false; default: true)
%       fi          "inf" flag (test for infs, default: true)
%       fn          "nan" flag (test for nans, default: true)
%
% Output fields:
%
%       o           output
%
% Note: this is a MEX (compiled) function

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

% argument check
if nargin < 4 || ...
   ~islogical(fn) || ...
    numel(fn) ~= 1
    fn = true;
end
if nargin < 3 || ...
   ~islogical(fi) || ...
    numel(fi) ~= 1
    fi = true;
end
if nargin < 2 || ...
   ~islogical(fa) || ...
    numel(fa) ~= 1
    fa = true;
end
if nargin < 1 || ...
   ~isnumeric(V)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% for non single/double input
if ~isa(V, 'double') || ...
   ~isa(V, 'single')

    % return scalar output
    if fa
        o = false;
    else
        o = false(size(V));
    end
    return;
end

% scalar output
if fa

    % what tests
    o = false;

    % inf-test
    if fi
        o = o || any(isinf(V(:)));
    end

    % nan-test
    if fn
        o = o || any(isnan(V(:)));
    end

% full output
else

    % no test
    if ~fi && ...
       ~fn

        % return false
        o = false(size(V));

    % at least inf-test
    elseif fi

        % also nan-test
        if fn
            o = isinf(V) | isnan(V);

        % only inf-test
        else
            o = isinf(V);
        end

    % only nan-test
    else
        o = isnan(V);
    end
end
