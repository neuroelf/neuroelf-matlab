function pscb = meanpsctrans(b, mb, mt)
% meanpsctrans  - perform PSC transformation after regression
%
% FORMAT:       pscb = meanpsctrans(b, mb [, mt])
%
% Input fields:
%
%       b           beta images
%       mb          mean confound beta image
%       mt          mask threshold (defaut: 0.25)
%
% Output fields:
%
%       pscb        PSC transformed beta image (masked)

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
if nargin < 2 || ...
   ~isnumeric(b) || ...
   ~isnumeric(mb)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
if nargin < 3 || ...
   ~isa(mt, 'double') || ...
    numel(mt) ~= 1 || ...
    isinf(mt) || ...
    isnan(mt) || ...
    mt < 0 || ...
    mt > 0.5
    mt = 0.25;
end
sb = size(b);
if prod(sb(1:end-1)) ~= numel(mb)
    error( ...
        'neuroelf:BadArgument', ...
        'Last dimension of b must be image dimension.' ...
    );
end

% reshape mean beta to make sure and get repmat argument
if numel(sb) == 2
    mb = reshape(mb, sb(1), 1);
    rma = [1, sb(2)];
else
    mb = reshape(mb, sb(1:end-1));
    rma = [ones(1, numel(sb) - 1), sb(end)];
end

% mask mean beta
mb = double(mb);
mmb = mt * mean(mb(mb > 0));
mb(mb < mmb) = Inf;

% compute PSC
pscb = repmat(100 ./ mb, rma) .* double(b);

% remove illegal values
pscb(isinf(pscb) | isnan(pscb)) = 0;
