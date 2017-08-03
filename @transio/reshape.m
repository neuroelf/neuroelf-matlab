function A = reshape(A, varargin)
% transio::reshape  - overloaded method

% Version:  v0.9c
% Build:    11052601
% Date:     May-10 2011, 11:49 PM EST
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

% check all arguments
if nargin < 2 || ...
    numel(struct(A)) ~= 1
    error( ...
        'transio:BadArgument', ...
        'Bad or missing argument.' ...
    );
end
nvargin = nargin - 1;
hasempty = 0;
wasempty = 0;
psize = 1;
if nargin == 2 && ...
    iscell(varargin{1}) && ...
   ~isempty(varargin{1})
    szin = varargin{1};
else
    szin = varargin;
end
for ac = 1:nvargin
    if (~isa(szin{ac}, 'double') && ...
        ~isa(szin{ac}, 'int32') && ...
        ~isa(szin{ac}, 'uint32')) || ...
        numel(szin{ac}) ~= max(size(szin{ac})) || ...
        any(isinf(szin{ac}) | isnan(szin{ac}) | szin{ac} < 0 | szin{ac} ~= fix(szin{ac}))
        error( ...
            'transio:BadArgument', ...
            'Bad size input argument %d.', ...
            ac + 1 ...
        );
    end
    if isempty(szin{ac})
        hasempty = hasempty + 1;
        wasempty = ac;
    else
        psize = psize * prod(double(szin{ac}));
        szin{ac} = szin{ac}(:)';
    end
end
sA = struct(A);
osz = prod(sA.DataDims);
if hasempty > 1
    error( ...
        'transio:BadArgument', ...
        'Too many empty size arguments.' ...
    );
elseif hasempty == 1
    szin{wasempty} = osz / psize;
    if szin{wasempty} ~= fix(szin{wasempty})
        error( ...
            'transio:BadArgument', ...
            'Cannot resolve [] size to integer dimension.' ...
        );
    end
elseif psize ~= osz
    error( ...
        'transio:BadArgument', ...
        'Number of elements must not change.' ...
    );
end

% set new size
nsize = cat(2, szin{:});
if iscell(sA.FileName) && ...
    nsize(end) ~= numel(sA.FileName)
    error( ...
        'transio:BadArgument', ...
        'Multi-file transio must be MxNx...xF in size, F = number of files.' ...
    );
end
sA.DataDims = nsize;
if numel(sA.DataDims) == 1
    sA.DataDims(2) = 1;
end

% clear any buffer
sA.IOBuffer = {{}, []};

% return as object
A = transio(0, 'makeobject', sA);
