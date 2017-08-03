function n = numel(obj, varargin)
% transio::numel  - return the number of underlying elements
%
% FORMAT:       N = numel(tio)
%
% Input fields:
%
%       tio         transio object
%
% Output fields:
%
%       N           number of elements

% Version:  v0.9a
% Build:    11052602
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

% return number of elements
if nargin == 1
    n = prod(obj.DataDims);
elseif nargin > 2
    error( ...
        'transio:BadArgument', ...
        'Indexing for numel only with (up to) one index.' ...
    );
else
    if ischar(varargin{1}) && ...
        strcmp(varargin{1}, ':')
        n = 1;
    elseif isa(varargin{1}, 'double') && ...
       ~any(isinf(varargin{1}(:)) | isnan(varargin{1}(:)) | varargin{1}(:) < 1 | varargin{1}(:) ~= fix(varargin{1}(:))) && ...
        all(varargin{1}(:) <= numel(obj.IOOffset))
        n = numel(varargin{1});
    elseif isa(varargin{1}, 'logical')
        n = sum(varargin{1}(:));
    else
        error( ...
            'transio:BadArgument', ...
            'Invalid numel-indexing expression.' ...
        );
    end
    if n < 1
        error( ...
            'transio:BadArgument', ...
            'Empty numel indexing not supported.' ...
        );
    end
    n = 1;
end
