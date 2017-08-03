function m = using(S, m)
% neuroelf::using  - make methods available in caller as function handles
%
% FORMAT:       using(neuroelf, 'METHOD')
%               using(neuroelf, {'METHOD', 'METHOD2', ...})
%
% No output fields.

% Version:  v0.9d
% Build:    14082113
% Date:     Aug-21 2014, 1:13 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2014, Jochen Weber
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

% valid
if nargin < 2
    m = 'all';
elseif (~iscell(m) && ...
    ~ischar(m)) || ...
    isempty(m)
    error( ...
        'neuroelf:BadArgument', ...
        'Bad or missing argument.' ...
    );
end

% get as struct and methods
neos = struct(S);
meth = neos(1).meth;

% for single char argument
if ischar(m)

    % set as single cell
    m = {m(:)'};

    % if all methods are requested
    if strcmpi(m{1}, 'all')

        % return methods if requested
        if nargout > 0
            m = neos.mets;
            return;
        end
        
        % otherwise, set method names as requested
        m = fieldnames(meth);
    end
end
for mc = 1:numel(m)
    if isfield(meth, m{mc})
        mh = meth.(m{mc}){1};
        assignin('caller', m{mc}, mh);
        m{mc} = mh;
    end
end
