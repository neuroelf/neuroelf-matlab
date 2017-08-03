function p = subsasgn(p, S, V)
% xprogress::subsref  - support for .Visible = assignments
%
% FORMAT:       pbar.Visible = VISFLAG;
%
% Input fields:
%
%       pbar        xprogress object

% Version:  v0.9a
% Build:    11062312
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
    numel(p) ~= 1 || ...
   ~isstruct(S) || ...
   ~isfield(S, 'type') || ...
   ~isfield(S, 'subs') || ...
    isempty(S) || ...
   ~strcmp(S(1).type, '.')
    error( ...
        'xprogress:BadCall', ...
        'Invalid call to xprogress::subsasgn.' ...
    );
end

% allow cell with one char call type
ssubs = S(1).subs;
if iscell(ssubs) && ...
    numel(ssubs) == 1
    ssubs = ssubs{1};
end
if ~ischar(ssubs) || ...
    isempty(ssubs)
    error( ...
        'xprogress:BadCall', ...
        'Invalid call to xprogress:subsasgn.' ...
    );
end

% what call
switch (lower(ssubs(:)'))
    case {'progress'}
        xprogress(p, V);
    case {'visible'}
        xprogress(p, 'setvisible', V);
    otherwise
        error( ...
            'xprogress:BadCall', ...
            'Invalid call to xprogress::subsasgn.' ...
        );
end
