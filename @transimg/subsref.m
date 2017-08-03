function fld = subsref(ti, S)
% transimg::subsref  - allow read-only access to properties

% Version:  v0.9c
% Build:    11050217
% Date:     May-2 2011, 5:21 AM EST
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

% check S
if ~isstruct(S) || ...
    isempty(S) || ...
   ~isfield(S, 'type') || ...
   ~isfield(S, 'subs') || ...
   ((~strcmp(S(1).type, '.') || ...
     ~any(strcmp(S(1).subs, ...
        {'Background', 'Handle', 'Height', 'IsRendered', 'Layer', 'Rendered', 'Width'}))) && ...
    ~strcmp(S(1).type, '()'))
    error( ...
        'transimg:InvalidCall', ...
        'Invalid subsref call to transimg object.' ...
    );
end

% global variables for storage
global tiobj ...
       tiobjlup;

% allow indexing to multi-objects
if strcmp(S(1).type, '()')
    try
        ti = subsref(struct(ti), S(1));
        S(1) = [];
    catch ne_eo;
        rethrow(ne_eo);
    end
    if isempty(S)
        fld = class(ti, 'transimg');
        return;
    end
end

% check arguments
try
    lup = find(tiobjlup == ti.L);
    if numel(lup) ~= 1
        error( ...
            'transimg:ObjectRemoved', ...
            'Object removed from global storage.' ...
        );
    end
catch ne_eo;
    rethrow(ne_eo);
end

% first get field
try
    fld = subsref(tiobj(lup), S(1));

    % then pass on if needed
    if numel(S) > 1
        fld = subsref(fld, S(2:end));
    end

% error handling
catch ne_eo;
    rethrow(ne_eo);
end
